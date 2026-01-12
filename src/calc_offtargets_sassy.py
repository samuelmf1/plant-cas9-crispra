#!/usr/bin/env python3
"""
Sassy Off-Target Calculator
===========================

This module executes `sassy` to identify off-target sites for CRISPR guides
and aggregates the results, optionally filtering by gene family regions.

It handles:
1.  Sequence trimming and orientation normalization.
2.  Input validation.
3.  Execution of the `sassy` search tool.
4.  Robust parsing of sassy output (handling potential format variations).
5.  Aggregation of on-target and off-target counts.

Usage:
    python calc_offtargets_sassy.py -i <input.tsv> -g <genome.fa> -o <output.tsv>

Dependencies:
    - pandas
    - biopython
    - sassy (command line tool)
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Set, TextIO

import pandas as pd
from Bio.Seq import Seq

# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Calculate off-targets using sassy search.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV with 'rs3_context'")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA (can be .gz)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    parser.add_argument("--mismatch_intolerance", type=int, default=2, help="Mismatch tolerance (k) for sassy (default: 2)")
    parser.add_argument("--family_id", action="store_true", help="Enable family ID filtering/classification")
    parser.add_argument("--family_flank", type=int, default=0, help="Flank size to extend family regions (bp)")
    parser.add_argument("--sassy_log_file", help="Path to save raw sassy output for debugging")
    return parser.parse_args()

def get_trimmed_seq(seq: str, orientation: str, strand: str) -> str:
    """
    Extracts the 23bp core sequence (Spacer + PAM) and normalizes it to the
    Top Strand orientation required by Sassy.

    Args:
        seq: The full 30bp rs3_context sequence.
        orientation: Guide orientation ('FWD' or 'RVS').
        strand: Gene strand ('+' or '-').

    Returns:
        The normalized 23bp sequence string.
        Returns original unique sequence if length < 23.
    """
    if not isinstance(seq, str) or len(seq) < 23:
        return str(seq)
    
    # rs3_context format: [4bp prefix] [20bp Spacer] [3bp PAM] [3bp suffix]
    # We strip the 4bp prefix and 3bp suffix to get the 23bp core.
    core = seq[4:-3]
    
    # Normalization Logic:
    # Sassy requires Top Strand query to find hits reliably.
    # (+ FWD) -> core is Top. Keep.
    # (+ RVS) -> core is Bottom. Need RC.
    # (- FWD) -> core is Bottom. Need RC.
    # (- RVS) -> core is Top. Keep.
    
    needs_rc = (strand == '+' and orientation == 'RVS') or \
               (strand == '-' and orientation == 'FWD')
        
    if needs_rc:
        return str(Seq(core).reverse_complement())
        
    return core

def parse_chrom_from_text_id(text_id: str) -> str:
    """
    Extracts chromosome name from the FASTA header ID.
    
    Expected format examples:
        - dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF
        - chromosome:TAIR10:1
    """
    # Regex to extract '1' after chromosome:TAIR10: or similar
    m = re.search(r'chromosome:[^:]+:([^:]+):', text_id)
    if m:
        return m.group(1)
    
    # Fallback: simple split if format differs
    if ':' in text_id:
        return text_id.split(':')[0]
        
    return text_id

def load_family_regions(df: pd.DataFrame, flank: int = 0) -> Dict[str, List[Tuple[str, int, int]]]:
    """
    Loads family regions into a dictionary for fast lookup.

    Args:
        df: DataFrame containing 'family_id', 'chrom', 'start', 'end'.
        flank: Number of base pairs to extend the region on both sides.

    Returns:
        Dict mapping family_id -> List of (chrom, start, end) tuples.
    """
    fam_map = defaultdict(list)
    if 'family_id' not in df.columns:
        logger.warning("--family_id requested but 'family_id' column missing. Skipping family logic.")
        return fam_map
    
    # Filter for valid regions
    required_cols = ['family_id', 'chrom', 'start', 'end']
    if missing := [c for c in required_cols if c not in df.columns]:
        logger.warning(f"Missing columns for family logic: {missing}")
        return fam_map

    valid_mask = df[required_cols].notna().all(axis=1)
    subset = df.loc[valid_mask, required_cols]
    
    for row in subset.itertuples(index=False):
        fid = row.family_id
        chrom = str(row.chrom)
        try:
            s, e = int(row.start), int(row.end)
            if s > e: s, e = e, s
            # Apply flank
            s = max(0, s - flank)
            e = e + flank
            fam_map[fid].append((chrom, s, e))
        except ValueError:
            continue
            
    return fam_map

def check_overlap(chrom: str, start: int, end: int, regions: List[Tuple[str, int, int]]) -> bool:
    """
    Checks if a target region overlaps with any predefined family regions.
    
    Args:
        chrom: Target chromosome.
        start: Target start coordinate.
        end: Target end coordinate.
        regions: List of (chrom, start, end) tuples for the family.
    """
    target_chrom = str(chrom)
    for r_chrom, r_start, r_end in regions:
        if target_chrom == str(r_chrom):
            # Overlap condition: max(start1, start2) <= min(end1, end2)
            if max(start, r_start) <= min(end, r_end):
                return True
    return False

def run_sassy_and_count(
    cmd: List[str], 
    pat_to_indices: Dict[str, List[int]], 
    id_to_pat: Dict[int, str], 
    df_fams: Optional[pd.Series], 
    fam_regions: Dict[str, List[Tuple[str, int, int]]], 
    use_family: bool, 
    sassy_log_file: Optional[str] = None
) -> Tuple[Dict, Dict, Dict]:
    """
    Executes sassy subprocess, streams output, and aggregates counts.

    Returns:
        Tuple of (counts_total, counts_in_family, counts_out_family)
        Each is a dict of dicts: {row_idx: {mismatch_count: hits}}
    """
    # Initialize count structures: row_idx -> mismatch_k -> count
    counts_total = defaultdict(lambda: defaultdict(int))
    counts_in = defaultdict(lambda: defaultdict(int))
    counts_out = defaultdict(lambda: defaultdict(int))
    
    logger.info(f"Running Sassy: {' '.join(cmd)}")
    logger.info("Streaming results...")

    # Pre-fetch family IDs for efficiency
    row_fams = df_fams.to_dict() if use_family and df_fams is not None else {}
    
    # Open optional debug log
    log_handle: Optional[TextIO] = None
    if sassy_log_file:
        try:
            log_handle = open(sassy_log_file, 'w')
        except OSError as e:
            logger.warning(f"Could not open sassy log file {sassy_log_file}: {e}")

    try:
        proc = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True, 
            bufsize=100000
        )
        
        col_map = {}
        header_parsed = False
        
        if not proc.stdout:
            raise RuntimeError("Failed to access sassy stdout.")

        for line in proc.stdout:
            line = line.strip()
            if not line: continue
            
            parts = line.split('\t')
            
            # 1. Parse Header
            if not header_parsed:
                if parts[0] == "pat_id":
                    col_map = {c: i for i, c in enumerate(parts)}
                    header_parsed = True
                    if log_handle: 
                        log_handle.write(line + "\tsequence\n")
                    continue
            
            # 2. Parse Hit
            try:
                # Default indices based on standard sassy output
                # pat_id(0), text_id(1), cost(2), strand(3), start(4), end(5)
                p_idx = col_map.get('pat_id', 0)
                c_idx = col_map.get('cost', 2)
                t_idx = col_map.get('text_id', 1)
                s_idx = col_map.get('start', 4)
                e_idx = col_map.get('end', 5)
                
                # Heuristic: Detect 9-column shifted output format
                # If column 3 looks like an integer cost, use shifted indices
                if len(parts) >= 9:
                    try:
                        int(parts[3])
                        # Apply shift: text_id(2), cost(3), strand(4), start(5), end(6)
                        c_idx, t_idx, s_idx, e_idx = 3, 2, 5, 6
                    except ValueError:
                        pass # Not the shifted format

                # Parse Pattern ID (handle "pattern 1")
                p_val = parts[p_idx]
                if p_val.startswith("pattern"):
                    pat_id = int(p_val.split()[-1])
                else:
                    pat_id = int(p_val)
                
                # Retrieve Sequence
                seq = id_to_pat.get(pat_id)
                if log_handle:
                    log_handle.write(line + f"\t{seq or 'UNKNOWN'}\n")
                
                if not seq: continue
                
                cost = int(parts[c_idx])
                row_indices = pat_to_indices[seq]

                # 3. Update Counts
                if not use_family:
                    for r in row_indices:
                        counts_total[r][cost] += 1
                else:
                    # Optimize family checks by grouping rows
                    text_id = parts[t_idx]
                    chrom = parse_chrom_from_text_id(text_id)
                    start = int(parts[s_idx])
                    end = int(parts[e_idx])
                    
                    fams_to_check = set()
                    rows_by_fam = defaultdict(list)
                    
                    for r in row_indices:
                        counts_total[r][cost] += 1
                        fid = row_fams.get(r)
                        if fid:
                            rows_by_fam[fid].append(r)
                            fams_to_check.add(fid)
                        else:
                            # Rows without family ID don't get in/out increments
                            # per standard logic "in vs out" applies only when family is defined.
                            pass

                    # Perform spatial check once per unique family involved
                    for fid in fams_to_check:
                        is_in = check_overlap(chrom, start, end, fam_regions.get(fid, []))
                        for r in rows_by_fam[fid]:
                            if is_in:
                                counts_in[r][cost] += 1
                            else:
                                counts_out[r][cost] += 1
                                
            except (ValueError, IndexError):
                # Skip malformed lines
                continue
        
        # Ensure subprocess finished correctly
        proc.wait()
        if proc.returncode != 0:
            stderr = proc.stderr.read() if proc.stderr else "Unknown error"
            logger.error(f"Sassy failed: {stderr}")
            sys.exit(1)
            
    finally:
        if log_handle: log_handle.close()
        
    return counts_total, counts_in, counts_out

def main():
    args = parse_args()
    
    logger.info(f"Reading input: {args.input}")
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        logger.error(f"Failed to read input file: {e}")
        sys.exit(1)
    
    # Validation
    required_cols = ['rs3_context', 'Orientation']
    if not all(c in df.columns for c in required_cols):
        logger.error(f"Input must contain columns: {required_cols}")
        sys.exit(1)
        
    if 'strand' not in df.columns:
        logger.warning("'strand' column missing. Defaulting to '+'. Off-target discovery for negative strand guides may be inaccurate.")
    
    # 1. Prepare Sequences
    logger.info("Preparing sequences...")
    trim_func = lambda row: get_trimmed_seq(
        row['rs3_context'], 
        row['Orientation'], 
        str(row.get('strand', '+'))
    )
    df['__trimmed_seq'] = df.apply(trim_func, axis=1)
    
    unique_patterns = [p for p in df['__trimmed_seq'].unique() if isinstance(p, str) and len(p) > 0]
    logger.info(f"Found {len(unique_patterns)} unique patterns from {len(df)} rows.")
    
    # Maps for lookup
    # Sassy pattern IDs are 1-based
    pat_to_id = {pat: i+1 for i, pat in enumerate(unique_patterns)}
    id_to_pat = {i+1: pat for i, pat in enumerate(unique_patterns)}
    
    pat_to_indices = defaultdict(list)
    for idx, seq in zip(df.index, df['__trimmed_seq']):
        if seq in pat_to_id:
            pat_to_indices[seq].append(idx)
            
    # 2. Load Family Regions
    use_family = args.family_id
    fam_regions = {}
    df_fams = None
    if use_family:
        if 'family_id' in df.columns:
            logger.info(f"Loading family regions (flank={args.family_flank}bp)...")
            fam_regions = load_family_regions(df, flank=args.family_flank)
            df_fams = df['family_id']
        else:
            logger.warning("--family_id specified but column missing. Disabling family logic.")
            use_family = False

    # 3. Create Temporary Pattern File
    tmp_pat_path = "sassy_patterns.tmp"
    with open(tmp_pat_path, 'w') as f:
        for pat in unique_patterns:
            f.write(f"{pat}\n")
            
    # 4. Run Sassy
    cmd = ["sassy", "search", "-k", str(args.mismatch_intolerance), "-l", tmp_pat_path, args.genome]
    
    try:
        counts_total, counts_in, counts_out = run_sassy_and_count(
            cmd, pat_to_indices, id_to_pat, df_fams, fam_regions, use_family, args.sassy_log_file
        )
    finally:
        if os.path.exists(tmp_pat_path):
            os.remove(tmp_pat_path)
            
    # 5. Update DataFrame
    logger.info("Aggregating results...")
    
    max_k = args.mismatch_intolerance
    
    # Use dictionary comprehension for bulk column creation
    new_cols = {}
    for k in range(max_k + 1):
        # Total counts
        new_cols[f'offtarget_{k}'] = [counts_total[i][k] for i in df.index]
        
        if use_family:
            new_cols[f'offtarget_{k}_in_family'] = [counts_in[i][k] for i in df.index]
            new_cols[f'offtarget_{k}_out_family'] = [counts_out[i][k] for i in df.index]
            
    # Assign all new columns at once (efficient)
    df = df.assign(**new_cols)

    # 6. Quality Checks
    zero_hits = (df['offtarget_0'] == 0).sum()
    if zero_hits > 0:
        logger.warning(f"{zero_hits} guides have 0 perfect matches! Verify genome/guide compatibility.")
        
    # Cleanup
    if '__trimmed_seq' in df.columns:
        df.drop(columns=['__trimmed_seq'], inplace=True)
        
    logger.info(f"Saving results to {args.output}")
    df.to_csv(args.output, sep='\t', index=False)
    logger.info("Done.")

if __name__ == "__main__":
    main()
