#!/usr/bin/env python3

import os
import re
import sys
import gzip
import argparse
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from typing import Tuple, List, Dict

class GenomeSequenceExtractor:
    """Loads genome and extracts 5'->3' promoter windows relative to TSS."""
    def __init__(self, genome_file: str):
        if not os.path.exists(genome_file):
            raise FileNotFoundError(f"Genome file not found: {genome_file}")
        
        open_func = gzip.open if genome_file.endswith(".gz") else open
        with open_func(genome_file, "rt") as handle:
            self.genome = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

        self.chrom_map = {}
        for name in self.genome.keys():
            base = re.sub(r"^Chr|^chr", "", name)
            for prefix in ("", "chr", "Chr"):
                self.chrom_map[f"{prefix}{base}"] = name

    def _resolve_chrom(self, chrom: str) -> str:
        c = str(chrom)
        if c in self.genome: return c
        if c in self.chrom_map: return self.chrom_map[c]
        raise ValueError(f"Chromosome {chrom} not found.")

    def extract_promoter(self, chrom: str, tss: int, strand: str, far_bp: int, close_bp: int) -> Tuple[str, int, int]:
        """Returns (Sequence, Start_Coord, End_Coord) relative to genomic Top strand."""
        chrom = self._resolve_chrom(chrom)
        seq_record = self.genome[chrom]
        
        if strand == "+":
            # Genomic Top is Coding. Promoter is upstream (lower coordinates).
            p_start = tss - far_bp
            p_end = tss - close_bp - 1
            p_seq = seq_record.seq[p_start - 1 : p_end]
        else:
            # Genomic Bottom is Coding. Promoter is upstream (higher coordinates).
            p_start = tss + close_bp + 1
            p_end = tss + far_bp
            # RC to get 5'->3' relative to transcription direction
            p_seq = seq_record.seq[p_start - 1 : p_end].reverse_complement()
        
        return str(p_seq), p_start, p_end

def find_guides(seq: str) -> List[Dict]:
    """
    Scans the coding-strand sequence for SpCas9 guides.
    Captures 30-mer context for Rule Set 3.
    """
    results = []
    # We start at index 4 and end 3bp early to allow for the RS3 context window
    for i in range(4, len(seq) - 22 - 3):
        window = seq[i:i+23].upper()
        
        # Capture 30-mer context: [4bp prefix][23bp guide][3bp suffix]
        context_30m = seq[i-4 : i+23+3].upper()
        
        if window.endswith("GG"):
            results.append({
                "Spacer": window[:20], "PAM": window[20:], 
                "Orientation": "FWD", "targets": "nc",
                "GC": round((window[:20].count("G")+window[:20].count("C"))/20, 2),
                "polyT": "TTTT" in window[:20],
                "rs3_context": context_30m,
                "relative_pos": i
            })
            
        if window.startswith("CC"):
            # For RVS, we RC the 30-mer context so RS3 sees a standard NGG
            # Top RVS context: [3bp suffix_rc] [Spacer_rc] [PAM_rc] [4bp prefix_rc]
            # To get this, we need [i-3 : i+27] on Top strand
            # Check bounds
            if i - 3 < 0 or i + 27 > len(seq):
                continue

            raw_rvs = seq[i-3 : i+27].upper()
            rc_context = str(Seq(raw_rvs).reverse_complement())
            
            results.append({
                "Spacer": rc_context[4:24], "PAM": rc_context[24:27], 
                "Orientation": "RVS", "targets": "c",
                "GC": round((rc_context[4:24].count("G")+rc_context[4:24].count("C"))/20, 2),
                "polyT": "TTTT" in rc_context[4:24],
                "rs3_context": rc_context,
                "relative_pos": i
            })
    return results

def main():
    parser = argparse.ArgumentParser(description="Extract promoters and design guides.")
    parser.add_argument("--genome", help="Genome FASTA")
    parser.add_argument("--input", help="CSV/TSV with 'chrom', 'tss', 'strand'")
    parser.add_argument("--output_prefix", help="Output prefix")
    parser.add_argument("--guides", action="store_true", help="Extract and output guides")
    parser.add_argument("--far-bp", type=int, default=250)
    parser.add_argument("--close-bp", type=int, default=25)
    args = parser.parse_args()

    if not args.genome or not args.input:
        print("[INFO] Run with --help for options. Running internal tests...")
        return

    extractor = GenomeSequenceExtractor(args.genome)
    df = pd.read_csv(args.input, sep=None, engine="python")
    lowered_cols = {c.lower(): c for c in df.columns}
    
    chr_col = lowered_cols.get('chrom') or lowered_cols.get('chr') or lowered_cols.get('chromosome')
    tss_col = lowered_cols.get('tss') or lowered_cols.get('start') or lowered_cols.get('tx_start')
    strand_col = lowered_cols.get('strand')
    
    meta_keys = ['start', 'end', 'family_id', 'gene']
    found_meta = {k: lowered_cols[k] for k in meta_keys if k in lowered_cols}

    all_guides = []
    with open(f"{args.output_prefix}.fa", "w") as fa_out:
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
            try:
                strand = str(row[strand_col]).strip()
                p_seq, p_start, p_end = extractor.extract_promoter(
                    row[chr_col], int(row[tss_col]), strand, args.far_bp, args.close_bp
                )

                # Header includes coordinates and gene metadata
                header_parts = [f"{row[chr_col]}:{p_start}-{p_end}({strand})"]
                for k in ['gene', 'family_id']:
                    if k in found_meta: header_parts.append(f"{k}:{row[found_meta[k]]}")
                fa_out.write(f">{ '|'.join(header_parts) }\n{p_seq}\n")
                
                if args.guides:
                    guides = find_guides(p_seq)
                    for g in guides:
                        rel_pos = g['relative_pos']
                        if strand == '+':
                            # Top strand extraction: start + rel_pos
                            # p_start is genomic coord of start of window
                            g_coord = p_start + rel_pos
                        else:
                            # Bottom strand extraction (RC)
                            # p_seq starts from p_end (genomic high) going backwards
                            g_coord = p_end - rel_pos

                        metadata = {k: row[v] for k, v in found_meta.items()}
                        g.update({
                            'chrom': row[chr_col], 
                            'strand': strand,
                            'tss_coord': row[tss_col], 
                            'coord': g_coord,
                            **metadata
                        })
                        all_guides.append(g)
            except Exception:
                continue

    if args.guides and all_guides:
        out_df = pd.DataFrame(all_guides)
        order = ['gene', 'strand', 'targets', 'Orientation', 'Spacer', 'PAM']
        cols = [c for c in order if c in out_df.columns] + [c for c in out_df.columns if c not in order]
        out_df[cols].to_csv(f"{args.output_prefix}_guides.tsv", sep="\t", index=False)
        print(f"\n[DONE] Saved to {args.output_prefix}_guides.tsv")
    elif args.guides:
        print("\n[WARN] No guides searched for. Specify --guides to do so.")


if __name__ == "__main__":
    main()
