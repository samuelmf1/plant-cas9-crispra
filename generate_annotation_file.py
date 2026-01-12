#!/usr/bin/env python3
"""
Generate Annotation File
========================

This script aggregates various gene annotation files into a single CSV for downstream analysis.
It merges family data, HDR/DDR classification, pseudogene status, and expression metrics.

Usage:
    python generate_annotation_file.py --data-dir data --output data/annotation.csv

Dependencies:
    - pandas
    - numpy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Set, Optional

import pandas as pd
import numpy as np

# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)

# Default Filenames
FILES = {
    "main": "Final_Combined_250625.txt",
    "family": "Ara_family.txt",
    "hdr": "HDR_total.txt",
    "ddr": "DNA_repair_total.txt",
    "pseudogenes": "Pseudogenes.txt",
    "hypothetical": "Hypothetical_proteins.txt"
}

def check_file(path: Path) -> Path:
    """Verifies that a file exists."""
    if not path.exists():
        logger.error(f"Required file not found: {path}")
        sys.exit(1)
    return path

def load_gene_set(path: Path, col: str = "Gene", **kwargs) -> Set[str]:
    """Loads a unique set of genes from a file."""
    try:
        df = pd.read_table(path, **kwargs)
        if df.empty:
            logger.warning(f"File {path.name} is empty.")
            return set()
        
        # Handle cases where column name might be missing or different in raw file
        # 'names' kwarg usually handles this, but we verify:
        if col not in df.columns:
            # Fallback: use first column
            genes = df.iloc[:, 0].unique()
        else:
            genes = df[col].unique()
            
        logger.info(f"Loaded {len(genes)} unique genes from {path.name}")
        return set(genes)
    except Exception as e:
        logger.error(f"Failed to load {path.name}: {e}")
        sys.exit(1)

def process_annotations(data_dir: Path, output_file: Path) -> None:
    """
    Main processing logic to merge annotation data.
    
    Args:
        data_dir (Path): Directory containing source text files.
        output_file (Path): Path to save the final CSV.
    """
    # 1. Load Main Data
    main_path = check_file(data_dir / FILES["main"])
    logger.info(f"Loading main annotation file: {main_path}")
    main_df = pd.read_table(main_path)
    
    initial_count = len(main_df)
    logger.info(f"Initial gene count: {initial_count}")

    # 2. Load and Merge Family Data
    fam_path = check_file(data_dir / FILES["family"])
    logger.info(f"Merging family data from: {fam_path}")
    fam_df = pd.read_table(fam_path).rename(columns={'Gene_ID': 'Gene'})
    
    main_df = main_df.merge(fam_df, on='Gene', how='left')
    
    # 3. Load Gene Sets (HDR, DDR, Pseudogenes)
    hdr_genes = load_gene_set(data_dir / FILES["hdr"], names=['Gene'])
    ddr_genes = load_gene_set(data_dir / FILES["ddr"], names=['Gene'])
    
    # Pseudogenes file has a header line to skip? Original code: skiprows=1, names=['Gene']
    # Let's verify file format. Assuming original code was correct.
    pseudo_genes = load_gene_set(data_dir / FILES["pseudogenes"], skiprows=1, names=['Gene'])

    # 4. Annotate Gene Types
    main_df['is_hdr_gene'] = main_df['Gene'].isin(hdr_genes)
    main_df['is_ddr_gene'] = main_df['Gene'].isin(ddr_genes)
    main_df['is_repair_gene'] = main_df['is_hdr_gene'] | main_df['is_ddr_gene']
    main_df['is_pseudogene'] = main_df['Gene'].isin(pseudo_genes)
    
    # 5. Hypothetical Proteins
    hypo_path = check_file(data_dir / FILES["hypothetical"])
    logger.info(f"Loading hypothetical proteins from: {hypo_path}")
    hypo_df = pd.read_table(hypo_path).rename(columns={'ID': 'Gene'})
    hypo_df['is_hypothetical_protein'] = True
    
    main_df = main_df.merge(hypo_df, on='Gene', how='left')
    main_df['is_hypothetical_protein'] = main_df['is_hypothetical_protein'].fillna(False)

    # 6. Fill Missing Values & Metrics
    # Using Infinity for easier filtering logic (original intent), though NaNs might be safer for stats.
    # Maintaining original logic for consistency.
    main_df['AA_length'] = main_df['AA_length'].fillna(np.inf)
    main_df['Max_expression_score'] = main_df['Max_expression_score'].fillna(np.inf)

    # 7. Apply Filtering Flag
    # Criteria: (Short AA < 100 AND Zero Expression) OR Pseudogene OR Very Short AA < 10
    filter_mask = (
        ((main_df['AA_length'] < 100) & (main_df['Max_expression_score'] == 0)) |
        (main_df['is_pseudogene']) |
        (main_df['AA_length'] < 10)
    )
    main_df['filt_shortnoexpAA_pseudogenes'] = filter_mask

    filtered_count = filter_mask.sum()
    logger.info(f"Flagged {filtered_count} genes as 'filt_shortnoexpAA_pseudogenes' ({filtered_count/len(main_df):.1%})")

    # 8. Save
    logger.info(f"Saving annotations to {output_file}")
    main_df.to_csv(output_file, index=False)
    logger.info("Done.")

def main():
    parser = argparse.ArgumentParser(description="Generate composite gene annotation file.")
    parser.add_argument(
        "--data-dir", 
        type=Path, 
        default=Path("data"), 
        help="Directory containing source annotation text files."
    )
    parser.add_argument(
        "--output", 
        type=Path, 
        default=None, 
        help="Output CSV path. Defaults to <data-dir>/annotation.csv"
    )

    args = parser.parse_args()

    # Determine output path
    if args.output is None:
        args.output = args.data_dir / "annotation.csv"

    # Validate Input Directory
    if not args.data_dir.is_dir():
        logger.error(f"Data directory does not exist: {args.data_dir}")
        sys.exit(1)

    process_annotations(args.data_dir, args.output)

if __name__ == "__main__":
    main()