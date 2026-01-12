#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
from rs3.seq import predict_seq

def run_prioritization(input_file, output_file, rs3_tracr='Hsu2013'):
    # Read TSV
    df = pd.read_csv(input_file, sep='\t')
    
    # 1. Validation & Hard Cleaning
    required = ['Spacer', 'targets', 'GC', 'tss_coord', 'rs3_context']
    for col in required:
        if col not in df.columns:
            sys.exit(f"[ERROR] Missing column '{col}'.")

    initial_len = len(df)
    
    # Drop any row where essential sequences are NaN
    df = df.dropna(subset=['rs3_context', 'Spacer'])
    
    # Ensure rs3_context is exactly 30nt AND Spacer is exactly 20nt
    # This prevents the Bio.SeqUtils IndexError
    df = df[
        (df['rs3_context'].str.len() == 30) & 
        (df['Spacer'].str.len() == 20)
    ].copy()
    
    # Ensure sequences contain only valid IUPAC characters (A, T, C, G)
    # rs3 and melting temp calculations will crash on 'N' or weird chars
    df = df[df['rs3_context'].str.contains('^[ATCGatcg]+$', regex=True)]
    
    dropped = initial_len - len(df)
    if dropped > 0:
        print(f"[WARN] Dropped {dropped} rows due to empty sequences, invalid lengths, or non-ATCG characters.")

    if df.empty:
        sys.exit("[ERROR] No valid sequences remaining to process.")

    # 2. Execute Rule Set 3
    print(f"[INFO] Calculating Rule Set 3 scores for {len(df):,} guides...")
    try:
        # Rule Set 3 handles the biochemical efficiency
        df['rs3_score'] = predict_seq(df['rs3_context'].tolist(), sequence_tracr=rs3_tracr)
    except Exception as e:
        sys.exit(f"[ERROR] RS3 Prediction failed: {e}")
    
    # 3. Apply CRISPR-Act3.0 Rice Rules
    # Feature A: High-Activity Zone (-25 to -200 bp)
    # Priority: 1. Exact coord (if available from extractor), 2. Midpoint of promoter window
    if 'tss_coord' in df.columns and 'coord' in df.columns:
        df['dist_to_tss'] = (df['tss_coord'] - df['coord']).abs()
        df['is_best_zone'] = df['dist_to_tss'].between(25, 200)
    elif 'prom_start' in df.columns and 'prom_end' in df.columns:
        # Use midpoint of the promoter window to estimate distance if exact guide pos isn't mapped
        df['dist_to_tss'] = (df['tss_coord'] - ((df['prom_start'] + df['prom_end'])/2)).abs()
        df['is_best_zone'] = df['dist_to_tss'].between(25, 200)
    else:
        df['is_best_zone'] = True 

    # Feature B: Non-coding Preference (Top Priority in Rice Act3.0 paper)
    # Extractor outputs 'nc' for noncoding/target-strand
    df['is_noncoding'] = df['targets'].isin(['noncoding', 'nc'])

    # Feature C: GC Content (45% - 60% per paper)
    df['is_best_gc'] = df['GC'].between(0.45, 0.60)

    # 4. Off-Target Safety Check
    # Check for any columns ending in '_out_family' and sum them.
    # If sum > 0, then is_offtarget_safe = False
    
    out_fam_cols = [c for c in df.columns if c.endswith('_out_family')]
    if out_fam_cols:
        # Sum across these columns for each row
        # FillNa with 0 just in case
        total_out = df[out_fam_cols].fillna(0).sum(axis=1)
        df['num_out_family'] = total_out
        df['is_offtarget_safe'] = total_out == 0
    else:
        # If columns missing, assume safe or handle as warning?
        # Assuming safe if logic wasn't run or columns not passed (e.g. no families defined)
        df['is_offtarget_safe'] = True

    # 5. Composite Ranking Score (0 to 100)
    # Weighted by Act3.0 priorities: Strand (50), Position (30), GC (20)
    df['plant_act_score'] = (
        (df['is_noncoding'].astype(int) * 50) + 
        (df['is_best_zone'].astype(int) * 30) + 
        (df['is_best_gc'].astype(int) * 20)
    )

    # Final Sort: Safety first, then Biology, then Biochemistry
    df = df.sort_values(by=['is_offtarget_safe', 'plant_act_score', 'rs3_score'], ascending=False)

    # 6. Tiering
    df['Tier'] = "Tier 3 (Average)"
    df.loc[df['plant_act_score'] >= 50, 'Tier'] = "Tier 2 (Good)"
    df.loc[df['plant_act_score'] >= 80, 'Tier'] = "Tier 1 (Elite)"
    
    # Downgrade unsafe guides
    df.loc[~df['is_offtarget_safe'], 'Tier'] = "Tier 4 (High Off-Target Risk)"

    # Downgrade polyT guides (Functional Failure) - Top Priority Exclusion
    if 'polyT' in df.columns:
        df.loc[df['polyT'] == True, 'Tier'] = "Tier 5 (PolyT)"

    # Clean output columns
    # Save all columns
    df.to_csv(output_file, sep='\t', index=False)
    print(f"[SUCCESS] Prioritization complete.")
    print(f"  - Elite (Tier 1): {len(df[df['Tier'].str.contains('Tier 1')]):,}")
    print(f"  - Good (Tier 2):  {len(df[df['Tier'].str.contains('Tier 2')]):,}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prioritize guides based on CRISPR-Act3.0 and RS3.")
    parser.add_argument("--input", required=True, help="Output from extractor_guides.tsv")
    parser.add_argument("--output", required=True, help="Filename for prioritized results")
    parser.add_argument("--rs3_tracr", default='Hsu2013', help="Tracr model for RS3 score (e.g., Hsu2013, Chen2013)")
    args = parser.parse_args()
    run_prioritization(args.input, args.output, args.rs3_tracr)