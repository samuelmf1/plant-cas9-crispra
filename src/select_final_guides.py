#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

def select_n_guides(gene_df, n, min_dist):
    """
    Selects up to n guides from sorted gene_df such that all selected guides
    are at least min_dist apart.
    """
    selected = []
    # gene_df is assumed sorted by priority
    for idx, row in gene_df.iterrows():
        if len(selected) >= n:
            break
        
        pos = row['coord']
        # Check distance constraint
        conflict = False
        for s in selected:
            # Distance constraint only strictly applies if same orientation (FWD/FWD or RVS/RVS)
            # If orientations differ (FWD/RVS), overlap is allowed.
            if s.get('Orientation') == row.get('Orientation'):
                if abs(s['coord'] - pos) <= min_dist:
                    conflict = True
                    break
        
        if not conflict:
            selected.append(row)
            
    return selected

def main():
    parser = argparse.ArgumentParser(description="Select best n guides per gene with relaxation.")
    parser.add_argument("--input", required=True, help="Prioritized guides TSV")
    parser.add_argument("--annotation", required=True, help="Original annotation CSV")
    parser.add_argument("--output_prefix", required=True, help="Output prefix")
    parser.add_argument("--mismatch_intolerance", type=int, default=2, help="Number of mismatch_intolerance allowed")
    parser.add_argument("--n_guides_pergene", type=int, default=3, help="Number of guides to select per gene")
    args = parser.parse_args()
    
    print(f"[INFO] Reading guides from {args.input}...")
    try:
        guides_df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read input guides: {e}")

    print(f"[INFO] Reading annotation from {args.annotation}...")
    try:
        annot_df = pd.read_csv(args.annotation, sep=None, engine='python')
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read annotation file: {e}")
        
    # Validation
    if 'gene' not in guides_df.columns:
         sys.exit("[ERROR] Input guides TSV missing 'gene' column.")
    
    # Identify gene column in annotation (case insensitive)
    # The snippet shows 'Gene', 'Chr', etc.
    annot_cols_lower = {c.lower(): c for c in annot_df.columns}
    if 'gene' not in annot_cols_lower:
        sys.exit(f"[ERROR] Annotation file missing 'Gene' column. Found: {annot_df.columns}")
    
    gene_col_name = annot_cols_lower['gene']
    
    # Filter unsafe guides
    if 'is_offtarget_safe' not in guides_df.columns:
        print("[WARN] 'is_offtarget_safe' column missing. Assuming all safe.")
        guides_df['is_offtarget_safe'] = True
        
    safe_mask = guides_df['is_offtarget_safe'] == True
    valid_guides = guides_df[safe_mask].copy()
    
    # Sort
    # Sort
    if 'Tier' in valid_guides.columns:
        # Exclude Tier 5 (PolyT) guides
        valid_guides = valid_guides[valid_guides['Tier'] != "Tier 5 (PolyT)"]
        valid_guides.sort_values(by=['Tier', 'plant_act_score', 'rs3_score'], ascending=[True, False, False], inplace=True)
    
    results = []
    gene_stats = {}
    
    grouped = valid_guides.groupby('gene')
    
    # Process ALL genes from annotation
    all_annot_genes = annot_df[gene_col_name].unique()
    print(f"[INFO] Processing {len(all_annot_genes)} genes from annotation...")
    
    for gene in all_annot_genes:
        # Default stats
        stats = {
            'gene': gene,
            'selection_status': 'No guides found',
            'n_guides_selected': 0,
            'relaxation_level': -1 # Indicator for none
        }
        
        if gene in grouped.groups:
            g_df = grouped.get_group(gene)
            
            # Selection logic
            best_selection = []
            status = "Success (Strict)"
            dist_used = 22
            
            # Try 22bp
            selection = select_n_guides(g_df, args.n_guides_pergene, 22)
            
            if len(selection) == args.n_guides_pergene:
                best_selection = selection
            else:
                # Relax
                found_full = False
                for dist in range(21, 4, -1):
                    selection = select_n_guides(g_df, args.n_guides_pergene, dist)
                    if len(selection) == args.n_guides_pergene:
                        best_selection = selection
                        status = f"Success (Relaxed >{dist:02d}bp)"
                        dist_used = dist
                        found_full = True
                        break
                
                if not found_full:
                    # Fallback to minimal distance (4bp)
                    selection = select_n_guides(g_df, args.n_guides_pergene, 4)
                    best_selection = selection
                    dist_used = 4
                    
                    if len(selection) == args.n_guides_pergene:
                        status = f"Success (Relaxed >04bp)"
                    else:
                        status = f"Partial ({len(selection)}/{args.n_guides_pergene})"
            
            # Save selected guides
            infamily_hits = 0
            for row in best_selection:
                r_dict = row.to_dict()
                r_dict['selection_status'] = status
                r_dict['min_dist_used'] = dist_used
                results.append(r_dict)
                
                # Check in-family off-targets (0, 1, 2 mismatch_intolerance)
                # Columns: offtarget_0_in_family, offtarget_1_in_family, etc.
                self_col = 'offtarget_0_in_family'
                for k in range(args.mismatch_intolerance + 1):
                    col = f'offtarget_{k}_in_family'
                    if col in row:
                        infamily_hits += row[col] if col != self_col else row[col] - 1

            stats['selection_status'] = status
            stats['n_guides_selected'] = len(best_selection)
            stats['relaxation_level'] = dist_used
            stats['has_infamily_offtargets'] = (infamily_hits > 0)
            stats['total_infamily_offtargets'] = infamily_hits
        
        gene_stats[gene] = stats

    # 1. Output Selected Guides
    if results:
        out_guides_df = pd.DataFrame(results)
        guides_filename = f"{args.output_prefix}_selected_guides.tsv"
        out_guides_df.to_csv(guides_filename, sep='\t', index=False)
        print(f"[SUCCESS] Saved selected guides to {guides_filename}")
    else:
        print("[WARN] No guides selected for any gene.")

    # 2. Augment Annotation File
    stats_df_list = list(gene_stats.values())
    stats_df = pd.DataFrame(stats_df_list)
    
    # Merge stats back to annotation
    stats_df.rename(columns={'gene': gene_col_name}, inplace=True)
    
    merged_annot = pd.merge(annot_df, stats_df, on=gene_col_name, how='left')
    
    annot_out_filename = f"{args.output_prefix}_annotation_with_stats.csv"
    merged_annot.to_csv(annot_out_filename, index=False)
    print(f"[SUCCESS] Saved annotated stats to {annot_out_filename}")
    
    # 3. EDA & Plots
    print("[INFO] Generating EDA plots...")
    
    plots_dir = "plots"
    os.makedirs(plots_dir, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    
    # 3. EDA & Plots
    print("[INFO] Generating EDA plots...")
    
    plots_dir = "plots"
    # Clean up old plots if needed (optional, but good practice if starting fresh)
    # shutil.rmtree(plots_dir, ignore_errors=True)
    os.makedirs(plots_dir, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)

    # --- 1. Filtering Logic ---
    print("[INFO] Applying boolean filters for EDA...")
    
    # Identify filter columns
    filter_cols = [c for c in merged_annot.columns if c.startswith("filt_") or c.startswith("filter_")]
    
    # Calculate initial stats
    total_genes = len(merged_annot)
    
    if filter_cols:
        # Determine genes to remove (True in any filter col)
        # Assuming True means "Filter Out" based on "filt_shortnoexpAA_pseudogenes" logic
        genes_to_remove_mask = merged_annot[filter_cols].any(axis=1)
        n_removed = genes_to_remove_mask.sum()
        n_remaining = total_genes - n_removed
        
        # Plot 1: Gene Removal Stats
        plt.figure(figsize=(8, 6))
        stats_df = pd.DataFrame({
            'Status': ['Removed', 'Remaining'],
            'Count': [n_removed, n_remaining]
        })
        stats_df['Percentage'] = (stats_df['Count'] / total_genes) * 100
        
        ax = sns.barplot(data=stats_df, x='Status', y='Count', hue='Status', palette=['red', 'green'], legend=False)
        plt.title(f"Gene Filtering Impact (Total: {total_genes})")
        plt.ylabel("Gene Count")
        
        # Add labels
        for i, row in stats_df.iterrows():
            ax.text(i, row['Count'], f"{row['Count']}\n({row['Percentage']:.1f}%)", 
                    ha='center', va='bottom')
            
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "eda_1_gene_filtering.svg"), format='svg', dpi=300)
        plt.close()
        
        # Apply filter for subsequent plots
        eda_df = merged_annot[~genes_to_remove_mask].copy()
        print(f"[INFO] Removed {n_removed} genes. EDA will proceed with {n_remaining} genes.")
    else:
        print("[INFO] No filter columns found. Proceeding with all genes.")
        eda_df = merged_annot.copy()

    if len(eda_df) == 0:
        print("[WARN] No genes remaining after filtering. Skipping plots.")
        return

    # --- 2. Chromosomal Position Plot ---
    # X=Chr, Y=Position (Start), Color=N_guides
    if 'Chr' in eda_df.columns and 'Start' in eda_df.columns:
        plt.figure(figsize=(12, 6))
        # Ensure Chr is sorted nicely if possible
        eda_df['Chr'] = eda_df['Chr'].astype(str)
        chr_order = sorted(eda_df['Chr'].unique())
        
        sns.scatterplot(data=eda_df, x='Chr', y='Start', hue='n_guides_selected', 
                        palette='viridis', alpha=0.6, s=20)
        
        plt.title("Gene Distribution by Chromosome and Position")
        plt.legend(title='Guides Selected', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "eda_2_chromosomal_dist.svg"), format='svg', dpi=300)
        plt.close()

    # --- 3. Guide Count Bar Chart ---
    plt.figure(figsize=(10, 6))
    ax = sns.countplot(data=eda_df, x='n_guides_selected', hue='n_guides_selected', palette='viridis', legend=False)
    plt.title("Detailed Guide Count per Gene")
    plt.xlabel("Number of Guides Selected")
    plt.ylabel("Gene Count")
    
    # Add count labels at top of bars
    for container in ax.containers:
        ax.bar_label(container)
        
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "eda_3_guide_counts.svg"), format='svg', dpi=300)
    plt.close()

    # --- 4. Failure Analysis (Low Guide Count) ---
    # For number of guides < n_guides_pergene, plot counts of flags
    failed_df = eda_df[eda_df['n_guides_selected'] < args.n_guides_pergene]
    
    if len(failed_df) > 0:
        flags = ['is_hdr_gene', 'is_ddr_gene', 'is_repair_gene', 'is_pseudogene', 'is_hypothetical_protein']
        valid_flags = [f for f in flags if f in failed_df.columns]
        
        if valid_flags:
            flag_counts = failed_df[valid_flags].sum().reset_index()
            flag_counts.columns = ['Feature', 'Count']
            
            plt.figure(figsize=(10, 6))
            sns.barplot(data=flag_counts, x='Feature', y='Count', hue='Feature', palette='Reds_d', legend=False)
            plt.title(f"Features of Failed Genes (< {args.n_guides_pergene} guides)")
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "eda_4_failure_analysis.svg"), format='svg', dpi=300)
            plt.close()

    # --- 5. Comparisons (All vs Selected) ---
    if results:
        sel_guides_df = pd.DataFrame(results)
        
        # A. RS3 Score Comparison
        if 'rs3_score' in valid_guides.columns and 'rs3_score' in sel_guides_df.columns:
            try:
                plt.figure(figsize=(10, 6))
                
                # Filter valid_guides to only include those from the remaining genes in eda_df
                # This ensures we are comparing the relevant pool
                gene_ids_remaining = set(eda_df['Gene'])
                pool_filtered = valid_guides[valid_guides['gene'].isin(gene_ids_remaining)]
                
                pool_scores = pool_filtered[['rs3_score']].copy()
                pool_scores['Group'] = 'All Candidates'
                
                sel_scores = sel_guides_df[['rs3_score']].copy()
                sel_scores['Group'] = 'Selected'
                
                combined = pd.concat([pool_scores, sel_scores], ignore_index=True)
                
                sns.histplot(data=combined, x='rs3_score', hue='Group', kde=True, 
                             stat="density", common_norm=False, palette='tab10')
                plt.title("RS3 Score Distribution: Candidates vs Selected")
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, "eda_5_rs3_comparison.svg"), format='svg', dpi=300)
                plt.close()
            except Exception as e:
                print(f"[WARN] Failed to plot RS3 comparison: {e}")

        # B. Off-target Heatmap: out_family offtarget counts
        # We want to compare the distribution of off-target counts between All and Selected.
        # Heatmap: Rows = Off-target Count, Cols = Group, Color = Frequency (normalized col-wise)
        
        # Calculate 'num_out_family' for pool if missing (summing _out_family cols)
        if 'num_out_family' not in valid_guides.columns:
             out_fam_cols = [c for c in valid_guides.columns if c.endswith('_out_family')]
             if out_fam_cols:
                 valid_guides['num_out_family'] = valid_guides[out_fam_cols].fillna(0).sum(axis=1)
             else:
                 valid_guides['num_out_family'] = 0 # Default if no info
                 
        if 'num_out_family' in valid_guides.columns and 'num_out_family' in sel_guides_df.columns:
            try:
                # Prepare data
                pool_counts = valid_guides['num_out_family'].reset_index(drop=True)
                sel_counts = sel_guides_df['num_out_family'].reset_index(drop=True)
                
                # Create a DataFrame for heatmap
                # Truncate high counts for readability if needed (e.g. 5+)
                MAX_OFF = 5
                
                def bucket_off(s):
                    return s.apply(lambda x: f"{int(x)}" if x < MAX_OFF else f"{MAX_OFF}+")
                
                pool_dist = bucket_off(pool_counts).value_counts(normalize=True).rename("All Candidates")
                sel_dist = bucket_off(sel_counts).value_counts(normalize=True).rename("Selected")
                
                # Merge and fill missing categories with 0
                heatmap_df = pd.concat([pool_dist, sel_dist], axis=1).fillna(0)
                
                # Sort index naturally (0, 1, 2, ..., 5+)
                # Helper for sorting strings "0", "1", ... "5+"
                def sort_key(idx):
                    return [int(x.replace('+', '')) for x in idx]
                    
                heatmap_df = heatmap_df.sort_index(key=sort_key)
                
                plt.figure(figsize=(6, 6))
                sns.heatmap(heatmap_df, annot=True, fmt=".1%", cmap="Blues", cbar_kws={'label': 'Percentage of Group'})
                plt.title("Off-Target Count Distribution")
                plt.ylabel("Number of Out-Family Off-Targets")
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, "eda_5b_offtarget_heatmap.svg"), format='svg', dpi=300)
                plt.close()
            except Exception as e:
                print(f"[WARN] Failed to plot Off-Target Heatmap: {e}")

    # --- 6. Orientation Heatmap ---
    # Orientation (FWD/RVS) vs # Guides per Gene
    # We want a heatmap where:
    # Rows = N guides (0, 1, 2...)
    # Cols = Orientation Mix? Or just raw counts?
    # "Orientation (FWD/RVS) vs number of guides per gene in a heatmap"
    # Maybe:
    # X axis: # Guides Selected (0, 1, 2, 3)
    # Y axis: Orientation Breakdown (e.g. 3F, 2F1R, 1F2R, 3R)?
    # Or simply: Count of FWD guides vs Count of RVS guides?
    
    # Interpretation: Heatmap where X = Count of FWD guides, Y = Count of RVS guides, Color = Count of Genes.
    # This neatly summarizes the mix.
    if results:
         try:
            # We need to aggregate by gene to get n_FWD and n_RVS
            sel_guides_df = pd.DataFrame(results)
            
            # Pivot to get counts per gene
            # gene | FWD | RVS
            ori_counts = sel_guides_df.pivot_table(index='gene', columns='Orientation', values='coord', aggfunc='count', fill_value=0)
            
            # Ensure columns exist
            if 'FWD' not in ori_counts.columns: ori_counts['FWD'] = 0
            if 'RVS' not in ori_counts.columns: ori_counts['RVS'] = 0
            
            # Merge with eda_df to include genes with 0 guides
            # eda_df has 'Gene'
            ori_merged = eda_df[['Gene']].merge(ori_counts, left_on='Gene', right_index=True, how='left').fillna(0)
            
            # Group by (FWD, RVS) to get gene counts
            heatmap_data = ori_merged.groupby(['FWD', 'RVS']).size().unstack(fill_value=0)
            
            plt.figure(figsize=(8, 6))
            sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='YlGnBu')
            plt.title("Gene Count by Guide Orientation Mix")
            plt.xlabel("Number of RVS Guides")
            plt.ylabel("Number of FWD Guides")
            plt.gca().invert_yaxis() # 0 at bottom
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "eda_6_orientation_heatmap.svg"), format='svg', dpi=300)
            plt.close()
            
         except Exception as e:
                print(f"[WARN] Failed to plot Orientation Heatmap: {e}")

    print(f"[SUCCESS] Plots saved to {plots_dir}/")

if __name__ == "__main__":
    main()
