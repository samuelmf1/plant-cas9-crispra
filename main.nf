nextflow.enable.dsl=2

params.test_subset = false
params.annotation = null
params.output_prefix = 'result'
params.n_guides_pergene = 3

process PROMOTER_GUIDES {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path annotation
    path genome
    
    output:
    path "*_guides.tsv", emit: guides
    path "*.fa", emit: fasta
    
    script:
    """
    python3 ${projectDir}/src/extract_promoter_seqs.py \\
        --input ${annotation} \\
        --genome ${genome} \\
        --output_prefix ${params.output_prefix} \\
        --far-bp ${params.far_bp} \\
        --close-bp ${params.close_bp} \\
        --guides
    """
}

process SPLIT {
    input:
    path guides
    
    output:
    path "chunk_*.tsv", emit: chunks
    
    script:
    """
    # Create subset if requested
    if [ "${params.test_subset}" = "true" ]; then
        head -n 1001 ${guides} > subset.tsv
        INPUT_FILE=subset.tsv
    else
        INPUT_FILE=${guides}
    fi

    # Save header
    head -n 1 \$INPUT_FILE > header.txt
    
    # Calculate lines per chunk
    TOTAL_LINES=\$(tail -n +2 \$INPUT_FILE | wc -l)
    if [ \$TOTAL_LINES -gt 0 ]; then
        LINES_PER_CHUNK=\$(( (\$TOTAL_LINES + ${params.chunks} - 1) / ${params.chunks} ))
        
        # Split body
        tail -n +2 \$INPUT_FILE | split -l \$LINES_PER_CHUNK -a 3 -d - chunk_body_
        
        # Add header to each chunk
        for f in chunk_body_*; do
            cat header.txt \$f > chunk_\${f#chunk_body_}.tsv
        done
    else
        # Handle empty case
        cp \$INPUT_FILE chunk_000.tsv
    fi
    """
}

process OFFTARGETS {
    publishDir "${params.outdir}/work", mode: 'copy', pattern: '*_sassy.out'
    
    input:
    path guides
    path genome
    
    output:
    path "*_offtargets.tsv", emit: offtargets
    path "*_sassy.out", emit: sassy_logs
    
    script:
    """
    python3 ${projectDir}/src/calc_offtargets_sassy.py \\
        -i ${guides} \\
        -g ${genome} \\
        -o ${guides.baseName}_offtargets.tsv \\
        --family_id \\
        --sassy_log_file ${guides.baseName}_sassy.out \\
        --mismatch_intolerance ${params.mismatch_intolerance} \\
        --family_flank ${params.family_flank}
    """
}

process RANK_GUIDES {
    // Run in parallel on chunks
    publishDir "${params.outdir}/work", mode: 'copy', pattern: '*_ranked.tsv'
    
    input:
    path offtargets
    
    output:
    path "*_ranked.tsv", emit: prioritized
    
    script:
    """
    python3 ${projectDir}/src/guide_prioritization.py \\
        --input ${offtargets} \\
        --output ${offtargets.baseName}_ranked.tsv \\
        --rs3_tracr ${params.rs3_tracr}
    """
}

process MERGE {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path ranked_files
    
    output:
    path "${params.output_prefix}_prioritized_guides.tsv", emit: merged
    
    script:
    """
    # Get first file for header
    # Ensure consistent ordering if needed, but cat *.tsv usually generic
    # ranked_files is a list in Nextflow. 
    # We can rely on ls or just loop. Use ls *.tsv
    
    first_file=\$(ls *_ranked.tsv | head -n 1)
    head -n 1 \$first_file > ${params.output_prefix}_prioritized_guides.tsv
    
    # Concatenate all bodies
    # Use -q to suppress headers if using tail? No, tail -q suppresses file name headers.
    
    for f in *_ranked.tsv; do
        tail -n +2 \$f >> ${params.output_prefix}_prioritized_guides.tsv
    done
    """
}

process SELECT_GUIDES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path prioritized_guides
    path annotation_file
    
    output:
    path "*_selected_guides.tsv", emit: selected
    path "*_annotation_with_stats.csv", emit: annotated_stats
    path "plots/*.svg", emit: plots
    
    script:
    """
    python3 ${projectDir}/src/select_final_guides.py \\
        --input ${prioritized_guides} \\
        --annotation ${annotation_file} \\
        --output_prefix ${params.output_prefix} \\
        --mismatch_intolerance ${params.mismatch_intolerance} \\
        --n_guides_pergene ${params.n_guides_pergene}
    """
}

workflow {
    if (!params.genome) {
        error "Please specify --genome path/to/genome.fa.gz"
    }
    if (!params.annotation) {
        error "Please specify --annotation path/to/annotation.csv"
    }
    
    genome_ch = Channel.fromPath(params.genome)
    annotation_ch = Channel.fromPath(params.annotation)
    
    PROMOTER_GUIDES(annotation_ch, genome_ch)
    
    SPLIT(PROMOTER_GUIDES.out.guides)
    
    chunks_ch = SPLIT.out.chunks.flatten()
    
    OFFTARGETS(chunks_ch, genome_ch.first())
    
    // Pipe calculated chunks directly to RANK_GUIDES
    RANK_GUIDES(OFFTARGETS.out.offtargets)
    
    // Merge the ranked chunks
    MERGE(RANK_GUIDES.out.prioritized.collect())
    
    SELECT_GUIDES(MERGE.out.merged, annotation_ch.first())
}
