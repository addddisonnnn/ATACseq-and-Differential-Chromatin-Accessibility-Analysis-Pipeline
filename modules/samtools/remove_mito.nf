#!/usr/bin/env nextflow

process REMOVE_MITO {
    container 'ghcr.io/bf528/samtools:latest'
    publishDir "${params.outdir}/no_mito", mode: 'copy'
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(condition), val(replicate), path(bam), path(bai)
    
    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}_noMT.bam"), path("${sample_id}_noMT.bam.bai"), emit: bam
    
    script:
    """
    # Debug: Show input files
    echo "=== DEBUG: REMOVE_MITO ==="
    echo "sample_id: ${sample_id}"
    echo "bam file: ${bam}"
    echo "bai file: ${bai}"
    
    # Remove mitochondrial reads (chrM or MT depending on reference)
    samtools idxstats "${bam}" | cut -f1 | grep -v -E 'chrM|MT' > keep_chroms.txt
    
    # Check what chromosomes we're keeping
    echo "Chromosomes to keep:"
    cat keep_chroms.txt
    
    # Create a safe list of chromosomes
    chrom_list=\$(cat keep_chroms.txt | tr '\\n' ' ')
    
    # Extract non-mitochondrial reads
    samtools view -@ ${task.cpus} -b -h "${bam}" \${chrom_list} > "${sample_id}_noMT.bam"
    
    # Index the output BAM
    samtools index "${sample_id}_noMT.bam"
    
    echo "Successfully created ${sample_id}_noMT.bam"
    """
    
    stub:
    """
    touch "${sample_id}_noMT.bam"
    touch "${sample_id}_noMT.bam.bai"
    """
}