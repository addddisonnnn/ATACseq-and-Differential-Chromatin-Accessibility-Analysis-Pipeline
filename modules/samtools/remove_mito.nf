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
    # Remove mitochondrial reads (chrM or MT depending on reference)
    samtools idxstats ${bam} | cut -f1 | grep -v -E 'chrM|MT' > keep_chroms.txt
    
    samtools view -@ ${task.cpus} -b -h ${bam} \$(cat keep_chroms.txt | tr '\\n' ' ') \\
        > ${sample_id}_noMT.bam
    
    samtools index ${sample_id}_noMT.bam
    """
    
    stub:
    """
    touch ${sample_id}_noMT.bam
    touch ${sample_id}_noMT.bam.bai
    """
}