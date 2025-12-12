#!/usr/bin/env nextflow

process TRIMMOMATIC {
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(condition), val(replicate), path(reads)
    
    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads
    
    script:
    """
    trimmomatic SE -threads ${task.cpus} \\
        ${reads} \\
        ${sample_id}_trimmed.fastq.gz \\
        ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
    
    stub:
    """
    touch ${sample_id}_trimmed.fastq.gz
    """
}