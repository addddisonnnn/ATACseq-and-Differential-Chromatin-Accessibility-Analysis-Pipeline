#!/usr/bin/env nextflow

process DOWNLOAD_SRA {
    container 'ghcr.io/bf528/sratools:latest'
    publishDir "${params.outdir}/fastq", mode: 'copy'
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(condition), val(replicate), val(srr)
    
    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}.fastq.gz"), emit: fastq
    
    script:
    """
    fasterq-dump ${srr} -O . -e ${task.cpus}
    mv ${srr}.fastq ${sample_id}.fastq
    gzip ${sample_id}.fastq
    """
    
    stub:
    """
    touch ${sample_id}.fastq.gz
    """
}