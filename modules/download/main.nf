#!/usr/bin/env nextflow

process DOWNLOAD_SRA {
    container 'ghcr.io/bf528/sratools:latest'
    publishDir "${params.outdir}/fastq", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(srr), val(cell_type), val(condition), val(replicate)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}.fastq.gz"), emit: reads

    script:
    """
    fasterq-dump ${srr} -O . -e ${task.cpus}
    mv ${srr}.fastq ${sample}.fastq
    gzip ${sample}.fastq
    """

    stub:
    """
    touch ${sample}.fastq.gz
    """
}