#!/usr/bin/env nextflow

process TRIMMOMATIC {
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(reads)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample}_trimming_report.txt"), emit: log

    script:
    """
    trimmomatic SE -threads ${task.cpus} \
        ${reads} \
        ${sample}_trimmed.fastq.gz \
        ILLUMINACLIP:\$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> ${sample}_trimming_report.txt
    """

    stub:
    """
    touch ${sample}_trimmed.fastq.gz
    touch ${sample}_trimming_report.txt
    """
}