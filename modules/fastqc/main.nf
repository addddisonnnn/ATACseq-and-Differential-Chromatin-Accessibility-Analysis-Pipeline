#!/usr/bin/env nextflow

process FASTQC {
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    label 'process_low'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(reads)

    output:
    path("*_fastqc.{zip,html}"), emit: fastqc_results
    path("*_fastqc.zip"), emit: fastqc_zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """

    stub:
    """
    touch ${sample}_fastqc.html
    touch ${sample}_fastqc.zip
    """
}