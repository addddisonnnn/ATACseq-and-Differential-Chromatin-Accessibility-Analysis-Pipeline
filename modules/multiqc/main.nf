#!/usr/bin/env nextflow

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    label 'process_low'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data/"), emit: data

    script:
    """
    multiqc .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}