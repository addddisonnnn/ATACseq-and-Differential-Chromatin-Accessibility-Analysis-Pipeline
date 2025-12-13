#!/usr/bin/env nextflow

process DIFF_ANALYSIS {
    container 'ghcr.io/bf528/pandas:latest'
    publishDir "${params.outdir}/differential", mode: 'copy'
    label 'process_low'

    input:
    tuple val(cell_type), val(samples), val(conditions), val(replicates), path(counts)

    output:
    tuple val(cell_type), path("${cell_type}_diff_peaks.bed"), emit: diff_peaks
    path("${cell_type}_diff_results.txt"), emit: results

    script:
    """
    diff_analysis.py \
        --counts ${counts} \
        --samples ${samples.join(',')} \
        --conditions ${conditions.join(',')} \
        --output ${cell_type}_diff_results.txt \
        --bed ${cell_type}_diff_peaks.bed
    """

    stub:
    """
    touch ${cell_type}_diff_peaks.bed
    touch ${cell_type}_diff_results.txt
    """
}