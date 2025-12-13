#!/usr/bin/env nextflow

process MACS2_CALLPEAK {
    container 'ghcr.io/bf528/macs3:latest'
    publishDir "${params.outdir}/peaks", mode: 'copy'
    label 'process_high'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_peaks.narrowPeak"), emit: peaks
    path("${sample}_peaks.xls"), emit: stats

    script:
    """
    macs3 callpeak \
        -t ${bam} \
        -f BAM \
        -g ${params.macs2_gsize} \
        -n ${sample} \
        --nomodel \
        --shift -100 \
        --extsize 200 \
        --keep-dup all \
        -q ${params.qvalue_cutoff}
    """

    stub:
    """
    touch ${sample}_peaks.narrowPeak
    touch ${sample}_peaks.xls
    """
}