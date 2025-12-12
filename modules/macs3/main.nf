#!/usr/bin/env nextflow

process MACS3_CALLPEAK {
    container 'ghcr.io/bf528/macs3:latest'
    publishDir "${params.outdir}/peaks", mode: 'copy'
    label 'process_medium'
    
    input:
    tuple val(condition), path(bams)
    
    output:
    tuple val(condition), path("${condition}_peaks.narrowPeak"), emit: narrowPeak
    path("${condition}_*"), emit: all_files
    
    script:
    def bam_files = bams.findAll { it.toString().endsWith('.bam') }
    def bam_list = bam_files.collect{it.toString()}.join(' ')
    """
    macs3 callpeak \\
        -t ${bam_list} \\
        -f BAM \\
        -g ${params.macs_gsize} \\
        -n ${condition} \\
        -q ${params.macs_qvalue} \\
        --nomodel \\
        --shift -75 \\
        --extsize 150 \\
        --keep-dup all
    """
    
    stub:
    """
    touch ${condition}_peaks.narrowPeak
    touch ${condition}_peaks.xls
    touch ${condition}_summits.bed
    """
}