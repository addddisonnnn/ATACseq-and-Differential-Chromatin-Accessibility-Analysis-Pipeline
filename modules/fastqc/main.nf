#!/usr/bin/env nextflow

process FASTQC {
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    label 'process_low'
    
    input:
    tuple val(sample_id), val(condition), val(replicate), path(reads)
    
    output:
    path("*_fastqc.zip"), emit: fastqc_results, optional: true
    
    script:
    """
    mkdir -p tmp_fastqc
    
    fastqc -t ${task.cpus} \
        -d tmp_fastqc \
        --noextract \
        ${reads} || echo "FastQC completed with warnings"
    
    rm -rf tmp_fastqc
    """
}