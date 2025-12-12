#!/usr/bin/env nextflow

process ANNOTATE_PEAKS {
    container 'ghcr.io/bf528/homer:latest'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    label 'process_medium'
    
    input:
    path(peaks)
    path(gtf)
    
    output:
    path("*_annotated.txt"), emit: annotated
    
    script:
    """
    for peak_file in ${peaks}; do
        base=\${peak_file%.bed}
        annotatePeaks.pl \$peak_file \\
            ${params.genome} \\
            -gtf ${gtf} \\
            > \${base}_annotated.txt
    done
    """
    
    stub:
    """
    for peak_file in ${peaks}; do
        touch \${peak_file%.bed}_annotated.txt
    done
    """
}