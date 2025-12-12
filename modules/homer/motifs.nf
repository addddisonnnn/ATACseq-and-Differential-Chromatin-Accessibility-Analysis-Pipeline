#!/usr/bin/env nextflow

process FIND_MOTIFS {
    container 'ghcr.io/bf528/homer:latest'
    publishDir "${params.outdir}/motifs", mode: 'copy'
    label 'process_high'
    
    input:
    path(peaks)
    path(genome)
    
    output:
    path("motif_*"), emit: motifs
    
    script:
    """
    for peak_file in ${peaks}; do
        base=\${peak_file%.bed}
        output_dir="motif_\${base}"
        findMotifsGenome.pl \$peak_file \\
            ${genome} \\
            \$output_dir \\
            -size 200 \\
            -mask \\
            -p ${task.cpus}
    done
    """
    
    stub:
    """
    for peak_file in ${peaks}; do
        mkdir -p motif_\${peak_file%.bed}
        touch motif_\${peak_file%.bed}/homerResults.html
    done
    """
}