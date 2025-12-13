#!/usr/bin/env nextflow

process ANNOTATE_PEAKS {
    container 'ghcr.io/bf528/homer:latest'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    label 'process_low'

    input:
    tuple val(cell_type), path(peaks)

    output:
    path("${cell_type}_annotated.txt"), emit: annotations

    script:
    """
    annotatePeaks.pl ${peaks} mm10 > ${cell_type}_annotated.txt
    """

    stub:
    """
    touch ${cell_type}_annotated.txt
    """
}

process MOTIF_ANALYSIS {
    container 'ghcr.io/bf528/homer:latest'
    publishDir "${params.outdir}/motifs", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(cell_type), path(peaks)

    output:
    path("${cell_type}_motifs/"), emit: motif_results

    script:
    """
    findMotifsGenome.pl ${peaks} mm10 ${cell_type}_motifs/ \
        -size 200 -mask
    """

    stub:
    """
    mkdir -p ${cell_type}_motifs
    touch ${cell_type}_motifs/homerResults.html
    """
}