#!/usr/bin/env nextflow

process FEATURE_COUNTS {
    container 'ghcr.io/bf528/featurecounts:latest'
    publishDir "${params.outdir}/counts_rna", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), 
          path("${sample}_counts.txt"), emit: counts
    path("${sample}_counts.txt.summary"), emit: summary

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${params.genome_gtf} \
        -o ${sample}_counts.txt \
        -t exon \
        -g gene_id \
        ${bam}
    """

    stub:
    """
    touch ${sample}_counts.txt
    touch ${sample}_counts.txt.summary
    """
}