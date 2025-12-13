#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir "${params.outdir}/aligned", mode: 'copy'
    label 'process_high'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(reads)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    path("${sample}_alignment_stats.txt"), emit: log

    script:
    """
    bowtie2 -p ${task.cpus} \
        --very-sensitive \
        -x ${params.genome_index} \
        -U ${reads} \
        2> ${sample}_alignment_stats.txt | \
    samtools view -@ ${task.cpus} -bS - | \
    samtools sort -@ ${task.cpus} -o ${sample}.bam -
    
    samtools index -@ ${task.cpus} ${sample}.bam
    """

    stub:
    """
    touch ${sample}.bam
    touch ${sample}.bam.bai
    touch ${sample}_alignment_stats.txt
    """
}