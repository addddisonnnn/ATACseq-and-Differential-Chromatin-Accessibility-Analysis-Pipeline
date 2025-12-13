#!/usr/bin/env nextflow

process FILTER_ALIGNMENTS {
    container 'ghcr.io/bf528/samtools:latest'
    publishDir "${params.outdir}/filtered", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_filtered.bam"), path("${sample}_filtered.bam.bai"), emit: filtered_bam
    path("${sample}_filter_stats.txt"), emit: stats

    script:
    """
    # Remove MT, low quality, duplicates
    samtools view -@ ${task.cpus} -h -q 30 -F 1804 ${bam} | \
        grep -v 'chrM' | \
        samtools view -@ ${task.cpus} -b - | \
        samtools sort -@ ${task.cpus} -o ${sample}_temp.bam -
    
    # Remove blacklist regions if available
    if [ -f "${params.blacklist}" ]; then
        bedtools intersect -v -abam ${sample}_temp.bam -b ${params.blacklist} > ${sample}_filtered.bam
        rm ${sample}_temp.bam
    else
        mv ${sample}_temp.bam ${sample}_filtered.bam
    fi
    
    samtools index -@ ${task.cpus} ${sample}_filtered.bam
    
    # Generate stats
    samtools flagstat ${sample}_filtered.bam > ${sample}_filter_stats.txt
    """

    stub:
    """
    touch ${sample}_filtered.bam
    touch ${sample}_filtered.bam.bai
    touch ${sample}_filter_stats.txt
    """
}