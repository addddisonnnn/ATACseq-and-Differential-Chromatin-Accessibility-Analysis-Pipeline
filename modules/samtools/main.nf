#!/usr/bin/env nextflow

process FILTER_ALIGNMENTS {
    container 'ghcr.io/bf528/bedtools_samtools:latest'
    publishDir "${params.outdir}/filtered", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_filtered.bam"), path("${sample}_filtered.bam.bai"), emit: filtered_bam
    path("${sample}_filter_stats.txt"), emit: stats

    script:
    """
    # Remove MT, low quality, and duplicates
    # F1804 = unmapped (4) + not primary (256) + PCR duplicate (1024) + supplementary (2048)
    samtools view -@ ${task.cpus} -h -q 20 -F 1804 ${bam} | \
        awk '\$3 != "chrM" && \$3 != "MT"' | \
        samtools view -@ ${task.cpus} -b - | \
        samtools sort -@ ${task.cpus} -o ${sample}_filtered.bam -
    
    samtools index -@ ${task.cpus} ${sample}_filtered.bam
    
    # Generate stats
    samtools flagstat ${sample}_filtered.bam > ${sample}_filter_stats.txt
    echo "" >> ${sample}_filter_stats.txt
    echo "Reads per chromosome (top 10):" >> ${sample}_filter_stats.txt
    samtools idxstats ${sample}_filtered.bam | sort -k3 -nr | head -10 >> ${sample}_filter_stats.txt
    """

    stub:
    """
    touch ${sample}_filtered.bam
    touch ${sample}_filtered.bam.bai
    touch ${sample}_filter_stats.txt
    """
}