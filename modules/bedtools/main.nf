#!/usr/bin/env nextflow

process MERGE_PEAKS {
    container 'ghcr.io/bf528/bedtools_samtools:latest'
    publishDir "${params.outdir}/merged_peaks", mode: 'copy'
    label 'process_low'

    input:
    tuple val(group_id), path(peaks)

    output:
    tuple val(group_id), path("${group_id}_merged.bed"), emit: merged_peaks

    script:
    """
    cat ${peaks} | \
        cut -f1-3 | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - > ${group_id}_merged.bed
    """

    stub:
    """
    touch ${group_id}_merged.bed
    """
}

process COUNT_PEAKS {
    container 'ghcr.io/bf528/bedtools_samtools:latest'
    publishDir "${params.outdir}/counts", mode: 'copy'
    label 'process_low'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(peaks)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_counts.txt"), emit: counts

    script:
    """
    # Check if peaks file is empty or has no peaks
    if [ ! -s ${peaks} ] || [ \$(wc -l < ${peaks}) -eq 0 ]; then
        echo "WARNING: Empty peaks file for ${sample}, creating dummy counts"
        # Create a dummy count file with zeros
        echo -e "chr1\t1\t1000\t0" > ${sample}_counts.txt
    else
        bedtools coverage -a ${peaks} -b ${bam} -counts > ${sample}_counts.txt
        
        # Check if coverage output is empty
        if [ ! -s ${sample}_counts.txt ]; then
            echo "WARNING: bedtools coverage produced empty output for ${sample}"
            echo -e "chr1\t1\t1000\t0" > ${sample}_counts.txt
        fi
    fi
    """

    stub:
    """
    touch ${sample}_counts.txt
    """
}