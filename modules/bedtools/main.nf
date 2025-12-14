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
    # Set temp directory to current work directory
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR
    
    echo "Processing group: ${group_id}"
    echo "Number of peak files: \$(echo ${peaks} | wc -w)"
    
    # Merge all peak files
    cat ${peaks} | \
        awk 'BEGIN{OFS="\t"} NF >= 3 {print \$1, \$2, \$3}' | \
        sort -T \$TMPDIR -k1,1 -k2,2n | \
        bedtools merge -i - > ${group_id}_merged.bed
    
    # Report result
    num_merged=\$(wc -l < ${group_id}_merged.bed)
    echo "Merged to \$num_merged regions for ${group_id}"
    
    # Cleanup
    rm -rf \$TMPDIR
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