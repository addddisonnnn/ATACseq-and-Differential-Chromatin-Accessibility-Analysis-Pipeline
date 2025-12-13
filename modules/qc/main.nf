#!/usr/bin/env nextflow

process FRIP_SCORE {
    container 'ghcr.io/bf528/bedtools_samtools:latest'
    publishDir "${params.outdir}/qc", mode: 'copy'
    label 'process_low'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai), path(peaks)

    output:
    path("${sample}_FRiP.txt"), emit: frip_stats

    script:
    """
    # Total reads
    total=\$(samtools view -c ${bam})
    
    # Reads in peaks
    inpeaks=\$(bedtools intersect -a ${bam} -b ${peaks} -u | samtools view -c)
    
    # Calculate FRiP
    frip=\$(echo "scale=4; \$inpeaks / \$total" | bc)
    
    echo -e "Sample\tTotal_Reads\tReads_in_Peaks\tFRiP" > ${sample}_FRiP.txt
    echo -e "${sample}\t\$total\t\$inpeaks\t\$frip" >> ${sample}_FRiP.txt
    """

    stub:
    """
    echo -e "Sample\tTotal_Reads\tReads_in_Peaks\tFRiP\n${sample}\t1000000\t500000\t0.5" > ${sample}_FRiP.txt
    """
}