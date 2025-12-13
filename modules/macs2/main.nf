#!/usr/bin/env nextflow

process MACS2_CALLPEAK {
    container 'ghcr.io/bf528/macs3:latest'
    publishDir "${params.outdir}/peaks", mode: 'copy', pattern: "*"
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_peaks.narrowPeak"), emit: peaks
    path("${sample}_peaks.xls"), emit: stats
    path("${sample}_macs2.log"), emit: log

    script:
    """
    # Check number of reads first
    total_reads=\$(samtools view -c ${bam})
    echo "Total reads in BAM: \$total_reads" > ${sample}_macs2.log
    
    # MACS3 for ATAC-seq - REMOVED -B flag to avoid numpy bug
    macs3 callpeak \
        -t ${bam} \
        -f BAM \
        -g ${params.macs2_gsize} \
        -n ${sample} \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --keep-dup all \
        --call-summits \
        -q 0.05 \
        2>&1 | tee -a ${sample}_macs2.log
    
    # Check if peaks were called
    if [ ! -f ${sample}_peaks.narrowPeak ]; then
        echo "ERROR: MACS3 did not produce narrowPeak file" >> ${sample}_macs2.log
        touch ${sample}_peaks.narrowPeak
        echo "# No peaks called - MACS3 failed" > ${sample}_peaks.xls
    elif [ ! -s ${sample}_peaks.narrowPeak ]; then
        echo "WARNING: MACS3 produced empty narrowPeak file" >> ${sample}_macs2.log
    else
        num_peaks=\$(wc -l < ${sample}_peaks.narrowPeak)
        echo "SUCCESS: Called \$num_peaks peaks" >> ${sample}_macs2.log
    fi
    
    # Ensure xls file exists
    if [ ! -f ${sample}_peaks.xls ]; then
        echo "# No peaks called" > ${sample}_peaks.xls
    fi
    """

    stub:
    """
    touch ${sample}_peaks.narrowPeak
    touch ${sample}_peaks.xls
    touch ${sample}_macs2.log
    """
}