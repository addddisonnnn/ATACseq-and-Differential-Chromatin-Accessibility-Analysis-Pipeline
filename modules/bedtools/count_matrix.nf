#!/usr/bin/env nextflow

process CREATE_COUNT_MATRIX {
    container 'ghcr.io/bf528/bedtools_samtools:latest'
    publishDir "${params.outdir}/counts", mode: 'copy'
    label 'process_high'
    
    input:
    path(peaks)
    path(bams)
    
    output:
    path("consensus_peaks.bed"), emit: consensus
    path("count_matrix.txt"), emit: counts
    path("sample_info.txt"), emit: sample_info
    
    script:
    """
    // Merge all peaks to create consensus peakset
    cat *.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > consensus_peaks.bed
    
    // Count reads in peaks for each BAM file
    for bam in *.bam; do
        sample=\${bam%_noMT.bam}
        bedtools coverage -a consensus_peaks.bed -b \$bam -counts | cut -f4 > \${sample}_counts.txt
    done
    
    // Create sample info file
    for bam in *.bam; do
        sample=\${bam%_noMT.bam}
        condition=\$(echo \$sample | grep -o 'cDC[12]_[WK][TO]')
        echo -e "\${sample}\\t\${condition}" >> sample_info.txt
    done
    
    // Combine into count matrix
    paste consensus_peaks.bed *_counts.txt > count_matrix.txt
    """
    
    stub:
    """
    touch consensus_peaks.bed
    touch count_matrix.txt
    touch sample_info.txt
    """
}