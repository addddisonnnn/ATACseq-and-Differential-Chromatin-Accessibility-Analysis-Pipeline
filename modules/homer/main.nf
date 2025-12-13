#!/usr/bin/env nextflow

process HOMER_CALLPEAK {
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/peaks", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}_peaks.narrowPeak"), emit: peaks
    path("${sample}_peaks.txt"), emit: stats
    path("${sample}_homer.log"), emit: log

    script:
    """
    # Create tag directory from BAM
    makeTagDirectory ${sample}_tags ${bam} -format sam 2>&1 | tee ${sample}_homer.log
    
    # Find peaks using HOMER for ATAC-seq
    # -style factor: for sharp peaks (ATAC-seq, ChIP-seq TFs)
    # -size 150: expected fragment size
    # -minDist 150: minimum distance between peaks
    findPeaks ${sample}_tags -style factor -size 150 -minDist 150 -o auto 2>&1 | tee -a ${sample}_homer.log
    
    # Convert HOMER peaks to narrowPeak format
    # HOMER format: PeakID chr start end strand Normalized Tag Count focus ratio
    pos2bed.pl ${sample}_tags/peaks.txt > ${sample}_peaks.bed
    
    # Convert to narrowPeak (chr start end name score strand signalValue pValue qValue peak)
    awk 'BEGIN{OFS="\t"} NR>1 {
        print \$2, \$3, \$4, \$1, \$6, \$5, \$6, -1, -1, int((\$4-\$3)/2)
    }' ${sample}_tags/peaks.txt > ${sample}_peaks.narrowPeak
    
    # Copy stats
    cp ${sample}_tags/peaks.txt ${sample}_peaks.txt
    
    # Count peaks
    num_peaks=\$(wc -l < ${sample}_peaks.narrowPeak)
    echo "SUCCESS: Called \$num_peaks peaks with HOMER" >> ${sample}_homer.log
    
    # Cleanup
    rm -rf ${sample}_tags
    """

    stub:
    """
    touch ${sample}_peaks.narrowPeak
    touch ${sample}_peaks.txt
    touch ${sample}_homer.log
    """
}

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
    if [ -s ${peaks} ]; then
        annotatePeaks.pl ${peaks} mm10 > ${cell_type}_annotated.txt
    else
        echo "# No peaks to annotate" > ${cell_type}_annotated.txt
    fi
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
    mkdir -p ${cell_type}_motifs
    
    if [ -s ${peaks} ] && [ \$(wc -l < ${peaks}) -gt 10 ]; then
        findMotifsGenome.pl ${peaks} mm10 ${cell_type}_motifs/ \
            -size 200 -mask
    else
        echo "<html><body>No significant peaks for motif analysis</body></html>" > ${cell_type}_motifs/homerResults.html
    fi
    """

    stub:
    """
    mkdir -p ${cell_type}_motifs
    touch ${cell_type}_motifs/homerResults.html
    """
}