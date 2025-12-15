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
    findPeaks ${sample}_tags -style factor -size 150 -minDist 150 -o auto 2>&1 | tee -a ${sample}_homer.log
    
    # Convert HOMER peaks to narrowPeak format
    grep -v "^#" ${sample}_tags/peaks.txt | \
    awk 'BEGIN{OFS="\t"} 
         \$1 != "PeakID" && \$1 != "" && NF >= 6 {
             chr = \$2
             start = \$3
             end = \$4
             name = \$1
             score = (\$6 > 1000) ? 1000 : int(\$6)
             strand = \$5
             signal = (NF >= 7) ? \$7 : \$6
             summit = int((end - start) / 2)
             
             print chr, start, end, name, score, strand, signal, -1, -1, summit
         }' > ${sample}_peaks.narrowPeak
    
    cp ${sample}_tags/peaks.txt ${sample}_peaks.txt
    
    num_peaks=\$(grep -v "^#" ${sample}_peaks.narrowPeak | wc -l)
    echo "SUCCESS: Called \$num_peaks peaks with HOMER" >> ${sample}_homer.log
    
    if [ "\$num_peaks" -eq 0 ]; then
        echo "WARNING: No peaks in final narrowPeak file!" >> ${sample}_homer.log
    fi
    
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
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    label 'process_low'

    input:
    tuple val(cell_type), path(peaks)

    output:
    path("${cell_type}_annotated.txt"), emit: annotations

    script:
    """
    if [ -s ${peaks} ] && [ \$(wc -l < ${peaks}) -gt 0 ]; then
        # Use genome fasta instead of pre-installed genome
        annotatePeaks.pl ${peaks} ${params.genome_fasta} -gtf ${params.genome_gtf} > ${cell_type}_annotated.txt
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
    container 'ghcr.io/bf528/homer_samtools:latest'
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
        # Use genome fasta instead of pre-installed genome
        findMotifsGenome.pl ${peaks} ${params.genome_fasta} ${cell_type}_motifs/ \
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