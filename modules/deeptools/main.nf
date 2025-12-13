#!/usr/bin/env nextflow

process BIGWIG_COVERAGE {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bam), path(bai)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path("${sample}.bw"), emit: bigwig

    script:
    """
    bamCoverage -b ${bam} -o ${sample}.bw \
        -p ${task.cpus} \
        --normalizeUsing CPM \
        --binSize 10
    """

    stub:
    """
    touch ${sample}.bw
    """
}

process COMPUTE_MATRIX {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/deeptools", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(cell_type), val(condition), path(bigwigs)

    output:
    tuple val(cell_type), val(condition), path("${cell_type}_${condition}_matrix.gz"), emit: matrix

    script:
    """
    computeMatrix reference-point \
        -S ${bigwigs.join(' ')} \
        -R ${params.genome_gtf} \
        --referencePoint TSS \
        -b 3000 -a 3000 \
        -p ${task.cpus} \
        --skipZeros \
        -o ${cell_type}_${condition}_matrix.gz
    """

    stub:
    """
    touch ${cell_type}_${condition}_matrix.gz
    """
}

process PLOT_HEATMAP {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/figures", mode: 'copy'
    label 'process_low'

    input:
    tuple val(cell_type), val(condition), path(matrix)

    output:
    path("${cell_type}_${condition}_heatmap.png"), emit: heatmap
    path("${cell_type}_${condition}_profile.png"), emit: profile

    script:
    """
    plotHeatmap -m ${matrix} \
        -o ${cell_type}_${condition}_heatmap.png \
        --colorMap RdBu \
        --whatToShow 'heatmap and colorbar'
    
    plotProfile -m ${matrix} \
        -o ${cell_type}_${condition}_profile.png \
        --perGroup \
        --plotTitle "${cell_type} ${condition} TSS Profile"
    """

    stub:
    """
    touch ${cell_type}_${condition}_heatmap.png
    touch ${cell_type}_${condition}_profile.png
    """
}

process TSS_ENRICHMENT {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/qc", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(bigwig)

    output:
    path("${sample}_TSS_enrichment.png"), emit: plot
    path("${sample}_TSS_matrix.gz"), emit: matrix

    script:
    """
    computeMatrix reference-point \
        -S ${bigwig} \
        -R ${params.genome_gtf} \
        --referencePoint TSS \
        -b 2000 -a 2000 \
        -p ${task.cpus} \
        -o ${sample}_TSS_matrix.gz
    
    plotProfile -m ${sample}_TSS_matrix.gz \
        -o ${sample}_TSS_enrichment.png \
        --perGroup
    """

    stub:
    """
    touch ${sample}_TSS_enrichment.png
    touch ${sample}_TSS_matrix.gz
    """
}