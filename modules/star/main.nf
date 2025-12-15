#!/usr/bin/env nextflow

process STAR_INDEX {
    container 'ghcr.io/bf528/star:latest'
    publishDir "${params.outdir}/star_index", mode: 'copy'
    label 'process_high'

    input:
    path(fasta)
    path(gtf)

    output:
    path("star_index/"), emit: index

    script:
    """
    mkdir -p star_index
    
    STAR --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang 99 \
        --runThreadN ${task.cpus}
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SA
    """
}

process STAR_ALIGN {
    container 'ghcr.io/bf528/star:latest'
    publishDir "${params.outdir}/star_aligned", mode: 'copy'
    label 'process_high'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(reads)
    path(index)

    output:
    tuple val(sample), val(cell_type), val(condition), val(replicate), 
          path("${sample}_Aligned.sortedByCoord.out.bam"), 
          path("${sample}_Aligned.sortedByCoord.out.bam.bai"), emit: bam
    path("${sample}_Log.final.out"), emit: log
    path("${sample}_ReadsPerGene.out.tab"), emit: counts

    script:
    """
    STAR --runMode alignReads \
        --genomeDir ${index} \
        --readFilesIn ${reads} \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --outFileNamePrefix ${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMunmapped Within \
        --outSAMattributes Standard
    
    # Index BAM file
    samtools index -@ ${task.cpus} ${sample}_Aligned.sortedByCoord.out.bam
    """

    stub:
    """
    touch ${sample}_Aligned.sortedByCoord.out.bam
    touch ${sample}_Aligned.sortedByCoord.out.bam.bai
    touch ${sample}_Log.final.out
    touch ${sample}_ReadsPerGene.out.tab
    """
}