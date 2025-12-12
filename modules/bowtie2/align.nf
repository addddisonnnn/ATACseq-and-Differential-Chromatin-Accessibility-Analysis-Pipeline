#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir "${params.outdir}/aligned", mode: 'copy'
    label 'process_high'
    
    input:
    tuple val(sample_id), val(condition), val(replicate), path(reads), path(index)
    
    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    path("${sample_id}_align.log"), emit: log
    
    script:
    """
    # Align with bowtie2
    bowtie2 -p ${task.cpus} \
        --very-sensitive \
        -X 2000 \
        -x ${index}/genome \
        -U ${reads} \
        2> ${sample_id}_align.log | \
        samtools view -@ ${task.cpus} -bS -q 30 - > ${sample_id}_unsorted.bam
    
    # Check if alignment succeeded
    if [ ! -s ${sample_id}_unsorted.bam ]; then
        echo "ERROR: Alignment failed or produced empty BAM"
        exit 1
    fi
    
    # Sort BAM
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam ${sample_id}_unsorted.bam
    
    # Index BAM
    samtools index ${sample_id}.bam
    
    # Clean up
    rm ${sample_id}_unsorted.bam
    """
}