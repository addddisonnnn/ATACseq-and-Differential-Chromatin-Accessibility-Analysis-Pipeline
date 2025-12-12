#!/usr/bin/env nextflow

process BOWTIE2_BUILD {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir "${params.outdir}/reference", mode: 'copy'
    label 'process_high'
    
    input:
    path(genome_fasta)
    
    output:
    path("genome_index"), emit: index
    
    script:
    """
    mkdir -p genome_index
    bowtie2-build --threads ${task.cpus} ${genome_fasta} genome_index/genome
    """
    
    stub:
    """
    mkdir -p genome_index
    touch genome_index/genome.1.bt2
    touch genome_index/genome.2.bt2
    touch genome_index/genome.3.bt2
    touch genome_index/genome.4.bt2
    touch genome_index/genome.rev.1.bt2
    touch genome_index/genome.rev.2.bt2
    """
}