#!/usr/bin/env nextflow

process FASTQC {
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    label 'process_low'

    input:
    tuple val(sample), val(cell_type), val(condition), val(replicate), path(reads)

    output:
    path("*_fastqc.{zip,html}"), emit: fastqc_results
    path("*_fastqc.zip"), emit: fastqc_zip

    script:
    """
    # Set temporary directory to current work directory
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR
    
    # Also set Java temp directory
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\$TMPDIR"
    
    fastqc -t ${task.cpus} --dir \$TMPDIR ${reads}
    """

    stub:
    """
    touch ${sample}_fastqc.html
    touch ${sample}_fastqc.zip
    """
}