#!/usr/bin/env nextflow

include { DOWNLOAD_SRA } from './modules/sra/main.nf'
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'

workflow {
    // Create channel from samplesheet with SRR numbers
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            tuple(
                row.sample_id,
                row.condition,
                row.replicate,
                row.srr
            )
        }
        .set { samples_ch }
    
    // 1. Download data from SRA
    DOWNLOAD_SRA(samples_ch)

    // 2. QC on raw reads
    FASTQC_RAW(DOWNLOAD_SRA.out.fastq)

    // 3. Trim adapters (Nextera for ATAC-seq)
    TRIMMOMATIC(DOWNLOAD_SRA.out.fastq)
}

