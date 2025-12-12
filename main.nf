#!/usr/bin/env nextflow

include { DOWNLOAD_SRA } from './modules/sra/main.nf'

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
}

