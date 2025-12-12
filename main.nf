#!/usr/bin/env nextflow

include {DOWNLOAD_SRA} from './modules/download_sra'

workflow {
    // Create channel from CSV with ATAC-seq sample information
    Channel.fromPath(params.atac_samples)
    | splitCsv(header: true)
    | map { row -> tuple(row.sample_id, row.srr_id) }
    | set { atac_ch }

    // Download ATAC-seq files from SRA
    DOWNLOAD_SRA(atac_ch)
}

