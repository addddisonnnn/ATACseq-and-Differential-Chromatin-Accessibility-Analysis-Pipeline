#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules - REUSING existing modules where possible
include { DOWNLOAD_SRA } from './modules/download/main.nf'
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/trimmomatic/main.nf'
include { STAR_INDEX } from './modules/star/main.nf'
include { STAR_ALIGN } from './modules/star/main.nf'
include { DESEQ2_RNA } from './modules/deseq2/main_rna.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

workflow {
    // Parse samplesheet
    Channel
        .fromPath(params.rna_samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            tuple(
                row.sample,
                row.srr,
                row.cell_type,
                row.condition,
                row.replicate
            )
        }
        .set { samples_ch }

    // 1. Download reads
    DOWNLOAD_SRA(samples_ch)
    
    // 2. QC raw reads
    FASTQC_RAW(DOWNLOAD_SRA.out.reads)
    
    // 3. Trim adapters
    TRIMMOMATIC(DOWNLOAD_SRA.out.reads)
    
    // 4. QC trimmed reads
    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads)
    
    // 5. Build STAR index (only once)
    STAR_INDEX(
        file(params.genome_fasta),
        file(params.genome_gtf)
    )
    
    // 6. Align with STAR (includes read counting)
    STAR_ALIGN(
        TRIMMOMATIC.out.trimmed_reads,
        STAR_INDEX.out.index
    )
    
    // 7. Group by cell type for differential expression
    STAR_ALIGN.out.counts
        .map { sample, cell_type, condition, replicate, counts ->
            tuple(cell_type, sample, condition, replicate, counts)
        }
        .groupTuple()
        .set { deseq_input }
    
    // 8. Differential expression analysis
    DESEQ2_RNA(deseq_input)
    
    // 9. MultiQC report
    FASTQC_RAW.out.fastqc_zip
        .mix(FASTQC_TRIMMED.out.fastqc_zip)
        .mix(STAR_ALIGN.out.log)
        .collect()
        .set { multiqc_input }
    
    MULTIQC(multiqc_input)
}