#!/usr/bin/env nextflow

include { DOWNLOAD_SRA } from './modules/sra/main.nf'
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/trimmomatic/main.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc/main.nf'
include { BOWTIE2_BUILD } from './modules/bowtie2/build.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2/align.nf'
include { REMOVE_MITO } from './modules/samtools/remove_mito.nf'


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

    // 4. QC on trimmed reads
    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads)

    // 5. Build Bowtie2 index
    genome_fasta = Channel.fromPath(params.genome)
    BOWTIE2_BUILD(genome_fasta)

    // 6. Align with ATAC-seq specific parameters
    TRIMMOMATIC.out.trimmed_reads
        .combine(BOWTIE2_BUILD.out.index)
        .set { reads_with_index }
    BOWTIE2_ALIGN(reads_with_index)

    // 7. Remove mitochondrial reads
    REMOVE_MITO(BOWTIE2_ALIGN.out.bam)
}
