#!/usr/bin/env nextflow

include { DOWNLOAD_SRA } from './modules/sra/main.nf'
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/trimmomatic/main.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc/main.nf'
include { BOWTIE2_BUILD } from './modules/bowtie2/build.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2/align.nf'
include { REMOVE_MITO } from './modules/samtools/remove_mito.nf'
include { MACS3_CALLPEAK } from './modules/macs3/main.nf'
include { CREATE_COUNT_MATRIX } from './modules/bedtools/count_matrix.nf'
include { DIFFERENTIAL_ACCESSIBILITY } from './modules/deseq2/main.nf'
include { ANNOTATE_PEAKS } from './modules/homer/annotate.nf'
include { FIND_MOTIFS } from './modules/homer/motifs.nf'

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

    // 8. Peak calling with MACS3 - group by condition
    REMOVE_MITO.out.bam
        .map { sample_id, condition, replicate, bam, bai -> 
            tuple(condition, bam) 
        }
        .groupTuple()
        .set { grouped_bams }

    MACS3_CALLPEAK(grouped_bams)
    
    // 9. Create count matrix for differential analysis
    // First, prepare peaks channel - collect all peaks files
    MACS3_CALLPEAK.out.narrowPeak
        .collect()
        .set { all_peaks }

    // Second, prepare BAMs channel - we need to map to include sample info
    REMOVE_MITO.out.bam
        .map { sample_id, condition, replicate, bam, bai -> 
            tuple(sample_id, condition, replicate, bam) 
        }
        .collect()
        .set { all_bams }

    CREATE_COUNT_MATRIX(all_peaks, all_bams)

    // 10. Differential accessibility analysis
    DIFFERENTIAL_ACCESSIBILITY(CREATE_COUNT_MATRIX.out.counts)

    // 11. Annotate differentially accessible peaks
    ANNOTATE_PEAKS(
        DIFFERENTIAL_ACCESSIBILITY.out.significant_peaks,
        Channel.fromPath(params.gtf)
    )

    // 12. Motif finding on differentially accessible peaks
    FIND_MOTIFS(
        DIFFERENTIAL_ACCESSIBILITY.out.significant_peaks,
        genome_fasta
    )
}
