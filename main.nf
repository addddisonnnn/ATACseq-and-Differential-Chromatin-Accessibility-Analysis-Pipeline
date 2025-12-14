#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { DOWNLOAD_SRA } from './modules/download/main.nf'
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/trimmomatic/main.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2/main.nf'
include { FILTER_ALIGNMENTS } from './modules/samtools/main.nf'
include { HOMER_CALLPEAK } from './modules/homer/main.nf'
include { MERGE_PEAKS } from './modules/bedtools/main.nf'
include { COUNT_PEAKS } from './modules/bedtools/main.nf'
include { DIFF_ANALYSIS } from './modules/deseq2/main.nf'
include { ANNOTATE_PEAKS } from './modules/homer/main.nf'
include { MOTIF_ANALYSIS } from './modules/homer/main.nf'
include { BIGWIG_COVERAGE } from './modules/deeptools/main.nf'
include { COMPUTE_MATRIX } from './modules/deeptools/main.nf'
include { PLOT_HEATMAP } from './modules/deeptools/main.nf'
include { TSS_ENRICHMENT } from './modules/deeptools/main.nf'
include { FRIP_SCORE } from './modules/qc/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

workflow {
    // Parse samplesheet
    Channel
        .fromPath(params.samplesheet)
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

    // 1-5. Download through filtering
    DOWNLOAD_SRA(samples_ch)
    FASTQC_RAW(DOWNLOAD_SRA.out.reads)
    TRIMMOMATIC(DOWNLOAD_SRA.out.reads)
    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads)
    BOWTIE2_ALIGN(TRIMMOMATIC.out.trimmed_reads)
    FILTER_ALIGNMENTS(BOWTIE2_ALIGN.out.bam)

    // 6. Call peaks per sample using HOMER
    HOMER_CALLPEAK(FILTER_ALIGNMENTS.out.filtered_bam)

    // 7. Group ALL peaks by cell_type only (not condition) for merging
    HOMER_CALLPEAK.out.peaks
        .map { sample, cell_type, condition, replicate, peaks -> 
            // Group by cell_type ONLY
            tuple(cell_type, peaks)
        }
        .groupTuple()
        .set { grouped_peaks }

    // Merge all peaks for each cell type (WT + KO together)
    MERGE_PEAKS(grouped_peaks)

    // 8. Count reads in merged peaks for all samples
    // Match each sample with its cell_type's merged peaks
    FILTER_ALIGNMENTS.out.filtered_bam
        .map { sample, cell_type, condition, replicate, bam, bai ->
            tuple(cell_type, sample, condition, replicate, bam, bai)
        }
        .combine(MERGE_PEAKS.out.merged_peaks, by: 0)
        .map { cell_type, sample, condition, replicate, bam, bai, peaks ->
            tuple(sample, cell_type, condition, replicate, bam, peaks)
        }
        .set { count_input }

    COUNT_PEAKS(count_input)

    // 9. Differential accessibility analysis (per cell type)
    COUNT_PEAKS.out.counts
        .map { sample, cell_type, condition, replicate, counts ->
            tuple(cell_type, sample, condition, replicate, counts)
        }
        .groupTuple()
        .set { diff_input }

    DIFF_ANALYSIS(diff_input)

    // 10. Annotate differentially accessible peaks
    ANNOTATE_PEAKS(DIFF_ANALYSIS.out.diff_peaks)

    // 11. Motif analysis on differential peaks
    MOTIF_ANALYSIS(DIFF_ANALYSIS.out.diff_peaks)

    // Generate bigWig files for visualization
    BIGWIG_COVERAGE(FILTER_ALIGNMENTS.out.filtered_bam)

    // Group bigwigs by cell type and condition
    BIGWIG_COVERAGE.out.bigwig
        .map { sample, cell_type, condition, replicate, bigwig ->
            tuple(cell_type, condition, bigwig)
        }
        .groupTuple(by: [0, 1])
        .set { grouped_bigwigs }

    // Compute matrix for heatmaps
    COMPUTE_MATRIX(grouped_bigwigs)

    // Plot heatmaps
    PLOT_HEATMAP(COMPUTE_MATRIX.out.matrix)

    // TSS enrichment per sample
    TSS_ENRICHMENT(BIGWIG_COVERAGE.out.bigwig)

    // Calculate FRiP scores
    FRIP_SCORE(
        FILTER_ALIGNMENTS.out.filtered_bam
            .join(HOMER_CALLPEAK.out.peaks, by: [0,1,2,3])
    )

    // Collect all QC outputs
    FASTQC_RAW.out.fastqc_zip
        .mix(FASTQC_TRIMMED.out.fastqc_zip)
        .mix(BOWTIE2_ALIGN.out.log)
        .mix(FRIP_SCORE.out.frip_stats)
        .collect()
        .set { multiqc_input }
    
    MULTIQC(multiqc_input)
}
/*
CONCEPTUAL MAP
    PARAMS.SAMPLESHEET
        │
    DOWNLOAD_SRA ───────────┐
        │               FASTQC_RAW
    TRIMMOMATIC ────────────┐
        │               FASTQC_TRIMMED
    BOWTIE2_ALIGN
        │
    FILTERED_ALIGNMENTS
        │
    MACS2_CALLPEAK ───────────────────────────────┐
        │                                         │
    MERGE_PEAKS                                   │
        │                                         │
    FILTERED_ALIGNMENTS ─────────┬────────────────┤
        │                        │                │
    COUNT_PEAKS           BIGWIG_COVERAGE    FRIP_SCORE
        │                        │                │
    DIFF_ANALYSIS         COMPUTE_MATRIX       MULTIQC
        │                        │
    ANNOTATE_PEAKS         PLOT_HEATMAP

┌──┬──┐ ╭-╮
├──┼──┤ ╰─╯
└──┴──┘

*/