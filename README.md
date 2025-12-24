# ATACseq and Differential Chromatin Accessibility Analysis Pipeline
Implementing a bioinformatics pipeline and performing analysis on a publication focused on chromatin accessibility as assayed by ATACseq.

## Project Overview
This repository implements a comprehensive bioinformatics pipeline for analyzing chromatin accessibility data from ATAC-seq experiments, specifically replicating the analysis from "The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity" (De Sá Fernandes et al., 2024, Cell Reports 43:114308).
## Scientific Background
The study investigates how the epigenetic regulator HDAC1 controls dendritic cell (DC) development and function. Using ATAC-seq, the authors demonstrated that HDAC1 deletion in hematopoietic progenitors alters chromatin accessibility in DC subsets, impairing cDC2 and pDC development while enhancing anti-tumor immunity through increased cDC1 activation. This pipeline reproduces their ATAC-seq analysis to validate these epigenetic findings.
## ATAC-seq Methodology
ATAC-seq (Assay for Transposase-Accessible Chromatin with sequencing) profiles genome-wide chromatin accessibility using a hyperactive Tn5 transposase that inserts sequencing adapters into accessible chromatin regions. This technique identifies open chromatin regions, transcription factor binding sites, and nucleosome positions with single-nucleotide resolution, making it ideal for studying epigenetic regulation in immune cell development.

## Pipeline Workflow
The pipeline implements the following analysis steps in Nextflow:
1. Data Acquisition: Download ATAC-seq data from SRA (GSE266584)
2. Quality Control: FastQC analysis of raw and trimmed reads
3. Read Processing: Adapter trimming and quality filtering with Trimmomatic
4. Alignment: Genome alignment with Bowtie2 (mm10/GRCm38)
5. Peak Calling: Identify open chromatin regions with HOMER
6. Differential Analysis: Detect differentially accessible regions (DARs) with DESeq2
7. Peak Annotation: Annotate DARs to genomic features using HOMER
8. Motif Analysis: Discover enriched transcription factor motifs
9. Visualization: Generate coverage tracks and quality metrics with deepTools

## Prerequisites
### Software Requirements
1.  **Nextflow:** (version 22.04 or later recommended)
2.  **Container Runtime:** Singularity (recommended for HPC/SGE) or Docker.
3.  **Reference Genome:** Bowtie2 index files and Salmon index for the mouse genome (mm10/GRCm38).

### System Requirements
1. Memory: 32+ GB RAM recommended for peak calling steps
2. Storage: 50+ GB free disk space
3. CPUs: 8+ cores for efficient parallel processing
## Installation and Setup

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/addddisonnnn/ATACseq-and-Differential-Chromatin-Accessibility-Analysis-Pipeline.git
    cd ATACseq-and-Differential-Chromatin-Accessibility-Analysis-Pipeline 
    ```
2.  **Configuration:** Edit `nextflow.config` to specify your execution environment.
```groovy
// For SGE cluster execution (recommended):
process.executor = 'sge'
process.clusterOptions = '-P bf528 -l h_rt=24:00:00'

// For local execution:
// process.executor = 'local'
// process.cpus = 8
// process.memory = '32.GB'

// Set reference genome paths
params {
    genome_fasta = "/path/to/reference/mm10/genome.fa"
    bowtie2_index = "/path/to/bowtie2/mm10_index/mm10"
    genome_gtf = "/path/to/annotation/gencode.vM25.annotation.gtf"
    chrom_sizes = "/path/to/mm10.chrom.sizes"
}

```
3.  **Reference Data:** Ensure the following paths are correctly set in your configuration file:
```bash
# Download mm10 reference
mkdir -p reference
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# Build Bowtie2 index
bowtie2-build mm10.fa mm10

# Download annotation (optional, for visualization)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz
```
4. **Prepare Sample Sheet**
Edit `atac_samples.csv` with your sample information:
```csv
sample,srr,cell_type,condition,replicate
cDC1_WT_1,SRR28895190,cDC1,WT,1
cDC2_WT_1,SRR28895189,cDC2,WT,1
cDC1_KO_1,SRR28895188,cDC1,KO,1
cDC2_KO_1,SRR28895187,cDC2,KO,1
cDC1_WT_2,SRR28895186,cDC1,WT,2
cDC2_WT_2,SRR28895185,cDC2,WT,2
cDC1_KO_2,SRR28895184,cDC1,KO,2
cDC2_KO_2,SRR28895183,cDC2,KO,2
```

## Running the Pipeline
### Basic Execution
```bash
# Run with default parameters
nextflow run main.nf -profile singularity

# Run with specific output directory
nextflow run main.nf -profile singularity --outdir results/hdac1_analysis

# Resume interrupted run
nextflow run main.nf -profile singularity -resume
```
## Execution Profiles
The pipeline supports multiple execution profiles:
```bash
# SGE cluster with Singularity (recommended for HPC)
nextflow run main.nf -profile cluster

# Local execution with Singularity
nextflow run main.nf -profile singularity

# Local execution with Conda environments
nextflow run main.nf -profile conda

# Test run with minimal data
nextflow run main.nf -profile test
```
## Monitor Executing
```bash
# View execution timeline
nextflow log -f timeline

# Check resource usage
nextflow log -f resources

# Generate execution report
nextflow run main.nf -profile singularity -with-report execution_report.html
```
## Project Structure
```
ATACseq-and-Differential-Chromatin-Accessibility-Analysis-Pipeline
├── nextflow.config         # Configuration
├── main.nf                 # Main workflow
├── atac_samples.csv        # Sample sheet with SRR IDs
├── bin/
│   └── download_references # Download reference
├── modules/
│   ├── sra/main.nf         # Download FASTQ
│   ├── fastqc/main.nf      # QC
│   ├── trimmomatic/main.nf # Trimming
│   ├── bowtie2/
│   │   ├── build.nf        # Index building
│   │   └── align.nf        # Alignment
│   ├── samtools/remove_mito.nf
│   ├── macs3/main.nf       # Peak calling
│   ├── bedtools/count_matrix.nf
│   ├── deseq2/main.nf      # Differential analysis
│   └── homer/
│       ├── annotate.nf     # Peak annotation
│       └── motifs.nf       # Motif finding
├── reference/
│   ├── bowtie2_index/      # Indexed genome
│   ├── genome.fa           # Reference Gneome (mm10/GRCm38)
│   ├── genes.gtf           # GTF Annotation
│   └── genome.chrom.sizes  # Chromosome size
└── results/                # Output directory
```
## Output Structure
```text
results/
├── fastqc/
│   ├── raw/           # Raw read quality reports
│   └── trimmed/       # Trimmed read quality reports
├── alignment/
│   ├── bam/           # Filtered alignment files
│   ├── bigwig/        # Normalized coverage tracks
│   └── stats/         # Alignment statistics
├── peaks/
│   ├── per_sample/    # Individual peak calls
│   ├── merged/        # Merged peak sets
│   └── counts/        # Peak count matrices
├── differential/
│   ├── deseq2/        # Differential accessibility results
│   ├── annotated/     # Peak annotations
│   └── motifs/        # Motif enrichment results
├── qc/
│   ├── multiqc/       # MultiQC summary report
│   ├── tss_enrichment/# TSS enrichment plots
│   └── frip/          # FRiP score calculations
└── figures/
    ├── heatmaps/      # Coverage heatmaps
    ├── browser_tracks/# Genome browser tracks
    └── correlation/   # ATAC-RNA integration plots
```

## Key Results and Interpretation

### Quality Control Metrics
The pipeline generates comprehensive QC reports including:

- FastQC reports: Per-base quality scores, GC content, adapter contamination
- Alignment statistics: Mapping rates, read distribution, mitochondrial content
- ATAC-specific metrics: TSS enrichment scores (>7 indicates high quality), FRiP scores (14-19% expected range)
### Differential Accessibility Analysis
Output includes:
- DAR lists: Differentially accessible regions for cDC1 and cDC2 (adjusted p-value < 0.15, |log2FC| > 0.5)
- Heatmap and siganl plots: Visual representation of accessibility changes
- Volcano plots: Significance vs. fold change distributions
### Functional Analysis
- Motif enrichment: Transcription factor binding sites enriched in DARs (PU.1, SPIB, IRF8 expected)
- Genomic annotation: Distribution of DARs across promoters, exons, introns, intergenic regions
- Integration with RNA-seq: Correlation plots linking chromatin accessibility to gene expression

## Reproducing Publication Figures
To reproduce specific figures from De Sá Fernandes et al. (2024):
```bash
# Generate Figure 6A/B: Summary histograms of DARs
python scripts/replicate_figures.py --figure 6AB

# Generate Figure 6C/E: ATAC-RNA correlation plots
python scripts/replicate_figures.py --figure 6CE

# Generate Figure 6D/F: Genome browser tracks
python scripts/visualize_loci.py --locus Maged1 --locus Spib
```
## Troubleshooting
### Common Issues and Solutions
1. Memory errors during peak calling:
```bash
# Increase memory allocation in nextflow.config
process {
    withName: 'HOMER_CALLPEAK' {
        memory = '64.GB'
        time = '24.h'
    }
}
```
2. Slow SRA downloads:
```bash
# Use prefetch with fasterq-dump
prefetch SRR28895190 --max-size 100G
fasterq-dump SRR28895190 --split-files
```
3. Missing reference files:
```bash
# Run the reference download script
bash bin/download_references.sh
```
### Debug Mode
```bash
# Run with debug output
nextflow run main.nf -profile singularity -debug

# Test individual processes
nextflow run main.nf -profile singularity -process HOMER_CALLPEAK
```

## Citation
De Sá Fernandes, C., Novoszel, P., Gastaldi, T., Krauß, D., Lang, M., Rica, R., Kutschat, A. P., Holcmann, M., Ellmeier, W., Seruggia, D., Strobl, H., & Sibilia, M. (2024). The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity. Cell reports, 43(6), 114308. https://doi.org/10.1016/j.celrep.2024.114308

```bibtex
@article{desafernandes2024hdac1,
  title={The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity},
  author={De S{\'a} Fernandes, Cristiano and Novoszel, Peter and Gastaldi, Tommaso and Krau{\ss}, Dagmar and Pone, Ilenia and Zikeli, Sabeth and Wache, Harald and Tonn, Theresa and Schlederer, Michael and Kenner, Lukas and others},
  journal={Cell Reports},
  volume={43},
  number={6},
  pages={114308},
  year={2024},
  publisher={Elsevier},
  doi={10.1016/j.celrep.2024.114308}
}
```
## Support and Contributing
- Issues: Report bugs or problems via GitHub Issues
- Questions: Use GitHub Discussions for usage questions
- Contributions: Submit pull requests for improvements or bug fixes

## License
This pipeline is released under the MIT License. See `LICENSE` for details.
