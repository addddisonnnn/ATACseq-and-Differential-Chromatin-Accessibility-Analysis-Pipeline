# ATACseq-and-Differential-Chromatin-Accessibility-Analysis-Pipeline
Implementing a bioinformatics pipeline and an analysis on a publication focused on chromatin accessiblity as assayed by ATACseq.

## Project Overview

### Understanding HDAC1's Role in Dendritic Cell Development and Cancer Immunity
Dendritic Cells (DCs) are the immune system's 'conductors', they conduct immune cells, such as T cells, what to attack. The paper discovered how a protein called HDAC1 can act as a 'volume control' for these DCs conductors, and when HDAC1 is removed, suprisingly, the immune system gets better at fighting tumors. 

### ATAC-Seq
ATAC-Seq which stands for Assay for Transposase-Accessible Chromatin with high-throughput sequencing is an established technique for identifying regions of open and closed chromatin via a transposase enzyme, called Tn5. This enzyme inserts sequencing adapters into accessible spots and creates a library of tagged genomic regions/fragments. Finally, these are amplified and sequenced. 

## Prerequisites

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
│   ├── genome.fa           # Mouse Reference Gneome
│   └── genes.gtf           # Mouse GTF Annotation
└── results/                # Output directory
```

## Workflow
1. Download ATAC-seq files (SRA Toolkit)
2. Download reference files (wget)
3. Trim adapters - Nextera for ATAC-seq (Trimmomatic)

## Key Results

## Citation

De Sá Fernandes, C., Novoszel, P., Gastaldi, T., Krauß, D., Lang, M., Rica, R., Kutschat, A. P., Holcmann, M., Ellmeier, W., Seruggia, D., Strobl, H., & Sibilia, M. (2024). The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity. Cell reports, 43(6), 114308. https://doi.org/10.1016/j.celrep.2024.114308

Documenting how long it takes it takes for each step/module to run
```
[20/26395f] DOWNLOAD_SRA (7) [100%] 8 of 8 ✔
Completed at: 11-Dec-2025 19:28:56
Duration    : 10m 26s
CPU hours   : 7.9
Succeeded   : 8
```
```
[9a/da0531] FASTQC_RAW (6)   [100%] 8 of 8 ✔
Completed at: 11-Dec-2025 19:36:45
Duration    : 3m 21s
CPU hours   : 8.8 (90% cached)
Succeeded   : 8
Cached      : 8
``` 
```
[76/95e46d] TRIMMOMATIC (3)  [100%] 8 of 8 ✔
Completed at: 11-Dec-2025 19:54:30
Duration    : 9m 16s
CPU hours   : 16.2 (54.2% cached)
Succeeded   : 8
Cached      : 16
```
```
[5c/deb8b1] FASTQC_TRIMMED (2) [100%] 8 of 8, ignored: 6 ✔
[45/1802ab] NOTE: Process `FASTQC_TRIMMED (4)` terminated with an error exit status (1) -- Error is ignored
[e4/e8cb6b] NOTE: Process `FASTQC_TRIMMED (6)` terminated with an error exit status (1) -- Error is ignored
[30/fc7fd2] NOTE: Process `FASTQC_TRIMMED (1)` terminated with an error exit status (1) -- Error is ignored
Completed at: 11-Dec-2025 20:15:34
Duration    : 3m 21s
CPU hours   : 17.0 (95.7% cached, 3% failed)
Succeeded   : 2
Cached      : 24
Ignored     : 6
Failed      : 6
```

```
[03/7a2b52] BOWTIE2_BUILD (1)  [100%] 1 of 1 ✔
Completed at: 11-Dec-2025 21:03:30
Duration    : 43m 29s
CPU hours   : 28.1 (58.5% cached)
Succeeded   : 7
Cached      : 26
```