#!/bin/bash
# Download and prepare mouse reference genome files for mm10

mkdir -p reference
cd reference

module load samtools
module load bowtie2

echo "Downloading mouse reference genome (GRCm38/mm10)..."
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
mv GRCm38.primary_assembly.genome.fa genome.fa

echo "Downloading mouse GTF annotation..."
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
gunzip gencode.vM25.primary_assembly.annotation.gtf.gz
mv gencode.vM25.primary_assembly.annotation.gtf genes.gtf

echo "Creating chromosome sizes file..."
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > genome.chrom.sizes

echo "Downloading ENCODE blacklist regions for mm10..."
wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz
gunzip mm10-blacklist.v2.bed.gz

echo "Building Bowtie2 index (this will take ~30 minutes)..."
mkdir -p bowtie2_index
bowtie2-build genome.fa bowtie2_index/GRCm38

echo "Reference files ready!"