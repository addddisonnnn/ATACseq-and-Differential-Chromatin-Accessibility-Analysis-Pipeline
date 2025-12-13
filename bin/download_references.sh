#!/bin/bash
# Download and prepare mouse reference genome files

mkdir -p reference
cd reference

module load samtools
module load bowtie2

#echo "Downloading mouse reference genome (GRCm39/mm39)..."
#wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
#gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
#mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa

#echo "Downloading mouse GTF annotation..."
#wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
#gunzip Mus_musculus.GRCm39.110.gtf.gz
#mv Mus_musculus.GRCm39.110.gtf genes.gtf

#echo "Creating chromosome sizes file..."
#samtools faidx genome.fa
#cut -f1,2 genome.fa.fai > genome.chrom.sizes

echo "Downloading ENCODE blacklist regions for mm39..."
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
gunzip mm39-blacklist.v2.bed.gz

echo "Building Bowtie2 index (this will take ~30 minutes)..."
mkdir -p bowtie2_index
bowtie2-build genome.fa bowtie2_index/GRCm39

echo "Reference files ready!"