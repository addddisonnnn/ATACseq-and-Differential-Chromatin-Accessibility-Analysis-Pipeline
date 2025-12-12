#!/bin/bash
# Download mouse reference genome and GTF

cd reference

echo "Downloading mouse reference genome (GRCm39/mm39)..."
wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa

echo "Downloading mouse GTF annotation..."
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
gunzip Mus_musculus.GRCm39.110.gtf.gz
mv Mus_musculus.GRCm39.110.gtf genes.gtf

echo "Reference files downloaded!"