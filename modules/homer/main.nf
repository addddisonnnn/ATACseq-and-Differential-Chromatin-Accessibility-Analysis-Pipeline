#!/usr/bin/env nextflow

process ANNOTATE_PEAKS {
    container 'ghcr.io/bf528/homer:latest'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    label 'process_low'

    input:
    tuple val(cell_type), path(peaks)

    output:
    path("${cell_type}_annotated.txt"), emit: annotations

    script:
    """
    annotatePeaks.pl ${peaks} mm10 > ${cell_type}_annotated.txt
    """

    stub:
    """
    touch ${cell_type}_annotated.txt
    """
}

process MOTIF_ANALYSIS {
    container 'ghcr.io/bf528/homer:latest'
    publishDir "${params.outdir}/motifs", mode: 'copy'
    label 'process_high'  // Changed to high since this is memory intensive
    
    // Add genome FASTA as input
    input:
    tuple val(cell_type), path(peaks)
    path(fasta)
    
    output:
    path("${cell_type}_motifs/"), emit: motif_results
    
    script:
    """
    echo "=== Preparing mm39 genome for HOMER ==="
    
    # Check if prepareGenome.pl exists
    if [ -f /opt/conda/share/homer/.//prepareGenome.pl ]; then
        echo "Using prepareGenome.pl..."
        # Prepare genome for HOMER
        perl /opt/conda/share/homer/.//prepareGenome.pl \\
            "${fasta}" \\
            prepared_mm39
        
        GENOME="prepared_mm39"
    else
        echo "prepareGenome.pl not found, creating custom genome directory..."
        # Create custom genome directory structure
        mkdir -p homer_genome/mm39
        cp "${fasta}" homer_genome/mm39/genome.fa
        
        # Create minimal configuration
        echo "name mm39" > homer_genome/mm39/genome.config
        echo "organism Mouse" >> homer_genome/mm39/genome.config
        echo "assembly GRCm39" >> homer_genome/mm39/genome.config
        
        GENOME="homer_genome/mm39"
    fi
    
    echo "=== Running motif analysis ==="
    echo "Using genome: \$GENOME"
    echo "Peaks file: ${peaks}"
    echo "Output directory: ${cell_type}_motifs/"
    
    # Run motif analysis
    findMotifsGenome.pl "${peaks}" "\$GENOME" "${cell_type}_motifs/" \\
        -size 200 \\
        -mask
    
    echo "=== Motif analysis complete ==="
    """
    
    stub:
    """
    mkdir -p ${cell_type}_motifs
    touch ${cell_type}_motifs/homerResults.html
    touch ${cell_type}_motifs/homerResults.txt
    """
}