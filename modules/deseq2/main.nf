#!/usr/bin/env nextflow

process DIFFERENTIAL_ACCESSIBILITY {
    container 'ghcr.io/bf528/bedtools_samtools:latest'
    publishDir "${params.outdir}/differential", mode: 'copy'
    label 'process_high'
    
    input:
    path(count_matrix)
    
    output:
    path("significant_up_peaks.bed"), emit: significant_peaks
    path("significant_down_peaks.bed")
    path("all_results.txt"), emit: results
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    
    // Read count matrix
    counts = pd.read_csv('${count_matrix}', sep='\\t', header=None)
    
    // Extract peak coordinates (first 3 columns)
    coords = counts.iloc[:, 0:3]
    coords.columns = ['chr', 'start', 'end']
    
    // Extract count data (remaining columns)
    count_data = counts.iloc[:, 3:]
    
    // Simple fold change analysis (replace with DESeq2 in R for production)
    // Assuming columns are: cDC1_WT_1, cDC1_WT_2, cDC1_KO_1, cDC1_KO_2, 
    //                        cDC2_WT_1, cDC2_WT_2, cDC2_KO_1, cDC2_KO_2
    
    // Calculate means for WT and KO
    wt_mean = count_data.iloc[:, [0,1,4,5]].mean(axis=1)
    ko_mean = count_data.iloc[:, [2,3,6,7]].mean(axis=1)
    
    // Calculate log2 fold change
    log2fc = np.log2((ko_mean + 1) / (wt_mean + 1))
    
    // Filter for significant peaks (|log2FC| > 1)
    results = pd.concat([coords, wt_mean, ko_mean, log2fc], axis=1)
    results.columns = ['chr', 'start', 'end', 'wt_mean', 'ko_mean', 'log2fc']
    results.to_csv('all_results.txt', sep='\\t', index=False)
    
    // Get significant peaks
    up_peaks = results[results['log2fc'] > 1]
    down_peaks = results[results['log2fc'] < -1]
    
    up_peaks[['chr', 'start', 'end']].to_csv('significant_up_peaks.bed', 
                                              sep='\\t', index=False, header=False)
    down_peaks[['chr', 'start', 'end']].to_csv('significant_down_peaks.bed', 
                                                sep='\\t', index=False, header=False)
    """
    
    stub:
    """
    touch significant_up_peaks.bed
    touch significant_down_peaks.bed
    touch all_results.txt
    """
}