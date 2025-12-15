#!/usr/bin/env nextflow

process DESEQ2_RNA {
    container 'ghcr.io/bf528/biopython:latest'
    publishDir "${params.outdir}/deseq2_rna", mode: 'copy'
    label 'process_low'

    input:
    tuple val(cell_type), val(samples), val(conditions), val(replicates), path(counts)

    output:
    tuple val(cell_type), path("${cell_type}_differential_genes.txt"), emit: diff_genes
    path("${cell_type}_deseq2_results.txt"), emit: results

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    from scipy import stats
    from scipy.stats import false_discovery_control
    
    # Parse inputs
    samples = "${samples}".split(',')
    conditions = "${conditions}".split(',')
    count_files = "${counts}".split()
    
    # Load count data from featureCounts or STAR output
    count_dict = {}
    for sample, count_file in zip(samples, count_files):
        # Try featureCounts format first
        try:
            df = pd.read_csv(count_file, sep='\t', comment='#', skiprows=1)
            if 'Geneid' in df.columns:
                count_dict[sample] = df.set_index('Geneid').iloc[:, 5]  # 6th column is counts
            else:
                raise ValueError("Not featureCounts format")
        except:
            # Try STAR ReadsPerGene format
            try:
                df = pd.read_csv(count_file, sep='\t', header=None, 
                               names=['gene', 'unstranded', 'forward', 'reverse'])
                # Skip first 4 lines (N_unmapped, etc.)
                df = df[4:].copy()
                count_dict[sample] = df.set_index('gene')['unstranded'].astype(int)
            except:
                print(f"ERROR: Could not parse {count_file}")
                continue
    
    # Create count matrix
    count_matrix = pd.DataFrame(count_dict)
    
    # Filter low count genes
    count_matrix = count_matrix.loc[(count_matrix >= 10).sum(axis=1) >= 3]
    
    # Separate by condition
    wt_samples = [s for s, c in zip(samples, conditions) if c == 'WT']
    ko_samples = [s for s, c in zip(samples, conditions) if c == 'KO']
    
    print(f"WT samples: {wt_samples}")
    print(f"KO samples: {ko_samples}")
    
    # Calculate means
    wt_mean = count_matrix[wt_samples].mean(axis=1)
    ko_mean = count_matrix[ko_samples].mean(axis=1)
    
    # Log2 fold change
    log2fc = np.log2((ko_mean + 1) / (wt_mean + 1))
    
    # T-test
    pvals = []
    for gene in count_matrix.index:
        wt_vals = count_matrix.loc[gene, wt_samples].values
        ko_vals = count_matrix.loc[gene, ko_samples].values
        
        if wt_vals.std() == 0 and ko_vals.std() == 0:
            pvals.append(1.0)
        else:
            _, pval = stats.ttest_ind(wt_vals, ko_vals)
            pvals.append(pval)
    
    # Create results dataframe
    results = pd.DataFrame({
        'gene': count_matrix.index,
        'wt_mean': wt_mean.values,
        'ko_mean': ko_mean.values,
        'log2FoldChange': log2fc.values,
        'pvalue': pvals
    })
    
    # FDR correction
    results['padj'] = false_discovery_control(results['pvalue'].values)
    
    # Save results
    results.to_csv('${cell_type}_differential_genes.txt', sep='\t', index=False)
    
    # Summary
    n_up = ((results['padj'] < 0.05) & (results['log2FoldChange'] > 1)).sum()
    n_down = ((results['padj'] < 0.05) & (results['log2FoldChange'] < -1)).sum()
    
    with open('${cell_type}_deseq2_results.txt', 'w') as f:
        f.write(f"Cell type: ${cell_type}\\n")
        f.write(f"Total genes: {len(results)}\\n")
        f.write(f"Upregulated (padj<0.05, log2FC>1): {n_up}\\n")
        f.write(f"Downregulated (padj<0.05, log2FC<-1): {n_down}\\n")
    
    print(f"Analysis complete: {n_up} up, {n_down} down")
    """

    stub:
    """
    touch ${cell_type}_differential_genes.txt
    touch ${cell_type}_deseq2_results.txt
    """
}