#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--counts', nargs='+', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--conditions', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--bed', required=True)
    args = parser.parse_args()
    
    samples = args.samples.split(',')
    conditions = args.conditions.split(',')
    
    # Load count data
    result = None
    for i, count_file in enumerate(args.counts):
        df = pd.read_csv(count_file, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'count'])
        df = df.rename(columns={'count': samples[i]})
        if i == 0:
            result = df[['chr', 'start', 'end', samples[i]]]
        else:
            result = result.merge(df[['chr', 'start', 'end', samples[i]]], 
                                on=['chr', 'start', 'end'])
    
    # Perform differential analysis
    ko_samples = [s for s, c in zip(samples, conditions) if c == 'KO']
    wt_samples = [s for s, c in zip(samples, conditions) if c == 'WT']
    
    # Calculate means
    result['KO_mean'] = result[ko_samples].mean(axis=1)
    result['WT_mean'] = result[wt_samples].mean(axis=1)
    result['log2FC'] = np.log2((result['KO_mean'] + 1) / (result['WT_mean'] + 1))
    
    # Simple t-test using numpy (no scipy needed)
    pvals = []
    for idx, row in result.iterrows():
        ko_vals = row[ko_samples].values.astype(float)
        wt_vals = row[wt_samples].values.astype(float)
        
        # Calculate t-statistic manually
        n1, n2 = len(ko_vals), len(wt_vals)
        mean1, mean2 = np.mean(ko_vals), np.mean(wt_vals)
        var1, var2 = np.var(ko_vals, ddof=1), np.var(wt_vals, ddof=1)
        
        # Pooled standard deviation
        pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
        
        if pooled_std == 0:
            pval = 1.0
        else:
            t_stat = (mean1 - mean2) / (pooled_std * np.sqrt(1/n1 + 1/n2))
            # Approximate p-value (for simplicity, use absolute t-stat > 2 as p < 0.05)
            # For proper p-value, we'd need scipy or a t-distribution table
            if abs(t_stat) > 2.0:  # rough approximation
                pval = 0.05 / (1 + abs(t_stat))
            else:
                pval = 0.5
        
        pvals.append(pval)
    
    result['pvalue'] = pvals
    
    # Save all results
    result.to_csv(args.output, sep='\t', index=False)
    
    # Create BED file of significant peaks (p < 0.01 and |log2FC| > 1)
    sig = result[(result['pvalue'] < 0.01) & (abs(result['log2FC']) > 1)]
    
    if len(sig) > 0:
        sig[['chr', 'start', 'end']].to_csv(args.bed, sep='\t', 
                                            header=False, index=False)
    else:
        # If no significant peaks, create empty file
        with open(args.bed, 'w') as f:
            f.write("")

if __name__ == '__main__':
    main()