#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy import stats

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
    dfs = []
    for i, count_file in enumerate(args.counts):
        df = pd.read_csv(count_file, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'count'])
        df = df.rename(columns={'count': samples[i]})
        if i == 0:
            result = df[['chr', 'start', 'end', samples[i]]]
        else:
            result = result.merge(df[['chr', 'start', 'end', samples[i]]], 
                                on=['chr', 'start', 'end'])
    
    # Perform differential analysis (simple t-test)
    ko_samples = [s for s, c in zip(samples, conditions) if c == 'KO']
    wt_samples = [s for s, c in zip(samples, conditions) if c == 'WT']
    
    result['KO_mean'] = result[ko_samples].mean(axis=1)
    result['WT_mean'] = result[wt_samples].mean(axis=1)
    result['log2FC'] = np.log2((result['KO_mean'] + 1) / (result['WT_mean'] + 1))
    
    pvals = []
    for idx, row in result.iterrows():
        ko_vals = row[ko_samples].values
        wt_vals = row[wt_samples].values
        _, pval = stats.ttest_ind(ko_vals, wt_vals)
        pvals.append(pval)
    
    result['pvalue'] = pvals
    result.to_csv(args.output, sep='\t', index=False)
    
    # Create BED file of significant peaks
    sig = result[result['pvalue'] < 0.01]
    sig[['chr', 'start', 'end']].to_csv(args.bed, sep='\t', 
                                        header=False, index=False)

if __name__ == '__main__':
    main()