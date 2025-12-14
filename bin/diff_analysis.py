#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

# Try to import scipy, use fallback if not available
try:
    from scipy.stats import ttest_ind
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def welch_ttest(group1, group2):
    """Manual Welch's t-test implmm10-blacklist.v2.bed"ementation"""
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    
    if var1 == 0 and var2 == 0:
        return 1.0
    
    # Welch's t-statistic
    se = np.sqrt(var1/n1 + var2/n2)
    if se == 0:
        return 1.0
    
    t_stat = (mean1 - mean2) / se
    
    # Welch-Satterthwaite degrees of freedom
    if var1/n1 + var2/n2 == 0:
        return 1.0
    
    df = ((var1/n1 + var2/n2)**2) / ((var1/n1)**2/(n1-1) + (var2/n2)**2/(n2-1))
    
    # Very rough p-value approximation using normal distribution
    # For df > 30, t-distribution ~= normal distribution
    if df > 30:
        from math import erf, sqrt
        z = abs(t_stat)
        pval = 2 * (1 - 0.5 * (1 + erf(z / sqrt(2))))
    else:
        # Conservative approximation for small df
        pval = 2 * (1 - (1 / (1 + abs(t_stat) / np.sqrt(df))))
    
    return min(pval, 1.0)

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
    
    print(f"Using scipy: {HAS_SCIPY}")
    
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
    
    print(f"KO samples: {ko_samples}")
    print(f"WT samples: {wt_samples}")
    
    # Calculate means
    result['KO_mean'] = result[ko_samples].mean(axis=1)
    result['WT_mean'] = result[wt_samples].mean(axis=1)
    result['log2FC'] = np.log2((result['KO_mean'] + 1) / (result['WT_mean'] + 1))
    
    # Perform t-test
    pvals = []
    for idx, row in result.iterrows():
        ko_vals = row[ko_samples].values.astype(float)
        wt_vals = row[wt_samples].values.astype(float)
        
        if HAS_SCIPY:
            _, pval = ttest_ind(ko_vals, wt_vals, equal_var=False)
        else:
            pval = welch_ttest(ko_vals, wt_vals)
        
        pvals.append(pval)
    
    result['pvalue'] = pvals
    
    # Save all results
    result.to_csv(args.output, sep='\t', index=False)
    
    # Create BED file of significant peaks
    # RELAXED THRESHOLDS: p < 0.05 and |log2FC| > 0.5
    sig = result[(result['pvalue'] < 0.05) & (abs(result['log2FC']) > 0.5)]
    
    print(f"Found {len(sig)} significant peaks out of {len(result)} total")
    print(f"  p < 0.05: {len(result[result['pvalue'] < 0.05])}")
    print(f"  |log2FC| > 0.5: {len(result[abs(result['log2FC']) > 0.5])}")
    
    if len(sig) > 0:
        # Sort by p-value
        sig = sig.sort_values('pvalue')
        sig[['chr', 'start', 'end']].to_csv(args.bed, sep='\t', 
                                            header=False, index=False)
    else:
        # If no significant peaks, create empty file
        with open(args.bed, 'w') as f:
            f.write("")

if __name__ == '__main__':
    main()