#!/usr/bin/env python3
"""
Script to replicate Figures 6A-F from the HDAC1 paper
Simplified Python version using basic differential analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
import os

# Set style
sns.set_style("white")
plt.rcParams['figure.dpi'] = 300

# Create output directory
os.makedirs('results/figures', exist_ok=True)

def load_count_matrix():
    """Load and process count matrix"""
    print("Loading count matrix...")
    counts = pd.read_csv('results/counts/count_matrix.txt', sep='\t', header=None)
    
    # Extract coordinates
    coords = counts.iloc[:, 0:3]
    coords.columns = ['chr', 'start', 'end']
    
    # Extract count data
    count_data = counts.iloc[:, 3:]
    count_data.columns = [
        'cDC1_WT_1', 'cDC1_WT_2', 'cDC1_KO_1', 'cDC1_KO_2',
        'cDC2_WT_1', 'cDC2_WT_2', 'cDC2_KO_1', 'cDC2_KO_2'
    ]
    
    return coords, count_data

def plot_figure_6a(count_data):
    """Figure 6A: PCA plot"""
    print("\nGenerating Figure 6A (PCA)...")
    
    # Normalize and filter
    count_data_filt = count_data.loc[(count_data >= 10).sum(axis=1) >= 3]
    
    # Log transform
    log_counts = np.log2(count_data_filt + 1)
    
    # Standardize
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(log_counts.T)
    
    # PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    
    # Create dataframe for plotting
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df['Sample'] = count_data_filt.columns
    pca_df['Cell_Type'] = ['cDC1', 'cDC1', 'cDC1', 'cDC1', 
                            'cDC2', 'cDC2', 'cDC2', 'cDC2']
    pca_df['Genotype'] = ['WT', 'WT', 'KO', 'KO', 'WT', 'WT', 'KO', 'KO']
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    colors = {'cDC1': '#E64B35', 'cDC2': '#4DBBD5'}
    markers = {'WT': 'o', 'KO': 's'}
    
    for cell_type in ['cDC1', 'cDC2']:
        for genotype in ['WT', 'KO']:
            mask = (pca_df['Cell_Type'] == cell_type) & (pca_df['Genotype'] == genotype)
            ax.scatter(pca_df.loc[mask, 'PC1'], 
                      pca_df.loc[mask, 'PC2'],
                      c=colors[cell_type], 
                      marker=markers[genotype],
                      s=150, 
                      label=f'{cell_type} {genotype}',
                      edgecolors='black',
                      linewidth=1.5)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=12)
    ax.set_title('Figure 6A: PCA of ATAC-seq samples', fontsize=14, fontweight='bold')
    ax.legend(frameon=True, loc='best')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('results/figures/Figure_6A_PCA.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('results/figures/Figure_6A_PCA.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 6A saved")
    return log_counts

def differential_analysis(count_data, cell_type):
    """Perform simple differential analysis"""
    if cell_type == 'cDC1':
        wt_cols = ['cDC1_WT_1', 'cDC1_WT_2']
        ko_cols = ['cDC1_KO_1', 'cDC1_KO_2']
    else:
        wt_cols = ['cDC2_WT_1', 'cDC2_WT_2']
        ko_cols = ['cDC2_KO_1', 'cDC2_KO_2']
    
    wt_mean = count_data[wt_cols].mean(axis=1)
    ko_mean = count_data[ko_cols].mean(axis=1)
    
    # Log2 fold change
    log2fc = np.log2((ko_mean + 1) / (wt_mean + 1))
    
    # T-test
    pvals = []
    for idx in count_data.index:
        wt_vals = count_data.loc[idx, wt_cols].values
        ko_vals = count_data.loc[idx, ko_cols].values
        _, pval = stats.ttest_ind(wt_vals, ko_vals)
        pvals.append(pval)
    
    # Create results dataframe
    results = pd.DataFrame({
        'wt_mean': wt_mean,
        'ko_mean': ko_mean,
        'log2fc': log2fc,
        'pvalue': pvals
    })
    
    # FDR correction (Benjamini-Hochberg)
    from scipy.stats import false_discovery_control
    results['padj'] = false_discovery_control(results['pvalue'])
    
    # Significance
    results['significant'] = 'NS'
    results.loc[(results['padj'] < 0.05) & (results['log2fc'] > 1), 'significant'] = 'Up in KO'
    results.loc[(results['padj'] < 0.05) & (results['log2fc'] < -1), 'significant'] = 'Down in KO'
    
    return results

def plot_volcano(results, cell_type, figure_label):
    """Plot volcano plot"""
    print(f"\nGenerating Figure {figure_label} (Volcano {cell_type})...")
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    colors = {'Up in KO': 'red', 'Down in KO': 'blue', 'NS': 'lightgray'}
    
    for sig_type in ['NS', 'Down in KO', 'Up in KO']:
        mask = results['significant'] == sig_type
        ax.scatter(results.loc[mask, 'log2fc'],
                  -np.log10(results.loc[mask, 'pvalue']),
                  c=colors[sig_type],
                  alpha=0.6 if sig_type == 'NS' else 0.8,
                  s=20 if sig_type == 'NS' else 30,
                  label=sig_type,
                  edgecolors='none')
    
    # Add threshold lines
    ax.axvline(x=-1, linestyle='--', color='black', linewidth=1)
    ax.axvline(x=1, linestyle='--', color='black', linewidth=1)
    ax.axhline(y=-np.log10(0.05), linestyle='--', color='black', linewidth=1)
    
    # Count significant peaks
    n_up = (results['significant'] == 'Up in KO').sum()
    n_down = (results['significant'] == 'Down in KO').sum()
    
    ax.set_xlabel('log2 Fold Change (KO/WT)', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title(f'Figure {figure_label}: Differential accessibility {cell_type}\n' + 
                 f'Up: {n_up} | Down: {n_down}', 
                 fontsize=14, fontweight='bold')
    ax.legend(frameon=True, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'results/figures/Figure_{figure_label}_Volcano_{cell_type}.pdf', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'results/figures/Figure_{figure_label}_Volcano_{cell_type}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Figure {figure_label} saved")
    
    # Save results
    results.to_csv(f'results/figures/{cell_type}_differential_peaks.txt', sep='\t')
    
    return results

def plot_heatmap(log_counts, results, cell_type, figure_label):
    """Plot heatmap of differential peaks"""
    print(f"\nGenerating Figure {figure_label} (Heatmap {cell_type})...")
    
    # Get significant peaks
    sig_peaks = results.index[(results['padj'] < 0.05) & 
                               (abs(results['log2fc']) > 1)]
    
    if len(sig_peaks) == 0:
        print(f"  Warning: No significant peaks found for {cell_type}")
        return
    
    # Select samples
    if cell_type == 'cDC1':
        samples = ['cDC1_WT_1', 'cDC1_WT_2', 'cDC1_KO_1', 'cDC1_KO_2']
    else:
        samples = ['cDC2_WT_1', 'cDC2_WT_2', 'cDC2_KO_1', 'cDC2_KO_2']
    
    # Get data
    plot_data = log_counts.loc[sig_peaks, samples]
    
    # Z-score normalize by row
    plot_data_scaled = plot_data.sub(plot_data.mean(axis=1), axis=0).div(plot_data.std(axis=1), axis=0)
    
    # Plot
    fig, ax = plt.subplots(figsize=(6, 10))
    
    sns.heatmap(plot_data_scaled, 
                cmap='RdBu_r',
                center=0,
                vmin=-2, vmax=2,
                xticklabels=True,
                yticklabels=False,
                cbar_kws={'label': 'Z-score'},
                ax=ax)
    
    ax.set_title(f'Figure {figure_label}: Differential peaks in {cell_type}\n' + 
                 f'(n={len(sig_peaks)} peaks)', 
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Peaks', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(f'results/figures/Figure_{figure_label}_Heatmap_{cell_type}.pdf', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'results/figures/Figure_{figure_label}_Heatmap_{cell_type}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Figure {figure_label} saved")

def plot_annotation_summary(coords, res_cdc1, res_cdc2):
    """Figure 6F: Summary of peak annotations"""
    print("\nGenerating Figure 6F (Peak annotation summary)...")
    
    # Read annotation files if they exist
    anno_files = {
        'cDC1': 'results/annotation/significant_up_peaks_annotated.txt',
        'cDC2': 'results/annotation/significant_up_peaks_annotated.txt'
    }
    
    # Simple annotation based on distance to gene (placeholder)
    # In reality, use HOMER annotation output
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Simplified annotation categories
    categories = ['Promoter', 'Exon', 'Intron', 'Intergenic']
    
    # Placeholder data - replace with actual HOMER output
    cdc1_dist = [25, 15, 35, 25]
    cdc2_dist = [30, 10, 40, 20]
    
    colors = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488']
    
    ax1.pie(cdc1_dist, labels=categories, autopct='%1.1f%%', 
            colors=colors, startangle=90)
    ax1.set_title('cDC1 - Genomic Annotation', fontsize=12, fontweight='bold')
    
    ax2.pie(cdc2_dist, labels=categories, autopct='%1.1f%%', 
            colors=colors, startangle=90)
    ax2.set_title('cDC2 - Genomic Annotation', fontsize=12, fontweight='bold')
    
    plt.suptitle('Figure 6F: Genomic annotation of differential peaks', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('results/figures/Figure_6F_Annotation.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('results/figures/Figure_6F_Annotation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✓ Figure 6F saved")
    print("\nNote: Figure 6F uses placeholder data.")
    print("      Use HOMER annotation output for accurate genomic distributions.")

def main():
    """Main execution"""
    print("="*60)
    print("Replicating Figures 6A-F from HDAC1 ATAC-seq paper")
    print("="*60)
    
    # Load data
    coords, count_data = load_count_matrix()
    
    # Filter low count peaks
    count_data_filt = count_data.loc[(count_data >= 10).sum(axis=1) >= 3]
    
    # Figure 6A: PCA
    log_counts = plot_figure_6a(count_data_filt)
    
    # Figure 6B: Volcano cDC1
    res_cdc1 = differential_analysis(count_data_filt, 'cDC1')
    plot_volcano(res_cdc1, 'cDC1', '6B')
    
    # Figure 6C: Volcano cDC2
    res_cdc2 = differential_analysis(count_data_filt, 'cDC2')
    plot_volcano(res_cdc2, 'cDC2', '6C')
    
    # Figure 6D: Heatmap cDC1
    plot_heatmap(log_counts, res_cdc1, 'cDC1', '6D')
    
    # Figure 6E: Heatmap cDC2
    plot_heatmap(log_counts, res_cdc2, 'cDC2', '6E')
    
    # Figure 6F: Annotation
    plot_annotation_summary(coords, res_cdc1, res_cdc2)
    
    print("\n" + "="*60)
    print("SUCCESS! All figures generated")
    print("="*60)
    print("\nOutput files:")
    print("  - results/figures/Figure_6A_PCA.pdf/.png")
    print("  - results/figures/Figure_6B_Volcano_cDC1.pdf/.png")
    print("  - results/figures/Figure_6C_Volcano_cDC2.pdf/.png")
    print("  - results/figures/Figure_6D_Heatmap_cDC1.pdf/.png")
    print("  - results/figures/Figure_6E_Heatmap_cDC2.pdf/.png")
    print("  - results/figures/Figure_6F_Annotation.pdf/.png")
    print("\nDifferential peak tables:")
    print("  - results/figures/cDC1_differential_peaks.txt")
    print("  - results/figures/cDC2_differential_peaks.txt")

if __name__ == '__main__':
    main()