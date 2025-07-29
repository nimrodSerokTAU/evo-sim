#!/usr/bin/env python3
"""
Plot Tree Comparison Results

Creates plots from tree comparison benchmark data.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
INPUT_FILE = "./tree_comparison/tree_comparison.csv"
OUTPUT_DIR = "./tree_comparison/plots"

def create_plots(df):
    """Create tree comparison plots"""
    
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Plot 1: Indel dominance across trees and algorithms
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Tree Comparison: Indel vs Substitution Timing Analysis', fontsize=16, fontweight='bold')
    
    # Subplot 1: Indel ratio by algorithm across all trees
    ax1 = axes[0, 0]
    algorithms = df['algorithm'].unique()
    ratios = [df[df['algorithm'] == alg]['indel_ratio'].values for alg in algorithms]
    
    bp = ax1.boxplot(ratios, labels=algorithms, patch_artist=True)
    colors = ['lightblue', 'lightgreen', 'lightcoral']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    ax1.axhline(0.5, color='red', linestyle='--', alpha=0.7, label='50% threshold')
    ax1.set_ylabel('Indel Time Fraction')
    ax1.set_title('Rate Limiting Analysis by Algorithm\n(Across All Trees)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Subplot 2: Indel ratio vs number of taxa
    ax2 = axes[0, 1]
    for alg in algorithms:
        alg_data = df[df['algorithm'] == alg]
        if len(alg_data) > 0:
            ax2.scatter(alg_data['num_taxa'], alg_data['indel_ratio'], 
                       label=alg, alpha=0.7, s=50)
    
    ax2.axhline(0.5, color='red', linestyle='--', alpha=0.7)
    ax2.set_xlabel('Number of Taxa')
    ax2.set_ylabel('Indel Time Fraction')
    ax2.set_title('Indel Dominance vs Tree Size')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Speedup analysis
    ax3 = axes[1, 0]
    speedup_data = df[df['algorithm'] != 'naive']  # Exclude naive from speedup plot
    
    if len(speedup_data) > 0:
        for alg in ['list', 'tree']:
            alg_data = speedup_data[speedup_data['algorithm'] == alg]
            if len(alg_data) > 0:
                ax3.scatter(alg_data['num_taxa'], alg_data['speedup_vs_naive'], 
                           label=f'{alg} vs naive', alpha=0.7, s=50)
    
    ax3.set_xlabel('Number of Taxa')
    ax3.set_ylabel('Speedup vs Naive Algorithm')
    ax3.set_title('Algorithm Speedup vs Tree Size')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    # Subplot 4: Tree-specific performance
    ax4 = axes[1, 1]
    
    # Show top trees by indel dominance
    tree_means = df.groupby('tree_file')['indel_ratio'].mean().sort_values(ascending=False)
    top_trees = tree_means.head(10)  # Top 10 trees
    
    if len(top_trees) > 0:
        y_pos = range(len(top_trees))
        bars = ax4.barh(y_pos, top_trees.values, alpha=0.7)
        
        # Color bars based on dominance level
        for i, bar in enumerate(bars):
            if top_trees.values[i] > 0.7:
                bar.set_color('red')
            elif top_trees.values[i] > 0.5:
                bar.set_color('orange')
            else:
                bar.set_color('lightblue')
        
        ax4.set_yticks(y_pos)
        ax4.set_yticklabels([name.split('.')[0][:15] for name in top_trees.index], fontsize=8)
        ax4.set_xlabel('Mean Indel Time Fraction')
        ax4.set_title('Top Trees by Indel Dominance')
        ax4.axvline(0.5, color='red', linestyle='--', alpha=0.7)
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "tree_comparison_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Detailed timing comparison
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle('Detailed Timing Analysis Across Trees', fontsize=16, fontweight='bold')
    
    # Absolute timing
    ax1 = axes[0]
    ax1.scatter(df['indel_time'], df['substitution_time'], 
               c=df['num_taxa'], alpha=0.6, s=50, cmap='viridis')
    
    max_time = max(df['indel_time'].max(), df['substitution_time'].max())
    ax1.plot([0, max_time], [0, max_time], 'r--', alpha=0.7, label='Equal time')
    
    ax1.set_xlabel('Indel Time (seconds)')
    ax1.set_ylabel('Substitution Time (seconds)')
    ax1.set_title('Absolute Timing Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(ax1.collections[0], ax=ax1)
    cbar.set_label('Number of Taxa')
    
    # Performance by tree complexity
    ax2 = axes[1]
    
    # Bin trees by number of taxa
    df['taxa_bin'] = pd.cut(df['num_taxa'], bins=5, labels=['Very Small', 'Small', 'Medium', 'Large', 'Very Large'])
    
    bin_data = []
    bin_labels = []
    
    for bin_label in df['taxa_bin'].cat.categories:
        bin_subset = df[df['taxa_bin'] == bin_label]
        if len(bin_subset) > 0:
            bin_data.append(bin_subset['indel_ratio'].values)
            bin_labels.append(f'{bin_label}\n(n={len(bin_subset)})')
    
    if bin_data:
        bp = ax2.boxplot(bin_data, labels=bin_labels, patch_artist=True)
        
        colors = plt.cm.Set3(range(len(bin_data)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
    
    ax2.axhline(0.5, color='red', linestyle='--', alpha=0.7)
    ax2.set_ylabel('Indel Time Fraction')
    ax2.set_title('Performance by Tree Complexity')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "tree_timing_details.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_dir}")

def main():
    if not Path(INPUT_FILE).exists():
        print(f"Input file not found: {INPUT_FILE}")
        print("Run tree_comparison_benchmark.py first")
        return
    
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded data for {len(df)} runs across {df['tree_file'].nunique()} trees")
    create_plots(df)

if __name__ == "__main__":
    main()