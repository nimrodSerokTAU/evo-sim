#!/usr/bin/env python3
"""
Plot Algorithm Comparison Results

Creates plots from indel algorithm comparison data.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
INPUT_FILE = "./algorithm_comparison/algorithm_comparison.csv"
OUTPUT_DIR = "./algorithm_comparison/plots"

def create_plots(df):
    """Create comparison plots"""
    
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Plot 1: Indel time and ratio comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Indel timing by algorithm
    ax1 = axes[0]
    df['total_indel_rate'] = df['insertion_rate'] + df['deletion_rate']
    
    for alg in ['naive', 'list', 'tree']:
        alg_data = df[df['algorithm'] == alg]
        if len(alg_data) > 0:
            ax1.scatter(alg_data['total_indel_rate'], alg_data['indel_time'], 
                       label=alg, alpha=0.7, s=50)
    
    ax1.set_xlabel('Total Indel Rate')
    ax1.set_ylabel('Indel Time (seconds)')
    ax1.set_title('Indel Simulation Time by Algorithm')
    ax1.legend()
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    
    # Indel ratio comparison
    ax2 = axes[1]
    algorithms = df['algorithm'].unique()
    ratios = [df[df['algorithm'] == alg]['indel_ratio'].values for alg in algorithms]
    
    bp = ax2.boxplot(ratios, labels=algorithms, patch_artist=True)
    colors = ['lightblue', 'lightgreen', 'lightcoral']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    ax2.axhline(0.5, color='red', linestyle='--', alpha=0.7, label='50% threshold')
    ax2.set_ylabel('Indel Time Fraction')
    ax2.set_title('Rate Limiting Analysis by Algorithm')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "algorithm_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Speedup analysis
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    pivot = df.pivot_table(values='speedup_vs_naive', 
                          index='sequence_length', 
                          columns=['algorithm', 'total_indel_rate'], 
                          aggfunc='mean')
    
    if not pivot.empty:
        sns.heatmap(pivot, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax)
        ax.set_title('Speedup vs Naive Algorithm')
        ax.set_xlabel('Algorithm & Indel Rate')
        ax.set_ylabel('Sequence Length')
    
    plt.tight_layout()
    plt.savefig(output_dir / "speedup_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_dir}")

def main():
    if not Path(INPUT_FILE).exists():
        print(f"Input file not found: {INPUT_FILE}")
        print("Run indel_algorithm_comparison.py first")
        return
    
    df = pd.read_csv(INPUT_FILE)
    create_plots(df)

if __name__ == "__main__":
    main()