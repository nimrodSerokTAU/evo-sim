#!/usr/bin/env python3
"""
Plotting Script for Simulator Comparison Results
Creates comparison plots from the CSV output of simulator_comparison.py
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
from matplotlib import pyplot as plt


def create_comparison_plots(results_df: pd.DataFrame, output_dir: Optional[Path] = None):
    """Create comparison plots from simulator results DataFrame."""
    
    if output_dir is None:
        output_dir = Path("assets")
    
    output_dir.mkdir(exist_ok=True)
    
    # Time comparison plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, (branch_scale, group_data) in enumerate(results_df.groupby("branch_scale")):
        ax = axes[i]
        
        # Plot time comparisons
        if 'indel_list_time' in group_data.columns:
            mask = group_data['indel_list_time'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'indel_list_time'], 
                       'o-', label='Indel (list)', markersize=6, linewidth=2)
        
        if 'indel_tree_time' in group_data.columns:
            mask = group_data['indel_tree_time'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'indel_tree_time'], 
                       's-', label='Indel (tree)', markersize=6, linewidth=2)
        
        if 'alisim_time' in group_data.columns:
            mask = group_data['alisim_time'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'alisim_time'], 
                       '^-', label='AliSim', markersize=6, linewidth=2)
        
        if 'phastsim_time' in group_data.columns:
            mask = group_data['phastsim_time'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'phastsim_time'], 
                       'D-', label='PhastSim', markersize=6, linewidth=2)
        
        ax.set_title(f"Runtime Comparison (Scale: {branch_scale})", size=14)
        ax.set_xlabel("Insertion Rate", size=12)
        ax.set_ylabel("Runtime (seconds)", size=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')  # Log scale for better visualization
    
    plt.tight_layout()
    plt.savefig(output_dir / "simulator_runtime_comparison.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "simulator_runtime_comparison.svg", bbox_inches="tight")
    plt.show()
    
    # Memory comparison plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, (branch_scale, group_data) in enumerate(results_df.groupby("branch_scale")):
        ax = axes[i]
        
        # Plot memory comparisons
        if 'indel_list_memory' in group_data.columns:
            mask = group_data['indel_list_memory'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'indel_list_memory'], 
                       'o-', label='Indel (list)', markersize=6, linewidth=2)
        
        if 'indel_tree_memory' in group_data.columns:
            mask = group_data['indel_tree_memory'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'indel_tree_memory'], 
                       's-', label='Indel (tree)', markersize=6, linewidth=2)
        
        if 'alisim_memory' in group_data.columns:
            mask = group_data['alisim_memory'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'alisim_memory'], 
                       '^-', label='AliSim', markersize=6, linewidth=2)
        
        if 'phastsim_memory' in group_data.columns:
            mask = group_data['phastsim_memory'].notna()
            if mask.any():
                ax.plot(group_data.loc[mask, 'insertion_rate'], 
                       group_data.loc[mask, 'phastsim_memory'], 
                       'D-', label='PhastSim', markersize=6, linewidth=2)
        
        ax.set_title(f"Memory Usage Comparison (Scale: {branch_scale})", size=14)
        ax.set_xlabel("Insertion Rate", size=12)
        ax.set_ylabel("Max Memory (MB)", size=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "simulator_memory_comparison.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "simulator_memory_comparison.svg", bbox_inches="tight")
    plt.show()


def create_summary_statistics(results_df: pd.DataFrame) -> pd.DataFrame:
    """Create summary statistics from the results."""
    
    # Group by branch scale and calculate means
    summary_stats = []
    
    for branch_scale, group in results_df.groupby('branch_scale'):
        stats = {'branch_scale': branch_scale}
        
        # Calculate average times
        for sim_type in ['indel_list', 'indel_tree', 'alisim', 'phastsim']:
            time_col = f'{sim_type}_time'
            memory_col = f'{sim_type}_memory'
            success_col = f'{sim_type}_success'
            
            if time_col in group.columns:
                valid_times = group[time_col].dropna()
                if len(valid_times) > 0:
                    stats[f'{sim_type}_avg_time'] = valid_times.mean()
                    stats[f'{sim_type}_std_time'] = valid_times.std()
                    stats[f'{sim_type}_min_time'] = valid_times.min()
                    stats[f'{sim_type}_max_time'] = valid_times.max()
            
            if memory_col in group.columns:
                valid_memory = group[memory_col].dropna()
                if len(valid_memory) > 0:
                    stats[f'{sim_type}_avg_memory'] = valid_memory.mean()
                    stats[f'{sim_type}_std_memory'] = valid_memory.std()
                    stats[f'{sim_type}_max_memory'] = valid_memory.max()
            
            if success_col in group.columns:
                stats[f'{sim_type}_success_rate'] = group[success_col].sum() / len(group)
        
        summary_stats.append(stats)
    
    return pd.DataFrame(summary_stats)


def create_speedup_analysis(results_df: pd.DataFrame) -> pd.DataFrame:
    """Create speedup analysis comparing different simulators."""
    
    speedup_analysis = []
    
    for (branch_scale, ins_rate), group in results_df.groupby(['branch_scale', 'insertion_rate']):
        if len(group) != 1:
            continue  # Skip if multiple rows for same condition
        
        row = group.iloc[0]
        analysis = {
            'branch_scale': branch_scale,
            'insertion_rate': ins_rate,
            'deletion_rate': row.get('deletion_rate', None)
        }
        
        # Use indel_tree as baseline for speedup calculations
        baseline_time = row.get('indel_tree_time', None)
        
        if baseline_time and pd.notna(baseline_time):
            for sim_type in ['indel_list', 'alisim', 'phastsim']:
                time_col = f'{sim_type}_time'
                sim_time = row.get(time_col, None)
                
                if sim_time and pd.notna(sim_time):
                    speedup = sim_time / baseline_time
                    analysis[f'{sim_type}_vs_tree_ratio'] = speedup
                    analysis[f'{sim_type}_faster_than_tree'] = speedup < 1.0
        
        speedup_analysis.append(analysis)
    
    return pd.DataFrame(speedup_analysis)


def main():
    """Main function for the plotting script."""
    parser = argparse.ArgumentParser(
        description="Create comparison plots from simulator results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot results from default CSV file
  python plot_comparison_results.py
  
  # Plot results from specific CSV file
  python plot_comparison_results.py --input my_results.csv
  
  # Save plots to specific directory
  python plot_comparison_results.py --output plots/
  
  # Generate summary statistics only
  python plot_comparison_results.py --stats-only
        """
    )
    
    parser.add_argument(
        "--input", "-i",
        type=str,
        default="simulator_comparison.csv",
        help="Input CSV file with comparison results (default: simulator_comparison.csv)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="assets",
        help="Output directory for plots (default: assets/)"
    )
    
    parser.add_argument(
        "--stats-only",
        action="store_true",
        help="Only generate summary statistics, don't create plots"
    )
    
    parser.add_argument(
        "--no-display",
        action="store_true",
        help="Don't display plots interactively (only save to files)"
    )
    
    args = parser.parse_args()
    
    # Load results
    input_file = Path(args.input)
    if not input_file.exists():
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    try:
        results_df = pd.read_csv(input_file)
        print(f"Loaded {len(results_df)} results from {input_file}")
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        sys.exit(1)
    
    # Set up output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    # Generate summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    summary_stats = create_summary_statistics(results_df)
    print(summary_stats.round(4))
    
    # Save summary statistics
    summary_file = output_dir / "summary_statistics.csv"
    summary_stats.to_csv(summary_file, index=False)
    print(f"\nSummary statistics saved to: {summary_file}")
    
    # Generate speedup analysis
    speedup_df = create_speedup_analysis(results_df)
    if not speedup_df.empty:
        print("\n" + "="*60)
        print("SPEEDUP ANALYSIS (relative to indel_tree)")
        print("="*60)
        print(speedup_df.round(4))
        
        speedup_file = output_dir / "speedup_analysis.csv"
        speedup_df.to_csv(speedup_file, index=False)
        print(f"\nSpeedup analysis saved to: {speedup_file}")
    
    # Create plots unless stats-only mode
    if not args.stats_only:
        if args.no_display:
            # Set matplotlib to non-interactive backend
            import matplotlib
            matplotlib.use('Agg')
        
        print(f"\nCreating comparison plots in: {output_dir}")
        create_comparison_plots(results_df, output_dir)
        print("Plots created successfully!")
    
    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
