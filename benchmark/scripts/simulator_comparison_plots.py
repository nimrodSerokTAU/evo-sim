#!/usr/bin/env python3
"""
Plotting Script for Dual-Mode Simulator Comparison Results
Creates comparison plots from the CSV output of the dual-mode simulator_comparison.py
"""

from pathlib import Path
from typing import Optional

import pandas as pd
from matplotlib import pyplot as plt


def create_comparison_plots(results_df: pd.DataFrame, output_dir: Optional[Path] = None):
    """Create comparison plots from dual-mode simulator results DataFrame."""
    
    if output_dir is None:
        output_dir = Path("assets")
    
    output_dir.mkdir(exist_ok=True)
    
    # First, create a standalone legend
    create_standalone_legend(output_dir)
    
    # Time comparison plots - using your exact style
    for i, (branch_scale, group_data) in enumerate(results_df.groupby("branch_scale")):
        
        # Create figure for this branch scale
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Define styles for different simulators and modes
        styles = ['C0s-', 'C1D-', 'C2^-', 'C3v-']  # Different colors and markers
        
        plot_data = []
        labels = []
        
        # Collect data for our simulator (both modes)
        for mode_idx, mode in enumerate(['with_indels', 'without_indels']):
            mode_data = group_data[group_data['mode'] == mode]
            if 'our_sim_time' in mode_data.columns:
                mask = mode_data['our_sim_time'].notna()
                if mask.any():
                    plot_data.append({
                        'x': mode_data.loc[mask, 'insertion_rate'].values,
                        'y': mode_data.loc[mask, 'our_sim_time'].values,
                        'style': styles[mode_idx],
                        'label': f'Our Sim ({"with" if mode == "with_indels" else "no"} indels)'
                    })
        
        # Collect data for AliSim (both modes)
        for mode_idx, mode in enumerate(['with_indels', 'without_indels']):
            mode_data = group_data[group_data['mode'] == mode]
            if 'alisim_time' in mode_data.columns:
                mask = mode_data['alisim_time'].notna()
                if mask.any():
                    plot_data.append({
                        'x': mode_data.loc[mask, 'insertion_rate'].values,
                        'y': mode_data.loc[mask, 'alisim_time'].values,
                        'style': styles[mode_idx + 2],
                        'label': f'AliSim ({"with" if mode == "with_indels" else "no"} indels)'
                    })
        
        # Plot all data WITHOUT legend
        for data in plot_data:
            ax.plot(data['x'], data['y'], data['style'], markersize=9)
        
        # Apply your exact styling (no legend)
        ax.set_title(f"Multiplication factor: {branch_scale}", size=18)
        ax.set_xlabel("Insertion rate", size=16)
        ax.set_ylabel("Runtime (seconds)", size=16)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_yscale('log')  # Log scale for better visualization
        
        plt.tight_layout()
        plt.savefig(output_dir / f"simulator_runtime_comparison_factor_{branch_scale}.svg", 
                   bbox_inches="tight", dpi=300)
        plt.savefig(output_dir / f"simulator_runtime_comparison_factor_{branch_scale}.png", 
                   bbox_inches="tight", dpi=300)
        plt.close()
    
    # Memory comparison plots - same style, no legends
    for i, (branch_scale, group_data) in enumerate(results_df.groupby("branch_scale")):
        
        # Create figure for this branch scale
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Define styles for different simulators and modes
        styles = ['C0s-', 'C1D-', 'C2^-', 'C3v-']
        
        plot_data = []
        
        # Collect data for our simulator memory (both modes)
        for mode_idx, mode in enumerate(['with_indels', 'without_indels']):
            mode_data = group_data[group_data['mode'] == mode]
            if 'our_sim_memory' in mode_data.columns:
                mask = mode_data['our_sim_memory'].notna()
                if mask.any():
                    plot_data.append({
                        'x': mode_data.loc[mask, 'insertion_rate'].values,
                        'y': mode_data.loc[mask, 'our_sim_memory'].values,
                        'style': styles[mode_idx],
                        'label': f'Our Sim ({"with" if mode == "with_indels" else "no"} indels)'
                    })
        
        # Collect data for AliSim memory (both modes)
        for mode_idx, mode in enumerate(['with_indels', 'without_indels']):
            mode_data = group_data[group_data['mode'] == mode]
            if 'alisim_memory' in mode_data.columns:
                mask = mode_data['alisim_memory'].notna()
                if mask.any():
                    plot_data.append({
                        'x': mode_data.loc[mask, 'insertion_rate'].values,
                        'y': mode_data.loc[mask, 'alisim_memory'].values,
                        'style': styles[mode_idx + 2],
                        'label': f'AliSim ({"with" if mode == "with_indels" else "no"} indels)'
                    })
        
        # Plot all data WITHOUT legend
        for data in plot_data:
            ax.plot(data['x'], data['y'], data['style'], markersize=9)
        
        # Apply your exact styling (no legend)
        ax.set_title(f"Multiplication factor: {branch_scale}", size=18)
        ax.set_xlabel("Insertion rate", size=16)
        ax.set_ylabel("Max Memory (MB)", size=16)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_yscale('log')  # Log scale for better visualization

        plt.tight_layout()
        plt.savefig(output_dir / f"simulator_memory_comparison_factor_{branch_scale}.svg", 
                   bbox_inches="tight", dpi=300)
        plt.savefig(output_dir / f"simulator_memory_comparison_factor_{branch_scale}.png", 
                   bbox_inches="tight", dpi=300)
        plt.close()


def create_standalone_legend(output_dir: Path):
    """Create a standalone legend plot."""
    
    # Create a figure just for the legend with minimal size
    fig, ax = plt.subplots(figsize=(8, 1.5))
    
    # Define the same styles and labels
    styles = ['C0s-', 'C1D-', 'C2^-', 'C3v-']
    labels = ['Block tree (with indels)', 'Block tree (no indels)', 'AliSim (with indels)', 'AliSim (no indels)']
    
    # Create invisible plots just to generate legend entries
    for style, label in zip(styles, labels):
        ax.plot([], [], style, markersize=9, label=label)
    
    # Hide the axes completely
    ax.axis('off')
    
    # Create the legend with tight spacing
    legend = ax.legend(loc='center', ncol=2, frameon=True, fontsize=14, 
                      columnspacing=1.0, handletextpad=0.5)
    
    # Set tight layout with minimal padding
    plt.tight_layout(pad=0.1)
    
    # Save with minimal bounding box
    plt.savefig(output_dir / "simulator_comparison_legend.svg", bbox_inches="tight", 
               dpi=300, pad_inches=0.05)
    plt.savefig(output_dir / "simulator_comparison_legend.png", bbox_inches="tight", 
               dpi=300, pad_inches=0.05)
    plt.close()


def create_summary_statistics(results_df: pd.DataFrame) -> pd.DataFrame:
    """Create summary statistics from the results."""
    
    summary_stats = []
    
    for (branch_scale, mode), group in results_df.groupby(['branch_scale', 'mode']):
        stats = {'branch_scale': branch_scale, 'mode': mode}
        
        # Calculate average times
        for sim_type in ['our_sim', 'alisim']:
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


def main():
    """Main function for the plotting script."""
    
    input_file = Path("benchmark/assets/data/simulator_comparison_dual_mode.csv")
    plots_output_dir = Path("benchmark/assets/plots/vs_alisim/indels")
    data_output_dir = Path("benchmark/assets/data")
    
    # Check if input file exists
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        print("Please run the dual-mode simulator_comparison.py script first.")
        return
    
    # Load and process data
    try:
        results_df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    print(f"Loaded {len(results_df)} comparison results")
    print("Creating comparison plots...")
    
    # Create plots output directory
    plots_output_dir.mkdir(parents=True, exist_ok=True)
    data_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create plots
    create_comparison_plots(results_df, plots_output_dir)
    
    # Create summary statistics
    summary_stats = create_summary_statistics(results_df)
    summary_stats.to_csv(data_output_dir / "summary_statistics_indels.csv", index=False)
    
    print(f"Plots saved to {plots_output_dir}")
    print(f"Summary statistics saved to {data_output_dir / 'summary_statistics_indels.csv'}")


if __name__ == "__main__":
    main()