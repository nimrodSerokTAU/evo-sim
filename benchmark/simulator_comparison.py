#!/usr/bin/env python3
"""
Simulator Comparison Script
Compares the indel simulator (list/tree) with external simulators like AliSim and PhastSim.
Uses the same setup as switch_factor.py for consistency.
"""

import subprocess
import psutil
import os
import tempfile
from pathlib import Path
import sys
import time
from typing import Dict, List, Tuple, Optional

import pandas as pd
from matplotlib import pyplot as plt


def run_indel_simulator(tree_file: Path, indel_rate: Tuple[float, float], 
                       sim_type: str, root_seq_length: int = 10000, 
                       temp_dir: Optional[Path] = None) -> Dict:
    """Run our indel simulator using the CLI tool and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    # Build command using the CLI tool
    cmd = [
        "indel-simulator",
        "--type", sim_type,
        "--insertion_rate", str(indel_rate[0]),
        "--deletion_rate", str(indel_rate[1]),
        "--tree_file", str(tree_file),
        "--original_sequence_length", str(root_seq_length),
        "--output_type", "drop_output",  # Don't write files
        "--seed", "42",
        "--benchmark"
    ]
    
    # Measure memory and time
    start_time = time.time()
    process = psutil.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Monitor memory usage
    max_memory_mb = 0
    try:
        while process.is_running():
            try:
                mem_info = process.memory_info()
                current_memory_mb = mem_info.rss / 1024 / 1024  # Convert to MB
                max_memory_mb = max(max_memory_mb, current_memory_mb)
                time.sleep(0.01)  # Check every 10ms
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except psutil.NoSuchProcess:
        pass
    
    stdout, stderr = process.communicate()
    end_time = time.time()
    runtime_seconds = end_time - start_time
    
    # Parse benchmark results from stdout if available
    benchmark_time = None
    if stdout:
        stdout_text = stdout.decode()
        # Look for "Average runtime" in the benchmark output
        for line in stdout_text.split('\n'):
            if "Average runtime:" in line:
                try:
                    benchmark_time = float(line.split(':')[1].strip().replace('s', ''))
                except (ValueError, IndexError):
                    pass
    
    return {
        'runtime_seconds': benchmark_time if benchmark_time else runtime_seconds,
        'max_memory_mb': max_memory_mb,
        'returncode': process.returncode,
        'stdout': stdout.decode() if stdout else "",
        'stderr': stderr.decode() if stderr else ""
    }


def run_alisim(tree_file: Path, model_params: Dict, root_seq_length: int = 10000,
               temp_dir: Optional[Path] = None) -> Dict:
    """Run AliSim simulator and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    output_prefix = temp_dir / "alisim_output"
    
    # Create a model that prevents substitutions but allows indels
    ins_rate = model_params.get('insertion_rate', 0.01)
    del_rate = model_params.get('deletion_rate', 0.01)
    
    # Use invariable sites close to 1 to prevent substitutions
    # Format: model+I{proportion}+INDEL{ins_rate,del_rate}
    model_str = f"JC+I{{0.9999999999}}+INDEL{{{ins_rate},{del_rate}}}"
    
    # Build AliSim command
    cmd = [
        "alisim",
        "--root-seq", str(root_seq_length),
        "-t", str(tree_file),
        "-m", model_str,
        "--length", str(root_seq_length),
        "--output-format", "FASTA",
        str(output_prefix)
    ]
    
    start_time = time.time()
    process = psutil.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Monitor memory usage
    max_memory_mb = 0
    try:
        while process.is_running():
            try:
                mem_info = process.memory_info()
                current_memory_mb = mem_info.rss / 1024 / 1024
                max_memory_mb = max(max_memory_mb, current_memory_mb)
                time.sleep(0.01)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except psutil.NoSuchProcess:
        pass
    
    stdout, stderr = process.communicate()
    end_time = time.time()
    runtime_seconds = end_time - start_time
    
    # Clean up output files
    for output_file in temp_dir.glob("alisim_output*"):
        try:
            output_file.unlink()
        except FileNotFoundError:
            pass
    
    return {
        'runtime_seconds': runtime_seconds,
        'max_memory_mb': max_memory_mb,
        'returncode': process.returncode,
        'stdout': stdout.decode() if stdout else "",
        'stderr': stderr.decode() if stderr else ""
    }


def run_phastsim(tree_file: Path, model_params: Dict, root_seq_length: int = 10000,
                 temp_dir: Optional[Path] = None) -> Dict:
    """Run PhastSim simulator and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    output_file = temp_dir / "phastsim_output.fa"
    
    ins_rate = model_params.get('insertion_rate', 0.01)
    del_rate = model_params.get('deletion_rate', 0.01)
    
    # Build PhastSim command with no substitutions (adjust based on actual PhastSim parameters)
    cmd = [
        "phastsim",
        "--tree", str(tree_file),
        "--rootSeqLen", str(root_seq_length),
        "--insertionRate", str(ins_rate),
        "--deletionRate", str(del_rate),
        "--substitutionRate", "0.0",  # No substitutions
        "--output", str(output_file)
    ]
    
    start_time = time.time()
    process = psutil.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Monitor memory usage
    max_memory_mb = 0
    try:
        while process.is_running():
            try:
                mem_info = process.memory_info()
                current_memory_mb = mem_info.rss / 1024 / 1024
                max_memory_mb = max(max_memory_mb, current_memory_mb)
                time.sleep(0.01)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except psutil.NoSuchProcess:
        pass
    
    stdout, stderr = process.communicate()
    end_time = time.time()
    runtime_seconds = end_time - start_time
    
    # Clean up output file
    if output_file.exists():
        output_file.unlink()
    
    return {
        'runtime_seconds': runtime_seconds,
        'max_memory_mb': max_memory_mb,
        'returncode': process.returncode,
        'stdout': stdout.decode() if stdout else "",
        'stderr': stderr.decode() if stderr else ""
    }


def check_simulator_availability() -> Dict[str, bool]:
    """Check which external simulators are available."""
    simulators = {}
    
    # Check indel-simulator (our CLI tool)
    try:
        result = subprocess.run(['indel-simulator', '--help'], 
                              capture_output=True, timeout=10)
        simulators['indel-simulator'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['indel-simulator'] = False
    
    # Check AliSim
    try:
        result = subprocess.run(['alisim', '--version'], 
                              capture_output=True, timeout=10)
        simulators['alisim'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['alisim'] = False
    
    # Check PhastSim
    try:
        result = subprocess.run(['phastsim', '--help'], 
                              capture_output=True, timeout=10)
        simulators['phastsim'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['phastsim'] = False
    
    return simulators


def main():
    """Main comparison function following switch_factor.py pattern."""
    
    ROOT_SEQUENCE_LENGTH = 10000
    SCALED_TREES_PATH = Path.cwd() / "benchmark" / "scaled_trees"
    
    # Check available simulators
    available_simulators = check_simulator_availability()
    print("Available simulators:", available_simulators)
    
    if not available_simulators.get('indel-simulator', False):
        print("ERROR: indel-simulator CLI tool not found. Please install it first.")
        return None
    
    # Same rates as switch_factor.py
    rates = [i * 0.01 for i in range(1, 10)]
    rates = list(zip(rates, rates[::-1]))
    
    all_results = []
    
    # Create temporary directory for this run
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        
        for RATE_MULTIPLIER in [1, 5, 10]:  # Same as switch_factor.py
            print(f"RATE_MULTIPLIER: {RATE_MULTIPLIER}")
            
            # Use pre-scaled tree files
            scaled_tree_path = SCALED_TREES_PATH / f"scaled_{RATE_MULTIPLIER}.tree"
            
            if not scaled_tree_path.exists():
                print(f"WARNING: Scaled tree file not found: {scaled_tree_path}")
                print(f"Available files in {SCALED_TREES_PATH}:")
                if SCALED_TREES_PATH.exists():
                    for f in SCALED_TREES_PATH.glob("*.tree"):
                        print(f"  {f.name}")
                continue
            
            for indel_rate in rates:
                print(f"Testing rates: ins={indel_rate[0]:.2f}, del={indel_rate[1]:.2f}")
                
                result_row = {
                    'branch_scale': RATE_MULTIPLIER,
                    'insertion_rate': indel_rate[0],
                    'deletion_rate': indel_rate[1]
                }
                
                # Run our indel simulator (list version)
                try:
                    indel_list_result = run_indel_simulator(
                        scaled_tree_path, indel_rate, "list", 
                        ROOT_SEQUENCE_LENGTH, temp_dir
                    )
                    result_row.update({
                        'indel_list_time': indel_list_result['runtime_seconds'],
                        'indel_list_memory': indel_list_result['max_memory_mb'],
                        'indel_list_success': indel_list_result['returncode'] == 0
                    })
                    print(f"  Indel (list): {indel_list_result['runtime_seconds']:.3f}s, "
                          f"{indel_list_result['max_memory_mb']:.1f}MB")
                except Exception as e:
                    print(f"  Indel (list) failed: {e}")
                    result_row.update({
                        'indel_list_time': None,
                        'indel_list_memory': None,
                        'indel_list_success': False
                    })
                
                # Run our indel simulator (tree version)
                try:
                    indel_tree_result = run_indel_simulator(
                        scaled_tree_path, indel_rate, "tree", 
                        ROOT_SEQUENCE_LENGTH, temp_dir
                    )
                    result_row.update({
                        'indel_tree_time': indel_tree_result['runtime_seconds'],
                        'indel_tree_memory': indel_tree_result['max_memory_mb'],
                        'indel_tree_success': indel_tree_result['returncode'] == 0
                    })
                    print(f"  Indel (tree): {indel_tree_result['runtime_seconds']:.3f}s, "
                          f"{indel_tree_result['max_memory_mb']:.1f}MB")
                except Exception as e:
                    print(f"  Indel (tree) failed: {e}")
                    result_row.update({
                        'indel_tree_time': None,
                        'indel_tree_memory': None,
                        'indel_tree_success': False
                    })
                
                # Run AliSim if available
                if available_simulators.get('alisim', False):
                    try:
                        alisim_result = run_alisim(
                            scaled_tree_path, 
                            {'insertion_rate': indel_rate[0], 'deletion_rate': indel_rate[1]},
                            ROOT_SEQUENCE_LENGTH, temp_dir
                        )
                        result_row.update({
                            'alisim_time': alisim_result['runtime_seconds'],
                            'alisim_memory': alisim_result['max_memory_mb'],
                            'alisim_success': alisim_result['returncode'] == 0
                        })
                        print(f"  AliSim: {alisim_result['runtime_seconds']:.3f}s, "
                              f"{alisim_result['max_memory_mb']:.1f}MB")
                    except Exception as e:
                        print(f"  AliSim failed: {e}")
                        result_row.update({
                            'alisim_time': None,
                            'alisim_memory': None,
                            'alisim_success': False
                        })
                else:
                    result_row.update({
                        'alisim_time': None,
                        'alisim_memory': None,
                        'alisim_success': False
                    })
                
                # Run PhastSim if available
                if available_simulators.get('phastsim', False):
                    try:
                        phastsim_result = run_phastsim(
                            scaled_tree_path,
                            {'insertion_rate': indel_rate[0], 'deletion_rate': indel_rate[1]},
                            ROOT_SEQUENCE_LENGTH, temp_dir
                        )
                        result_row.update({
                            'phastsim_time': phastsim_result['runtime_seconds'],
                            'phastsim_memory': phastsim_result['max_memory_mb'],
                            'phastsim_success': phastsim_result['returncode'] == 0
                        })
                        print(f"  PhastSim: {phastsim_result['runtime_seconds']:.3f}s, "
                              f"{phastsim_result['max_memory_mb']:.1f}MB")
                    except Exception as e:
                        print(f"  PhastSim failed: {e}")
                        result_row.update({
                            'phastsim_time': None,
                            'phastsim_memory': None,
                            'phastsim_success': False
                        })
                else:
                    result_row.update({
                        'phastsim_time': None,
                        'phastsim_memory': None,
                        'phastsim_success': False
                    })
                
                all_results.append(result_row)
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(all_results)
    results_df.to_csv("simulator_comparison.csv", index=False)
    print(f"\nResults saved to simulator_comparison.csv")
    
    # Create comparison plots
    create_comparison_plots(results_df)
    
    return results_df


def create_comparison_plots(results_df: pd.DataFrame):
    """Create comparison plots similar to switch_factor.py."""
    
    # Time comparison plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, (branch_scale, group_data) in enumerate(results_df.groupby("branch_scale")):
        ax = axes[i]
        
        # Plot time comparisons
        if 'indel_list_time' in group_data.columns:
            ax.plot(group_data['insertion_rate'], group_data['indel_list_time'], 
                   'o-', label='Indel (list)', markersize=6)
        
        if 'indel_tree_time' in group_data.columns:
            ax.plot(group_data['insertion_rate'], group_data['indel_tree_time'], 
                   's-', label='Indel (tree)', markersize=6)
        
        if 'alisim_time' in group_data.columns and group_data['alisim_time'].notna().any():
            ax.plot(group_data['insertion_rate'], group_data['alisim_time'], 
                   '^-', label='AliSim', markersize=6)
        
        if 'phastsim_time' in group_data.columns and group_data['phastsim_time'].notna().any():
            ax.plot(group_data['insertion_rate'], group_data['phastsim_time'], 
                   'D-', label='PhastSim', markersize=6)
        
        ax.set_title(f"Runtime Comparison (Scale: {branch_scale})", size=14)
        ax.set_xlabel("Insertion Rate", size=12)
        ax.set_ylabel("Runtime (seconds)", size=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("assets/simulator_runtime_comparison.png", dpi=300, bbox_inches="tight")
    plt.savefig("assets/simulator_runtime_comparison.svg", bbox_inches="tight")
    plt.show()
    
    # Memory comparison plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, (branch_scale, group_data) in enumerate(results_df.groupby("branch_scale")):
        ax = axes[i]
        
        # Plot memory comparisons
        if 'indel_list_memory' in group_data.columns:
            ax.plot(group_data['insertion_rate'], group_data['indel_list_memory'], 
                   'o-', label='Indel (list)', markersize=6)
        
        if 'indel_tree_memory' in group_data.columns:
            ax.plot(group_data['insertion_rate'], group_data['indel_tree_memory'], 
                   's-', label='Indel (tree)', markersize=6)
        
        if 'alisim_memory' in group_data.columns and group_data['alisim_memory'].notna().any():
            ax.plot(group_data['insertion_rate'], group_data['alisim_memory'], 
                   '^-', label='AliSim', markersize=6)
        
        if 'phastsim_memory' in group_data.columns and group_data['phastsim_memory'].notna().any():
            ax.plot(group_data['insertion_rate'], group_data['phastsim_memory'], 
                   'D-', label='PhastSim', markersize=6)
        
        ax.set_title(f"Memory Usage Comparison (Scale: {branch_scale})", size=14)
        ax.set_xlabel("Insertion Rate", size=12)
        ax.set_ylabel("Max Memory (MB)", size=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("assets/simulator_memory_comparison.png", dpi=300, bbox_inches="tight")
    plt.savefig("assets/simulator_memory_comparison.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    results = main()
    print("\nComparison completed!")
    print(f"Shape of results: {results.shape}")
    print("\nSample results:")
    print(results.head())
