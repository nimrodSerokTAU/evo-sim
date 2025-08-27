#!/usr/bin/env python3
"""
Simulator Comparison Script - Dual Mode Version
Compares the indel simulator with external simulators like AliSim in two modes:
1. With substitutions (default AliSim, msa-simulator endpoint)
2. Without substitutions (AliSim +I{0.9999999}, indel-simulator endpoint)
"""

import subprocess
import psutil
import tempfile
from pathlib import Path
import sys
import time
from typing import Dict, Tuple, Optional

import pandas as pd

IQTREE_EXE_NAME = "iqtree3_intel"


def process_monitor(cmd: list[str]):
    # Measure memory and time
    start_time = time.time()
    process = psutil.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Monitor memory usage
    max_memory_mb = 0
    try:
        while process.poll() is None:
            try:
                mem_info = process.memory_info()
                current_memory_mb = mem_info.rss / 1024 / 1024  # Convert to MB
                max_memory_mb = max(max_memory_mb, current_memory_mb)
                time.sleep(0.01)  # Check every 10ms
                elapsed_time = time.time() - start_time
                if elapsed_time > 300:
                    process.kill()
                    break
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except psutil.NoSuchProcess:
        pass

    try:
        mem_info = process.memory_info()
        current_memory_mb = mem_info.rss / 1024 / 1024
        max_memory_mb = max(max_memory_mb, current_memory_mb)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass        

    
    stdout, stderr = process.communicate()
    end_time = time.time()
    runtime_seconds = end_time - start_time

    return (
        runtime_seconds,
        max_memory_mb,
        process.returncode,
        stdout.decode() if stdout else "",
        stderr.decode() if stderr else ""
    )


def run_our_simulator(tree_file: Path, indel_rate: Tuple[float, float], 
                      sim_type: str, root_seq_length: int = 10000, 
                      temp_dir: Optional[Path] = None, with_substitutions: bool = True) -> Dict:
    """Run our simulator using the appropriate CLI tool and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    # Choose the appropriate CLI tool based on mode
    if with_substitutions:
        # Mode 1: With substitutions - use msa-simulator (combined indel + substitution)
        cli_tool = "msa-simulator"
    else:
        # Mode 2: Without substitutions - use indel-simulator (indel only)
        cli_tool = "indel-simulator"
    
    cmd = [
        cli_tool,
        "--type", sim_type,
        "--insertion_rate", str(indel_rate[0]),
        "--deletion_rate", str(indel_rate[1]),
        "--tree_file", str(tree_file),
        "--original_sequence_length", str(root_seq_length),
        # "--output_type", "drop_output",
        "--keep_in_memory",
        "--seed", "420"
    ]
    
    runtime_seconds, max_memory_mb, returncode, stdout, stderr = process_monitor(cmd=cmd)
    # print(stdout, stderr)
    # Parse benchmark results from stdout if available
    benchmark_time = None
    if stdout:
        stdout_text = stdout
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
        'returncode': returncode,
        'stdout': stdout,
        'stderr': stderr
    }


def run_alisim(tree_file: Path, model_params: Dict, root_seq_length: int = 10000,
               temp_dir: Optional[Path] = None, with_substitutions: bool = True) -> Dict:
    """Run AliSim simulator and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    output_prefix = temp_dir / "alisim_output"
    
    ins_rate = model_params.get('insertion_rate', 0.01)
    del_rate = model_params.get('deletion_rate', 0.01)
    
    if with_substitutions:
        # Mode 1: With substitutions (default AliSim behavior)
        model_str = "JTT"
    else:
        # Mode 2: Without substitutions - use invariant sites
        model_str = "JTT+I{0.9999999}"
    
    # Build AliSim command (part of iqtree)
    cmd = [
        IQTREE_EXE_NAME,
        "--alisim", str(output_prefix),
        "-t", str(tree_file),
        "-m", model_str,
        "--length", str(root_seq_length),
        "--indel", f"{ins_rate},{del_rate}",
        "--indel-size", "POW{2.0/50},POW{2.0/50}",
        "--seed", "420"
    ]
    
    runtime_seconds, max_memory_mb, returncode, stdout, stderr = process_monitor(cmd=cmd)

    # Clean up output files
    for output_file in temp_dir.glob("alisim_output*"):
        try:
            output_file.unlink()
        except FileNotFoundError:
            pass
    
    return {
        'runtime_seconds': runtime_seconds,
        'max_memory_mb': max_memory_mb,
        'returncode': returncode,
        'stdout': stdout,
        'stderr': stderr
    }


def check_simulator_availability() -> Dict[str, bool]:
    """Check which external simulators are available."""
    simulators = {}
    
    # Check msa-simulator (our combined CLI tool)
    try:
        result = subprocess.run(['msa-simulator', '--help'], 
                              capture_output=True, timeout=10)
        simulators['msa-simulator'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['msa-simulator'] = False
    
    # Check indel-simulator (our indel-only CLI tool)
    try:
        result = subprocess.run(['indel-simulator', '--help'], 
                              capture_output=True, timeout=10)
        simulators['indel-simulator'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['indel-simulator'] = False
    
    # Check AliSim (part of IQ-TREE)
    try:
        result = subprocess.run([IQTREE_EXE_NAME, '--help'], 
                              capture_output=True, timeout=10)
        simulators['alisim'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['alisim'] = False
    
    return simulators


def main():
    """Main comparison function with dual-mode support."""
    
    ROOT_SEQUENCE_LENGTH = 10000
    SCALED_TREES_PATH = Path.cwd() / "benchmark" / "scaled_trees"
    
    # Check available simulators
    available_simulators = check_simulator_availability()
    print("Available simulators:", available_simulators)
    
    if not available_simulators.get('msa-simulator', False):
        print("WARNING: msa-simulator CLI tool not found. Mode 1 (with substitutions) will be skipped.")
    
    if not available_simulators.get('indel-simulator', False):
        print("WARNING: indel-simulator CLI tool not found. Mode 2 (without substitutions) will be skipped.")
    
    if not any([available_simulators.get('msa-simulator', False), 
                available_simulators.get('indel-simulator', False)]):
        print("ERROR: Neither msa-simulator nor indel-simulator CLI tools found. Please install them first.")
        return None
    
    # Same rates as switch_factor.py
    rates = [i * 0.01 for i in range(1, 10)]
    rates = list(zip(rates, rates[::-1]))

    all_results = []
    
    # Create temporary directory for this run
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        
        for RATE_MULTIPLIER in [1, 5, 10]:  # Same as switch_factor.py
            print(f"\nRATE_MULTIPLIER: {RATE_MULTIPLIER}")
            
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
                
                # Run both modes for each rate combination
                for mode_idx, (with_substitutions, mode_name) in enumerate([
                    (True, "with_substitutions"),
                    (False, "without_substitutions")
                ]):
                    print(f"  Mode: {mode_name}")
                    
                    result_row = {
                        'branch_scale': RATE_MULTIPLIER,
                        'insertion_rate': indel_rate[0],
                        'deletion_rate': indel_rate[1],
                        'mode': mode_name,
                        'with_substitutions': with_substitutions
                    }
                    
                    # Run our simulator (tree version) for the current mode
                    simulator_available = (
                        available_simulators.get('msa-simulator', False) if with_substitutions 
                        else available_simulators.get('indel-simulator', False)
                    )
                    
                    if simulator_available:
                        try:
                            our_sim_result = run_our_simulator(
                                scaled_tree_path, indel_rate, "tree", 
                                ROOT_SEQUENCE_LENGTH, temp_dir, with_substitutions
                            )
                            result_row.update({
                                'our_sim_time': our_sim_result['runtime_seconds'],
                                'our_sim_memory': our_sim_result['max_memory_mb'],
                                'our_sim_success': our_sim_result['returncode'] == 0
                            })
                            print(f"    Our simulator: {our_sim_result['runtime_seconds']:.3f}s, "
                                  f"{our_sim_result['max_memory_mb']:.1f}MB")
                        except Exception as e:
                            print(f"    Our simulator failed: {e}")
                            result_row.update({
                                'our_sim_time': None,
                                'our_sim_memory': None,
                                'our_sim_success': False
                            })
                    else:
                        result_row.update({
                            'our_sim_time': None,
                            'our_sim_memory': None,
                            'our_sim_success': False
                        })
                    
                    # Run AliSim if available
                    if available_simulators.get('alisim', False):
                        try:
                            alisim_result = run_alisim(
                                scaled_tree_path, 
                                {'insertion_rate': indel_rate[0], 'deletion_rate': indel_rate[1]},
                                ROOT_SEQUENCE_LENGTH, temp_dir, with_substitutions
                            )
                            result_row.update({
                                'alisim_time': alisim_result['runtime_seconds'],
                                'alisim_memory': alisim_result['max_memory_mb'],
                                'alisim_success': alisim_result['returncode'] == 0
                            })
                            print(f"    AliSim: {alisim_result['runtime_seconds']:.3f}s, "
                                  f"{alisim_result['max_memory_mb']:.1f}MB")
                        except Exception as e:
                            print(f"    AliSim failed: {e}")
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
                    
                    all_results.append(result_row)
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(all_results)
    results_df.to_csv("benchmark/assets/data/simulator_comparison_subs_dual_mode.csv", index=False)
    print(f"\nResults saved to simulator_comparison_subs_dual_mode.csv")
    
    # Print summary
    print(f"\nComparison completed!")
    print(f"Shape of results: {results_df.shape}")
    print(f"\nModes tested:")
    for mode in results_df['mode'].unique():
        count = len(results_df[results_df['mode'] == mode])
        print(f"  {mode}: {count} experiments")
    
    print("\nSample results:")
    print(results_df.head())
    
    return results_df


if __name__ == "__main__":
    results = main()
    if results is not None:
        print("\nDual-mode comparison completed successfully!")
    else:
        print("\nComparison failed!")
        sys.exit(1)