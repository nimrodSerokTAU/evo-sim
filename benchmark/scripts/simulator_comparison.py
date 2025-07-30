#!/usr/bin/env python3
"""
Simulator Comparison Script
Compares the indel simulator (list/tree) with external simulators like AliSim.
Uses the same setup as switch_factor.py for consistency.
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

def run_msa_simulator(tree_file: Path, indel_rate: Tuple[float, float], 
                       sim_type: str, root_seq_length: int = 10000, 
                       temp_dir: Optional[Path] = None) -> Dict:
    """Run our indel simulator using the CLI tool and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    # Build command using the CLI tool
    cmd = [
        "msa-simulator",
        "--type", sim_type,
        "--insertion_rate", str(indel_rate[0]),
        "--deletion_rate", str(indel_rate[1]),
        "--tree_file", str(tree_file),
        "--original_sequence_length", str(root_seq_length),
        "--output_type", "drop_output",  # Don't write files
        "--seed", "42", "--verbose"
    ]
    
    runtime_seconds, max_memory_mb, returncode, stdout, stderr   = process_monitor(cmd=cmd)

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
               temp_dir: Optional[Path] = None) -> Dict:
    """Run AliSim simulator and measure performance."""
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp())
    
    output_prefix = temp_dir / "alisim_output"
    
    # AliSim is part of IQ-TREE. For indel-only simulation, we need to use a model
    # that prevents substitutions but allows indels. The model syntax is:
    # --alisim <output_prefix> -t <tree> -m <model> --length <seq_length>
    ins_rate = model_params.get('insertion_rate', 0.01)
    del_rate = model_params.get('deletion_rate', 0.01)
    # Use invariable sites to prevent substitutions and specify indel rates
    # Format: model
    model_str = "JTT"
    
    # Build AliSim command (part of iqtree)
    cmd = [
        IQTREE_EXE_NAME,
        "--alisim", str(output_prefix),
        "-t", str(tree_file),
        "-m", model_str,
        "--length", str(root_seq_length),
        "--indel", f"{ins_rate},{del_rate}",
        "--indel-size", "POW{2.0/50},POW{2.0/50}",
        "--seed", "42"
    ]
    
    runtime_seconds, max_memory_mb, returncode, stdout, stderr   = process_monitor(cmd=cmd)

    
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
    
    # Check msa-simulator (our CLI tool)
    try:
        result = subprocess.run(['msa-simulator', '--help'], 
                              capture_output=True, timeout=10)
        simulators['msa-simulator'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['msa-simulator'] = False
    
    # Check AliSim (part of IQ-TREE)
    try:
        result = subprocess.run([IQTREE_EXE_NAME, '--help'], 
                              capture_output=True, timeout=10)
        simulators['alisim'] = result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        simulators['alisim'] = False
    
    
    return simulators


def main():
    """Main comparison function following switch_factor.py pattern."""
    
    ROOT_SEQUENCE_LENGTH = 2000
    SCALED_TREES_PATH = Path.cwd() / "benchmark" / "scaled_trees"
    
    # Check available simulators
    available_simulators = check_simulator_availability()
    print("Available simulators:", available_simulators)
    
    if not available_simulators.get('msa-simulator', False):
        print("ERROR: msa-simulator CLI tool not found. Please install it first.")
        return None
    
    # Same rates as switch_factor.py
    rates = [i * 0.01 for i in range(1, 10)]
    rates = list(zip(rates, rates[::-1]))
    # rates = [(1.,.5)]

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
                
                # Run our indel simulator (tree version)
                try:
                    indel_tree_result = run_msa_simulator(
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
                
                # Run our indel simulator (list version)
                try:
                    indel_list_result = run_msa_simulator(
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
                

                

                                
                all_results.append(result_row)
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(all_results)
    results_df.to_csv("benchmark/assets/data/simulator_comparison.csv", index=False)
    print(f"\nResults saved to simulator_comparison.csv")
    
    # Optionally create plots - run the "plot_comparison_results.py" script
    
    return results_df


if __name__ == "__main__":
    results = main()
    print("\nComparison completed!")
    print(f"Shape of results: {results.shape}")
    print("\nSample results:")
    print(results.head())
