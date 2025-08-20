#!/usr/bin/env python3
"""
Indel Algorithm Performance Comparison

Compares naive, list, and tree indel algorithms to demonstrate
that indel simulation is the rate limiting step.
"""

import subprocess
import pandas as pd
import tempfile
from pathlib import Path

# Configuration
OUTPUT_DIR = "benchmark/assets/data/"
ALGORITHMS = ['naive', 'list', 'tree']
# INDEL_RATES = [0.01, 0.05, 0.1, 0.5]  # Combined insertion + deletion rate
INSERTION_RATE = 0.03
DELETION_RATE = 0.09
SEQUENCE_LENGTHS = [100, 500, 1000, 5000]
BRANCH_LENGTHS = [0.01, 0.05, 0.1, 0.5]
NUM_SIMULATIONS = 3
TREE_FILE = "benchmark/scaled_trees/scaled_1.tree"  # Path to your tree file


def create_test_tree(branch_length=0.1):
    """Read tree file and set all branches to specified length"""
    try:
        with open(TREE_FILE, 'r') as f:
            tree_content = f.read().strip()
        
        # Replace all branch lengths with the new value
        # This regex finds patterns like :0.123 and replaces with :branch_length
        import re
        tree_with_new_branches = re.sub(r':\d*\.?\d+', f':{branch_length}', tree_content)
        return tree_with_new_branches
        
    except FileNotFoundError:
        print(f"Warning: Tree file {TREE_FILE} not found")
    
def run_simulation(algorithm, ins_rate, del_rate, seq_len, branch_len, tree_content):
    """Run single simulation and extract timing"""
    with tempfile.TemporaryDirectory() as temp_dir:
        tree_file = Path(temp_dir) / "tree.newick"
        tree_file.write_text(tree_content)
        cmd = [
            "msa-simulator",
            "--type", algorithm,
            "--insertion_rate", str(ins_rate),
            "--deletion_rate", str(del_rate),
            "--original_sequence_length", str(seq_len),
            "--tree_file", str(tree_file),
            "--number_of_simulations", str(NUM_SIMULATIONS),
            "--keep_in_memory",
            "--output_type", "drop_output",
            "--benchmark",
        ]
        # print(" ".join(cmd))
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                return None
            # Parse timing from output
            indel_time = sub_time = total_time = 0.0
            for line in result.stdout.split('\n'):
                if "Indel simulation:" in line:
                    indel_time = float(line.split(':')[1].split('s')[0].strip())
                elif "Substitution simulation:" in line:
                    sub_time = float(line.split(':')[1].split('s')[0].strip())
                elif "Total runtime:" in line:
                    total_time = float(line.split(':')[1].split('s')[0].strip())
            
            return {
                'algorithm': algorithm,
                'insertion_rate': ins_rate,
                'deletion_rate': del_rate,
                'sequence_length': seq_len,
                'branch_length': branch_len,
                'indel_time': indel_time,
                'substitution_time': sub_time,
                'total_time': total_time,
                'indel_ratio': indel_time / (indel_time + sub_time) if (indel_time + sub_time) > 0 else 0,
                'speedup_vs_naive': 1.0
            }
        except Exception as e:
            print(e)
            return None

def main():
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Indel Algorithm Performance Comparison")
    print("=" * 50)
    INDEL_RATES = INSERTION_RATE + DELETION_RATE
    results = []
    total_runs = len(ALGORITHMS) * len(SEQUENCE_LENGTHS) * len(BRANCH_LENGTHS)
    current_run = 0
    
    for algorithm in ALGORITHMS:
        print(f"\nTesting {algorithm} algorithm...")
        
        # for indel_rate in INDEL_RATES:
        for seq_len in SEQUENCE_LENGTHS:
            for branch_len in BRANCH_LENGTHS:
                current_run += 1
                ins_rate = INSERTION_RATE
                del_rate = DELETION_RATE
                print(f"  [{current_run:2d}/{total_runs}] "
                        f"indel={INDEL_RATES:.3f}, len={seq_len}, branch={branch_len:.3f}", 
                        end=" ... ")
                
                tree_content = create_test_tree(branch_len)
                result = run_simulation(algorithm, ins_rate, del_rate, seq_len, branch_len, tree_content)
                
                if result:
                    results.append(result)
                    print(f"✓ indel_ratio={result['indel_ratio']:.2f}")
                else:
                    print("✗")
    
    # Calculate speedups vs naive
    df = pd.DataFrame(results)
    for _, group in df.groupby(['insertion_rate', 'deletion_rate', 'sequence_length', 'branch_length']):
        naive_time = group[group['algorithm'] == 'naive']['indel_time'].iloc[0] if 'naive' in group['algorithm'].values else 1.0
        for idx, row in group.iterrows():
            df.at[idx, 'speedup_vs_naive'] = naive_time / row['indel_time'] if row['indel_time'] > 0 else 1.0
    
    # Save results
    output_file = output_dir / "naive_list_tree_comparison.csv"
    df.to_csv(output_file, index=False)
    
    # Print summary
    print(f"\n{'='*50}")
    print("RESULTS SUMMARY")
    print(f"{'='*50}")
    print(f"Total successful runs: {len(df)}")
    
    for alg in ALGORITHMS:
        alg_data = df[df['algorithm'] == alg]
        if len(alg_data) > 0:
            print(f"  {alg:>6}: avg_indel_time={alg_data['indel_time'].mean():.4f}s, "
                  f"avg_ratio={alg_data['indel_ratio'].mean():.3f}, "
                  f"avg_speedup={alg_data['speedup_vs_naive'].mean():.2f}x")
    
    indel_dominant = (df['indel_ratio'] > 0.5).mean() * 100
    print(f"\nIndel is rate limiting in {indel_dominant:.1f}% of cases")
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()