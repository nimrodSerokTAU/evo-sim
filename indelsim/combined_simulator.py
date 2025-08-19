#!/usr/bin/env python3
"""
Combined Indel and Substitution Simulator CLI Tool

This tool combines indel and substitution simulations in a two-step process:
1. First performs indel simulation to create the MSA template
2. Then performs substitution simulation on the resulting MSA length
3. Finally replaces X placeholders with actual amino acid sequences
"""

import argparse
import sys
import os
import pathlib
import traceback
from typing import List, Dict, Any, Tuple
import time
from datetime import datetime
import numpy as np

# Import existing CLI classes to reuse their functionality
from indelsim.indel_simulator import IndelSimulatorCLI, TEMP_FILE_NAME as TEMP_INDEL_FILE
from indelsim.substitution_simulator import SubstitutionSimulatorCLI, TEMP_FILE_NAME as TEMP_SUBS_FILE


class CombinedSimulatorCLI:
    """Command-line interface for the combined indel and substitution simulator."""
    
    def __init__(self):
        self.indel_cli = IndelSimulatorCLI()
        self.substitution_cli = SubstitutionSimulatorCLI()
        self.parser = self._create_parser()
    
    def _create_parser(self) -> argparse.ArgumentParser:
        """Create combined argument parser using both existing parsers."""
        parser = argparse.ArgumentParser(
            description="Combined indel and substitution simulation along phylogenetic trees",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Basic combined simulation
  python combined_simulator.py --type list --insertion_rate 0.01 --deletion_rate 0.01 \\
                               --substitution_rate 1.0 --tree_file tree.newick \\
                               --output_directory ./results
  
  # Advanced combined simulation
  python combined_simulator.py --type tree --insertion_rate 0.03 --deletion_rate 0.09 \\
                               --substitution_rate 2.0 --algorithm matrix \\
                               --insertion_length_distribution_parameter 2.0 \\
                               --deletion_length_distribution_parameter 2.0 \\
                               --original_sequence_length 500 \\
                               --number_of_simulations 10 --tree_file tree.newick \\
                               --output_directory ./combined_results --benchmark --verbose
            """
        )
        
        # Get arguments from both existing parsers
        indel_parser = self.indel_cli.parser
        substitution_parser = self.substitution_cli.parser
        
        # Add all indel arguments
        for action in indel_parser._actions:
            if action.dest not in ['help', 'output_directory', 'number_of_simulations', 'seed', 
                                   'output_type', 'verbose', 'benchmark', 'tree_file', 'original_sequence_length',
                                   'keep_in_memory']:
                parser.add_argument(*action.option_strings, **{
                    'type': action.type,
                    'default': action.default,
                    'help': action.help,
                    'choices': action.choices,
                    'required': action.required
                })
        
        # Add substitution-specific arguments (avoiding duplicates)
        for action in substitution_parser._actions:
            if action.dest in ['substitution_rate', 'algorithm']:
                parser.add_argument(*action.option_strings, **{
                    'type': action.type,
                    'default': action.default,
                    'help': action.help,
                    'choices': action.choices,
                    'required': action.required
                })
        
        # Add common arguments (from either parser, avoiding duplicates)
        common_args = ['tree_file', 'original_sequence_length', 'number_of_simulations', 
                      'seed', 'output_type', 'output_directory', 'verbose', 'benchmark', 'keep_in_memory']
        
        for action in indel_parser._actions:
            if action.dest in common_args:
                required = action.dest == 'tree_file'
                default_dir = "./combined_simulation_results" if action.dest == 'output_directory' else action.default
                
                kwargs = {'default': default_dir, 'help': action.help, 'required': required}
                if hasattr(action, 'choices') and action.choices:
                    kwargs['choices'] = action.choices
                if isinstance(action, argparse._StoreTrueAction):
                    kwargs['action'] = 'store_true'
                elif isinstance(action, argparse._StoreFalseAction):
                    kwargs['action'] = 'store_false'
                elif action.type is not None:
                    kwargs['type'] = action.type
                parser.add_argument(*action.option_strings, **kwargs)
        
        return parser
    
    def _init_output_file(self, args: argparse.Namespace) -> None:
        """Save simulation results to files."""
        args.output_directory = pathlib.Path(args.output_directory)
        if args.output_type == "drop_output":
            return
        
        # Create output directory
        args.output_directory.mkdir(parents=True, exist_ok=True)
        with open(args.output_directory / TEMP_INDEL_FILE, 'w') as f:
            f.write("")
        with open(args.output_directory / TEMP_SUBS_FILE, 'w') as f:
            f.write("")



    
    def _run_indel_simulation(self, args: argparse.Namespace, sim_num: int) -> Tuple[Dict[str, Any], int]:
        """
        Run indel simulation using existing IndelSimulatorCLI and return results with MSA length.
        """
        if args.verbose:
            print(f"  Step 1: Running indel simulation...")
        
        start_time = time.perf_counter()
        
        # Use the existing indel simulator method
        indel_result = self.indel_cli._run_single_simulation(args, sim_num)
        
        end_time = time.perf_counter()
        indel_runtime = end_time - start_time
        
        # Extract MSA length from the indel result
        msa_length = indel_result["msa"]._msa_length
        
        if args.verbose:
            print(f"    Indel simulation completed in {indel_runtime:.3f} seconds")
            print(f"    Final MSA length: {msa_length}")
        
        return indel_result, msa_length
    
    def _run_substitution_simulation(self, args: argparse.Namespace, msa_length: int, sim_num: int) -> Dict[str, Any]:
        """
        Run substitution simulation using existing SubstitutionSimulatorCLI.
        """
        if args.verbose:
            print(f"  Step 2: Running substitution simulation on MSA length {msa_length}...")
        
        # Create modified args for substitution simulation with the correct sequence length
        sub_args = argparse.Namespace(**vars(args))
        sub_args.output_directory = pathlib.Path(sub_args.output_directory)
        sub_args.original_sequence_length = msa_length
        
        start_time = time.perf_counter()
        
        # Use the existing substitution simulator method
        substitution_result = self.substitution_cli._run_single_simulation(sub_args, sim_num)
        
        end_time = time.perf_counter()
        substitution_runtime = end_time - start_time
        
        if args.verbose:
            print(f"    Substitution simulation completed in {substitution_runtime:.3f} seconds")
        
        return substitution_result
    
    def _merge_simulations(self, indel_result: Dict[str, Any], substitution_result: Dict[str, Any]) -> Dict[str, str]:
        """
        Merge indel template with substitution sequences by replacing X placeholders.
        """
        # Get the aligned sequences from indel simulation
        indel_sequences = indel_result["msa"]._aligned_sequences
        id_to_name = indel_result["msa"]._id_to_name
        
        # Get the substitution sequences  
        substitution_sequences = substitution_result["msa"]
        # Replace X placeholders with actual amino acids
        merged_sequences = {}
        
        for seq_id, template_seq in indel_sequences.items():
            # print(template_seq)
            if seq_id in substitution_sequences:
                substitution_seq = (substitution_sequences[seq_id])
                for idx,c in enumerate(template_seq):
                    if c=="-":
                        substitution_seq[idx] = '-'
                merged_sequences[id_to_name[seq_id]] = ''.join(substitution_seq)
            else:
                # Keep original template if no substitution sequence available
                merged_sequences[id_to_name[seq_id]] = template_seq
        
        return merged_sequences
    
    def _run_single_simulation(self, args: argparse.Namespace, sim_num: int) -> Dict[str, Any]:
        """Run a single combined simulation and return results."""
        if args.verbose:
            print(f"Running combined simulation {sim_num + 1}/{args.number_of_simulations}...")
        
        total_start_time = time.perf_counter()
        
        # Step 1: Run indel simulation
        indel_result, msa_length = self._run_indel_simulation(args, sim_num)
        
        # Step 2: Run substitution simulation
        substitution_result = self._run_substitution_simulation(args, msa_length, sim_num)
        
        # Step 3: Merge simulations
        # if args.verbose:
        #     print(f"  Step 3: Merging indel template with substitution sequences...")
        
        # merge_start_time = time.perf_counter()
        merged_sequences = None
        if args.keep_in_memory:
            merged_sequences = self._merge_simulations(indel_result, substitution_result)
        # merge_time = time.perf_counter() - merge_start_time
        
        total_end_time = time.perf_counter()
        total_runtime = total_end_time - total_start_time
        
        # Collect combined results
        results = {
            "simulation_number": sim_num + 1,
            "total_runtime_seconds": total_runtime,
            "indel_runtime_seconds": indel_result["runtime_seconds"],
            "substitution_runtime_seconds": substitution_result["runtime_seconds"],
            "indel_config": indel_result["config"],
            "substitution_config": substitution_result["config"],
            "final_msa": merged_sequences,
            "indel_events": indel_result.get("events", []),
            "indel_type": indel_result["simulation_type"],
            "substitution_algorithm": substitution_result["algorithm"]
        }
        
        if args.verbose:
            print(f"  Total simulation completed in {total_runtime:.3f} seconds")
            print(f"    Indel time: {indel_result['runtime_seconds']:.3f}s")
            print(f"    Substitution time: {substitution_result['runtime_seconds']:.3f}s")
        
        return results
    
    
    def _save_multiple_files(self, result: List[Dict[str, Any]], args: argparse.Namespace, output_dir: pathlib.Path) -> None:
        """Save each simulation to a separate file."""
        sim_num = result["simulation_number"]
        
        output_msa_path = output_dir / f"combined_sim_{sim_num:04d}.fasta"
        (output_dir / TEMP_SUBS_FILE).rename(output_msa_path)
        
        if args.keep_in_memory:
            self._write_fasta(result["final_msa"], output_msa_path)

        if args.verbose:
            print(f"Saved simulation {sim_num} to {output_msa_path}")

    def _save_single_file(self, result: List[Dict[str, Any]], args: argparse.Namespace, output_dir: pathlib.Path) -> None:
        """Save all simulations to a single file."""

        temp_path = output_dir / TEMP_SUBS_FILE
        

        with open(temp_path, 'a') as f:
            if args.keep_in_memory:
                msa = result["final_msa"]
                for species_name, sequence in msa.items():
                    f.write(f">{species_name}\n")
                    f.write(''.join(sequence))
                    f.write("\n")

            f.write(f"# Combined Simulation {result['simulation_number']}\n")
            f.write(f"# Total Runtime: {result['total_runtime_seconds']:.3f}s\n")
            f.write(f"# Indel Runtime: {result['indel_runtime_seconds']:.3f}s\n")
            f.write(f"# Substitution Runtime: {result['substitution_runtime_seconds']:.3f}s\n")
            f.write(f"# Indel Type: {result['indel_type']}\n")
            f.write(f"# Substitution Algorithm: {result['substitution_algorithm']}\n")
            f.write(f"# Insertion Rate: {result['indel_config']['insertion_rate']}\n")
            f.write(f"# Deletion Rate: {result['indel_config']['deletion_rate']}\n")
            f.write(f"# Substitution Rate: {result['substitution_config']['substitution_rate']}\n")
            f.write(f"# Original Sequence Length: {result['indel_config']['original_sequence_length']}\n")
            f.write(f"# Seed: {result['indel_config']['seed']}\n")
            
            # Write MSA sequences
            f.write("\n\n")
        
        if args.verbose:
            print(f"Saved simulation to {temp_path}")

                
    
    def _write_fasta(self, msa: Dict[str, str], filename: pathlib.Path, result: Dict[str, Any]) -> None:
        """Write MSA in FASTA format with metadata."""
        with open(filename, 'w') as f:
            # Write metadata as comments
            f.write(f"# Combined Simulation {result['simulation_number']}\n")
            f.write(f"# Total Runtime: {result['total_runtime_seconds']:.3f}s\n")
            f.write(f"# Indel Runtime: {result['indel_runtime_seconds']:.3f}s\n")
            f.write(f"# Substitution Runtime: {result['substitution_runtime_seconds']:.3f}s\n")
            f.write(f"# Merge Runtime: {result['merge_runtime_seconds']:.3f}s\n")
            f.write(f"# Indel Type: {result['indel_type']}\n")
            f.write(f"# Substitution Algorithm: {result['substitution_algorithm']}\n")
            
            # Write sequences
            for species_id, sequence in msa.items():
                f.write(f">{species_id}\n")
                f.write(sequence)
                f.write("\n")
    
    def _print_benchmark_results(self, results: List[Dict[str, Any]], args: argparse.Namespace) -> None:
        """Print benchmarking statistics."""
        total_runtimes = [r["total_runtime_seconds"] for r in results]
        indel_runtimes = [r["indel_runtime_seconds"] for r in results]
        substitution_runtimes = [r["substitution_runtime_seconds"] for r in results]
        
        print("\n" + "="*60)
        print("COMBINED SIMULATION BENCHMARK RESULTS")
        print("="*60)
        print(f"Indel Algorithm: {args.type}")
        print(f"Substitution Algorithm: {args.algorithm}")
        print(f"Number of simulations: {len(results)}")
        print(f"Original sequence length: {args.original_sequence_length}")
        print(f"Insertion rate: {args.insertion_rate}")
        print(f"Deletion rate: {args.deletion_rate}")
        print(f"Substitution rate: {args.substitution_rate}")
        # if len(results) > 0:
        #     final_length = len(list(results[0]["final_msa"].values())[0])
        #     print(f"Average final MSA length: {final_length}")
        print("-"*60)
        print("TIMING BREAKDOWN:")
        print(f"Total runtime: {sum(total_runtimes):.3f}s")
        if sum(total_runtimes) > 0:
            print(f"  Indel simulation: {sum(indel_runtimes):.3f}s ({sum(indel_runtimes)/sum(total_runtimes)*100:.1f}%)")
            print(f"  Substitution simulation: {sum(substitution_runtimes):.3f}s ({sum(substitution_runtimes)/sum(total_runtimes)*100:.1f}%)")
        print("-"*60)
        print("AVERAGE TIMES PER SIMULATION:")
        print(f"Total average: {sum(total_runtimes)/len(total_runtimes):.3f}s")
        print(f"  Indel average: {sum(indel_runtimes)/len(indel_runtimes):.3f}s")
        print(f"  Substitution average: {sum(substitution_runtimes)/len(substitution_runtimes):.3f}s")
        print("-"*60)
        print(f"Min total runtime: {min(total_runtimes):.3f}s")
        print(f"Max total runtime: {max(total_runtimes):.3f}s")
        if len(total_runtimes) > 1:
            import statistics
            print(f"Std deviation: {statistics.stdev(total_runtimes):.3f}s")
        print("="*60)
    
    def run(self) -> None:
        """Main entry point for the CLI."""
        args = self.parser.parse_args()
        
        # Validate arguments
        if not os.path.exists(args.tree_file):
            print(f"Error: Tree file '{args.tree_file}' not found.", file=sys.stderr)
            sys.exit(1)
        
        if args.number_of_simulations < 1:
            print("Error: Number of simulations must be at least 1.", file=sys.stderr)
            sys.exit(1)
        
        if args.verbose:
            print("Starting combined indel and substitution simulations...")
            print(f"Indel algorithm: {args.type}")
            print(f"Substitution algorithm: {args.algorithm}")
            print(f"Number of simulations: {args.number_of_simulations}")
            print(f"Tree file: {args.tree_file}")
            print(f"Output directory: {args.output_directory}")
            print()
        
        # Run simulations
        results = []
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self._init_output_file(args)
        combined_file_path = args.output_directory / f"combined_simulations_{timestamp}.fasta"

        for i in range(args.number_of_simulations):
            result = self._run_single_simulation(args, i)
            results.append(result)
            if args.output_type == "multiple_files":
                self._save_multiple_files(result, args, args.output_directory)
                self._init_output_file(args)

            elif args.output_type == "single_file":
                self._save_single_file(result, args, args.output_directory)
                (args.output_directory / TEMP_INDEL_FILE).unlink(missing_ok=True)
        
        if args.output_type == "single_file":
            (args.output_directory / TEMP_SUBS_FILE).rename(combined_file_path)
        (args.output_directory / TEMP_SUBS_FILE).unlink(missing_ok=True)
        (args.output_directory / TEMP_INDEL_FILE).unlink(missing_ok=True)

        # Print benchmark results if requested
        if args.benchmark:
            self._print_benchmark_results(results, args)
        
        if args.verbose:
            print(f"\nCompleted {len(results)} combined simulations successfully!")

def main():
    cli = CombinedSimulatorCLI()
    cli.run()

if __name__ == "__main__":
    main()

