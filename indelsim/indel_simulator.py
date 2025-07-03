#!/usr/bin/env python3
"""
CLI Indel Simulator Tool
A command-line interface for simulating indel events along phylogenetic trees
using different data structures (naive, list, tree) for performance optimization.

This tool implements the algorithms described in:
"Efficient algorithms for simulating sequences along a phylogenetic tree"
by Elya Wygoda, Asher Moshe, Nimrod Serok, Edo Dotan, Noa Ecker, Omer Israeli, 
Itsik Pe'er, and Tal Pupko.

Authors: Elya Wygoda, Asher Moshe, Nimrod Serok, Edo Dotan, Noa Ecker, 
         Omer Israeli, Itsik Pe'er, Tal Pupko
Contact: talp@tauex.tau.ac.il
License: Academic Free License v. 3.0
Repository: https://github.com/nimrodSerokTAU/evo-sim

Licensed under the Academic Free License version 3.0
"""

import argparse
import sys
import os
import pathlib
from typing import List, Dict, Any
import time
from datetime import datetime

from indelsim.classes.simulation import Simulation
from indelsim.classes.sim_config import SimConfiguration
from indelsim.enums import SimulationTypes


class IndelSimulatorCLI:
    """Command-line interface for the indel simulator."""
    
    def __init__(self):
        self.parser = self._create_parser()
    
    def _create_parser(self) -> argparse.ArgumentParser:
        """Create and configure the argument parser."""
        parser = argparse.ArgumentParser(
            description="Simulate indel events along phylogenetic trees",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Basic simulation with block list approach
  python indel_simulator.py --type list --insertion_rate 0.01 --deletion_rate 0.01 
                           --tree_file tree.newick --output_directory ./results
  
  # Advanced simulation with custom length distributions
  python indel_simulator.py --type tree --insertion_rate 0.03 --deletion_rate 0.09
                           --insertion_length_distribution zipf --insertion_length_distribution_parameter 2.0
                           --deletion_length_distribution zipf --deletion_length_distribution_parameter 2.0
                           --tree_file tree.newick --number_of_simulations 100
                           --output_type single_file --output_directory ./sim_results
            """
        )
        
        # Required arguments
        parser.add_argument(
            "--type",
            choices=["naive", "list", "tree"],
            required=True,
            help="Simulation algorithm type: naive (O(k*n)), list (O(k*b)), or tree (O(k*log(b)))"
        )
        
        parser.add_argument(
            "--insertion_rate",
            type=float,
            required=True,
            help="Insertion rate per site per unit time"
        )
        
        parser.add_argument(
            "--deletion_rate",
            type=float,
            required=True,
            help="Deletion rate per site per unit time"
        )
        
        parser.add_argument(
            "--tree_file",
            type=str,
            required=True,
            help="Path to Newick format phylogenetic tree file"
        )
        
        # Length distribution arguments
        parser.add_argument(
            "--insertion_length_distribution_parameter",
            type=float,
            default=2.0,
            help="Parameter for insertion length distribution (default: 2.0)"
        )
        
        parser.add_argument(
            "--insertion_length_truncation",
            type=int,
            default=50,
            help="Maximum insertion length (default: 50)"
        )
        
        
        parser.add_argument(
            "--deletion_length_distribution_parameter",
            type=float,
            default=2.0,
            help="Parameter for deletion length distribution (default: 2.0)"
        )
        
        parser.add_argument(
            "--deletion_length_truncation",
            type=int,
            default=50,
            help="Maximum deletion length (default: 50)"
        )
        
        # Simulation parameters
        parser.add_argument(
            "--original_sequence_length",
            type=int,
            default=1000,
            help="Length of the root sequence (default: 1000)"
        )
        
        parser.add_argument(
            "--deletion_extra_edge_length",
            type=int,
            default=49,
            help="Extra positions before sequence start for deletion events (default: 49)"
        )
        
        parser.add_argument(
            "--number_of_simulations",
            type=int,
            default=1,
            help="Number of independent simulation runs (default: 1)"
        )
        
        parser.add_argument(
            "--seed",
            type=int,
            default=42,
            help="Random seed for reproducibility (default: 42)"
        )
        
        # Output options
        parser.add_argument(
            "--output_type",
            choices=["drop_output", "multiple_files", "single_file"],
            default="single_file",
            help="Output format: drop_output (no files), multiple_files (one fasta per sim), single_file (combined fastas)"
        )
        
        parser.add_argument(
            "--output_directory",
            type=str,
            default="./simulation_results",
            help="Directory to save simulation results (default: ./simulation_results)"
        )
        
        parser.add_argument(
            "--verbose",
            action="store_true",
            help="Enable verbose output"
        )
        
        parser.add_argument(
            "--benchmark",
            action="store_true",
            help="Run benchmarking and report performance statistics"
        )
        
        return parser
    
    def _validate_args(self, args: argparse.Namespace) -> None:
        """Validate command-line arguments."""
        # Check if tree file exists
        if not os.path.exists(args.tree_file):
            raise FileNotFoundError(f"Tree file not found: {args.tree_file}")
        
        # Validate rates
        if args.insertion_rate < 0 or args.deletion_rate < 0:
            raise ValueError("Insertion and deletion rates must be non-negative")
        
        # Validate length distribution parameters
        if args.insertion_length_distribution_parameter <= 0:
            raise ValueError("Insertion length distribution parameter must be positive")
        
        if args.deletion_length_distribution_parameter <= 0:
            raise ValueError("Deletion length distribution parameter must be positive")
        
        # Validate truncation values
        if args.insertion_length_truncation <= 0 or args.deletion_length_truncation <= 0:
            raise ValueError("Length truncation values must be positive")
        
        # Validate sequence length
        if args.original_sequence_length <= 0:
            raise ValueError("Original sequence length must be positive")
        
        # Validate number of simulations
        if args.number_of_simulations <= 0:
            raise ValueError("Number of simulations must be positive")
    
    def _create_sim_config(self, args: argparse.Namespace) -> SimConfiguration:
        """Create simulation configuration from arguments."""
        return SimConfiguration(
            original_sequence_length=args.original_sequence_length,
            indel_length_alpha=args.insertion_length_distribution_parameter,  # Using insertion param as default
            indel_truncated_length=args.insertion_length_truncation,
            deletion_extra_edge_length=args.deletion_extra_edge_length,
            rate_ins=args.insertion_rate,
            rate_del=args.deletion_rate,
            seed=args.seed
        )
    
    def _get_simulation_type(self, sim_type: str) -> SimulationTypes:
        """Convert string simulation type to enum."""
        type_mapping = {
            "naive": SimulationTypes.NAIVE,
            "list": SimulationTypes.BLOCK_LIST,
            "tree": SimulationTypes.BLOCK_TREE
        }
        return type_mapping[sim_type]
    
    def _run_single_simulation(self, args: argparse.Namespace, sim_num: int) -> Dict[str, Any]:
        """Run a single simulation and return results."""
        if args.verbose:
            print(f"Running simulation {sim_num + 1}/{args.number_of_simulations}...")
        
        # Create configuration with incremented seed for each simulation
        config = self._create_sim_config(args)
        config.random_seed = args.seed + sim_num
        
        start_time = time.time()
        
        # Create and run simulation
        simulation = Simulation(args.tree_file, config)
        
        # Choose simulation method based on type
        sim_type = args.type
        if sim_type == "naive":
            simulation.msa_from_naive()
        elif sim_type == "list":
            simulation.msa_from_blocklist()
        elif sim_type == "tree":
            simulation.msa_from_blocktree()
        
        end_time = time.time()
        runtime = end_time - start_time
        
        # Collect results
        results = {
            "simulation_number": sim_num + 1,
            "runtime_seconds": runtime,
            "simulation_type": sim_type,
            "config": {
                "insertion_rate": args.insertion_rate,
                "deletion_rate": args.deletion_rate,
                "original_sequence_length": args.original_sequence_length,
                "seed": config.random_seed
            },
            "msa": simulation.msa,
            "events": [str(node) for node in simulation.get_events()]
        }
        
        if args.verbose:
            print(f"  Completed in {runtime:.3f} seconds")
        
        return results
    
    def _save_results(self, results: List[Dict[str, Any]], args: argparse.Namespace) -> None:
        """Save simulation results to files."""
        if args.output_type == "drop_output":
            return
        
        # Create output directory
        output_dir = pathlib.Path(args.output_directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if args.output_type == "multiple_files":
            self._save_multiple_files(results, args, output_dir)
        elif args.output_type == "single_file":
            self._save_single_file(results, args, output_dir)
    
    def _save_multiple_files(self, results: List[Dict[str, Any]], args: argparse.Namespace, output_dir: pathlib.Path) -> None:
        """Save each simulation to a separate file."""
        for result in results:
            sim_num = result["simulation_number"]
            
            filename = output_dir / f"simulation_{sim_num:04d}.fasta"
            self._write_fasta(result["msa"], filename)
            
            if args.verbose:
                print(f"Saved simulation {sim_num} to {filename}")
    
    def _save_single_file(self, results: List[Dict[str, Any]], args: argparse.Namespace, output_dir: pathlib.Path) -> None:
        """Save all simulations to a single file."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        filename = output_dir / f"combined_simulations_{timestamp}.fasta"
        with open(filename, 'w') as f:
            for i, result in enumerate(results):
                f.write(f"# Simulation {result['simulation_number']}\n")
                f.write(f"# Runtime: {result['runtime_seconds']:.3f}s\n")
                f.write(f"# Simulation type: {result['simulation_type']}\n")
                f.write(f"# Insertion rate: {result['config']['insertion_rate']}\n")
                f.write(f"# Deletion rate: {result['config']['deletion_rate']}\n")
                f.write(f"# Seed: {result['config']['seed']}\n")
                if isinstance(result["msa"], str):
                    f.write(result["msa"])
                else:
                    f.write(str(result["msa"]))
                f.write("\n\n")
        
        if args.verbose:
            print(f"Saved {len(results)} simulations to {filename}")
    
    def _write_fasta(self, msa, filename: pathlib.Path) -> None:
        """Write MSA in FASTA format."""
        with open(filename, 'w') as f:
            if isinstance(msa, str):
                f.write(msa)
            else:
                f.write(str(msa))
    
    def _print_benchmark_results(self, results: List[Dict[str, Any]], args: argparse.Namespace) -> None:
        """Print benchmarking statistics."""
        runtimes = [r["runtime_seconds"] for r in results]
        
        print("\n" + "="*50)
        print("BENCHMARK RESULTS")
        print("="*50)
        print(f"Simulation type: {args.type}")
        print(f"Number of simulations: {len(results)}")
        print(f"Original sequence length: {args.original_sequence_length}")
        print(f"Insertion rate: {args.insertion_rate}")
        print(f"Deletion rate: {args.deletion_rate}")
        print("-"*50)
        print(f"Total runtime: {sum(runtimes):.3f}s")
        print(f"Average runtime: {sum(runtimes)/len(runtimes):.3f}s")
        print(f"Min runtime: {min(runtimes):.3f}s")
        print(f"Max runtime: {max(runtimes):.3f}s")
        if len(runtimes) > 1:
            import statistics
            print(f"Std deviation: {statistics.stdev(runtimes):.3f}s")
        print("="*50)
    
    def run(self) -> None:
        """Main entry point for the CLI."""
        try:
            args = self.parser.parse_args()
            self._validate_args(args)
            
            if args.verbose:
                print(f"Starting indel simulation with {args.type} algorithm...")
                print(f"Tree file: {args.tree_file}")
                print(f"Number of simulations: {args.number_of_simulations}")
                print(f"Output directory: {args.output_directory}")
            
            # Run simulations
            results = []
            total_start_time = time.time()
            
            for i in range(args.number_of_simulations):
                result = self._run_single_simulation(args, i)
                results.append(result)
            
            total_end_time = time.time()
            
            # Save results
            self._save_results(results, args)
            
            # Print benchmark results if requested
            if args.benchmark or args.verbose:
                self._print_benchmark_results(results, args)
            
            if args.verbose:
                print(f"\nAll simulations completed in {total_end_time - total_start_time:.3f}s")
                print(f"Results saved to: {args.output_directory}")
        
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)


def main():
    """Main function."""
    simulator = IndelSimulatorCLI()
    simulator.run()


if __name__ == "__main__":
    main()