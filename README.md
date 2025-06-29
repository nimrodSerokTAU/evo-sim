# Indel Simulator CLI Tool

A high-performance command-line tool for simulating insertion and deletion (indel) events along phylogenetic trees. This tool implements three different algorithms with varying computational complexities for handling indel events efficiently.

Based on the research paper: **"Efficient algorithms for simulating sequences along a phylogenetic tree"** by Elya Wygoda, Asher Moshe, Nimrod Serok, Edo Dotan, Noa Ecker, Omer Israeli, Itsik Pe'er, and Tal Pupko.

## Citation

If you use this tool in your research, please cite:

```
Wygoda, E., Moshe, A., Serok, N., Dotan, E., Ecker, N., Israeli, O., Pe'er, I., & Pupko, T. 
Efficient algorithms for simulating sequences along a phylogenetic tree. 
```

## Features

- **Multiple Algorithms**: Choose from naive O(k×n), block list O(k×b), or block tree O(k×log(b)) implementations
- **Flexible Parameters**: Customize insertion/deletion rates, length distributions, and sequence parameters
- **Performance Benchmarking**: Built-in timing and performance analysis
- **Reproducible Results**: Seed-based random number generation

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Install from Source

```bash
git clone https://github.com/nimrodSerokTAU/evo-sim.git
cd evo-sim
pip install -e .
```

## Quick Start

### Basic Usage

```bash
python indel_simulator.py \
    --type list \
    --insertion_rate 0.01 \
    --deletion_rate 0.01 \
    --tree_file tree.newick \
    --output_directory ./results
```

### Advanced Usage

```bash
python indel_simulator.py \
    --type tree \
    --insertion_rate 0.03 \
    --deletion_rate 0.09 \
    --insertion_length_distribution zipf \
    --insertion_length_distribution_parameter 2.0 \
    --insertion_length_truncation 50 \
    --deletion_length_distribution zipf \
    --deletion_length_distribution_parameter 2.0 \
    --deletion_length_truncation 50 \
    --original_sequence_length 1000 \
    --tree_file tree.newick \
    --number_of_simulations 100 \
    --output_type single_file \
    --output_directory ./results \
    --benchmark \
    --verbose
```

## Command Line Arguments

### Required Arguments

- `--type {naive,list,tree}`: Simulation algorithm type
- `--insertion_rate FLOAT`: Insertion rate per site per unit time
- `--deletion_rate FLOAT`: Deletion rate per site per unit time
- `--tree_file PATH`: Path to Newick format phylogenetic tree file

### Length Distribution Arguments

- `--insertion_length_distribution {zipf,geometric}`: Distribution for insertion lengths (default: zipf)
- `--insertion_length_distribution_parameter FLOAT`: Parameter for insertion length distribution (default: 2.0)
- `--insertion_length_truncation INT`: Maximum insertion length (default: 50)
- `--deletion_length_distribution {zipf,geometric}`: Distribution for deletion lengths (default: zipf)
- `--deletion_length_distribution_parameter FLOAT`: Parameter for deletion length distribution (default: 2.0)
- `--deletion_length_truncation INT`: Maximum deletion length (default: 50)

### Simulation Parameters

- `--original_sequence_length INT`: Length of the root sequence (default: 1000)
- `--deletion_extra_edge_length INT`: Extra positions before sequence start for deletion events (default: 49)
- `--switch_factor FLOAT`: Threshold for switching between list and tree methods in hybrid mode (default: 100.0)
- `--number_of_simulations INT`: Number of independent simulation runs (default: 1)
- `--seed INT`: Random seed for reproducibility (default: 42)

### Output Options

- `--output_type {drop_output,multiple_files,single_file}`: Output format (default: single_file)
- `--output_directory PATH`: Directory to save simulation results (default: ./simulation_results)
- `--output_format {fasta}`: Output file format (always FASTA)
- `--verbose`: Enable verbose output
- `--benchmark`: Run benchmarking and report performance statistics

## Algorithm Comparison

| Algorithm | Time Complexity | Memory Usage | Best Use Case |
|-----------|----------------|--------------|---------------|
| Naive     | O(k×n)         | O(n)         | Small sequences, few events |
| Block List| O(k×b)         | O(b)         | Medium sequences, moderate events |
| Block Tree| O(k×log(b))    | O(b)         | Large sequences, many events |

Where:
- k = number of indel events
- n = sequence length  
- b = number of blocks (≤ min(k,n))

## Output Formats

### FASTA Format
```
# Simulation 1
# Runtime: 0.123s
# Simulation type: tree
# Insertion rate: 0.01
# Deletion rate: 0.01
# Seed: 42
>1
ATGCGATCGATCG--ATCGATCG
>2
ATGC--TCGATCGATCGATCG
```

The simulator outputs multiple sequence alignments in FASTA format with header comments containing simulation metadata.

## Examples

See the example scripts in the repository:
- `example_usage.sh`: Basic usage examples
- `performance_comparison.py`: Algorithm performance comparison
- `batch_simulations.sh`: Batch processing with different parameters

## Performance Tips

1. **For small sequences (< 500 bp)**: Use `--type naive` for simplicity
2. **For medium sequences (500-5000 bp)**: Use `--type list` for good performance
3. **For large sequences (> 5000 bp) or many events**: Use `--type tree` for best performance
4. **For benchmarking**: Use `--benchmark` flag to get detailed timing information
5. **For reproducibility**: Always set `--seed` to a fixed value

## Tree File Format

The tool accepts phylogenetic trees in Newick format:

```
((A:0.1,B:0.1):0.05,(C:0.08,D:0.12):0.05);
```

Branch lengths represent evolutionary time in substitutions per site.

## Integration with SpartaABC

This simulator can be integrated with SpartaABC for Approximate Bayesian Computation-based inference of indel parameters. Visit [https://spartaabc.tau.ac.il/](https://spartaabc.tau.ac.il/) for more information.

## Contributing

We welcome contributions! Please see our contributing guidelines and submit pull requests via GitHub.

## License

This project is licensed under the Academic Free License v. 3.0. See the LICENSE file for details.

## Contact

For questions or support, please contact:
- Tal Pupko: talp@tauex.tau.ac.il
- GitHub Issues: [https://github.com/nimrodSerokTAU/evo-sim/issues](https://github.com/nimrodSerokTAU/evo-sim/issues)

## Acknowledgments

- The Shmunis School of Biomedicine and Cancer Research, Tel Aviv University
- The Henry and Marilyn Taub Faculty of Computer Science, Technion
- Department of Computer Science, Columbia University

This research was supported by the Israel Science Foundation (ISF) [2818/21 to T.P.].