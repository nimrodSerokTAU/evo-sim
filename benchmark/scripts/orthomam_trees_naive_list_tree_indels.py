
import timeit
from pathlib import Path
import re

import pandas as pd

from indelsim.classes.sim_config import SimConfiguration
from indelsim.classes.simulation import Simulation


TREES_PATH = Path("benchmark/TRUE_TREES").resolve()

if not TREES_PATH.exists():
    print(f"Missing tree data at {TREES_PATH}")
    exit(1)


def scale_tree(tree_path: str, scale_factor: float, overwrite: bool=False):
    # Read the tree
    tree_path = Path(tree_path)

    with open(tree_path, 'r') as f:
        newick_str = f.read().strip()
    
    # Function to scale branch lengths in Newick string
    def scale_length(match):
        length = float(match.group(1))
        scaled = length * scale_factor
        return f":{scaled:.10f}"
    
    # Scale all branch lengths using regex
    scaled_newick = re.sub(r':([0-9.]+)', scale_length, newick_str)

    scaled_tree_path = (tree_path.parent / f"scaled_{scale_factor}.tree") if not overwrite else tree_path
    # Write the scaled tree directly
    with open(scaled_tree_path, 'w') as f:
        f.write(scaled_newick)
    return scaled_tree_path


def refresh_sim(tree_file, root_length):

    sim_config = SimConfiguration(original_sequence_length=root_length, indel_length_alpha=2.0,
                                indel_truncated_length=50,
                                rate_ins=0.03, rate_del=0.09,
                                deletion_extra_edge_length=50,
                                seed=111)
    sim = Simulation(input_tree=str(tree_file), config=sim_config)
    return sim

def run_simulation(sim: Simulation, type: str):
    if type == "Block list":
        sim.msa_from_blocklist()
    if type == "Block tree":
        sim.msa_from_blocktree()
    if type == "Naive":
        sim.msa_from_naive()

    sim.msa.compute_msa()

df_labels = ["tree", "Root length", "Method", "Time"]

algorithm_types = ["Block list", "Block tree",  "Naive"]


measured_times: list[list[float]] = []

for root_length in [50, 100, 500, 1000]:
    for tree_path in list(TREES_PATH.iterdir()):
        sim = refresh_sim(tree_path, root_length)
        for simulator_type in algorithm_types:
            simulation_time = timeit.timeit(lambda: run_simulation(sim, simulator_type), number=1)
            sim.msa.clear()
            row_data = [tree_path, root_length, simulator_type, simulation_time]
            measured_times.append(row_data)
        
    
    print(f"Done with root length: {root_length}")

times_df = pd.DataFrame(measured_times, columns=df_labels)

times_df.to_csv("benchmark/assets/data/algorithms_comparison_indels_orthomam.csv", index=False)


