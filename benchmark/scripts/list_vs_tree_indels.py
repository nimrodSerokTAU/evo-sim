
import timeit
import re
from pathlib import Path

import pandas as pd

from indelsim.classes.sim_config import SimConfiguration
from indelsim.classes.simulation import Simulation



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


ROOT_SEQUENCE_LENGTH = 10000

TREES_PATH = Path.cwd() / "TRUE_TREES"
time_winners = []
measured_times: list[list[float]] = []



def refresh_sim(tree_file, indel_rate):
    sim_config = SimConfiguration(original_sequence_length=ROOT_SEQUENCE_LENGTH, indel_length_alpha=2.0,
                                indel_truncated_length=50,
                                rate_ins=indel_rate[0], rate_del=indel_rate[1],
                                deletion_extra_edge_length=50,
                                seed=420)
    sim = Simulation(input_tree=str(tree_file), config=sim_config)

    return sim

def run_simulation(sim: Simulation, type: str):
    print(type)
    if type == "list":
        sim.msa_from_blocklist()
    if type == "tree":
        sim.msa_from_blocktree()

    sim.msa.compute_msa()

branches_dfs = []


rates = [i*0.01 for i in range(1,10)]
rates = list(zip(rates, rates[::-1]))

#TODO: Change here to allow this benchmark as a script and take care of the 
# changes that were made to the MSA construction.
for RATE_MULTIPLIER in [1,5,10]:#[0.5, 1,2,4,8,16]:
    print("RATE_MULTIPLIER", RATE_MULTIPLIER)
    for indel_rate in rates:
        factor_compare = {"branch_scale": [], 
                    "blocklist_time": [], "blocktree_time": [],
                    "true": [], "factor": [],
                    "insertion_rate": [indel_rate[0]], "deletion_rate": [indel_rate[1]]
                    }
        scaled_tree_path = scale_tree("benchmark/scaled_trees/test_tree.txt", RATE_MULTIPLIER)
        sim = refresh_sim(scaled_tree_path, indel_rate)
        factor_compare["true"].append(len(sim.sim_nodes[1].list_of_events))
        factor_compare["factor"].append(sim.sim_nodes[1].hybrid_factor)
        for simulator_type in ["list", "tree"]:
            simulation_time = timeit.timeit(lambda: run_simulation(sim, simulator_type), number=3)/3.0
            factor_compare[f"block{simulator_type}_time"].append(simulation_time)

        factor_compare["branch_scale"].append(RATE_MULTIPLIER)
        comparison_df = pd.DataFrame(factor_compare)
        branches_dfs.append(comparison_df)

        print(indel_rate)


branches_df = pd.concat(branches_dfs)
# branches_df[branches_df["blocklist_time"] > branches_df["blocktree_time"]]
branches_df = branches_df.reset_index(drop=True)


branches_df.to_csv("benchmark/assets/data/list_vs_tree_indel_only.csv")

