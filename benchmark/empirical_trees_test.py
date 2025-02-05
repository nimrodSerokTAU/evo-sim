
import timeit
from pathlib import Path


import sys
print(Path.cwd())
print(sys.path)
sys.path.append(str(Path.cwd()))
sys.path.append(str(Path.cwd() / "classes"))


from sim_config import SimConfiguration
from simulation import Simulation
from msa import Msa

def argmin(a):
    return min(range(len(a)), key=lambda x : a[x])

ROOT_SEQUENCE_LENGTH = 1000

TREES_PATH = Path("benchmark/TRUE_TREES")
# RATE_MULTIPLIER = 15
time_winners = []

for tree_path in list(TREES_PATH.iterdir())[:100]:#[0.5, 1,2,4,8,16]:

    sim_config = SimConfiguration(original_sequence_length=ROOT_SEQUENCE_LENGTH, indel_length_alpha=1.5,
                                indel_truncated_length=5,
                                rate_ins=0.03, rate_del=0.09,
                                deletion_extra_edge_length=5, seed=420)

    sim = Simulation(input_tree=str(tree_path), config=sim_config)
    print(f"Events are ready!")

    hybrid_time = timeit.timeit(sim.msa_from_hybrid, number=1)
    print("Hybrid time","is", hybrid_time, "sec")
    hybrid_msa = sim.msa.msa_str_rep()

    blocktree_time = timeit.timeit(sim.msa_from_blocktree, number=1)
    print("Blocktree time","is", blocktree_time, "sec")
    blocktree_msa = sim.msa.msa_str_rep()
    # with open("blocktree_msa.fasta", 'w') as f:
    #     f.write(blocktree_msa)

    blocklist_time = timeit.timeit(sim.msa_from_blocklist, number=1)
    print("Blocklist time","is", blocklist_time, "sec")
    blocklist_msa = sim.msa.msa_str_rep()
    # with open("blocklist_msa.fasta", 'w') as f:
    #     f.write(blocklist_msa)

    naive_time = timeit.timeit(sim.msa_from_naive, number=1)
    print("Naive_time","is", naive_time, "sec")
    # with open("naive_msa.fasta", 'w') as f:
    #     f.write(sim.msa)
    time_labels = ["hybrid_time", "blocktree_time", "blocklist_time", "naive_time"]

    time_measures = [hybrid_time, blocktree_time, blocklist_time, naive_time]

    min_time = argmin(time_measures)
    time_winners.append(min_time)

    if blocktree_time < blocklist_time:
        print(f"tree is faster than list on the tree {tree_path}")


for i in range(4):
    print(time_labels[i], time_winners.count(i))





