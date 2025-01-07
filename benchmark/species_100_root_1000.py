
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

ROOT_SEQUENCE_LENGTH = 1000
# RATE_MULTIPLIER = 15

for RATE_MULTIPLIER in range(1,50):

    sim_config = SimConfiguration(original_sequence_length=ROOT_SEQUENCE_LENGTH, indel_length_alpha=1.5,
                                indel_truncated_length=5,
                                rate_ins=0.03 * RATE_MULTIPLIER, rate_del=0.09 * RATE_MULTIPLIER,
                                deletion_extra_edge_length=5, seed=42)

    sim = Simulation(input_tree="benchmark/normalbranches_nLeaves10.treefile", config=sim_config)
    print("events are ready!")

    blocktree_time = timeit.timeit(sim.msa_from_blocktree, number=3)
    print("blocktree_time","is", blocktree_time, "sec")

    blocklist_time = timeit.timeit(sim.msa_from_blocklist, number=3)
    print("blocklist_time","is", blocklist_time, "sec")

    if blocktree_time < blocklist_time:
        print(f"tree is faster than list once the rate multiplier is {RATE_MULTIPLIER}")
# blocklist_msa = sim.msa.msa_str_rep()
# with open("blocklist_msa.fasta", 'w') as f:
#     f.write(blocklist_msa)


# with open("blocktree_msa.fasta", 'w') as f:
#     f.write(sim.msa.msa_str_rep())

# naive_time = timeit.timeit(sim.msa_from_naive, number=10)
# print("naive_time","is", naive_time, "sec")

# with open("naive_msa.fasta", 'w') as f:
#     f.write(sim.msa)






