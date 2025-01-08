
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

for RATE_MULTIPLIER in range(1,16):#[0.5, 1,2,4,8,16]:

    sim_config = SimConfiguration(original_sequence_length=ROOT_SEQUENCE_LENGTH, indel_length_alpha=1.5,
                                indel_truncated_length=5,
                                rate_ins=0.03 * RATE_MULTIPLIER, rate_del=0.09 * RATE_MULTIPLIER,
                                deletion_extra_edge_length=5, seed=420)

    sim = Simulation(input_tree="benchmark/normalbranches_nLeaves100.treefile", config=sim_config)
    print(f"Events for MULTIPLIER {RATE_MULTIPLIER} are ready!")

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

    # if blocktree_time < blocklist_time:
    #     print(f"tree is faster than list once the rate multiplier is {RATE_MULTIPLIER}")








