import random

import numpy as np

from simulation.classes.indel_event import IndelEvent
from simulation.classes.seq_node import SequenceNode
from simulation.classes.sim_config import SimConfiguration
from simulation.constants import NUCLEOTIDES
from utils import calc_inserted_seq, calc_inserted_seq_len_and_place


class SimulatedNode:
    indel_length_alpha: float
    indel_truncated_length: int
    list_of_options: list[str]
    seq_node: SequenceNode

    def __init__(self, config: SimConfiguration):
        self.is_use_nucleotides = config.is_use_nucleotides
        if config.is_use_nucleotides:
            self.list_of_options = NUCLEOTIDES
        else:
            self.list_of_options = ['A']
        self.indel_length_alpha = config.indel_length_alpha
        self.indel_truncated_length = config.indel_truncated_length
        self.seq_node = SequenceNode(config.original_sequence_length)

    # def create_original_sequence(self, seq: str | None, seq_len: str | None) -> list[str]:
    #     if seq is not None:
    #         return list(seq)
    #     else:
    #         if seq_len is None:
    #             seq_len = 1
    #         return self.create_rand_seq(seq_len)
    #
    # def create_rand_seq(self, seq_length: int) -> list[str]:
    #
    #     fair_prob: list[float] = [1 / len(self.list_of_options)] * len(self.list_of_options)
    #     seq: list[str] = np.random.choice(self.list_of_options, seq_length, p=fair_prob)
    #     return seq
    #
    # def calc_is_indel(self) -> bool:
    #     rand_res: float = np.random.default_rng().uniform()
    #     return rand_res < self.indel_per_sub_ratio
    #
    # def calc_efficient_gillespie(self):
    #     current_time: float = 0
    #     while current_time < self.new_seq_bl:
    #         current_time += np.random.default_rng().exponential(1.0 / len(self.new_seq), 1)
    #         if current_time < self.new_seq_bl:
    #             is_indel: bool = self.calc_is_indel()
    #             if is_indel:
    #                 self.create_indel_event()
    #             else:
    #                 c_inx: int = random.randint(0, len(self.new_seq) - 1)
    #                 self.new_seq[c_inx] = self.get_substitution(self.new_seq[c_inx])
    #
    # def get_substitution(self, current_code: str) -> str:
    #     options: list[str] = NUCLEOTIDES.copy()
    #     options.remove(current_code)
    #     inx: int = random.randint(0, 2)
    #     return options[inx]

    def create_insertion_event(self, optional_event_override: IndelEvent | None):
        if optional_event_override is None:
            event = calc_inserted_seq_len_and_place(self.seq_node.my_length, self.indel_length_alpha)
        else:
            event = optional_event_override
        self.seq_node.calculate_insertion_event(event)

    def create_deletion_event(self, optional_event_override: IndelEvent | None):
        if optional_event_override is None:
            event = calc_inserted_seq_len_and_place(self.seq_node.my_length, self.indel_length_alpha)
        else:
            event = optional_event_override
        self.seq_node.calculate_deletion_event(event)



    # def simulate_node_events(self, sim_config: SimConfiguration):
    #     pass


    # def task_1(alpha: float):
    #     seq: list[str] = ['A'] * 1000
    #     sim_iterations: int = 100000
    #     seq_lengths: list[int] = []
    #     for i in range(sim_iterations):
    #         s = sim_indels(seq=seq.copy(), alpha=alpha, indels_count=20)
    #         seq_lengths.append(len(s))
    #     avg: float = sum(seq_lengths) / sim_iterations
    #     stdev: float = statistics.stdev(seq_lengths)

