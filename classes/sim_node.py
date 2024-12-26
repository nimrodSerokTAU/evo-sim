import numpy as np
import random as rnd

from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree
from classes.seq_node_naive import SequenceNodeNaive
from classes.sim_config import SimConfiguration
from utils import calc_trunc_zipf


class SimulatedNode:
    id: int
    branch_length: float
    list_of_events: list[IndelEvent]
    seq_node_naive: SequenceNodeNaive
    seq_node_as_list: SequenceNodeAsList
    seq_node_as_tree: SequenceNodeAsTree

    def __init__(self, node_id: int, branch_length: float, config: SimConfiguration, father_seq_length: int):

        self.id = node_id
        self.branch_length = branch_length
        self.list_of_events = self.create_events(config, father_seq_length)

    def create_events(self, config: SimConfiguration, father_seq_length: int) -> list[IndelEvent]:
        events: list[IndelEvent] = []
        current_time: float = 0
        current_running_length: int = father_seq_length
        while True:
            total_rate_across_entire_sequence = config.rate_ins * (current_running_length + 1) + config.rate_del * current_running_length
            event_time = np.random.exponential(1.0/ total_rate_across_entire_sequence)
            current_time += event_time
            if current_time > self.branch_length:
                break
            insertion_prob = config.rate_ins * (current_running_length + 1) / total_rate_across_entire_sequence
            is_insert = insertion_prob < rnd.uniform(0, 1)
            if is_insert:
                event: IndelEvent = insertion_event(config, current_running_length)
                events.append(event)
                current_running_length += event.length
            else:
                event: IndelEvent | None = deletion_event(config, current_running_length)
                if event is not None:
                    events.append(event)
                    current_running_length -= event.length
        return events


def insertion_event(config: SimConfiguration, current_running_length: int) -> IndelEvent:
    place: int = rnd.randint(0, current_running_length)
    insertion_size: int = calc_trunc_zipf(config.indel_length_alpha, config.indel_truncated_length)
    return IndelEvent(True, place, insertion_size)

def deletion_event(config: SimConfiguration, current_running_length: int) -> IndelEvent | None:
    start: int = -config.deletion_extra_edge_length
    place: int = rnd.randint(start, current_running_length - 1)
    deletion_size: int = calc_trunc_zipf(config.indel_length_alpha, config.indel_truncated_length)
    if place + deletion_size > 0:
        return IndelEvent(True, place, deletion_size)
    return None

