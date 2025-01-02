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

    def __init__(self, node_id: int, parent_id: int, branch_length: float, config: SimConfiguration, father_seq_length: int):

        self.id = node_id
        self.parent_id = parent_id

        self.branch_length = branch_length
        self.length_of_sequence_before = father_seq_length
        self.length_of_sequence_after_events = -1
        rnd.seed(config.random_seed)
        np.random.seed(config.random_seed)
        self.list_of_events = self.create_events(config, father_seq_length)
        config.random_seed += 1

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
        self.length_of_sequence_after_events = current_running_length
        return events
    
    def apply_events_with_tree(self):
        self.seq_node_as_tree = SequenceNodeAsTree(self.id, self.length_of_sequence_before)
        
        for event in self.list_of_events:
            self.seq_node_as_tree.calculate_event(event)
        

    
    def __repr__(self):
        return str(self.id) + "\n" + "\n".join(str(x) for x in self.list_of_events)



def insertion_event(config: SimConfiguration, current_running_length: int) -> IndelEvent:
    place: int = rnd.randint(0, current_running_length)
    insertion_size: int = calc_trunc_zipf(config.indel_length_alpha, config.indel_truncated_length)
    return IndelEvent(True, place, insertion_size)

def deletion_event(config: SimConfiguration, current_running_length: int) -> IndelEvent | None:
    start: int = -config.deletion_extra_edge_length
    place: int = rnd.randint(start, current_running_length - 1)
    deletion_size: int = calc_trunc_zipf(config.indel_length_alpha, config.indel_truncated_length)
    if place + deletion_size > 0:
        return IndelEvent(False, place, deletion_size)
    return None

