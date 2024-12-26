from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree
from classes.seq_node_naive import SequenceNodeNaive


class SimulatedNode:
    id: int
    branch_length: float
    list_of_events: list[IndelEvent]
    seq_node_naive: SequenceNodeNaive
    seq_node_as_list: SequenceNodeAsList
    seq_node_as_tree: SequenceNodeAsTree

    def __init__(self, node_id: int, branch_length: float):

        self.id = node_id
        self.branch_length = branch_length
        self.list_of_events = self.create_events()

    def create_events(self) -> list[IndelEvent]:
        current_running_length: int
        # TODO: here calculate events according to branch length and current_running_length
        pass

