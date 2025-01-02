from pathlib import Path
import json
from ete3 import Tree, TreeNode, PhyloNode

import sys
print(Path.cwd())
print(sys.path)
sys.path.append(str(Path.cwd()))

from sim_config import SimConfiguration
from sim_node import SimulatedNode
from super_sequence import SuperSequence
from sequence import Sequence
from msa import Msa

from block import Block

class Simulation:
    tree: Tree
    sim_nodes: list[SimulatedNode]
    config: SimConfiguration
    nodes_to_align: set[int]
    msa: Msa

    def __init__(self, input_tree_path: Path):
        self.newick = '(A:1,(B:1,(E:1,D:1):0.5):0.5);'  # TODO: Nimrod, read the file from the path and get it from there
        self.tree = Tree(self.newick)
        self.config = SimConfiguration(original_sequence_length=50, indel_length_alpha=1.5,
                                       indel_truncated_length=5, rate_ins=0.1, rate_del=0.1,
                                       deletion_extra_edge_length=5, seed=4)
        self.nodes_to_align = set()
        self.sim_nodes = [None]
        node: TreeNode = None
        for idx, node in enumerate(self.tree.traverse("preorder")):
            node.add_features(id=idx)  # Assigning an ID based on index
            if node.id == 0:
                node.add_features(sequence_length=self.config.original_sequence_length)  # Assigning root length to root node
                continue # nothing to simulated in root node
            if node.is_leaf():
                self.nodes_to_align.add(node.id)

            simulatedNode = SimulatedNode(node.id, node.up.id, node.dist, self.config, node.up.sequence_length)
            simulatedNode.apply_events_with_tree()
            node.add_features(sequence_length=simulatedNode.length_of_sequence_after_events)


            self.sim_nodes.append(simulatedNode)


    def generate_alignement_from_blocks(self):
        super_seq = SuperSequence(self.sim_nodes[1].length_of_sequence_before, 7)
        sequences = []
        parent_seq = Sequence(super_seq, False, 0)
        parent_seq.init_root_seq()
        sequences.append(parent_seq)

        sequences_to_save = []
        for node in self.sim_nodes[1:]:
            current_seq = Sequence(super_seq, node.id in self.nodes_to_align, node.id)
            blocks = node.seq_node_as_tree.blocks_iterator()

            current_seq.generate_sequence(blocks, sequences[node.parent_id])
            sequences.append(current_seq)
            if node.id in self.nodes_to_align:
                sequences_to_save.append(current_seq)

        self.msa = Msa(super_seq)
        self.msa.compute_msa(sequences_to_save)


    def get_events(self):
        return self.sim_nodes
    
    def __repr__(self):
        return "\n$".join(f"{node}" for node in self.sim_nodes)

events = Simulation("")

events.generate_alignement_from_blocks()

print(events)
print(events.msa)