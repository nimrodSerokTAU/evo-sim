from pathlib import Path
from ete3 import Tree

from classes.sim_config import SimConfiguration
from classes.sim_node import SimulatedNode


class Simulation:
    tree: Tree
    sim_nodes: list[SimulatedNode]
    config: SimConfiguration

    def __init__(self, input_tree_path: Path):
        self.newick = '(A:1,(B:1,(E:1,D:1):0.5):0.5);'  # TODO: Nimrod, read the file from the path and get it from there
        self.tree = Tree(self.newick)
        # for node in self.tree.get_tree_root().get_root_children():
        #     self.sim_nodes.append(SimulatedNode(node.id: node.dist(), config, node.father_seq_length))  # TODO: Elya, here get it from the tree package
