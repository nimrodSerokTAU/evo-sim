from pathlib import Path
from ete3 import Tree

from classes.sim_node import SimulatedNode


class Simulation:
    tree: Tree
    sim_nodes: list[SimulatedNode]

    # TODO: here we should use a package that consider the internal nodes as well
    def __init__(self, input_tree_path: Path):
        self.newick = '(A:1,(B:1,(E:1,D:1):0.5):0.5);'
        self.tree = Tree(self.newick)
        for node in self.tree.iter_prepostorder():
            pass
