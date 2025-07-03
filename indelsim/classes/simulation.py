from pathlib import Path
from ete3 import Tree, TreeNode


from indelsim.classes.sim_config import SimConfiguration
from indelsim.classes.sim_node import SimulatedNode
from indelsim.classes.super_sequence import SuperSequence
from indelsim.classes.sequence import Sequence
from indelsim.classes.msa import Msa
from indelsim.utils import calc_msa_from_naive_nodes
from indelsim.classes.seq_node_as_list import SequenceNodeAsList
from indelsim.classes.seq_node_as_tree import SequenceNodeAsTree
from indelsim.classes.seq_node_naive import SequenceNodeNaive

class Simulation:
    tree: Tree
    sim_nodes: list[SimulatedNode]
    config: SimConfiguration
    nodes_to_align: set[int]
    msa: Msa | str

    def __init__(self, input_tree: Path|str, config: SimConfiguration):
        self.tree = Tree(input_tree)
        self.config = config
        self.nodes_to_align = set()
        self.nodes_to_align.add(0)
        self.sim_nodes = [None]
        node: TreeNode = None
        for idx, node in enumerate(self.tree.traverse("preorder")):
            node.add_features(id=idx)  # Assigning an ID based on index
            if node.id == 0:
                node.add_features(sequence_length=self.config.original_sequence_length)  # Assigning root length to root node
                continue # nothing to simulated in root node
            # if node.is_leaf():
            self.nodes_to_align.add(node.id)

            simulatedNode = SimulatedNode(node.id, node.up.id, node.dist, self.config, node.up.sequence_length)
            node.add_features(sequence_length=simulatedNode.length_of_sequence_after_events)


            self.sim_nodes.append(simulatedNode)


    def msa_from_blocklist(self):
        super_seq = SuperSequence(self.sim_nodes[1].length_of_sequence_before, len(self.nodes_to_align))
        sequences = []
        parent_seq = Sequence(super_seq, True, 0)
        parent_seq.init_root_seq()
        sequences.append(parent_seq)
        sequences_to_save = []
        for node in self.sim_nodes[1:]:
            # print("node id:", node.id)
            node.seq_node_as_list = SequenceNodeAsList(node.id, node.length_of_sequence_before)

            for event in node.list_of_events:
                node.seq_node_as_list.calculate_event(event)
            # print("done with events!")
            current_seq = Sequence(super_seq, node.id in self.nodes_to_align, node.id)
            blocks = node.seq_node_as_list.blocks_iterator()
            current_seq.generate_sequence(blocks, sequences[node.parent_id])
            # print("parent length",len(sequences[node.parent_id]._sequence))
            # print("child length", len(current_seq._sequence))

            sequences.append(current_seq)
            if node.id in self.nodes_to_align:
                sequences_to_save.append(current_seq)

        self.msa = Msa(super_seq)
        self.msa.compute_msa(sequences_to_save)
    
    def msa_from_blocktree(self):
        super_seq = SuperSequence(self.sim_nodes[1].length_of_sequence_before, len(self.nodes_to_align))
        parent_seq = Sequence(super_seq, True, 0)
        parent_seq.init_root_seq()
        sequences = [parent_seq]

        sequences_to_save = []
        for node in self.sim_nodes[1:]:
            node.seq_node_as_list = SequenceNodeAsTree(node.id, node.length_of_sequence_before)
            for event in node.list_of_events:
                node.seq_node_as_list.calculate_event(event)

            current_seq = Sequence(super_seq, node.id in self.nodes_to_align, node.id)
            blocks = node.seq_node_as_list.blocks_iterator()
            current_seq.generate_sequence(blocks, sequences[node.parent_id])

            sequences.append(current_seq)
            if node.id in self.nodes_to_align:
                sequences_to_save.append(current_seq)

        self.msa = Msa(super_seq)
        self.msa.compute_msa(sequences_to_save)

    def msa_from_naive(self):
        original_sequence_length: int = self.config.original_sequence_length

        root = SequenceNodeNaive(seq_id=0, original_sequence=[i for i in range(original_sequence_length)])
        sequences = [root.seq]

        sequences_to_save = []
        ids_to_save = []
        parent_ids_to_save = [-1]
        for node in self.sim_nodes[1:]:
            node.seq_node_naive = SequenceNodeNaive(node.id, sequences[node.parent_id])
        
            for event in node.list_of_events:
                node.seq_node_naive.calculate_event(event)
            current_seq = node.seq_node_naive.seq
            sequences.append(current_seq)

            if node.id in self.nodes_to_align:
                sequences_to_save.append(current_seq)
                ids_to_save.append(node.id)
                parent_ids_to_save.append(node.parent_id)

        msa: list[list[int]] = calc_msa_from_naive_nodes(sequences, parent_ids_to_save)

        msa_str = ""
        for id, seq in zip(ids_to_save, msa[1:]):
            seq_str = f">{id}\n" + "".join(["X" if (site != -1) else "-" for site in seq])
            msa_str += f"{seq_str}\n"
        
        self.msa = msa_str

    def msa_from_hybrid(self):
        super_seq = SuperSequence(self.sim_nodes[1].length_of_sequence_before, len(self.nodes_to_align))
        parent_seq = Sequence(super_seq, True, 0)
        parent_seq.init_root_seq()
        sequences = [parent_seq]

        sequences_to_save = []
        for node in self.sim_nodes[1:]:
            if node.hybrid_switch:
                node.seq_node_as_list = SequenceNodeAsTree(node.id, node.length_of_sequence_before)
            else:
                node.seq_node_as_list = SequenceNodeAsList(node.id, node.length_of_sequence_before)
 
            for event in node.list_of_events:
                node.seq_node_as_list.calculate_event(event)

            current_seq = Sequence(super_seq, node.id in self.nodes_to_align, node.id)
            blocks = node.seq_node_as_list.blocks_iterator()
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
    

