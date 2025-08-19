from __future__ import annotations

from llist import sllistnode #https://ajakubek.github.io/python-llist/index.html#llist.sllistnode

from indelsim.classes.super_sequence import SuperSequence
from indelsim.classes.block import Block



class Sequence:
    _super_seq: SuperSequence
    _is_save_sequence: bool
    _node_id: int
    _sequence: list[sllistnode]
    _number_of_children: int

    def __init__(self, super_seq: SuperSequence, is_save_seq: bool, node_id: int, number_of_children: int):
        self._super_seq = super_seq
        self._is_save_sequence = is_save_seq
        self._node_id = node_id
        self._number_of_children = number_of_children
        self._sequence = []

    def init_root_seq(self) -> None:
        super_seq_iterator = self._super_seq.get_iterator()

        for node in super_seq_iterator:
            self._sequence.append(node)
            if self._is_save_sequence:
                self._super_seq.reference_position(node)
    
    def generate_sequence(self, blocks: Iterator[Block], parent_seq: Sequence):
        """
        Generate a sequence based on a blocklist and parent sequence.
        
        Args:
            blocks: iterator of blocks containing position, length, and insertion info
            parent_seq: Parent sequence object containing the base sequence
        """
        
        #insert anchor site
        self._sequence.append(parent_seq._sequence[0])
        for block in blocks:
            self.apply_block(block, parent_seq)

    def apply_block(self, block: Block, parent_seq: Sequence):
        random_pos = self._super_seq.get_num_inserted_positions()

        position = block.index_in_predecessor
        length = block.copy_sites_count
        insertion = block.inserted_seq_count
        if position == -1:
            length = 0
        if length == 0 and insertion == 0:
            return
        if position == 0 and length == 0:
            position = -1
        # Copy parent sequence elements
        i = 0
        for i in range(length):
            if self._is_save_sequence:
                self._super_seq.reference_position(parent_seq._sequence[position + i + 1])
            self._sequence.append(parent_seq._sequence[position + i + 1])
        
        # Get iterator position
        super_seq_iterator = parent_seq._sequence[position + i + 1]
        if position > 0 and length == 0:
            super_seq_iterator = parent_seq._sequence[position + i]
        if position == -1:
            super_seq_iterator = parent_seq._sequence[0]
        if insertion == 0:
            return

        # Handle insertions
        for i in range(insertion):
            super_seq_iterator = self._super_seq.insert_item_at_position(
                super_seq_iterator, 
                random_pos,
                self._is_save_sequence
            )
            self._sequence.append(super_seq_iterator)
            # super_seq_iterator = super_seq_iterator.next
            random_pos = self._super_seq.increment_num_inserted_positions()
    
        if self._is_save_sequence:
            self._super_seq.increment_leaf_num()



    def get_super_sequence(self) -> SuperSequence:
        return self._super_seq
    
    def get_ref_to_super_sequence(self, pos) -> sllistnode:
        return self._sequence[pos]
    
    def get_sequence_node_id(self) -> int:
        return self._node_id

    def __repr__(self):
        return "·".join([str(node()['position']) for node in self._sequence])
    
    def __len__(self) -> int:
        return len(self._sequence)
    
    def __getitem__(self, index) -> sllistnode:
        return self._sequence[index]