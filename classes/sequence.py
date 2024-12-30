from __future__ import annotations

from llist import sllistnode #https://ajakubek.github.io/python-llist/index.html#llist.sllistnode

from super_sequence import SuperSequence
from block import Block



class Sequence:
    _super_seq: SuperSequence
    _is_save_sequence: bool
    _node_id: int
    _sequence: list[sllistnode]

    def __init__(self, super_seq: SuperSequence, is_save_seq: bool, node_id: int):
        self._super_seq = super_seq
        self._is_save_sequence = is_save_seq
        self._node_id = node_id
        self._sequence = []

    def init_root_seq(self) -> None:
        super_seq_iterator = self._super_seq.get_iterator()

        for node in super_seq_iterator:
            print(node)
            self._sequence.append(node)
            if self._is_save_sequence:
                self._super_seq.reference_position(node)

    def generate_sequence(self, blocks: Iterator[Block], parent_seq: Sequence):
        """
        Generate a sequence based on a blocklist and parent sequence.
        
        Args:
            blocklist: List of blocks containing position, length, and insertion info
            parent_seq: Parent sequence object containing the base sequence
        """
        
        random_pos = self._super_seq.get_num_inserted_positions()
        
        for block in blocks:
            position = block.index_in_predecessor
            length = block.copy_sites_count
            insertion = block.inserted_seq_count
            
            # Special case for single element
            if position == 0 and length == 1 and insertion == 0:
                self._sequence.append(parent_seq._sequence[0])
                continue
                
            # Adjust position or length based on starting position
            if position != 0:
                position -= 1
            else:
                length -= 1
            # Copy sequence elements
            for i in range(length):
                if self._is_save_sequence:
                    self._super_seq.reference_position(parent_seq._sequence[position + i])
                self._sequence.append(parent_seq._sequence[position + i])
            
            # Get iterator position
            super_seq_iterator = parent_seq._sequence[position]
            if self._sequence:
                super_seq_iterator = parent_seq._sequence[position + length - 1]
                super_seq_iterator = super_seq_iterator.next if (super_seq_iterator.next is not None) else super_seq_iterator

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

    
super_seq = SuperSequence(50,3)

root_seq = Sequence(super_seq, True, 0)
root_seq.init_root_seq()
print(root_seq)

child_seq = Sequence(super_seq, False, 1)
blocks = [Block(0,10,5),Block(15,35,4)]
child_seq.generate_sequence(blocks, root_seq)

print(child_seq)
grandchild_seq = Sequence(super_seq, True, 1)
blocks = [Block(0,5,2),Block(5,15,0),Block(30,24,1)]
grandchild_seq.generate_sequence(blocks, child_seq)
print(grandchild_seq)
print(super_seq)