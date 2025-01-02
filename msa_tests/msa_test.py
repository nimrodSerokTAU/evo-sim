import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd()))
sys.path.append(str(pathlib.Path.cwd() / "classes"))

print(sys.path)

from super_sequence import SuperSequence
from block import Block
from sequence import Sequence
from msa import Msa
from indel_event import IndelEvent
from seq_node_as_list import SequenceNodeAsList

def root_to_leaf_alignment():
    super_seq = SuperSequence(20,3)

    root_seq = Sequence(super_seq, True, 0)
    root_seq.init_root_seq()

    child_seq = Sequence(super_seq, True, 1)
    blocks = [Block(-1,0,1),Block(0,10,5),Block(10,10,4)]
    child_seq.generate_sequence(blocks, root_seq)

    grandchild_seq = Sequence(super_seq, True, 2)
    blocks = [Block(0,5,2),Block(5,15,0),Block(23,7,1)]
    grandchild_seq.generate_sequence(blocks, child_seq)


    msa = Msa(super_seq)
    msa.compute_msa([root_seq, child_seq, grandchild_seq])


    return (msa.__repr__())


def test_root_to_leaf_alignment():
    true_msa = """>0
                  -XXXX--XXXXXX-----XXXXXXXXXX-----
                  >1
                  XXXXX--XXXXXXXXXXXXXXXXXXXXXXXXX-
                  >2
                  XXXXXXXXXXXXXXXXXXXXXX---XXXXXXXX
                  """
    true_msa = true_msa.replace(" ","")
    assert true_msa == root_to_leaf_alignment()

def full_tree_alignment():
    original_sequence_length: int = 100
    # organism_a = SequenceNodeAsTree(seq_id=0, original_sequence_length=original_sequence_length)
    organism_b = SequenceNodeAsList(seq_id=1, original_sequence_length=original_sequence_length)
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=4, place=119))
    organism_b.calculate_event(IndelEvent(is_insertion=False, length=4, place=3))

    organism_c = SequenceNodeAsList(seq_id=2, original_sequence_length=original_sequence_length)
    organism_c.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    organism_c.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    organism_c.calculate_event(IndelEvent(is_insertion=False, length=20, place=10))

    organism_d = SequenceNodeAsList(seq_id=3, original_sequence_length=organism_b.get_length())
    organism_d.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    organism_d.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    organism_d.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))

    super_seq = SuperSequence(original_sequence_length, 4)
    root_seq = Sequence(super_seq, is_save_seq=True, node_id=0)
    root_seq.init_root_seq()

    b_seq = Sequence(super_seq, is_save_seq=True, node_id=organism_b.id)
    print(organism_b.blocks_iterator())
    b_seq.generate_sequence(organism_b.blocks_iterator(), root_seq)

    c_seq = Sequence(super_seq, is_save_seq=True, node_id=organism_c.id)
    print(organism_c.blocks_iterator())
    c_seq.generate_sequence(organism_c.blocks_iterator(), root_seq)

    d_seq = Sequence(super_seq, is_save_seq=True, node_id=organism_d.id)
    print(organism_d.blocks_iterator())
    d_seq.generate_sequence(organism_d.blocks_iterator(), b_seq)

    msa = Msa(super_seq)
    msa.compute_msa([root_seq, b_seq , d_seq,c_seq])

    print(msa)


# def test_full_tree_alignment():
    
#     true_msa = true_msa.replace(" ","")
#     assert true_msa == root_to_leaf_alignment()


full_tree_alignment()