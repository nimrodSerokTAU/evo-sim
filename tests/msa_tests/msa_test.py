import re

from indelsim import utils

from indelsim.classes.super_sequence import SuperSequence
from indelsim.classes.block import Block
from indelsim.classes.sequence import Sequence
from indelsim.classes.msa import Msa
from indelsim.classes.indel_event import IndelEvent
from indelsim.classes.seq_node_as_list import SequenceNodeAsList
from indelsim.classes.seq_node_naive import SequenceNodeNaive


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

def msa_from_blocks_4_nodes():
    original_sequence_length: int = 100

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
    b_seq.generate_sequence(organism_b.blocks_iterator(), root_seq)

    c_seq = Sequence(super_seq, is_save_seq=True, node_id=organism_c.id)
    c_seq.generate_sequence(organism_c.blocks_iterator(), root_seq)

    d_seq = Sequence(super_seq, is_save_seq=True, node_id=organism_d.id)
    d_seq.generate_sequence(organism_d.blocks_iterator(), b_seq)

    msa = Msa(super_seq)
    msa.compute_msa([root_seq, b_seq, c_seq, d_seq])

    return msa.msa_str_rep()


def msa_from_naive_4_nodes():
    original_sequence_length: int = 100
    organism_a = SequenceNodeNaive(seq_id=0, original_sequence=[i for i in range(original_sequence_length)])
    organism_b = SequenceNodeNaive(seq_id=1, original_sequence=organism_a.seq)
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    organism_b.calculate_event(IndelEvent(is_insertion=True, length=4, place=119))
    organism_b.calculate_event(IndelEvent(is_insertion=False, length=4, place=3))
    organism_c = SequenceNodeNaive(seq_id=2, original_sequence=organism_a.seq)
    organism_c.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    organism_c.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    organism_c.calculate_event(IndelEvent(is_insertion=False, length=20, place=10))
    organism_d = SequenceNodeNaive(seq_id=3, original_sequence=organism_b.seq)
    organism_d.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    organism_d.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    organism_d.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    msa: list[list[int]] = utils.calc_msa_from_naive_nodes([organism_a.seq, organism_b.seq, organism_c.seq,
                                                            organism_d.seq], [-1, 0, 0, 1])
    res: list[str] = utils.get_msa_as_str_list(msa, 3)

    orgnism_ids: list[int] = [organism_a.id, organism_b.id, organism_c.id, organism_d.id]
    msa = []
    for id,seq in zip(orgnism_ids , res):
        seq = re.sub(r"\s*\-1\,*\s*", "-", seq)
        seq = re.sub(r"\d+\,*\s*", "X", seq)
        msa.append(f">{id}\n{seq}")

    msa: str = "\n".join(msa) + "\n"
    return msa


def test_msas_4_nodes():
    blocks_msa = msa_from_blocks_4_nodes()
    naive_msa = msa_from_naive_4_nodes()
    assert blocks_msa == naive_msa
    

