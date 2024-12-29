import sys, pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd()))

from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree



def test_list():
    new_organism = SequenceNodeAsList(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]
    

def test_tree():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=95, place=0))

    res = new_organism.get_block_list_str()

    print(res["blocks"])
    return res["blocks"]


assert test_list() == test_tree()
