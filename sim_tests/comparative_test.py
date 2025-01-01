import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd()))

from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree


# Test 1

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

    res = new_organism.get_blocklist_str()

    print(res["blocks"])
    return res["blocks"]


assert test_list() == test_tree()

# Test 2
def test_insertion_at_0_list():
    new_organism = SequenceNodeAsList(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=0))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]


def test_insertion_at_0_tree():
    new_organism = SequenceNodeAsTree(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=0))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]




assert test_insertion_at_0_list() == test_insertion_at_0_tree()
