import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd()))



from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree


# Test 1

def list_test_case_a():
    new_organism = SequenceNodeAsList(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]
    

def tree_test_case_a():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=95, place=0))

    res = new_organism.get_blocklist_str()

    print(res["blocks"])
    return res["blocks"]

def test_case_a():
    list_res = list_test_case_a()
    tree_res = tree_test_case_a()
    assert list_res == tree_res

# Test 2
def insertion_at_0_list():
    new_organism = SequenceNodeAsList(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=0))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]


def insertion_at_0_tree():
    new_organism = SequenceNodeAsTree(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=0))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]

def test_case_b():
    list_res = insertion_at_0_list()
    tree_res = insertion_at_0_tree()
    assert list_res == tree_res


def deletion_at_0_list():
    new_organism = SequenceNodeAsList(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=False, length=5, place=0))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]

def deletion_at_0_tree():
    new_organism = SequenceNodeAsTree(seq_id=0 ,original_sequence_length=100)

    new_organism.calculate_event(IndelEvent(is_insertion=False, length=5, place=0))
    # new_organism.calculate_event(IndelEvent(is_insertion=False, length=91, place=0))

    res = new_organism.get_blocklist_str()
    print(res["blocks"])
    return res["blocks"]

def test_case_c():
    list_res = deletion_at_0_list()
    tree_res = deletion_at_0_tree()
    assert list_res == tree_res


def test_zero_length_block_case_list():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=20, place=10))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 10, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 84}


def test_zero_length_block_case_tree():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=20, place=10))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 10, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 84}


def test_2_insertions_at_zero_list():
    new_organism = SequenceNodeAsList(seq_id=3, original_sequence_length=76)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    res = new_organism.get_dto()
    print(res)
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 24',
            'predecessor index: 0, #copied sites: 76, inserted len: 0'],
        'length': 100}


def test_2_insertions_at_zero_tree():
    new_organism = SequenceNodeAsTree(seq_id=3, original_sequence_length=76)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    res = new_organism.get_clean_dto()
    print(res)
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 24',
            'predecessor index: 0, #copied sites: 76, inserted len: 0'],
        'length': 100}