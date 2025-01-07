import sys
import pathlib
import random
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd()))



from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree
from classes.seq_node_naive import SequenceNodeNaive
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
    

def test_insertion_at_boundary_event_list():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=11)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=3))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=5))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=1, place=10))
    res = new_organism.get_dto()
    print(res)
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 3, inserted len: 7',
            'predecessor index: 4, #copied sites: 7, inserted len: 0'],
        'length': 17}


def test_insertion_at_boundary_event_tree():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=11)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=3))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=5))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=1, place=10))
    res = new_organism.get_clean_dto()
    print(res)
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 3, inserted len: 7',
            'predecessor index: 4, #copied sites: 7, inserted len: 0'],
        'length': 17}
    
def event_creator(previous_length: int) -> tuple[IndelEvent, int]:
    event_type = random.choice([True, False])
    if previous_length == 0 :
        event_type = True

    event_size = random.randint(1,20)

    event_position = random.randint(0, previous_length -1 + event_type)

    if not event_type and (event_size + event_position > previous_length): # handle deletion bleeding out
        event_size = previous_length - event_position

    if event_type:
        previous_length += event_size
    else:
        previous_length -= event_size

    final_event = IndelEvent(event_type, event_position, event_size)
    return final_event, previous_length

def ordinal(n: int):
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    return str(n) + suffix

def test_random_events_tree_vs_list():
    seed = random.randint(1, 2147483647)
    # seed = 841039037
    random.seed(seed)
    current_sequence_length = 100
    blocklist = SequenceNodeAsList(0, current_sequence_length)
    blocktree = SequenceNodeAsTree(0, current_sequence_length)
    print(f"the seed for this run was: {seed}")
    for event_number in range(10000):
        if event_number == 181:
            stop = True
        if current_sequence_length == 0:
            break
        current_event, current_sequence_length = event_creator(current_sequence_length)
        print(f"the {ordinal(event_number)} event is: ", current_event)
        blocklist.calculate_event(current_event)
        blocktree.calculate_event(current_event)
        if blocklist.get_dto() != blocktree.get_clean_dto():
            print(blocklist.blocks_iterator())
            print(blocktree.blocks_iterator())
            print(f"the seed for this run was: {seed}")

        assert blocklist.get_dto() == blocktree.get_clean_dto()
    print("Sequence length is 0, halting\n")
    assert True



test_random_events_tree_vs_list()
