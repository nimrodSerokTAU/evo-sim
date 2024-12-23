from classes.avl_tree import AVLTree
from classes.block import Block
from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree
from classes.sim_config import SimConfiguration

basic_config: SimConfiguration = SimConfiguration(
    original_sequence=None, original_sequence_length=100, substitutions_per_site_rate=1, indel_length_alpha=1.7,
    indel_truncated_length=50, is_use_nucleotides=True, indel_per_sub_ratio=0.01)


def test_insertion_including_inside_copied_and_at_end():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}


def test_deletion_case_start_on_copy_start_and_contained():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=35))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 33, #copied sites: 2, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 114}


def test_deletion_case_start_on_copy_mid_and_contained():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=6, place=10))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 10, inserted len: 0',
            'predecessor index: 16, #copied sites: 14, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 111}


def test_deletion_case_start_on_copy_mid_and_not_contained():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=9, place=37))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 6',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 108}


def test_deletion_case_start_on_copy_mid_and_continue_to_next_block():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=17, place=37))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 0',
            'predecessor index: 37, #copied sites: 63, inserted len: 0'],
        'length': 100}


def test_deletion_case_start_on_insert_contained():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=4, place=43))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 8',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 113}


def test_deletion_case_start_on_insert_continue_to_next():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=19, place=42))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 2',
            'predecessor index: 44, #copied sites: 56, inserted len: 0'],
        'length': 98}


def test_insertion_inside_insertion():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=3, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=4, place=41))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=1, place=48))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=1, place=54))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 7',
            'predecessor index: 35, #copied sites: 1, inserted len: 1',
            'predecessor index: 36, #copied sites: 5, inserted len: 1',
            'predecessor index: 41, #copied sites: 59, inserted len: 0'],
        'length': 114}


def test_deletion_case_of_deleting_several_blocks():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=3, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=4, place=41))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=1, place=48))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=1, place=54))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=21, place=37))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 0',
            'predecessor index: 44, #copied sites: 56, inserted len: 0'],
        'length': 93}


def test_deletion_case_on_copied_more_than_end():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=17, place=111))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 59, inserted len: 0'],
        'length': 111}


def test_deletion_case_on_inserted_more_than_end():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=8, place=117))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=17, place=120))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 3'],
        'length': 120}


def test_deletion_before_start_not_affecting():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=-5))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_deletion_before_start_affecting():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=11, place=-5))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 6, #copied sites: 24, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 111}


def test_deletion_after_end_not_affecting():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=127))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_before_start():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=10, place=-3))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 0, inserted len: 7',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 126}


def test_insertion_before_start_not_affecting():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=-4))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_after_end_not_affecting():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=119))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_at_start():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=3, place=42))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 0, inserted len: 12',
            'predecessor index: 0, #copied sites: 30, inserted len: 8',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 120}


def test_deletion_at_start_of_copied():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 3, #copied sites: 27, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 102}


def test_deletion_at_start_insertion_only():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 0, inserted len: 9',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 114}


def test_deletion_at_end_of_copied():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 104}


def test_deletion_after_end_of_copied():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=105))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 105}


def test_deletion_at_end_of_inserted():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=2, place=118))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 1'],
        'length': 118}


def test_deletion_after_end_of_inserted():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=2, place=119))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}


def test_deletion_at_end_of_copied_plus_inserted():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=2, place=117))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_at_end_of_inserted():
    new_node = SequenceNodeAsList(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_node.calculate_event(IndelEvent(is_insertion=True, length=4, place=119))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 6'],
        'length': 123}


############################################################################# avl #####################


def avl_insertion_including_inside_copied_and_at_end():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}


def test_tree_search():
    tree = AVLTree(Block(30, 5, 0))
    tree.insert_block(Block(4, 8, 2))
    tree.insert_block(Block(0, 1, 0))
    tree.insert_block(Block(12, 18, 5))
    tree.insert_block(Block(47, 48, 20))
    tree.insert_block(Block(95, 6, 0))
    node, seq_length_with_block = tree.search(tree.root, 37, 0)
    assert seq_length_with_block == 39

def test_avl_insertion_elya_a():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 0, predecessor index: 0, #copied sites: 30, inserted len: 5, length_under_including: 105',
            'id: 1, predecessor index: 30, #copied sites: 70, inserted len: 0, length_under_including: 70'],
        'length': 105}


def test_avl_deletion_elya_b():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 0, predecessor index: 0, #copied sites: 30, inserted len: 5, length_under_including: 35',
            'id: 1, predecessor index: 30, #copied sites: 5, inserted len: 0, length_under_including: 93',
            'id: 2, predecessor index: 47, #copied sites: 53, inserted len: 0, length_under_including: 53'],
        'length': 93}


def test_avl_deletion_elya_c():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 0, predecessor index: 0, #copied sites: 12, inserted len: 2, length_under_including: 37',
            'id: 3, predecessor index: 12, #copied sites: 18, inserted len: 5, length_under_including: 23',
            'id: 1, predecessor index: 30, #copied sites: 5, inserted len: 0, length_under_including: 95',
            'id: 2, predecessor index: 47, #copied sites: 53, inserted len: 0, length_under_including: 53'],
        'length': 95}


def test_avl_deletion_elya_d():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=20, place=90))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 0, predecessor index: 0, #copied sites: 12, inserted len: 2, length_under_including: 37',
            'id: 3, predecessor index: 12, #copied sites: 18, inserted len: 5, length_under_including: 23',
            'id: 1, predecessor index: 30, #copied sites: 5, inserted len: 0, length_under_including: 115',
            'id: 2, predecessor index: 47, #copied sites: 48, inserted len: 20, length_under_including: 73',
            'id: 4, predecessor index: 95, #copied sites: 5, inserted len: 0, length_under_including: 5'],
        'length': 115}


def test_avl_deletion_elya_e():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=20, place=90))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=4, place=0))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 5, predecessor index: 4, #copied sites: 8, inserted len: 2, length_under_including: 10',
            'id: 3, predecessor index: 12, #copied sites: 18, inserted len: 5, length_under_including: 33',
            'id: 1, predecessor index: 30, #copied sites: 5, inserted len: 0, length_under_including: 111',
            'id: 2, predecessor index: 47, #copied sites: 48, inserted len: 20, length_under_including: 73',
            'id: 4, predecessor index: 95, #copied sites: 5, inserted len: 0, length_under_including: 5'],
        'length': 111}