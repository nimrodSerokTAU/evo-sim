from classes.avl_node import AVLNode
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
    new_organism = SequenceNodeAsList(original_sequence_length=100)
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


def test_deletion_case_start_on_copy_start_and_contained():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=35))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 33, #copied sites: 2, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 114}


def test_deletion_case_start_on_copy_mid_and_contained():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=6, place=10))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 10, inserted len: 0',
            'predecessor index: 16, #copied sites: 14, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 111}


def test_deletion_case_start_on_copy_mid_and_not_contained():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=9, place=37))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 6',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 108}


def test_deletion_case_start_on_copy_mid_and_continue_to_next_block():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=37))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 0',
            'predecessor index: 37, #copied sites: 63, inserted len: 0'],
        'length': 100}


def test_deletion_case_start_on_insert_contained():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=4, place=43))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 8',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 113}


def test_deletion_case_start_on_insert_continue_to_next():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=19, place=42))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 2',
            'predecessor index: 44, #copied sites: 56, inserted len: 0'],
        'length': 98}


def test_insertion_inside_insertion():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=4, place=41))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=48))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=54))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 7',
            'predecessor index: 35, #copied sites: 1, inserted len: 1',
            'predecessor index: 36, #copied sites: 5, inserted len: 1',
            'predecessor index: 41, #copied sites: 59, inserted len: 0'],
        'length': 114}


def test_deletion_case_of_deleting_several_blocks():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=4, place=41))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=48))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=54))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=21, place=37))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 0',
            'predecessor index: 44, #copied sites: 56, inserted len: 0'],
        'length': 93}


def test_deletion_case_on_copied_more_than_end():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=111))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 59, inserted len: 0'],
        'length': 111}


def test_deletion_case_on_inserted_more_than_end():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=8, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=120))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 3'],
        'length': 120}


def test_deletion_before_start_not_affecting():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=-5))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_deletion_before_start_affecting():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=11, place=-5))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 6, #copied sites: 24, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 111}


def test_deletion_after_end_not_affecting():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=127))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_before_start():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=10, place=-3))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 0, inserted len: 7',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 126}


def test_insertion_before_start_not_affecting():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=-4))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_after_end_not_affecting():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=119))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_at_start():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=42))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 0, inserted len: 12',
            'predecessor index: 0, #copied sites: 30, inserted len: 8',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 120}


def test_deletion_at_start_of_copied():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 3, #copied sites: 27, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 102}


def test_deletion_at_start_insertion_only():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 0, inserted len: 9',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 114}


def test_deletion_at_end_of_copied():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 104}


def test_deletion_after_end_of_copied():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=105))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 105}


def test_deletion_at_end_of_inserted():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=2, place=118))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 1'],
        'length': 118}


def test_deletion_after_end_of_inserted():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=2, place=119))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}


def test_deletion_at_end_of_copied_plus_inserted():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=2, place=117))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_insertion_at_end_of_inserted():
    new_organism = SequenceNodeAsList(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=4, place=119))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 6'],
        'length': 123}


################################################### avl #############################################


def test_tree_search_case_a():
    tree = AVLTree(Block(30, 5, 0))
    tree.insert_block(Block(4, 8, 2))
    tree.insert_block(Block(0, 1, 0))
    tree.insert_block(Block(12, 18, 5))
    tree.insert_block(Block(47, 48, 20))
    tree.insert_block(Block(95, 6, 0))
    node, seq_length_with_block_0 = tree.search_for_insert(tree.root, 0, 0)
    node, seq_length_with_block_9 = tree.search_for_insert(tree.root, 9, 0)
    ode, seq_length_with_block_20 = tree.search_for_insert(tree.root, 20, 0)
    node, seq_length_with_block_37 = tree.search_for_insert(tree.root, 37, 0)
    node, seq_length_with_block_40 = tree.search_for_insert(tree.root, 40, 0)
    node, seq_length_with_block_108 = tree.search_for_insert(tree.root, 108, 0)
    res = {'a': seq_length_with_block_0, 'b': seq_length_with_block_9, 'c': seq_length_with_block_20,
           'd': seq_length_with_block_37, 'e': seq_length_with_block_40, 'f': seq_length_with_block_108}
    assert res == {'a': 1, 'b': 11, 'c': 34, 'd': 39, 'e': 107, 'f': 113}


def test_tree_search_case_b():
    tree = AVLTree(Block(30, 5, 7))
    tree.insert_block(Block(0, 30, 5))
    tree.insert_block(Block(36, 5, 1))
    tree.insert_block(Block(35, 1, 1))
    tree.insert_block(Block(41, 59, 0))
    node, seq_length_with_block_1 = tree.search_for_insert(tree.root, 1, 0)
    node, seq_length_with_block_37 = tree.search_for_insert(tree.root, 37, 0)
    ode, seq_length_with_block_48 = tree.search_for_insert(tree.root, 48, 0)
    node, seq_length_with_block_51 = tree.search_for_insert(tree.root, 51, 0)
    node, seq_length_with_block_110 = tree.search_for_insert(tree.root, 110, 0)
    res = {'a': seq_length_with_block_1, 'b': seq_length_with_block_37, 'c': seq_length_with_block_48,
           'd': seq_length_with_block_51, 'e': seq_length_with_block_110}
    assert res == {'a': 35, 'b': 47, 'c': 49, 'd': 55, 'e': 114}


def test_tree_search_case_c():
    tree = AVLTree(Block(30, 5, 7))
    tree.insert_block(Block(3, 15, 5))
    tree.insert_block(Block(0, 1, 1))
    tree.insert_block(Block(35, 3, 0))
    tree.insert_block(Block(58, 5, 1))
    tree.insert_block(Block(41, 1, 2))
    tree.insert_block(Block(54, 3, 4))
    tree.insert_block(Block(67, 11, 2))
    tree.insert_block(Block(80, 7, 0))
    tree.insert_block(Block(90, 4, 3))
    tree.insert_block(Block(94, 8, 2))
    tree.insert_block(Block(110, 3, 6))
    tree.insert_block(Block(89, 1, 2))
    node_2, seq_length_with_block_1 = tree.search_for_insert(tree.root, 1, 0)
    node_1, seq_length_with_block_11 = tree.search_for_insert(tree.root, 11, 0)
    node_0, seq_length_with_block_32 = tree.search_for_insert(tree.root, 32, 0)
    node_3, seq_length_with_block_35 = tree.search_for_insert(tree.root, 35, 0)
    node_5, seq_length_with_block_39 = tree.search_for_insert(tree.root, 39, 0)
    node_6, seq_length_with_block_44 = tree.search_for_insert(tree.root, 44, 0)
    node_4, seq_length_with_block_50 = tree.search_for_insert(tree.root, 50, 0)
    node_7, seq_length_with_block_54 = tree.search_for_insert(tree.root, 54, 0)
    node_8, seq_length_with_block_67 = tree.search_for_insert(tree.root, 67, 0)
    node_12, seq_length_with_block_74 = tree.search_for_insert(tree.root, 74, 0)
    node_9, seq_length_with_block_77 = tree.search_for_insert(tree.root, 77, 0)
    node_10, seq_length_with_block_84 = tree.search_for_insert(tree.root, 84, 0)
    node_11, seq_length_with_block_94 = tree.search_for_insert(tree.root, 94, 0)
    res = {'a': seq_length_with_block_1, 'b': seq_length_with_block_11, 'c': seq_length_with_block_32,
           'd': seq_length_with_block_35, 'e': seq_length_with_block_39, 'f': seq_length_with_block_44,
           'g': seq_length_with_block_50, 'h': seq_length_with_block_54, 'i': seq_length_with_block_67,
           'j': seq_length_with_block_74, 'k': seq_length_with_block_77, 'l': seq_length_with_block_84,
           'm': seq_length_with_block_94}
    # res = {'a': calc(node_2, True, 22), 'b': calc(node_1, True, 37), 'c': calc(node_0, False, 22),
    #        'd': calc(node_3, True, 66), 'e': calc(node_5, True, 47), 'f': calc(node_6, False, 37),
    #        'g': calc(node_4, False, 47), 'h': calc(node_7, False, 0), 'i': calc(node_8, True, 83),
    #        'j': calc(node_12, False, 73), 'k': calc(node_9, False, 66), 'l': calc(node_10, False, 83),
    #        'm': calc(node_11, False, 93)}
    assert res == {'a': 2, 'b': 22, 'c': 34, 'd': 37, 'e': 40, 'f': 47, 'g': 53, 'h': 66, 'i': 73, 'j': 76, 'k': 83,
                   'l': 93, 'm': 102}


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


def test_avl_deletion_elya_f():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=20, place=90))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=4, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=10))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 5, predecessor index: 4, #copied sites: 8, inserted len: 2, length_under_including: 16',
            'id: 6, predecessor index: 29, #copied sites: 1, inserted len: 5, length_under_including: 6',
            'id: 1, predecessor index: 30, #copied sites: 5, inserted len: 0, length_under_including: 94',
            'id: 2, predecessor index: 47, #copied sites: 48, inserted len: 20, length_under_including: 73',
            'id: 4, predecessor index: 95, #copied sites: 5, inserted len: 0, length_under_including: 5'],
        'length': 94}


def test_avl_insertion_including_inside_copied_and_at_end():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}
    
def test_avl_deletion_case_start_on_copy_start_and_contained():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=35))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 33, #copied sites: 2, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 114}


def test_avl_deletion_case_start_on_copy_mid_and_contained():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=6, place=10))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 10, inserted len: 0',
            'predecessor index: 16, #copied sites: 14, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 111}


def test_avl_deletion_case_start_on_copy_mid_and_not_contained():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=9, place=37))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 6',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 108}


def test_avl_deletion_case_start_on_copy_mid_and_continue_to_next_block():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=37))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 0',
            'predecessor index: 37, #copied sites: 63, inserted len: 0'],
        'length': 100}


def test_avl_deletion_case_start_on_insert_contained():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=4, place=43))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 8',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 113}


def test_avl_deletion_case_start_on_insert_continue_to_next():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=19, place=42))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 2',
            'predecessor index: 44, #copied sites: 56, inserted len: 0'],
        'length': 98}


def test_avl_insertion_inside_insertion():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=4, place=41))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=48))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=54))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 7',
            'predecessor index: 35, #copied sites: 1, inserted len: 1',
            'predecessor index: 36, #copied sites: 5, inserted len: 1',
            'predecessor index: 41, #copied sites: 59, inserted len: 0'],
        'length': 114}


def test_avl_deletion_case_of_deleting_several_blocks():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=4, place=41))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=48))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=1, place=54))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=21, place=37))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 2, inserted len: 0',
            'predecessor index: 44, #copied sites: 56, inserted len: 0'],
        'length': 93}


def test_avl_deletion_case_on_copied_more_than_end():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=111))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 59, inserted len: 0'],
        'length': 111}


def test_avl_deletion_case_on_inserted_more_than_end():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=8, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=17, place=120))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 3'],
        'length': 120}


def test_avl_deletion_before_start_not_affecting():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=-5))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_avl_deletion_before_start_affecting():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=11, place=-5))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 6, #copied sites: 24, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 111}


def test_avl_deletion_after_end_not_affecting():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=127))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_avl_insertion_before_start():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=10, place=-3))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 7',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 126}


def test_avl_insertion_before_start_not_affecting():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=-4))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_avl_insertion_after_end_not_affecting():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=119))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_avl_insertion_at_start():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=42))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 12',
            'predecessor index: 0, #copied sites: 30, inserted len: 8',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 120}


def test_avl_deletion_at_start_of_copied():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 3, #copied sites: 27, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 102}


def test_avl_deletion_at_start_insertion_only():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 9',  # TODO: consider overriding when exporting
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 114}


def test_avl_deletion_at_end_of_copied():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 104}


def test_avl_deletion_after_end_of_copied():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=105))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 105}


def test_avl_deletion_at_end_of_inserted():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=2, place=118))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 1'],
        'length': 118}


def test_avl_deletion_after_end_of_inserted():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=2, place=119))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}


def test_avl_deletion_at_end_of_copied_plus_inserted():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=2, place=117))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 0'],
        'length': 117}


def test_avl_insertion_at_end_of_inserted():
    new_organism = SequenceNodeAsTree(original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=4, place=119))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 6'],
        'length': 123}