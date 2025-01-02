import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd()))

import utils
from classes.avl_tree import AVLTree
from classes.block import Block
from classes.indel_event import IndelEvent
from classes.seq_node_as_list import SequenceNodeAsList
from classes.seq_node_as_tree import SequenceNodeAsTree
from classes.seq_node_naive import SequenceNodeNaive
from classes.sim_config import SimConfiguration
from classes.super_sequence import SuperSequence
from classes.sequence import Sequence
from classes.msa import Msa

basic_config: SimConfiguration = SimConfiguration(
    original_sequence_length=100, indel_length_alpha=1.7, indel_truncated_length=50, deletion_extra_edge_length=49,
    rate_ins=0.01, rate_del=0.01, seed=1)


def test_insertion_including_inside_copied_and_at_end():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=10, place=-3))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 7',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 126}


def test_insertion_before_start_not_affecting():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=3, place=42))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 12',
            'predecessor index: 0, #copied sites: 30, inserted len: 8',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 120}


def test_deletion_at_start_of_copied():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 3, #copied sites: 27, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 102}


def test_deletion_at_start_insertion_only():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=0))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: -1, #copied sites: 0, inserted len: 9',
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 114}


def test_deletion_at_end_of_copied():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 104}


def test_deletion_after_end_of_copied():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=105))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 105}


def test_deletion_at_end_of_inserted():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
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


def test_list_case_el_a():
    new_organism = SequenceNodeAsList(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 12, inserted len: 2',
            'predecessor index: 12, #copied sites: 18, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 0',
            'predecessor index: 47, #copied sites: 53, inserted len: 0'],
        'length': 95}

################################################### avl #############################################


def test_tree_search_case_a():
    tree = AVLTree(Block(30, 5, 0))
    tree.insert_block(Block(4, 8, 2))
    tree.insert_block(Block(0, 1, 0))
    tree.insert_block(Block(12, 18, 5))
    tree.insert_block(Block(47, 48, 20))
    tree.insert_block(Block(95, 6, 0))
    node, seq_length_with_block_0 = tree.search(tree.root, 0)
    node, seq_length_with_block_9 = tree.search(tree.root, 9)
    ode, seq_length_with_block_20 = tree.search(tree.root, 20)
    node, seq_length_with_block_37 = tree.search(tree.root, 37)
    node, seq_length_with_block_40 = tree.search(tree.root, 40)
    node, seq_length_with_block_108 = tree.search(tree.root, 108)
    res = {'a': seq_length_with_block_0, 'b': seq_length_with_block_9, 'c': seq_length_with_block_20,
           'd': seq_length_with_block_37, 'e': seq_length_with_block_40, 'f': seq_length_with_block_108}
    assert res == {'a': 0, 'b': 8, 'c': 9, 'd': 3, 'e': 1, 'f': 1}


def test_tree_search_case_b():
    tree = AVLTree(Block(30, 5, 7))
    tree.insert_block(Block(0, 30, 5))
    tree.insert_block(Block(36, 5, 1))
    tree.insert_block(Block(35, 1, 1))
    tree.insert_block(Block(41, 59, 0))
    node, seq_length_with_block_1 = tree.search(tree.root, 1)
    node, seq_length_with_block_37 = tree.search(tree.root, 37)
    ode, seq_length_with_block_48 = tree.search(tree.root, 48)
    node, seq_length_with_block_51 = tree.search(tree.root, 51)
    node, seq_length_with_block_110 = tree.search(tree.root, 110)
    res = {'a': seq_length_with_block_1, 'b': seq_length_with_block_37, 'c': seq_length_with_block_48,
           'd': seq_length_with_block_51, 'e': seq_length_with_block_110}
    assert res == {'a': 1, 'b': 2, 'c': 1, 'd': 2, 'e': 55}


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
    node_2, seq_length_with_block_1 = tree.search(tree.root, 1)
    node_1, seq_length_with_block_11 = tree.search(tree.root, 11)
    node_0, seq_length_with_block_32 = tree.search(tree.root, 32)
    node_3, seq_length_with_block_35 = tree.search(tree.root, 35)
    node_5, seq_length_with_block_39 = tree.search(tree.root, 39)
    node_6, seq_length_with_block_44 = tree.search(tree.root, 44)
    node_4, seq_length_with_block_50 = tree.search(tree.root, 50)
    node_7, seq_length_with_block_54 = tree.search(tree.root, 54)
    node_8, seq_length_with_block_67 = tree.search(tree.root, 67)
    node_12, seq_length_with_block_74 = tree.search(tree.root, 74)
    node_9, seq_length_with_block_77 = tree.search(tree.root, 77)
    node_10, seq_length_with_block_84 = tree.search(tree.root, 84)
    node_11, seq_length_with_block_94 = tree.search(tree.root, 94)
    res = {'a': seq_length_with_block_1, 'b': seq_length_with_block_11, 'c': seq_length_with_block_32,
           'd': seq_length_with_block_35, 'e': seq_length_with_block_39, 'f': seq_length_with_block_44,
           'g': seq_length_with_block_50, 'h': seq_length_with_block_54, 'i': seq_length_with_block_67,
           'j': seq_length_with_block_74, 'k': seq_length_with_block_77, 'l': seq_length_with_block_84,
           'm': seq_length_with_block_94}
    assert res == {'a': 1, 'b': 9, 'c': 10, 'd': 1, 'e': 2, 'f': 4, 'g': 3, 'h': 1, 'i': 1, 'j': 1, 'k': 1, 'l': 1, 'm': 1}


def test_avl_insertion_elya_a():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    res = new_organism.get_dto()
    assert res == {
        'blocks': [
            'id: 0, predecessor index: 0, #copied sites: 30, inserted len: 5, length_under_including: 105',
            'id: 1, predecessor index: 30, #copied sites: 70, inserted len: 0, length_under_including: 70'],
        'length': 105}


def test_avl_deletion_elya_b():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 3, #copied sites: 27, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 102}


def test_avl_deletion_at_start_insertion_only():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 104}


def test_avl_deletion_after_end_of_copied():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=3, place=105))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 105}


def test_avl_deletion_at_end_of_inserted():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
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


def test_avl_case_el_a():
    new_organism = SequenceNodeAsTree(seq_id=0, original_sequence_length=100)
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=False, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=12))
    res = new_organism.get_clean_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 12, inserted len: 2',
            'predecessor index: 12, #copied sites: 18, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 0',
            'predecessor index: 47, #copied sites: 53, inserted len: 0'],
        'length': 95}


################################################### naive #############################################


def test_naive_insertion_including_inside_copied_and_at_end_regular_dto():
    original_sequence_length: int = 100
    new_organism = SequenceNodeNaive(seq_id=0, original_sequence=[i for i in range(original_sequence_length)])
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    res = new_organism.get_dto()
    assert res == {
        'seq': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                100, 101, 102, 103, 104,
                30, 31, 32, 33, 34,
                105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
                85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
                117, 118],
        'length': 119}


def test_naive_insertion_including_inside_copied_and_at_end_t():
    original_sequence_length: int = 100
    new_organism = SequenceNodeNaive(seq_id=0, original_sequence=[i for i in range(original_sequence_length)])
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=12, place=40))
    new_organism.calculate_event(IndelEvent(is_insertion=True, length=2, place=117))
    res = new_organism.get_block_dto_from_single_branch(original_sequence_length)
    assert res == {
         'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 5, inserted len: 12',
            'predecessor index: 35, #copied sites: 65, inserted len: 2'],
        'length': 119}


def test_naive_super_sequence():
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
    import re
    msa = []
    for seq in res:
        seq = re.sub("\s*\-1\,*\s*", "-", seq)
        seq = re.sub("\d+\,*\s*", "X", seq)
        msa.append(seq)
        # seq = seq.split(",")
        # seq
    print(">\n" + "\n>\n".join(msa))
    assert res == [
        ' -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 30, 31, 32, 33, 34, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, -1, -1, -1, -1, -1, -1',
        ' -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  2, -1, -1, -1, -1,  7,  8,  9,  -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,100,101,102,103, -1, -1, -1, -1, -1,104, 30, 31, 32, 33, 34,105,106,107,108,109,110,111,112,113,114,115,116, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,117,118,119,120,121,122',
        ' -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 100,101,102,103,104, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 30, 31, 32, 33, 34, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, -1, -1, -1, -1, -1, -1, -1',
        '131,132,133,134,135,136,137,138,139,  0,  1,  2, -1, -1, -1, -1,  7,  8,  9,  -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,100,101,102,103,123,124,125,126,127,104, 30, 31, 32, 33, 34,105,106,107,108,109,110,111,112,113,114,115,116, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,117,118,119,120,121,122',
    ]