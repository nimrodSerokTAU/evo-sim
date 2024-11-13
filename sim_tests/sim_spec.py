from classes.indel_event import IndelEvent
from classes.seq_node import SequenceNode
from classes.sim_config import SimConfiguration

basic_config: SimConfiguration = SimConfiguration(
    original_sequence=None, original_sequence_length=100, substitutions_per_site_rate=1, indel_length_alpha=1.7,
    indel_truncated_length=50, is_use_nucleotides=True, indel_per_sub_ratio=0.01)


def test_insertion_including_inside_copied_and_at_end():
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=0))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 3, #copied sites: 27, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 102}


def test_deletion_at_start_insertion_only():
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=104))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 69, inserted len: 0'],
        'length': 104}


def test_deletion_after_end_of_copied():
    new_node = SequenceNode(original_sequence_length=100)
    new_node.calculate_event(IndelEvent(is_insertion=True, length=5, place=30))
    new_node.calculate_event(IndelEvent(is_insertion=False, length=3, place=105))
    res = new_node.get_dto()
    assert res == {
        'blocks': [
            'predecessor index: 0, #copied sites: 30, inserted len: 5',
            'predecessor index: 30, #copied sites: 70, inserted len: 0'],
        'length': 105}


def test_deletion_at_end_of_inserted():
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
    new_node = SequenceNode(original_sequence_length=100)
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
