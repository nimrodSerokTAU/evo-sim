import random as rnd

from utils import calc_trunc_zipf


class BlockListItem:
    index_in_predecessor: int
    copy_sites_count: int
    inserted_seq_count: int

    def __init__(self, index_in_predecessor: int = 0, copy_sites_count: int = 0, inserted_seq_count: int = 0):
        self.index_in_predecessor = index_in_predecessor
        self.copy_sites_count = copy_sites_count
        self.inserted_seq_count = inserted_seq_count

    def update_insert_count(self, event_length: int):
        self.inserted_seq_count = event_length

    def inc_insert_count(self, event_length: int):
        self.inserted_seq_count += event_length

    def update_copy_sites_count(self, copied_length: int):
        self.copy_sites_count = copied_length

    def inc_copy_sites_count(self, copied_length: int):
        self.copy_sites_count += copied_length

    def get_dto_str(self) -> str:
        return f"predecessor index: {self.index_in_predecessor}, #copied sites: {self.copy_sites_count}, inserted len: {self.inserted_seq_count}"

    def is_redundant(self):
        return self.copy_sites_count == 0 and self.inserted_seq_count == 0
