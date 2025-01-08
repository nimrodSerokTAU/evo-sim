from typing_extensions import Self
from classes.block import Block

class AVLNode:
    id: int
    bl: Block
    left: Self
    right: Self
    father: Self
    height: int
    length_under_including: int

    def __init__(self, node_id: int, bl: Block):
        self.id = node_id
        self.bl = bl
        self.left = None
        self.right = None
        self.father = None
        self.height = 1
        self.length_under_including = self.calc_len_under_includes_from_bl()

    def set_len_under_includes(self, length_including: int):
        self.length_under_including = length_including

    def calc_len_under_includes_from_bl(self) -> int:
        return self.bl.inserted_seq_count + self.bl.copy_sites_count

    def update_on_same_location(self, copy_sites_count: int | None, inserted_seq_count: int | None):
        if copy_sites_count is not None:
            self.bl.copy_sites_count = copy_sites_count
        if inserted_seq_count is not None:
            self.bl.update_insert_count(inserted_seq_count)
        self.update_length_under_including_recursive()

    def inc_on_same_location(self, delta_copy_sites_count: int | None, delta_inserted_count: int | None):
        if delta_copy_sites_count is not None:
            self.bl.inc_copy_sites_count(delta_copy_sites_count)
        if delta_inserted_count is not None:
            self.bl.inc_insert_count(delta_inserted_count)
        self.update_length_under_including_recursive()

    def update_length_under_including(self):
        self.length_under_including = self.bl.copy_sites_count + self.bl.inserted_seq_count
        if self.left is not None:
            self.length_under_including += self.left.length_under_including
        if self.right is not None:
            self.length_under_including += self.right.length_under_including

    def update_length_under_including_recursive(self):
        self.update_length_under_including()
        if self.father is not None:
            self.father.update_length_under_including_recursive()

    def set_a_father(self, father: Self):
        self.father = father

    def update_key_to_insert_only(self):
        self.bl.index_in_predecessor = -1

    def is_redundant(self) -> bool:
        return self.bl.is_redundant()

    def get_my_own_length(self) -> int:
        return self.bl.copy_sites_count + self.bl.inserted_seq_count

    def get_dto_str(self) -> str:
        return f"id: {self.id}, predecessor index: {self.bl.index_in_predecessor}, #copied sites: {self.bl.copy_sites_count}, " + \
            f"inserted len: {self.bl.inserted_seq_count}, length_under_including: {self.length_under_including}"

    def get_block_str(self) -> str:
        return f"{self.bl.index_in_predecessor}|{self.bl.copy_sites_count}|{self.bl.inserted_seq_count}"


    def get_clean_dto_str(self) -> str:
        return f"predecessor index: {self.bl.index_in_predecessor}, #copied sites: {self.bl.copy_sites_count}, inserted len: {self.bl.inserted_seq_count}"

    def debug_node_structure(self) -> bool:
        k: int = self.bl.inserted_seq_count + self.bl.copy_sites_count
        if self.left is not None:
            if self.left.bl.index_in_predecessor > self.bl.index_in_predecessor:
                return False
            k += self.left.length_under_including
            if not self.left.debug_node_structure():
                return False
        if self.right is not None:
            if self.right.bl.index_in_predecessor < self.bl.index_in_predecessor:
                return False
            k += self.right.length_under_including
            if not self.right.debug_node_structure():
                return False
        return self.length_under_including == k
