class Block:
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

    def get_my_length(self) -> int:
        return self.copy_sites_count + self.inserted_seq_count

    def get_dto_str(self) -> str:
        index_in_predecessor: int = -1 if self.index_in_predecessor == 0 and self.copy_sites_count == 0 else self.index_in_predecessor
        return f"predecessor index: {index_in_predecessor}, #copied sites: {self.copy_sites_count}, inserted len: {self.inserted_seq_count}"

    def get_block_str(self) -> str:
        index_in_predecessor: int = -1 if self.index_in_predecessor == 0 and self.copy_sites_count == 0 else self.index_in_predecessor
        return f"{index_in_predecessor}|{self.copy_sites_count}|{self.inserted_seq_count}"


    def is_redundant(self):
        return self.copy_sites_count == 0 and self.inserted_seq_count == 0

    def __repr__(self):
        return f"{self.index_in_predecessor}|{self.copy_sites_count}|{self.inserted_seq_count}"
