from indelsim.classes.indel_event import IndelEvent


class SequenceNodeNaive:
    id: int
    seq: list[int]
    max_count: int

    def __init__(self, seq_id: int, original_sequence: list[int]):
        self.id = seq_id
        self.seq: list[int] = original_sequence.copy()
        self.max_count = max(original_sequence)

    def calculate_event(self, event: IndelEvent):
        if event.length < 0 or event.place > self.get_length():
            return
        if event.is_insertion:
            self.calculate_insertion_event(event)
        else:
            self.calculate_deletion_event(event)

    def calculate_insertion_event(self, event: IndelEvent):
        insertion: list[int] = self.calc_inserted_seq(event.length)
        self.seq[event.place:event.place] = insertion

    def calculate_deletion_event(self, event: IndelEvent):
        if event.place + event.length > 0:
            del self.seq[max(0, event.place):event.place + event.length]

    def calc_inserted_seq(self, length) -> list[int]:
        insertion: list[int] = []
        for i in range(length):
            self.max_count += 1
            insertion.append(self.max_count)
        return insertion

    def get_length(self) -> int:
        return len(self.seq)

    def get_dto(self) -> dict:
        return {'seq': self.seq, 'length': self.get_length()}

    def get_block_dto_from_single_branch(self, orig_seq_count: int) -> dict:
        res: list[str] = []
        prev_val: int = -1
        block_start: int = self.seq[0]
        copied_sites: int  = 0
        inserted_len: int = 0
        is_copied: bool = True
        for i in self.seq:
            if i == prev_val + 1:
                if is_copied:
                    copied_sites += 1
                else:
                    inserted_len += 1
            elif i < orig_seq_count:
                res.append(f'predecessor index: {block_start}, #copied sites: {copied_sites}, inserted len: {inserted_len}')
                copied_sites = 1
                inserted_len = 0
                block_start = i
                is_copied = True
            elif is_copied:
                is_copied = False
                inserted_len += 1
            else:
                res.append(
                    f'predecessor index: {block_start}, #copied sites: {copied_sites}, inserted len: {inserted_len}')
                copied_sites = 0
                inserted_len = 0
                block_start = 0
            prev_val = i
        res.append(
            f'predecessor index: {block_start}, #copied sites: {copied_sites}, inserted len: {inserted_len}')
        return {'blocks': res, 'length': self.get_length()}

