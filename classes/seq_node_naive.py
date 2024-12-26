from classes.indel_event import IndelEvent


class SequenceNodeNaive:
    id: int
    seq: list[int]
    max_count: int

    def __init__(self, seq_id: int, original_sequence_length: int):
        self.id = seq_id
        self.seq: list[int] = [i for i in range(original_sequence_length)]
        self.max_count = original_sequence_length
        # for i in range(original_sequence_length):
        #     self.seq.append(i)

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
            insertion.append(self.max_count)
            self.max_count += 1
        return insertion

    def get_length(self) -> int:
        return len(self.seq)

    def get_dto(self) -> dict:
        return {'seq': self.seq, 'length': self.get_length()}

