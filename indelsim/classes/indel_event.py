class IndelEvent:
    __slots__ = ['is_insertion', 'place', 'length']

    is_insertion: bool
    place: int
    length: int

    def __init__(self, is_insertion: bool, place: int, length: int):
        self.is_insertion = is_insertion
        self.place = place
        self.length = length
        # Support for deletions and insertions that start before the beginning of the sequence, 
        # is it necessary for insertions?
        if place < 0: 
            self.length = length - abs(place)
            self.place = 0


    def __repr__(self):
        return f"IndelEvent(is_insertion={(self.is_insertion)}, length={self.length}, place={self.place})"