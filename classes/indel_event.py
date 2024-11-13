class IndelEvent:
    is_insertion: bool
    place: int
    length: int

    def __init__(self, is_insertion: bool, place: int, length: int):
        self.is_insertion = is_insertion
        self.place = place
        self.length = length
        if place < 0:
            self.length = length - abs(place)
            self.place = 0
