from enum import Enum


class SimulationTypes(Enum):
    NAIVE = 0,
    BLOCK_LIST = 1,
    BLOCK_TREE = 2,


class EventSubTypes(Enum):
    INSERTION_AT_START = 0
    INSERTION_AT_START_ADDITION = 1
    INSERTION_INSIDE_COPIED = 2
    INSERTION_INSIDE_INSERTED = 3
    INSERTION_AT_END = 4
    DELETION_INSIDE_COPIED_CONTAINED_AT_MID = 5
    DELETION_INSIDE_COPIED_CONTAINED_AT_START = 6
    DELETION_INSIDE_COPIED_UNCONTAINED = 7
    DELETION_OF_COPIED = 8
    DELETION_ALL_COPIED_UNCONTAINED = 9
    DELETION_ALL_COPIED_UNCONTAINED_AT_START = 10
    DELETION_INSIDE_INSERTED_CONTAINED = 11
    DELETION_INSIDE_INSERTED_UNCONTAINED = 12
    DELETION_OF_INSERTED = 13
    OUT_OF_SEQUENCE = 14


