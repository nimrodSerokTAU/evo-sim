from enum import Enum, IntEnum


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


class AminoAcid(IntEnum):
    A = 0   # Alanine
    R = 1   # Arginine
    N = 2   # Asparagine
    D = 3   # Aspartic acid
    C = 4   # Cysteine
    Q = 5   # Glutamine
    E = 6   # Glutamic acid
    G = 7   # Glycine
    H = 8   # Histidine
    I = 9   # Isoleucine
    L = 10  # Leucine
    K = 11  # Lysine
    M = 12  # Methionine
    F = 13  # Phenylalanine
    P = 14  # Proline
    S = 15  # Serine
    T = 16  # Threonine
    W = 17  # Tryptophan
    Y = 18  # Tyrosine
    V = 19  # Valine

def amino_acid_to_index(aa: str) -> int:
    """Convert amino acid character to index (0-19)."""
    try:
        return AminoAcid[aa].value
    except KeyError:
        raise ValueError(f"Invalid amino acid: {aa}")

def index_to_amino_acid(index: int) -> str:
    """Convert index (0-19) to amino acid character."""
    try:
        return AminoAcid(index).name
    except ValueError:
        raise ValueError(f"Invalid amino acid index: {index}")
    
PROTEIN_ALPHABET = [aa.name for aa in AminoAcid]