class SimConfiguration:
    original_sequence: str | None
    original_sequence_length: int | None
    substitutions_per_site_rate: float
    indel_length_alpha: float
    indel_truncated_length: int
    is_use_nucleotides: bool
    indel_per_sub_ratio: float

    def __init__(self, original_sequence: str | None, original_sequence_length: int | None,
                 substitutions_per_site_rate: float, indel_length_alpha: float, indel_truncated_length: int,
                 is_use_nucleotides: bool, indel_per_sub_ratio: float):
        self.original_sequence = original_sequence
        self.original_sequence_length = original_sequence_length
        self.substitutions_per_site_rate = substitutions_per_site_rate
        self.indel_length_alpha = indel_length_alpha
        self.indel_truncated_length = indel_truncated_length
        self.is_use_nucleotides = is_use_nucleotides
        self.indel_per_sub_ratio = indel_per_sub_ratio
