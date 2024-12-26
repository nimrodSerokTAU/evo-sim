class SimConfiguration:
    original_sequence_length: int
    indel_length_alpha: float
    indel_truncated_length: int
    rate_ins: float
    rate_del: float
    deletion_extra_edge_length: int

    def __init__(self, original_sequence_length: int, indel_length_alpha: float, indel_truncated_length: int,
                 rate_ins: float, rate_del: float, deletion_extra_edge_length: int):

        self.rate_ins = rate_ins
        self.rate_del = rate_del
        self.original_sequence_length = original_sequence_length
        self.indel_length_alpha = indel_length_alpha
        self.indel_truncated_length = indel_truncated_length
        self.deletion_extra_edge_length = deletion_extra_edge_length
