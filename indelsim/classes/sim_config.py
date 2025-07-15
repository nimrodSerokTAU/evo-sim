class SimConfiguration:
    original_sequence_length: int
    indel_length_alpha: float
    indel_truncated_length: int
    rate_ins: float
    rate_del: float
    deletion_extra_edge_length: int
    random_seed: int
    # Substitution parameters
    substitution_rate: float
    enable_substitutions: bool
    substitution_model: str
    substitution_algorithm: str

    def __init__(self, original_sequence_length: int, indel_length_alpha: float, indel_truncated_length: int,
                 rate_ins: float, rate_del: float, deletion_extra_edge_length: int, seed: int,
                 substitution_rate: float = 1.0, enable_substitutions: bool = False, 
                 substitution_model: str = 'jtt', substitution_algorithm: str = 'gillespie'):

        self.rate_ins = rate_ins
        self.rate_del = rate_del
        self.original_sequence_length = original_sequence_length
        self.indel_length_alpha = indel_length_alpha
        self.indel_truncated_length = indel_truncated_length
        self.deletion_extra_edge_length = deletion_extra_edge_length
        self.random_seed = seed
        
        # Substitution parameters with validation
        if substitution_rate < 0:
            raise ValueError("substitution_rate must be non-negative")
        if substitution_model not in ['jtt']:
            raise ValueError(f"substitution_model must be 'jtt', got '{substitution_model}'")
        if substitution_algorithm not in ['gillespie', 'matrix']:
            raise ValueError(f"substitution_algorithm must be 'gillespie' or 'matrix', got '{substitution_algorithm}'")
            
        self.substitution_rate = substitution_rate
        self.enable_substitutions = enable_substitutions
        self.substitution_model = substitution_model
        self.substitution_algorithm = substitution_algorithm