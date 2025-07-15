#!/usr/bin/env python3
"""
CLI tool to simulate substitutions in a sequence.

"""

import numpy as np
from typing import List, Optional
from indelsim.classes.jtt import JTTModel
from indelsim.enums import AminoAcid, amino_acid_to_index, index_to_amino_acid

# Constants for validation
MAX_BRANCH_LENGTH = 1000.0
MIN_SUBSTITUTION_RATE = 1e-10

class SubstitutionEvolver:
    """
    Evolver for amino acid substitutions using the JTT model.
    
    This class simulates substitution-only evolution along phylogenetic branches,
    maintaining constant sequence length while applying amino acid changes based
    on the JTT substitution model using the Gillespie algorithm.
    """
    
    def __init__(self, substitution_rate: float = 1.0, seed: Optional[int] = None):
        if substitution_rate < MIN_SUBSTITUTION_RATE:
            raise ValueError(f"substitution_rate must be >= {MIN_SUBSTITUTION_RATE}")
            
        self.substitution_rate = substitution_rate
        self.jtt_model = JTTModel()
        self.jtt_model.compute_model()

        if seed is not None:
            self.rng = np.random.default_rng(seed)
        else:
            self.rng = np.random.default_rng()
    
    def evolve_branch_substitutions_gillespie(
        self,
        sequence: List[int],
        branch_length: float
    ) -> List[int]:
        """
        Evolve a sequence along a branch using substitution-only evolution.
        
        Args:
            sequence: List of amino acid indices (0-19)
            branch_length: Length of the phylogenetic branch
            
        Returns:
            Evolved sequence with substitutions applied
        """
        if branch_length < 0:
            raise ValueError("branch_length must be non-negative")
        if branch_length > MAX_BRANCH_LENGTH:
            raise ValueError(f"branch_length {branch_length} exceeds maximum {MAX_BRANCH_LENGTH}")
        if sequence is None or len(sequence) == 0:
            raise ValueError("sequence cannot be empty")

        # local copies for speed
        rng = self.rng
        Q = self.jtt_model.rate_matrix * self.substitution_rate

        seq = sequence.copy()
        L = len(seq)
        # Exit rates per site
        exit_r = np.asarray([-Q[aa, aa] for aa in seq], dtype=float)
        total_r = exit_r.sum()

        t = 0.0
        while t < branch_length and total_r > 0.0:
            # 1. waiting time
            t += rng.exponential(1.0 / total_r)
            if t >= branch_length:
                break

            # 2. choose site
            site = rng.choice(L, p=exit_r / total_r)

            old_aa = seq[site]
            if exit_r[site] == 0.0:  # should never happen if Q is valid
                continue

            # 3. choose new aa from chain
            probs = Q[old_aa].copy()
            probs[old_aa] = 0.0 # remove diagonal
            probs /= exit_r[site]  # normalise
            new_aa = rng.choice(20, p=probs)

            # 4. update seq and rates
            seq[site] = new_aa

            new_exit = -Q[new_aa, new_aa]
            total_r += new_exit - exit_r[site]
            exit_r[site] = new_exit

        return seq
    
    def evolve_branch_substitutions_jtt(
        self,
        sequence: List[int],
        branch_length: float
    ) -> List[int]:
        """
        Evolve a protein sequence along a single tree branch under the
        Jones–Taylor–Thornton (JTT) substitution model using matrix exponentiation.

        Args:
            sequence: List of amino acid indices (0-19)
            branch_length: Length of the phylogenetic branch
            
        Returns:
            Evolved sequence with substitutions applied
        """
        if branch_length < 0:
            raise ValueError("branch_length must be non-negative")
        if branch_length > MAX_BRANCH_LENGTH:
            raise ValueError(f"branch_length {branch_length} exceeds maximum {MAX_BRANCH_LENGTH}")
        if sequence is None or len(sequence) == 0:
            raise ValueError("sequence cannot be empty")

        effective_time = branch_length * self.substitution_rate

        # Vectorized approach using cumulative probabilities
        P_t = self.jtt_model.transition_probability(effective_time) # (20, 20)
        probs = P_t[np.asarray(sequence, dtype=np.uint8)]

        cumprob = np.cumsum(probs, axis=1)
        u = self.rng.random((len(sequence), 1))

        evolved_sequence = np.argmax(u < cumprob, axis=1).astype(np.uint8)

        return evolved_sequence.tolist()
    
    def evolve_sequence_chars(self, sequence_chars: List[str], branch_length: float) -> List[str]:
        """
        Convenience method for evolving character sequences.
        
        Args:
            sequence_chars: List of amino acid characters
            branch_length: Length of the phylogenetic branch
            
        Returns:
            Evolved sequence as list of amino acid characters
        """
        # Convert to indices
        indices = [amino_acid_to_index(aa) for aa in sequence_chars]
        
        # Evolve using Gillespie algorithm
        evolved_indices = self.evolve_branch_substitutions_gillespie(indices, branch_length)
        
        # Convert back to characters
        return [index_to_amino_acid(idx) for idx in evolved_indices]
