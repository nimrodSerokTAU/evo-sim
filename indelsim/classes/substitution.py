#!/usr/bin/env python3
"""
Implements both a Gillespie-style CTMC sampler (`evolve_branch_substitutions_gillespie`)
and a matrix-exponential sampler (`evolve_branch_substitutions_jtt`)
under the JTT model.
"""

from __future__ import annotations
from typing import List, Optional


import numpy as np
from indelsim.classes.jtt import get_jtt_model
from indelsim.enums import amino_acid_to_index, index_to_amino_acid

# Constants for validation
MAX_BRANCH_LENGTH = 1000.0
MIN_SUBSTITUTION_RATE = 1e-10

class SubstitutionEvolver:
    """
    Evolves amino-acid sequences along a branch under the JTT model.
    
    Parameters
    ----------
    substitution_rate
        Scalar multiplier applied to JTT's normalised rate matrix.
    seed
        If given, a NumPy PCG64 RNG is created with this seed;
        otherwise a fresh unpredictable RNG is used.
    """
    
    def __init__(self, substitution_rate: float = 1.0, seed: Optional[int] = None):
        if substitution_rate < MIN_SUBSTITUTION_RATE:
            raise ValueError(f"substitution_rate must be >= {MIN_SUBSTITUTION_RATE}")
            
        self.substitution_rate = float(substitution_rate)
        self.jtt_model = get_jtt_model()

        if seed is not None:
            self.rng = np.random.default_rng(seed)
        else:
            self.rng = np.random.default_rng()
    
    def evolve_branch_substitutions_gillespie(
        self,
        sequence: List[int],
        branch_length: float
    ) -> np.ndarray:
        """
        Simulate substitutions along a branch with the Gillespie algorithm.

        Parameters
        ----------
        sequence
            List of length L with integer amino-acid codes 0-19.
        branch_length
            Time units for this branch (must be non-negative).

        Returns
        -------
        list[int]
            New list (length L) after substitutions.
        """
        self._validate_inputs(sequence, branch_length)

        Q: np.ndarray = self.jtt_model.rate_matrix * self.substitution_rate
        L: int = len(sequence)
        seq: list[int] = sequence.copy()

        # Pre-compute exit rates λ_i = −Q_ii
        exit_rates = -Q[seq, seq]
        total_rate: float = exit_rates.sum()

        t: float = 0.0
        rng = self.rng

        while t < branch_length and total_rate > 0.0:
            # 1. waiting time to next event
            t += rng.exponential(1.0 / total_rate)
            if t >= branch_length:
                break

            # 2. choose site, proportional to exit rate
            site = rng.choice(L, p=exit_rates / total_rate)

            old_aa = seq[site]

            # 3. draw target residue conditional on leaving old_aa
            probs = Q[old_aa].copy()
            probs[old_aa] = 0.0
            probs /= exit_rates[site]          # normalise row
            new_aa = rng.choice(20, p=probs)

            # 4. update sequence and rates
            seq[site] = new_aa
            new_exit = -Q[new_aa, new_aa]
            total_rate += new_exit - exit_rates[site]
            exit_rates[site] = new_exit

        return seq
    
    def evolve_branch_substitutions_jtt(
        self,
        sequence: np.ndarray,
        branch_length: float
    ) -> np.ndarray:
        """
        Simulate substitutions via pre-computed transition matrix exp(Q t).

        Faster for long sequences / many sites, but incurs one matrix
        exponential per distinct branch length.

        Returns a fresh list with evolved residues.
        """

        self._validate_inputs(sequence, branch_length)

        eff_time = branch_length * self.substitution_rate
        P_t: np.ndarray = self.jtt_model.transition_probability(eff_time)
        probs = P_t[sequence]          # shape (L, 20)
        cumprob = np.cumsum(probs, axis=1)
        u = self.rng.random((len(sequence), 1))
        evolved = np.argmax(u < cumprob, axis=1).astype(np.uint8)
        return evolved
    
    def evolve_sequence_chars(self, sequence_chars: List[str], branch_length: float) -> List[str]:
        """Same as gillespie method, but with ACDE.. chars"""
        indices = [amino_acid_to_index(aa) for aa in sequence_chars]
        evolved_indices = self.evolve_branch_substitutions_gillespie(indices, branch_length)
        return [index_to_amino_acid(idx) for idx in evolved_indices]

    def _validate_inputs(self, sequence: np.ndarray, branch_length: float) -> None:
        if branch_length < 0.0:
            raise ValueError("branch_length must be non-negative")
        if branch_length > MAX_BRANCH_LENGTH:
            raise ValueError(f"branch_length {branch_length} exceeds maximum {MAX_BRANCH_LENGTH}")
        if sequence.size == 0:
            raise ValueError("sequence cannot be empty")
