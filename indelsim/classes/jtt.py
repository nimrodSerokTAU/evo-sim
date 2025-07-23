"""
Enhanced JTT (Jones-Taylor-Thornton) Substitution Model Implementation

This module provides a robust, object-oriented implementation of the JTT model
for calculating amino acid substitution probabilities over evolutionary time.

# Simple usage
model = JTTModel()
P_t = model.transition_probability(1.0)  # Auto-computes model if needed

# With custom configuration
config = JTTConfig(ZERO_TOLERANCE=1e-10, CACHE_SIZE=256)
model = JTTModel(config)

# Validate model properties
validation_results = model.validate_model_properties()
all_passed = all(validation_results.values())

"""

import numpy as np
from typing import Optional, Tuple, Dict, Any
from dataclasses import dataclass
from scipy.linalg import expm
import warnings
from functools import lru_cache

# -------------------------------------------------------------------------- #
# Global singleton helper                                                    #
# -------------------------------------------------------------------------- #

__all__ = ["JTTConfig", "JTTModel", "JTTModelError", "get_jtt_model"]

# Cache a single computed model so multiple SubstitutionEvolver instances
# don’t redo the expensive eigendecomposition.
_GLOBAL_MODEL: "JTTModel | None" = None


def get_jtt_model() -> "JTTModel":
    """Return a globally shared, pre-computed JTTModel instance."""
    global _GLOBAL_MODEL
    if _GLOBAL_MODEL is None:
        _GLOBAL_MODEL = JTTModel()
        _GLOBAL_MODEL.compute_model()
    return _GLOBAL_MODEL

# -------------------------------------------------------------------------- #

@dataclass(frozen=True)
class JTTConfig:
    """Configuration parameters for the JTT model."""
    
    # Number of amino acids
    NUM_AMINO_ACIDS: int = 20
    
    # Numerical tolerances
    ZERO_TOLERANCE: float = 1e-8
    BALANCE_TOLERANCE: float = 1e-8
    RECONSTRUCTION_TOLERANCE: float = 1e-10
    
    # Time limits for calculations
    MIN_TIME: float = 1e-10
    MAX_TIME: float = 1000.0
    
    # Cache size for transition matrices
    CACHE_SIZE: int = 128


class JTTModelError(Exception):
    """Base exception for JTT model errors."""
    pass


class JTTModel:
    """
    Jones-Taylor-Thornton substitution model for amino acid evolution.
    
    This class encapsulates the JTT model, providing methods to compute
    transition probabilities and validate model properties.

    # TODO: share model across simulations. 
    """
    
    def __init__(self, config: Optional[JTTConfig] = None):
        """
        Initialize the JTT model.
        
        Args:
            config: Configuration parameters (uses defaults if None)
        """
        self.config = config or JTTConfig()
        self._is_computed = False
        self._cache = {}
        
        # Model components (initialized when computed)
        self._Q: Optional[np.ndarray] = None
        self._eigenvalues: Optional[np.ndarray] = None
        self._Y: Optional[np.ndarray] = None
        self._Y_inv: Optional[np.ndarray] = None
        self._pi: Optional[np.ndarray] = None
        self._v_matrix: Optional[np.ndarray] = None
        
        # Initialize model data
        self._initialize_data()
    
    def _initialize_data(self) -> None:
        """Initialize the empirical data for the JTT model."""
        # Amino acid equilibrium frequencies
        self._pi_data = np.array([
            0.076748, 0.051691, 0.042645, 0.051544, 0.019803,
            0.040752, 0.061830, 0.073152, 0.022944, 0.053761,
            0.091904, 0.058676, 0.023826, 0.040126, 0.050901,
            0.068765, 0.058565, 0.014261, 0.032102, 0.066005
        ])
        
        # Empirical substitution counts (upper triangle)
        self._tri_counts = np.array([
            58, 54, 45, 81, 16, 528, 56, 113, 34, 10, 57, 310, 86, 49, 9,
            105, 29, 58, 767, 5, 323, 179, 137, 81, 130, 59, 26, 119, 27,
            328, 391, 112, 69, 597, 26, 23, 36, 22, 47, 11, 17, 9, 12, 6,
            16, 30, 38, 12, 7, 23, 72, 9, 6, 56, 229, 35, 646, 263, 26,
            7, 292, 181, 27, 45, 21, 14, 54, 44, 30, 15, 31, 43, 18, 14,
            33, 479, 388, 65, 15, 5, 10, 4, 78, 4, 5, 5, 40, 89, 248, 4,
            43, 194, 74, 15, 15, 14, 164, 18, 24, 115, 10, 102, 21, 16,
            17, 378, 101, 503, 59, 223, 53, 30, 201, 73, 40, 59, 47, 29,
            92, 285, 475, 64, 232, 38, 42, 51, 32, 33, 46, 245, 25, 103,
            226, 12, 118, 477, 9, 126, 8, 4, 115, 18, 10, 55, 8, 9, 52,
            10, 24, 53, 6, 35, 12, 11, 20, 70, 46, 209, 24, 7, 8, 573,
            32, 24, 8, 18, 536, 10, 63, 21, 71, 298, 17, 16, 31, 62, 20,
            45, 47, 11, 961, 180, 14, 323, 62, 23, 38, 112, 25, 16
        ], dtype=float)
        
        # Validate data integrity
        self._validate_data()
    
    def _validate_data(self) -> None:
        """Validate the empirical data."""
        n = self.config.NUM_AMINO_ACIDS
        expected_tri_count = n * (n - 1) // 2
        
        if len(self._pi_data) != n:
            raise JTTModelError(f"Expected {n} amino acid frequencies, got {len(self._pi_data)}")
        
        if len(self._tri_counts) != expected_tri_count:
            raise JTTModelError(f"Expected {expected_tri_count} substitution counts, got {len(self._tri_counts)}")
        
        if not np.isclose(self._pi_data.sum(), 1.0):
            raise JTTModelError(f"Amino acid frequencies sum to {self._pi_data.sum():.6f}, expected 1.0")
        
        if np.any(self._pi_data <= 0):
            raise JTTModelError("All amino acid frequencies must be positive")
        
        if np.any(self._tri_counts < 0):
            raise JTTModelError("All substitution counts must be non-negative")
    
    def compute_model(self) -> None:
        """
        Compute the JTT model components.
        
        This method builds the rate matrix Q, performs eigendecomposition,
        and validates the resulting model.
        """
        if self._is_computed:
            return
        
        try:
            # Build the symmetric substitution matrix S
            S = self._build_substitution_matrix()
            
            # Construct the rate matrix Q
            Q = self._build_rate_matrix(S)
            
            # Scale Q to evolutionary time units
            Q_scaled = self._scale_rate_matrix(Q)
            
            # Perform eigendecomposition
            eigenvalues, Y, Y_inv, v_matrix = self._eigendecomposition(Q_scaled)
            
            # Validate the decomposition
            self._validate_decomposition(Q_scaled, eigenvalues, Y, Y_inv)
            
            # Store results
            self._Q = Q_scaled
            self._eigenvalues = eigenvalues
            self._Y = Y
            self._Y_inv = Y_inv
            self._pi = self._pi_data.copy()
            self._v_matrix = v_matrix
            
            self._is_computed = True
            
        except Exception as e:
            raise JTTModelError(f"Failed to compute JTT model: {str(e)}") from e
    
    def _build_substitution_matrix(self) -> np.ndarray:
        """Build the symmetric substitution matrix S from empirical counts."""
        n = self.config.NUM_AMINO_ACIDS
        S = np.zeros((n, n))
        
        k = 0
        for i in range(n):
            for j in range(i):
                S[i, j] = S[j, i] = self._tri_counts[k]
                k += 1
        
        return S
    
    def _build_rate_matrix(self, S: np.ndarray) -> np.ndarray:
        """Build the rate matrix Q from the substitution matrix S."""
        # Q_ij = S_ij * π_j for i ≠ j
        Q = S * self._pi_data[np.newaxis, :]
        
        # Set diagonal to zero initially
        np.fill_diagonal(Q, 0.0)
        
        # Set diagonal elements to ensure row sums are zero
        Q -= np.diag(Q.sum(axis=1))
        
        return Q
    
    def _scale_rate_matrix(self, Q: np.ndarray) -> np.ndarray:
        """Scale the rate matrix Q to evolutionary time units."""
        # Calculate the average substitution rate
        rate = -(self._pi_data * np.diag(Q)).sum()
        
        if rate <= 0:
            raise JTTModelError(f"Invalid substitution rate: {rate}")
        
        # Scale Q so that the average rate is 1.0
        Q_scaled = Q / rate
        
        return Q_scaled
    
    def _eigendecomposition(self, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Perform eigendecomposition of the rate matrix Q."""
        # Use symmetrization for numerical stability
        sqrt_pi = np.sqrt(self._pi_data)
        
        # Avoid division by very small numbers
        min_sqrt_pi = np.min(sqrt_pi)
        if min_sqrt_pi < 1e-10:
            warnings.warn(f"Very small amino acid frequency detected: {min_sqrt_pi**2:.2e}")
        
        # Build transformation matrices
        T = np.diag(sqrt_pi)
        T_inv = np.diag(1.0 / sqrt_pi)
        
        # Symmetrize: M = T * Q * T^(-1)
        M = T @ Q @ T_inv
        
        # Eigendecomposition of symmetric matrix
        eigenvalues, v_matrix = np.linalg.eigh(M)
        
        # Sort eigenvalues and eigenvectors in descending order
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        v_matrix = v_matrix[:, idx]
        
        # Transform back to original space
        Y = T_inv @ v_matrix
        Y_inv = v_matrix.T @ T
        
        return eigenvalues, Y, Y_inv, v_matrix
    
    def _validate_decomposition(self, Q: np.ndarray, eigenvalues: np.ndarray, Y: np.ndarray, Y_inv: np.ndarray) -> None:
        """Validate the eigendecomposition results."""
        # Check reconstruction
        lambda_matrix = np.diag(eigenvalues)
        Q_reconstructed = Y @ lambda_matrix @ Y_inv
        
        if not np.allclose(Q, Q_reconstructed, atol=self.config.RECONSTRUCTION_TOLERANCE):
            raise JTTModelError("Eigendecomposition reconstruction failed")
        
        # Check eigenvalue structure
        zero_evals = np.abs(eigenvalues) < self.config.ZERO_TOLERANCE
        neg_evals = eigenvalues < -self.config.ZERO_TOLERANCE
        
        if np.sum(zero_evals) != 1:
            raise JTTModelError(f"Expected 1 zero eigenvalue, found {np.sum(zero_evals)}")
        
        if np.sum(neg_evals) != self.config.NUM_AMINO_ACIDS - 1:
            raise JTTModelError(f"Expected {self.config.NUM_AMINO_ACIDS - 1} negative eigenvalues, found {np.sum(neg_evals)}")
    

    @lru_cache(maxsize=128)
    def _transition_probability_cached(self, t_quantized: float) -> np.ndarray:
        """Internal cached method for transition probability calculation."""
        # Calculate P(t) = Y * exp(Λt) * Y^(-1)
        exp_lambda_t = np.diag(np.exp(self._eigenvalues * t_quantized))
        P_t = self._Y @ exp_lambda_t @ self._Y_inv
        
        # Validate result
        self._validate_transition_matrix(P_t)
        
        return P_t
    
    def transition_probability(self, t: float) -> np.ndarray:
        """
        Calculate the transition probability matrix P(t) = exp(Qt).
        
        Args:
            t: Evolutionary time (must be positive)
            
        Returns:
            Transition probability matrix P(t)
            
        Raises:
            JTTModelError: If time is invalid or model not computed
        """
        if not self._is_computed:
            self.compute_model()
        
        # Validate time parameter
        if not isinstance(t, (int, float)) or np.isnan(t) or np.isinf(t):
            raise JTTModelError(f"Invalid time value: {t}")
        
        if t < self.config.MIN_TIME:
            raise JTTModelError(f"Time {t} is too small (minimum: {self.config.MIN_TIME})")
        
        if t > self.config.MAX_TIME:
            warnings.warn(f"Large time value {t} may cause numerical issues")
        
        # Quantize time to 10 decimal places for cache efficiency
        t_quantized = round(t, 10)
        
        return self._transition_probability_cached(t_quantized)
    
    def _validate_transition_matrix(self, P: np.ndarray) -> None:
        """Validate the transition probability matrix."""
        # Check row sums
        row_sums = P.sum(axis=1)
        if not np.allclose(row_sums, 1.0, atol=self.config.ZERO_TOLERANCE):
            raise JTTModelError(f"Transition matrix rows don't sum to 1.0: {row_sums}")
        
        # Check bounds
        if np.any(P < -self.config.ZERO_TOLERANCE) or np.any(P > 1.0 + self.config.ZERO_TOLERANCE):
            raise JTTModelError("Transition matrix contains invalid probabilities")
    
    def validate_model_properties(self) -> Dict[str, bool]:
        """
        Validate all mathematical properties of the JTT model.
        
        Returns:
            Dictionary of validation results
        """
        if not self._is_computed:
            self.compute_model()
        
        results = {}
        
        # Test 1: Q row sums should be zero
        results['q_row_sums'] = np.allclose(self._Q.sum(axis=1), 0, atol=self.config.ZERO_TOLERANCE)
        
        # Test 2: Detailed balance
        left = self._pi[:, None] * self._Q
        right = left.T
        results['detailed_balance'] = np.allclose(left, right, atol=self.config.BALANCE_TOLERANCE)
        
        # Test 3: Eigenvalue structure
        zero_evals = np.abs(self._eigenvalues) < self.config.ZERO_TOLERANCE
        neg_evals = self._eigenvalues < -self.config.ZERO_TOLERANCE
        results['eigenvalue_structure'] = (np.sum(zero_evals) == 1 and np.sum(neg_evals) == 19)
        
        # Test 4: Eigenvector orthogonality
        results['orthogonality'] = np.allclose(self._v_matrix.T @ self._v_matrix, np.eye(20), atol=self.config.ZERO_TOLERANCE)
        
        # Test 5: Stationarity
        P_test = self.transition_probability(1.0)
        results['stationarity'] = np.allclose(self._pi @ P_test, self._pi, atol=self.config.ZERO_TOLERANCE)
        
        # Test 6: Comparison with matrix exponential
        P_direct = expm(self._Q * 1.0)
        results['expm_comparison'] = np.allclose(P_test, P_direct, atol=self.config.ZERO_TOLERANCE)
        
        return results
    
    @property
    def rate_matrix(self) -> np.ndarray:
        """Get the rate matrix Q."""
        if not self._is_computed:
            self.compute_model()
        return self._Q.copy()
    
    @property
    def equilibrium_frequencies(self) -> np.ndarray:
        """Get the normalised equilibrium amino acid frequencies."""
        return (self._pi_data / self._pi_data.sum()).copy()
    
    @property
    def eigenvalues(self) -> np.ndarray:
        """Get the eigenvalues of the rate matrix."""
        if not self._is_computed:
            self.compute_model()
        return self._eigenvalues.copy()
    
    def __repr__(self) -> str:
        """String representation of the JTT model."""
        status = "computed" if self._is_computed else "not computed"
        return f"JTTModel(status={status}, amino_acids={self.config.NUM_AMINO_ACIDS})"


# Convenience functions for backward compatibility
def compute_jtt_model() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Legacy function for backward compatibility.
    
    Returns:
        Tuple of (Q, eigenvalues, Y, Y_inv, pi)
    """
    model = JTTModel()
    model.compute_model()
    return (model._Q, model._eigenvalues, model._Y, model._Y_inv, model._pi)


def calculate_p_t(t: float, eigenvalues: np.ndarray, Y: np.ndarray, Y_inv: np.ndarray) -> np.ndarray:
    """
    Legacy function for backward compatibility.
    
    Args:
        t: Evolutionary time
        eigenvalues: Eigenvalues of the rate matrix
        Y: Eigenvector matrix
        Y_inv: Inverse eigenvector matrix
        
    Returns:
        Transition probability matrix P(t)
    """
    exp_lambda_t = np.diag(np.exp(eigenvalues * t))
    P_t = Y @ exp_lambda_t @ Y_inv
    return P_t


# Example usage and testing
if __name__ == "__main__":
    # Create and compute the JTT model
    model = JTTModel()
    print(f"JTT Model: {model}")
    
    # Compute model components
    model.compute_model()
    print(f"JTT Model: {model}")
    
    # Calculate transition probabilities
    t_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    
    print("\n--- Transition Probabilities ---")
    for t in t_values:
        P_t = model.transition_probability(t)
        print(f"t={t}: P(t) shape={P_t.shape}, min={P_t.min():.6f}, max={P_t.max():.6f}")
    
    # Validate model properties
    print("\n--- Model Validation ---")
    validation_results = model.validate_model_properties()
    for test_name, passed in validation_results.items():
        status = "PASS" if passed else "FAIL"
        print(f"{test_name}: {status}")
    
    # Display some model properties
    print(f"\n--- Model Properties ---")
    print(f"Equilibrium frequencies sum: {model.equilibrium_frequencies.sum():.6f}")
    print(f"Eigenvalues range: [{model.eigenvalues.min():.6f}, {model.eigenvalues.max():.6f}]")
    print(f"Rate matrix trace: {np.trace(model.rate_matrix):.6f}")


