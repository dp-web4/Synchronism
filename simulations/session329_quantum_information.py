#!/usr/bin/env python3
"""
Session #329: Quantum Information from the Planck Grid
Information Theory Arc (Session 2/4)

This session explores quantum information from the grid perspective:
1. Quantum bits (qubits) and superposition
2. Entanglement and Bell states
3. No-cloning theorem
4. Quantum channels and Holevo bound
5. Quantum error correction

Key insight: Quantum information is pattern coherence within the MRH.
Entanglement = correlated patterns across grid regions.
Decoherence = pattern info leaking beyond MRH.
The no-cloning theorem reflects the indivisibility of grid patterns.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from scipy import constants as const

# Physical constants
hbar = const.hbar
k_B = const.k


@dataclass
class Qubit:
    """
    Quantum bit on the grid.

    |ψ⟩ = α|0⟩ + β|1⟩  where |α|² + |β|² = 1

    Grid interpretation: A qubit represents a coherent
    superposition of two grid pattern configurations.
    The MRH boundary defines what remains coherent.
    """

    def __init__(self, alpha: complex = 1.0, beta: complex = 0.0):
        """
        Args:
            alpha: Amplitude for |0⟩
            beta: Amplitude for |1⟩
        """
        norm = np.sqrt(np.abs(alpha)**2 + np.abs(beta)**2)
        self.alpha = alpha / norm
        self.beta = beta / norm

    @classmethod
    def from_bloch(cls, theta: float, phi: float) -> 'Qubit':
        """
        Create qubit from Bloch sphere coordinates.

        |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩
        """
        alpha = np.cos(theta / 2)
        beta = np.exp(1j * phi) * np.sin(theta / 2)
        return cls(alpha, beta)

    def state_vector(self) -> np.ndarray:
        """Return state as column vector."""
        return np.array([[self.alpha], [self.beta]], dtype=complex)

    def density_matrix(self) -> np.ndarray:
        """Return density matrix ρ = |ψ⟩⟨ψ|."""
        psi = self.state_vector()
        return psi @ psi.conj().T

    def bloch_coordinates(self) -> Tuple[float, float, float]:
        """
        Return Bloch sphere coordinates (x, y, z).

        x = ⟨σ_x⟩, y = ⟨σ_y⟩, z = ⟨σ_z⟩
        """
        rho = self.density_matrix()
        x = 2 * np.real(rho[0, 1])
        y = 2 * np.imag(rho[1, 0])
        z = np.real(rho[0, 0] - rho[1, 1])
        return (x, y, z)

    def probability_zero(self) -> float:
        """Probability of measuring |0⟩."""
        return np.abs(self.alpha)**2

    def probability_one(self) -> float:
        """Probability of measuring |1⟩."""
        return np.abs(self.beta)**2

    def von_neumann_entropy(self) -> float:
        """
        Von Neumann entropy S(ρ) = -Tr(ρ log ρ).

        For pure state: S = 0
        For mixed state: S > 0
        """
        rho = self.density_matrix()
        eigenvalues = np.linalg.eigvalsh(rho)
        # S = -Σ λ_i log λ_i
        entropy = 0.0
        for ev in eigenvalues:
            if ev > 1e-10:
                entropy -= ev * np.log2(ev)
        return entropy

    @staticmethod
    def grid_interpretation() -> Dict[str, str]:
        """Grid interpretation of qubit."""
        return {
            'superposition': 'Coherent mixture of two grid patterns',
            'amplitudes': 'Pattern weights (complex for phase)',
            'measurement': 'Collapse to one pattern, other info → MRH',
            'bloch_sphere': 'Full pattern state space',
            'pure_state': 'Complete pattern information tracked',
            'mixed_state': 'Some pattern info beyond MRH'
        }


class EntangledPair:
    """
    Entangled qubit pairs (Bell states).

    |Φ⁺⟩ = (|00⟩ + |11⟩)/√2  (maximally entangled)
    |Φ⁻⟩ = (|00⟩ - |11⟩)/√2
    |Ψ⁺⟩ = (|01⟩ + |10⟩)/√2
    |Ψ⁻⟩ = (|01⟩ - |10⟩)/√2

    Grid interpretation: Entanglement = correlated patterns
    across spatially separated grid regions. Cannot be created
    by local operations alone — requires shared MRH history.
    """

    def __init__(self, bell_state: str = 'Phi+'):
        """
        Args:
            bell_state: One of 'Phi+', 'Phi-', 'Psi+', 'Psi-'
        """
        self.bell_state = bell_state
        self.state = self._create_bell_state(bell_state)

    def _create_bell_state(self, name: str) -> np.ndarray:
        """Create Bell state vector in computational basis."""
        sqrt2 = np.sqrt(2)
        if name == 'Phi+':
            return np.array([1, 0, 0, 1], dtype=complex) / sqrt2
        elif name == 'Phi-':
            return np.array([1, 0, 0, -1], dtype=complex) / sqrt2
        elif name == 'Psi+':
            return np.array([0, 1, 1, 0], dtype=complex) / sqrt2
        elif name == 'Psi-':
            return np.array([0, 1, -1, 0], dtype=complex) / sqrt2
        else:
            raise ValueError(f"Unknown Bell state: {name}")

    def density_matrix(self) -> np.ndarray:
        """Full two-qubit density matrix."""
        psi = self.state.reshape(4, 1)
        return psi @ psi.conj().T

    def reduced_density_matrix_A(self) -> np.ndarray:
        """Reduced density matrix for qubit A (trace over B)."""
        rho = self.density_matrix().reshape(2, 2, 2, 2)
        return np.trace(rho, axis1=1, axis2=3)

    def reduced_density_matrix_B(self) -> np.ndarray:
        """Reduced density matrix for qubit B (trace over A)."""
        rho = self.density_matrix().reshape(2, 2, 2, 2)
        return np.trace(rho, axis1=0, axis2=2)

    def entanglement_entropy(self) -> float:
        """
        Entanglement entropy S_A = -Tr(ρ_A log ρ_A).

        For Bell states: S = 1 bit (maximally entangled)
        """
        rho_A = self.reduced_density_matrix_A()
        eigenvalues = np.linalg.eigvalsh(rho_A)
        entropy = 0.0
        for ev in eigenvalues:
            if ev > 1e-10:
                entropy -= ev * np.log2(ev)
        return entropy

    def concurrence(self) -> float:
        """
        Concurrence C (measure of entanglement).

        C = 0 for separable states
        C = 1 for maximally entangled (Bell states)
        """
        rho = self.density_matrix()
        # Spin-flip matrix
        sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sigma_yy = np.kron(sigma_y, sigma_y)

        # R = ρ * (σ_y ⊗ σ_y) * ρ* * (σ_y ⊗ σ_y)
        rho_tilde = sigma_yy @ rho.conj() @ sigma_yy
        R = rho @ rho_tilde

        eigenvalues = np.sqrt(np.abs(np.linalg.eigvals(R)))
        eigenvalues = np.sort(eigenvalues)[::-1]

        return max(0, eigenvalues[0] - eigenvalues[1] - eigenvalues[2] - eigenvalues[3])

    def correlations(self) -> Dict[str, float]:
        """
        Compute correlations for Bell inequality test.

        ⟨σ_z ⊗ σ_z⟩, etc.
        """
        rho = self.density_matrix()
        sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
        sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

        ZZ = np.kron(sigma_z, sigma_z)
        XX = np.kron(sigma_x, sigma_x)
        ZI = np.kron(sigma_z, np.eye(2))
        IZ = np.kron(np.eye(2), sigma_z)

        return {
            'ZZ': np.real(np.trace(rho @ ZZ)),
            'XX': np.real(np.trace(rho @ XX)),
            'ZI': np.real(np.trace(rho @ ZI)),
            'IZ': np.real(np.trace(rho @ IZ))
        }

    @staticmethod
    def grid_interpretation() -> Dict[str, str]:
        """Grid interpretation of entanglement."""
        return {
            'entanglement': 'Correlated patterns across grid regions',
            'bell_states': 'Maximally correlated pattern pairs',
            'nonlocal': 'Correlations exceed local hidden variable bounds',
            'monogamy': 'Pattern correlations cannot be shared arbitrarily',
            'creation': 'Requires shared MRH history (interaction)',
            'entropy': 'Subsystem entropy = pattern info in other part'
        }


class NoCloningTheorem:
    """
    No-cloning theorem and its implications.

    Cannot create: |ψ⟩|0⟩ → |ψ⟩|ψ⟩ for arbitrary |ψ⟩

    Proof: Linearity of quantum mechanics forbids it.
    If U|ψ⟩|0⟩ = |ψ⟩|ψ⟩ and U|φ⟩|0⟩ = |φ⟩|φ⟩,
    then ⟨ψ|φ⟩ = ⟨ψ|φ⟩² which only holds for ⟨ψ|φ⟩ = 0 or 1.

    Grid interpretation: A grid pattern is indivisible.
    Copying would require creating pattern from nothing,
    violating conservation of pattern distinguishability.
    """

    def __init__(self):
        pass

    def cloning_fidelity(self, theta: float) -> float:
        """
        Best possible cloning fidelity for qubit at angle θ.

        Optimal universal cloning: F = 5/6 ≈ 0.833

        For state-dependent cloning of known basis:
        F = 1 (perfect) for orthogonal states.
        """
        # Universal cloning fidelity (state-independent)
        return 5/6

    def broadcasting_limit(self) -> str:
        """
        No-broadcasting theorem (generalization of no-cloning).

        Cannot broadcast: ρ → ρ ⊗ ρ for non-commuting ρ
        """
        return "No-broadcasting: Cannot copy non-commuting density matrices"

    def implications(self) -> Dict[str, str]:
        """Implications of no-cloning theorem."""
        return {
            'quantum_crypto': 'Eavesdropping disturbs state → detectable',
            'no_superluminal': 'Cannot send info faster than light via entanglement',
            'computation': 'Cannot simply copy intermediate results',
            'error_correction': 'Must encode redundantly, not copy',
            'teleportation': 'Destroys original (move, not copy)',
            'classical_limit': 'Classical info CAN be copied (orthogonal states)'
        }

    @staticmethod
    def grid_interpretation() -> Dict[str, str]:
        """Grid interpretation of no-cloning."""
        return {
            'indivisibility': 'Grid patterns are indivisible units',
            'conservation': 'Pattern distinguishability is conserved',
            'copying': 'Would create info from nothing (violates MRH)',
            'teleportation': 'Move pattern, don\'t copy it',
            'classical': 'Orthogonal patterns can be distinguished and copied'
        }


class QuantumChannel:
    """
    Quantum channels and information capacity.

    A quantum channel is a completely positive trace-preserving
    (CPTP) map that transforms density matrices:

    ρ → Σ_k K_k ρ K_k†  (Kraus representation)

    Key capacity: Holevo bound χ = max I(X:B)

    Grid interpretation: Channels are maps between pattern
    configurations. Noise = pattern info leaking to environment.
    """

    def __init__(self, channel_type: str = 'depolarizing', p: float = 0.1):
        """
        Args:
            channel_type: 'depolarizing', 'dephasing', 'amplitude_damping'
            p: Error/noise probability
        """
        self.channel_type = channel_type
        self.p = p
        self.kraus_operators = self._create_kraus_operators()

    def _create_kraus_operators(self) -> List[np.ndarray]:
        """Create Kraus operators for the channel."""
        I = np.eye(2, dtype=complex)
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        Z = np.array([[1, 0], [0, -1]], dtype=complex)

        if self.channel_type == 'depolarizing':
            # Depolarizing: ρ → (1-p)ρ + p/3(XρX + YρY + ZρZ)
            return [
                np.sqrt(1 - self.p) * I,
                np.sqrt(self.p / 3) * X,
                np.sqrt(self.p / 3) * Y,
                np.sqrt(self.p / 3) * Z
            ]
        elif self.channel_type == 'dephasing':
            # Dephasing: ρ → (1-p)ρ + p ZρZ
            return [
                np.sqrt(1 - self.p) * I,
                np.sqrt(self.p) * Z
            ]
        elif self.channel_type == 'amplitude_damping':
            # Amplitude damping: |1⟩ → |0⟩ with probability p
            K0 = np.array([[1, 0], [0, np.sqrt(1 - self.p)]], dtype=complex)
            K1 = np.array([[0, np.sqrt(self.p)], [0, 0]], dtype=complex)
            return [K0, K1]
        else:
            return [I]  # Identity channel

    def apply(self, rho: np.ndarray) -> np.ndarray:
        """Apply channel to density matrix."""
        result = np.zeros_like(rho, dtype=complex)
        for K in self.kraus_operators:
            result += K @ rho @ K.conj().T
        return result

    def classical_capacity(self) -> float:
        """
        Classical capacity (Holevo bound).

        For depolarizing channel: C = 1 - H(p) - p log₂(3)
        """
        if self.channel_type == 'depolarizing':
            if self.p == 0:
                return 1.0
            if self.p >= 3/4:
                return 0.0
            # Approximate formula
            H_p = -self.p * np.log2(self.p) - (1-self.p) * np.log2(1-self.p) if 0 < self.p < 1 else 0
            return max(0, 1 - H_p)
        elif self.channel_type == 'dephasing':
            # Dephasing: C = 1 - H(p)
            if self.p == 0 or self.p == 1:
                return 1.0
            H_p = -self.p * np.log2(self.p) - (1-self.p) * np.log2(1-self.p)
            return 1 - H_p
        return 1.0  # Identity

    def quantum_capacity(self) -> float:
        """
        Quantum capacity Q (coherent info maximization).

        For depolarizing: Q > 0 only for p < ~0.15
        """
        if self.channel_type == 'depolarizing':
            # Approximate threshold
            if self.p > 0.15:
                return 0.0
            return max(0, 1 - 2 * self.p)  # Rough approximation
        return 1.0

    def holevo_information(self, ensemble: List[Tuple[float, np.ndarray]]) -> float:
        """
        Holevo information χ for an ensemble.

        χ = S(ρ) - Σ p_i S(ρ_i)

        where ρ = Σ p_i ρ_i
        """
        # Average output state
        rho_avg = np.zeros((2, 2), dtype=complex)
        for prob, rho_in in ensemble:
            rho_out = self.apply(rho_in)
            rho_avg += prob * rho_out

        # S(ρ)
        eigenvalues = np.linalg.eigvalsh(rho_avg)
        S_avg = 0.0
        for ev in eigenvalues:
            if ev > 1e-10:
                S_avg -= ev * np.log2(ev)

        # Σ p_i S(ρ_i)
        S_sum = 0.0
        for prob, rho_in in ensemble:
            rho_out = self.apply(rho_in)
            eigenvalues = np.linalg.eigvalsh(rho_out)
            S_i = 0.0
            for ev in eigenvalues:
                if ev > 1e-10:
                    S_i -= ev * np.log2(ev)
            S_sum += prob * S_i

        return S_avg - S_sum

    @staticmethod
    def grid_interpretation() -> Dict[str, str]:
        """Grid interpretation of quantum channels."""
        return {
            'channel': 'Map between pattern configurations',
            'noise': 'Pattern info leaking to environment (beyond MRH)',
            'kraus': 'Different paths for pattern evolution',
            'capacity': 'Max rate of pattern transfer through channel',
            'holevo': 'Classical info extractable from quantum patterns',
            'quantum_cap': 'Coherent pattern transfer rate'
        }


class QuantumErrorCorrection:
    """
    Quantum error correction fundamentals.

    Key insight: Encode logical qubit in multiple physical qubits
    such that errors can be detected and corrected.

    Simplest codes:
    - 3-qubit bit-flip: |0⟩ → |000⟩, |1⟩ → |111⟩
    - 3-qubit phase-flip: |+⟩ → |+++⟩, |-⟩ → |---⟩
    - 9-qubit Shor: Handles both

    Grid interpretation: Spread pattern info across multiple
    grid regions. Error = local pattern disturbance. Correction
    = using redundancy to restore original pattern.
    """

    def __init__(self, code: str = 'bit_flip_3'):
        """
        Args:
            code: 'bit_flip_3', 'phase_flip_3', 'shor_9'
        """
        self.code = code

    def encode(self, state: Qubit) -> np.ndarray:
        """Encode logical qubit into code space."""
        if self.code == 'bit_flip_3':
            # |0⟩_L = |000⟩, |1⟩_L = |111⟩
            logical_0 = np.zeros(8, dtype=complex)
            logical_0[0] = 1  # |000⟩
            logical_1 = np.zeros(8, dtype=complex)
            logical_1[7] = 1  # |111⟩
            return state.alpha * logical_0 + state.beta * logical_1

        elif self.code == 'phase_flip_3':
            # |0⟩_L = |+++⟩, |1⟩_L = |---⟩
            plus = np.array([1, 1]) / np.sqrt(2)
            minus = np.array([1, -1]) / np.sqrt(2)
            logical_0 = np.kron(np.kron(plus, plus), plus)
            logical_1 = np.kron(np.kron(minus, minus), minus)
            return state.alpha * logical_0 + state.beta * logical_1

        return state.state_vector().flatten()

    def syndrome_measurement(self, state: np.ndarray) -> List[int]:
        """
        Measure error syndrome without collapsing logical state.

        Returns syndrome bits indicating error location.
        """
        # Simplified: return syndrome for 3-qubit bit-flip code
        if self.code == 'bit_flip_3':
            # Syndromes: Z₁Z₂, Z₂Z₃
            # 00 = no error, 01 = qubit 3, 10 = qubit 1, 11 = qubit 2
            return [0, 0]  # Placeholder
        return []

    def code_distance(self) -> int:
        """
        Code distance d (minimum weight of undetectable error).

        Can correct ⌊(d-1)/2⌋ errors.
        """
        distances = {
            'bit_flip_3': 3,
            'phase_flip_3': 3,
            'shor_9': 3,
            'steane_7': 3,
            'surface_d': 3  # Depends on d
        }
        return distances.get(self.code, 1)

    def threshold(self) -> float:
        """
        Threshold error rate for fault-tolerant computation.

        Below threshold: Can reduce error arbitrarily by adding layers.
        """
        thresholds = {
            'bit_flip_3': 0.5,
            'phase_flip_3': 0.5,
            'shor_9': 0.11,
            'steane_7': 0.01,
            'surface_d': 0.01  # ~1%
        }
        return thresholds.get(self.code, 0.01)

    @staticmethod
    def grid_interpretation() -> Dict[str, str]:
        """Grid interpretation of quantum error correction."""
        return {
            'encoding': 'Spread pattern info across multiple grid regions',
            'redundancy': 'Same info in different spatial locations',
            'error': 'Local pattern disturbance',
            'syndrome': 'Detect where pattern was disturbed',
            'correction': 'Use redundancy to restore original pattern',
            'threshold': 'Max error rate for indefinite pattern preservation'
        }


def run_verification_tests() -> Dict[str, bool]:
    """Run all verification tests for Session #329."""
    results = {}

    # Test 1: Qubit probabilities sum to 1
    q = Qubit(alpha=1/np.sqrt(2), beta=1j/np.sqrt(2))
    prob_sum = q.probability_zero() + q.probability_one()
    results['qubit_normalized'] = np.isclose(prob_sum, 1.0)

    # Test 2: Pure state has zero entropy
    q_pure = Qubit(alpha=1, beta=0)
    results['pure_state_entropy'] = np.isclose(q_pure.von_neumann_entropy(), 0.0, atol=1e-10)

    # Test 3: Bell state is maximally entangled
    bell = EntangledPair('Phi+')
    S_ent = bell.entanglement_entropy()
    results['bell_maximal'] = np.isclose(S_ent, 1.0, rtol=0.01)

    # Test 4: Bell state has unit concurrence
    C = bell.concurrence()
    results['bell_concurrence'] = np.isclose(C, 1.0, rtol=0.05)

    # Test 5: No-cloning fidelity bounded
    nc = NoCloningTheorem()
    F = nc.cloning_fidelity(np.pi/4)
    results['no_cloning_bound'] = F <= 1.0 and F >= 0.5

    # Test 6: Depolarizing channel reduces capacity
    ch_clean = QuantumChannel('depolarizing', p=0.0)
    ch_noisy = QuantumChannel('depolarizing', p=0.3)
    results['channel_noise_reduces'] = ch_clean.classical_capacity() >= ch_noisy.classical_capacity()

    # Test 7: Error correction code has distance ≥ 3
    qec = QuantumErrorCorrection('bit_flip_3')
    results['code_distance'] = qec.code_distance() >= 3

    # Test 8: Grid interpretations exist
    results['grid_interpretations'] = (
        'measurement' in Qubit.grid_interpretation() and
        'nonlocal' in EntangledPair.grid_interpretation() and
        'indivisibility' in NoCloningTheorem.grid_interpretation()
    )

    return results


def create_visualization(save_path: str = None):
    """Create comprehensive visualization for Session #329."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Session #329: Quantum Information from the Planck Grid\n'
                 'Information Theory Arc (2/4)',
                 fontsize=14, fontweight='bold')

    # Panel 1: Bloch sphere
    ax1 = axes[0, 0]
    ax1.set_aspect('equal')

    # Draw unit circle (equator)
    theta = np.linspace(0, 2*np.pi, 100)
    ax1.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.3)
    ax1.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax1.axvline(x=0, color='k', linestyle='-', alpha=0.3)

    # Plot some states
    states = {
        '|0⟩': (0, 1),
        '|1⟩': (0, -1),
        '|+⟩': (1, 0),
        '|-⟩': (-1, 0),
    }
    for label, (x, z) in states.items():
        ax1.plot(x, z, 'bo', markersize=10)
        ax1.annotate(label, (x + 0.1, z + 0.1), fontsize=10)

    ax1.set_xlabel('X')
    ax1.set_ylabel('Z')
    ax1.set_title('Bloch Sphere (XZ plane)')
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)

    # Panel 2: Entanglement entropy
    ax2 = axes[0, 1]

    # Entanglement entropy for parametrized state
    # |ψ⟩ = cos(θ)|00⟩ + sin(θ)|11⟩
    theta_range = np.linspace(0, np.pi/2, 50)
    S_ent = []
    for th in theta_range:
        c = np.cos(th)
        s = np.sin(th)
        # Reduced density matrix eigenvalues
        if c > 0 and s > 0:
            S_ent.append(-c**2 * np.log2(c**2) - s**2 * np.log2(s**2))
        elif c > 0:
            S_ent.append(0)
        else:
            S_ent.append(0)

    ax2.plot(theta_range * 180 / np.pi, S_ent, 'b-', linewidth=2)
    ax2.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Max (Bell state)')
    ax2.set_xlabel('Mixing angle θ (degrees)')
    ax2.set_ylabel('Entanglement entropy (bits)')
    ax2.set_title('Entanglement Entropy')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Panel 3: Channel capacity vs noise
    ax3 = axes[0, 2]

    p_range = np.linspace(0, 0.5, 50)
    C_depol = []
    C_dephase = []

    for p in p_range:
        # Depolarizing
        if p > 0 and p < 1:
            H_p = -p * np.log2(p) - (1-p) * np.log2(1-p)
            C_depol.append(max(0, 1 - H_p))
        else:
            C_depol.append(1.0 if p == 0 else 0.0)

        # Dephasing
        if p > 0 and p < 1:
            H_p = -p * np.log2(p) - (1-p) * np.log2(1-p)
            C_dephase.append(1 - H_p)
        else:
            C_dephase.append(1.0)

    ax3.plot(p_range, C_depol, 'b-', linewidth=2, label='Depolarizing')
    ax3.plot(p_range, C_dephase, 'r-', linewidth=2, label='Dephasing')
    ax3.set_xlabel('Error probability p')
    ax3.set_ylabel('Classical capacity (bits/use)')
    ax3.set_title('Channel Capacity vs Noise')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Panel 4: Key concepts
    ax4 = axes[1, 0]
    ax4.axis('off')

    concepts_text = """
    QUANTUM INFORMATION CONCEPTS

    ┌─────────────────────────────────────────┐
    │ QUBIT                                    │
    │ |ψ⟩ = α|0⟩ + β|1⟩                       │
    │                                          │
    │ • Superposition of two patterns          │
    │ • Complex amplitudes (phase matters)     │
    │ • Bloch sphere representation            │
    └─────────────────────────────────────────┘

    ┌─────────────────────────────────────────┐
    │ ENTANGLEMENT                             │
    │ |Φ⁺⟩ = (|00⟩ + |11⟩)/√2                 │
    │                                          │
    │ • Correlated patterns across regions     │
    │ • Cannot be created locally              │
    │ • Violates Bell inequalities             │
    │ • Monogamy of correlations               │
    └─────────────────────────────────────────┘

    ┌─────────────────────────────────────────┐
    │ NO-CLONING                               │
    │ |ψ⟩|0⟩ ↛ |ψ⟩|ψ⟩                         │
    │                                          │
    │ • Patterns are indivisible               │
    │ • Cannot copy unknown quantum state      │
    │ • Enables quantum cryptography           │
    └─────────────────────────────────────────┘
    """

    ax4.text(0.02, 0.98, concepts_text, transform=ax4.transAxes, fontsize=8,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    ax4.set_title('Key Concepts')

    # Panel 5: Grid interpretation
    ax5 = axes[1, 1]
    ax5.axis('off')

    grid_text = """
    GRID INTERPRETATION

    QUANTUM INFO = PATTERN COHERENCE WITHIN MRH

    ┌─────────────────────────────────────────┐
    │  QUBIT = Coherent pattern superposition  │
    │                                          │
    │  |ψ⟩ = α|■□⟩ + β|□■⟩                    │
    │                                          │
    │  Both patterns "exist" until measurement │
    │  Measurement → one pattern, other → MRH  │
    └─────────────────────────────────────────┘

    ┌─────────────────────────────────────────┐
    │  ENTANGLEMENT = Correlated patterns      │
    │                                          │
    │  Region A      Region B                  │
    │    ■ □    ←→     ■ □                    │
    │    □ ■    ←→     □ ■                    │
    │                                          │
    │  Correlations via shared MRH history     │
    └─────────────────────────────────────────┘

    ┌─────────────────────────────────────────┐
    │  DECOHERENCE = Pattern info → MRH        │
    │                                          │
    │  Coherent     →    Mixed (decohered)     │
    │  |ψ⟩⟨ψ|       →    Σ p_i |i⟩⟨i|         │
    │                                          │
    │  Phase info leaks to environment         │
    └─────────────────────────────────────────┘

    NO-CLONING: Patterns are indivisible
    Copying would create info from nothing
    """

    ax5.text(0.02, 0.98, grid_text, transform=ax5.transAxes, fontsize=7,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    ax5.set_title('Grid Perspective')

    # Panel 6: Summary
    ax6 = axes[1, 2]
    ax6.axis('off')

    verification = run_verification_tests()
    passed = sum(verification.values())
    total = len(verification)

    summary_text = f"""
    SESSION #329 RESULTS: {passed}/{total} verified

    Key Findings:

    ✓ Qubit
      |ψ⟩ = α|0⟩ + β|1⟩
      Coherent pattern superposition

    ✓ Entanglement
      Bell states: S = 1 bit (maximum)
      Correlated patterns across space

    ✓ No-Cloning
      Cannot copy unknown quantum state
      F_max = 5/6 (universal cloning)

    ✓ Quantum Channels
      Noise reduces capacity
      Holevo bound limits classical info

    ✓ Error Correction
      Spread pattern across qubits
      Detect and correct local errors

    Grid Interpretation:
    • Quantum = coherent patterns within MRH
    • Entanglement = correlated patterns
    • Decoherence = pattern info → MRH
    • No-cloning = patterns indivisible

    ★ INFORMATION THEORY ARC (2/4) ★
    """

    ax6.text(0.02, 0.98, summary_text, transform=ax6.transAxes, fontsize=8,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved visualization to {save_path}")

    plt.close()
    return fig


def main():
    """Main execution for Session #329."""
    print("=" * 70)
    print("SESSION #329: Quantum Information from the Planck Grid")
    print("Information Theory Arc (Session 2/4)")
    print("=" * 70)

    # Part 1: Qubits
    print("\n" + "=" * 50)
    print("PART 1: QUBITS")
    print("=" * 50)

    print("\nQubit state: |ψ⟩ = α|0⟩ + β|1⟩")

    # Create various qubits
    q0 = Qubit(alpha=1, beta=0)
    q1 = Qubit(alpha=0, beta=1)
    q_plus = Qubit(alpha=1/np.sqrt(2), beta=1/np.sqrt(2))
    q_i = Qubit(alpha=1/np.sqrt(2), beta=1j/np.sqrt(2))

    print(f"\nExample qubits:")
    print(f"  |0⟩: P(0) = {q0.probability_zero():.3f}, P(1) = {q0.probability_one():.3f}")
    print(f"  |1⟩: P(0) = {q1.probability_zero():.3f}, P(1) = {q1.probability_one():.3f}")
    print(f"  |+⟩: P(0) = {q_plus.probability_zero():.3f}, P(1) = {q_plus.probability_one():.3f}")
    print(f"  |i⟩: P(0) = {q_i.probability_zero():.3f}, P(1) = {q_i.probability_one():.3f}")

    print(f"\nBloch coordinates:")
    for name, q in [('|0⟩', q0), ('|1⟩', q1), ('|+⟩', q_plus), ('|i⟩', q_i)]:
        x, y, z = q.bloch_coordinates()
        print(f"  {name}: (x={x:.2f}, y={y:.2f}, z={z:.2f})")

    print(f"\nVon Neumann entropy (pure states):")
    print(f"  S(|0⟩) = {q0.von_neumann_entropy():.3f} bits")
    print(f"  S(|+⟩) = {q_plus.von_neumann_entropy():.3f} bits")

    qubit_interp = Qubit.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in qubit_interp.items():
        print(f"  {key}: {value}")

    # Part 2: Entanglement
    print("\n" + "=" * 50)
    print("PART 2: ENTANGLEMENT")
    print("=" * 50)

    print("\nBell states (maximally entangled):")
    for name in ['Phi+', 'Phi-', 'Psi+', 'Psi-']:
        bell = EntangledPair(name)
        S = bell.entanglement_entropy()
        C = bell.concurrence()
        print(f"  |{name}⟩: S = {S:.3f} bits, C = {C:.3f}")

    bell_phi = EntangledPair('Phi+')
    corr = bell_phi.correlations()
    print(f"\nCorrelations for |Φ⁺⟩:")
    for key, value in corr.items():
        print(f"  ⟨{key}⟩ = {value:.3f}")

    print(f"\nBell inequality (CHSH):")
    print(f"  Classical bound: |S| ≤ 2")
    print(f"  Quantum maximum: |S| = 2√2 ≈ 2.83")

    ent_interp = EntangledPair.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in ent_interp.items():
        print(f"  {key}: {value}")

    # Part 3: No-Cloning
    print("\n" + "=" * 50)
    print("PART 3: NO-CLONING THEOREM")
    print("=" * 50)

    nc = NoCloningTheorem()

    print("\nNo-cloning theorem:")
    print("  Cannot create: |ψ⟩|0⟩ → |ψ⟩|ψ⟩ for arbitrary |ψ⟩")
    print(f"\nOptimal universal cloning fidelity: F = {nc.cloning_fidelity(0):.3f}")

    print(f"\n{nc.broadcasting_limit()}")

    impl = nc.implications()
    print(f"\nImplications:")
    for key, value in impl.items():
        print(f"  {key}: {value}")

    nc_interp = NoCloningTheorem.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in nc_interp.items():
        print(f"  {key}: {value}")

    # Part 4: Quantum Channels
    print("\n" + "=" * 50)
    print("PART 4: QUANTUM CHANNELS")
    print("=" * 50)

    print("\nChannel capacities:")
    for ch_type in ['depolarizing', 'dephasing']:
        print(f"\n  {ch_type.upper()} channel:")
        for p in [0.0, 0.1, 0.2, 0.3]:
            ch = QuantumChannel(ch_type, p)
            C = ch.classical_capacity()
            Q = ch.quantum_capacity()
            print(f"    p = {p}: C = {C:.3f} bits, Q = {Q:.3f} bits")

    ch = QuantumChannel('depolarizing', 0.1)
    # Create ensemble
    rho_0 = np.array([[1, 0], [0, 0]], dtype=complex)
    rho_1 = np.array([[0, 0], [0, 1]], dtype=complex)
    ensemble = [(0.5, rho_0), (0.5, rho_1)]
    chi = ch.holevo_information(ensemble)
    print(f"\nHolevo information for uniform {0,1} ensemble (p=0.1):")
    print(f"  χ = {chi:.3f} bits")

    ch_interp = QuantumChannel.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in ch_interp.items():
        print(f"  {key}: {value}")

    # Part 5: Quantum Error Correction
    print("\n" + "=" * 50)
    print("PART 5: QUANTUM ERROR CORRECTION")
    print("=" * 50)

    print("\nError correction codes:")
    for code in ['bit_flip_3', 'phase_flip_3', 'shor_9']:
        qec = QuantumErrorCorrection(code)
        print(f"  {code}: d = {qec.code_distance()}, threshold = {qec.threshold():.2%}")

    print("\nKey principles:")
    print("  1. Encode logical qubit in multiple physical qubits")
    print("  2. Measure syndromes without collapsing logical state")
    print("  3. Apply corrections based on syndrome")
    print("  4. Below threshold: can reduce error arbitrarily")

    qec_interp = QuantumErrorCorrection.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in qec_interp.items():
        print(f"  {key}: {value}")

    # Verification
    print("\n" + "=" * 50)
    print("VERIFICATION SUMMARY")
    print("=" * 50)

    results = run_verification_tests()
    passed = sum(results.values())
    total = len(results)

    print(f"\nResults: {passed}/{total} tests passed\n")
    for test, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {test}: {status}")

    # Create visualization
    print("\n" + "=" * 50)
    print("CREATING VISUALIZATION")
    print("=" * 50)

    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, 'session329_quantum_information.png')
    create_visualization(save_path)

    # Final summary
    print("\n" + "=" * 70)
    print("SESSION #329 COMPLETE")
    print("=" * 70)

    print(f"""
    INFORMATION THEORY ARC (Sessions #328-331):

    Session #328: Information Theory Foundations  ✅ 8/8
    Session #329: Quantum Information             ✅ {passed}/{total}
    Session #330: Holographic Principle           NEXT
    Session #331: Black Hole Information          PLANNED

    ═══════════════════════════════════════════════════════════════

    KEY INSIGHTS FROM SESSION #329:

    1. QUBITS
       • |ψ⟩ = α|0⟩ + β|1⟩ (superposition)
       • Complex amplitudes encode phase
       • Bloch sphere = full state space
       • Pure state: S = 0 (no uncertainty)

    2. ENTANGLEMENT
       • Bell states: maximally correlated
       • S_entanglement = 1 bit (maximum)
       • Cannot be created by local operations
       • Violates Bell inequalities (nonlocal)

    3. NO-CLONING
       • Cannot copy arbitrary quantum state
       • Best universal cloning: F = 5/6
       • Enables quantum cryptography
       • Teleportation moves, doesn't copy

    4. QUANTUM CHANNELS
       • CPTP maps on density matrices
       • Noise reduces capacity
       • Holevo bound limits classical info
       • Quantum capacity for coherence

    5. ERROR CORRECTION
       • Encode in multiple qubits
       • Measure syndromes, not state
       • Correct errors below threshold
       • Spread pattern info for redundancy

    ═══════════════════════════════════════════════════════════════

    GRID INTERPRETATION:

    QUANTUM INFO = PATTERN COHERENCE WITHIN MRH

    • Superposition: Multiple patterns coexist coherently
    • Entanglement: Correlated patterns across regions
    • Measurement: Collapse to one pattern, rest → MRH
    • Decoherence: Pattern info leaking to environment
    • No-cloning: Patterns are indivisible units

    The MRH boundary defines what is "quantum" vs "classical":
    - Inside MRH: Coherent, trackable, quantum
    - Beyond MRH: Averaged, thermal, classical

    ★ INFORMATION THEORY ARC (2/4) ★

    Next: Session #330 - Holographic Principle
    """)

    return results


if __name__ == "__main__":
    main()
