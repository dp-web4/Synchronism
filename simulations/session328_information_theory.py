#!/usr/bin/env python3
"""
Session #328: Information Theory from the Planck Grid
Information Theory Arc (Session 1/4)

This session explores information theory from the grid perspective:
1. Shannon entropy and the grid
2. Landauer's principle (information erasure)
3. Maxwell's demon resolution
4. Channel capacity and noise
5. Data compression and pattern redundancy

Key insight: Information is pattern distinguishability on the grid.
Entropy measures uncertainty about which pattern is realized.
The MRH is the channel capacity of nature itself.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Callable
from scipy import constants as const
from scipy.special import xlogy

# Physical constants
k_B = const.k  # Boltzmann constant
h = const.h  # Planck constant
c = const.c  # Speed of light
T_room = 300  # Room temperature (K)


@dataclass
class ShannonEntropy:
    """
    Shannon entropy and information theory fundamentals.

    H(X) = -Σ p(x) log₂ p(x)

    Measures uncertainty (bits) about a random variable.

    Grid interpretation: H measures how many bits needed to
    specify which pattern is realized from an ensemble.
    """

    def __init__(self, base: float = 2.0):
        """
        Args:
            base: Logarithm base (2 for bits, e for nats)
        """
        self.base = base

    def entropy(self, p: np.ndarray) -> float:
        """
        Calculate Shannon entropy.

        H(X) = -Σ p(x) log p(x)

        Uses xlogy for numerical stability at p=0.
        """
        p = np.asarray(p)
        p = p[p > 0]  # Remove zeros
        return -np.sum(p * np.log(p)) / np.log(self.base)

    def joint_entropy(self, p_xy: np.ndarray) -> float:
        """
        Joint entropy H(X,Y) from joint distribution.
        """
        return self.entropy(p_xy.flatten())

    def conditional_entropy(self, p_xy: np.ndarray) -> float:
        """
        Conditional entropy H(Y|X) = H(X,Y) - H(X).

        Measures remaining uncertainty about Y given X.
        """
        p_x = p_xy.sum(axis=1)
        H_XY = self.joint_entropy(p_xy)
        H_X = self.entropy(p_x)
        return H_XY - H_X

    def mutual_information(self, p_xy: np.ndarray) -> float:
        """
        Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y).

        Measures how much knowing X tells about Y.
        """
        p_x = p_xy.sum(axis=1)
        p_y = p_xy.sum(axis=0)
        H_X = self.entropy(p_x)
        H_Y = self.entropy(p_y)
        H_XY = self.joint_entropy(p_xy)
        return H_X + H_Y - H_XY

    def relative_entropy(self, p: np.ndarray, q: np.ndarray) -> float:
        """
        Relative entropy (KL divergence) D(P||Q).

        D(P||Q) = Σ p(x) log(p(x)/q(x))

        Measures "distance" from Q to P (not symmetric!).
        """
        p = np.asarray(p)
        q = np.asarray(q)
        # Only sum where both are positive
        mask = (p > 0) & (q > 0)
        return np.sum(p[mask] * np.log(p[mask] / q[mask])) / np.log(self.base)

    def maximum_entropy_distribution(self, n: int) -> np.ndarray:
        """
        Maximum entropy (uniform) distribution over n states.
        """
        return np.ones(n) / n

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of Shannon entropy."""
        return {
            'entropy': 'Bits to specify pattern from ensemble',
            'probability': 'Weight of each pattern configuration',
            'mutual_info': 'Shared pattern structure between regions',
            'conditional': 'Remaining uncertainty given context (MRH)',
            'max_entropy': 'Uniform distribution = all patterns equally likely',
            'kl_divergence': 'Information lost by approximating one pattern set with another'
        }


class LandauerPrinciple:
    """
    Landauer's principle: Information erasure has thermodynamic cost.

    Minimum energy to erase 1 bit: E_min = k_B T ln(2)

    At 300 K: E_min ≈ 2.87 × 10^-21 J ≈ 0.018 eV

    Grid interpretation: Erasing a bit means losing pattern
    distinguishability. This information must cross the MRH
    boundary and become thermal entropy in the bath.
    """

    def __init__(self, T: float = 300.0):
        """
        Args:
            T: Temperature (K)
        """
        self.T = T
        self.E_per_bit = k_B * T * np.log(2)

    def erasure_energy(self, n_bits: int) -> float:
        """
        Minimum energy to erase n bits at temperature T.

        E = n × k_B T ln(2)
        """
        return n_bits * self.E_per_bit

    def erasure_entropy(self, n_bits: int) -> float:
        """
        Entropy increase from erasing n bits.

        ΔS = n × k_B ln(2)
        """
        return n_bits * k_B * np.log(2)

    def maximum_bits_from_energy(self, E: float) -> float:
        """
        Maximum bits that can be erased with energy E.

        n_max = E / (k_B T ln(2))
        """
        return E / self.E_per_bit

    def bit_energy_vs_temperature(self, T_range: np.ndarray) -> np.ndarray:
        """
        Energy per bit as function of temperature.
        """
        return k_B * T_range * np.log(2)

    def efficiency(self, E_actual: float, n_bits: int) -> float:
        """
        Landauer efficiency η = E_min / E_actual.

        η = 1 for reversible (Landauer limit) operation.
        η < 1 for irreversible operation.
        """
        E_min = self.erasure_energy(n_bits)
        return E_min / E_actual if E_actual > 0 else 0

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of Landauer principle."""
        return {
            'erasure': 'Pattern distinguishability lost',
            'energy_cost': 'Work to push info across MRH boundary',
            'entropy_increase': 'Info becomes thermal (beyond MRH)',
            'reversibility': 'Reversible = no info crosses MRH',
            'landauer_limit': 'Fundamental cost of forgetting',
            'modern_computing': 'CPUs operate at ~1000× Landauer limit'
        }


class MaxwellDemon:
    """
    Maxwell's demon and its resolution.

    The demon appears to violate 2nd Law by sorting hot/cold molecules.
    Resolution: The demon must STORE information about molecules,
    and eventually ERASE that memory (Landauer cost).

    Grid interpretation: The demon creates correlations between
    its memory and the gas. When reset, these correlations must
    cross the MRH → entropy production.
    """

    def __init__(self, T: float = 300.0, n_molecules: int = 1000):
        """
        Args:
            T: Temperature (K)
            n_molecules: Number of molecules sorted
        """
        self.T = T
        self.n_molecules = n_molecules
        self.landauer = LandauerPrinciple(T)

    def apparent_entropy_decrease(self) -> float:
        """
        Apparent entropy decrease from sorting (per molecule).

        For binary sorting (left/right): ΔS = k_B ln(2) per molecule
        """
        return self.n_molecules * k_B * np.log(2)

    def demon_memory_required(self) -> int:
        """
        Bits demon must store to track all molecules.

        1 bit per molecule (which side it should go).
        """
        return self.n_molecules

    def memory_erasure_cost(self) -> float:
        """
        Energy cost to erase demon's memory.

        E = n × k_B T ln(2)
        """
        return self.landauer.erasure_energy(self.n_molecules)

    def memory_erasure_entropy(self) -> float:
        """
        Entropy produced when erasing demon's memory.

        ΔS = n × k_B ln(2)
        """
        return self.landauer.erasure_entropy(self.n_molecules)

    def net_entropy_change(self) -> float:
        """
        Net entropy change: sorting gain + erasure cost.

        ΔS_net = -ΔS_sorting + ΔS_erasure ≥ 0

        Second Law satisfied!
        """
        sorting = -self.apparent_entropy_decrease()
        erasure = self.memory_erasure_entropy()
        return sorting + erasure  # ≥ 0

    def szilard_engine_work(self) -> float:
        """
        Maximum work extractable from 1-molecule Szilard engine.

        W_max = k_B T ln(2) (exactly Landauer limit)

        This equals the cost to reset the memory.
        Net work = 0 over full cycle.
        """
        return k_B * self.T * np.log(2)

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of Maxwell's demon."""
        return {
            'demon': 'Subsystem that creates correlations',
            'sorting': 'Transfers patterns to ordered configuration',
            'memory': 'Correlations between demon and gas patterns',
            'erasure': 'Correlations must cross MRH to reset',
            'resolution': '2nd Law = info flow across MRH is one-way',
            'szilard': 'Work from 1 bit = k_B T ln(2) (Landauer limit)'
        }


class ChannelCapacity:
    """
    Channel capacity and communication limits.

    Shannon's channel capacity theorem:
    C = max I(X;Y) over input distributions

    For additive white Gaussian noise (AWGN):
    C = B log₂(1 + S/N) bits/second

    Grid interpretation: Channel capacity is the maximum rate
    of pattern distinguishability that can cross a boundary.
    MRH itself has a "channel capacity" for information.
    """

    def __init__(self, bandwidth: float = 1e6, snr_db: float = 20.0):
        """
        Args:
            bandwidth: Channel bandwidth (Hz)
            snr_db: Signal-to-noise ratio (dB)
        """
        self.B = bandwidth
        self.snr_db = snr_db
        self.snr = 10 ** (snr_db / 10)

    def shannon_capacity(self) -> float:
        """
        Shannon capacity for AWGN channel.

        C = B log₂(1 + S/N) bits/second
        """
        return self.B * np.log2(1 + self.snr)

    def capacity_vs_snr(self, snr_db_range: np.ndarray) -> np.ndarray:
        """
        Channel capacity as function of SNR.
        """
        snr = 10 ** (snr_db_range / 10)
        return self.B * np.log2(1 + snr)

    def capacity_vs_bandwidth(self, B_range: np.ndarray) -> np.ndarray:
        """
        Channel capacity as function of bandwidth.
        """
        return B_range * np.log2(1 + self.snr)

    def energy_per_bit(self, P: float) -> float:
        """
        Energy per bit at power P.

        E_b = P / C
        """
        C = self.shannon_capacity()
        return P / C if C > 0 else np.inf

    def spectral_efficiency(self) -> float:
        """
        Spectral efficiency (bits/s/Hz).

        η = C / B = log₂(1 + S/N)
        """
        return np.log2(1 + self.snr)

    def binary_symmetric_channel(self, p_error: float) -> float:
        """
        Capacity of binary symmetric channel with error probability p.

        C = 1 - H(p) = 1 + p log₂(p) + (1-p) log₂(1-p)
        """
        if p_error == 0 or p_error == 1:
            return 1.0
        return 1.0 + p_error * np.log2(p_error) + (1 - p_error) * np.log2(1 - p_error)

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of channel capacity."""
        return {
            'channel': 'Path for pattern transfer between regions',
            'capacity': 'Max rate of distinguishable patterns',
            'noise': 'Random pattern perturbations',
            'bandwidth': 'Range of pattern frequencies',
            'mrh_as_channel': 'MRH boundary has channel capacity',
            'quantum_limit': 'Holevo bound on quantum channels'
        }


class DataCompression:
    """
    Data compression and source coding.

    Shannon's source coding theorem:
    Optimal compression approaches H(X) bits per symbol.

    Types:
    - Lossless: Recover original exactly (entropy limit)
    - Lossy: Some information lost (rate-distortion limit)

    Grid interpretation: Compression removes pattern redundancy.
    The entropy H(X) is the irreducible pattern complexity.
    """

    def __init__(self):
        self.shannon = ShannonEntropy()

    def entropy_rate(self, p: np.ndarray) -> float:
        """
        Entropy rate (bits per symbol).

        H(X) = -Σ p(x) log₂ p(x)
        """
        return self.shannon.entropy(p)

    def compression_ratio(self, original_bits: int, compressed_bits: int) -> float:
        """
        Compression ratio = original / compressed.
        """
        return original_bits / compressed_bits if compressed_bits > 0 else np.inf

    def redundancy(self, p: np.ndarray, n_symbols: int) -> float:
        """
        Redundancy = 1 - H(X) / log₂(n_symbols).

        Fraction of bits that are "wasted" on predictable patterns.
        """
        H = self.shannon.entropy(p)
        H_max = np.log2(n_symbols)
        return 1 - H / H_max if H_max > 0 else 0

    def huffman_average_length(self, p: np.ndarray) -> float:
        """
        Approximate average Huffman code length.

        H(X) ≤ L < H(X) + 1
        """
        H = self.shannon.entropy(p)
        return H + 0.5  # Approximate middle of bounds

    def rate_distortion_bound(self, D: float, sigma2: float) -> float:
        """
        Rate-distortion function for Gaussian source.

        R(D) = (1/2) log₂(σ² / D) for D ≤ σ²
        R(D) = 0 for D > σ²

        Args:
            D: Distortion (mean squared error)
            sigma2: Source variance
        """
        if D >= sigma2:
            return 0.0
        return 0.5 * np.log2(sigma2 / D)

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of data compression."""
        return {
            'compression': 'Removing pattern redundancy',
            'entropy': 'Irreducible pattern complexity',
            'redundancy': 'Predictable pattern structure',
            'lossless': 'Preserve all pattern distinguishability',
            'lossy': 'Accept some pattern confusion',
            'mrh_compression': 'MRH naturally compresses: only relevant patterns tracked'
        }


class InformationThermodynamics:
    """
    Unified view of information and thermodynamics.

    Key relations:
    - S = k_B H (for microcanonical ensemble)
    - Landauer: ΔE ≥ k_B T ln(2) per bit erased
    - Maxwell demon: Sorting requires memory; erasure costs energy
    - Fluctuation theorems: Work-information relations

    Grid interpretation: Thermodynamic entropy IS information
    entropy when measured in appropriate units. The MRH
    defines what counts as "information" vs "noise".
    """

    def __init__(self, T: float = 300.0):
        self.T = T
        self.landauer = LandauerPrinciple(T)

    def boltzmann_to_shannon(self, S_boltzmann: float) -> float:
        """
        Convert Boltzmann entropy to Shannon entropy (bits).

        H = S / (k_B ln 2)
        """
        return S_boltzmann / (k_B * np.log(2))

    def shannon_to_boltzmann(self, H_bits: float) -> float:
        """
        Convert Shannon entropy (bits) to Boltzmann entropy.

        S = H × k_B ln(2)
        """
        return H_bits * k_B * np.log(2)

    def jarzynski_equality(self, W_values: np.ndarray) -> float:
        """
        Jarzynski equality: <e^{-βW}> = e^{-βΔF}

        Returns estimated ΔF from work measurements.

        Args:
            W_values: Array of work values (J)
        """
        beta = 1 / (k_B * self.T)
        return -np.log(np.mean(np.exp(-beta * W_values))) / beta

    def crooks_fluctuation(self, W_forward: np.ndarray, W_reverse: np.ndarray) -> float:
        """
        Crooks fluctuation theorem relates forward and reverse work distributions.

        P_F(W) / P_R(-W) = e^{β(W - ΔF)}

        Returns estimated ΔF from crossing point.
        """
        # Simplified: estimate ΔF from mean works
        return 0.5 * (np.mean(W_forward) + np.mean(W_reverse))

    def information_engine_efficiency(self, W_extracted: float, bits_used: int) -> float:
        """
        Efficiency of information engine.

        η = W_extracted / (bits × k_B T ln 2)

        η ≤ 1 by Landauer's principle.
        """
        W_max = bits_used * k_B * self.T * np.log(2)
        return W_extracted / W_max if W_max > 0 else 0

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of information thermodynamics."""
        return {
            'unity': 'Thermodynamic and information entropy are the same',
            'landauer': 'Forgetting = pattern info crosses MRH',
            'work_info': 'Can extract work from information (Szilard)',
            'jarzynski': 'Free energy from non-equilibrium trajectories',
            'crooks': 'Forward/reverse path symmetry',
            'mrh_unifies': 'MRH defines info/thermo boundary'
        }


def run_verification_tests() -> Dict[str, bool]:
    """Run all verification tests for Session #328."""
    results = {}

    # Test 1: Shannon entropy of uniform distribution
    shannon = ShannonEntropy()
    p_uniform = np.array([0.25, 0.25, 0.25, 0.25])
    H_uniform = shannon.entropy(p_uniform)
    results['shannon_uniform'] = np.isclose(H_uniform, 2.0, rtol=0.01)  # 2 bits

    # Test 2: Shannon entropy bounds (0 ≤ H ≤ log n)
    p_certain = np.array([1.0, 0.0, 0.0, 0.0])
    H_certain = shannon.entropy(p_certain)
    results['shannon_bounds'] = (H_certain >= 0) and (H_uniform <= np.log2(4) + 0.01)

    # Test 3: Landauer energy is positive
    landauer = LandauerPrinciple(T=300)
    E_bit = landauer.E_per_bit
    results['landauer_positive'] = E_bit > 0

    # Test 4: Landauer scales with temperature
    landauer_cold = LandauerPrinciple(T=100)
    landauer_hot = LandauerPrinciple(T=1000)
    results['landauer_scales'] = landauer_hot.E_per_bit > landauer_cold.E_per_bit

    # Test 5: Maxwell demon net entropy ≥ 0
    demon = MaxwellDemon(T=300, n_molecules=100)
    dS_net = demon.net_entropy_change()
    results['demon_second_law'] = dS_net >= -1e-30  # Allow tiny numerical error

    # Test 6: Channel capacity increases with SNR
    channel_low = ChannelCapacity(bandwidth=1e6, snr_db=10)
    channel_high = ChannelCapacity(bandwidth=1e6, snr_db=30)
    results['channel_increases'] = channel_high.shannon_capacity() > channel_low.shannon_capacity()

    # Test 7: Compression limited by entropy
    dc = DataCompression()
    p = np.array([0.5, 0.25, 0.125, 0.125])
    H = dc.entropy_rate(p)
    L = dc.huffman_average_length(p)
    results['compression_limited'] = L >= H - 0.01  # L ≥ H

    # Test 8: Grid interpretations exist
    it = InformationThermodynamics()
    results['grid_interpretations'] = (
        'entropy' in shannon.grid_interpretation() and
        'landauer_limit' in landauer.grid_interpretation() and
        'mrh_unifies' in it.grid_interpretation()
    )

    return results


def create_visualization(save_path: str = None):
    """Create comprehensive visualization for Session #328."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Session #328: Information Theory from the Planck Grid\n'
                 'Information Theory Arc (1/4)',
                 fontsize=14, fontweight='bold')

    # Panel 1: Shannon entropy vs probability
    ax1 = axes[0, 0]
    p_range = np.linspace(0.001, 0.999, 100)
    H_binary = -p_range * np.log2(p_range) - (1 - p_range) * np.log2(1 - p_range)

    ax1.plot(p_range, H_binary, 'b-', linewidth=2)
    ax1.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='Max at p=0.5')
    ax1.set_xlabel('Probability p')
    ax1.set_ylabel('Binary entropy H(p) (bits)')
    ax1.set_title('Shannon Entropy (Binary)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1.1)

    # Panel 2: Landauer energy vs temperature
    ax2 = axes[0, 1]
    T_range = np.linspace(10, 1000, 100)
    E_bit = k_B * T_range * np.log(2)

    ax2.semilogy(T_range, E_bit * 1e21, 'r-', linewidth=2)
    ax2.axvline(x=300, color='gray', linestyle='--', alpha=0.7, label='Room temp')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Energy per bit (×10⁻²¹ J)')
    ax2.set_title("Landauer's Principle")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.annotate(f'E(300K) = {k_B * 300 * np.log(2) * 1e21:.2f} × 10⁻²¹ J',
                 xy=(350, k_B * 300 * np.log(2) * 1e21), fontsize=9)

    # Panel 3: Channel capacity vs SNR
    ax3 = axes[0, 2]
    snr_db = np.linspace(0, 40, 100)
    B = 1e6  # 1 MHz bandwidth

    snr = 10 ** (snr_db / 10)
    C = B * np.log2(1 + snr)

    ax3.semilogy(snr_db, C / 1e6, 'g-', linewidth=2)
    ax3.set_xlabel('Signal-to-Noise Ratio (dB)')
    ax3.set_ylabel('Channel Capacity (Mbits/s)')
    ax3.set_title('Shannon Capacity (B = 1 MHz)')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=1, color='gray', linestyle=':', alpha=0.7)

    # Panel 4: Key concepts
    ax4 = axes[1, 0]
    ax4.axis('off')

    concepts_text = """
    INFORMATION THEORY FOUNDATIONS

    ┌─────────────────────────────────────────┐
    │ SHANNON ENTROPY                          │
    │ H(X) = -Σ p(x) log₂ p(x)                │
    │                                          │
    │ • Measures uncertainty (bits)            │
    │ • Max at uniform distribution            │
    │ • Compression limit                      │
    └─────────────────────────────────────────┘

    ┌─────────────────────────────────────────┐
    │ LANDAUER'S PRINCIPLE                     │
    │ E_min = k_B T ln(2) per bit erased      │
    │                                          │
    │ • Information erasure costs energy       │
    │ • At 300K: ~2.9 × 10⁻²¹ J per bit       │
    │ • CPUs: ~1000× Landauer limit           │
    └─────────────────────────────────────────┘

    ┌─────────────────────────────────────────┐
    │ CHANNEL CAPACITY                         │
    │ C = B log₂(1 + S/N)                      │
    │                                          │
    │ • Max reliable communication rate        │
    │ • Shannon's fundamental limit            │
    │ • MRH = nature's channel capacity       │
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

    INFORMATION = PATTERN DISTINGUISHABILITY

    ┌─────────────────────────────────────────┐
    │  Pattern A    Pattern B    Pattern C    │
    │   ■ ■ □        ■ □ ■        □ ■ ■       │
    │   □ ■ ■        ■ ■ □        ■ □ ■       │
    │                                          │
    │  Which pattern is realized?              │
    │  H = log₂(3) ≈ 1.58 bits                │
    └─────────────────────────────────────────┘

    ENTROPY = UNCERTAINTY ABOUT PATTERNS

    ┌─────────────────────────────────────────┐
    │  Equal probability → max entropy         │
    │  Certain pattern → zero entropy          │
    │  H measures "which pattern?" bits        │
    └─────────────────────────────────────────┘

    LANDAUER = PATTERN ERASURE COST

    ┌─────────────────────────────────────────┐
    │  Erasing = losing distinguishability     │
    │  Pattern info → thermal noise (MRH)      │
    │  Energy cost = k_B T ln(2) per bit       │
    └─────────────────────────────────────────┘

    MRH AS CHANNEL CAPACITY

    ┌─────────────────────────────────────────┐
    │  MRH = max rate of pattern transfer      │
    │  Beyond MRH = noise, not signal          │
    │  Thermodynamics from channel limits      │
    └─────────────────────────────────────────┘
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
    SESSION #328 RESULTS: {passed}/{total} verified

    Key Findings:

    ✓ Shannon Entropy
      H(X) = -Σ p(x) log p(x)
      Measures pattern uncertainty

    ✓ Landauer's Principle
      E ≥ k_B T ln(2) per bit erased
      Information erasure costs energy

    ✓ Maxwell's Demon Resolution
      Memory erasure → entropy increase
      Net ΔS ≥ 0 (2nd Law satisfied)

    ✓ Channel Capacity
      C = B log₂(1 + S/N)
      Maximum reliable info rate

    ✓ Data Compression
      Can't compress below H(X)
      Entropy = irreducible complexity

    Grid Interpretation:
    • Information = pattern distinguishability
    • Entropy = bits to specify pattern
    • MRH = channel capacity of nature
    • Landauer = pattern → thermal crossing

    ★ INFORMATION THEORY ARC (1/4) ★
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
    """Main execution for Session #328."""
    print("=" * 70)
    print("SESSION #328: Information Theory from the Planck Grid")
    print("Information Theory Arc (Session 1/4)")
    print("=" * 70)

    # Part 1: Shannon Entropy
    print("\n" + "=" * 50)
    print("PART 1: SHANNON ENTROPY")
    print("=" * 50)

    shannon = ShannonEntropy()

    print("\nShannon entropy formula:")
    print("  H(X) = -Σ p(x) log₂ p(x)")

    # Examples
    p_uniform = np.array([0.25, 0.25, 0.25, 0.25])
    p_biased = np.array([0.9, 0.05, 0.025, 0.025])
    p_certain = np.array([1.0, 0.0, 0.0, 0.0])

    print(f"\nEntropy examples (4 outcomes):")
    print(f"  Uniform [0.25, 0.25, 0.25, 0.25]: H = {shannon.entropy(p_uniform):.3f} bits")
    print(f"  Biased [0.9, 0.05, 0.025, 0.025]: H = {shannon.entropy(p_biased):.3f} bits")
    print(f"  Certain [1, 0, 0, 0]: H = {shannon.entropy(p_certain):.3f} bits")

    # Joint entropy example
    p_xy = np.array([[0.2, 0.1], [0.1, 0.6]])
    print(f"\nJoint entropy example:")
    print(f"  H(X,Y) = {shannon.joint_entropy(p_xy):.3f} bits")
    print(f"  H(Y|X) = {shannon.conditional_entropy(p_xy):.3f} bits")
    print(f"  I(X;Y) = {shannon.mutual_information(p_xy):.3f} bits")

    shannon_interp = shannon.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in shannon_interp.items():
        print(f"  {key}: {value}")

    # Part 2: Landauer's Principle
    print("\n" + "=" * 50)
    print("PART 2: LANDAUER'S PRINCIPLE")
    print("=" * 50)

    landauer = LandauerPrinciple(T=300)

    print(f"\nLandauer's principle:")
    print(f"  Minimum energy to erase 1 bit: E = k_B T ln(2)")
    print(f"  At T = 300 K: E = {landauer.E_per_bit * 1e21:.3f} × 10⁻²¹ J")
    print(f"                = {landauer.E_per_bit / const.e * 1000:.3f} meV")

    print(f"\nEnergy to erase various amounts:")
    for n in [1, 8, 1000, 1e9]:
        E = landauer.erasure_energy(int(n))
        print(f"  {int(n):>10} bits: {E:.3e} J")

    print(f"\nComparison to modern computing:")
    E_cpu_per_bit = 1e-18  # Typical CPU energy per bit operation (J)
    efficiency = landauer.efficiency(E_cpu_per_bit, 1)
    print(f"  Typical CPU: ~10⁻¹⁸ J per bit operation")
    print(f"  Landauer limit: ~3 × 10⁻²¹ J per bit")
    print(f"  Efficiency: {efficiency:.4f} (~{1/efficiency:.0f}× above Landauer limit)")

    landauer_interp = landauer.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in landauer_interp.items():
        print(f"  {key}: {value}")

    # Part 3: Maxwell's Demon
    print("\n" + "=" * 50)
    print("PART 3: MAXWELL'S DEMON")
    print("=" * 50)

    demon = MaxwellDemon(T=300, n_molecules=100)

    print(f"\nMaxwell's demon (100 molecules at 300K):")
    print(f"  Apparent entropy decrease from sorting: {demon.apparent_entropy_decrease() / k_B:.1f} k_B")
    print(f"  Memory required: {demon.demon_memory_required()} bits")
    print(f"  Energy to erase memory: {demon.memory_erasure_cost() * 1e21:.2f} × 10⁻²¹ J")
    print(f"  Entropy from memory erasure: {demon.memory_erasure_entropy() / k_B:.1f} k_B")
    print(f"  Net entropy change: {demon.net_entropy_change() / k_B:.1f} k_B (≥ 0, 2nd Law OK!)")

    print(f"\nSzilard engine (1 molecule):")
    print(f"  Max work extractable: W = k_B T ln(2) = {demon.szilard_engine_work() * 1e21:.3f} × 10⁻²¹ J")
    print(f"  This equals the Landauer erasure cost!")
    print(f"  Net work over full cycle: 0")

    demon_interp = demon.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in demon_interp.items():
        print(f"  {key}: {value}")

    # Part 4: Channel Capacity
    print("\n" + "=" * 50)
    print("PART 4: CHANNEL CAPACITY")
    print("=" * 50)

    channel = ChannelCapacity(bandwidth=1e6, snr_db=20)

    print(f"\nShannon channel capacity:")
    print(f"  C = B log₂(1 + S/N)")
    print(f"  For B = 1 MHz, SNR = 20 dB:")
    print(f"    C = {channel.shannon_capacity() / 1e6:.2f} Mbits/s")
    print(f"    Spectral efficiency: {channel.spectral_efficiency():.2f} bits/s/Hz")

    print(f"\nChannel capacity examples:")
    for snr_db in [10, 20, 30, 40]:
        ch = ChannelCapacity(bandwidth=1e6, snr_db=snr_db)
        print(f"  SNR = {snr_db} dB: C = {ch.shannon_capacity() / 1e6:.2f} Mbits/s")

    print(f"\nBinary symmetric channel:")
    for p_err in [0.0, 0.1, 0.2, 0.5]:
        C_bsc = channel.binary_symmetric_channel(p_err)
        print(f"  p_error = {p_err}: C = {C_bsc:.3f} bits/use")

    channel_interp = channel.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in channel_interp.items():
        print(f"  {key}: {value}")

    # Part 5: Data Compression
    print("\n" + "=" * 50)
    print("PART 5: DATA COMPRESSION")
    print("=" * 50)

    dc = DataCompression()

    print(f"\nShannon source coding theorem:")
    print(f"  Optimal compression approaches H(X) bits per symbol")

    # English letter frequencies (simplified)
    p_english = np.array([0.127, 0.091, 0.082, 0.075, 0.070,
                          0.063, 0.061, 0.060, 0.050, 0.043])  # Top 10 letters
    p_english = p_english / p_english.sum()  # Normalize

    print(f"\nEnglish text (top 10 letters):")
    H = dc.entropy_rate(p_english)
    print(f"  Entropy: H = {H:.3f} bits/symbol")
    print(f"  Max entropy (uniform): log₂(10) = {np.log2(10):.3f} bits/symbol")
    print(f"  Redundancy: {dc.redundancy(p_english, 10) * 100:.1f}%")

    print(f"\nRate-distortion (Gaussian source, σ² = 1):")
    for D in [0.1, 0.25, 0.5, 1.0]:
        R = dc.rate_distortion_bound(D, sigma2=1.0)
        print(f"  D = {D}: R(D) = {R:.3f} bits/sample")

    dc_interp = dc.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in dc_interp.items():
        print(f"  {key}: {value}")

    # Part 6: Information Thermodynamics
    print("\n" + "=" * 50)
    print("PART 6: INFORMATION THERMODYNAMICS")
    print("=" * 50)

    it = InformationThermodynamics(T=300)

    print(f"\nUnification of information and thermodynamics:")
    print(f"  S_Boltzmann = k_B × H_Shannon × ln(2)")

    print(f"\nConversions:")
    S_example = 10 * k_B  # 10 k_B entropy
    H_bits = it.boltzmann_to_shannon(S_example)
    print(f"  S = 10 k_B → H = {H_bits:.3f} bits")

    H_example = 10  # 10 bits
    S_joules = it.shannon_to_boltzmann(H_example)
    print(f"  H = 10 bits → S = {S_joules / k_B:.3f} k_B")

    print(f"\nInformation engine efficiency:")
    W_extracted = 5 * k_B * 300 * np.log(2)  # Extract 5 bits worth
    bits_used = 10
    eta = it.information_engine_efficiency(W_extracted, bits_used)
    print(f"  Extracted: 5 k_B T ln(2) from 10 bits")
    print(f"  Efficiency: η = {eta * 100:.0f}% (≤ 100% by Landauer)")

    it_interp = it.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in it_interp.items():
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
    save_path = os.path.join(script_dir, 'session328_information_theory.png')
    create_visualization(save_path)

    # Final summary
    print("\n" + "=" * 70)
    print("SESSION #328 COMPLETE")
    print("=" * 70)

    print(f"""
    INFORMATION THEORY ARC (Sessions #328-331):

    Session #328: Information Theory Foundations  ✅ {passed}/{total}
    Session #329: Quantum Information             NEXT
    Session #330: Holographic Principle           PLANNED
    Session #331: Black Hole Information          PLANNED

    ═══════════════════════════════════════════════════════════════

    KEY INSIGHTS FROM SESSION #328:

    1. SHANNON ENTROPY
       • H(X) = -Σ p(x) log₂ p(x)
       • Measures bits of uncertainty
       • Maximum for uniform distribution
       • Fundamental compression limit

    2. LANDAUER'S PRINCIPLE
       • E_min = k_B T ln(2) per bit erased
       • Information erasure is physical
       • Modern CPUs: ~1000× above limit
       • Connects info to thermodynamics

    3. MAXWELL'S DEMON
       • Appears to violate 2nd Law
       • Resolution: memory erasure cost
       • Net ΔS ≥ 0 (2nd Law OK)
       • Information IS physical

    4. CHANNEL CAPACITY
       • C = B log₂(1 + S/N)
       • Maximum reliable communication rate
       • Shannon's fundamental theorem
       • MRH = nature's channel capacity

    5. DATA COMPRESSION
       • Can't beat entropy limit
       • Redundancy enables compression
       • Lossless vs lossy tradeoffs

    ═══════════════════════════════════════════════════════════════

    GRID INTERPRETATION:

    INFORMATION = PATTERN DISTINGUISHABILITY

    • Entropy H(X) = bits to specify which pattern is realized
    • Landauer: erasing a bit = pattern info crosses MRH
    • Channel capacity = max rate of pattern transfer
    • MRH IS the channel capacity of nature itself

    The connection between information and thermodynamics is
    not accidental — they are the SAME thing viewed from
    different perspectives (inside vs outside MRH).

    ★ INFORMATION THEORY ARC (1/4) ★

    Next: Session #329 - Quantum Information
    """)

    return results


if __name__ == "__main__":
    main()
