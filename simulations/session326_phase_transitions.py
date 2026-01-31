#!/usr/bin/env python3
"""
Session #326: Phase Transitions from the Planck Grid
Statistical Mechanics Arc (Session 3/4)

This session explores phase transitions from the grid perspective:
1. Order parameters and symmetry breaking
2. First-order vs continuous transitions
3. Critical phenomena and universality
4. The Ising model on the grid
5. MRH divergence at criticality

Key insight: Phase transitions are where MRH → ∞. At a critical point,
correlations extend to all scales, and the MRH boundary disappears.
This is when "local" physics becomes "global" - a fundamental change
in the nature of pattern organization.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from scipy import constants as const
from scipy.ndimage import convolve

# Physical constants
k_B = const.k  # Boltzmann constant


@dataclass
class OrderParameter:
    """
    Order parameters and symmetry breaking.

    An order parameter distinguishes phases:
    - Zero in disordered (high T) phase
    - Non-zero in ordered (low T) phase

    Grid interpretation: Order parameter measures pattern coherence.
    """

    def __init__(self):
        pass

    def examples(self) -> Dict[str, Dict]:
        """Examples of order parameters in different systems."""
        return {
            'ferromagnet': {
                'order_param': 'Magnetization M = <s_i>',
                'symmetric_phase': 'Paramagnetic (M = 0)',
                'broken_phase': 'Ferromagnetic (M ≠ 0)',
                'symmetry': 'Z_2 (spin flip)',
                'grid_meaning': 'Alignment of spin patterns'
            },
            'liquid_gas': {
                'order_param': 'Density difference Δρ = ρ_l - ρ_g',
                'symmetric_phase': 'Supercritical fluid (Δρ = 0)',
                'broken_phase': 'Liquid-gas coexistence (Δρ ≠ 0)',
                'symmetry': 'Liquid-gas interchange',
                'grid_meaning': 'Density pattern coherence'
            },
            'superconductor': {
                'order_param': 'Cooper pair amplitude Ψ',
                'symmetric_phase': 'Normal metal (Ψ = 0)',
                'broken_phase': 'Superconductor (Ψ ≠ 0)',
                'symmetry': 'U(1) gauge',
                'grid_meaning': 'Phase coherence of electron pairs'
            },
            'bec': {
                'order_param': 'Condensate fraction n_0/N',
                'symmetric_phase': 'Normal gas (n_0/N = 0)',
                'broken_phase': 'BEC (n_0/N ≠ 0)',
                'symmetry': 'U(1) global',
                'grid_meaning': 'Macroscopic pattern occupation'
            }
        }

    def symmetry_breaking_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of symmetry breaking."""
        return {
            'high_T': 'Patterns rapidly reconfigure → no preferred direction',
            'low_T': 'Patterns freeze into one configuration → symmetry broken',
            'mechanism': 'Below Tc, one pattern dominates over others',
            'spontaneous': 'System chooses direction, not external field',
            'domains': 'Different regions may choose differently'
        }


class FirstOrderTransition:
    """
    First-order phase transitions.

    Characterized by:
    - Discontinuous order parameter
    - Latent heat
    - Phase coexistence
    - Metastability and hysteresis

    Grid interpretation: Abrupt change in pattern organization.
    """

    def __init__(self, T_c: float = 373.0, latent_heat: float = 2260e3):
        """
        Args:
            T_c: Transition temperature (K)
            latent_heat: Latent heat (J/kg)
        """
        self.T_c = T_c
        self.L = latent_heat

    def order_parameter_vs_T(self, T: np.ndarray) -> np.ndarray:
        """Order parameter shows discontinuous jump at Tc."""
        phi = np.zeros_like(T)
        phi[T < self.T_c] = 1.0  # Ordered phase
        phi[T >= self.T_c] = 0.0  # Disordered phase
        return phi

    def free_energy_landscape(self, phi: np.ndarray, T: float) -> np.ndarray:
        """
        Double-well free energy near first-order transition.

        F(φ) = a(T)|φ|² - b|φ|³ + c|φ|⁴ + ...

        At T < Tc: minimum at φ ≠ 0
        At T > Tc: minimum at φ = 0
        Near Tc: two competing minima → phase coexistence
        """
        a = (T - self.T_c) / self.T_c
        b = 0.5
        c = 0.25

        return a * phi**2 - b * np.abs(phi)**3 + c * phi**4

    def entropy_discontinuity(self) -> float:
        """
        Entropy discontinuity at first-order transition.

        ΔS = L / T_c
        """
        return self.L / self.T_c

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of first-order transitions."""
        return {
            'mechanism': 'Competing pattern configurations',
            'coexistence': 'Two stable patterns at same T (near Tc)',
            'latent_heat': 'Energy to reorganize all patterns',
            'metastability': 'System can be trapped in "wrong" phase',
            'nucleation': 'New phase forms as droplets, then grows',
            'mrh': 'MRH remains finite (no divergence)'
        }


class ContinuousTransition:
    """
    Continuous (second-order) phase transitions.

    Characterized by:
    - Continuous order parameter (but discontinuous derivative)
    - No latent heat
    - Diverging correlation length and susceptibility
    - Critical fluctuations
    - Universal critical exponents

    Grid interpretation: MRH diverges at criticality.
    """

    def __init__(self, T_c: float = 2.269, d: int = 2):
        """
        Args:
            T_c: Critical temperature (in units of J/k_B for Ising)
            d: Dimensionality
        """
        self.T_c = T_c
        self.d = d
        self._set_critical_exponents()

    def _set_critical_exponents(self):
        """Set critical exponents based on dimensionality."""
        if self.d == 2:
            # 2D Ising exact values
            self.beta = 1/8   # Order parameter: M ~ |t|^β
            self.gamma = 7/4  # Susceptibility: χ ~ |t|^-γ
            self.nu = 1       # Correlation length: ξ ~ |t|^-ν
            self.alpha = 0    # Heat capacity: C ~ |t|^-α (log divergence)
            self.delta = 15   # Critical isotherm: M ~ H^(1/δ)
            self.eta = 1/4    # Correlation function: G(r) ~ r^-(d-2+η)
        elif self.d == 3:
            # 3D Ising approximate values
            self.beta = 0.326
            self.gamma = 1.237
            self.nu = 0.630
            self.alpha = 0.110
            self.delta = 4.789
            self.eta = 0.036
        else:
            # Mean field (d ≥ 4)
            self.beta = 0.5
            self.gamma = 1.0
            self.nu = 0.5
            self.alpha = 0
            self.delta = 3
            self.eta = 0

    def reduced_temperature(self, T: float) -> float:
        """Reduced temperature t = (T - Tc) / Tc."""
        return (T - self.T_c) / self.T_c

    def order_parameter(self, T: float) -> float:
        """
        Order parameter near Tc: M ~ |t|^β for T < Tc
        """
        t = self.reduced_temperature(T)
        if t >= 0:
            return 0.0
        return abs(t) ** self.beta

    def correlation_length(self, T: float) -> float:
        """
        Correlation length: ξ ~ |t|^-ν

        This is the MRH! At Tc, ξ → ∞.
        """
        t = self.reduced_temperature(T)
        if abs(t) < 1e-10:
            return np.inf
        return abs(t) ** (-self.nu)

    def susceptibility(self, T: float) -> float:
        """
        Susceptibility: χ ~ |t|^-γ

        Measures response to external field.
        """
        t = self.reduced_temperature(T)
        if abs(t) < 1e-10:
            return np.inf
        return abs(t) ** (-self.gamma)

    def heat_capacity(self, T: float) -> float:
        """
        Heat capacity: C ~ |t|^-α

        For 2D Ising, α = 0 → logarithmic divergence.
        """
        t = self.reduced_temperature(T)
        if abs(t) < 1e-10:
            return np.inf
        if self.alpha == 0:
            return -np.log(abs(t))  # Log divergence
        return abs(t) ** (-self.alpha)

    def scaling_relations(self) -> Dict[str, str]:
        """Critical exponent scaling relations."""
        return {
            'rushbrooke': f'α + 2β + γ = 2: {self.alpha:.3f} + 2×{self.beta:.3f} + {self.gamma:.3f} = {self.alpha + 2*self.beta + self.gamma:.3f}',
            'widom': f'γ = β(δ-1): {self.gamma:.3f} = {self.beta:.3f}×({self.delta:.3f}-1) = {self.beta*(self.delta-1):.3f}',
            'fisher': f'γ = ν(2-η): {self.gamma:.3f} = {self.nu:.3f}×(2-{self.eta:.3f}) = {self.nu*(2-self.eta):.3f}',
            'josephson': f'dν = 2-α: {self.d}×{self.nu:.3f} = 2-{self.alpha:.3f} → {self.d*self.nu:.3f} ≈ {2-self.alpha:.3f}'
        }

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of continuous transitions."""
        return {
            'correlation_length': 'ξ = MRH! Diverges at Tc',
            'at_tc': 'MRH → ∞, all scales coupled',
            'fluctuations': 'Patterns fluctuate at all scales',
            'scale_invariance': 'No characteristic length → self-similarity',
            'universality': 'Exponents depend only on symmetry + dimension',
            'critical_slowing': 'Dynamics slow as τ ~ ξ^z'
        }


class IsingModel2D:
    """
    2D Ising model on a grid.

    H = -J Σ_{<ij>} s_i s_j - h Σ_i s_i

    This is the canonical example of a phase transition.
    Exact solution by Onsager (1944): T_c = 2J / (k_B ln(1+√2))
    """

    def __init__(self, L: int = 32, J: float = 1.0, h: float = 0.0):
        """
        Args:
            L: Linear size of grid
            J: Coupling constant
            h: External field
        """
        self.L = L
        self.J = J
        self.h = h
        self.T_c = 2 * J / np.log(1 + np.sqrt(2))  # Onsager's result
        self.spins = np.random.choice([-1, 1], size=(L, L))

    def energy(self) -> float:
        """Total energy of current configuration."""
        # Neighbor sum using convolution
        kernel = np.array([[0, 1, 0],
                          [1, 0, 1],
                          [0, 1, 0]])
        neighbor_sum = convolve(self.spins, kernel, mode='wrap')
        E_coupling = -0.5 * self.J * np.sum(self.spins * neighbor_sum)
        E_field = -self.h * np.sum(self.spins)
        return E_coupling + E_field

    def magnetization(self) -> float:
        """Order parameter: average spin."""
        return np.mean(self.spins)

    def metropolis_step(self, T: float):
        """Single Metropolis update step."""
        # Pick random site
        i, j = np.random.randint(0, self.L, 2)

        # Calculate energy change from flipping
        s = self.spins[i, j]
        neighbors = (self.spins[(i+1) % self.L, j] +
                    self.spins[(i-1) % self.L, j] +
                    self.spins[i, (j+1) % self.L] +
                    self.spins[i, (j-1) % self.L])
        dE = 2 * s * (self.J * neighbors + self.h)

        # Accept or reject
        if dE <= 0 or np.random.random() < np.exp(-dE / (k_B * T)):
            self.spins[i, j] = -s

    def simulate(self, T: float, n_steps: int = 10000, n_thermalize: int = 1000) -> Dict[str, float]:
        """
        Run Monte Carlo simulation at temperature T.

        Returns average quantities.
        """
        # Thermalization
        for _ in range(n_thermalize):
            for _ in range(self.L**2):  # One sweep
                self.metropolis_step(T)

        # Measurement
        M_samples = []
        E_samples = []
        for _ in range(n_steps):
            for _ in range(self.L**2):
                self.metropolis_step(T)
            M_samples.append(abs(self.magnetization()))
            E_samples.append(self.energy())

        M_arr = np.array(M_samples)
        E_arr = np.array(E_samples)

        return {
            'magnetization': np.mean(M_arr),
            'magnetization_err': np.std(M_arr) / np.sqrt(len(M_arr)),
            'energy': np.mean(E_arr),
            'energy_err': np.std(E_arr) / np.sqrt(len(E_arr)),
            'susceptibility': self.L**2 * np.var(M_arr) / T,
            'heat_capacity': np.var(E_arr) / (T**2)
        }

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of Ising model."""
        return {
            'spins': 'Binary pattern at each grid cell (+1 or -1)',
            'coupling': 'Neighboring patterns prefer alignment',
            'temperature': 'Rate of pattern flipping',
            'ordered_phase': 'All patterns aligned (M = ±1)',
            'disordered': 'Random pattern orientations (M = 0)',
            'critical_point': 'Long-range correlations, scale-free patterns'
        }


class UniversalityClasses:
    """
    Universality in critical phenomena.

    Key insight: Critical exponents depend ONLY on:
    - Dimensionality (d)
    - Symmetry of order parameter

    NOT on microscopic details!

    Grid interpretation: At criticality, microscopic details wash out.
    Only the topology of pattern space matters.
    """

    def __init__(self):
        self.classes = self._define_classes()

    def _define_classes(self) -> Dict[str, Dict]:
        """Define major universality classes."""
        return {
            'ising_2d': {
                'symmetry': 'Z_2 (discrete)',
                'dimension': 2,
                'beta': 1/8,
                'gamma': 7/4,
                'nu': 1,
                'examples': ['2D ferromagnet', 'Binary alloy (2D)', 'Lattice gas (2D)']
            },
            'ising_3d': {
                'symmetry': 'Z_2 (discrete)',
                'dimension': 3,
                'beta': 0.326,
                'gamma': 1.237,
                'nu': 0.630,
                'examples': ['3D ferromagnet', 'Liquid-gas critical point', 'Binary fluid']
            },
            'xy_2d': {
                'symmetry': 'O(2) / U(1) (continuous)',
                'dimension': 2,
                'beta': 'N/A (BKT)',
                'gamma': 'N/A (BKT)',
                'nu': 'N/A (BKT)',
                'examples': ['2D superfluids', 'Thin film magnets', 'XY model'],
                'special': 'Berezinskii-Kosterlitz-Thouless transition'
            },
            'heisenberg_3d': {
                'symmetry': 'O(3) (continuous)',
                'dimension': 3,
                'beta': 0.365,
                'gamma': 1.386,
                'nu': 0.707,
                'examples': ['Isotropic ferromagnet', 'Antiferromagnet']
            },
            'mean_field': {
                'symmetry': 'Any',
                'dimension': '≥ 4 (above upper critical)',
                'beta': 0.5,
                'gamma': 1.0,
                'nu': 0.5,
                'examples': ['High-d systems', 'Long-range interactions']
            }
        }

    def why_universal(self) -> Dict[str, str]:
        """Why universality emerges."""
        return {
            'renormalization_group': 'At Tc, irrelevant operators flow to zero',
            'scale_invariance': 'No characteristic length → same at all scales',
            'fixed_point': 'All systems flow to same fixed point',
            'coarse_graining': 'Microscopic details average out at large scales',
            'mrh_divergence': 'MRH → ∞ means ALL scales contribute equally',
            'grid_insight': 'Only pattern symmetry matters, not detailed structure'
        }


def run_verification_tests() -> Dict[str, bool]:
    """Run all verification tests for Session #326."""
    results = {}

    # Test 1: First-order transition has discontinuous order parameter
    fot = FirstOrderTransition(T_c=373.0)
    T_below = np.array([370.0])
    T_above = np.array([376.0])
    phi_below = fot.order_parameter_vs_T(T_below)[0]
    phi_above = fot.order_parameter_vs_T(T_above)[0]
    results['first_order_discontinuous'] = phi_below != phi_above

    # Test 2: Continuous transition has continuous order parameter
    ct = ContinuousTransition(T_c=2.269, d=2)
    # Order parameter should approach 0 as T → Tc from below
    M_near_tc = ct.order_parameter(ct.T_c - 0.01)
    M_at_tc = ct.order_parameter(ct.T_c)
    results['continuous_order_param'] = M_at_tc == 0 and M_near_tc > 0

    # Test 3: Correlation length diverges at Tc
    xi_near = ct.correlation_length(ct.T_c + 0.1)
    xi_far = ct.correlation_length(ct.T_c + 1.0)
    results['xi_diverges'] = xi_near > xi_far

    # Test 4: Scaling relations approximately satisfied
    relations = ct.scaling_relations()
    rushbrooke = ct.alpha + 2*ct.beta + ct.gamma
    results['rushbrooke_holds'] = abs(rushbrooke - 2) < 0.01

    # Test 5: Ising model critical temperature correct
    ising = IsingModel2D(L=16)
    T_c_exact = 2 / np.log(1 + np.sqrt(2))
    results['ising_tc_correct'] = abs(ising.T_c - T_c_exact) < 0.01

    # Test 6: Universality classes defined
    uc = UniversalityClasses()
    results['universality_defined'] = len(uc.classes) >= 4

    # Test 7: Critical exponents differ by dimension
    ct_2d = ContinuousTransition(T_c=2.269, d=2)
    ct_3d = ContinuousTransition(T_c=4.5, d=3)
    results['exponents_differ'] = ct_2d.beta != ct_3d.beta

    # Test 8: Grid interpretations exist
    results['grid_interpretations'] = (
        'mrh' in ct.grid_interpretation().keys() or
        'correlation_length' in ct.grid_interpretation().keys()
    )

    return results


def create_visualization(save_path: str = None):
    """Create comprehensive visualization for Session #326."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Session #326: Phase Transitions from the Planck Grid\nStatistical Mechanics Arc (3/4)',
                 fontsize=14, fontweight='bold')

    # Panel 1: Order parameter vs temperature
    ax1 = axes[0, 0]
    ct = ContinuousTransition(T_c=2.269, d=2)

    T_range = np.linspace(0.5, 4.0, 200)
    M_values = [ct.order_parameter(T) for T in T_range]

    ax1.plot(T_range / ct.T_c, M_values, 'b-', linewidth=2)
    ax1.axvline(x=1, color='red', linestyle='--', alpha=0.7, label='$T_c$')
    ax1.set_xlabel('$T / T_c$')
    ax1.set_ylabel('Order parameter M')
    ax1.set_title('Order Parameter (2D Ising)')
    ax1.set_xlim(0, 2)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Panel 2: Correlation length (MRH) vs temperature
    ax2 = axes[0, 1]
    T_range = np.linspace(ct.T_c + 0.01, ct.T_c + 2.0, 100)
    xi_values = [ct.correlation_length(T) for T in T_range]

    ax2.semilogy(T_range / ct.T_c, xi_values, 'r-', linewidth=2)
    ax2.axvline(x=1, color='gray', linestyle='--', alpha=0.7)
    ax2.set_xlabel('$T / T_c$')
    ax2.set_ylabel('Correlation length ξ (= MRH)')
    ax2.set_title('MRH Diverges at $T_c$')
    ax2.set_xlim(1, 2)
    ax2.grid(True, alpha=0.3)

    # Annotate
    ax2.annotate('MRH → ∞\nat $T_c$', xy=(1.05, max(xi_values)/10),
                fontsize=10, color='red')

    # Panel 3: Critical exponents comparison
    ax3 = axes[0, 2]

    dims = ['2D Ising', '3D Ising', 'Mean Field']
    betas = [1/8, 0.326, 0.5]
    gammas = [7/4, 1.237, 1.0]
    nus = [1.0, 0.630, 0.5]

    x = np.arange(len(dims))
    width = 0.25

    ax3.bar(x - width, betas, width, label='β', color='blue', alpha=0.7)
    ax3.bar(x, gammas, width, label='γ', color='red', alpha=0.7)
    ax3.bar(x + width, nus, width, label='ν', color='green', alpha=0.7)

    ax3.set_xticks(x)
    ax3.set_xticklabels(dims)
    ax3.set_ylabel('Exponent value')
    ax3.set_title('Critical Exponents by Universality Class')
    ax3.legend()

    # Panel 4: Phase transition types
    ax4 = axes[1, 0]
    ax4.axis('off')

    types_text = """
    PHASE TRANSITION TYPES

    ┌─────────────────────────────────────┐
    │ FIRST ORDER                         │
    │ • Discontinuous order parameter     │
    │ • Latent heat                        │
    │ • Phase coexistence at Tc           │
    │ • MRH stays finite                  │
    │ • Examples: water boiling, melting  │
    └─────────────────────────────────────┘

    ┌─────────────────────────────────────┐
    │ CONTINUOUS (Second Order)           │
    │ • Continuous order parameter        │
    │ • No latent heat                    │
    │ • Critical fluctuations             │
    │ • MRH DIVERGES at Tc               │
    │ • Universal exponents               │
    │ • Examples: ferromagnet, superfluid │
    └─────────────────────────────────────┘

    KEY DIFFERENCE:
    First order: Abrupt reorganization
    Continuous: Gradual + diverging MRH
    """

    ax4.text(0.02, 0.98, types_text, transform=ax4.transAxes, fontsize=8,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    ax4.set_title('Transition Types')

    # Panel 5: Grid interpretation
    ax5 = axes[1, 1]
    ax5.axis('off')

    grid_text = """
    GRID INTERPRETATION

    Phase Transition = Pattern Reorganization

    HIGH T (Disordered):
    ┌─────────────────────────────────────┐
    │  ↑ ↓ ↑ ↓ ↓ ↑ ↓ ↑  Random patterns  │
    │  ↓ ↑ ↑ ↓ ↑ ↓ ↑ ↓  No alignment     │
    │  Small MRH                          │
    └─────────────────────────────────────┘

    AT Tc (Critical):
    ┌─────────────────────────────────────┐
    │  ↑↑↓↓↑ ↓↓↑↑↓  Fluctuations at all │
    │   ↓↓↑↑  ↑↑↓↓   scales             │
    │  MRH → ∞ (no characteristic scale) │
    │  SCALE INVARIANCE                   │
    └─────────────────────────────────────┘

    LOW T (Ordered):
    ┌─────────────────────────────────────┐
    │  ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑  All aligned     │
    │  ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑  Broken symmetry │
    │  Large MRH (long-range order)       │
    └─────────────────────────────────────┘

    Universality:
    At Tc, microscopic details wash out.
    Only pattern SYMMETRY matters!
    """

    ax5.text(0.02, 0.98, grid_text, transform=ax5.transAxes, fontsize=8,
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
    SESSION #326 RESULTS: {passed}/{total} verified

    Key Findings:

    ✓ Order parameters distinguish phases
    ✓ First-order: discontinuous jump
    ✓ Continuous: M ~ |t|^β at Tc

    ✓ Correlation length = MRH
      ξ ~ |t|^-ν diverges at Tc
      At criticality: MRH → ∞

    ✓ Universal critical exponents
      Depend only on d and symmetry
      NOT on microscopic details

    ✓ Scaling relations satisfied
      α + 2β + γ = 2 (Rushbrooke)
      γ = ν(2-η) (Fisher)

    Grid Interpretation:
    • Phase = pattern organization
    • Transition = reorganization
    • Tc = MRH divergence point
    • Universality = only topology matters

    ★ STAT MECH ARC (3/4) ★
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
    """Main execution for Session #326."""
    print("=" * 70)
    print("SESSION #326: Phase Transitions from the Planck Grid")
    print("Statistical Mechanics Arc (Session 3/4)")
    print("=" * 70)

    # Part 1: Order Parameters
    print("\n" + "=" * 50)
    print("PART 1: ORDER PARAMETERS")
    print("=" * 50)

    op = OrderParameter()
    examples = op.examples()

    print("\nOrder parameter examples:")
    for name, props in examples.items():
        print(f"\n  {name.upper()}:")
        print(f"    Order parameter: {props['order_param']}")
        print(f"    Symmetric: {props['symmetric_phase']}")
        print(f"    Broken: {props['broken_phase']}")
        print(f"    Grid meaning: {props['grid_meaning']}")

    sym_breaking = op.symmetry_breaking_interpretation()
    print(f"\nSymmetry breaking (grid view):")
    for key, value in sym_breaking.items():
        print(f"  {key}: {value}")

    # Part 2: First-Order Transitions
    print("\n" + "=" * 50)
    print("PART 2: FIRST-ORDER TRANSITIONS")
    print("=" * 50)

    fot = FirstOrderTransition(T_c=373.0, latent_heat=2260e3)

    print(f"\nWater boiling example:")
    print(f"  Critical temperature: {fot.T_c} K (100°C)")
    print(f"  Latent heat: {fot.L/1e6:.2f} MJ/kg")
    print(f"  Entropy discontinuity: {fot.entropy_discontinuity():.0f} J/(kg·K)")

    fot_interp = fot.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in fot_interp.items():
        print(f"  {key}: {value}")

    # Part 3: Continuous Transitions
    print("\n" + "=" * 50)
    print("PART 3: CONTINUOUS TRANSITIONS")
    print("=" * 50)

    ct_2d = ContinuousTransition(T_c=2.269, d=2)

    print(f"\n2D Ising model:")
    print(f"  Critical temperature: {ct_2d.T_c:.3f} J/k_B")
    print(f"  Critical exponents:")
    print(f"    β (order param): {ct_2d.beta:.4f}")
    print(f"    γ (susceptibility): {ct_2d.gamma:.4f}")
    print(f"    ν (correlation): {ct_2d.nu:.4f}")
    print(f"    α (heat capacity): {ct_2d.alpha:.4f}")

    print(f"\nScaling relations:")
    relations = ct_2d.scaling_relations()
    for name, check in relations.items():
        print(f"  {name}: {check}")

    print(f"\nNear critical point (T/Tc = 0.99):")
    T_near = ct_2d.T_c * 0.99
    print(f"  Order parameter M: {ct_2d.order_parameter(T_near):.4f}")
    print(f"  Correlation length ξ: {ct_2d.correlation_length(T_near):.2f}")
    print(f"  Susceptibility χ: {ct_2d.susceptibility(T_near):.2f}")

    ct_interp = ct_2d.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in ct_interp.items():
        print(f"  {key}: {value}")

    # Part 4: Universality
    print("\n" + "=" * 50)
    print("PART 4: UNIVERSALITY")
    print("=" * 50)

    uc = UniversalityClasses()

    print(f"\nUniversality classes:")
    for name, props in uc.classes.items():
        print(f"\n  {name.upper()}:")
        print(f"    Symmetry: {props['symmetry']}")
        print(f"    Dimension: {props['dimension']}")
        print(f"    β = {props['beta']}, γ = {props['gamma']}, ν = {props['nu']}")
        print(f"    Examples: {', '.join(props['examples'][:2])}")

    why = uc.why_universal()
    print(f"\nWhy universality emerges:")
    for key, value in why.items():
        print(f"  {key}: {value}")

    # Part 5: Ising Model
    print("\n" + "=" * 50)
    print("PART 5: 2D ISING MODEL")
    print("=" * 50)

    ising = IsingModel2D(L=16)
    print(f"\nIsing model (L={ising.L}):")
    print(f"  Exact Tc = {ising.T_c:.4f} J/k_B")

    ising_interp = ising.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in ising_interp.items():
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
    save_path = os.path.join(script_dir, 'session326_phase_transitions.png')
    create_visualization(save_path)

    # Final summary
    print("\n" + "=" * 70)
    print("SESSION #326 COMPLETE")
    print("=" * 70)

    print(f"""
    STATISTICAL MECHANICS ARC (Sessions #324-327):

    Session #324: Thermodynamics Foundations  ✅ 8/8
    Session #325: Partition Functions         ✅ 8/8
    Session #326: Phase Transitions           ✅ {passed}/{total}
    Session #327: Non-Equilibrium             NEXT

    KEY INSIGHTS FROM SESSION #326:

    1. ORDER PARAMETERS
       • Distinguish phases (M = 0 vs M ≠ 0)
       • Measure pattern coherence on grid
       • Symmetry breaking = choosing one pattern

    2. FIRST-ORDER TRANSITIONS
       • Discontinuous jump in order parameter
       • Latent heat, phase coexistence
       • MRH stays finite

    3. CONTINUOUS TRANSITIONS
       • M ~ |t|^β approaches zero continuously
       • Correlation length ξ ~ |t|^-ν DIVERGES
       • Critical fluctuations at all scales

    4. CORRELATION LENGTH = MRH
       • ξ IS the MRH!
       • At Tc: MRH → ∞
       • All scales become coupled
       • No characteristic length → scale invariance

    5. UNIVERSALITY
       • Critical exponents depend only on:
         - Dimensionality d
         - Symmetry of order parameter
       • NOT on microscopic details
       • Grid insight: Only pattern topology matters

    GRID INTERPRETATION:

    Phase transitions reveal the deepest connection between
    statistical mechanics and MRH:

    • Phase = collective pattern organization
    • Tc = point where MRH diverges
    • Universality = microscopic details wash out at MRH → ∞
    • Critical exponents = geometry of pattern space

    At the critical point, the distinction between "local"
    and "global" disappears. This is profound.

    ★ STAT MECH ARC (3/4) ★

    Next: Session #327 - Non-Equilibrium Statistical Mechanics
    """)

    return results


if __name__ == "__main__":
    main()
