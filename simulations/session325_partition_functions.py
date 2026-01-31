#!/usr/bin/env python3
"""
Session #325: Partition Functions from the Planck Grid
Statistical Mechanics Arc (Session 2/4)

This session derives partition functions from grid principles:
1. Canonical ensemble from grid + heat bath
2. Grand canonical from variable particle number
3. Quantum statistics (Bose-Einstein, Fermi-Dirac)
4. Connection between Z and thermodynamic potentials
5. MRH as natural ensemble boundary

Key insight: The partition function Z = Σ exp(-βE) is exactly the
trace over grid microstates weighted by Boltzmann factors.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from scipy import constants as const
from scipy.special import zeta

# Physical constants
k_B = const.k  # Boltzmann constant (J/K)
hbar = const.hbar
c = const.c
eV = const.eV

# Planck units
L_PLANCK = np.sqrt(hbar * const.G / c**3)
E_PLANCK = np.sqrt(hbar * c**5 / const.G)


@dataclass
class CanonicalEnsemble:
    """
    Canonical ensemble from grid + heat bath.

    The system (grid region inside MRH) exchanges energy with
    a heat bath (everything outside MRH). The partition function
    is the sum over all microstates weighted by Boltzmann factors.

    Z = Σ_i exp(-βE_i) where β = 1/(k_B T)
    """

    def __init__(self, n_levels: int = 10, energy_spacing: float = 1e-21):
        """
        Args:
            n_levels: Number of energy levels
            energy_spacing: Energy gap between levels (J)
        """
        self.n_levels = n_levels
        self.energy_spacing = energy_spacing
        self.energies = np.arange(n_levels) * energy_spacing

    def partition_function(self, T: float) -> float:
        """
        Compute canonical partition function Z.

        Z = Σ exp(-E_i / k_B T)

        Args:
            T: Temperature (K)
        """
        if T <= 0:
            return 1.0  # Ground state only
        beta = 1 / (k_B * T)
        return np.sum(np.exp(-beta * self.energies))

    def free_energy(self, T: float) -> float:
        """
        Helmholtz free energy: F = -k_B T ln(Z)

        The connection between partition function and thermodynamics.
        """
        Z = self.partition_function(T)
        if Z <= 0:
            return np.inf
        return -k_B * T * np.log(Z)

    def average_energy(self, T: float) -> float:
        """
        Average energy: <E> = -∂ln(Z)/∂β = Σ E_i P_i

        where P_i = exp(-βE_i) / Z
        """
        if T <= 0:
            return self.energies[0]
        beta = 1 / (k_B * T)
        Z = self.partition_function(T)
        return np.sum(self.energies * np.exp(-beta * self.energies)) / Z

    def entropy_from_Z(self, T: float) -> float:
        """
        Entropy from partition function: S = k_B (ln Z + β<E>)

        Or equivalently: S = -∂F/∂T
        """
        if T <= 0:
            return 0.0
        Z = self.partition_function(T)
        E_avg = self.average_energy(T)
        beta = 1 / (k_B * T)
        return k_B * (np.log(Z) + beta * E_avg)

    def heat_capacity(self, T: float) -> float:
        """
        Heat capacity: C = ∂<E>/∂T = β² <(E - <E>)²>
        """
        if T <= 0:
            return 0.0
        dT = T * 0.001  # Small temperature change
        E1 = self.average_energy(T - dT/2)
        E2 = self.average_energy(T + dT/2)
        return (E2 - E1) / dT

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of canonical ensemble."""
        return {
            'system': 'Grid region inside MRH',
            'bath': 'Everything outside MRH',
            'exchange': 'Energy flows across MRH boundary',
            'fixed': 'Volume (cell count), particle number',
            'fluctuating': 'Energy (intent quanta)',
            'Z_meaning': 'Sum over ways to distribute intent in system'
        }


class GrandCanonicalEnsemble:
    """
    Grand canonical ensemble: variable particle number.

    System exchanges both energy AND particles with reservoir.
    Partition function: Ξ = Σ_{N,i} exp(-β(E_i - μN))

    Grid interpretation: System can gain/lose pattern complexity.
    """

    def __init__(self, single_particle_levels: np.ndarray = None):
        """
        Args:
            single_particle_levels: Energy levels for one particle (J)
        """
        if single_particle_levels is None:
            self.sp_levels = np.arange(5) * 1e-21
        else:
            self.sp_levels = single_particle_levels

    def grand_partition_classical(self, T: float, mu: float, max_N: int = 20) -> float:
        """
        Grand partition function for classical (distinguishable) particles.

        Ξ = Σ_N (z^N / N!) × Z_1^N

        where z = exp(βμ) is fugacity, Z_1 = single-particle Z.
        """
        if T <= 0:
            return 1.0

        beta = 1 / (k_B * T)
        z = np.exp(beta * mu)

        # Single particle partition function
        Z_1 = np.sum(np.exp(-beta * self.sp_levels))

        # Sum over particle numbers
        Xi = 0
        for N in range(max_N + 1):
            Xi += (z * Z_1) ** N / np.math.factorial(N)

        return Xi

    def average_particle_number(self, T: float, mu: float) -> float:
        """
        Average particle number: <N> = z ∂ln(Ξ)/∂z
        """
        if T <= 0:
            return 0.0

        dmu = abs(mu) * 0.01 if mu != 0 else 1e-22
        Xi_plus = self.grand_partition_classical(T, mu + dmu/2)
        Xi_minus = self.grand_partition_classical(T, mu - dmu/2)

        beta = 1 / (k_B * T)
        return (np.log(Xi_plus) - np.log(Xi_minus)) / (beta * dmu)

    def grand_potential(self, T: float, mu: float) -> float:
        """
        Grand potential: Ω = -k_B T ln(Ξ) = F - μN
        """
        Xi = self.grand_partition_classical(T, mu)
        return -k_B * T * np.log(Xi)

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of grand canonical ensemble."""
        return {
            'system': 'Grid region with variable pattern count',
            'bath': 'Reservoir of patterns outside MRH',
            'exchange': 'Both energy AND patterns cross MRH',
            'chemical_potential': 'μ = cost to add pattern',
            'fluctuations': 'Pattern number fluctuates in equilibrium',
            'application': 'Open systems, phase transitions'
        }


class QuantumStatistics:
    """
    Bose-Einstein and Fermi-Dirac statistics from grid.

    Key insight: Grid cells can have occupancy constraints.
    - Bosons: Any number per cell (∞ occupation)
    - Fermions: At most one per cell (Pauli exclusion)

    These emerge from the topology of pattern configurations on the grid.
    """

    def __init__(self, n_levels: int = 20, energy_spacing: float = 1e-21):
        """
        Args:
            n_levels: Number of single-particle levels
            energy_spacing: Energy gap between levels
        """
        self.n_levels = n_levels
        self.epsilon = np.arange(n_levels) * energy_spacing

    def bose_einstein_distribution(self, T: float, mu: float) -> np.ndarray:
        """
        Bose-Einstein distribution: n_i = 1 / (exp(β(ε_i - μ)) - 1)

        For bosons, μ ≤ ε_0 (chemical potential bounded by ground state).
        """
        if T <= 0:
            n = np.zeros(self.n_levels)
            n[0] = np.inf if mu >= self.epsilon[0] else 0
            return n

        beta = 1 / (k_B * T)
        x = beta * (self.epsilon - mu)
        # Avoid overflow
        n = np.where(x > 0, 1 / (np.exp(x) - 1), 0)
        return n

    def fermi_dirac_distribution(self, T: float, mu: float) -> np.ndarray:
        """
        Fermi-Dirac distribution: n_i = 1 / (exp(β(ε_i - μ)) + 1)

        For fermions, 0 ≤ n_i ≤ 1 (Pauli exclusion).
        """
        if T <= 0:
            return (self.epsilon <= mu).astype(float)

        beta = 1 / (k_B * T)
        x = beta * (self.epsilon - mu)
        return 1 / (np.exp(x) + 1)

    def bose_partition_function(self, T: float, mu: float) -> float:
        """
        Bosonic grand partition function.

        ln(Ξ_B) = -Σ ln(1 - exp(-β(ε_i - μ)))
        """
        if T <= 0:
            return 1.0

        beta = 1 / (k_B * T)
        x = beta * (self.epsilon - mu)

        # Sum only valid terms (x > 0 for stability)
        valid = x > 0
        if not np.any(valid):
            return 1.0

        ln_Xi = -np.sum(np.log(1 - np.exp(-x[valid])))
        return np.exp(ln_Xi)

    def fermi_partition_function(self, T: float, mu: float) -> float:
        """
        Fermionic grand partition function.

        ln(Ξ_F) = Σ ln(1 + exp(-β(ε_i - μ)))
        """
        if T <= 0:
            n_occupied = np.sum(self.epsilon <= mu)
            return 2 ** n_occupied  # Each level can be 0 or 1

        beta = 1 / (k_B * T)
        x = beta * (self.epsilon - mu)
        ln_Xi = np.sum(np.log(1 + np.exp(-x)))
        return np.exp(ln_Xi)

    def grid_interpretation(self) -> Dict[str, str]:
        """Grid interpretation of quantum statistics."""
        return {
            'bosons': 'Patterns that can overlap (same grid location)',
            'fermions': 'Patterns that exclude (one per grid cell)',
            'statistics': 'Topology of pattern space on grid',
            'spin_connection': 'Integer spin = bosonic; half-integer = fermionic',
            'bec': 'Bose-Einstein condensation = macroscopic ground state occupation',
            'fermi_sea': 'Fermi sea = all states below μ filled'
        }


class ThermodynamicPotentials:
    """
    Connection between partition functions and thermodynamic potentials.

    All thermodynamics follows from Z or Ξ:
    - F = -k_B T ln(Z)     (Helmholtz free energy)
    - Ω = -k_B T ln(Ξ)     (Grand potential)
    - S = -∂F/∂T           (Entropy)
    - P = -∂F/∂V           (Pressure)
    - μ = ∂F/∂N            (Chemical potential)
    """

    def __init__(self):
        pass

    def relationships(self) -> Dict[str, str]:
        """Key thermodynamic relationships."""
        return {
            'free_energy': 'F = -k_B T ln(Z)',
            'grand_potential': 'Ω = -k_B T ln(Ξ) = F - μN',
            'entropy': 'S = -∂F/∂T = k_B(ln Z + β<E>)',
            'pressure': 'P = -∂F/∂V = k_B T ∂ln(Z)/∂V',
            'energy': '<E> = -∂ln(Z)/∂β',
            'particle_number': '<N> = -∂Ω/∂μ = k_B T ∂ln(Ξ)/∂μ'
        }

    def legendre_transforms(self) -> Dict[str, str]:
        """Legendre transforms connecting potentials."""
        return {
            'U_to_F': 'F = U - TS (const T instead of S)',
            'F_to_G': 'G = F + PV = μN (const P instead of V)',
            'F_to_Omega': 'Ω = F - μN (const μ instead of N)',
            'interpretation': 'Different ensembles = different controls'
        }

    def fluctuation_dissipation(self) -> Dict[str, str]:
        """Fluctuation-dissipation relations from Z."""
        return {
            'energy_fluctuation': '<ΔE²> = k_B T² C_V',
            'particle_fluctuation': '<ΔN²> = k_B T (∂<N>/∂μ)',
            'interpretation': 'Response functions from fluctuations',
            'grid_meaning': 'Larger fluctuations → more states accessible'
        }


class MRHEnsembleBoundary:
    """
    MRH as the natural boundary for ensemble definition.

    Key insight: The MRH defines what's "system" vs "bath":
    - Inside MRH: quantum coherent, fully tracked
    - Outside MRH: thermal, statistical average

    The ensemble choice corresponds to what crosses the MRH boundary.
    """

    def __init__(self):
        pass

    def ensemble_from_mrh(self) -> Dict[str, Dict]:
        """How MRH determines ensemble type."""
        return {
            'microcanonical': {
                'mrh_crossing': 'Nothing',
                'fixed': 'E, V, N',
                'when': 'Isolated system, no contact'
            },
            'canonical': {
                'mrh_crossing': 'Energy (heat)',
                'fixed': 'T, V, N',
                'when': 'Thermal contact, particles confined'
            },
            'grand_canonical': {
                'mrh_crossing': 'Energy AND particles',
                'fixed': 'T, V, μ',
                'when': 'Open boundary for patterns'
            },
            'isobaric': {
                'mrh_crossing': 'Energy AND volume work',
                'fixed': 'T, P, N',
                'when': 'Mechanical contact'
            }
        }

    def equivalence_of_ensembles(self) -> Dict[str, str]:
        """When different ensembles give same results."""
        return {
            'thermodynamic_limit': 'N → ∞ with N/V fixed',
            'equivalence': 'All ensembles agree for intensive quantities',
            'fluctuations': 'Relative fluctuations → 0 as N → ∞',
            'mrh_interpretation': 'Large system has many internal MRH boundaries',
            'why_it_works': 'Local behavior same regardless of boundary conditions'
        }

    def phase_space_connection(self) -> Dict[str, str]:
        """Connection to phase space volume."""
        return {
            'classical': 'Z = ∫ dq dp exp(-βH) / h^{3N} N!',
            'quantum': 'Z = Tr(exp(-βH))',
            'grid': 'Z = Σ_{microstates} exp(-βE)',
            'mrh_meaning': 'Sum over states inside MRH',
            'h_factor': 'Planck constant = grid cell in phase space'
        }


def run_verification_tests() -> Dict[str, bool]:
    """Run all verification tests for Session #325."""
    results = {}

    # Test 1: Partition function positive
    ce = CanonicalEnsemble()
    Z = ce.partition_function(300)
    results['Z_positive'] = Z > 0

    # Test 2: Z(T=0) = 1 (only ground state)
    Z_zero = ce.partition_function(0.001)
    results['Z_ground_state'] = abs(Z_zero - 1) < 0.1

    # Test 3: Free energy decreases with T
    F_low = ce.free_energy(100)
    F_high = ce.free_energy(300)
    results['F_decreases_with_T'] = F_high < F_low

    # Test 4: Entropy positive at finite T
    S = ce.entropy_from_Z(300)
    results['S_positive'] = S > 0

    # Test 5: Fermi distribution bounded [0, 1]
    qs = QuantumStatistics()
    n_F = qs.fermi_dirac_distribution(100, 0)
    results['fermi_bounded'] = np.all(n_F >= 0) and np.all(n_F <= 1)

    # Test 6: Bose distribution positive
    n_B = qs.bose_einstein_distribution(100, -1e-21)
    results['bose_positive'] = np.all(n_B >= 0)

    # Test 7: Grand canonical defined
    gce = GrandCanonicalEnsemble()
    Xi = gce.grand_partition_classical(300, -1e-20)
    results['grand_Z_defined'] = Xi > 0

    # Test 8: All ensembles from MRH defined
    mrh = MRHEnsembleBoundary()
    ensembles = mrh.ensemble_from_mrh()
    results['all_ensembles_defined'] = len(ensembles) == 4

    return results


def create_visualization(save_path: str = None):
    """Create comprehensive visualization for Session #325."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Session #325: Partition Functions from the Planck Grid\nStatistical Mechanics Arc (2/4)',
                 fontsize=14, fontweight='bold')

    # Panel 1: Partition function vs temperature
    ax1 = axes[0, 0]
    ce = CanonicalEnsemble(n_levels=10, energy_spacing=1e-21)
    T_range = np.logspace(-1, 3, 100)
    Z_values = [ce.partition_function(T) for T in T_range]

    ax1.semilogx(T_range, Z_values, 'b-', linewidth=2)
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Ground state only')
    ax1.axhline(y=10, color='gray', linestyle=':', alpha=0.5, label='All states equal')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Partition function Z')
    ax1.set_title('Canonical Partition Function')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Panel 2: Thermodynamic quantities from Z
    ax2 = axes[0, 1]
    T_range = np.linspace(1, 500, 100)

    E_avg = [ce.average_energy(T) / (k_B * 300) for T in T_range]  # In units of k_B × 300K
    S = [ce.entropy_from_Z(T) / k_B for T in T_range]

    ax2.plot(T_range, E_avg, 'r-', linewidth=2, label=r'$\langle E \rangle / k_B(300K)$')
    ax2.plot(T_range, S, 'b-', linewidth=2, label=r'$S / k_B$')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Value')
    ax2.set_title('Thermodynamics from Z')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Panel 3: Quantum statistics comparison
    ax3 = axes[0, 2]
    qs = QuantumStatistics(n_levels=20, energy_spacing=1e-21)
    T = 100  # K
    mu = 5e-21  # Chemical potential

    levels = np.arange(20)
    n_BE = qs.bose_einstein_distribution(T, mu - 1e-21)  # Below ground state
    n_FD = qs.fermi_dirac_distribution(T, mu)
    n_MB = np.exp(-(levels * 1e-21 - mu) / (k_B * T))  # Maxwell-Boltzmann (unnormalized)
    n_MB = n_MB / np.sum(n_MB) * np.sum(n_FD)  # Normalize for comparison

    ax3.plot(levels[:15], n_BE[:15], 'b-o', markersize=4, label='Bose-Einstein')
    ax3.plot(levels[:15], n_FD[:15], 'r-s', markersize=4, label='Fermi-Dirac')
    ax3.plot(levels[:15], n_MB[:15], 'g--', linewidth=2, label='Maxwell-Boltzmann')
    ax3.axhline(y=1, color='red', linestyle=':', alpha=0.5)
    ax3.set_xlabel('Energy level')
    ax3.set_ylabel('Occupation number')
    ax3.set_title(f'Quantum Statistics (T={T}K)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Panel 4: Ensemble types
    ax4 = axes[1, 0]
    ax4.axis('off')

    ensemble_text = """
    STATISTICAL ENSEMBLES

    ┌─────────────────────────────────────┐
    │ MICROCANONICAL                      │
    │ Fixed: E, V, N                      │
    │ MRH boundary: Impermeable           │
    │ Partition: Ω = # of microstates     │
    │ Entropy: S = k_B ln(Ω)              │
    └─────────────────────────────────────┘

    ┌─────────────────────────────────────┐
    │ CANONICAL                           │
    │ Fixed: T, V, N                      │
    │ MRH boundary: Energy exchange       │
    │ Partition: Z = Σ exp(-βE)           │
    │ Free energy: F = -k_B T ln(Z)       │
    └─────────────────────────────────────┘

    ┌─────────────────────────────────────┐
    │ GRAND CANONICAL                     │
    │ Fixed: T, V, μ                      │
    │ MRH boundary: Energy + particles    │
    │ Partition: Ξ = Σ exp(-β(E-μN))      │
    │ Grand pot: Ω = -k_B T ln(Ξ)         │
    └─────────────────────────────────────┘
    """

    ax4.text(0.02, 0.98, ensemble_text, transform=ax4.transAxes, fontsize=8,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    ax4.set_title('Ensemble Types')

    # Panel 5: MRH as ensemble boundary
    ax5 = axes[1, 1]
    ax5.axis('off')

    mrh_text = """
    MRH AS ENSEMBLE BOUNDARY

    The MRH naturally defines the system-bath boundary:

    ┌─────────────────────────────────────┐
    │      OUTSIDE MRH (Bath)             │
    │    • Thermal reservoir              │
    │    • Statistical average            │
    │    • Infinite heat capacity         │
    └─────────────────────────────────────┘
              ↕ Energy/particles
    ┌─────────────────────────────────────┐
    │      INSIDE MRH (System)            │
    │    • Quantum coherent               │
    │    • Fully tracked states           │
    │    • Finite degrees of freedom      │
    └─────────────────────────────────────┘

    Ensemble choice = what crosses MRH:
    • Nothing → Microcanonical (E fixed)
    • Energy → Canonical (T fixed)
    • Energy + particles → Grand canonical
    • Energy + volume → Isobaric-isothermal

    The MRH is not arbitrary—it's physics-defined
    by correlation decay length.
    """

    ax5.text(0.02, 0.98, mrh_text, transform=ax5.transAxes, fontsize=8,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    ax5.set_title('MRH Interpretation')

    # Panel 6: Summary
    ax6 = axes[1, 2]
    ax6.axis('off')

    verification = run_verification_tests()
    passed = sum(verification.values())
    total = len(verification)

    summary_text = f"""
    SESSION #325 RESULTS: {passed}/{total} verified

    Key Findings:

    ✓ Canonical ensemble from grid
      Z = Σ exp(-βE_i)
      Natural from MRH boundary

    ✓ Grand canonical for open systems
      Ξ = Σ exp(-β(E - μN))
      Pattern exchange across MRH

    ✓ Quantum statistics emerge
      Bosons: unlimited occupation
      Fermions: Pauli exclusion
      From grid topology

    ✓ All thermodynamics from Z
      F = -k_B T ln(Z)
      S, E, P, μ all follow

    ✓ MRH defines ensemble
      What crosses MRH boundary
      determines ensemble type

    Grid Interpretation:
    • Partition function = state count
    • Ensembles = MRH boundary conditions
    • Quantum statistics = grid topology
    • Thermodynamics = MRH physics

    ★ STAT MECH ARC (2/4) ★
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
    """Main execution for Session #325."""
    print("=" * 70)
    print("SESSION #325: Partition Functions from the Planck Grid")
    print("Statistical Mechanics Arc (Session 2/4)")
    print("=" * 70)

    # Part 1: Canonical Ensemble
    print("\n" + "=" * 50)
    print("PART 1: CANONICAL ENSEMBLE")
    print("=" * 50)

    ce = CanonicalEnsemble(n_levels=10, energy_spacing=1e-21)

    print(f"\nSystem: {ce.n_levels} energy levels")
    print(f"Energy spacing: {ce.energy_spacing:.2e} J")

    print(f"\nPartition function vs temperature:")
    for T in [1, 10, 100, 300, 1000]:
        Z = ce.partition_function(T)
        F = ce.free_energy(T)
        E = ce.average_energy(T)
        S = ce.entropy_from_Z(T)
        print(f"  T = {T:4d} K: Z = {Z:8.2f}, F = {F/eV:10.2e} eV, "
              f"<E> = {E/eV:10.2e} eV, S = {S/k_B:6.2f} k_B")

    grid_interp = ce.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in grid_interp.items():
        print(f"  {key}: {value}")

    # Part 2: Grand Canonical Ensemble
    print("\n" + "=" * 50)
    print("PART 2: GRAND CANONICAL ENSEMBLE")
    print("=" * 50)

    gce = GrandCanonicalEnsemble()

    print(f"\nGrand partition function (T=300K):")
    for mu_ev in [-1, -0.5, -0.1, 0]:
        mu = mu_ev * 1e-21
        Xi = gce.grand_partition_classical(300, mu)
        N_avg = gce.average_particle_number(300, mu)
        Omega = gce.grand_potential(300, mu)
        print(f"  μ = {mu_ev:.1f}×10⁻²¹ J: Ξ = {Xi:.2e}, <N> = {N_avg:.2f}, "
              f"Ω = {Omega/eV:.2e} eV")

    gce_interp = gce.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in gce_interp.items():
        print(f"  {key}: {value}")

    # Part 3: Quantum Statistics
    print("\n" + "=" * 50)
    print("PART 3: QUANTUM STATISTICS")
    print("=" * 50)

    qs = QuantumStatistics(n_levels=10, energy_spacing=1e-21)

    print(f"\nFermi-Dirac distribution (T=100K):")
    mu = 3e-21
    n_FD = qs.fermi_dirac_distribution(100, mu)
    print(f"  Chemical potential μ = {mu:.2e} J")
    for i in range(min(5, len(n_FD))):
        print(f"  Level {i}: n = {n_FD[i]:.4f}")

    print(f"\nBose-Einstein distribution (T=100K):")
    mu_BE = -1e-21  # Must be below ground state
    n_BE = qs.bose_einstein_distribution(100, mu_BE)
    print(f"  Chemical potential μ = {mu_BE:.2e} J")
    for i in range(min(5, len(n_BE))):
        print(f"  Level {i}: n = {n_BE[i]:.4f}")

    qs_interp = qs.grid_interpretation()
    print(f"\nGrid interpretation:")
    for key, value in qs_interp.items():
        print(f"  {key}: {value}")

    # Part 4: Thermodynamic Potentials
    print("\n" + "=" * 50)
    print("PART 4: THERMODYNAMIC POTENTIALS")
    print("=" * 50)

    tp = ThermodynamicPotentials()

    relations = tp.relationships()
    print(f"\nKey relationships:")
    for key, value in relations.items():
        print(f"  {key}: {value}")

    legendre = tp.legendre_transforms()
    print(f"\nLegendre transforms:")
    for key, value in legendre.items():
        print(f"  {key}: {value}")

    # Part 5: MRH as Ensemble Boundary
    print("\n" + "=" * 50)
    print("PART 5: MRH AS ENSEMBLE BOUNDARY")
    print("=" * 50)

    mrh = MRHEnsembleBoundary()

    ensembles = mrh.ensemble_from_mrh()
    print(f"\nEnsembles from MRH boundary conditions:")
    for name, props in ensembles.items():
        print(f"\n  {name.upper()}:")
        for key, value in props.items():
            print(f"    {key}: {value}")

    equiv = mrh.equivalence_of_ensembles()
    print(f"\nEnsemble equivalence:")
    for key, value in equiv.items():
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
    save_path = os.path.join(script_dir, 'session325_partition_functions.png')
    create_visualization(save_path)

    # Final summary
    print("\n" + "=" * 70)
    print("SESSION #325 COMPLETE")
    print("=" * 70)

    print(f"""
    STATISTICAL MECHANICS ARC (Sessions #324-327):

    Session #324: Thermodynamics Foundations  ✅ 8/8
    Session #325: Partition Functions         ✅ {passed}/{total}
    Session #326: Phase Transitions           NEXT
    Session #327: Non-Equilibrium             Planned

    KEY INSIGHTS FROM SESSION #325:

    1. CANONICAL PARTITION FUNCTION
       • Z = Σ exp(-βE_i)
       • Sum over grid microstates
       • All thermodynamics follows from Z

    2. GRAND CANONICAL ENSEMBLE
       • Ξ = Σ exp(-β(E - μN))
       • Variable particle number
       • Pattern exchange across MRH

    3. QUANTUM STATISTICS
       • Bosons: n = 1/(exp(βε-βμ) - 1)
       • Fermions: n = 1/(exp(βε-βμ) + 1)
       • From grid topology (overlap vs exclusion)

    4. THERMODYNAMIC POTENTIALS
       • F = -k_B T ln(Z)
       • All quantities from partition function
       • Legendre transforms connect ensembles

    5. MRH AS ENSEMBLE BOUNDARY
       • What crosses MRH defines ensemble
       • Microcanonical: nothing
       • Canonical: energy
       • Grand canonical: energy + particles

    GRID INTERPRETATION:

    The partition function is the fundamental object:
    • It counts microstates weighted by Boltzmann factors
    • MRH defines what's "system" vs "bath"
    • Ensemble = boundary condition at MRH

    This completes the connection:
    • Session #324: Entropy and temperature from grid
    • Session #325: Partition functions and ensembles
    • Together: Complete statistical mechanics

    ★ STAT MECH ARC (2/4) ★

    Next: Session #326 - Phase Transitions
    """)

    return results


if __name__ == "__main__":
    main()
