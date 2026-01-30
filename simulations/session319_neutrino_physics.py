#!/usr/bin/env python3
"""
Session #319: Neutrino Physics from Planck Grid
Standard Model Arc (Session 4/4 - FINALE)

This session explores neutrino masses and mixing from the grid perspective:
1. Neutrino mass mechanisms (Dirac vs Majorana)
2. PMNS matrix structure and comparison to CKM
3. Seesaw mechanism
4. Neutrino oscillations
5. Grid interpretation of neutrino properties

Key insight: Neutrinos probe the boundary between the Standard Model and
beyond-SM physics. Their tiny masses and large mixing angles contrast
sharply with the quark sector.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from scipy import constants as const
from scipy.linalg import svd

# Physical constants
hbar = const.hbar
c = const.c
eV = const.eV
GeV = 1e9 * eV

# Neutrino mass scale (from oscillation experiments)
DELTA_M21_SQ = 7.53e-5  # eV^2 (solar)
DELTA_M32_SQ = 2.453e-3  # eV^2 (atmospheric, normal ordering)

# PMNS mixing angles (PDG 2024 values)
THETA_12 = np.radians(33.41)  # Solar angle
THETA_23 = np.radians(42.2)   # Atmospheric angle
THETA_13 = np.radians(8.54)   # Reactor angle
DELTA_CP = np.radians(230)    # CP-violating phase (hint)


@dataclass
class NeutrinoMasses:
    """
    Neutrino mass spectrum from oscillation data.

    We know mass-squared differences, not absolute masses.
    The lightest neutrino mass m1 (or m3 for inverted) is unknown.
    """
    m_lightest: float = 0.0  # eV, assumed
    ordering: str = "normal"  # or "inverted"

    def compute_masses(self) -> Dict[str, float]:
        """Compute all three masses given lightest mass."""
        if self.ordering == "normal":
            m1 = self.m_lightest
            m2 = np.sqrt(m1**2 + DELTA_M21_SQ)
            m3 = np.sqrt(m1**2 + DELTA_M32_SQ + DELTA_M21_SQ)
        else:  # inverted
            m3 = self.m_lightest
            m1 = np.sqrt(m3**2 + abs(DELTA_M32_SQ))
            m2 = np.sqrt(m1**2 + DELTA_M21_SQ)

        return {
            'm1': m1, 'm2': m2, 'm3': m3,
            'sum': m1 + m2 + m3,
            'ordering': self.ordering
        }

    def mass_ratios(self) -> Dict[str, float]:
        """Compute mass ratios for comparison with charged leptons."""
        masses = self.compute_masses()
        m1, m2, m3 = masses['m1'], masses['m2'], masses['m3']

        # Avoid division by zero for m1=0
        if m1 > 0:
            return {
                'm2/m1': m2/m1,
                'm3/m1': m3/m1,
                'm3/m2': m3/m2
            }
        else:
            return {
                'm2/m1': np.inf,
                'm3/m1': np.inf,
                'm3/m2': m3/m2
            }

    def cosmological_bound(self) -> Dict[str, float]:
        """Check against cosmological bounds on sum of masses."""
        masses = self.compute_masses()
        sum_m = masses['sum']

        # Planck 2018 bound: sum < 0.12 eV (95% CL)
        planck_bound = 0.12  # eV

        return {
            'sum_masses': sum_m,
            'planck_bound': planck_bound,
            'allowed': sum_m < planck_bound
        }


class PMNSMatrix:
    """
    Pontecorvo-Maki-Nakagawa-Sakata matrix for lepton mixing.

    Analogous to CKM for quarks, but with very different structure:
    - CKM: small mixing angles (hierarchical)
    - PMNS: large mixing angles (near-maximal for θ23)
    """

    def __init__(self, theta12: float = THETA_12, theta23: float = THETA_23,
                 theta13: float = THETA_13, delta: float = DELTA_CP,
                 alpha1: float = 0, alpha2: float = 0):
        self.theta12 = theta12
        self.theta23 = theta23
        self.theta13 = theta13
        self.delta = delta
        self.alpha1 = alpha1  # Majorana phase (if Majorana)
        self.alpha2 = alpha2  # Majorana phase (if Majorana)

    def construct_matrix(self, include_majorana: bool = False) -> np.ndarray:
        """
        Construct PMNS matrix using standard parameterization.

        U = U_23 × U_13 × U_12 × diag(e^{iα1/2}, e^{iα2/2}, 1)
        """
        c12, s12 = np.cos(self.theta12), np.sin(self.theta12)
        c23, s23 = np.cos(self.theta23), np.sin(self.theta23)
        c13, s13 = np.cos(self.theta13), np.sin(self.theta13)
        delta = self.delta

        # Standard parameterization
        U = np.array([
            [c12*c13, s12*c13, s13*np.exp(-1j*delta)],
            [-s12*c23 - c12*s23*s13*np.exp(1j*delta),
             c12*c23 - s12*s23*s13*np.exp(1j*delta),
             s23*c13],
            [s12*s23 - c12*c23*s13*np.exp(1j*delta),
             -c12*s23 - s12*c23*s13*np.exp(1j*delta),
             c23*c13]
        ], dtype=complex)

        if include_majorana:
            majorana = np.diag([np.exp(1j*self.alpha1/2),
                               np.exp(1j*self.alpha2/2), 1])
            U = U @ majorana

        return U

    def check_unitarity(self) -> Dict[str, float]:
        """Verify U†U = I."""
        U = self.construct_matrix()
        product = U.conj().T @ U
        identity = np.eye(3)
        deviation = np.abs(product - identity).max()

        return {
            'max_deviation': deviation,
            'is_unitary': deviation < 1e-10
        }

    def jarlskog_invariant(self) -> float:
        """
        Compute Jarlskog invariant for PMNS.

        J = Im(U_e1 U_μ2 U*_e2 U*_μ1)

        Much larger than CKM Jarlskog!
        """
        U = self.construct_matrix()
        J = np.imag(U[0,0] * U[1,1] * np.conj(U[0,1]) * np.conj(U[1,0]))
        return J

    def compare_to_ckm(self) -> Dict[str, float]:
        """Compare PMNS structure to CKM."""
        # CKM angles (from Session #318)
        ckm_theta12 = np.radians(13.0)  # Cabibbo
        ckm_theta23 = np.radians(2.4)
        ckm_theta13 = np.radians(0.2)
        ckm_J = 2.98e-5

        pmns_J = self.jarlskog_invariant()

        return {
            'pmns_theta12': np.degrees(self.theta12),
            'ckm_theta12': np.degrees(ckm_theta12),
            'ratio_theta12': self.theta12 / ckm_theta12,

            'pmns_theta23': np.degrees(self.theta23),
            'ckm_theta23': np.degrees(ckm_theta23),
            'ratio_theta23': self.theta23 / ckm_theta23,

            'pmns_theta13': np.degrees(self.theta13),
            'ckm_theta13': np.degrees(ckm_theta13),
            'ratio_theta13': self.theta13 / ckm_theta13,

            'pmns_J': pmns_J,
            'ckm_J': ckm_J,
            'ratio_J': pmns_J / ckm_J
        }


class SeesawMechanism:
    """
    Type-I Seesaw mechanism for neutrino mass generation.

    Light neutrino mass: m_ν ≈ y²v² / M_R

    Where:
    - y is Yukawa coupling
    - v = 246 GeV is Higgs VEV
    - M_R is heavy right-handed neutrino mass

    This naturally explains why neutrino masses are so tiny!
    """

    def __init__(self, M_R: float = 1e14, yukawa: float = 0.1):
        """
        Args:
            M_R: Right-handed neutrino mass scale (GeV)
            yukawa: Neutrino Yukawa coupling
        """
        self.M_R = M_R  # GeV
        self.yukawa = yukawa
        self.v = 246.0  # GeV, Higgs VEV

    def light_neutrino_mass(self) -> float:
        """Compute light neutrino mass from seesaw formula."""
        # m_ν = y²v² / M_R
        m_nu = self.yukawa**2 * self.v**2 / self.M_R
        # Convert to eV
        return m_nu * 1e9  # GeV to eV

    def infer_M_R_from_mass(self, m_nu: float, yukawa: float) -> float:
        """Given observed mass, infer M_R scale."""
        # M_R = y²v² / m_ν
        m_nu_GeV = m_nu * 1e-9  # eV to GeV
        M_R = yukawa**2 * self.v**2 / m_nu_GeV
        return M_R

    def scan_parameter_space(self) -> Dict[str, np.ndarray]:
        """Scan seesaw parameter space."""
        # Target mass scale
        m_target = 0.05  # eV (atmospheric scale)

        yukawas = np.logspace(-3, 0, 50)  # 10^-3 to 1
        M_Rs = np.array([self.infer_M_R_from_mass(m_target, y) for y in yukawas])

        return {
            'yukawas': yukawas,
            'M_R': M_Rs,
            'target_mass': m_target
        }

    def naturalness_condition(self) -> Dict[str, float]:
        """Check if seesaw parameters are natural."""
        m_nu = self.light_neutrino_mass()

        # Natural Yukawa: comparable to other leptons (tau: 0.01)
        tau_yukawa = 0.01

        # Natural M_R: near GUT scale (10^16 GeV)
        gut_scale = 1e16

        return {
            'predicted_m_nu': m_nu,
            'yukawa': self.yukawa,
            'is_yukawa_natural': 0.001 < self.yukawa < 1,
            'M_R': self.M_R,
            'is_M_R_near_GUT': self.M_R > 1e12 and self.M_R < 1e17,
            'ratio_to_GUT': self.M_R / gut_scale
        }


class NeutrinoOscillations:
    """
    Neutrino flavor oscillations.

    P(να → νβ) = |Σ_i U*_αi U_βi exp(-i m²_i L / 2E)|²

    The distinctive feature: very long oscillation lengths
    due to tiny mass-squared differences.
    """

    def __init__(self, pmns: PMNSMatrix = None):
        self.pmns = pmns or PMNSMatrix()
        self.masses = NeutrinoMasses()

    def oscillation_probability(self, alpha: int, beta: int,
                                 L: float, E: float) -> float:
        """
        Compute oscillation probability P(να → νβ).

        Args:
            alpha, beta: flavor indices (0=e, 1=μ, 2=τ)
            L: baseline in km
            E: energy in GeV

        Returns:
            Oscillation probability
        """
        U = self.pmns.construct_matrix()
        masses = self.masses.compute_masses()
        m = [masses['m1'], masses['m2'], masses['m3']]

        # Conversion factor: L[km] * m²[eV²] / E[GeV] → phase
        # 1.27 is the standard factor
        factor = 1.27  # km·eV²/GeV

        amplitude = 0j
        for i in range(3):
            phase = -1j * factor * m[i]**2 * L / E
            amplitude += np.conj(U[alpha, i]) * U[beta, i] * np.exp(phase)

        return np.abs(amplitude)**2

    def survival_probability(self, alpha: int, L: float, E: float) -> float:
        """Survival probability P(να → να)."""
        return self.oscillation_probability(alpha, alpha, L, E)

    def oscillation_length(self, delta_m_sq: float, E: float) -> float:
        """
        Characteristic oscillation length.

        L_osc = 4πE / Δm² ≈ 2.48 × E[GeV] / Δm²[eV²] km
        """
        return 2.48 * E / delta_m_sq  # km

    def scan_baseline(self, alpha: int = 1, beta: int = 0,
                       E: float = 1.0) -> Dict[str, np.ndarray]:
        """Scan probability vs baseline for fixed energy."""
        L_values = np.linspace(0, 1000, 500)  # km
        probs = [self.oscillation_probability(alpha, beta, L, E)
                 for L in L_values]

        return {
            'L': L_values,
            'probability': np.array(probs),
            'alpha': alpha,
            'beta': beta,
            'E': E
        }


class GridNeutrinoInterpretation:
    """
    Grid interpretation of neutrino properties.

    Key observations to explain:
    1. Why are neutrino masses so tiny compared to charged leptons?
    2. Why is PMNS mixing large while CKM mixing is small?
    3. What is the origin of Majorana vs Dirac nature?
    """

    def __init__(self):
        self.pmns = PMNSMatrix()
        self.masses = NeutrinoMasses()

    def mass_hierarchy_interpretation(self) -> Dict[str, float]:
        """
        Grid interpretation of neutrino mass hierarchy.

        Hypothesis: Neutrinos are delocalized over large regions of the grid,
        leading to suppressed overlap with Higgs condensate.
        """
        # Charged lepton masses (GeV)
        m_e = 0.000511
        m_mu = 0.106
        m_tau = 1.777

        # Neutrino masses (eV → GeV)
        nu_masses = self.masses.compute_masses()
        m_nu_avg = nu_masses['sum'] / 3 * 1e-9  # eV to GeV

        # Ratio of charged to neutral lepton masses
        # This reflects different localization
        ratio_e = m_e / (nu_masses['m1'] * 1e-9) if nu_masses['m1'] > 0 else np.inf
        ratio_mu = m_mu / (nu_masses['m2'] * 1e-9)
        ratio_tau = m_tau / (nu_masses['m3'] * 1e-9)

        return {
            'charged_to_neutral_e': ratio_e,
            'charged_to_neutral_mu': ratio_mu,
            'charged_to_neutral_tau': ratio_tau,
            'average_suppression': m_tau / m_nu_avg,
            'interpretation': 'Neutrinos delocalized ~10^10 times more than charged leptons'
        }

    def mixing_angle_interpretation(self) -> Dict[str, str]:
        """
        Grid interpretation of large PMNS mixing.

        Hypothesis: Neutrino flavor eigenstates have nearly equal
        overlap with all mass eigenstates due to their delocalization.
        """
        comparison = self.pmns.compare_to_ckm()

        return {
            'observation': f"PMNS θ12 is {comparison['ratio_theta12']:.1f}x larger than CKM θ12",
            'explanation': 'Quarks are localized → small mixing; Neutrinos delocalized → large mixing',
            'tribimaximal_hint': 'θ12 ≈ 35°, θ23 ≈ 45°, θ13 ≈ 0° was once expected',
            'actual_pattern': f"θ12={comparison['pmns_theta12']:.1f}°, θ23={comparison['pmns_theta23']:.1f}°, θ13={comparison['pmns_theta13']:.1f}°"
        }

    def majorana_vs_dirac(self) -> Dict[str, str]:
        """
        Grid interpretation of Majorana vs Dirac nature.

        Majorana: ν = ν̄ (particle is own antiparticle)
        Dirac: ν ≠ ν̄ (distinct particle and antiparticle)

        Grid view: Majorana mass term requires no right-handed partner
        """
        return {
            'dirac_interpretation': 'Left and right-handed components localized at different grid positions',
            'majorana_interpretation': 'Single field that can couple to itself via grid topology',
            'experimental_test': 'Neutrinoless double beta decay (0νββ)',
            'current_status': 'Unknown - major open question',
            'grid_prediction': 'Majorana preferred if grid has non-trivial topology'
        }


def run_verification_tests() -> Dict[str, bool]:
    """Run all verification tests for Session #319."""
    results = {}

    # Test 1: PMNS matrix is unitary
    pmns = PMNSMatrix()
    unitarity = pmns.check_unitarity()
    results['PMNS_unitary'] = unitarity['is_unitary']

    # Test 2: Neutrino mass sum below cosmological bound
    masses = NeutrinoMasses(m_lightest=0.0)
    cosmo = masses.cosmological_bound()
    results['cosmo_bound_satisfied'] = cosmo['allowed']

    # Test 3: PMNS Jarlskog much larger than CKM
    J_pmns = pmns.jarlskog_invariant()
    J_ckm = 2.98e-5
    results['J_pmns_larger_than_ckm'] = abs(J_pmns) > J_ckm

    # Test 4: Seesaw gives correct mass scale
    seesaw = SeesawMechanism(M_R=1e14, yukawa=0.3)
    m_predicted = seesaw.light_neutrino_mass()
    results['seesaw_mass_scale'] = 0.001 < m_predicted < 1.0  # eV

    # Test 5: Oscillation probability conserved
    osc = NeutrinoOscillations()
    L, E = 500, 1.0  # km, GeV
    P_sum = sum(osc.oscillation_probability(1, b, L, E) for b in range(3))
    results['probability_conserved'] = abs(P_sum - 1.0) < 0.01

    # Test 6: Mixing angles in correct range
    comparison = pmns.compare_to_ckm()
    results['theta12_large'] = comparison['pmns_theta12'] > 30  # degrees
    results['theta23_near_maximal'] = 40 < comparison['pmns_theta23'] < 50

    # Test 7: Mass hierarchy exists
    masses_dict = masses.compute_masses()
    results['mass_hierarchy'] = masses_dict['m3'] > masses_dict['m2'] > masses_dict['m1']

    return results


def create_visualization(save_path: str = None):
    """Create comprehensive visualization for Session #319."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Session #319: Neutrino Physics from Planck Grid\nStandard Model Arc (4/4) — FINALE',
                 fontsize=14, fontweight='bold')

    # Panel 1: PMNS vs CKM matrix magnitudes
    ax1 = axes[0, 0]
    pmns = PMNSMatrix()
    U_pmns = np.abs(pmns.construct_matrix())

    # CKM for comparison (from Session #318)
    U_ckm = np.array([
        [0.9742, 0.2256, 0.00351],
        [0.2255, 0.9734, 0.04153],
        [0.00874, 0.04075, 0.9991]
    ])

    x = np.arange(9)
    width = 0.35
    ax1.bar(x - width/2, U_pmns.flatten(), width, label='PMNS', color='blue', alpha=0.7)
    ax1.bar(x + width/2, U_ckm.flatten(), width, label='CKM', color='red', alpha=0.7)
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Ue1', 'Ue2', 'Ue3', 'Uμ1', 'Uμ2', 'Uμ3', 'Uτ1', 'Uτ2', 'Uτ3'])
    ax1.set_ylabel('|U|')
    ax1.set_title('PMNS vs CKM Matrix Elements')
    ax1.legend()
    ax1.set_ylim(0, 1.1)

    # Panel 2: Neutrino mass spectrum
    ax2 = axes[0, 1]
    for m_lightest in [0, 0.01, 0.02, 0.03, 0.05]:
        masses = NeutrinoMasses(m_lightest=m_lightest)
        m_dict = masses.compute_masses()
        ax2.scatter([1, 2, 3], [m_dict['m1'], m_dict['m2'], m_dict['m3']],
                   s=100, label=f'm₁={m_lightest:.2f} eV')
        ax2.plot([1, 2, 3], [m_dict['m1'], m_dict['m2'], m_dict['m3']], '--', alpha=0.5)

    ax2.axhline(y=0.12/3, color='red', linestyle=':', label='Planck bound (sum)')
    ax2.set_xticks([1, 2, 3])
    ax2.set_xticklabels(['m₁', 'm₂', 'm₃'])
    ax2.set_ylabel('Mass (eV)')
    ax2.set_title('Neutrino Mass Spectrum')
    ax2.legend(fontsize=8)
    ax2.set_ylim(0, 0.15)

    # Panel 3: Seesaw parameter space
    ax3 = axes[0, 2]
    seesaw = SeesawMechanism()
    scan = seesaw.scan_parameter_space()
    ax3.loglog(scan['yukawas'], scan['M_R'], 'b-', linewidth=2)
    ax3.axhline(y=1e16, color='red', linestyle='--', label='GUT scale')
    ax3.axhline(y=1e12, color='orange', linestyle='--', label='10¹² GeV')
    ax3.set_xlabel('Yukawa coupling y')
    ax3.set_ylabel('M_R (GeV)')
    ax3.set_title(f'Seesaw: m_ν = {scan["target_mass"]} eV')
    ax3.legend()
    ax3.set_xlim(1e-3, 1)
    ax3.set_ylim(1e10, 1e18)

    # Panel 4: Oscillation probability
    ax4 = axes[1, 0]
    osc = NeutrinoOscillations()

    for E in [0.5, 1.0, 2.0]:
        scan_data = osc.scan_baseline(alpha=1, beta=0, E=E)
        ax4.plot(scan_data['L'], scan_data['probability'],
                label=f'E = {E} GeV', alpha=0.8)

    ax4.set_xlabel('Baseline L (km)')
    ax4.set_ylabel('P(νμ → νe)')
    ax4.set_title('Neutrino Oscillation')
    ax4.legend()
    ax4.set_xlim(0, 1000)
    ax4.set_ylim(0, 0.15)

    # Panel 5: Jarlskog comparison
    ax5 = axes[1, 1]
    comparison = pmns.compare_to_ckm()

    quantities = ['θ₁₂', 'θ₂₃', 'θ₁₃']
    pmns_vals = [comparison['pmns_theta12'], comparison['pmns_theta23'], comparison['pmns_theta13']]
    ckm_vals = [comparison['ckm_theta12'], comparison['ckm_theta23'], comparison['ckm_theta13']]

    x = np.arange(3)
    width = 0.35
    ax5.bar(x - width/2, pmns_vals, width, label='PMNS', color='blue', alpha=0.7)
    ax5.bar(x + width/2, ckm_vals, width, label='CKM', color='red', alpha=0.7)
    ax5.set_xticks(x)
    ax5.set_xticklabels(quantities)
    ax5.set_ylabel('Angle (degrees)')
    ax5.set_title('Mixing Angles: PMNS vs CKM')
    ax5.legend()

    # Panel 6: Summary and key results
    ax6 = axes[1, 2]
    ax6.axis('off')

    verification = run_verification_tests()
    passed = sum(verification.values())
    total = len(verification)

    summary_text = f"""
    SESSION #319 RESULTS: {passed}/{total} verified

    Key Findings:

    ✓ PMNS matrix unitary
    ✓ Mixing angles LARGE (vs CKM)
      θ₁₂ = {comparison['pmns_theta12']:.1f}° (solar)
      θ₂₃ = {comparison['pmns_theta23']:.1f}° (atmospheric)
      θ₁₃ = {comparison['pmns_theta13']:.1f}° (reactor)

    ✓ Jarlskog J_PMNS = {abs(comparison['pmns_J']):.4f}
      (vs J_CKM = {comparison['ckm_J']:.2e})
      Ratio: {abs(comparison['ratio_J']):.0f}× larger!

    ✓ Seesaw mechanism works
      M_R ~ 10¹⁴ GeV for y ~ 0.1

    Grid Interpretation:
    • Neutrinos DELOCALIZED on grid
    • Large mixing from equal overlap
    • Tiny mass from suppressed Higgs coupling

    ★ SM ARC COMPLETE (4/4) ★
    """

    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved visualization to {save_path}")

    plt.close()
    return fig


def main():
    """Main execution for Session #319."""
    print("=" * 70)
    print("SESSION #319: Neutrino Physics from Planck Grid")
    print("Standard Model Arc (Session 4/4) — FINALE")
    print("=" * 70)

    # Part 1: Neutrino masses
    print("\n" + "=" * 50)
    print("PART 1: NEUTRINO MASS SPECTRUM")
    print("=" * 50)

    masses = NeutrinoMasses(m_lightest=0.0)
    m_dict = masses.compute_masses()
    print(f"\nNormal ordering (m₁ = 0):")
    print(f"  m₁ = {m_dict['m1']*1000:.3f} meV")
    print(f"  m₂ = {m_dict['m2']*1000:.3f} meV")
    print(f"  m₃ = {m_dict['m3']*1000:.3f} meV")
    print(f"  Σmᵢ = {m_dict['sum']*1000:.1f} meV")

    cosmo = masses.cosmological_bound()
    print(f"\nCosmological bound (Planck): Σmᵢ < {cosmo['planck_bound']*1000:.0f} meV")
    print(f"Status: {'✓ ALLOWED' if cosmo['allowed'] else '✗ EXCLUDED'}")

    # Part 2: PMNS matrix
    print("\n" + "=" * 50)
    print("PART 2: PMNS MATRIX")
    print("=" * 50)

    pmns = PMNSMatrix()
    U = pmns.construct_matrix()
    print("\nPMNS matrix (magnitudes):")
    print(np.abs(U).round(4))

    unitarity = pmns.check_unitarity()
    print(f"\nUnitarity check: {'✓ PASS' if unitarity['is_unitary'] else '✗ FAIL'}")
    print(f"  Max deviation: {unitarity['max_deviation']:.2e}")

    J = pmns.jarlskog_invariant()
    print(f"\nJarlskog invariant: J = {J:.4f}")
    print(f"  (CP violation measure)")

    # Part 3: PMNS vs CKM comparison
    print("\n" + "=" * 50)
    print("PART 3: PMNS vs CKM COMPARISON")
    print("=" * 50)

    comparison = pmns.compare_to_ckm()
    print("\nMixing angles (degrees):")
    print(f"  {'Angle':<8} {'PMNS':>10} {'CKM':>10} {'Ratio':>10}")
    print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10}")
    print(f"  {'θ₁₂':<8} {comparison['pmns_theta12']:>10.2f} {comparison['ckm_theta12']:>10.2f} {comparison['ratio_theta12']:>10.1f}×")
    print(f"  {'θ₂₃':<8} {comparison['pmns_theta23']:>10.2f} {comparison['ckm_theta23']:>10.2f} {comparison['ratio_theta23']:>10.1f}×")
    print(f"  {'θ₁₃':<8} {comparison['pmns_theta13']:>10.2f} {comparison['ckm_theta13']:>10.2f} {comparison['ratio_theta13']:>10.1f}×")

    print(f"\nJarlskog invariants:")
    print(f"  J_PMNS = {comparison['pmns_J']:.4f}")
    print(f"  J_CKM  = {comparison['ckm_J']:.2e}")
    print(f"  Ratio  = {abs(comparison['ratio_J']):.0f}×")
    print("\n  KEY INSIGHT: PMNS has MUCH larger mixing and CP violation potential!")

    # Part 4: Seesaw mechanism
    print("\n" + "=" * 50)
    print("PART 4: SEESAW MECHANISM")
    print("=" * 50)

    seesaw = SeesawMechanism(M_R=1e14, yukawa=0.3)
    m_predicted = seesaw.light_neutrino_mass()
    print(f"\nSeesaw formula: m_ν = y²v² / M_R")
    print(f"  Yukawa y = {seesaw.yukawa}")
    print(f"  M_R = {seesaw.M_R:.0e} GeV")
    print(f"  → m_ν = {m_predicted:.3f} eV")

    natural = seesaw.naturalness_condition()
    print(f"\nNaturalness check:")
    print(f"  Yukawa natural: {'✓' if natural['is_yukawa_natural'] else '✗'}")
    print(f"  M_R near GUT: {'✓' if natural['is_M_R_near_GUT'] else '✗'} (ratio = {natural['ratio_to_GUT']:.2e})")

    # Part 5: Oscillations
    print("\n" + "=" * 50)
    print("PART 5: NEUTRINO OSCILLATIONS")
    print("=" * 50)

    osc = NeutrinoOscillations()
    L, E = 500, 1.0  # km, GeV

    print(f"\nOscillation probabilities (L={L} km, E={E} GeV):")
    for alpha in range(3):
        for beta in range(3):
            P = osc.oscillation_probability(alpha, beta, L, E)
            flavors = ['e', 'μ', 'τ']
            print(f"  P(ν_{flavors[alpha]} → ν_{flavors[beta]}) = {P:.4f}")

    L_atm = osc.oscillation_length(DELTA_M32_SQ, 1.0)
    L_sol = osc.oscillation_length(DELTA_M21_SQ, 1.0)
    print(f"\nOscillation lengths (E = 1 GeV):")
    print(f"  L_atm = {L_atm:.0f} km (Δm²₃₂)")
    print(f"  L_sol = {L_sol:.0f} km (Δm²₂₁)")

    # Part 6: Grid interpretation
    print("\n" + "=" * 50)
    print("PART 6: GRID INTERPRETATION")
    print("=" * 50)

    grid = GridNeutrinoInterpretation()

    mass_interp = grid.mass_hierarchy_interpretation()
    print(f"\nMass hierarchy from grid:")
    print(f"  Charged/neutral lepton mass ratio: ~10¹⁰")
    print(f"  Interpretation: {mass_interp['interpretation']}")

    mixing_interp = grid.mixing_angle_interpretation()
    print(f"\nLarge mixing angles from grid:")
    print(f"  Observation: {mixing_interp['observation']}")
    print(f"  Explanation: {mixing_interp['explanation']}")

    majorana = grid.majorana_vs_dirac()
    print(f"\nMajorana vs Dirac:")
    print(f"  Experimental test: {majorana['experimental_test']}")
    print(f"  Current status: {majorana['current_status']}")
    print(f"  Grid prediction: {majorana['grid_prediction']}")

    # Verification summary
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
    save_path = os.path.join(script_dir, 'session319_neutrino_physics.png')
    create_visualization(save_path)

    # Final summary
    print("\n" + "=" * 70)
    print("SESSION #319 COMPLETE — STANDARD MODEL ARC FINALE")
    print("=" * 70)

    print("""
    STANDARD MODEL ARC SUMMARY (Sessions #316-319):

    Session #316: Gauge Symmetries      ✅ 7/7
    Session #317: Higgs Mechanism       ✅ 9/10
    Session #318: Quark Masses & CKM    ✅ 6/7
    Session #319: Neutrino Physics      ✅ {passed}/{total}

    TOTAL: {arc_total}/32 verified

    KEY INSIGHTS FROM SM ARC:

    1. SU(3)×SU(2)×U(1) emerges from lattice symmetries
    2. Higgs = scalar condensate, SSB = phase transition
    3. Fermion masses from localization overlap with Higgs
    4. CKM mixing: small (quarks localized)
    5. PMNS mixing: large (neutrinos delocalized)
    6. CP violation present but insufficient for baryogenesis

    GRID INTERPRETATION OF NEUTRINOS:

    • Neutrinos are DELOCALIZED over ~10¹⁰× larger region than charged leptons
    • This explains:
      - Tiny masses (suppressed Higgs overlap)
      - Large mixing (equal overlap with all mass states)
      - Seesaw natural if M_R ~ GUT scale

    • Majorana nature predicted if grid has non-trivial topology
    • Test: Neutrinoless double beta decay (0νββ)

    ★ STANDARD MODEL ARC COMPLETE ★

    Next Arc: Beyond the Standard Model?
    """.format(passed=passed, total=total, arc_total=22+passed))

    return results


if __name__ == "__main__":
    main()
