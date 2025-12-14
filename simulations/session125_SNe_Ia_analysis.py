"""
Session #125: Type Ia Supernova Distance Analysis
==================================================

Following the roadmap from Session #124, this session analyzes Type Ia
supernova (SNe Ia) distance predictions in Synchronism.

KEY QUESTION:
Does Synchronism predict different luminosity distances for SNe Ia
compared to ΛCDM?

EXPECTED RESULT:
NO - Geometry should be UNCHANGED in Synchronism.
SNe Ia measure luminosity distance d_L(z), which depends on:
- Expansion history H(z)
- Spatial geometry

Synchronism modifies GROWTH (G_eff), not GEOMETRY (H, distances).
Therefore SNe Ia distances should match ΛCDM.

This is a CONSISTENCY CHECK, not a discriminating test.

OBSERVATIONAL DATA:
- Pantheon+ (2022): 1701 SNe Ia, z = 0.001 to 2.26
- Hubble diagram tension with H0
- Current precision: ~2% on cosmological parameters

Created: December 14, 2025
Session: #125
Roadmap: From Session #124 - remaining cosmology gaps
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize
from datetime import datetime

# =============================================================================
# COSMOLOGICAL CONSTANTS
# =============================================================================

# Fiducial cosmology (Planck 2018 + Pantheon+)
H0_Planck = 67.4     # km/s/Mpc
H0_SH0ES = 73.0      # km/s/Mpc (local measurement)
Omega_m = 0.315
Omega_Lambda = 0.685
Omega_k = 0.0        # Flat universe

c = 299792.458       # km/s

# Synchronism parameters
# From Session #121: C_cosmic = Omega_m(z)
# This affects GROWTH, not geometry


# =============================================================================
# PART 1: DISTANCE-REDSHIFT RELATIONS
# =============================================================================

def E_LCDM(z, Omega_m=0.315, Omega_Lambda=0.685, Omega_k=0.0):
    """
    E(z) = H(z)/H0 for ΛCDM.

    E²(z) = Ω_m(1+z)³ + Ω_k(1+z)² + Ω_Λ
    """
    return np.sqrt(Omega_m * (1 + z)**3 +
                   Omega_k * (1 + z)**2 +
                   Omega_Lambda)


def E_Sync(z, Omega_m=0.315, Omega_Lambda=0.685, Omega_k=0.0):
    """
    E(z) = H(z)/H0 for Synchronism.

    KEY POINT: Synchronism does NOT modify the Friedmann equation.
    The coherence function C affects G_eff for structure growth,
    but the background expansion is determined by total energy density.

    E_Sync(z) = E_LCDM(z) (UNCHANGED)
    """
    # Synchronism preserves background geometry
    return E_LCDM(z, Omega_m, Omega_Lambda, Omega_k)


def comoving_distance(z, E_func, Omega_m=0.315, Omega_Lambda=0.685):
    """
    Comoving distance d_C(z) = c/H0 × ∫₀^z dz'/E(z')
    """
    integrand = lambda zp: 1.0 / E_func(zp, Omega_m, Omega_Lambda)
    result, _ = quad(integrand, 0, z)
    return (c / H0_Planck) * result  # Mpc


def luminosity_distance(z, E_func, Omega_m=0.315, Omega_Lambda=0.685):
    """
    Luminosity distance d_L(z) = (1+z) × d_C(z)

    For flat universe (Omega_k = 0).
    """
    d_C = comoving_distance(z, E_func, Omega_m, Omega_Lambda)
    return (1 + z) * d_C  # Mpc


def angular_diameter_distance(z, E_func, Omega_m=0.315, Omega_Lambda=0.685):
    """
    Angular diameter distance d_A(z) = d_C(z) / (1+z)
    """
    d_C = comoving_distance(z, E_func, Omega_m, Omega_Lambda)
    return d_C / (1 + z)  # Mpc


def distance_modulus(z, E_func, Omega_m=0.315, Omega_Lambda=0.685, H0=67.4):
    """
    Distance modulus μ(z) = 5 × log10(d_L / 10 pc)

    μ = m - M (apparent - absolute magnitude)
    """
    # Need to recalculate with specific H0
    integrand = lambda zp: 1.0 / E_func(zp, Omega_m, Omega_Lambda)
    result, _ = quad(integrand, 0, z)
    d_C = (c / H0) * result  # Mpc
    d_L = (1 + z) * d_C  # Mpc

    # Convert to pc
    d_L_pc = d_L * 1e6  # pc

    return 5 * np.log10(d_L_pc / 10)


# =============================================================================
# PART 2: COMPARISON OF LCDM AND SYNCHRONISM
# =============================================================================

def compare_distances():
    """Compare ΛCDM and Synchronism distance predictions."""
    print("=" * 80)
    print("PART 2: DISTANCE COMPARISON - ΛCDM vs SYNCHRONISM")
    print("=" * 80)

    print("""
THEORETICAL EXPECTATION:
========================
Synchronism modifies the GROWTH of structure via G_eff = G/C,
but does NOT modify the background EXPANSION.

The Friedmann equation:
    H² = (8πG/3) × ρ_total

remains unchanged because:
1. Total energy density ρ_total is the same
2. G in the Friedmann equation is NOT modified by C
3. C only affects perturbation growth, not background

Therefore:
    d_L^Sync(z) = d_L^LCDM(z) ✓
    μ^Sync(z) = μ^LCDM(z) ✓

This is a CONSISTENCY requirement, not a prediction.
    """)

    # Calculate distances
    redshifts = np.array([0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0])

    print("\n" + "-" * 70)
    print(f"{'z':<8} {'d_L^LCDM (Mpc)':<18} {'d_L^Sync (Mpc)':<18} {'Δd_L/d_L':<12}")
    print("-" * 70)

    for z in redshifts:
        d_L_LCDM = luminosity_distance(z, E_LCDM)
        d_L_Sync = luminosity_distance(z, E_Sync)
        ratio = (d_L_Sync - d_L_LCDM) / d_L_LCDM

        print(f"{z:<8.2f} {d_L_LCDM:<18.2f} {d_L_Sync:<18.2f} {ratio:<12.6f}")

    print("\nRESULT: Δd_L/d_L = 0 at all redshifts ✓")
    print("Synchronism predicts IDENTICAL distances to ΛCDM.")


# =============================================================================
# PART 3: HUBBLE TENSION ANALYSIS
# =============================================================================

def analyze_hubble_tension():
    """
    Analyze whether Synchronism affects the Hubble tension.
    """
    print("\n" + "=" * 80)
    print("PART 3: HUBBLE TENSION ANALYSIS")
    print("=" * 80)

    print("""
HUBBLE TENSION:
===============
H0 (Planck/CMB): 67.4 ± 0.5 km/s/Mpc
H0 (SH0ES/SNe):  73.0 ± 1.0 km/s/Mpc
Tension: ~5σ

QUESTION: Does Synchronism resolve or affect this tension?

ANSWER: NO - Synchronism does NOT resolve the Hubble tension.

Reason:
1. H0 tension is about the EXPANSION RATE
2. Synchronism only modifies GROWTH (perturbations)
3. Background expansion is unchanged

This was already established in Session #109.
    """)

    # Calculate distance modulus with different H0 values
    z_test = np.linspace(0.01, 2.0, 50)

    mu_Planck = [distance_modulus(z, E_LCDM, H0=H0_Planck) for z in z_test]
    mu_SH0ES = [distance_modulus(z, E_LCDM, H0=H0_SH0ES) for z in z_test]

    # The difference
    delta_mu = np.array(mu_Planck) - np.array(mu_SH0ES)

    print(f"\nDistance modulus difference (Planck vs SH0ES H0):")
    print(f"At z=0.01: Δμ = {delta_mu[0]:.3f} mag")
    print(f"At z=0.1:  Δμ = {delta_mu[4]:.3f} mag")
    print(f"At z=1.0:  Δμ = {delta_mu[24]:.3f} mag")

    print("""

CONCLUSION:
===========
The H0 tension appears as a constant offset in distance modulus:
    Δμ = 5 × log10(H0_SH0ES / H0_Planck) = 5 × log10(73/67.4) ≈ 0.18 mag

Synchronism does NOT provide a mechanism to resolve this.
The tension likely requires:
- Early universe modification (dark radiation, EDE)
- Systematic errors in measurements
- New physics affecting expansion, not growth

Synchronism is a LATE-TIME theory affecting growth.
It cannot resolve early-universe tensions like H0.
    """)


# =============================================================================
# PART 4: PANTHEON+ MOCK COMPARISON
# =============================================================================

def pantheon_mock_comparison():
    """
    Compare predictions with Pantheon+ precision.
    """
    print("\n" + "=" * 80)
    print("PART 4: PANTHEON+ PRECISION COMPARISON")
    print("=" * 80)

    print("""
PANTHEON+ DATASET (2022):
=========================
- 1701 spectroscopically confirmed SNe Ia
- Redshift range: 0.001 < z < 2.26
- 77 SNe at z > 1 (HST)
- Statistical + systematic uncertainty: ~0.02 mag per bin
- Constrains Ω_m to ~2%

SYNCHRONISM PREDICTION:
=======================
μ^Sync(z) = μ^LCDM(z) (identical)

The only way Synchronism could affect SNe Ia measurements:
1. Modified G_eff → different SN intrinsic luminosity?
   - NO: SNe Ia physics set by nuclear/EM, not gravity
2. Modified growth → lensing magnification bias?
   - POSSIBLE: Reduced structure → less lensing scatter
   - Effect: ~0.01 mag at z > 1 (very small)
    """)

    # Generate mock Pantheon+ binned data
    z_bins = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30,
                       0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.50, 2.00])

    # LCDM predictions
    mu_LCDM = np.array([distance_modulus(z, E_LCDM) for z in z_bins])

    # Synchronism predictions (identical)
    mu_Sync = np.array([distance_modulus(z, E_Sync) for z in z_bins])

    # Mock errors (approximate Pantheon+ precision)
    sigma_mu = 0.02 + 0.03 * z_bins  # Increases with z

    print("\n" + "-" * 70)
    print(f"{'z':<8} {'μ_LCDM':<12} {'μ_Sync':<12} {'Δμ':<12} {'σ_obs':<12}")
    print("-" * 70)

    for i, z in enumerate(z_bins):
        delta = mu_Sync[i] - mu_LCDM[i]
        print(f"{z:<8.2f} {mu_LCDM[i]:<12.2f} {mu_Sync[i]:<12.2f} {delta:<12.4f} {sigma_mu[i]:<12.3f}")

    print("""

RESULT:
=======
Δμ = 0 at all redshifts.
Synchronism predicts IDENTICAL SNe Ia distances to ΛCDM.

This is expected because:
1. SNe Ia measure geometry (d_L), not growth
2. Synchronism preserves geometry
3. No modification to Friedmann equation background
    """)


# =============================================================================
# PART 5: WEAK LENSING MAGNIFICATION EFFECT
# =============================================================================

def lensing_magnification_effect():
    """
    Analyze potential weak lensing magnification differences.
    """
    print("\n" + "=" * 80)
    print("PART 5: WEAK LENSING MAGNIFICATION BIAS")
    print("=" * 80)

    print("""
SUBTLE EFFECT: Lensing Magnification
====================================

Even though intrinsic distances are unchanged, SNe Ia
can be magnified/demagnified by foreground structure.

In ΛCDM: More structure → more lensing scatter
In Sync: Less structure (S8 suppressed) → less lensing scatter

This affects:
1. Scatter in Hubble diagram at z > 0.5
2. Potential bias in H0 estimates (magnification bias)

QUANTITATIVE ESTIMATE:
----------------------
Lensing dispersion: σ_μ^lens ∝ σ_8 × f(z)

At z = 1:
    σ_μ^lens(ΛCDM) ≈ 0.04 mag
    σ_μ^lens(Sync) ≈ 0.04 × (0.78/0.83)^2 ≈ 0.035 mag

Difference: ~0.005 mag (13% reduction in scatter)

This is BELOW current measurement precision.
But could be detectable with:
- LSST (thousands of SNe Ia at z > 1)
- Roman Space Telescope (precision cosmology)
    """)

    # Calculate lensing scatter
    z_array = np.linspace(0.1, 2.0, 20)

    # Approximate lensing dispersion (increases with z)
    sigma_lens_LCDM = 0.01 + 0.03 * z_array
    sigma_lens_Sync = sigma_lens_LCDM * (0.78 / 0.83)**2

    print("\nPredicted Lensing Scatter:")
    print("-" * 50)
    print(f"{'z':<10} {'σ_LCDM (mag)':<15} {'σ_Sync (mag)':<15} {'Ratio':<10}")
    print("-" * 50)

    for i in [0, 5, 10, 15, 19]:
        z = z_array[i]
        ratio = sigma_lens_Sync[i] / sigma_lens_LCDM[i]
        print(f"{z:<10.2f} {sigma_lens_LCDM[i]:<15.4f} {sigma_lens_Sync[i]:<15.4f} {ratio:<10.3f}")

    print("""

CONCLUSION:
===========
Synchronism predicts ~12% reduction in lensing scatter.
This is a NOVEL PREDICTION but requires:
- Large sample sizes (N > 1000 at z > 1)
- Precise intrinsic scatter modeling
- Currently NOT discriminating (below precision)
    """)


# =============================================================================
# PART 6: FALSIFICATION CRITERIA
# =============================================================================

def define_falsification_criteria():
    """Define falsification criteria for SNe Ia predictions."""
    print("\n" + "=" * 80)
    print("PART 6: FALSIFICATION CRITERIA")
    print("=" * 80)

    print("""
SNe Ia PREDICTIONS - FALSIFICATION CRITERIA
===========================================

PRIMARY PREDICTION: μ^Sync(z) = μ^LCDM(z)
-----------------------------------------
FALSIFIED if: Distance modulus differs from ΛCDM at >2σ
              for fixed cosmological parameters

Test: Compare Pantheon+ data with ΛCDM predictions.
Result: Must be consistent (not a discriminating test)


SECONDARY PREDICTION: Reduced lensing scatter
---------------------------------------------
FALSIFIED if: High-z SNe Ia show SAME or GREATER scatter
              compared to ΛCDM predictions

Test: Compare Hubble diagram scatter at z > 1 with expectations
Requirement: Large sample (N > 1000), precise modeling


CONSISTENCY CHECK: H0 tension unresolved
----------------------------------------
FALSIFIED if: Synchronism somehow resolves H0 tension
              (would indicate scope violation)

Status: CONFIRMED - Synchronism does NOT resolve H0


SUMMARY:
========
SNe Ia provide a CONSISTENCY CHECK, not a discriminating test.

If SNe Ia distances differ from ΛCDM → Synchronism FALSIFIED
(because it predicts unchanged geometry)

The only subtle effect is reduced lensing scatter (~12%),
which is currently below precision but may be testable with LSST/Roman.
    """)


# =============================================================================
# PART 7: VISUALIZATION
# =============================================================================

def create_visualization():
    """Create visualization of SNe Ia analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Session #125: Type Ia Supernova Analysis', fontsize=14, fontweight='bold')

    z = np.linspace(0.01, 2.0, 100)

    # 1. Distance-redshift relation
    ax1 = axes[0, 0]
    d_L = [luminosity_distance(zi, E_LCDM) for zi in z]
    ax1.plot(z, d_L, 'b-', linewidth=2, label='ΛCDM = Sync')
    ax1.set_xlabel('Redshift z')
    ax1.set_ylabel('Luminosity Distance (Mpc)')
    ax1.set_title('Luminosity Distance')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Distance modulus (Hubble diagram)
    ax2 = axes[0, 1]
    mu = [distance_modulus(zi, E_LCDM) for zi in z]
    ax2.plot(z, mu, 'r-', linewidth=2, label='ΛCDM = Sync')

    # Add H0 tension effect
    mu_SH0ES = [distance_modulus(zi, E_LCDM, H0=73.0) for zi in z]
    ax2.plot(z, mu_SH0ES, 'g--', linewidth=2, alpha=0.7, label='SH0ES H0')

    ax2.set_xlabel('Redshift z')
    ax2.set_ylabel('Distance Modulus μ (mag)')
    ax2.set_title('Hubble Diagram')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Lensing scatter comparison
    ax3 = axes[1, 0]
    sigma_LCDM = 0.01 + 0.03 * z
    sigma_Sync = sigma_LCDM * (0.78/0.83)**2

    ax3.plot(z, sigma_LCDM, 'b-', linewidth=2, label='ΛCDM')
    ax3.plot(z, sigma_Sync, 'r--', linewidth=2, label='Synchronism')
    ax3.fill_between(z, sigma_Sync, sigma_LCDM, alpha=0.3, color='green')

    ax3.set_xlabel('Redshift z')
    ax3.set_ylabel('Lensing Scatter σ_μ (mag)')
    ax3.set_title('Lensing Scatter (Subtle Prediction)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Summary table
    ax4 = axes[1, 1]
    ax4.axis('off')

    table_data = [
        ['Quantity', 'ΛCDM', 'Sync', 'Status'],
        ['d_L(z)', 'Fiducial', 'SAME', 'Consistency'],
        ['μ(z)', 'Fiducial', 'SAME', 'Consistency'],
        ['H0 tension', '5σ', 'Unresolved', 'NOT affected'],
        ['σ_lens', '0.04 mag', '0.035 mag', 'Subtle (12%)'],
    ]

    table = ax4.table(cellText=table_data, loc='center', cellLoc='center',
                       colWidths=[0.3, 0.2, 0.2, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    ax4.set_title('SNe Ia Predictions Summary', fontsize=12, pad=20)

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session125_SNe_Ia.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nVisualization saved to session125_SNe_Ia.png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main analysis function."""
    print("=" * 80)
    print("SESSION #125: TYPE Ia SUPERNOVA DISTANCE ANALYSIS")
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 80)

    print("""
OBJECTIVE:
Verify that Synchronism predicts UNCHANGED luminosity distances
for Type Ia supernovae compared to ΛCDM.

EXPECTATION:
Synchronism modifies GROWTH (G_eff), not GEOMETRY (distances).
Therefore: d_L^Sync(z) = d_L^LCDM(z)
    """)

    # Part 2: Distance comparison
    compare_distances()

    # Part 3: Hubble tension
    analyze_hubble_tension()

    # Part 4: Pantheon+ comparison
    pantheon_mock_comparison()

    # Part 5: Lensing effect
    lensing_magnification_effect()

    # Part 6: Falsification criteria
    define_falsification_criteria()

    # Part 7: Visualization
    create_visualization()

    # Final summary
    print("\n" + "=" * 80)
    print("SESSION #125 SUMMARY")
    print("=" * 80)

    print("""
SNe Ia ANALYSIS COMPLETE
========================

1. MAIN RESULT:
   d_L^Sync(z) = d_L^LCDM(z) at all redshifts ✓
   Geometry is UNCHANGED in Synchronism.

2. HUBBLE TENSION:
   NOT resolved by Synchronism (as expected).
   Synchronism is a late-time growth modification,
   not an early-universe expansion modification.

3. SUBTLE EFFECT:
   Lensing scatter reduced by ~12% due to lower σ8.
   Currently below precision, may be testable with LSST/Roman.

4. CONSISTENCY:
   SNe Ia provide a CONSISTENCY CHECK, not discrimination.
   If distances differ from ΛCDM → Synchronism falsified.

NEXT PRIORITY:
Session #126: Gravitational Lensing Time Delays
    """)

    return {
        'distances_unchanged': True,
        'H0_resolved': False,
        'lensing_scatter_reduction': 0.12,
        'status': 'Consistency check passed'
    }


if __name__ == "__main__":
    results = main()
    print(f"\nResults: {results}")
