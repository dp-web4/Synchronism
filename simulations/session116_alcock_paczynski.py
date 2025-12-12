"""
Session #116: Alcock-Paczynski Test Predictions for Synchronism
===============================================================

The Alcock-Paczynski (AP) test is a geometric probe that measures the ratio
of angular to radial distances. This provides a test of expansion history
independent of structure growth.

Key question: Does Synchronism predict any deviation from ΛCDM in the AP test?

Physics:
- AP parameter: F_AP = D_A(z) × H(z) / c
- In ΛCDM: F_AP depends only on H(z) and geometry
- In Synchronism: H(z) matches ΛCDM exactly (by construction)
- Therefore: AP should be UNCHANGED

This session verifies this prediction and explores any subtleties.

Created: December 12, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Cosmological parameters (Planck 2018)
H0 = 70.0  # km/s/Mpc
Omega_m = 0.3
Omega_Lambda = 0.7
c = 299792.458  # km/s


def H_LCDM(z):
    """Hubble parameter H(z) in ΛCDM in km/s/Mpc."""
    return H0 * np.sqrt(Omega_m * (1+z)**3 + Omega_Lambda)


def H_Sync(z):
    """
    Hubble parameter in Synchronism.

    KEY POINT: H(z) is UNCHANGED from ΛCDM.

    This is because:
    1. Modified Friedmann: H² = (8πG/3C) × ρ_m
    2. With C_cosmic = Ω_m(z), this gives exact ΛCDM H(z)
    3. The modification affects GROWTH, not EXPANSION
    """
    return H_LCDM(z)  # Identical by construction


def comoving_distance_LCDM(z):
    """
    Comoving distance D_C(z) in ΛCDM.

    D_C = c × ∫_0^z dz'/H(z')
    """
    def integrand(z_prime):
        return 1 / H_LCDM(z_prime)

    result, _ = integrate.quad(integrand, 0, z)
    return c * result  # Mpc


def comoving_distance_Sync(z):
    """
    Comoving distance in Synchronism.

    Since H(z) is unchanged, D_C is also unchanged.
    """
    def integrand(z_prime):
        return 1 / H_Sync(z_prime)

    result, _ = integrate.quad(integrand, 0, z)
    return c * result  # Mpc


def angular_diameter_distance_LCDM(z):
    """
    Angular diameter distance D_A(z) in ΛCDM.

    D_A = D_C / (1 + z)   [flat universe]
    """
    D_C = comoving_distance_LCDM(z)
    return D_C / (1 + z)


def angular_diameter_distance_Sync(z):
    """
    Angular diameter distance in Synchronism.

    Since D_C is unchanged and geometry is unchanged,
    D_A is also unchanged.
    """
    D_C = comoving_distance_Sync(z)
    return D_C / (1 + z)


def F_AP_LCDM(z):
    """
    Alcock-Paczynski parameter F_AP in ΛCDM.

    F_AP = D_A(z) × H(z) / c

    This quantity is dimensionless and depends on cosmology.
    """
    D_A = angular_diameter_distance_LCDM(z)
    H_z = H_LCDM(z)
    return D_A * H_z / c


def F_AP_Sync(z):
    """
    Alcock-Paczynski parameter in Synchronism.

    Since both D_A and H(z) are unchanged, F_AP is unchanged.
    """
    D_A = angular_diameter_distance_Sync(z)
    H_z = H_Sync(z)
    return D_A * H_z / c


def epsilon_parallel(z):
    """
    Radial stretch factor for DESI BAO.

    ε_∥ = H_fid(z) / H(z)

    For Synchronism vs ΛCDM: ε_∥ = 1.0 exactly.
    """
    return H_LCDM(z) / H_Sync(z)


def epsilon_perp(z):
    """
    Transverse stretch factor for DESI BAO.

    ε_⊥ = D_A(z) / D_A_fid(z)

    For Synchronism vs ΛCDM: ε_⊥ = 1.0 exactly.
    """
    return angular_diameter_distance_Sync(z) / angular_diameter_distance_LCDM(z)


def alpha_parallel(z, r_d_fid=147.09, r_d=147.09):
    """
    DESI α_∥ parameter.

    α_∥ = (D_H(z) / r_d) / (D_H_fid(z) / r_d_fid)

    where D_H = c / H(z)
    """
    D_H_Sync = c / H_Sync(z)
    D_H_LCDM = c / H_LCDM(z)
    return (D_H_Sync / r_d) / (D_H_LCDM / r_d_fid)


def alpha_perp(z, r_d_fid=147.09, r_d=147.09):
    """
    DESI α_⊥ parameter.

    α_⊥ = (D_M(z) / r_d) / (D_M_fid(z) / r_d_fid)

    where D_M = (1+z) × D_A(z) = D_C(z)
    """
    D_M_Sync = comoving_distance_Sync(z)
    D_M_LCDM = comoving_distance_LCDM(z)
    return (D_M_Sync / r_d) / (D_M_LCDM / r_d_fid)


def analyze_AP_predictions():
    """
    Analyze Synchronism predictions for Alcock-Paczynski test.
    """
    print("=" * 70)
    print("Session #116: Alcock-Paczynski Test Predictions for Synchronism")
    print("=" * 70)

    print("\n" + "-" * 70)
    print("1. THE ALCOCK-PACZYNSKI TEST")
    print("-" * 70)

    print("""
The AP test measures the ratio of angular to radial sizes of objects
that are intrinsically isotropic (like BAO or void stacks).

Any anisotropy in the measured correlation function indicates
incorrect cosmological parameters.

AP Parameter: F_AP = D_A(z) × H(z) / c

This is PURELY GEOMETRIC - it tests the expansion history, not growth.
    """)

    print("\n" + "-" * 70)
    print("2. SYNCHRONISM PREDICTION: F_AP UNCHANGED")
    print("-" * 70)

    print("""
KEY THEORETICAL RESULT:

In Synchronism, the modified Friedmann equation is:

    H² = (8πG/3C) × ρ_m

With C_cosmic = Ω_m(z), this EXACTLY reproduces ΛCDM expansion:

    H(z)_Sync = H(z)_LCDM

Therefore:
- Comoving distance D_C(z) is UNCHANGED
- Angular diameter distance D_A(z) is UNCHANGED
- AP parameter F_AP is UNCHANGED

Synchronism modifies GROWTH, not GEOMETRY.
    """)

    print("\n" + "-" * 70)
    print("3. NUMERICAL VERIFICATION")
    print("-" * 70)

    print("\n| z | H_LCDM | H_Sync | D_A_LCDM | D_A_Sync | F_AP_LCDM | F_AP_Sync |")
    print("|---|--------|--------|----------|----------|-----------|-----------|")

    for z in [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5]:
        H_L = H_LCDM(z)
        H_S = H_Sync(z)
        D_A_L = angular_diameter_distance_LCDM(z)
        D_A_S = angular_diameter_distance_Sync(z)
        F_L = F_AP_LCDM(z)
        F_S = F_AP_Sync(z)
        print(f"| {z:.1f} | {H_L:.1f} | {H_S:.1f} | {D_A_L:.0f} | {D_A_S:.0f} | {F_L:.4f} | {F_S:.4f} |")

    print("\nDeviation from ΛCDM:")
    print("\n| z | ΔH/H | ΔD_A/D_A | ΔF_AP/F_AP |")
    print("|---|------|----------|------------|")

    for z in [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5]:
        dH = (H_Sync(z) - H_LCDM(z)) / H_LCDM(z) * 100
        dD = (angular_diameter_distance_Sync(z) - angular_diameter_distance_LCDM(z)) / angular_diameter_distance_LCDM(z) * 100
        dF = (F_AP_Sync(z) - F_AP_LCDM(z)) / F_AP_LCDM(z) * 100
        print(f"| {z:.1f} | {dH:.4f}% | {dD:.4f}% | {dF:.4f}% |")

    print("\n" + "-" * 70)
    print("4. DESI BAO + AP PREDICTIONS")
    print("-" * 70)

    print("""
DESI measures BAO in both radial and transverse directions:

- α_∥ = (D_H / r_d) / (D_H_fid / r_d_fid)  [radial]
- α_⊥ = (D_M / r_d) / (D_M_fid / r_d_fid)  [transverse]

For Synchronism with fiducial = ΛCDM:
    """)

    print("\n| z | α_∥ | α_⊥ | ε_∥ | ε_⊥ |")
    print("|---|-----|-----|-----|-----|")

    for z in [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
        a_par = alpha_parallel(z)
        a_perp = alpha_perp(z)
        e_par = epsilon_parallel(z)
        e_perp = epsilon_perp(z)
        print(f"| {z:.1f} | {a_par:.4f} | {a_perp:.4f} | {e_par:.4f} | {e_perp:.4f} |")

    print("""
Result: ALL parameters = 1.0000 exactly.

This is a NULL PREDICTION: Synchronism predicts PERFECT agreement
with ΛCDM for all geometric/AP measurements.
    """)

    print("\n" + "-" * 70)
    print("5. WHY THIS IS IMPORTANT")
    print("-" * 70)

    print("""
The AP test is a CONSISTENCY CHECK for Synchronism:

1. If DESI BAO + AP shows deviation from ΛCDM:
   → Does NOT support Synchronism
   → Synchronism predicts EXACT agreement

2. If DESI BAO + AP agrees with ΛCDM:
   → CONSISTENT with Synchronism
   → But also consistent with ΛCDM (not discriminating)

The discriminating tests are:
- fσ8 (growth rate) - Session #103
- S8 (structure amplitude) - Session #102
- ISW effect - Session #104

These measure GROWTH, which IS modified by Synchronism.
BAO/AP measures GEOMETRY, which is UNCHANGED.
    """)

    print("\n" + "-" * 70)
    print("6. COMPARISON TO OTHER MODIFIED GRAVITY THEORIES")
    print("-" * 70)

    print("""
Different theories modify different things:

Theory        | H(z) | Growth | AP | fσ8
--------------|------|--------|----|---------
ΛCDM          | Std  | Std    | Std| Std
f(R)          | Mod  | Mod    | Mod| Mod
DGP           | Mod  | Mod    | Mod| Mod
Quintessence  | Mod  | Std    | Mod| Std
Synchronism   | Std  | Mod    | Std| Mod
              ↑                    ↑
           UNCHANGED            MODIFIED

KEY DISTINCTION:
Synchronism is UNIQUE in having:
- UNCHANGED expansion (H, D_A, BAO, AP all match ΛCDM)
- MODIFIED growth (fσ8, S8, σ8 all suppressed)

This makes it clearly distinguishable from f(R), DGP, etc.
    """)

    print("\n" + "-" * 70)
    print("7. FALSIFICATION CRITERIA")
    print("-" * 70)

    print("""
To FALSIFY Synchronism with AP/BAO:

1. If α_∥ ≠ 1.0 at any z:
   → Synchronism is FALSIFIED (predicts exact 1.0)

2. If α_⊥ ≠ 1.0 at any z:
   → Synchronism is FALSIFIED (predicts exact 1.0)

3. If F_AP deviates from ΛCDM:
   → Synchronism is FALSIFIED

Current DESI precision: ~1% on α parameters
Expected precision by 2030: ~0.3%

If DESI finds ANY deviation from ΛCDM in BAO/AP, Synchronism is ruled out.
    """)

    print("\n" + "-" * 70)
    print("8. INTEGRATION WITH PREVIOUS SESSIONS")
    print("-" * 70)

    print("""
Complete prediction picture (Sessions #102-116):

Category     | Observable      | Sync/ΛCDM | Discriminating?
-------------|-----------------|-----------|----------------
GROWTH       | S8              | 0.94      | YES (3.5σ current)
             | fσ8 (z=0.5)     | 0.88      | YES (3.3σ current)
             | σ8              | 0.94      | YES
             | γ (growth idx)  | 0.73      | YES
GEOMETRY     | H(z)            | 1.00      | NO
             | D_A(z)          | 1.00      | NO
             | F_AP            | 1.00      | NO
             | BAO scale       | 1.00      | NO
CROSS        | ISW             | 1.23      | MARGINAL
             | ISW/κg ratio    | 1.31      | YES (by 2030)

Synchronism makes a CLEAN SPLIT:
- Growth: MODIFIED (suppressed by ~6-12%)
- Geometry: UNCHANGED (identical to ΛCDM)

This is the UNIQUE SIGNATURE of the theory.
    """)

    print("\n" + "-" * 70)
    print("9. SUMMARY")
    print("-" * 70)

    print("""
Session #116 Key Findings:

1. THEORETICAL: AP test is UNCHANGED from ΛCDM
   - H(z) identical by construction (C_cosmic = Ω_m(z))
   - D_A(z) identical (geometry unchanged)
   - F_AP = D_A × H / c identical

2. NUMERICAL: α_∥ = α_⊥ = 1.0000 exactly
   - No deviation at any redshift
   - BAO scale preserved

3. DISCRIMINATING POWER: ZERO (for Synchronism vs ΛCDM)
   - AP is a CONSISTENCY CHECK
   - Any deviation would FALSIFY Synchronism

4. UNIQUE SIGNATURE CONFIRMED:
   - Synchronism: Growth modified, Geometry unchanged
   - f(R)/DGP: Both modified
   - This is the key discriminator

5. OBSERVATIONAL STRATEGY:
   - Use BAO/AP to RULE OUT Synchronism (if deviations found)
   - Use fσ8/S8 to CONFIRM Synchronism (growth suppression)
   - Combined analysis maximizes power

BOTTOM LINE: The AP test does not discriminate Synchronism from ΛCDM,
but any deviation would falsify the theory. This is a stringent
consistency requirement that Synchronism must pass.
    """)

    return {
        'F_AP': 'UNCHANGED',
        'alpha_parallel': 1.0,
        'alpha_perp': 1.0,
        'discriminating': False,
        'falsification': 'Any BAO/AP deviation rules out Synchronism'
    }


def create_visualization():
    """Create visualization of AP predictions."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    z = np.linspace(0.1, 3.0, 100)

    # 1. H(z) comparison
    ax1 = axes[0, 0]
    H_L = [H_LCDM(zi) for zi in z]
    H_S = [H_Sync(zi) for zi in z]

    ax1.plot(z, H_L, 'b-', linewidth=2, label=r'$H(z)$ ΛCDM')
    ax1.plot(z, H_S, 'r--', linewidth=2, label=r'$H(z)$ Synchronism')
    ax1.set_xlabel('Redshift z', fontsize=12)
    ax1.set_ylabel('H(z) [km/s/Mpc]', fontsize=12)
    ax1.set_title('Hubble Parameter: IDENTICAL', fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. D_A(z) comparison
    ax2 = axes[0, 1]
    D_A_L = [angular_diameter_distance_LCDM(zi) for zi in z]
    D_A_S = [angular_diameter_distance_Sync(zi) for zi in z]

    ax2.plot(z, D_A_L, 'b-', linewidth=2, label=r'$D_A(z)$ ΛCDM')
    ax2.plot(z, D_A_S, 'r--', linewidth=2, label=r'$D_A(z)$ Synchronism')
    ax2.set_xlabel('Redshift z', fontsize=12)
    ax2.set_ylabel(r'$D_A(z)$ [Mpc]', fontsize=12)
    ax2.set_title('Angular Diameter Distance: IDENTICAL', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. F_AP comparison
    ax3 = axes[1, 0]
    F_L = [F_AP_LCDM(zi) for zi in z]
    F_S = [F_AP_Sync(zi) for zi in z]

    ax3.plot(z, F_L, 'b-', linewidth=2, label=r'$F_{AP}$ ΛCDM')
    ax3.plot(z, F_S, 'r--', linewidth=2, label=r'$F_{AP}$ Synchronism')
    ax3.set_xlabel('Redshift z', fontsize=12)
    ax3.set_ylabel(r'$F_{AP} = D_A H / c$', fontsize=12)
    ax3.set_title('AP Parameter: IDENTICAL', fontsize=12)
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Ratio plot (showing they're identical)
    ax4 = axes[1, 1]
    ratio_H = [H_Sync(zi) / H_LCDM(zi) for zi in z]
    ratio_DA = [angular_diameter_distance_Sync(zi) / angular_diameter_distance_LCDM(zi) for zi in z]
    ratio_F = [F_AP_Sync(zi) / F_AP_LCDM(zi) for zi in z]

    ax4.plot(z, ratio_H, 'b-', linewidth=2, label=r'$H_{Sync}/H_{ΛCDM}$')
    ax4.plot(z, ratio_DA, 'g--', linewidth=2, label=r'$D_{A,Sync}/D_{A,ΛCDM}$')
    ax4.plot(z, ratio_F, 'r:', linewidth=2, label=r'$F_{AP,Sync}/F_{AP,ΛCDM}$')
    ax4.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Redshift z', fontsize=12)
    ax4.set_ylabel('Ratio', fontsize=12)
    ax4.set_title('All Ratios = 1.0 (Unchanged)', fontsize=12)
    ax4.set_ylim(0.98, 1.02)
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # Add annotation
    ax4.annotate('Synchronism: Geometry UNCHANGED\nGrowth MODIFIED',
                xy=(1.5, 0.99), fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session116_alcock_paczynski.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nVisualization saved to session116_alcock_paczynski.png")


if __name__ == "__main__":
    results = analyze_AP_predictions()
    create_visualization()

    print("\n" + "=" * 70)
    print("SESSION #116 COMPLETE")
    print("=" * 70)
    print(f"\nResults: {results}")
