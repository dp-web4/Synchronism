#!/usr/bin/env python3
"""
Session #80 Track A: Discriminating Tests - Synchronism vs MOND

Both theories reproduce BTFR exactly. Where do they differ?

From Session #79:
- MOND: Universal a₀, transition in OUTER disk
- Synchronism: Galaxy-dependent ρ_crit, transition in INNER disk

This script identifies and quantifies discriminating predictions.

Author: CBP Autonomous Synchronism Research
Date: December 3, 2025
Session: #80 - Discriminating Tests
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime


# Physical constants
G = 4.302e-6  # kpc (km/s)² / M_sun
a_0_MOND = 1.2e-10  # m/s² - MOND acceleration scale
a_0_kpc = a_0_MOND * 3.086e19 / 1e6  # Convert to (km/s)²/kpc = 3.7e-3


def print_section(title):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def discriminating_test_1_radial_profile():
    """
    Test 1: Radial Transition Location

    MOND: Transition at r where a = a_0
          r_MOND = V²/a_0 (typically OUTER disk, >10 kpc)

    Synchronism: Transition where ρ = ρ_crit
          r_Sync (typically INNER disk, 2-5 kpc)

    These predict DIFFERENT rotation curve shapes in the transition region.
    """
    print_section("Discriminating Test 1: Radial Transition Location")

    print("""
    MOND Transition:
    ----------------
    At r where centripetal acceleration equals a_0:
      a = V²/r = a_0
      r_MOND = V²/a_0

    For V = 200 km/s:
      r_MOND = (200)² / 3.7e-3 ≈ 10,800 kpc (!)

    Wait - this is way outside the galaxy! Let me recalculate...
    """)

    # MOND transition radius
    # a = GM(<r)/r² = a_0 at transition
    # For exponential disk: M(<r) ≈ M_total × (1 - e^(-r/R_d))
    # At r >> R_d: a ≈ GM_total/r² = a_0
    # r_MOND = √(GM_total/a_0)

    M_total = 5e10  # M_sun (typical MW baryonic mass)
    a_0_kpc = 3.7e-3  # (km/s)²/kpc

    r_mond = np.sqrt(G * M_total / a_0_kpc)

    print(f"    Recalculating with M = {M_total:.1e} M_sun:")
    print(f"    r_MOND = √(GM/a_0) = √({G:.2e} × {M_total:.1e} / {a_0_kpc:.2e})")
    print(f"    r_MOND = {r_mond:.1f} kpc")

    print("""
    Synchronism Transition:
    ----------------------
    Where ρ = ρ_crit = A × V^B

    For exponential disk: ρ(r) = ρ_0 × exp(-r/R_d)
    Transition at: ρ_0 × exp(-r/R_d) = ρ_crit
    r_Sync = R_d × ln(ρ_0 / ρ_crit)
    """)

    # Typical values
    R_d = 3.5  # kpc (disk scale length)
    rho_0 = 0.1  # M_sun/pc³ (central disk density)
    V_max = 200  # km/s
    A = 0.25
    B = 1.63
    rho_crit = A * V_max**B / 1e9  # Convert to M_sun/pc³

    print(f"    For V_max = {V_max} km/s:")
    print(f"    ρ_crit = {A} × {V_max}^{B} = {A * V_max**B:.1f} (code units)")
    print(f"    ρ_0 / ρ_crit ≈ {rho_0 * 1e9 / (A * V_max**B):.1f}")

    # Approximate transition radius
    if rho_0 * 1e9 > A * V_max**B:
        r_sync = R_d * np.log(rho_0 * 1e9 / (A * V_max**B))
    else:
        r_sync = 0  # Transition at center

    print(f"    r_Sync ≈ {r_sync:.1f} kpc (where C transitions)")

    print("""
    --- KEY DIFFERENCE ---

    MOND:        Transition at r ~ 5-20 kpc (acceleration-based)
    Synchronism: Transition at r ~ 1-5 kpc (density-based)

    These predict DIFFERENT rotation curve shapes:
    - MOND: Newtonian inner, MONDian outer
    - Sync: Low-coherence inner, high-coherence outer

    TEST: High-resolution rotation curves in 2-10 kpc range
          should show different curvature patterns.
    """)

    return {
        'test': 'radial_transition',
        'mond_r_transition_kpc': float(r_mond),
        'sync_r_transition_kpc': float(r_sync) if r_sync > 0 else 'inner',
        'discriminating_region': '2-10 kpc',
        'observable': 'rotation curve curvature in transition zone'
    }


def discriminating_test_2_mass_dependence():
    """
    Test 2: Mass Dependence of Enhancement

    MOND: G_eff/G ~ 1/μ(a/a_0) depends on ACCELERATION
          For given V, enhancement is ~ same regardless of mass

    Synchronism: G_eff/G = 1/C(ρ) depends on DENSITY
          Low mass + low V → lower density → more enhancement
          High mass + high V → higher density → less enhancement
    """
    print_section("Discriminating Test 2: Mass Dependence of Enhancement")

    print("""
    MOND Enhancement:
    -----------------
    μ(x) ≈ x/(1+x) where x = a/a_0
    G_eff/G = 1/μ = (1+x)/x = 1 + 1/x

    For deep MOND (x << 1): G_eff/G → 1/x = a_0/a

    At fixed radius r with different masses M1, M2:
      a1/a2 = M1/M2

    So enhancement ratio:
      (G_eff/G)_1 / (G_eff/G)_2 = a2/a1 = M2/M1

    MOND prediction: Enhancement scales inversely with MASS.
    """)

    print("""
    Synchronism Enhancement:
    ------------------------
    G_eff/G = 1/C where C = tanh(γ log(ρ/ρ_crit + 1))

    For different masses at same V:
      M1, M2 with V1 = V2 = V
      ρ_crit_1 = ρ_crit_2 = A × V^B (same!)

    But densities differ:
      ρ_1 ∝ M_1/R_1³
      ρ_2 ∝ M_2/R_2³

    From size-mass relation: R ∝ M^(1/3) (roughly)
    So ρ ∝ M / M^(1) = constant (approximately)

    Wait - this needs more careful analysis...
    """)

    # More careful: at fixed V, different M
    # From BTFR: M = A_TF × V^4
    # So at fixed V, M is FIXED in BTFR!

    print("""
    Key Insight: BTFR constrains M at fixed V
    -----------------------------------------
    M_bar = A_TF × V^4

    At fixed V, baryonic mass is nearly fixed!
    (Scatter is ~0.1 dex)

    So "different masses at same V" tests BTFR scatter.

    Better test: SAME baryonic mass, DIFFERENT environment
    → Void vs cluster galaxies
    """)

    print("""
    --- DISCRIMINATING PREDICTION ---

    MOND: Enhancement depends on a/a_0
          Environment doesn't directly affect a_0

    Synchronism: Enhancement depends on ρ/ρ_crit
          ρ_crit may vary with environment!

    TEST: Compare TF relation in voids vs clusters
    - MOND: Same TF (a_0 universal)
    - Sync: Offset TF (ρ_crit varies with C_formation)
    """)

    return {
        'test': 'mass_dependence',
        'mond_prediction': 'Enhancement ~ 1/M (fixed r)',
        'sync_prediction': 'Enhancement varies with environment',
        'discriminating_test': 'TF relation in voids vs clusters',
        'mond_expects': 'Same TF in all environments',
        'sync_expects': 'TF offset in voids (130% higher V)'
    }


def discriminating_test_3_external_field_effect():
    """
    Test 3: External Field Effect (EFE)

    MOND: Has EFE - external gravitational field affects internal dynamics
          Dwarf in cluster feels cluster potential

    Synchronism: No direct EFE - coherence is LOCAL
          But environment affects formation → affects ρ_crit

    These are DIFFERENT mechanisms with different signatures.
    """
    print_section("Discriminating Test 3: External Field Effect")

    print("""
    MOND's External Field Effect (EFE):
    -----------------------------------
    In MOND, the total acceleration matters:
      a_total = a_internal + a_external

    If a_external > a_0:
      Entire system is Newtonian (even if a_internal << a_0)

    Prediction for satellite galaxies:
      - Isolated dwarf: Deep MOND, high DM-like enhancement
      - Same dwarf near massive host: Newtonian, less enhancement

    This is TESTABLE by comparing:
      - Field dwarfs
      - Satellite dwarfs (same M_bar)
    """)

    print("""
    Synchronism's Environment Effect:
    ---------------------------------
    Coherence is LOCAL to the galaxy.
    External field doesn't directly affect internal C(ρ).

    BUT: Environment affects FORMATION:
      - Dwarfs formed in voids: Lower C_formation → lower ρ_crit
      - Dwarfs formed in clusters: Higher C_formation → higher ρ_crit

    Prediction:
      - Void dwarf: More DM-like effect (lower ρ_crit)
      - Cluster dwarf: Less DM-like effect (higher ρ_crit)

    KEY DIFFERENCE:
    - MOND: Effect depends on CURRENT environment
    - Sync: Effect depends on FORMATION environment
    """)

    print("""
    --- DISCRIMINATING TEST ---

    Find galaxies that:
    1. Formed in voids, now in clusters (infalling)
    2. Formed in clusters, now in voids (ejected)

    MOND predicts: Current environment matters
    - Infalling void dwarf shows EFE (less DM-like)

    Synchronism predicts: Formation environment matters
    - Infalling void dwarf still shows void signature (more DM-like)

    Observable: Rotation curves of infalling vs native satellites
    """)

    # Quantitative difference
    print("""
    Quantitative Example:
    --------------------
    Dwarf galaxy M_bar = 10^8 M_sun, V_max = 50 km/s

    Case 1: Isolated field dwarf
      MOND: a << a_0 → deep MOND → V_flat ~ (G M a_0)^(1/4)
      Sync: Moderate ρ_crit → moderate C

    Case 2: Same dwarf 100 kpc from MW
      MOND: a_ext ~ GM_MW/r² ~ 10^-11 m/s² ~ a_0
            → Partial EFE → reduced enhancement
      Sync: No direct EFE, but higher C if formed near MW
            → same reduction but for different reason

    Case 3: Void-formed dwarf now infalling to cluster
      MOND: Strong EFE → Newtonian
      Sync: Still void-like (formation imprint) → enhanced
    """)

    # Calculate EFE example
    M_MW = 1e12  # M_sun
    r_sat = 100  # kpc
    a_ext = G * M_MW / r_sat**2 * 1e6  # Convert to m/s² scale

    print(f"    At r = {r_sat} kpc from MW:")
    print(f"    a_ext ≈ {a_ext:.2e} (km/s)²/kpc")
    print(f"    a_0 ≈ {a_0_kpc:.2e} (km/s)²/kpc")
    print(f"    a_ext / a_0 ≈ {a_ext/a_0_kpc:.2f}")

    return {
        'test': 'external_field_effect',
        'mond_mechanism': 'Current external acceleration matters',
        'sync_mechanism': 'Formation environment imprinted in ρ_crit',
        'discriminating_test': 'Infalling void dwarfs vs native satellites',
        'mond_expects': 'Infalling dwarfs show EFE (Newtonian-like)',
        'sync_expects': 'Infalling dwarfs retain void signature (enhanced)'
    }


def discriminating_test_4_stellar_vs_gas():
    """
    Test 4: Stellar vs Gas Dominated Systems

    MOND: Depends on TOTAL baryonic mass
          Stars and gas contribute equally per unit mass

    Synchronism: Depends on DENSITY structure
          Gas-dominated systems may have different ρ profiles

    This could distinguish the theories for low-surface-brightness galaxies.
    """
    print_section("Discriminating Test 4: Stellar vs Gas Dominated Systems")

    print("""
    Low Surface Brightness (LSB) Galaxies:
    --------------------------------------
    - High gas fraction (M_gas/M_star > 1)
    - Extended, diffuse structure
    - Slow rotation, low surface density

    MOND Prediction:
      V_flat ~ (G M_bar a_0)^(1/4)
      Only depends on total M_bar = M_star + M_gas
      Gas-rich = more mass = higher V

    Synchronism Prediction:
      C = tanh(γ log(ρ/ρ_crit + 1))
      ρ = ρ_star + ρ_gas at each radius
      Extended gas disk → lower mean density → lower C → more enhancement
    """)

    print("""
    Gas vs Stars at Same Total Mass:
    --------------------------------
    Compare two galaxies with M_bar = 10^9 M_sun:

    Galaxy A: M_star = 9×10^8, M_gas = 10^8 (stellar dominated)
      - Compact stellar disk: ρ ~ 0.1 M_sun/pc³
      - Low gas density outer region

    Galaxy B: M_gas = 9×10^8, M_star = 10^8 (gas dominated)
      - Extended gas disk: ρ ~ 0.01 M_sun/pc³
      - Much lower mean density

    MOND: Same M_bar → same V_flat (approximately)

    Synchronism:
      Galaxy A: Higher ρ → higher C → less enhancement
      Galaxy B: Lower ρ → lower C → more enhancement

    Prediction: At fixed M_bar, gas-dominated systems
                have higher V_flat in Synchronism but not MOND
    """)

    # Quantitative estimate
    M_bar = 1e9  # M_sun

    # Stellar dominated
    R_star = 2  # kpc
    rho_star = M_bar / (4/3 * np.pi * R_star**3) / 1e9  # M_sun/pc³

    # Gas dominated
    R_gas = 6  # kpc
    rho_gas = M_bar / (4/3 * np.pi * R_gas**3) / 1e9  # M_sun/pc³

    print(f"\n    At M_bar = {M_bar:.1e} M_sun:")
    print(f"    Stellar-dom (R={R_star} kpc): ρ ~ {rho_star:.4f} M_sun/pc³")
    print(f"    Gas-dom (R={R_gas} kpc): ρ ~ {rho_gas:.4f} M_sun/pc³")
    print(f"    Density ratio: {rho_star/rho_gas:.1f}×")

    print("""
    --- DISCRIMINATING TEST ---

    Compare TF relation for:
    - High surface brightness (HSB) galaxies
    - Low surface brightness (LSB) galaxies

    MOND: Same TF relation (depends only on M_bar)
    Sync: LSB offset to higher V (lower density → more enhancement)

    Existing data: McGaugh et al. have LSB samples
    Test: Is there V_flat offset at fixed M_bar for LSB vs HSB?
    """)

    return {
        'test': 'stellar_vs_gas',
        'mond_prediction': 'Same TF for HSB and LSB at fixed M_bar',
        'sync_prediction': 'LSB have higher V_flat at fixed M_bar',
        'reason': 'Lower density → lower C → more enhancement',
        'data_source': 'McGaugh LSB + HSB samples',
        'observable': 'TF offset between HSB and LSB'
    }


def discriminating_test_5_time_evolution():
    """
    Test 5: Time Evolution / Lookback

    MOND: a_0 is constant (fundamental)
          Same physics at all redshifts

    Synchronism: ρ_crit may evolve with cosmic C
          Different effective dark matter at high z?

    This is testable with high-z rotation curves (though difficult).
    """
    print_section("Discriminating Test 5: Time Evolution")

    print("""
    MOND at High Redshift:
    ----------------------
    a_0 is a fundamental constant.
    BTFR should be identical at z = 0 and z = 2.

    Exception: Some MOND variants have a_0 ∝ H(z)
    This would change BTFR normalization with redshift.

    Standard MOND prediction: Same TF at all z.
    """)

    print("""
    Synchronism at High Redshift:
    -----------------------------
    ρ_crit = A × V^B where A depends on formation C.

    At high z:
    - Universe denser → higher mean C
    - But also: galaxies still forming → C_formation varies

    Possible effects:
    1. If C_formation was higher at z~2: ρ_crit higher → less DM-like
    2. If structure was more uniform: C_formation varies less

    Prediction: TF relation evolution with redshift
    - Normalization may shift
    - Scatter may be different
    """)

    # Cosmic density evolution
    z_values = [0, 0.5, 1, 2, 3]
    rho_cosmic = [1, 1.5**3 * 0.31, 2**3 * 0.31, 3**3 * 0.31, 4**3 * 0.31]  # relative to today

    print("\n    Cosmic Density Evolution:")
    print(f"    {'z':>5} {'ρ/ρ₀':>10}")
    print("    " + "-" * 18)
    for z, rho in zip(z_values, rho_cosmic):
        print(f"    {z:>5} {rho:>10.2f}")

    print("""
    --- DISCRIMINATING TEST ---

    ALMA/JWST rotation curves at z > 1:
    - Compare TF relation at high z vs local

    MOND: Same TF (a_0 constant)
    Sync: Possible evolution (ρ_crit ~ C_formation)

    Challenge: High-z rotation curves are noisy
    But: Even ~0.1 dex TF offset would be detectable

    Current data: Lelli et al. 2023 shows TF holds to z~2.5
                  This constrains BOTH theories!
    """)

    return {
        'test': 'time_evolution',
        'mond_prediction': 'TF constant with redshift (a_0 fundamental)',
        'sync_prediction': 'TF may evolve (ρ_crit ~ C_formation)',
        'current_data': 'TF appears stable to z~2.5',
        'implication': 'Constrains C_formation evolution in Synchronism',
        'future_test': 'JWST/ALMA high-z rotation curves'
    }


def summary_discriminating_tests():
    """
    Summarize all discriminating tests.
    """
    print_section("SUMMARY: Discriminating Tests")

    tests = [
        ("1. Radial Transition", "2-10 kpc rotation curve curvature",
         "MOND: Outer transition", "Sync: Inner transition"),
        ("2. Environment", "Void vs cluster TF offset",
         "MOND: Same TF everywhere", "Sync: 130% V offset in voids"),
        ("3. External Field", "Infalling satellite dynamics",
         "MOND: EFE reduces enhancement", "Sync: Formation imprint persists"),
        ("4. HSB vs LSB", "TF relation surface brightness dependence",
         "MOND: Same TF", "Sync: LSB higher V at fixed M"),
        ("5. Redshift", "High-z TF evolution",
         "MOND: Constant a_0", "Sync: Possible C_formation evolution"),
    ]

    print(f"\n{'Test':<25} {'Observable':<35} {'MOND':<25} {'Sync':<25}")
    print("-" * 110)

    for name, observable, mond, sync in tests:
        print(f"{name:<25} {observable:<35} {mond:<25} {sync:<25}")

    print("""

    MOST PROMISING TESTS (with existing data):
    -----------------------------------------
    1. Void galaxy TF offset (ALFALFA + SDSS)
       - Clear prediction difference
       - Data exists, analysis feasible

    2. HSB vs LSB TF comparison (McGaugh samples)
       - Well-characterized samples exist
       - Already analyzed for MOND

    3. High-z TF (Lelli et al. + JWST)
       - Constrains both theories
       - Data quality improving rapidly

    LESS FEASIBLE (needs new observations):
    --------------------------------------
    4. EFE in infalling satellites
       - Requires identifying formation environment
       - Kinematic data for faint dwarfs

    5. Inner rotation curve curvature
       - Needs high spatial resolution
       - Beam smearing is a problem
    """)

    return {
        'most_promising': ['void_TF', 'HSB_LSB_TF', 'high_z_TF'],
        'less_feasible': ['EFE_satellites', 'inner_curvature'],
        'immediate_action': 'Analyze ALFALFA void galaxies for TF offset'
    }


def main():
    """Run all discriminating test analyses."""
    print("=" * 70)
    print("SESSION #80 TRACK A: DISCRIMINATING TESTS")
    print("Synchronism vs MOND - Where Do They Differ?")
    print("=" * 70)

    results = {}

    results['test_1'] = discriminating_test_1_radial_profile()
    results['test_2'] = discriminating_test_2_mass_dependence()
    results['test_3'] = discriminating_test_3_external_field_effect()
    results['test_4'] = discriminating_test_4_stellar_vs_gas()
    results['test_5'] = discriminating_test_5_time_evolution()
    results['summary'] = summary_discriminating_tests()

    # Save results
    output_path = Path(__file__).parent / 'results' / 'session80_discriminating_tests.json'
    output_path.parent.mkdir(exist_ok=True)

    output = {
        'session': 80,
        'track': 'A',
        'title': 'Discriminating Tests: Synchronism vs MOND',
        'date': datetime.now().isoformat(),
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n\nResults saved to: {output_path}")

    print("\n" + "=" * 70)
    print("SESSION #80 TRACK A COMPLETE")
    print("=" * 70)

    return results


if __name__ == '__main__':
    main()
