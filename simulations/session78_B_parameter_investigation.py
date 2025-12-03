#!/usr/bin/env python3
"""
Session #78 Track A: Investigating the B Parameter Discrepancy

The critical finding from Session #77:
- B_derived = 0.5 (from Jeans + R_half ∝ V^0.75 scaling)
- B_empirical = 1.62 (from SPARC rotation curve fitting)

This ~3x discrepancy is significant. This script investigates:
1. What galaxy scaling would give B = 1.62?
2. Is the R_half ∝ V^0.75 assumption correct?
3. What alternative physical derivations could yield B ≈ 1.6?

Author: CBP Autonomous Synchronism Research
Date: December 3, 2025
Session: #78 - B Parameter Investigation
"""

import numpy as np
import json
from datetime import datetime
from scipy import stats

print("=" * 70)
print("Session #78 Track A: B Parameter Discrepancy Investigation")
print("=" * 70)


# Physical constants
G = 4.3e-3  # kpc³ / (M_sun × (km/s)² × Gyr²) - convenient units
G_SI = 6.674e-11  # m³/(kg·s²)
M_sun = 1.989e30  # kg
pc = 3.086e16  # m


def analyze_scaling_for_B():
    """
    Analyze what galaxy scaling would give B = 1.62.
    """
    print("\n" + "=" * 70)
    print("PART 1: What Galaxy Scaling Gives B = 1.62?")
    print("=" * 70)

    print("""
DERIVATION REVIEW:

From Jeans criterion:
    ρ_crit = V² / (G α² R_half²)

If R_half = R_0 × V^δ, then:
    ρ_crit = V² / (G α² R_0² V^(2δ))
           = V^(2 - 2δ) / (G α² R_0²)
           = A × V^B

where B = 2 - 2δ

Current derivation: δ = 0.75 → B = 2 - 1.5 = 0.5
Empirical: B = 1.62 → δ = (2 - B) / 2 = (2 - 1.62) / 2 = 0.19

CRITICAL FINDING:
To get B = 1.62, we need R_half ∝ V^0.19 (almost constant!)
This is VERY different from R_half ∝ V^0.75
""")

    # Plot what this means
    v_range = np.linspace(30, 300, 100)

    # Different δ values
    delta_values = [0.75, 0.5, 0.19, 0.0]
    B_values = [2 - 2*d for d in delta_values]

    print("δ (R ∝ V^δ)  →  B (ρ_crit ∝ V^B)")
    print("-" * 40)
    for d, b in zip(delta_values, B_values):
        marker = " ← DERIVED" if abs(b - 0.5) < 0.1 else ""
        marker = " ← EMPIRICAL" if abs(b - 1.62) < 0.1 else marker
        print(f"δ = {d:.2f}  →  B = {b:.2f}{marker}")

    return {
        'derived_delta': 0.75,
        'derived_B': 0.5,
        'empirical_B': 1.62,
        'required_delta': 0.19,
        'interpretation': 'Empirical B requires nearly constant R_half vs V'
    }


def test_real_galaxy_scaling():
    """
    Test actual galaxy scaling relations from SPARC-like data.
    """
    print("\n" + "=" * 70)
    print("PART 2: Real Galaxy Scaling Relations")
    print("=" * 70)

    # SPARC-like galaxy data (approximate values)
    # Format: name, V_flat (km/s), R_eff (kpc), Type
    galaxies = [
        # Dwarfs
        ("DDO 154", 47, 1.5, "dwarf"),
        ("NGC 2366", 55, 2.0, "dwarf"),
        ("DDO 168", 60, 1.3, "dwarf"),
        ("IC 2574", 65, 6.5, "dwarf"),
        ("NGC 1569", 50, 0.8, "dwarf"),

        # Intermediate
        ("NGC 3109", 67, 2.5, "irregular"),
        ("NGC 2403", 136, 3.9, "spiral"),
        ("NGC 7793", 108, 2.2, "spiral"),
        ("NGC 3521", 212, 3.0, "spiral"),
        ("NGC 3031", 215, 3.2, "spiral"),

        # Massive
        ("NGC 7331", 250, 4.5, "spiral"),
        ("UGC 2953", 260, 5.5, "spiral"),
        ("NGC 4736", 155, 1.8, "spiral"),
        ("NGC 5055", 200, 3.8, "spiral"),
        ("NGC 6946", 175, 3.5, "spiral"),
    ]

    v_vals = np.array([g[1] for g in galaxies])
    r_vals = np.array([g[2] for g in galaxies])

    # Fit R = R_0 × V^δ
    log_v = np.log10(v_vals)
    log_r = np.log10(r_vals)

    slope, intercept, r_value, p_value, std_err = stats.linregress(log_v, log_r)

    delta_fit = slope
    R_0_fit = 10**intercept

    print(f"\nFitting R_half = R_0 × V^δ to {len(galaxies)} galaxies:")
    print(f"  δ (fitted) = {delta_fit:.3f} ± {std_err:.3f}")
    print(f"  R_0 (fitted) = {R_0_fit:.3f} kpc/(km/s)^δ")
    print(f"  R² = {r_value**2:.3f}")
    print(f"  Implied B = 2 - 2×{delta_fit:.3f} = {2 - 2*delta_fit:.3f}")

    # Compare to expectations
    print("\n" + "-" * 50)
    print("COMPARISON:")
    print("-" * 50)
    print(f"  Derived (Session #53): δ = 0.75 → B = 0.5")
    print(f"  This fit: δ = {delta_fit:.2f} → B = {2 - 2*delta_fit:.2f}")
    print(f"  Empirical target: B = 1.62 → δ = 0.19")

    return {
        'delta_fit': float(delta_fit),
        'implied_B': float(2 - 2*delta_fit),
        'r_squared': float(r_value**2),
        'n_galaxies': len(galaxies)
    }


def alternative_derivations():
    """
    Explore alternative physical derivations that could give B ≈ 1.6.
    """
    print("\n" + "=" * 70)
    print("PART 3: Alternative Physical Derivations for B")
    print("=" * 70)

    print("""
HYPOTHESIS 1: Dynamical Timescale Scaling
══════════════════════════════════════════════════════════════════════

If coherence depends on the DYNAMICAL TIME rather than Jeans length:
    τ_dyn = R / V

And the coherence condition is:
    ρ_crit = f(τ_dyn) = f(R/V)

For ρ_crit ∝ V^B, we need to know how R/V scales.

If R ∝ V^δ:
    R/V ∝ V^(δ-1)

For Tully-Fisher disk scaling (M ∝ V^4, ρ ∝ M/R³, Σ ~ const):
    R² ∝ M ∝ V^4 → R ∝ V^2 → δ = 2
    τ_dyn = R/V ∝ V^1 (dynamical time grows with V)

If ρ_crit ∝ τ_dyn^(-2) (decoherence rate²):
    ρ_crit ∝ (V/R)² = V^(2-2δ) = V^(2-4) = V^(-2)

This gives B = -2, wrong sign!


HYPOTHESIS 2: Surface Density Scaling
══════════════════════════════════════════════════════════════════════

Galaxy surface brightness is relatively constant (Freeman's law):
    Σ = M / R² ≈ constant for disk galaxies

If coherence depends on SURFACE density rather than volume density:
    Σ_crit = Σ_0 × V^B_Σ

And volume density: ρ ∝ Σ / h where h = disk scale height

If h ∝ R ∝ V^δ:
    ρ = Σ / h ∝ V^(B_Σ - δ)

For B = 1.62 with δ ≈ 0.5:
    B_Σ = B + δ = 1.62 + 0.5 = 2.12

This means: Σ_crit ∝ V^2.12 ≈ V²

Interpretation: Critical SURFACE density scales as V²!
    Σ_crit ∝ V² is related to centripetal acceleration:
    a_c = V²/R ∝ V^(2-δ)

This connects to MOND's a_0 acceleration scale!


HYPOTHESIS 3: Baryonic Tully-Fisher Scaling
══════════════════════════════════════════════════════════════════════

The baryonic Tully-Fisher relation: M_bar ∝ V^4 (exactly!)

If the coherence scale tracks baryonic mass:
    ρ_crit ∝ M_bar / R³ ∝ V^4 / V^(3δ) = V^(4 - 3δ)

For B = 1.62:
    4 - 3δ = 1.62
    δ = (4 - 1.62) / 3 = 0.79

This gives δ ≈ 0.8, close to the Session #53 assumption of 0.75!

WAIT - this is almost exactly what Session #53 derived!
The discrepancy might be in the FORMULA, not the scaling.
""")

    # Test BTFR-based formula
    delta_btfr = (4 - 1.62) / 3
    print(f"\nBTFR implies δ = {delta_btfr:.2f}")
    print(f"Session #53 assumed δ = 0.75")
    print(f"Difference: {abs(delta_btfr - 0.75):.2f}")

    return {
        'hypothesis_1': 'Dynamical time - gives B = -2 (wrong)',
        'hypothesis_2': 'Surface density - gives Σ_crit ∝ V², connects to MOND',
        'hypothesis_3': 'BTFR - gives δ ≈ 0.79, close to 0.75'
    }


def reexamine_jeans_derivation():
    """
    Re-examine the Jeans derivation more carefully.
    """
    print("\n" + "=" * 70)
    print("PART 4: Re-examining the Jeans Derivation")
    print("=" * 70)

    print("""
ORIGINAL DERIVATION (Session #53):
══════════════════════════════════════════════════════════════════════

Jeans criterion: λ_J = α × R_half
    where λ_J = V / √(G ρ)

Solving: ρ_crit = V² / (G α² R²)

With R = R_0 × V^δ:
    ρ_crit = V^(2-2δ) / (G α² R_0²)


THE ISSUE: What is R in this formula?
══════════════════════════════════════════════════════════════════════

Session #53 used R = R_half (stellar half-light radius).

But for ROTATION CURVES, what matters is the TOTAL MATTER distribution,
not just stars. The relevant scale might be:

1. R_disk (gas + stars disk scale length)
2. R_vir (virial radius of total halo)
3. R_DM (dark matter scale radius, if present)

If R_relevant ≠ R_half, the scaling could be different.


ALTERNATIVE: Use MOND-like reasoning
══════════════════════════════════════════════════════════════════════

In MOND, the transition occurs at acceleration a_0:
    a_MOND = a_0 ≈ 1.2 × 10^-10 m/s²

The centripetal acceleration is:
    a_c = V² / R

Setting a_c = a_crit:
    a_crit = V² / R = ρ_crit × G × R (approximately)

    Therefore: ρ_crit × R² ∝ V²

If R ∝ V^δ:
    ρ_crit ∝ V² / R² = V^(2-2δ)

Same formula! The issue is what δ actually is.


NEW APPROACH: What if δ is NOT constant?
══════════════════════════════════════════════════════════════════════

Maybe R/V relationship depends on galaxy mass/type:
    δ_dwarf ≠ δ_spiral ≠ δ_massive

If there's a systematic trend:
    δ_eff(V) = δ_0 + ε × log(V)

Then the effective B would be:
    B_eff = 2 - 2δ_eff(V) = 2 - 2δ_0 - 2ε × log(V)

This creates a LOG TERM in the ρ_crit formula!

Actually, this might explain the tanh form:
    tanh(γ × log(ρ/ρ_crit)) captures a log-scale transition
""")

    return {
        'issue': 'R in Jeans formula - which scale?',
        'possibilities': ['R_half (stars)', 'R_disk (total)', 'R_vir (halo)'],
        'insight': 'δ may not be constant across galaxy types'
    }


def numerical_investigation():
    """
    Numerically investigate what gives B = 1.62.
    """
    print("\n" + "=" * 70)
    print("PART 5: What Physical Formula Gives B ≈ 1.6?")
    print("=" * 70)

    # The empirical formula is: ρ_crit = A × V^B with B = 1.62

    print("""
SEARCH FOR PHYSICAL FORMULAE:
══════════════════════════════════════════════════════════════════════

We want ρ_crit ∝ V^1.62. What combination of known scalings gives this?

Known galaxy scaling relations:
    - Tully-Fisher: M ∝ V^4 (exact in MOND)
    - Size-velocity: R ∝ V^(0.3 to 1.0) depending on type
    - Mass-concentration: c ∝ M^(-0.1) (from NFW fits)
    - Faber-Jackson (ellipticals): L ∝ σ^4

Let's test combinations:
""")

    # Test different formula combinations
    formulas = [
        ("V²/R² with R ∝ V^0.19", lambda v: v**2 / (v**0.19)**2, 2 - 2*0.19),
        ("V²/R² with R ∝ V^0.75", lambda v: v**2 / (v**0.75)**2, 2 - 2*0.75),
        ("M/R³ with M ∝ V^4, R ∝ V^0.79", lambda v: v**4 / (v**0.79)**3, 4 - 3*0.79),
        ("V^4/R^3 with R ∝ V^1.5", lambda v: v**4 / (v**1.5)**3, 4 - 4.5),
        ("Surface density Σ ∝ V^2 / R with R ∝ V^0.4", lambda v: v**2 / v**0.4, 2 - 0.4),
    ]

    print(f"{'Formula':<50} {'Implied B':<12}")
    print("-" * 65)

    for name, func, B_implied in formulas:
        marker = " ← CLOSE!" if abs(B_implied - 1.62) < 0.2 else ""
        print(f"{name:<50} {B_implied:<12.2f}{marker}")

    # The key insight
    print("""

KEY FINDING:
══════════════════════════════════════════════════════════════════════

Surface density Σ ∝ V^2/R with R ∝ V^0.4 gives B = 1.6!

This means:
    Σ_crit ∝ V^2 / R ∝ V^2 / V^0.4 = V^1.6

Physical interpretation:
    - V^2 is the centripetal acceleration × R
    - Dividing by R gives surface density scale
    - R ∝ V^0.4 is MUCH flatter than R ∝ V^0.75

For LATE-TYPE SPIRALS specifically:
    The R-V relation is flatter than the general sample!

This suggests the B = 1.62 is calibrated to a specific galaxy population.
""")

    return {
        'finding': 'Σ_crit ∝ V^2/R with R ∝ V^0.4 gives B ≈ 1.6',
        'interpretation': 'Surface density interpretation, flat R-V relation'
    }


def summary():
    """
    Summarize findings.
    """
    print("\n" + "=" * 70)
    print("SESSION #78 TRACK A SUMMARY")
    print("=" * 70)

    summary_text = """
B PARAMETER DISCREPANCY: What We Learned
══════════════════════════════════════════════════════════════════════

DERIVED (Session #53): B = 0.5
    From: ρ_crit = V² / (G α² R²) with R ∝ V^0.75
    Gives: B = 2 - 2×0.75 = 0.5

EMPIRICAL: B = 1.62
    From: SPARC rotation curve fitting
    Implies: R ∝ V^0.19 (nearly constant R!)

THE GAP: B differs by ~3x (0.5 vs 1.62)


POSSIBLE EXPLANATIONS:
══════════════════════════════════════════════════════════════════════

1. WRONG R-V SCALING ASSUMED
   Session #53 used R ∝ V^0.75 (general galaxy population)
   But SPARC galaxies might have flatter R ∝ V^0.4 relation

2. WRONG SCALE (R) USED
   Half-light radius R_half ≠ relevant scale for coherence
   Should use total matter radius or coherence length

3. SURFACE DENSITY vs VOLUME DENSITY
   If coherence depends on Σ rather than ρ:
   Σ_crit ∝ V^2/R ∝ V^1.6 for R ∝ V^0.4

4. GALAXY POPULATION BIAS
   The empirical B is calibrated to SPARC late-type spirals
   The derived B might apply to different populations


NEXT STEPS:
══════════════════════════════════════════════════════════════════════

1. Check R-V scaling for SPARC galaxies specifically
2. Test if B varies with galaxy type (dwarfs vs spirals)
3. Investigate if surface density formulation works better
4. Consider variable B(type) rather than universal B
"""
    print(summary_text)

    return {
        'B_derived': 0.5,
        'B_empirical': 1.62,
        'explanation_1': 'Wrong R-V scaling (0.75 vs needed 0.19-0.4)',
        'explanation_2': 'Wrong scale R used (R_half vs coherence scale)',
        'explanation_3': 'Surface vs volume density',
        'explanation_4': 'Galaxy population bias'
    }


# Run all analyses
if __name__ == '__main__':
    scaling_result = analyze_scaling_for_B()
    galaxy_result = test_real_galaxy_scaling()
    alt_result = alternative_derivations()
    jeans_result = reexamine_jeans_derivation()
    numerical_result = numerical_investigation()
    summary_result = summary()

    # Save results
    results = {
        'session': 78,
        'track': 'A',
        'title': 'B Parameter Discrepancy Investigation',
        'date': datetime.now().isoformat(),

        'findings': {
            'scaling_analysis': scaling_result,
            'galaxy_fit': galaxy_result,
            'alternative_derivations': alt_result,
            'jeans_reexamination': jeans_result,
            'numerical_investigation': numerical_result,
            'summary': summary_result
        }
    }

    output_path = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/results/session78_B_investigation.json'
    import os
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {output_path}")
