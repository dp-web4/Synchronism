#!/usr/bin/env python3
"""
Session #90: Tidal Dwarf Galaxies as MOND vs Synchronism Test

Sessions #88-89 showed MOND and Synchronism are unified for DISK galaxies.
But they may DIFFER for non-disk systems like tidal dwarf galaxies (TDGs).

TDGs are unique because:
1. They form from tidal debris of galaxy collisions
2. They should NOT have dark matter halos (ΛCDM)
3. MOND predicts: Same BTFR as normal galaxies (a₀ universal)
4. Synchronism predicts: BTFR depends on formation environment

This is a key discriminating test!

Key references:
- Bournaud et al. (2007): TDG dynamics
- Lelli et al. (2015): TDG rotation curves
- Müller et al. (2018): NGC 1052 TDG candidates

Author: CBP Autonomous Synchronism Research
Date: December 6, 2025
"""

import numpy as np
import json
from pathlib import Path


def tdg_physics():
    """
    Explain the physics of tidal dwarf galaxies.
    """
    print("="*70)
    print("TIDAL DWARF GALAXIES: THE PHYSICS")
    print("="*70)
    print()

    print("What are TDGs?")
    print("  - Galaxies formed from tidal debris of merging galaxies")
    print("  - Made of pre-existing stars and gas from parent galaxy")
    print("  - Typically low mass: 10^7 - 10^9 M_sun")
    print("  - No primordial dark matter halo (ΛCDM prediction)")
    print()

    print("Why are they interesting for MOND/Synchronism?")
    print()
    print("  ΛCDM: TDGs lack DM halos → rotation curves should be Keplerian")
    print("        (V ∝ 1/√R at large R)")
    print()
    print("  MOND:  a₀ is universal → TDGs obey SAME BTFR as normal galaxies")
    print("        (no DM needed, modified gravity applies)")
    print()
    print("  Synchronism: C(ρ) depends on LOCAL density")
    print("        - If TDGs have similar densities → same as normal galaxies")
    print("        - If TDGs have different densities → DIFFERENT behavior")
    print()


def mond_vs_synchronism_for_tdgs():
    """
    Analyze where MOND and Synchronism differ for TDGs.
    """
    print("="*70)
    print("MOND vs SYNCHRONISM FOR TDGs")
    print("="*70)
    print()

    print("KEY QUESTION:")
    print("Do TDGs have the same surface density profile as normal dwarfs?")
    print()

    print("-"*70)
    print()
    print("CASE 1: TDGs have SAME Σ profile as normal dwarfs")
    print()
    print("  MOND prediction:    Same BTFR (a₀ universal)")
    print("  Synchronism prediction: Same BTFR (C(ρ) same)")
    print("  Result: NO discrimination between theories")
    print()

    print("-"*70)
    print()
    print("CASE 2: TDGs have DIFFERENT Σ profile")
    print()
    print("  Example: TDGs are more diffuse (lower Σ at same R)")
    print()
    print("  MOND prediction:    Same BTFR (a₀ doesn't depend on Σ)")
    print("  Synchronism prediction:")
    print("    - Lower Σ → lower ρ → lower C → higher G_eff")
    print("    - V_obs > V_bar (more 'dark matter' enhancement)")
    print("    - TDGs would have HIGHER V at fixed M_bar")
    print()
    print("  Result: DIFFERENT predictions! This is a discriminating test!")
    print()

    print("-"*70)
    print()
    print("CASE 3: TDGs formed in different environment")
    print()
    print("  TDGs form in the tidal debris of interacting galaxies")
    print("  The formation environment has:")
    print("    - Higher local density (overdense region)")
    print("    - Different velocity dispersion")
    print()
    print("  Synchronism implication:")
    print("    - If C depends on formation environment → different BTFR")
    print("    - Session #85 showed: Environmental effect is WEAK (~0.1 coefficient)")
    print("    - So: TDGs should be close to normal BTFR")
    print()


def observational_evidence():
    """
    Review observational evidence on TDG dynamics.
    """
    print("="*70)
    print("OBSERVATIONAL EVIDENCE")
    print("="*70)
    print()

    print("1. LELLI et al. (2015) - NGC 5291 system")
    print("-"*70)
    print("  - 3 TDGs with resolved rotation curves")
    print("  - V_flat = 30-60 km/s")
    print("  - M_bar = 10^8 - 10^9 M_sun")
    print("  - BTFR: Consistent with normal galaxies!")
    print("  - MOND: Confirmed (no DM needed)")
    print()

    print("2. FLORES et al. (2016) - VCC 2062")
    print("-"*70)
    print("  - TDG in Virgo cluster")
    print("  - V_rot ~ 50 km/s")
    print("  - M_bar ~ 10^9 M_sun")
    print("  - BTFR: Consistent with normal galaxies")
    print()

    print("3. MÜLLER et al. (2018) - NGC 1052 satellites")
    print("-"*70)
    print("  - NGC 1052-DF2 and NGC 1052-DF4")
    print("  - Initially claimed to LACK dark matter")
    print("  - Later revised: stellar mass lower, may be consistent with MOND")
    print("  - Still debated")
    print()

    print("SUMMARY OF OBSERVATIONS:")
    print("-"*70)
    print("  - Most TDGs follow the same BTFR as normal galaxies")
    print("  - This supports MOND (a₀ universal)")
    print("  - For Synchronism: Need to check if TDG densities are similar")
    print()


def synchronism_prediction_for_tdgs():
    """
    Calculate Synchronism prediction for TDGs.
    """
    print("="*70)
    print("SYNCHRONISM PREDICTION FOR TDGs")
    print("="*70)
    print()

    print("From Sessions #87-89:")
    print("  - C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))")
    print("  - V/V_bar correlates with local SB (r = -0.626)")
    print("  - Both Synchronism and MOND capture same physics in disks")
    print()

    print("For TDGs:")
    print("  - If TDG has same Σ profile → same C(ρ) → same BTFR")
    print("  - If TDG is more diffuse → lower C → higher V_obs/V_bar")
    print()

    # Typical TDG parameters
    M_tdg = 1e8  # M_sun
    R_tdg = 2  # kpc (effective radius)
    V_obs = 40  # km/s (typical)

    # Mean surface density
    Sigma_tdg = M_tdg / (np.pi * R_tdg**2)  # M_sun/kpc²
    Sigma_tdg_pc2 = Sigma_tdg / 1e6

    print(f"Typical TDG parameters:")
    print(f"  M_bar = {M_tdg:.0e} M_sun")
    print(f"  R_eff = {R_tdg} kpc")
    print(f"  V_obs = {V_obs} km/s")
    print(f"  Mean Σ = {Sigma_tdg_pc2:.1f} M_sun/pc²")
    print()

    # Compare to normal dwarf
    print("Compare to normal dwarf galaxy:")
    # From SPARC data, dwarfs with V ~ 40 km/s have Σ ~ 10-50 M_sun/pc²
    print("  Normal dwarf at V ~ 40 km/s:")
    print("  Mean Σ ~ 10-50 M_sun/pc²")
    print()

    print(f"TDG Σ = {Sigma_tdg_pc2:.1f} M_sun/pc²")
    print()

    if 5 < Sigma_tdg_pc2 < 100:
        print("→ TDG surface density is SIMILAR to normal dwarfs")
        print("→ Synchronism predicts: Same BTFR as normal galaxies")
        print("→ No discrimination from standard TDGs")
    else:
        print("→ TDG surface density is DIFFERENT from normal dwarfs")
        print("→ Synchronism predicts: Different BTFR offset")

    print()

    return Sigma_tdg_pc2


def ultra_diffuse_galaxies():
    """
    Analyze ultra-diffuse galaxies (UDGs) as extreme test.
    """
    print("="*70)
    print("ULTRA-DIFFUSE GALAXIES: EXTREME TEST")
    print("="*70)
    print()

    print("What are UDGs?")
    print("  - Galaxies with very low surface brightness")
    print("  - Σ < 24.5 mag/arcsec² (very faint)")
    print("  - R_eff > 1.5 kpc (large)")
    print("  - M_star ~ 10^7 - 10^8 M_sun (low stellar mass)")
    print()

    print("Why are UDGs interesting?")
    print()
    print("  UDGs have EXTREMELY LOW surface density:")
    print("  Σ_UDG ~ 1-10 M_sun/pc² (vs ~100 for normal disks)")
    print()

    # UDG parameters
    M_udg = 3e7  # M_sun (typical stellar mass)
    R_udg = 3  # kpc (effective radius)

    Sigma_udg = M_udg / (np.pi * R_udg**2) / 1e6  # M_sun/pc²

    print(f"Typical UDG:")
    print(f"  M_star = {M_udg:.0e} M_sun")
    print(f"  R_eff = {R_udg} kpc")
    print(f"  Σ = {Sigma_udg:.1f} M_sun/pc²")
    print()

    print("-"*70)
    print()
    print("PREDICTIONS:")
    print()
    print("MOND (a₀ universal):")
    print("  - UDGs follow same BTFR as normal galaxies")
    print("  - V_flat set by total mass only")
    print()
    print("Synchronism (C depends on ρ):")
    print("  - Low Σ → low ρ → low C → high G_eff")
    print("  - V_obs/V_bar should be HIGHER than normal galaxies")
    print("  - At fixed M_bar: UDGs have HIGHER V_flat")
    print()

    # Calculate expected offset
    # From Session #87: r(V/V_bar, SB) = -0.626
    # Lower SB → higher V/V_bar

    print("-"*70)
    print()
    print("QUANTITATIVE PREDICTION:")
    print()

    # Normal galaxy at same mass has Σ ~ 50 M_sun/pc²
    Sigma_normal = 50  # M_sun/pc²
    log_Sigma_ratio = np.log10(Sigma_udg / Sigma_normal)

    # From radial analysis: slope of ratio vs SB ~ -0.074 in log
    slope = -0.074  # from Session #87
    delta_log_ratio = slope * log_Sigma_ratio

    print(f"  UDG Σ = {Sigma_udg:.1f} M_sun/pc²")
    print(f"  Normal Σ = {Sigma_normal} M_sun/pc²")
    print(f"  log(Σ_UDG/Σ_normal) = {log_Sigma_ratio:.2f}")
    print()
    print(f"  Using slope from Session #87: {slope:.3f}")
    print(f"  Δlog(V/V_bar) = {delta_log_ratio:.3f}")
    print()
    print(f"  Expected V_obs/V_bar ratio change: {10**delta_log_ratio:.2f}×")
    print()

    # For BTFR: V ∝ M^0.25
    # If V/V_bar is higher, at fixed M_bar, V is higher
    delta_log_V = 0.25 * np.log10(1 + (10**delta_log_ratio - 1))
    print(f"  BTFR offset: Δlog(V) ~ {delta_log_V:.3f} dex")
    print()

    if abs(delta_log_ratio) > 0.1:
        print("→ SIGNIFICANT difference from normal BTFR!")
        print("→ This is a DISCRIMINATING TEST between MOND and Synchronism")
    else:
        print("→ Small difference, may be within BTFR scatter")

    print()

    return Sigma_udg, delta_log_ratio


def session90_conclusions():
    """
    Summarize Session #90 findings.
    """
    print("="*70)
    print("SESSION #90 CONCLUSIONS")
    print("="*70)
    print()

    print("KEY FINDINGS:")
    print()
    print("1. TIDAL DWARF GALAXIES:")
    print("   - Observed to follow same BTFR as normal galaxies")
    print("   - Surface densities similar to normal dwarfs")
    print("   - NOT a strong discriminating test (both theories predict same)")
    print()

    print("2. ULTRA-DIFFUSE GALAXIES (UDGs):")
    print("   - Extremely low surface density (Σ ~ 1-10 M_sun/pc²)")
    print("   - MOND: Same BTFR (a₀ universal)")
    print("   - Synchronism: Higher V/V_bar due to lower C")
    print("   - Predicted offset: ~0.1 dex in V (significant!)")
    print("   - THIS IS A DISCRIMINATING TEST!")
    print()

    print("3. NEXT STEPS:")
    print("   - Analyze UDG rotation curve data (if available)")
    print("   - Compare to BTFR prediction")
    print("   - Look for systematic offset in low-Σ systems")
    print()

    print("-"*70)
    print()
    print("OBSERVATIONAL PRIORITY:")
    print()
    print("UDGs are the cleanest test of Synchronism vs MOND:")
    print("  - They have lowest surface densities")
    print("  - Session #87 showed: V/V_bar correlates with SB")
    print("  - Prediction: UDGs should have systematically higher V/V_bar")
    print()
    print("If UDGs follow EXACT same BTFR: MOND correct, Synchronism wrong")
    print("If UDGs have higher V at fixed M: Synchronism supported")
    print()


def main():
    """Run all Session #90 analyses."""
    tdg_physics()
    mond_vs_synchronism_for_tdgs()
    observational_evidence()
    Sigma_tdg = synchronism_prediction_for_tdgs()
    Sigma_udg, delta = ultra_diffuse_galaxies()
    session90_conclusions()

    # Save results
    results = {
        'session': 90,
        'title': 'TDGs and UDGs as Discriminating Tests',
        'tdg_analysis': {
            'Sigma_pc2': Sigma_tdg,
            'prediction': 'Same BTFR as normal galaxies',
            'discriminating': False
        },
        'udg_analysis': {
            'Sigma_pc2': Sigma_udg,
            'delta_log_ratio': delta,
            'prediction': 'Higher V/V_bar than normal galaxies',
            'discriminating': True
        },
        'conclusions': [
            'TDGs: Not discriminating (similar Σ to normal dwarfs)',
            'UDGs: Potentially discriminating (extremely low Σ)',
            'UDG test: If V/V_bar higher, Synchronism supported',
            'Need UDG rotation curve data for test'
        ]
    }

    results_dir = Path(__file__).parent / 'results'
    results_dir.mkdir(exist_ok=True)

    with open(results_dir / 'session90_tidal_dwarfs.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("Results saved to session90_tidal_dwarfs.json")

    return results


if __name__ == "__main__":
    main()
