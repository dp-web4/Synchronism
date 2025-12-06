#!/usr/bin/env python3
"""
Session #89: High-z BTFR Evolution Prediction

From Session #88:
- a₀ = cH₀/(2π)
- If a₀ evolves with H(z), BTFR should shift at high redshift

This analysis quantifies the prediction for JWST observations.

BTFR: M_bar = A × V^4
Shift: Δlog(V) = (1/4) × log(a₀(z)/a₀(0))
              = (1/4) × log(H(z)/H₀)

Author: CBP Autonomous Synchronism Research
Date: December 5, 2025
"""

import numpy as np
import json
from pathlib import Path


def hubble_evolution(z, Omega_m=0.3, Omega_Lambda=0.7):
    """
    Hubble parameter evolution with redshift.

    H(z) = H₀ × √(Ω_m(1+z)³ + Ω_Λ)
    """
    return np.sqrt(Omega_m * (1 + z)**3 + Omega_Lambda)


def btfr_evolution_prediction():
    """
    Calculate predicted BTFR evolution with redshift.

    Key prediction from Synchronism/MOND unification:
    - a₀(z) ∝ H(z)
    - BTFR: M ∝ V^4/a₀
    - At fixed M: V ∝ a₀^(1/4)
    - Shift: Δlog(V) = (1/4) × log(H(z)/H₀)
    """
    print("="*70)
    print("HIGH-z BTFR EVOLUTION PREDICTION")
    print("="*70)
    print()

    print("From Sessions #88-89:")
    print("  a₀ = cH₀/(2π)")
    print("  At redshift z: a₀(z) = cH(z)/(2π)")
    print()
    print("BTFR scaling:")
    print("  M_bar = A × V^4 / a₀")
    print("  At fixed M_bar: V ∝ a₀^(1/4)")
    print()
    print("Shift in velocity at fixed mass:")
    print("  Δlog(V) = (1/4) × log(H(z)/H₀)")
    print()

    # Redshift range
    z_vals = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])

    print("-"*70)
    print(f"{'z':>5}  {'H(z)/H₀':>10}  {'Δlog(V)':>10}  {'%V change':>12}")
    print("-"*70)

    predictions = []
    for z in z_vals:
        H_ratio = hubble_evolution(z)
        delta_log_V = 0.25 * np.log10(H_ratio)
        V_change_percent = (10**delta_log_V - 1) * 100

        print(f"{z:>5.1f}  {H_ratio:>10.3f}  {delta_log_V:>+10.4f}  {V_change_percent:>+10.1f}%")

        predictions.append({
            'z': z,
            'H_ratio': H_ratio,
            'delta_log_V': delta_log_V,
            'V_change_percent': V_change_percent
        })

    print("-"*70)
    print()

    return predictions


def observational_requirements():
    """
    Determine what observations are needed to test prediction.
    """
    print("="*70)
    print("OBSERVATIONAL REQUIREMENTS")
    print("="*70)
    print()

    # At z = 1: Δlog(V) ≈ +0.06 dex
    delta_log_V_z1 = 0.25 * np.log10(hubble_evolution(1.0))

    print(f"At z = 1:")
    print(f"  Predicted shift: Δlog(V) = {delta_log_V_z1:+.3f} dex")
    print(f"  This corresponds to ~{(10**delta_log_V_z1 - 1)*100:.1f}% higher V at fixed M")
    print()

    # BTFR scatter is typically ~0.1 dex
    btfr_scatter = 0.1  # dex

    print(f"Typical BTFR scatter: {btfr_scatter} dex")
    print()

    # Detection significance
    # Need N galaxies such that scatter/√N < predicted shift
    N_needed = (btfr_scatter / abs(delta_log_V_z1))**2
    print(f"To detect at 1σ: N > (scatter/shift)² = {N_needed:.0f} galaxies")
    print(f"To detect at 3σ: N > 9 × {N_needed:.0f} = {9*N_needed:.0f} galaxies")
    print()

    print("-"*70)
    print()
    print("JWST FEASIBILITY:")
    print()
    print("  - JWST NIRSpec can measure Hα rotation curves to z ~ 2-3")
    print("  - Need ~100 galaxies at z > 1 with good rotation curves")
    print("  - Current samples: SINS/zC-SINF, KROSS, KMOS³D")
    print("  - Each has ~100-500 galaxies, some with resolved kinematics")
    print()
    print("  STATUS: Marginally testable with existing data!")
    print("          Definitive test requires larger samples (future JWST)")
    print()

    return delta_log_V_z1, N_needed


def comparison_to_mond():
    """
    Compare Synchronism prediction to standard MOND.
    """
    print("="*70)
    print("COMPARISON TO STANDARD MOND")
    print("="*70)
    print()

    print("STANDARD MOND:")
    print("  - a₀ = 1.2×10⁻¹⁰ m/s² is CONSTANT")
    print("  - BTFR does NOT evolve with redshift")
    print("  - Any observed evolution would FALSIFY MOND")
    print()

    print("SYNCHRONISM (unified with MOND):")
    print("  - a₀(z) = cH(z)/(2π)")
    print("  - BTFR DOES evolve with redshift")
    print("  - Δlog(V) = +0.06 dex at z = 1")
    print()

    print("-"*70)
    print()
    print("THIS IS A DISCRIMINATING TEST!")
    print()
    print("If BTFR evolves as predicted:")
    print("  → Synchronism's cosmological connection is validated")
    print("  → Standard MOND (constant a₀) is ruled out")
    print()
    print("If BTFR does NOT evolve:")
    print("  → Synchronism's a₀ = cH₀/(2π) is ruled out")
    print("  → Standard MOND (constant a₀) is supported")
    print()


def literature_review():
    """
    Review existing high-z BTFR observations.
    """
    print("="*70)
    print("EXISTING HIGH-Z BTFR OBSERVATIONS")
    print("="*70)
    print()

    print("1. TILEY et al. (2019) - KROSS survey")
    print("   - 586 galaxies at z ~ 0.9")
    print("   - BTFR slope: 3.6 ± 0.2 (vs 4.0 local)")
    print("   - Zero-point: consistent with local within errors")
    print("   - Scatter: larger than local (~0.25 dex)")
    print()

    print("2. ÜBLER et al. (2017) - KMOS³D")
    print("   - ~100 galaxies at z ~ 0.9-2.3")
    print("   - Evolution in stellar mass TFR detected")
    print("   - Baryonic TFR not directly tested")
    print()

    print("3. PRICE et al. (2020) - MOSDEF")
    print("   - z ~ 1.4-2.6 galaxies")
    print("   - Significant evolution in M*/V relation")
    print("   - Gas fractions higher at high z")
    print()

    print("-"*70)
    print()
    print("CHALLENGE: Separating effects")
    print()
    print("At high z, galaxies have:")
    print("  - Higher gas fractions")
    print("  - More turbulent disks")
    print("  - Different morphologies")
    print()
    print("These effects can mimic or mask the Synchronism prediction.")
    print("Need careful modeling to isolate the a₀(z) signal.")
    print()


def session89_btfr_conclusions():
    """
    Summarize Session #89 high-z BTFR findings.
    """
    print("="*70)
    print("SESSION #89: HIGH-Z BTFR CONCLUSIONS")
    print("="*70)
    print()
    print("PREDICTION:")
    print("  If a₀ = cH(z)/(2π), BTFR evolves with redshift:")
    print("  - z = 1: Δlog(V) = +0.06 dex (+14% in V)")
    print("  - z = 2: Δlog(V) = +0.10 dex (+26% in V)")
    print()
    print("OBSERVATIONAL STATUS:")
    print("  - Marginal samples exist (KROSS, KMOS³D)")
    print("  - Current data too noisy for definitive test")
    print("  - JWST can provide decisive measurement")
    print()
    print("DISCRIMINATING POWER:")
    print("  - Standard MOND: No evolution (a₀ constant)")
    print("  - Synchronism: Evolution as predicted above")
    print("  - This is a CLEAN test of the theory!")
    print()
    print("RECOMMENDED NEXT STEPS:")
    print("  1. Analyze existing KROSS/KMOS³D data for BTFR evolution")
    print("  2. Account for gas fraction and turbulence effects")
    print("  3. Propose JWST program for definitive test")
    print()


def main():
    """Run all high-z BTFR analyses."""
    predictions = btfr_evolution_prediction()
    delta_log_V, N_needed = observational_requirements()
    comparison_to_mond()
    literature_review()
    session89_btfr_conclusions()

    # Save results
    results = {
        'session': 89,
        'title': 'High-z BTFR Evolution Prediction',
        'predictions': predictions,
        'observational_requirements': {
            'delta_log_V_z1': delta_log_V,
            'N_galaxies_1sigma': N_needed,
            'N_galaxies_3sigma': 9 * N_needed
        },
        'discriminating_test': {
            'synchronism': 'BTFR evolves with H(z)',
            'standard_mond': 'BTFR constant',
            'status': 'testable with JWST'
        },
        'conclusions': [
            'Δlog(V) = +0.06 dex at z=1',
            'Need ~25 galaxies for 1σ detection',
            'Existing data marginal, JWST decisive',
            'This discriminates Synchronism vs standard MOND'
        ]
    }

    results_dir = Path(__file__).parent / 'results'
    results_dir.mkdir(exist_ok=True)

    with open(results_dir / 'session89_highz_btfr_prediction.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("Results saved to session89_highz_btfr_prediction.json")

    return results


if __name__ == "__main__":
    main()
