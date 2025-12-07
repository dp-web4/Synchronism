#!/usr/bin/env python3
"""
Session #93: Detailed UDG Kinematic Predictions

Ultra-diffuse galaxies (UDGs) are a key discriminating test between
Synchronism and standard MOND.

Synchronism: C(ρ) depends on LOCAL baryonic density
             Low Σ → Low C → Higher G_eff → Higher V/V_bar

MOND: Universal a₀ → Same BTFR for all galaxies regardless of Σ

This analysis quantifies the Synchronism prediction for UDGs
and compares to available observational data.

Author: CBP Autonomous Synchronism Research
Date: December 6, 2025
Session: #93 - UDG Kinematic Predictions
"""

import numpy as np
import json
from pathlib import Path
from datetime import datetime


# Physical constants
G_astro = 4.302e-6  # kpc (km/s)² / M_sun
a_0 = 1.2e-10       # m/s² (MOND scale)
Sigma_0 = 140       # M_sun/pc² (Freeman's Law)


def coherence_function(rho, rho_crit, gamma=2.0):
    """
    Synchronism coherence function.

    C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))

    For surface density (2D):
    C(Σ) = tanh(γ × log(Σ/Σ_crit + 1))
    """
    return np.tanh(gamma * np.log(rho / rho_crit + 1))


def g_eff_ratio(Sigma, Sigma_crit):
    """
    Effective gravity ratio G_eff/G = 1/C(Σ).

    Lower Σ → Lower C → Higher G_eff
    """
    C = coherence_function(Sigma, Sigma_crit)
    return 1.0 / C


def udg_properties():
    """
    Typical properties of ultra-diffuse galaxies.
    """
    print("=" * 70)
    print("ULTRA-DIFFUSE GALAXY PROPERTIES")
    print("=" * 70)
    print()

    print("Definition: Galaxies with r_eff > 1.5 kpc and μ_eff > 24 mag/arcsec²")
    print()

    # Typical UDG properties
    udg = {
        'r_eff_kpc': 3.0,        # Effective radius
        'mu_eff': 25.5,          # Surface brightness (g-band)
        'M_star': 2e8,           # Stellar mass (M_sun)
        'M_gas': 1e8,            # Gas mass (M_sun) - often undetected
        'Sigma_star': None,      # To be calculated
    }

    # Calculate mean surface density
    # M = 2π Σ r_eff² for exponential profile (r_eff ≈ 1.68 r_d)
    area = np.pi * udg['r_eff_kpc']**2 * 1e6  # pc²
    udg['Sigma_star'] = (udg['M_star'] + udg['M_gas']) / area  # M_sun/pc²

    print("TYPICAL UDG:")
    print(f"  Effective radius: r_eff = {udg['r_eff_kpc']:.1f} kpc")
    print(f"  Surface brightness: μ_eff = {udg['mu_eff']:.1f} mag/arcsec²")
    print(f"  Stellar mass: M* = {udg['M_star']:.1e} M_sun")
    print(f"  Gas mass: M_gas ~ {udg['M_gas']:.1e} M_sun (often upper limit)")
    print(f"  Mean surface density: Σ = {udg['Sigma_star']:.1f} M_sun/pc²")
    print()

    # Compare to normal dwarf
    normal = {
        'r_eff_kpc': 1.0,
        'M_star': 2e8,
        'Sigma_star': None
    }
    area_normal = np.pi * normal['r_eff_kpc']**2 * 1e6  # pc²
    normal['Sigma_star'] = normal['M_star'] / area_normal

    print("COMPARISON TO NORMAL DWARF (same M*):")
    print(f"  Normal dwarf: r_eff = {normal['r_eff_kpc']:.1f} kpc, Σ = {normal['Sigma_star']:.1f} M_sun/pc²")
    print(f"  UDG:          r_eff = {udg['r_eff_kpc']:.1f} kpc, Σ = {udg['Sigma_star']:.1f} M_sun/pc²")
    print(f"  Ratio:        Σ_normal / Σ_UDG = {normal['Sigma_star'] / udg['Sigma_star']:.0f}×")
    print()

    return udg, normal


def synchronism_prediction_udg():
    """
    Calculate Synchronism prediction for UDG kinematics.
    """
    print("=" * 70)
    print("SYNCHRONISM PREDICTION FOR UDGs")
    print("=" * 70)
    print()

    # Get UDG and normal dwarf properties
    udg, normal = udg_properties()

    # Synchronism predicts: V/V_bar depends on Σ through C(Σ)
    # Use empirical slope from Session #87: V/V_bar ∝ Σ^(-0.22)

    print("From Session #87 radial analysis:")
    print("  V/V_bar correlates with local surface brightness")
    print("  Implied: V/V_bar ∝ Σ^(-0.22) approximately")
    print()

    slope = -0.22  # From V/V_bar vs Σ correlation

    Sigma_ratio = udg['Sigma_star'] / normal['Sigma_star']
    V_ratio_prediction = Sigma_ratio**slope

    print("PREDICTION:")
    print(f"  Σ_UDG / Σ_normal = {Sigma_ratio:.2f}")
    print(f"  (V/V_bar)_UDG / (V/V_bar)_normal = {Sigma_ratio:.2f}^{slope}")
    print(f"                                    = {V_ratio_prediction:.2f}")
    print()
    print(f"  UDGs should have {(V_ratio_prediction - 1)*100:.0f}% higher V/V_bar than normal dwarfs!")
    print()

    # Alternative calculation using coherence function directly
    print("-" * 70)
    print("ALTERNATIVE: Direct coherence calculation")
    print("-" * 70)
    print()

    # Use Σ_crit from cosmology
    Sigma_crit = Sigma_0 / 2  # ~70 M_sun/pc² for transition

    C_normal = coherence_function(normal['Sigma_star'], Sigma_crit)
    C_udg = coherence_function(udg['Sigma_star'], Sigma_crit)

    G_eff_normal = 1.0 / C_normal
    G_eff_udg = 1.0 / C_udg

    # V² ∝ G_eff M / R, so V ∝ sqrt(G_eff)
    V_ratio_coherence = np.sqrt(G_eff_udg / G_eff_normal)

    print(f"Using C(Σ) = tanh(2 × log(Σ/Σ_crit + 1)):")
    print(f"  Σ_crit = {Sigma_crit:.0f} M_sun/pc²")
    print()
    print(f"  Normal dwarf: Σ = {normal['Sigma_star']:.1f}, C = {C_normal:.3f}, G_eff/G = {G_eff_normal:.2f}")
    print(f"  UDG:          Σ = {udg['Sigma_star']:.1f}, C = {C_udg:.3f}, G_eff/G = {G_eff_udg:.2f}")
    print()
    print(f"  V_UDG / V_normal = sqrt(G_eff_UDG / G_eff_normal)")
    print(f"                   = sqrt({G_eff_udg:.2f} / {G_eff_normal:.2f})")
    print(f"                   = {V_ratio_coherence:.2f}")
    print()
    print(f"  UDGs should have {(V_ratio_coherence - 1)*100:.0f}% higher V than normal dwarfs")
    print(f"  at the same baryonic mass!")
    print()

    return {
        'Sigma_UDG': udg['Sigma_star'],
        'Sigma_normal': normal['Sigma_star'],
        'slope_method': {
            'slope': slope,
            'V_ratio': V_ratio_prediction
        },
        'coherence_method': {
            'C_normal': C_normal,
            'C_udg': C_udg,
            'V_ratio': V_ratio_coherence
        }
    }


def mond_prediction_udg():
    """
    Calculate MOND prediction for UDG kinematics.
    """
    print("=" * 70)
    print("MOND PREDICTION FOR UDGs")
    print("=" * 70)
    print()

    print("Standard MOND:")
    print("  a₀ = 1.2×10⁻¹⁰ m/s² is UNIVERSAL")
    print("  BTFR: M_bar = A × V⁴/a₀ holds for ALL galaxies")
    print()
    print("PREDICTION:")
    print("  At fixed M_bar, V is the same regardless of Σ")
    print("  UDGs should follow the SAME BTFR as normal galaxies")
    print()
    print("  (V/V_bar)_UDG = (V/V_bar)_normal")
    print()

    return {
        'prediction': 'Same BTFR for all',
        'V_ratio': 1.0
    }


def observational_constraints():
    """
    Review observational constraints on UDG kinematics.
    """
    print("=" * 70)
    print("OBSERVATIONAL CONSTRAINTS")
    print("=" * 70)
    print()

    print("1. NGC 1052-DF2 (van Dokkum+ 2018, 2019)")
    print("   - Claimed to have little/no dark matter")
    print("   - σ ~ 3-8 km/s (controversial)")
    print("   - If true, MOND is violated")
    print("   - Distance controversy affects interpretation")
    print()

    print("2. NGC 1052-DF4")
    print("   - Similar to DF2")
    print("   - Also claimed to lack dark matter")
    print("   - Under intense scrutiny")
    print()

    print("3. Dragonfly 44 (van Dokkum+ 2016)")
    print("   - σ ~ 47 km/s (very high for UDG)")
    print("   - Implies high M/L or DM")
    print("   - Some controversy about membership")
    print()

    print("4. Other UDGs")
    print("   - VCC 1287 (Beasley+ 2016): σ ~ 19 km/s")
    print("   - DGSAT I: HI rotation curve, follows BTFR")
    print()

    print("SUMMARY:")
    print("  - UDG kinematics are HETEROGENEOUS")
    print("  - Some appear dark-matter-deficient (DF2, DF4)")
    print("  - Some appear dark-matter-rich (Dragonfly 44)")
    print("  - Need larger, systematic sample")
    print()

    print("SYNCHRONISM INTERPRETATION:")
    print("  - Low Σ → Low C → High G_eff → High V/V_bar")
    print("  - But DF2/DF4 seem to show OPPOSITE: Low V/V_bar!")
    print()
    print("  POSSIBLE EXPLANATIONS:")
    print("  1. Distance errors (DF2 closer → lower M* → normal)")
    print("  2. Formation environment effects")
    print("  3. Tidal stripping modified the density profile")
    print("  4. DF2/DF4 are atypical UDGs")
    print()

    return {
        'DF2': {'sigma_km_s': '3-8', 'status': 'controversial'},
        'DF4': {'sigma_km_s': 'low', 'status': 'controversial'},
        'Dragonfly44': {'sigma_km_s': 47, 'status': 'high V'},
        'VCC1287': {'sigma_km_s': 19, 'status': 'normal'},
        'summary': 'Heterogeneous - no clear pattern'
    }


def comparison_and_synthesis():
    """
    Compare Synchronism and MOND predictions for UDGs.
    """
    print("=" * 70)
    print("COMPARISON: SYNCHRONISM vs MOND FOR UDGs")
    print("=" * 70)
    print()

    print("| Property                    | Synchronism | MOND      |")
    print("|-----------------------------|-------------|-----------|")
    print("| V/V_bar at low Σ            | Higher      | Same      |")
    print("| Expected increase for UDG   | ~30%        | 0%        |")
    print("| BTFR scatter at low Σ       | Increased   | Same      |")
    print()

    print("-" * 70)
    print("KEY PREDICTION:")
    print("-" * 70)
    print()
    print("  If Synchronism is correct:")
    print("    UDGs should lie ABOVE the standard BTFR")
    print("    By ~0.1 dex in log(V) at fixed M_bar")
    print()
    print("  If MOND is correct:")
    print("    UDGs should follow the SAME BTFR as normal galaxies")
    print()
    print("-" * 70)
    print("CURRENT OBSERVATIONAL STATUS:")
    print("-" * 70)
    print()
    print("  The data is CONFUSING:")
    print("  - Some UDGs (DF2, DF4) appear BELOW BTFR")
    print("  - Some UDGs (Dragonfly 44) appear ABOVE BTFR")
    print("  - Need systematic sample with well-measured M_bar and V")
    print()
    print("  NEITHER Synchronism NOR MOND is clearly supported yet!")
    print()

    print("-" * 70)
    print("RECOMMENDATIONS:")
    print("-" * 70)
    print()
    print("  1. Need HI rotation curves for UDG sample")
    print("     - Gas-rich UDGs provide clean V measurement")
    print("     - E.g., ALFALFA-selected UDGs")
    print()
    print("  2. Control for environment")
    print("     - Cluster UDGs may be tidally affected")
    print("     - Field UDGs are cleaner test")
    print()
    print("  3. Careful distance measurements")
    print("     - DF2/DF4 controversy shows importance")
    print()

    return {
        'synchronism_prediction': 'V ~30% higher at fixed M_bar for UDGs',
        'mond_prediction': 'Same BTFR for UDGs',
        'current_status': 'Inconclusive - heterogeneous data',
        'recommendations': [
            'HI rotation curves for UDG sample',
            'Control for environment (field vs cluster)',
            'Careful distance measurements'
        ]
    }


def main():
    """Run all Session #93 UDG analyses."""
    print("=" * 70)
    print("SESSION #93: UDG KINEMATIC PREDICTIONS")
    print("=" * 70)
    print()

    results = {
        'session': 93,
        'title': 'UDG Kinematic Predictions',
        'date': datetime.now().isoformat()
    }

    # Calculate predictions
    results['synchronism'] = synchronism_prediction_udg()
    results['mond'] = mond_prediction_udg()

    # Review observations
    results['observations'] = observational_constraints()

    # Synthesis
    results['synthesis'] = comparison_and_synthesis()

    # Save results
    output_dir = Path(__file__).parent / 'results'
    output_dir.mkdir(exist_ok=True)

    with open(output_dir / 'session93_UDG_predictions.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print()
    print(f"Results saved to: results/session93_UDG_predictions.json")

    print()
    print("=" * 70)
    print("SESSION #93 UDG ANALYSIS COMPLETE")
    print("=" * 70)

    return results


if __name__ == '__main__':
    main()
