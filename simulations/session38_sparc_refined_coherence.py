#!/usr/bin/env python3
"""
Session #38: SPARC Dark Matter - Refined Coherence Formula

Test refined coherence formula that avoids saturation in high-density regions.

Current formula (Session #17):
    C_vis = (ρ_vis/ρ_0)^γ     (γ = 0.30)
    Problem: C → 1 as ρ_vis → ∞ (saturates)

Refined formula (this session):
    C_vis = 1 - exp(-(ρ_vis/ρ_crit)^γ)
    Advantage: Never fully saturates, always allows some DM

Dark matter formula (unchanged):
    ρ_DM = α(1 - C_vis) × ρ_vis^β     (β = 0.30)

Goal: Test if refined formula improves massive spiral fits while maintaining
      irregular galaxy success.

Author: CBP Autonomous Synchronism Research
Date: 2025-11-22
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json

# Import SPARC parsing utilities
sys.path.append(str(Path(__file__).parent))
from synchronism_real_sparc_validation import (
    load_sparc_galaxy,
    calculate_visible_density,
    fit_alpha_parameter
)

def coherence_original(rho_vis, rho_0=1.0, gamma=0.30):
    """
    Original coherence formula from Session #14-17.

    Problem: Saturates to 1 for high density.
    """
    return (rho_vis / rho_0) ** gamma


def coherence_refined(rho_vis, rho_crit=1.0, gamma=0.30):
    """
    Refined coherence formula that avoids full saturation.

    C_vis = 1 - exp(-(ρ_vis/ρ_crit)^γ)

    Properties:
    - C → 0 as ρ_vis → 0 (no coherence in vacuum)
    - C → 1 as ρ_vis → ∞ but NEVER reaches 1 exactly
    - Always leaves room for dark matter: (1 - C_vis) > 0

    Parameters:
    -----------
    rho_vis : float or ndarray
        Visible matter density (M_sun / pc^3)
    rho_crit : float
        Critical density scale (fit parameter)
    gamma : float
        Power law exponent (0.30 from theory)

    Returns:
    --------
    C_vis : float or ndarray
        Visible matter coherence [0, 1)
    """
    return 1.0 - np.exp(-(rho_vis / rho_crit) ** gamma)


def dark_matter_density(rho_vis, alpha, beta=0.30, coherence_func=coherence_original, **coherence_params):
    """
    Dark matter density from Synchronism coherence.

    ρ_DM = α(1 - C_vis) × ρ_vis^β

    Parameters:
    -----------
    rho_vis : ndarray
        Visible matter density profile
    alpha : float
        Normalization parameter (fit per galaxy)
    beta : float
        Scaling exponent (0.30 from theory)
    coherence_func : callable
        Coherence function (original or refined)
    coherence_params : dict
        Parameters for coherence function

    Returns:
    --------
    rho_DM : ndarray
        Dark matter density profile
    """
    C_vis = coherence_func(rho_vis, **coherence_params)
    return alpha * (1.0 - C_vis) * (rho_vis ** beta)


def rotation_curve_from_density(r, rho_vis, rho_DM):
    """
    Calculate rotation velocity from density profiles.

    v_rot^2 = G M(<r) / r

    where M(<r) = ∫[0 to r] 4πr'^2 (ρ_vis + ρ_DM) dr'
    """
    G = 4.302e-6  # Gravitational constant (kpc (km/s)^2 / M_sun)

    # Numerical integration (trapezoidal rule)
    r_kpc = r / 1000.0  # Convert pc to kpc

    # Calculate enclosed mass at each radius
    M_enc = np.zeros_like(r)
    for i in range(1, len(r)):
        # Integrate from 0 to r[i]
        r_int = r[:i+1] / 1000.0  # kpc
        rho_tot = (rho_vis[:i+1] + rho_DM[:i+1]) * 1e9  # Convert to M_sun/kpc^3
        integrand = 4 * np.pi * r_int**2 * rho_tot
        M_enc[i] = np.trapz(integrand, r_int)

    # Avoid division by zero
    M_enc[0] = M_enc[1] if len(M_enc) > 1 else 0

    # Calculate rotation velocity
    v_rot = np.sqrt(G * M_enc / np.maximum(r_kpc, 1e-10))

    return v_rot


def fit_galaxy_refined(galaxy_data, rho_crit_values=np.logspace(-2, 2, 50)):
    """
    Fit refined coherence model to galaxy rotation curve.

    Free parameters:
    - alpha: DM normalization (fit per galaxy)
    - rho_crit: Critical density scale for coherence

    Fixed parameters:
    - gamma = 0.30 (theory prediction)
    - beta = 0.30 (theory prediction)

    Returns:
    --------
    best_params : dict
        Best-fit parameters and statistics
    """
    r = galaxy_data['r']  # pc
    v_obs = galaxy_data['v']  # km/s
    dv = galaxy_data['dv']  # km/s

    # Calculate visible matter density
    rho_vis = calculate_visible_density(galaxy_data)

    best_chi2 = np.inf
    best_params = None

    # Grid search over rho_crit
    for rho_crit in rho_crit_values:
        # For each rho_crit, fit alpha
        C_vis = coherence_refined(rho_vis, rho_crit=rho_crit, gamma=0.30)

        # Fit alpha (same method as Session #17)
        result = fit_alpha_parameter(galaxy_data, rho_vis, C_vis, beta=0.30)

        if result['chi2'] < best_chi2:
            best_chi2 = result['chi2']
            best_params = {
                'alpha': result['alpha'],
                'rho_crit': rho_crit,
                'chi2': result['chi2'],
                'chi2_red': result['chi2_red'],
                'M_DM': result['M_DM'],
                'M_vis': result['M_vis'],
                'M_DM_M_vis': result['M_DM'] / result['M_vis'] if result['M_vis'] > 0 else 0
            }

    return best_params


def compare_models(galaxy_name, galaxy_data):
    """
    Compare original vs refined coherence for single galaxy.
    """
    r = galaxy_data['r']
    v_obs = galaxy_data['v']
    dv = galaxy_data['dv']

    rho_vis = calculate_visible_density(galaxy_data)

    # Original model
    C_orig = coherence_original(rho_vis, rho_0=1.0, gamma=0.30)
    result_orig = fit_alpha_parameter(galaxy_data, rho_vis, C_orig, beta=0.30)

    # Refined model
    result_refined = fit_galaxy_refined(galaxy_data)

    return {
        'galaxy': galaxy_name,
        'original': result_orig,
        'refined': result_refined,
        'improvement': result_orig['chi2_red'] - result_refined['chi2_red']
    }


def run_full_sparc_comparison(data_dir='sparc_real_data', max_galaxies=None):
    """
    Run comparison on all available SPARC galaxies.
    """
    data_path = Path(data_dir)
    if not data_path.exists():
        print(f"Error: {data_dir} not found")
        return None

    galaxy_files = sorted(data_path.glob('*.dat'))
    if max_galaxies:
        galaxy_files = galaxy_files[:max_galaxies]

    results = []

    print("\n" + "="*80)
    print("SESSION #38: REFINED COHERENCE FORMULA TEST")
    print("="*80)
    print(f"\nTesting on {len(galaxy_files)} SPARC galaxies")
    print("Original: C_vis = (ρ_vis/ρ_0)^0.30")
    print("Refined:  C_vis = 1 - exp(-(ρ_vis/ρ_crit)^0.30)")
    print("\n" + "-"*80)

    for i, galaxy_file in enumerate(galaxy_files, 1):
        galaxy_name = galaxy_file.stem

        try:
            galaxy_data = load_sparc_galaxy(str(galaxy_file))
            comparison = compare_models(galaxy_name, galaxy_data)
            results.append(comparison)

            if i % 10 == 0 or i == len(galaxy_files):
                print(f"Processed {i}/{len(galaxy_files)} galaxies...")

        except Exception as e:
            print(f"  Error processing {galaxy_name}: {e}")
            continue

    return results


def analyze_results(results):
    """
    Statistical analysis of original vs refined comparison.
    """
    if not results:
        print("No results to analyze")
        return

    print("\n" + "="*80)
    print("RESULTS ANALYSIS")
    print("="*80)

    # Extract metrics
    chi2_orig = np.array([r['original']['chi2_red'] for r in results])
    chi2_refined = np.array([r['refined']['chi2_red'] for r in results])
    improvements = chi2_orig - chi2_refined

    # Overall statistics
    print("\n### Overall Performance ###")
    print(f"Galaxies analyzed: {len(results)}")
    print(f"\nOriginal model:")
    print(f"  Median χ²_red: {np.median(chi2_orig):.2f}")
    print(f"  Mean χ²_red: {np.mean(chi2_orig):.2f} ± {np.std(chi2_orig):.2f}")
    print(f"\nRefined model:")
    print(f"  Median χ²_red: {np.median(chi2_refined):.2f}")
    print(f"  Mean χ²_red: {np.mean(chi2_refined):.2f} ± {np.std(chi2_refined):.2f}")

    # Improvement statistics
    print(f"\n### Improvement Statistics ###")
    print(f"Median improvement: {np.median(improvements):.3f}")
    print(f"Mean improvement: {np.mean(improvements):.3f} ± {np.std(improvements):.3f}")
    print(f"Galaxies improved: {np.sum(improvements > 0)} ({100*np.sum(improvements > 0)/len(results):.1f}%)")
    print(f"Galaxies worsened: {np.sum(improvements < 0)} ({100*np.sum(improvements < 0)/len(results):.1f}%)")

    # Success rate comparison
    excellent_orig = np.sum(chi2_orig < 2)
    good_orig = np.sum((chi2_orig >= 2) & (chi2_orig < 5))

    excellent_refined = np.sum(chi2_refined < 2)
    good_refined = np.sum((chi2_refined >= 2) & (chi2_refined < 5))

    print(f"\n### Success Rates ###")
    print("Original:")
    print(f"  Excellent (χ² < 2): {excellent_orig} ({100*excellent_orig/len(results):.1f}%)")
    print(f"  Good (2 ≤ χ² < 5): {good_orig} ({100*good_orig/len(results):.1f}%)")
    print(f"  Combined: {excellent_orig + good_orig} ({100*(excellent_orig + good_orig)/len(results):.1f}%)")

    print("\nRefined:")
    print(f"  Excellent (χ² < 2): {excellent_refined} ({100*excellent_refined/len(results):.1f}%)")
    print(f"  Good (2 ≤ χ² < 5): {good_refined} ({100*good_refined/len(results):.1f}%)")
    print(f"  Combined: {excellent_refined + good_refined} ({100*(excellent_refined + good_refined)/len(results):.1f}%)")

    # Top improvements
    print(f"\n### Top 10 Improvements ###")
    sorted_idx = np.argsort(improvements)[::-1]
    for i in range(min(10, len(results))):
        idx = sorted_idx[i]
        r = results[idx]
        print(f"{i+1}. {r['galaxy']}: Δχ² = {improvements[idx]:.2f} "
              f"({r['original']['chi2_red']:.2f} → {r['refined']['chi2_red']:.2f})")

    # Top worsenings
    print(f"\n### Top 10 Worsenings ###")
    for i in range(min(10, len(results))):
        idx = sorted_idx[-(i+1)]
        r = results[idx]
        print(f"{i+1}. {r['galaxy']}: Δχ² = {improvements[idx]:.2f} "
              f"({r['original']['chi2_red']:.2f} → {r['refined']['chi2_red']:.2f})")

    return {
        'chi2_orig': chi2_orig,
        'chi2_refined': chi2_refined,
        'improvements': improvements,
        'success_orig': (excellent_orig + good_orig) / len(results),
        'success_refined': (excellent_refined + good_refined) / len(results)
    }


def save_results(results, filename='session38_refined_coherence_results.json'):
    """Save results to JSON file."""
    output = {
        'session': 38,
        'date': '2025-11-22',
        'formula': 'C_vis = 1 - exp(-(ρ_vis/ρ_crit)^0.30)',
        'n_galaxies': len(results),
        'results': results
    }

    with open(filename, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n✓ Results saved to {filename}")


def main():
    """
    Main execution: Test refined coherence formula on SPARC galaxies.
    """
    print("\n" + "="*80)
    print("SYNCHRONISM SESSION #38")
    print("SPARC Dark Matter: Refined Coherence Formula Test")
    print("="*80)

    # Run comparison on all SPARC galaxies
    results = run_full_sparc_comparison()

    if results:
        # Analyze results
        stats = analyze_results(results)

        # Save results
        save_results(results)

        print("\n" + "="*80)
        print("SESSION #38 COMPLETE")
        print("="*80)
        print(f"\n✓ Tested refined coherence formula on {len(results)} galaxies")
        print(f"✓ Original model success: {stats['success_orig']*100:.1f}%")
        print(f"✓ Refined model success: {stats['success_refined']*100:.1f}%")

        if stats['success_refined'] > stats['success_orig']:
            print(f"\n🎯 IMPROVEMENT: +{(stats['success_refined'] - stats['success_orig'])*100:.1f}% success rate")
        elif stats['success_refined'] < stats['success_orig']:
            print(f"\n⚠ WARNING: -{(stats['success_orig'] - stats['success_refined'])*100:.1f}% success rate")
        else:
            print(f"\n→ No change in overall success rate")

    else:
        print("\n❌ No results generated - check data directory")


if __name__ == '__main__':
    main()
