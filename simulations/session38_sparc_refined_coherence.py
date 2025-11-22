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
from synchronism_sparc_validation import SynchronismPredictor, SPARCValidator, SPARCGalaxy
from synchronism_real_sparc_validation import RealSPARCLoader

class RefinedCoherencePredictor(SynchronismPredictor):
    """
    Extended predictor with refined coherence formula.

    C_vis = 1 - exp(-(ρ_vis/ρ_crit)^γ)

    Never saturates to 1, maintains quantum residual.
    """

    def __init__(self, rho_crit: float = 1.0, gamma: float = 0.30, beta: float = 0.30):
        """Initialize with critical density parameter."""
        super().__init__(gamma=gamma, beta=beta)
        self.rho_crit = rho_crit

    def compute_coherence(self, rho_vis: np.ndarray, rho_0=None) -> np.ndarray:
        """
        Compute refined coherence: C_vis = 1 - exp(-(ρ_vis/ρ_crit)^γ)

        Properties:
        - C → 0 as ρ_vis → 0 (no coherence in vacuum)
        - C → 1 as ρ_vis → ∞ but NEVER reaches 1 exactly
        - Always leaves room for dark matter: (1 - C_vis) > 0
        """
        # Use instance's rho_crit instead of rho_0
        C_vis = 1.0 - np.exp(-(rho_vis / self.rho_crit) ** self.gamma)

        # Ensure 0 ≤ C < 1
        C_vis = np.clip(C_vis, 0.0, 0.99999)

        return C_vis


def fit_galaxy_refined(galaxy: SPARCGalaxy, rho_crit_values=np.logspace(-2, 2, 30)):
    """
    Fit refined coherence model to galaxy with grid search over rho_crit.

    Returns best-fit parameters and chi2.
    """
    best_chi2 = np.inf
    best_alpha = None
    best_rho_crit = None
    best_result = None

    # Grid search over rho_crit
    for rho_crit in rho_crit_values:
        # Create predictor with this rho_crit
        predictor = RefinedCoherencePredictor(rho_crit=rho_crit)
        validator = SPARCValidator(predictor)

        # Fit alpha for this galaxy
        result = validator.fit_single_galaxy(galaxy)

        if result['chi2_red'] < best_chi2:
            best_chi2 = result['chi2_red']
            best_alpha = result['alpha_best']
            best_rho_crit = rho_crit
            best_result = result

    # Add rho_crit to result
    best_result['rho_crit'] = best_rho_crit

    return best_result


def compare_models(galaxy: SPARCGalaxy):
    """
    Compare original vs refined coherence for single galaxy.
    """
    # Original model
    predictor_orig = SynchronismPredictor()
    validator_orig = SPARCValidator(predictor_orig)
    result_orig = validator_orig.fit_single_galaxy(galaxy)

    # Refined model (grid search)
    result_refined = fit_galaxy_refined(galaxy)

    return {
        'name': galaxy.name,
        'original': result_orig,
        'refined': result_refined,
        'improvement': result_orig['chi2_red'] - result_refined['chi2_red']
    }


def run_full_sparc_comparison(max_galaxies=None):
    """
    Run comparison on all available SPARC galaxies.
    """
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=max_galaxies)

    if not galaxies:
        print("Error: No SPARC galaxies loaded")
        return None

    results = []

    print("\n" + "="*80)
    print("SESSION #38: REFINED COHERENCE FORMULA TEST")
    print("="*80)
    print(f"\nTesting on {len(galaxies)} SPARC galaxies")
    print("Original: C_vis = (ρ_vis/ρ_0)^0.30")
    print("Refined:  C_vis = 1 - exp(-(ρ_vis/ρ_crit)^0.30)")
    print("\n" + "-"*80)

    for i, galaxy in enumerate(galaxies, 1):
        try:
            comparison = compare_models(galaxy)
            results.append(comparison)

            if i % 10 == 0 or i == len(galaxies):
                print(f"Processed {i}/{len(galaxies)} galaxies...")

        except Exception as e:
            print(f"  Error processing {galaxy.name}: {e}")
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
    # Convert numpy types to Python types for JSON serialization
    serializable_results = []
    for r in results:
        serializable_results.append({
            'name': r['name'],
            'original_chi2_red': float(r['original']['chi2_red']),
            'original_alpha': float(r['original']['alpha_best']),
            'refined_chi2_red': float(r['refined']['chi2_red']),
            'refined_alpha': float(r['refined']['alpha_best']),
            'refined_rho_crit': float(r['refined']['rho_crit']),
            'improvement': float(r['improvement'])
        })

    output = {
        'session': 38,
        'date': '2025-11-22',
        'formula': 'C_vis = 1 - exp(-(ρ_vis/ρ_crit)^0.30)',
        'n_galaxies': len(results),
        'results': serializable_results
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
