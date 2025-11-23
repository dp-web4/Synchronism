#!/usr/bin/env python3
"""
Session #40: Extended ρ_crit Grid Search for Massive Galaxies

Motivation:
- Session #39 analysis revealed NGC massive spirals hitting ρ_crit = 100 M☉/pc² upper bound
- 63% of best-performing AND worst-performing NGC galaxies at bound
- Strong evidence that true optimal ρ_crit > 100 for massive systems

Changes from Session #38:
- Extended ρ_crit range: [0.01, 100] → [0.01, 10000] M☉/pc²
- Increased grid density: 30 → 40 points (better resolution)
- Target: NGC massive spirals that hit previous bound

Expected Impact:
- NGC galaxies: +5-10 pp additional success rate improvement
- UGC galaxies: +2-3 pp improvement (~30% at bound)
- F/DDO galaxies: ~0 pp (already optimal in low range)

Author: CBP Autonomous Synchronism Research
Date: 2025-11-22
Session: #40 - Empirical Validation Extension
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime

# Import SPARC validation utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    RefinedCoherencePredictor,
    fit_galaxy_refined,
    RealSPARCLoader,
    SPARCGalaxy
)

def fit_galaxy_extended(galaxy: SPARCGalaxy):
    """
    Fit refined coherence model with EXTENDED rho_crit range.

    Session #38: rho_crit ∈ [0.01, 100] M☉/pc² (30 points)
    Session #40: rho_crit ∈ [0.01, 10000] M☉/pc² (40 points)

    Returns comparison of Session #38 vs Session #40 fit.
    """
    # Session #38 grid (original)
    rho_crit_s38 = np.logspace(-2, 2, 30)  # [0.01, 100]
    result_s38 = fit_galaxy_refined(galaxy, rho_crit_values=rho_crit_s38)

    # Session #40 grid (extended)
    rho_crit_s40 = np.logspace(-2, 4, 40)  # [0.01, 10000]
    result_s40 = fit_galaxy_refined(galaxy, rho_crit_values=rho_crit_s40)

    return {
        'name': galaxy.name,
        'session38': result_s38,
        'session40': result_s40,
        'additional_improvement': result_s38['chi2_red'] - result_s40['chi2_red'],
        'rho_crit_changed': abs(result_s40['rho_crit'] - result_s38['rho_crit']) > 1.0,
        'hit_s38_bound': result_s38['rho_crit'] >= 95.0,  # Within 5% of bound
        'hit_s40_bound': result_s40['rho_crit'] >= 9500.0
    }


def identify_target_galaxies(session38_results_file='session38_refined_coherence_results.json'):
    """
    Identify NGC galaxies that hit Session #38 upper bound.

    These are priority targets for extended grid search.
    """
    with open(Path(__file__).parent / session38_results_file) as f:
        s38_data = json.load(f)

    targets = []
    for result in s38_data['results']:
        # Target galaxies with rho_crit close to upper bound (100)
        if result['refined_rho_crit'] >= 95.0:
            targets.append(result['name'])

    return targets


def run_extended_validation(target_only=True):
    """
    Run extended ρ_crit validation.

    If target_only=True: Only test galaxies that hit Session #38 bound
    If target_only=False: Retest all 175 SPARC galaxies with extended range
    """
    loader = RealSPARCLoader()
    all_galaxies = loader.load_all_galaxies()

    if target_only:
        # Get list of galaxies that hit bound
        target_names = identify_target_galaxies()
        galaxies = [g for g in all_galaxies if g.name in target_names]
        print(f"\nTargeting {len(galaxies)} galaxies that hit Session #38 bound")
    else:
        galaxies = all_galaxies
        print(f"\nRetesting all {len(galaxies)} SPARC galaxies with extended range")

    results = []

    print("\n" + "="*80)
    print("SESSION #40: EXTENDED ρ_crit GRID SEARCH")
    print("="*80)
    print("Session #38 range: ρ_crit ∈ [0.01, 100] M☉/pc² (30 points)")
    print("Session #40 range: ρ_crit ∈ [0.01, 10000] M☉/pc² (40 points)")
    print("\n" + "-"*80)

    for i, galaxy in enumerate(galaxies, 1):
        try:
            result = fit_galaxy_extended(galaxy)
            results.append(result)

            # Report significant improvements
            if result['additional_improvement'] > 1.0:
                print(f"  {galaxy.name}: Δχ² = {result['additional_improvement']:.2f} "
                      f"(ρ_crit: {result['session38']['rho_crit']:.1f} → "
                      f"{result['session40']['rho_crit']:.1f})")

            if i % 10 == 0 or i == len(galaxies):
                print(f"\nProcessed {i}/{len(galaxies)} galaxies...")

        except Exception as e:
            print(f"  Error processing {galaxy.name}: {e}")
            continue

    return results


def analyze_extended_results(results):
    """
    Analyze impact of extended ρ_crit range.
    """
    n_total = len(results)
    n_improved = sum(1 for r in results if r['additional_improvement'] > 0)
    n_significant = sum(1 for r in results if r['additional_improvement'] > 1.0)
    n_rho_changed = sum(1 for r in results if r['rho_crit_changed'])
    n_s38_bound = sum(1 for r in results if r['hit_s38_bound'])
    n_s40_bound = sum(1 for r in results if r['hit_s40_bound'])

    total_improvement = sum(r['additional_improvement'] for r in results)
    avg_improvement = total_improvement / n_total if n_total > 0 else 0

    # Success rate change (χ² < 5 threshold)
    s38_success = sum(1 for r in results if r['session38']['chi2_red'] < 5.0)
    s40_success = sum(1 for r in results if r['session40']['chi2_red'] < 5.0)

    print("\n" + "="*80)
    print("SESSION #40 EXTENDED GRID ANALYSIS")
    print("="*80)
    print(f"\nGalaxies tested: {n_total}")
    print(f"Hit Session #38 bound (ρ_crit ≥ 95): {n_s38_bound} ({100*n_s38_bound/n_total:.1f}%)")
    print(f"Hit Session #40 bound (ρ_crit ≥ 9500): {n_s40_bound} ({100*n_s40_bound/n_total:.1f}%)")
    print()
    print(f"Galaxies improved: {n_improved}/{n_total} ({100*n_improved/n_total:.1f}%)")
    print(f"Significant improvement (Δχ² > 1): {n_significant} ({100*n_significant/n_total:.1f}%)")
    print(f"ρ_crit changed substantially: {n_rho_changed} ({100*n_rho_changed/n_total:.1f}%)")
    print()
    print(f"Average additional improvement: Δχ² = {avg_improvement:.3f}")
    print(f"Total χ² reduction: {total_improvement:.1f}")
    print()
    print(f"Session #38 success rate (χ² < 5): {s38_success}/{n_total} ({100*s38_success/n_total:.1f}%)")
    print(f"Session #40 success rate (χ² < 5): {s40_success}/{n_total} ({100*s40_success/n_total:.1f}%)")
    print(f"Success rate improvement: +{100*(s40_success - s38_success)/n_total:.1f} pp")
    print()

    # Find biggest improvements
    top_improvements = sorted(results, key=lambda r: r['additional_improvement'], reverse=True)[:10]
    print("\nTop 10 Additional Improvements:")
    print("-" * 80)
    for r in top_improvements:
        print(f"  {r['name']:15s} Δχ² = {r['additional_improvement']:+7.2f}  "
              f"(ρ_crit: {r['session38']['rho_crit']:6.1f} → {r['session40']['rho_crit']:8.1f})")

    return {
        'n_total': n_total,
        'n_improved': n_improved,
        'n_significant': n_significant,
        'avg_improvement': avg_improvement,
        'success_rate_s38': 100 * s38_success / n_total,
        'success_rate_s40': 100 * s40_success / n_total,
        'success_improvement': 100 * (s40_success - s38_success) / n_total,
        'top_improvements': [{'name': r['name'], 'delta_chi2': r['additional_improvement']}
                            for r in top_improvements]
    }


def save_results(results, analysis, output_file='session40_extended_rho_crit_results.json'):
    """
    Save Session #40 extended grid results.
    """
    output_path = Path(__file__).parent / output_file

    data = {
        'session': 40,
        'date': datetime.now().strftime('%Y-%m-%d'),
        'description': 'Extended ρ_crit grid search [0.01, 10000] M☉/pc²',
        'rho_crit_range': {'min': 0.01, 'max': 10000, 'n_points': 40},
        'analysis': analysis,
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"\nResults saved to: {output_file}")
    return output_path


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Session #40 Extended ρ_crit Grid Search')
    parser.add_argument('--all', action='store_true',
                       help='Test all 175 galaxies (default: target bound-hitters only)')
    parser.add_argument('--limit', type=int, default=None,
                       help='Limit number of galaxies for testing')

    args = parser.parse_args()

    # Run extended validation
    results = run_extended_validation(target_only=not args.all)

    if args.limit:
        results = results[:args.limit]

    # Analyze results
    analysis = analyze_extended_results(results)

    # Save results
    save_results(results, analysis)

    print("\n" + "="*80)
    print("SESSION #40 COMPLETE")
    print("="*80)
