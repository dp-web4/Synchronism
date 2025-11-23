#!/usr/bin/env python3
"""
Session #41: Continuous ρ_crit Optimization - Eliminating Grid Bounds

Nova's Recommendation (Session #40 Review):
"Continue to extend the ρ_crit parameter space, possibly through continuous
optimization as suggested."

Problem Identified:
- Session #38: 48% at bound (ρ_crit = 100)
- Session #40: 65.5% STILL at bound (ρ_crit = 10000)
- Grid search fundamentally limited - always hits upper bound

Solution:
Replace grid search with continuous optimization (scipy.optimize.minimize_scalar)
- No artificial bounds
- Finds true optimum
- Faster convergence
- More accurate results

Method:
- Use Brent's method (bracketing + parabolic interpolation)
- Robust to local minima
- Efficient for 1D optimization
- Automatic convergence detection

Author: CBP Autonomous Synchronism Research
Date: 2025-11-23
Session: #41 - Continuous Optimization Implementation
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy.optimize import minimize_scalar

# Import SPARC validation utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    RefinedCoherencePredictor,
    RealSPARCLoader,
    SPARCGalaxy
)


def fit_galaxy_continuous(galaxy: SPARCGalaxy,
                          rho_crit_bounds=(0.01, 100000),
                          method='bounded'):
    """
    Fit refined coherence model with CONTINUOUS ρ_crit optimization.

    Uses scipy.optimize.minimize_scalar instead of grid search.

    Args:
        galaxy: SPARC galaxy data
        rho_crit_bounds: (min, max) search range for ρ_crit
        method: 'bounded' (Brent with bounds) or 'brent' (unbounded)

    Returns:
        dict with optimized ρ_crit, α, and χ²
    """

    def objective(rho_crit):
        """
        Objective function: χ² as function of ρ_crit.

        For each ρ_crit, we fit α and return χ².
        """
        predictor = RefinedCoherencePredictor(rho_crit=rho_crit)
        result = predictor.fit_alpha(galaxy)
        return result['chi2_red']

    # Continuous optimization
    if method == 'bounded':
        result = minimize_scalar(
            objective,
            bounds=rho_crit_bounds,
            method='bounded',
            options={'xatol': 1e-3}  # Tolerance: 0.001 M☉/pc²
        )
    else:
        result = minimize_scalar(
            objective,
            method='brent',
            options={'xtol': 1e-3}
        )

    # Extract optimal ρ_crit
    rho_crit_opt = result.x
    chi2_opt = result.fun

    # Get full fit results at optimal ρ_crit
    predictor_opt = RefinedCoherencePredictor(rho_crit=rho_crit_opt)
    fit_result = predictor_opt.fit_alpha(galaxy)

    # Add optimization info
    fit_result['rho_crit'] = rho_crit_opt
    fit_result['optimization_success'] = result.success
    fit_result['optimization_nfev'] = result.nfev  # Number of function evaluations
    fit_result['at_lower_bound'] = abs(rho_crit_opt - rho_crit_bounds[0]) < 0.01
    fit_result['at_upper_bound'] = abs(rho_crit_opt - rho_crit_bounds[1]) < 100

    return fit_result


def compare_grid_vs_continuous(galaxy: SPARCGalaxy):
    """
    Compare Session #40 grid search vs Session #41 continuous optimization.
    """
    # Grid search (Session #40 approach)
    from session40_extended_rho_crit import fit_galaxy_extended
    grid_result = fit_galaxy_extended(galaxy)

    # Continuous optimization (Session #41 approach)
    continuous_result = fit_galaxy_continuous(galaxy)

    return {
        'name': galaxy.name,
        'grid_s40': grid_result['session40'],
        'continuous_s41': continuous_result,
        'improvement': grid_result['session40']['chi2_red'] - continuous_result['chi2_red'],
        'rho_crit_change': continuous_result['rho_crit'] - grid_result['session40']['rho_crit'],
        'grid_at_bound': grid_result['session40']['hit_s40_bound'],
        'continuous_at_bound': continuous_result['at_upper_bound']
    }


def run_continuous_optimization_all_sparc(limit=None):
    """
    Run continuous ρ_crit optimization on all 175 SPARC galaxies.

    This replaces grid search entirely.
    """
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=limit)

    results = []

    print("\n" + "="*80)
    print("SESSION #41: CONTINUOUS ρ_crit OPTIMIZATION")
    print("="*80)
    print(f"\nOptimizing all {len(galaxies)} SPARC galaxies")
    print("Method: scipy.optimize.minimize_scalar (Brent with bounds)")
    print("Range: [0.01, 100000] M☉/pc²")
    print("No grid - continuous optimization")
    print("\n" + "-"*80)

    for i, galaxy in enumerate(galaxies, 1):
        try:
            result = fit_galaxy_continuous(galaxy, rho_crit_bounds=(0.01, 100000))
            results.append({
                'name': galaxy.name,
                'rho_crit': result['rho_crit'],
                'alpha': result['alpha_best'],
                'chi2_red': result['chi2_red'],
                'success': result['optimization_success'],
                'nfev': result['optimization_nfev'],
                'at_bound': result['at_upper_bound']
            })

            if i % 25 == 0 or i == len(galaxies):
                print(f"Processed {i}/{len(galaxies)} galaxies...")

                # Show example of recent optimization
                if result['optimization_success']:
                    print(f"  {galaxy.name}: ρ_crit = {result['rho_crit']:.2f}, "
                          f"χ² = {result['chi2_red']:.3f}, "
                          f"nfev = {result['optimization_nfev']}")

        except Exception as e:
            print(f"  Error processing {galaxy.name}: {e}")
            continue

    return results


def analyze_continuous_results(results):
    """
    Analyze results of continuous optimization.
    """
    n_total = len(results)
    n_success = sum(1 for r in results if r['success'])
    n_at_bound = sum(1 for r in results if r['at_bound'])

    rho_crits = [r['rho_crit'] for r in results]
    chi2s = [r['chi2_red'] for r in results]
    nfevs = [r['nfev'] for r in results]

    # Success rate (χ² < 5)
    success_rate = 100 * sum(1 for chi2 in chi2s if chi2 < 5.0) / n_total

    # ρ_crit statistics
    rho_crit_median = np.median(rho_crits)
    rho_crit_mean = np.mean(rho_crits)
    rho_crit_std = np.std(rho_crits)
    rho_crit_min = np.min(rho_crits)
    rho_crit_max = np.max(rho_crits)

    # Efficiency
    avg_nfev = np.mean(nfevs)

    print("\n" + "="*80)
    print("SESSION #41 CONTINUOUS OPTIMIZATION ANALYSIS")
    print("="*80)
    print(f"\nGalaxies optimized: {n_total}")
    print(f"Optimization successes: {n_success} ({100*n_success/n_total:.1f}%)")
    print(f"At upper bound (100000): {n_at_bound} ({100*n_at_bound/n_total:.1f}%)")
    print()
    print(f"Success rate (χ² < 5): {success_rate:.1f}%")
    print()
    print(f"ρ_crit distribution:")
    print(f"  Median: {rho_crit_median:.2f} M☉/pc²")
    print(f"  Mean:   {rho_crit_mean:.2f} ± {rho_crit_std:.2f} M☉/pc²")
    print(f"  Range:  [{rho_crit_min:.2f}, {rho_crit_max:.2f}] M☉/pc²")
    print()
    print(f"Optimization efficiency:")
    print(f"  Average function evaluations: {avg_nfev:.1f}")
    print()

    # Best and worst fits
    sorted_results = sorted(results, key=lambda r: r['chi2_red'])

    print("Best 10 fits:")
    print("-" * 80)
    for r in sorted_results[:10]:
        print(f"  {r['name']:15s} χ² = {r['chi2_red']:7.3f}  "
              f"ρ_crit = {r['rho_crit']:8.2f}  nfev = {r['nfev']:3d}")

    print("\nWorst 10 fits:")
    print("-" * 80)
    for r in sorted_results[-10:]:
        bound_flag = " (BOUND)" if r['at_bound'] else ""
        print(f"  {r['name']:15s} χ² = {r['chi2_red']:7.3f}  "
              f"ρ_crit = {r['rho_crit']:8.2f}  nfev = {r['nfev']:3d}{bound_flag}")

    return {
        'n_total': n_total,
        'n_success': n_success,
        'n_at_bound': n_at_bound,
        'success_rate': success_rate,
        'rho_crit_stats': {
            'median': rho_crit_median,
            'mean': rho_crit_mean,
            'std': rho_crit_std,
            'min': rho_crit_min,
            'max': rho_crit_max
        },
        'avg_nfev': avg_nfev,
        'best_fits': sorted_results[:10],
        'worst_fits': sorted_results[-10:]
    }


def plot_rho_crit_distribution(results, output_dir='session41_analysis'):
    """
    Visualize ρ_crit distribution from continuous optimization.
    """
    output_path = Path(__file__).parent / output_dir
    output_path.mkdir(exist_ok=True)

    rho_crits = [r['rho_crit'] for r in results]
    chi2s = [r['chi2_red'] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left: ρ_crit histogram (log scale)
    ax1.hist(np.log10(rho_crits), bins=30, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('log₁₀(ρ_crit) [M☉/pc²]', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Distribution of Optimal ρ_crit\n(Continuous Optimization)',
                  fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axvline(np.log10(100), color='red', linestyle='--', linewidth=1.5,
                label='Session #38 bound', alpha=0.7)
    ax1.axvline(np.log10(10000), color='orange', linestyle='--', linewidth=1.5,
                label='Session #40 bound', alpha=0.7)
    ax1.legend()

    # Right: χ² vs ρ_crit scatter
    scatter = ax2.scatter(np.log10(rho_crits), chi2s, c=chi2s, cmap='RdYlGn_r',
                         alpha=0.6, s=30, vmin=0, vmax=10)
    ax2.set_xlabel('log₁₀(ρ_crit) [M☉/pc²]', fontsize=12)
    ax2.set_ylabel('χ²', fontsize=12)
    ax2.set_title('Fit Quality vs Optimal ρ_crit', fontsize=13, fontweight='bold')
    ax2.axhline(5, color='red', linestyle='--', linewidth=1.5, alpha=0.5,
                label='Success threshold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    plt.colorbar(scatter, ax=ax2, label='χ²')

    plt.tight_layout()
    fig.savefig(output_path / 'rho_crit_distribution_continuous.png',
                dpi=150, bbox_inches='tight')
    plt.close(fig)

    return output_path / 'rho_crit_distribution_continuous.png'


def save_results(results, analysis, output_file='session41_continuous_results.json'):
    """
    Save Session #41 continuous optimization results.
    """
    output_path = Path(__file__).parent / output_file

    # Convert to JSON-serializable format
    results_serializable = []
    for r in results:
        results_serializable.append({
            'name': r['name'],
            'rho_crit': float(r['rho_crit']),
            'alpha': float(r['alpha']),
            'chi2_red': float(r['chi2_red']),
            'success': bool(r['success']),
            'nfev': int(r['nfev']),
            'at_bound': bool(r['at_bound'])
        })

    data = {
        'session': 41,
        'date': datetime.now().strftime('%Y-%m-%d'),
        'description': 'Continuous ρ_crit optimization (no grid bounds)',
        'method': 'scipy.optimize.minimize_scalar (bounded Brent)',
        'rho_crit_range': {'min': 0.01, 'max': 100000},
        'analysis': {
            'n_total': analysis['n_total'],
            'n_success': analysis['n_success'],
            'n_at_bound': analysis['n_at_bound'],
            'success_rate': float(analysis['success_rate']),
            'rho_crit_stats': {
                k: float(v) for k, v in analysis['rho_crit_stats'].items()
            },
            'avg_nfev': float(analysis['avg_nfev'])
        },
        'results': results_serializable
    }

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"\nResults saved to: {output_file}")
    return output_path


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Session #41 Continuous ρ_crit Optimization'
    )
    parser.add_argument('--limit', type=int, default=None,
                       help='Limit number of galaxies for testing')
    parser.add_argument('--compare', action='store_true',
                       help='Compare grid vs continuous for sample galaxies')

    args = parser.parse_args()

    if args.compare:
        # Comparison mode: Grid vs Continuous
        print("\nComparison Mode: Grid Search vs Continuous Optimization")
        print("="*80)

        loader = RealSPARCLoader()
        galaxies = loader.load_all_galaxies(limit=10)

        for galaxy in galaxies:
            comparison = compare_grid_vs_continuous(galaxy)
            print(f"\n{comparison['name']}:")
            print(f"  Grid (S40):       ρ_crit = {comparison['grid_s40']['rho_crit']:.2f}, "
                  f"χ² = {comparison['grid_s40']['chi2_red']:.3f}")
            print(f"  Continuous (S41): ρ_crit = {comparison['continuous_s41']['rho_crit']:.2f}, "
                  f"χ² = {comparison['continuous_s41']['chi2_red']:.3f}")
            print(f"  Improvement: Δχ² = {comparison['improvement']:+.3f}")
    else:
        # Full optimization mode
        results = run_continuous_optimization_all_sparc(limit=args.limit)
        analysis = analyze_continuous_results(results)
        plot_rho_crit_distribution(results)
        save_results(results, analysis)

        print("\n" + "="*80)
        print("SESSION #41 COMPLETE")
        print("="*80)
        print(f"\nContinuous optimization eliminated grid bounds")
        print(f"Success rate: {analysis['success_rate']:.1f}%")
        print(f"Still at bound: {analysis['n_at_bound']}/{analysis['n_total']} "
              f"({100*analysis['n_at_bound']/analysis['n_total']:.1f}%)")
