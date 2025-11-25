#!/usr/bin/env python3
"""
Session #44 Track A: Failure Population Analysis

Session #43 discovered that failure galaxies have significantly higher v_max:
- Success v_max median: 85.7 km/s
- Failure v_max median: 141.0 km/s
- T-test p-value: 6.9×10⁻⁵

This session investigates:
1. What makes high-v_max galaxies fail?
2. Is there a v_max threshold effect?
3. Are there other distinguishing properties?
4. Can we identify a sub-population where Synchronism excels?

Author: CBP Autonomous Synchronism Research
Date: 2025-11-24
Session: #44 - Failure Population Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy import stats

# Import utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    RealSPARCLoader,
    SPARCGalaxy
)
from session43_combined_predictor import (
    TanhCoherencePredictor,
    predict_rho_crit_from_v_max
)


def json_safe(obj):
    """Convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [json_safe(v) for v in obj]
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def load_virial_params():
    """Load virial scaling parameters from Session #42."""
    results_path = Path(__file__).parent / 'session42_virial_predictor_results.json'
    with open(results_path, 'r') as f:
        data = json.load(f)
    return data['virial_params']


def compute_galaxy_properties(galaxy: SPARCGalaxy) -> dict:
    """Compute various galaxy properties for analysis."""
    # Velocity properties
    v_obs = galaxy.v_obs
    v_max = np.max(v_obs)
    v_mean = np.mean(v_obs)
    v_flat = v_obs[-3:].mean() if len(v_obs) >= 3 else v_max  # Outer velocity

    # Density properties
    rho_vis = galaxy.total_baryonic_density()
    rho_max = np.max(rho_vis)
    rho_min = np.min(rho_vis[rho_vis > 0]) if np.any(rho_vis > 0) else 1e-10
    rho_mean = np.mean(rho_vis)

    # Compactness and spread
    compactness = rho_max / rho_mean if rho_mean > 0 else 0
    density_spread = np.log10(rho_max / rho_min) if rho_min > 0 else 0

    # Radial extent
    r_max = np.max(galaxy.radius)
    r_half = galaxy.radius[len(galaxy.radius) // 2]

    # Rotation curve shape
    v_rise = v_obs[len(v_obs) // 4] / v_max if v_max > 0 else 0  # How fast to rise
    v_decline = v_flat / v_max if v_max > 0 else 0  # How much decline (flat=1)

    # Number of data points
    n_points = len(galaxy.radius)

    return {
        'name': galaxy.name,
        'v_max': float(v_max),
        'v_mean': float(v_mean),
        'v_flat': float(v_flat),
        'rho_max': float(rho_max),
        'rho_mean': float(rho_mean),
        'rho_min': float(rho_min),
        'compactness': float(compactness),
        'density_spread': float(density_spread),
        'r_max': float(r_max),
        'r_half': float(r_half),
        'v_rise': float(v_rise),
        'v_decline': float(v_decline),
        'n_points': int(n_points)
    }


def run_best_predictor(galaxy: SPARCGalaxy, virial_params: dict) -> dict:
    """Run Session #43's best predictor (tanh + γ=2.0) on a galaxy."""
    v_max = np.max(galaxy.v_obs)
    rho_crit_pred = predict_rho_crit_from_v_max(v_max, virial_params)

    predictor = TanhCoherencePredictor(rho_crit=rho_crit_pred, gamma=2.0)
    fit_result = predictor.fit_alpha(galaxy)

    return {
        'alpha': float(fit_result['alpha_best']),
        'chi2_red': float(fit_result['chi2_red']),
        'success': fit_result['chi2_red'] < 5.0
    }


def analyze_failure_population():
    """Comprehensive analysis of why certain galaxies fail."""

    print("\n" + "="*80)
    print("SESSION #44 TRACK A: FAILURE POPULATION ANALYSIS")
    print("="*80)

    # Load data
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()
    virial_params = load_virial_params()

    print(f"\nLoaded {len(galaxies)} galaxies")
    print(f"Using best predictor: Tanh + γ=2.0 + virial ρ_crit")

    # Run analysis on all galaxies
    all_data = []

    for galaxy in galaxies:
        props = compute_galaxy_properties(galaxy)
        fit = run_best_predictor(galaxy, virial_params)

        all_data.append({
            **props,
            **fit
        })

    # Separate success and failure populations
    success = [d for d in all_data if d['success']]
    failure = [d for d in all_data if not d['success']]

    print(f"\nSuccess: {len(success)} galaxies ({100*len(success)/len(all_data):.1f}%)")
    print(f"Failure: {len(failure)} galaxies ({100*len(failure)/len(all_data):.1f}%)")

    # Property comparison
    print("\n" + "-"*80)
    print("PROPERTY COMPARISON: Success vs Failure")
    print("-"*80)

    properties = ['v_max', 'v_flat', 'v_decline', 'rho_max', 'compactness',
                  'density_spread', 'r_max', 'n_points']

    results = {}

    for prop in properties:
        s_vals = np.array([d[prop] for d in success])
        f_vals = np.array([d[prop] for d in failure])

        s_median = np.median(s_vals)
        f_median = np.median(f_vals)

        # Mann-Whitney U test (non-parametric)
        try:
            statistic, pval = stats.mannwhitneyu(s_vals, f_vals, alternative='two-sided')
        except:
            pval = 1.0

        sig = "*" if pval < 0.05 else ""
        sig = "**" if pval < 0.01 else sig
        sig = "***" if pval < 0.001 else sig

        print(f"  {prop:20s}: Success={s_median:10.2f}, Failure={f_median:10.2f}, p={pval:.2e} {sig}")

        results[prop] = {
            'success_median': float(s_median),
            'failure_median': float(f_median),
            'p_value': float(pval),
            'significant': bool(pval < 0.05)
        }

    # v_max threshold analysis
    print("\n" + "-"*80)
    print("V_MAX THRESHOLD ANALYSIS")
    print("-"*80)

    thresholds = [50, 75, 100, 125, 150, 175, 200]
    threshold_results = {}

    for thresh in thresholds:
        below = [d for d in all_data if d['v_max'] <= thresh]
        above = [d for d in all_data if d['v_max'] > thresh]

        if len(below) > 0:
            success_below = 100 * sum(1 for d in below if d['success']) / len(below)
        else:
            success_below = 0

        if len(above) > 0:
            success_above = 100 * sum(1 for d in above if d['success']) / len(above)
        else:
            success_above = 0

        threshold_results[thresh] = {
            'n_below': len(below),
            'n_above': len(above),
            'success_below': float(success_below),
            'success_above': float(success_above)
        }

        print(f"  v_max ≤ {thresh:3d}: {success_below:5.1f}% success ({len(below):3d} galaxies)")
        print(f"  v_max >  {thresh:3d}: {success_above:5.1f}% success ({len(above):3d} galaxies)")
        print()

    # Identify optimal v_max cutoff
    print("\n" + "-"*80)
    print("OPTIMAL V_MAX CUTOFF ANALYSIS")
    print("-"*80)

    v_max_values = np.array([d['v_max'] for d in all_data])
    v_max_range = np.arange(50, 220, 10)

    best_cutoff = 0
    best_success_below = 0

    cutoff_analysis = []

    for cutoff in v_max_range:
        below = [d for d in all_data if d['v_max'] <= cutoff]
        if len(below) >= 20:  # Require at least 20 galaxies
            success_rate = 100 * sum(1 for d in below if d['success']) / len(below)
            cutoff_analysis.append({
                'cutoff': float(cutoff),
                'n_galaxies': len(below),
                'success_rate': float(success_rate)
            })

            if success_rate > best_success_below:
                best_success_below = success_rate
                best_cutoff = cutoff

    print(f"\nBest cutoff: v_max ≤ {best_cutoff} km/s")

    below_best = [d for d in all_data if d['v_max'] <= best_cutoff]
    above_best = [d for d in all_data if d['v_max'] > best_cutoff]

    success_best_below = 100 * sum(1 for d in below_best if d['success']) / len(below_best)
    success_best_above = 100 * sum(1 for d in above_best if d['success']) / len(above_best) if above_best else 0

    print(f"  v_max ≤ {best_cutoff}: {success_best_below:.1f}% success ({len(below_best)} galaxies)")
    print(f"  v_max >  {best_cutoff}: {success_best_above:.1f}% success ({len(above_best)} galaxies)")

    # Worst failures analysis
    print("\n" + "-"*80)
    print("WORST FAILURE ANALYSIS (Top 10 by χ²)")
    print("-"*80)

    sorted_failures = sorted(failure, key=lambda x: x['chi2_red'], reverse=True)
    worst_10 = sorted_failures[:10]

    print(f"\n{'Galaxy':<15} {'v_max':>8} {'χ²':>10} {'compactness':>12}")
    print("-" * 50)

    for d in worst_10:
        print(f"{d['name']:<15} {d['v_max']:>8.1f} {d['chi2_red']:>10.1f} {d['compactness']:>12.1f}")

    # Success stories analysis
    print("\n" + "-"*80)
    print("BEST SUCCESS ANALYSIS (Top 10 by lowest χ²)")
    print("-"*80)

    sorted_success = sorted(success, key=lambda x: x['chi2_red'])
    best_10 = sorted_success[:10]

    print(f"\n{'Galaxy':<15} {'v_max':>8} {'χ²':>10} {'compactness':>12}")
    print("-" * 50)

    for d in best_10:
        print(f"{d['name']:<15} {d['v_max']:>8.1f} {d['chi2_red']:>10.2f} {d['compactness']:>12.1f}")

    # Save results
    output = {
        'session': 44,
        'track': 'A - Failure Population Analysis',
        'date': datetime.now().isoformat(),
        'summary': {
            'total_galaxies': len(all_data),
            'success_count': len(success),
            'failure_count': len(failure),
            'success_rate': 100 * len(success) / len(all_data)
        },
        'property_comparison': results,
        'threshold_analysis': threshold_results,
        'optimal_cutoff': {
            'v_max_cutoff': best_cutoff,
            'n_below': len(below_best),
            'n_above': len(above_best),
            'success_below': success_best_below,
            'success_above': success_best_above
        },
        'cutoff_curve': cutoff_analysis,
        'worst_failures': [
            {'name': d['name'], 'v_max': d['v_max'], 'chi2': d['chi2_red']}
            for d in worst_10
        ],
        'best_successes': [
            {'name': d['name'], 'v_max': d['v_max'], 'chi2': d['chi2_red']}
            for d in best_10
        ]
    }

    output_path = Path(__file__).parent / 'session44_failure_analysis_results.json'
    with open(output_path, 'w') as f:
        json.dump(json_safe(output), f, indent=2)

    print(f"\n\nResults saved to: {output_path}")

    return output


if __name__ == '__main__':
    results = analyze_failure_population()

    print("\n" + "="*80)
    print("SESSION #44 TRACK A COMPLETE")
    print("="*80)

    print(f"\nKey finding: Synchronism with v_max ≤ {results['optimal_cutoff']['v_max_cutoff']} km/s:")
    print(f"  Success rate: {results['optimal_cutoff']['success_below']:.1f}%")
    print(f"  Coverage: {results['optimal_cutoff']['n_below']} / {results['summary']['total_galaxies']} galaxies")
