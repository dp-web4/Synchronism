#!/usr/bin/env python3
"""
Session #43 Track A: γ Parameter Analysis - Physical Interpretation

Session #42 Results:
- Tanh coherence achieved 64.6% success with C = tanh(γ × log(ρ/ρ_c + 1))
- γ varies per galaxy (unlike fixed γ = 0.30 in baseline)
- Many galaxies hit γ bounds (108 at upper, 46 at lower)

Nova's Feedback:
- "The γ parameter, a crucial part of the Tanh function, lacks a clear physical interpretation"
- "Future sessions should focus on theoretical aspects"

This Track Investigates:
1. What is the distribution of optimal γ values?
2. Does γ correlate with galaxy properties (like ρ_crit does with v_max)?
3. Can we predict γ from galaxy properties?
4. Physical interpretation: What does γ represent?

Physical Hypotheses for γ:
- γ = "decoherence width" in log-density space
- γ ∝ galaxy compactness (ρ_vis distribution width)
- γ ∝ morphology (early vs late type)
- γ ∝ distance from bound (galaxies near bound may have γ issues)

Author: CBP Autonomous Synchronism Research
Date: 2025-11-24
Session: #43 - γ Parameter Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy import stats
from scipy.optimize import curve_fit

# Import SPARC validation utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    RealSPARCLoader,
    SPARCGalaxy
)


def load_session42_tanh_results():
    """Load Session #42 tanh coherence results."""
    results_path = Path(__file__).parent / 'session42_alternative_coherence_results.json'

    with open(results_path, 'r') as f:
        data = json.load(f)

    return data


def extract_tanh_gamma_values(s42_data):
    """Extract γ values from tanh results."""
    gamma_values = []
    rho_crit_values = []
    chi2_values = []
    galaxy_names = []
    at_upper = []
    at_lower = []

    for galaxy_result in s42_data['results']:
        tanh_result = galaxy_result['results']['tanh']

        gamma_values.append(tanh_result['params']['gamma'])
        rho_crit_values.append(tanh_result['params']['rho_crit'])
        chi2_values.append(tanh_result['chi2_red'])
        galaxy_names.append(galaxy_result['galaxy'])
        at_upper.append(tanh_result.get('gamma_at_upper', False))
        at_lower.append(tanh_result.get('gamma_at_lower', False))

    return {
        'gamma': np.array(gamma_values),
        'rho_crit': np.array(rho_crit_values),
        'chi2': np.array(chi2_values),
        'names': galaxy_names,
        'gamma_at_upper': np.array(at_upper),
        'gamma_at_lower': np.array(at_lower)
    }


def get_galaxy_properties():
    """Load SPARC galaxies and extract properties for correlation analysis."""
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()

    properties = {}
    for g in galaxies:
        v_max = np.max(g.v_obs)
        rho_vis = g.total_baryonic_density()  # It's a method!
        rho_vis_max = np.max(rho_vis)
        rho_vis_mean = np.mean(rho_vis)
        rho_vis_std = np.std(rho_vis)

        # Compactness proxy: ratio of max to mean density
        compactness = rho_vis_max / rho_vis_mean if rho_vis_mean > 0 else 0

        # Density spread: log10(max/min) - how many orders of magnitude
        rho_min = np.min(rho_vis[rho_vis > 0])
        density_spread = np.log10(rho_vis_max / rho_min) if rho_min > 0 else 0

        # Number of data points
        n_points = len(g.radius)

        # Radial extent
        r_max = np.max(g.radius)
        r_min = np.min(g.radius)

        properties[g.name] = {
            'v_max': v_max,
            'rho_vis_max': rho_vis_max,
            'rho_vis_mean': rho_vis_mean,
            'rho_vis_std': rho_vis_std,
            'compactness': compactness,
            'density_spread': density_spread,
            'n_points': n_points,
            'r_max': r_max,
            'r_min': r_min
        }

    return properties


def test_gamma_correlations(tanh_results, galaxy_props):
    """Test correlations between γ and galaxy properties."""

    # Match galaxies
    gamma_vals = []
    v_max_vals = []
    rho_max_vals = []
    compactness_vals = []
    density_spread_vals = []
    chi2_vals = []
    names = []

    for i, name in enumerate(tanh_results['names']):
        if name in galaxy_props:
            gamma_vals.append(tanh_results['gamma'][i])
            chi2_vals.append(tanh_results['chi2'][i])
            v_max_vals.append(galaxy_props[name]['v_max'])
            rho_max_vals.append(galaxy_props[name]['rho_vis_max'])
            compactness_vals.append(galaxy_props[name]['compactness'])
            density_spread_vals.append(galaxy_props[name]['density_spread'])
            names.append(name)

    gamma_vals = np.array(gamma_vals)
    v_max_vals = np.array(v_max_vals)
    rho_max_vals = np.array(rho_max_vals)
    compactness_vals = np.array(compactness_vals)
    density_spread_vals = np.array(density_spread_vals)
    chi2_vals = np.array(chi2_vals)

    # Filter out galaxies at γ bounds (unreliable for correlation analysis)
    not_at_bound = ~(tanh_results['gamma_at_upper'][:len(gamma_vals)] |
                     tanh_results['gamma_at_lower'][:len(gamma_vals)])

    print(f"\nγ Parameter Correlation Analysis")
    print(f"="*60)
    print(f"Total galaxies: {len(gamma_vals)}")
    print(f"At γ bounds (excluded): {(~not_at_bound).sum()}")
    print(f"Used for correlations: {not_at_bound.sum()}")

    correlations = {}

    # Test each property
    properties_to_test = [
        ('v_max', v_max_vals, 'Maximum velocity [km/s]'),
        ('rho_vis_max', rho_max_vals, 'Max visible density [M☉/pc²]'),
        ('compactness', compactness_vals, 'Compactness (ρ_max/ρ_mean)'),
        ('density_spread', density_spread_vals, 'Density spread [decades]'),
    ]

    print(f"\n### Spearman Correlations (γ vs property):\n")

    for prop_name, prop_vals, prop_label in properties_to_test:
        # Use only non-bound galaxies
        gamma_fit = gamma_vals[not_at_bound]
        prop_fit = prop_vals[not_at_bound]

        # Spearman correlation
        rho, p_val = stats.spearmanr(gamma_fit, prop_fit)

        correlations[prop_name] = {
            'spearman_rho': float(rho),
            'p_value': float(p_val),
            'label': prop_label
        }

        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(f"  γ vs {prop_label:40s}: ρ = {rho:+.3f} (p = {p_val:.2e}) {sig}")

    # Also test γ vs χ² (are higher γ associated with better/worse fits?)
    rho_chi2, p_chi2 = stats.spearmanr(gamma_vals[not_at_bound], chi2_vals[not_at_bound])
    print(f"\n  γ vs χ² (fit quality)                        : ρ = {rho_chi2:+.3f} (p = {p_chi2:.2e})")

    return correlations, {
        'gamma': gamma_vals,
        'v_max': v_max_vals,
        'compactness': compactness_vals,
        'density_spread': density_spread_vals,
        'chi2': chi2_vals,
        'not_at_bound': not_at_bound,
        'names': names
    }


def fit_gamma_predictor(data):
    """
    Attempt to fit γ = f(galaxy properties) predictor.

    Similar to ρ_crit = A × v_max^B from Session #41/42.
    """

    # Use only non-bound galaxies
    mask = data['not_at_bound']
    gamma = data['gamma'][mask]

    print(f"\n### γ Predictor Fitting")
    print(f"="*60)

    # Try: γ = A × v_max^B
    v_max = data['v_max'][mask]

    def power_law(x, A, B):
        return A * x**B

    try:
        popt_v, pcov_v = curve_fit(power_law, v_max, gamma, p0=[1.0, 0.0], maxfev=5000)
        gamma_pred_v = power_law(v_max, *popt_v)
        r2_v = 1 - np.sum((gamma - gamma_pred_v)**2) / np.sum((gamma - np.mean(gamma))**2)

        print(f"\n1. γ = A × v_max^B:")
        print(f"   A = {popt_v[0]:.4f}")
        print(f"   B = {popt_v[1]:.4f}")
        print(f"   R² = {r2_v:.3f}")
    except:
        popt_v = None
        r2_v = 0
        print(f"\n1. γ = A × v_max^B: FAILED TO FIT")

    # Try: γ = A × compactness^B
    compactness = data['compactness'][mask]

    try:
        popt_c, pcov_c = curve_fit(power_law, compactness, gamma, p0=[1.0, 0.0], maxfev=5000)
        gamma_pred_c = power_law(compactness, *popt_c)
        r2_c = 1 - np.sum((gamma - gamma_pred_c)**2) / np.sum((gamma - np.mean(gamma))**2)

        print(f"\n2. γ = A × compactness^B:")
        print(f"   A = {popt_c[0]:.4f}")
        print(f"   B = {popt_c[1]:.4f}")
        print(f"   R² = {r2_c:.3f}")
    except:
        popt_c = None
        r2_c = 0
        print(f"\n2. γ = A × compactness^B: FAILED TO FIT")

    # Try: γ = A × density_spread^B
    density_spread = data['density_spread'][mask]

    try:
        popt_d, pcov_d = curve_fit(power_law, density_spread, gamma, p0=[1.0, 0.0], maxfev=5000)
        gamma_pred_d = power_law(density_spread, *popt_d)
        r2_d = 1 - np.sum((gamma - gamma_pred_d)**2) / np.sum((gamma - np.mean(gamma))**2)

        print(f"\n3. γ = A × density_spread^B:")
        print(f"   A = {popt_d[0]:.4f}")
        print(f"   B = {popt_d[1]:.4f}")
        print(f"   R² = {r2_d:.3f}")
    except:
        popt_d = None
        r2_d = 0
        print(f"\n3. γ = A × density_spread^B: FAILED TO FIT")

    # Simple statistics
    print(f"\n### γ Distribution Statistics:")
    print(f"   Median γ: {np.median(gamma):.3f}")
    print(f"   Mean γ:   {np.mean(gamma):.3f} ± {np.std(gamma):.3f}")
    print(f"   Range:    [{np.min(gamma):.3f}, {np.max(gamma):.3f}]")

    return {
        'v_max': {'params': popt_v, 'r2': r2_v} if popt_v is not None else None,
        'compactness': {'params': popt_c, 'r2': r2_c} if popt_c is not None else None,
        'density_spread': {'params': popt_d, 'r2': r2_d} if popt_d is not None else None,
        'stats': {
            'median': float(np.median(gamma)),
            'mean': float(np.mean(gamma)),
            'std': float(np.std(gamma)),
            'min': float(np.min(gamma)),
            'max': float(np.max(gamma))
        }
    }


def analyze_failure_cases(tanh_results, galaxy_props):
    """
    Investigate why 35.4% of galaxies still fail (χ² > 5).

    Questions:
    - Are failures concentrated in certain galaxy types?
    - Do failures correlate with γ bounds?
    - Are failures more massive or less massive?
    """

    chi2 = tanh_results['chi2']
    gamma = tanh_results['gamma']
    names = tanh_results['names']

    success = chi2 < 5.0
    failure = chi2 >= 5.0

    print(f"\n### Failure Case Analysis")
    print(f"="*60)
    print(f"Successes (χ² < 5): {success.sum()} ({100*success.sum()/len(chi2):.1f}%)")
    print(f"Failures (χ² ≥ 5):  {failure.sum()} ({100*failure.sum()/len(chi2):.1f}%)")

    # Extract properties for success vs failure
    success_v_max = []
    failure_v_max = []
    success_gamma = []
    failure_gamma = []
    success_at_bound = []
    failure_at_bound = []

    for i, name in enumerate(names):
        if name in galaxy_props:
            at_bound = tanh_results['gamma_at_upper'][i] or tanh_results['gamma_at_lower'][i]

            if success[i]:
                success_v_max.append(galaxy_props[name]['v_max'])
                success_gamma.append(gamma[i])
                success_at_bound.append(at_bound)
            else:
                failure_v_max.append(galaxy_props[name]['v_max'])
                failure_gamma.append(gamma[i])
                failure_at_bound.append(at_bound)

    print(f"\n### v_max Distribution:")
    print(f"   Success median: {np.median(success_v_max):.1f} km/s")
    print(f"   Failure median: {np.median(failure_v_max):.1f} km/s")

    # T-test for v_max difference
    t_stat, t_pval = stats.ttest_ind(success_v_max, failure_v_max)
    print(f"   T-test p-value: {t_pval:.3e} {'*' if t_pval < 0.05 else ''}")

    print(f"\n### γ Distribution:")
    print(f"   Success median: {np.median(success_gamma):.3f}")
    print(f"   Failure median: {np.median(failure_gamma):.3f}")

    print(f"\n### Bound-Hitting:")
    print(f"   Success at γ bound: {sum(success_at_bound)} ({100*sum(success_at_bound)/len(success_at_bound):.1f}%)")
    print(f"   Failure at γ bound: {sum(failure_at_bound)} ({100*sum(failure_at_bound)/len(failure_at_bound):.1f}%)")

    # List worst failures
    print(f"\n### Worst 10 Failures:")
    sorted_idx = np.argsort(chi2)[::-1][:10]
    for rank, idx in enumerate(sorted_idx, 1):
        name = names[idx]
        v_max = galaxy_props.get(name, {}).get('v_max', 0)
        at_bound = tanh_results['gamma_at_upper'][idx] or tanh_results['gamma_at_lower'][idx]
        bound_str = " [γ BOUND]" if at_bound else ""
        print(f"   {rank:2d}. {name:15s}  χ² = {chi2[idx]:7.2f}  γ = {gamma[idx]:.3f}  v_max = {v_max:6.1f}{bound_str}")

    return {
        'success_v_max_median': float(np.median(success_v_max)),
        'failure_v_max_median': float(np.median(failure_v_max)),
        'v_max_ttest_pval': float(t_pval),
        'success_gamma_median': float(np.median(success_gamma)),
        'failure_gamma_median': float(np.median(failure_gamma)),
        'success_at_bound_frac': sum(success_at_bound) / len(success_at_bound),
        'failure_at_bound_frac': sum(failure_at_bound) / len(failure_at_bound)
    }


def run_gamma_analysis():
    """Main analysis function."""

    print("\n" + "="*80)
    print("SESSION #43 TRACK A: γ PARAMETER ANALYSIS")
    print("="*80)

    # Load Session #42 results
    s42_data = load_session42_tanh_results()
    print(f"\nLoaded Session #42 results: {s42_data['n_galaxies']} galaxies")

    # Extract tanh γ values
    tanh_results = extract_tanh_gamma_values(s42_data)
    print(f"Extracted γ values for {len(tanh_results['gamma'])} galaxies")

    # Load galaxy properties
    galaxy_props = get_galaxy_properties()
    print(f"Loaded properties for {len(galaxy_props)} SPARC galaxies")

    # Test correlations
    correlations, matched_data = test_gamma_correlations(tanh_results, galaxy_props)

    # Fit predictors
    predictors = fit_gamma_predictor(matched_data)

    # Analyze failures
    failure_analysis = analyze_failure_cases(tanh_results, galaxy_props)

    # Save results
    output_path = Path(__file__).parent / 'session43_gamma_analysis_results.json'

    save_data = {
        'session': 43,
        'track': 'A - γ Parameter Analysis',
        'date': datetime.now().isoformat(),
        'correlations': correlations,
        'predictors': {
            k: {'params': list(v['params']) if v and v['params'] is not None else None, 'r2': v['r2'] if v else None}
            for k, v in predictors.items() if k != 'stats'
        },
        'gamma_stats': predictors['stats'],
        'failure_analysis': failure_analysis
    }

    with open(output_path, 'w') as f:
        json.dump(save_data, f, indent=2)

    print(f"\n\nResults saved to: {output_path}")

    return save_data


if __name__ == '__main__':
    results = run_gamma_analysis()

    print("\n" + "="*80)
    print("SESSION #43 TRACK A COMPLETE")
    print("="*80)
