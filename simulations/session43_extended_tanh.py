#!/usr/bin/env python3
"""
Session #43 Track B: Extended γ Bounds for Tanh Coherence

Session #43 Track A Discovery:
- 154/175 galaxies (88%) are at γ bounds [0.10, 2.0]!
- Bound-hitting explains poor fits for many galaxies
- Failure cases (35.4%) correlate with massive galaxies (higher v_max)
- Current bounds [0.10, 2.0] are too restrictive

This Track Tests:
1. Extend γ bounds to [0.01, 5.0] (wider range)
2. Retest all 175 SPARC galaxies with tanh coherence
3. Compare success rate: old bounds vs extended bounds
4. If successful, retry combined predictor with new γ values

Hypothesis:
- Extended bounds will reduce bound-hitting from 88% to <30%
- Success rate will improve beyond 64.6%
- γ will become more predictable (not saturated at limits)

Author: CBP Autonomous Synchronism Research
Date: 2025-11-24
Session: #43 - Extended Tanh Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy.optimize import minimize, minimize_scalar
from typing import Dict

# Import SPARC validation utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    SynchronismPredictor,
    RealSPARCLoader,
    SPARCGalaxy
)


class ExtendedTanhPredictor(SynchronismPredictor):
    """
    Tanh coherence with extended γ bounds.

    C = tanh(γ × log(ρ/ρ_c + 1))

    Extended bounds: γ ∈ [0.01, 5.0] (was [0.10, 2.0])
    """

    def __init__(self, rho_crit: float = 1.0, gamma: float = 0.30, beta: float = 0.30):
        super().__init__(gamma=gamma, beta=beta)
        self.rho_crit = rho_crit
        self.tanh_gamma = gamma  # Rename to avoid confusion

    def compute_coherence(self, rho_vis: np.ndarray, rho_0=None) -> np.ndarray:
        """Compute tanh coherence with extended γ."""
        # Tanh of log form
        C_vis = np.tanh(self.tanh_gamma * np.log(rho_vis / self.rho_crit + 1.0))

        # Ensure 0 ≤ C < 1
        C_vis = np.clip(C_vis, 0.0, 0.99999)

        return C_vis


def fit_galaxy_extended_tanh(galaxy: SPARCGalaxy,
                              rho_crit_bounds=(0.01, 100000),
                              gamma_bounds=(0.01, 5.0)) -> Dict:
    """
    Fit tanh coherence with EXTENDED γ bounds.

    Args:
        galaxy: SPARC galaxy data
        rho_crit_bounds: bounds for ρ_crit [M☉/pc²]
        gamma_bounds: bounds for γ (EXTENDED from [0.10, 2.0] to [0.01, 5.0])

    Returns:
        Dictionary with optimized parameters and fit quality
    """

    def objective(params):
        """Objective: minimize χ²."""
        rho_crit, gamma = params

        predictor = ExtendedTanhPredictor(rho_crit=rho_crit, gamma=gamma)
        result = predictor.fit_alpha(galaxy)

        return result['chi2_red']

    # Initial guess (log-space center)
    x0 = [
        np.sqrt(rho_crit_bounds[0] * rho_crit_bounds[1]),
        np.sqrt(gamma_bounds[0] * gamma_bounds[1])
    ]

    # Bounds
    bounds = [rho_crit_bounds, gamma_bounds]

    # Optimize
    result = minimize(
        objective,
        x0=x0,
        bounds=bounds,
        method='L-BFGS-B',
        options={'ftol': 1e-6, 'maxiter': 100}
    )

    rho_crit_opt, gamma_opt = result.x
    chi2_opt = result.fun

    # Get full fit at optimal params
    predictor_opt = ExtendedTanhPredictor(rho_crit=rho_crit_opt, gamma=gamma_opt)
    fit_result = predictor_opt.fit_alpha(galaxy)

    # Check bounds
    rho_crit_at_lower = abs(rho_crit_opt - rho_crit_bounds[0]) < 0.01 * rho_crit_bounds[0]
    rho_crit_at_upper = abs(rho_crit_opt - rho_crit_bounds[1]) < 0.01 * rho_crit_bounds[1]
    gamma_at_lower = abs(gamma_opt - gamma_bounds[0]) < 0.01 * gamma_bounds[0]
    gamma_at_upper = abs(gamma_opt - gamma_bounds[1]) < 0.01 * gamma_bounds[1]

    return {
        'name': galaxy.name,
        'rho_crit': float(rho_crit_opt),
        'gamma': float(gamma_opt),
        'chi2_red': float(chi2_opt),
        'alpha': float(fit_result['alpha_best']),
        'nfev': int(result.nfev),
        'success': bool(result.success),
        'rho_crit_at_lower': bool(rho_crit_at_lower),
        'rho_crit_at_upper': bool(rho_crit_at_upper),
        'gamma_at_lower': bool(gamma_at_lower),
        'gamma_at_upper': bool(gamma_at_upper)
    }


def run_extended_tanh_analysis(limit=None):
    """
    Test extended tanh bounds on all SPARC galaxies.

    Compare with Session #42 results (γ bounds [0.10, 2.0]).
    """
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=limit)

    print("\n" + "="*80)
    print("SESSION #43 TRACK B: EXTENDED TANH ANALYSIS")
    print("="*80)
    print(f"\nTesting extended γ bounds on {len(galaxies)} SPARC galaxies")
    print(f"OLD bounds: γ ∈ [0.10, 2.0]")
    print(f"NEW bounds: γ ∈ [0.01, 5.0]")
    print("\n" + "-"*80)

    results = []

    for i, galaxy in enumerate(galaxies, 1):
        result = fit_galaxy_extended_tanh(galaxy)
        results.append(result)

        if i % 25 == 0:
            print(f"Processed {i}/{len(galaxies)} galaxies...")
            print(f"  {galaxy.name}: ρ_crit = {result['rho_crit']:.2f}, "
                  f"γ = {result['gamma']:.3f}, χ² = {result['chi2_red']:.2f}")

    # Analyze results
    print("\n" + "="*80)
    print("SESSION #43 EXTENDED TANH ANALYSIS")
    print("="*80)

    chi2_values = np.array([r['chi2_red'] for r in results])
    gamma_values = np.array([r['gamma'] for r in results])
    rho_crit_values = np.array([r['rho_crit'] for r in results])

    success_rate = 100.0 * np.sum(chi2_values < 5.0) / len(results)
    median_chi2 = np.median(chi2_values)

    # Bound-hitting
    gamma_at_lower = sum(1 for r in results if r['gamma_at_lower'])
    gamma_at_upper = sum(1 for r in results if r['gamma_at_upper'])
    rho_at_lower = sum(1 for r in results if r['rho_crit_at_lower'])
    rho_at_upper = sum(1 for r in results if r['rho_crit_at_upper'])

    print(f"\n### Extended Tanh Performance:")
    print(f"  Success rate (χ² < 5): {success_rate:.1f}% ({np.sum(chi2_values < 5.0)}/{len(results)})")
    print(f"  Median χ²: {median_chi2:.2f}")

    print(f"\n### γ Distribution (Extended [0.01, 5.0]):")
    print(f"  Median: {np.median(gamma_values):.3f}")
    print(f"  Mean:   {np.mean(gamma_values):.3f} ± {np.std(gamma_values):.3f}")
    print(f"  Range:  [{np.min(gamma_values):.3f}, {np.max(gamma_values):.3f}]")

    print(f"\n### Bound-Hitting (Extended Bounds):")
    print(f"  γ at lower (0.01): {gamma_at_lower} ({100*gamma_at_lower/len(results):.1f}%)")
    print(f"  γ at upper (5.0):  {gamma_at_upper} ({100*gamma_at_upper/len(results):.1f}%)")
    print(f"  ρ_crit at upper:   {rho_at_upper} ({100*rho_at_upper/len(results):.1f}%)")

    # Compare with Session #42
    print(f"\n### Comparison with Session #42 (γ ∈ [0.10, 2.0]):")
    print(f"  OLD success rate: 64.6%")
    print(f"  NEW success rate: {success_rate:.1f}%")
    print(f"  Improvement: {success_rate - 64.6:+.1f} pp")

    # Distribution of γ values
    print(f"\n### γ Value Categories:")
    very_low = np.sum(gamma_values < 0.10)
    low = np.sum((gamma_values >= 0.10) & (gamma_values < 0.50))
    medium = np.sum((gamma_values >= 0.50) & (gamma_values < 2.0))
    high = np.sum((gamma_values >= 2.0) & (gamma_values < 5.0))
    very_high = np.sum(gamma_values >= 5.0)

    print(f"  γ < 0.10 (very low):      {very_low} ({100*very_low/len(results):.1f}%) - NEEDED extended range!")
    print(f"  0.10 ≤ γ < 0.50 (low):    {low} ({100*low/len(results):.1f}%)")
    print(f"  0.50 ≤ γ < 2.0 (medium):  {medium} ({100*medium/len(results):.1f}%)")
    print(f"  2.0 ≤ γ < 5.0 (high):     {high} ({100*high/len(results):.1f}%) - NEEDED extended range!")
    print(f"  γ ≥ 5.0 (very high):      {very_high} ({100*very_high/len(results):.1f}%) - at bound")

    # Top 10 best fits
    print(f"\n### Top 10 Best Fits (Extended Tanh):")
    sorted_results = sorted(results, key=lambda r: r['chi2_red'])
    for i, r in enumerate(sorted_results[:10], 1):
        print(f"  {i:2d}. {r['name']:15s}  χ² = {r['chi2_red']:6.3f}  "
              f"γ = {r['gamma']:.3f}  ρ_crit = {r['rho_crit']:.2f}")

    # Worst 10 fits
    print(f"\n### Worst 10 Fits (Extended Tanh):")
    for i, r in enumerate(sorted_results[-10:], 1):
        bound_str = ""
        if r['gamma_at_lower'] or r['gamma_at_upper']:
            bound_str = " [γ BOUND]"
        if r['rho_crit_at_upper']:
            bound_str += " [ρ BOUND]"
        print(f"  {i:2d}. {r['name']:15s}  χ² = {r['chi2_red']:6.2f}  "
              f"γ = {r['gamma']:.3f}  ρ_crit = {r['rho_crit']:.2f}{bound_str}")

    # Save results
    output_path = Path(__file__).parent / 'session43_extended_tanh_results.json'

    save_data = {
        'session': 43,
        'track': 'B - Extended Tanh',
        'date': datetime.now().isoformat(),
        'bounds': {
            'gamma': [0.01, 5.0],
            'rho_crit': [0.01, 100000]
        },
        'n_galaxies': len(results),
        'summary': {
            'success_rate': float(success_rate),
            'median_chi2': float(median_chi2),
            'gamma_median': float(np.median(gamma_values)),
            'gamma_mean': float(np.mean(gamma_values)),
            'gamma_std': float(np.std(gamma_values)),
            'gamma_at_lower': gamma_at_lower,
            'gamma_at_upper': gamma_at_upper,
            'rho_crit_at_upper': rho_at_upper,
            'improvement_vs_s42': float(success_rate - 64.6)
        },
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(save_data, f, indent=2)

    print(f"\n\nResults saved to: {output_path}")

    return save_data


if __name__ == '__main__':
    results = run_extended_tanh_analysis(limit=None)

    print("\n" + "="*80)
    print("SESSION #43 TRACK B COMPLETE")
    print("="*80)
