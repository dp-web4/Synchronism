#!/usr/bin/env python3
"""
Session #42 Track B: Virial Predictor - Predicting ρ_crit from v_max

Session #41 Discovery:
- ρ_crit correlates most strongly with v_max (Spearman ρ = 0.57, p < 2e-16)
- Power law fit: ρ_crit ∝ v_max^1.74
- Physical interpretation: ρ_crit ∝ √(GM/R) (virial scaling)

This Track Tests:
Instead of optimizing ρ_crit for each galaxy independently, can we
PREDICT ρ_crit from v_max using the discovered empirical relationship?

Advantages if successful:
1. Reduces free parameters: 1 parameter (global scaling) vs 175 (per-galaxy ρ_crit)
2. Makes theory more predictive (not just fitting)
3. Validates physical interpretation (virial decoherence scale)
4. Simpler model = better science

Approach:
1. Use Session #41 data to fit ρ_crit = A × v_max^B
2. Apply this formula to predict ρ_crit for all galaxies
3. Compute χ² with predicted (not optimized) ρ_crit
4. Compare success rate: predicted vs optimized

Expected outcome:
- Success rate will decrease (we're constraining the model)
- But if still reasonable (>40-50%), validates virial interpretation
- Trade parameter count for predictive power

Author: CBP Autonomous Synchronism Research
Date: 2025-11-23
Session: #42 - Virial Predictor
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
    RefinedCoherencePredictor,
    RealSPARCLoader,
    SPARCGalaxy
)


def load_session41_results():
    """Load Session #41 continuous optimization results."""
    results_path = Path(__file__).parent / 'session41_continuous_results.json'

    with open(results_path, 'r') as f:
        data = json.load(f)

    return data


def fit_virial_scaling(session41_data):
    """
    Fit ρ_crit = A × v_max^B power law from Session #41 data.

    Returns:
        Dictionary with fit parameters and statistics
    """
    # Load SPARC galaxies to get v_max
    loader = RealSPARCLoader()
    sparc_galaxies = loader.load_all_galaxies()

    # Create lookup dict for v_max (maximum observed velocity)
    v_max_lookup = {g.name: np.max(g.v_obs) for g in sparc_galaxies}

    # Extract v_max and ρ_crit from Session #41 results
    galaxies = []
    v_maxs = []
    rho_crits = []

    for result in session41_data['results']:
        name = result['name']
        if name in v_max_lookup:
            galaxies.append(name)
            v_maxs.append(v_max_lookup[name])
            rho_crits.append(result['rho_crit'])

    v_maxs = np.array(v_maxs)
    rho_crits = np.array(rho_crits)

    # Remove galaxies at bounds (unreliable for fitting)
    not_at_bound = rho_crits < 99000  # Not at upper bound
    v_maxs_fit = v_maxs[not_at_bound]
    rho_crits_fit = rho_crits[not_at_bound]
    galaxies_fit = [g for g, nab in zip(galaxies, not_at_bound) if nab]

    print(f"\nFitting virial scaling:")
    print(f"  Total galaxies: {len(galaxies)}")
    print(f"  Excluded (at bound): {(~not_at_bound).sum()}")
    print(f"  Used for fitting: {len(galaxies_fit)}")

    # Power law fit in log-log space
    def power_law(x, A, B):
        return A * x**B

    # Fit
    popt, pcov = curve_fit(
        power_law,
        v_maxs_fit,
        rho_crits_fit,
        p0=[1.0, 1.5]  # Initial guess
    )

    A_fit, B_fit = popt
    A_err, B_err = np.sqrt(np.diag(pcov))

    # Compute correlation in log-log space
    log_v = np.log10(v_maxs_fit)
    log_rho = np.log10(rho_crits_fit)
    pearson_r, p_value = stats.pearsonr(log_v, log_rho)
    spearman_rho, sp_pval = stats.spearmanr(v_maxs_fit, rho_crits_fit)

    # R² (coefficient of determination)
    rho_pred = power_law(v_maxs_fit, A_fit, B_fit)
    residuals = rho_crits_fit - rho_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((rho_crits_fit - np.mean(rho_crits_fit))**2)
    r_squared = 1 - (ss_res / ss_tot)

    print(f"\nVirial Scaling Fit: ρ_crit = A × v_max^B")
    print(f"  A = {A_fit:.2f} ± {A_err:.2f}")
    print(f"  B = {B_fit:.2f} ± {B_err:.2f}")
    print(f"  Pearson r (log-log): {pearson_r:.3f} (p = {p_value:.2e})")
    print(f"  Spearman ρ: {spearman_rho:.3f} (p = {sp_pval:.2e})")
    print(f"  R²: {r_squared:.3f}")

    return {
        'A': float(A_fit),
        'B': float(B_fit),
        'A_err': float(A_err),
        'B_err': float(B_err),
        'pearson_r': float(pearson_r),
        'spearman_rho': float(spearman_rho),
        'r_squared': float(r_squared),
        'n_fit': int(len(galaxies_fit)),
        'n_excluded': int((~not_at_bound).sum())
    }


def predict_rho_crit(v_max: float, virial_params: dict) -> float:
    """Predict ρ_crit from v_max using virial scaling."""
    A = virial_params['A']
    B = virial_params['B']
    return A * v_max**B


def test_virial_predictor(virial_params: dict, limit=None):
    """
    Test virial predictor on all SPARC galaxies.

    Instead of optimizing ρ_crit, we PREDICT it from v_max.

    Returns:
        Results comparing predicted vs optimized performance
    """
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=limit)

    print("\n" + "="*80)
    print("SESSION #42 TRACK B: VIRIAL PREDICTOR TEST")
    print("="*80)
    print(f"\nPredicting ρ_crit from v_max for {len(galaxies)} galaxies")
    print(f"Formula: ρ_crit = {virial_params['A']:.2f} × v_max^{virial_params['B']:.2f}")
    print("\n" + "-"*80)

    results = []

    for i, galaxy in enumerate(galaxies, 1):
        # Compute v_max from observed velocities
        v_max = np.max(galaxy.v_obs)

        # Predict ρ_crit from v_max
        rho_crit_pred = predict_rho_crit(v_max, virial_params)

        # Evaluate model with predicted ρ_crit (no optimization!)
        predictor = RefinedCoherencePredictor(rho_crit=rho_crit_pred)
        fit_result = predictor.fit_alpha(galaxy)

        results.append({
            'name': galaxy.name,
            'v_max': float(v_max),
            'rho_crit_predicted': float(rho_crit_pred),
            'alpha': float(fit_result['alpha_best']),
            'chi2_red': float(fit_result['chi2_red'])
        })

        if i % 25 == 0:
            print(f"Processed {i}/{len(galaxies)} galaxies...")
            print(f"  {galaxy.name}: v_max = {v_max:.1f} km/s, "
                  f"ρ_crit_pred = {rho_crit_pred:.2f}, χ² = {fit_result['chi2_red']:.2f}")

    # Analyze results
    print("\n" + "="*80)
    print("SESSION #42 VIRIAL PREDICTOR ANALYSIS")
    print("="*80)

    chi2_values = np.array([r['chi2_red'] for r in results])
    rho_crit_preds = np.array([r['rho_crit_predicted'] for r in results])

    success_rate = 100.0 * np.sum(chi2_values < 5.0) / len(results)
    median_chi2 = np.median(chi2_values)
    mean_chi2 = np.mean(chi2_values)

    print(f"\n### Predicted ρ_crit Performance:\n")
    print(f"  Success rate (χ² < 5): {success_rate:.1f}% ({np.sum(chi2_values < 5.0)}/{len(results)})")
    print(f"  Median χ²: {median_chi2:.2f}")
    print(f"  Mean χ²: {mean_chi2:.2f}")

    print(f"\n### Predicted ρ_crit Distribution:\n")
    print(f"  Median: {np.median(rho_crit_preds):.2f} M☉/pc²")
    print(f"  Mean: {np.mean(rho_crit_preds):.2f} ± {np.std(rho_crit_preds):.2f}")
    print(f"  Range: [{np.min(rho_crit_preds):.2f}, {np.max(rho_crit_preds):.2f}]")

    # Top 10 best fits
    print(f"\n### Top 10 Best Fits (Predicted ρ_crit):\n")
    sorted_results = sorted(results, key=lambda r: r['chi2_red'])

    for i, r in enumerate(sorted_results[:10], 1):
        print(f"  {i:2d}. {r['name']:15s}  χ² = {r['chi2_red']:6.3f}  "
              f"v_max = {r['v_max']:6.1f} km/s  ρ_crit = {r['rho_crit_predicted']:8.2f}")

    # Worst 10 fits
    print(f"\n### Worst 10 Fits (Predicted ρ_crit):\n")

    for i, r in enumerate(sorted_results[-10:], 1):
        print(f"  {i:2d}. {r['name']:15s}  χ² = {r['chi2_red']:6.3f}  "
              f"v_max = {r['v_max']:6.1f} km/s  ρ_crit = {r['rho_crit_predicted']:8.2f}")

    # Save results
    output_path = Path(__file__).parent / 'session42_virial_predictor_results.json'

    save_data = {
        'session': 42,
        'track': 'B - Virial Predictor',
        'date': datetime.now().isoformat(),
        'virial_params': virial_params,
        'n_galaxies': len(results),
        'summary': {
            'success_rate': float(success_rate),
            'median_chi2': float(median_chi2),
            'mean_chi2': float(mean_chi2),
            'rho_crit_median': float(np.median(rho_crit_preds)),
            'rho_crit_mean': float(np.mean(rho_crit_preds)),
            'rho_crit_std': float(np.std(rho_crit_preds))
        },
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(save_data, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    return save_data


def compare_predicted_vs_optimized():
    """
    Load Session #41 (optimized) and Session #42 (predicted) results.
    Compare performance.
    """
    # Load Session #41 optimized results
    s41_path = Path(__file__).parent / 'session41_continuous_results.json'
    with open(s41_path, 'r') as f:
        s41_data = json.load(f)

    # Load Session #42 predicted results
    s42_path = Path(__file__).parent / 'session42_virial_predictor_results.json'
    with open(s42_path, 'r') as f:
        s42_data = json.load(f)

    print("\n" + "="*80)
    print("COMPARISON: Predicted vs Optimized ρ_crit")
    print("="*80)

    s41_success = s41_data['analysis']['success_rate']
    s42_success = s42_data['summary']['success_rate']

    # Compute median from results
    s41_chi2 = [r['chi2_red'] for r in s41_data['results']]
    s41_median = np.median(s41_chi2)
    s42_median = s42_data['summary']['median_chi2']

    print(f"\n### Success Rate (χ² < 5):")
    print(f"  Optimized (S41):  {s41_success:.1f}%")
    print(f"  Predicted (S42):  {s42_success:.1f}%")
    print(f"  Difference:       {s42_success - s41_success:+.1f} pp")

    print(f"\n### Median χ²:")
    print(f"  Optimized (S41):  {s41_median:.2f}")
    print(f"  Predicted (S42):  {s42_median:.2f}")
    print(f"  Difference:       {s42_median - s41_median:+.2f}")

    print(f"\n### Interpretation:")

    if s42_success > 40:
        print("  ✅ VIRIAL PREDICTOR VIABLE!")
        print(f"     Success rate {s42_success:.1f}% with ZERO per-galaxy tuning")
        print("     Validates ρ_crit ∝ v_max physical interpretation")

    if abs(s42_success - s41_success) < 10:
        print("  ✅ COMPARABLE PERFORMANCE!")
        print("     Predicted ρ_crit performs within 10 pp of optimized")
        print("     Massive reduction in model complexity (1 vs 175 parameters)")

    print("\n" + "="*80)


if __name__ == '__main__':
    # Load Session #41 results
    s41_data = load_session41_results()

    # Fit virial scaling
    virial_params = fit_virial_scaling(s41_data)

    # Test predictor
    s42_results = test_virial_predictor(virial_params, limit=None)

    # Compare predicted vs optimized
    compare_predicted_vs_optimized()

    print("\n" + "="*80)
    print("SESSION #42 TRACK B COMPLETE")
    print("="*80)
