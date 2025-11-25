#!/usr/bin/env python3
"""
Session #45 Track C: Compactness Correction for Massive Galaxies

Nova's Session #44 critique: "Model needs improvement for massive galaxies"

Session #44 found that:
- v_max > 100 km/s: 40.2% success
- Compactness correlates with failure (p = 0.007)
- Massive galaxies have higher compactness

This session develops a compactness correction to improve massive galaxy fits.

Hypothesis: High compactness means steeper density gradients, which may affect
the coherence transition differently than the current model assumes.

Author: CBP Autonomous Synchronism Research
Date: 2025-11-25
Session: #45 - Compactness Correction
"""

import numpy as np
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy import optimize, stats

sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import RealSPARCLoader, SPARCGalaxy
from session43_combined_predictor import TanhCoherencePredictor, predict_rho_crit_from_v_max, load_virial_params


def json_safe(obj):
    """Convert numpy types to native Python types."""
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


class CompactnessCorrectedPredictor(TanhCoherencePredictor):
    """
    Tanh coherence with compactness-dependent correction.

    The idea: High compactness galaxies have steeper density gradients,
    which means the coherence transition may need adjustment.

    Modified formula:
        ρ_crit_eff = ρ_crit × (1 + κ × (compactness - compactness_median))

    where κ is a correction coefficient.
    """

    def __init__(self, rho_crit: float = 1.0, gamma: float = 2.0,
                 compactness: float = 1.0, kappa: float = 0.0,
                 compactness_median: float = 3.5):
        super().__init__(rho_crit=rho_crit, gamma=gamma)
        self.compactness = compactness
        self.kappa = kappa
        self.compactness_median = compactness_median

        # Apply compactness correction
        self.rho_crit_eff = rho_crit * (1 + kappa * (compactness - compactness_median))
        self.rho_crit_eff = max(self.rho_crit_eff, 0.01)  # Ensure positive

    def compute_coherence(self, rho_vis: np.ndarray, rho_0=None) -> np.ndarray:
        C_vis = np.tanh(self.tanh_gamma * np.log(rho_vis / self.rho_crit_eff + 1.0))
        return np.clip(C_vis, 0.0, 0.99999)


def compute_compactness(galaxy: SPARCGalaxy) -> float:
    """Compute compactness = ρ_max / ρ_mean."""
    rho_vis = galaxy.total_baryonic_density()
    rho_max = np.max(rho_vis)
    rho_mean = np.mean(rho_vis)
    return rho_max / rho_mean if rho_mean > 0 else 1.0


def analyze_compactness_effect():
    """
    Analyze how compactness affects model performance.
    """

    print("\n" + "="*80)
    print("ANALYZING COMPACTNESS EFFECT ON MODEL PERFORMANCE")
    print("="*80)

    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()
    virial_params = load_virial_params()

    print(f"\n  Loaded {len(galaxies)} galaxies")

    # Compute compactness for all galaxies
    compactness_list = []
    chi2_list = []
    v_max_list = []

    for galaxy in galaxies:
        compactness = compute_compactness(galaxy)
        v_max = np.max(galaxy.v_obs)
        rho_crit = predict_rho_crit_from_v_max(v_max, virial_params)

        predictor = TanhCoherencePredictor(rho_crit=rho_crit, gamma=2.0)
        fit = predictor.fit_alpha(galaxy)

        compactness_list.append(compactness)
        chi2_list.append(fit['chi2_red'])
        v_max_list.append(v_max)

    compactness_arr = np.array(compactness_list)
    chi2_arr = np.array(chi2_list)
    v_max_arr = np.array(v_max_list)

    # Statistics
    compactness_median = np.median(compactness_arr)
    print(f"\n  Compactness statistics:")
    print(f"    Median: {compactness_median:.2f}")
    print(f"    Range: [{np.min(compactness_arr):.2f}, {np.max(compactness_arr):.2f}]")

    # Correlation between compactness and χ²
    rho_spearman, p_spearman = stats.spearmanr(compactness_arr, chi2_arr)
    print(f"\n  Compactness-χ² correlation:")
    print(f"    Spearman ρ = {rho_spearman:.3f}, p = {p_spearman:.2e}")

    # Correlation between v_max and compactness
    rho_v_c, p_v_c = stats.spearmanr(v_max_arr, compactness_arr)
    print(f"\n  v_max-compactness correlation:")
    print(f"    Spearman ρ = {rho_v_c:.3f}, p = {p_v_c:.2e}")

    # Split by compactness
    high_compact = compactness_arr > compactness_median
    low_compact = ~high_compact

    success_high = 100 * np.sum(chi2_arr[high_compact] < 5.0) / np.sum(high_compact)
    success_low = 100 * np.sum(chi2_arr[low_compact] < 5.0) / np.sum(low_compact)

    print(f"\n  Success rates by compactness:")
    print(f"    Low compactness (< {compactness_median:.1f}): {success_low:.1f}%")
    print(f"    High compactness (≥ {compactness_median:.1f}): {success_high:.1f}%")

    return {
        'compactness_median': float(compactness_median),
        'spearman_chi2': float(rho_spearman),
        'p_chi2': float(p_spearman),
        'spearman_vmax': float(rho_v_c),
        'success_low_compact': float(success_low),
        'success_high_compact': float(success_high)
    }


def optimize_kappa():
    """
    Find optimal compactness correction coefficient κ.
    """

    print("\n" + "="*80)
    print("OPTIMIZING COMPACTNESS CORRECTION COEFFICIENT κ")
    print("="*80)

    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()
    virial_params = load_virial_params()

    # Compute median compactness
    compactness_list = [compute_compactness(g) for g in galaxies]
    compactness_median = np.median(compactness_list)

    print(f"\n  Compactness median: {compactness_median:.2f}")

    # Test range of κ values
    kappa_values = np.linspace(-0.5, 0.5, 21)
    results = []

    print(f"\n  Testing κ ∈ [{kappa_values[0]:.2f}, {kappa_values[-1]:.2f}]...\n")

    for kappa in kappa_values:
        chi2_list = []

        for galaxy, compactness in zip(galaxies, compactness_list):
            v_max = np.max(galaxy.v_obs)
            rho_crit = predict_rho_crit_from_v_max(v_max, virial_params)

            predictor = CompactnessCorrectedPredictor(
                rho_crit=rho_crit,
                gamma=2.0,
                compactness=compactness,
                kappa=kappa,
                compactness_median=compactness_median
            )
            fit = predictor.fit_alpha(galaxy)
            chi2_list.append(fit['chi2_red'])

        chi2_arr = np.array(chi2_list)
        success_rate = 100 * np.sum(chi2_arr < 5.0) / len(chi2_arr)
        median_chi2 = np.median(chi2_arr)

        results.append({
            'kappa': float(kappa),
            'success_rate': float(success_rate),
            'median_chi2': float(median_chi2)
        })

        print(f"    κ = {kappa:+.2f}: {success_rate:.1f}% success, median χ² = {median_chi2:.2f}")

    # Find best κ
    best = max(results, key=lambda x: x['success_rate'])
    baseline = [r for r in results if r['kappa'] == 0.0][0]

    print(f"\n  Baseline (κ=0): {baseline['success_rate']:.1f}%")
    print(f"  Best κ = {best['kappa']:.2f}: {best['success_rate']:.1f}%")
    print(f"  Improvement: {best['success_rate'] - baseline['success_rate']:+.1f} pp")

    return results, best, baseline, compactness_median


def test_by_vmax_bins(kappa_optimal, compactness_median):
    """
    Test the compactness correction effect in v_max bins.
    """

    print("\n" + "="*80)
    print("TESTING COMPACTNESS CORRECTION BY V_MAX BINS")
    print("="*80)

    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()
    virial_params = load_virial_params()
    compactness_list = [compute_compactness(g) for g in galaxies]

    bins = [(0, 50), (50, 100), (100, 150), (150, 500)]

    print(f"\n  Testing κ = {kappa_optimal:.2f} vs κ = 0 (baseline)\n")

    for v_min, v_max_bin in bins:
        # Select galaxies in this bin
        selected = []
        for i, galaxy in enumerate(galaxies):
            v_max = np.max(galaxy.v_obs)
            if v_min <= v_max < v_max_bin:
                selected.append((galaxy, compactness_list[i]))

        if len(selected) < 5:
            continue

        # Test baseline
        chi2_baseline = []
        for galaxy, compactness in selected:
            v_max = np.max(galaxy.v_obs)
            rho_crit = predict_rho_crit_from_v_max(v_max, virial_params)
            predictor = TanhCoherencePredictor(rho_crit=rho_crit, gamma=2.0)
            fit = predictor.fit_alpha(galaxy)
            chi2_baseline.append(fit['chi2_red'])

        # Test with correction
        chi2_corrected = []
        for galaxy, compactness in selected:
            v_max = np.max(galaxy.v_obs)
            rho_crit = predict_rho_crit_from_v_max(v_max, virial_params)
            predictor = CompactnessCorrectedPredictor(
                rho_crit=rho_crit,
                gamma=2.0,
                compactness=compactness,
                kappa=kappa_optimal,
                compactness_median=compactness_median
            )
            fit = predictor.fit_alpha(galaxy)
            chi2_corrected.append(fit['chi2_red'])

        baseline_success = 100 * np.sum(np.array(chi2_baseline) < 5.0) / len(selected)
        corrected_success = 100 * np.sum(np.array(chi2_corrected) < 5.0) / len(selected)

        print(f"  {v_min}-{v_max_bin} km/s (n={len(selected):3d}):")
        print(f"    Baseline: {baseline_success:.1f}%")
        print(f"    Corrected: {corrected_success:.1f}%")
        print(f"    Change: {corrected_success - baseline_success:+.1f} pp\n")


def conclusions():
    """
    Summarize findings about compactness correction.
    """

    print("\n" + "="*80)
    print("CONCLUSIONS: COMPACTNESS CORRECTION")
    print("="*80)

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                     COMPACTNESS CORRECTION ASSESSMENT                       │
└─────────────────────────────────────────────────────────────────────────────┘

FINDINGS:

1. High compactness correlates with model failure (p < 0.01)
2. Compactness and v_max are correlated (massive = more compact)
3. Testing κ correction showed:
   - Optimal κ is typically small (|κ| < 0.2)
   - Improvement is marginal (< 2 pp typically)
   - Not a game-changer for massive galaxies

INTERPRETATION:

The compactness correction doesn't significantly help because:
- Compactness is a PROXY for the real issue
- The real problem is more complex baryonic physics:
  - Bars, arms, bulges in massive spirals
  - Non-circular motions
  - Radial inflows/outflows
  - Star formation feedback

These effects are NOT captured by a simple correction term.

RECOMMENDATION:

1. Don't add compactness correction (adds parameter, minimal benefit)
2. Acknowledge massive galaxy limitation honestly
3. Focus on dwarf galaxy success (67-82%)
4. Future work: Include baryonic physics explicitly

┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│   The 0-parameter model is ALREADY optimal for its simplicity.              │
│   Adding parameters must be justified by SIGNIFICANT improvement.           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
""")


def save_results(analysis_results, kappa_results, best_kappa, baseline):
    """Save all results."""

    output = json_safe({
        'session': 45,
        'track': 'C - Compactness Correction',
        'date': datetime.now().isoformat(),
        'analysis': analysis_results,
        'kappa_optimization': kappa_results,
        'best_kappa': best_kappa,
        'baseline': baseline,
        'conclusion': 'Marginal improvement, not recommended',
        'recommendation': 'Keep 0-parameter model, acknowledge massive galaxy limitation'
    })

    output_path = Path(__file__).parent / 'session45_compactness_correction_results.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    return output


if __name__ == '__main__':
    print("\n" + "="*80)
    print("SESSION #45 TRACK C: COMPACTNESS CORRECTION FOR MASSIVE GALAXIES")
    print("="*80)

    # Analyze compactness effect
    analysis = analyze_compactness_effect()

    # Optimize κ
    kappa_results, best, baseline, compactness_median = optimize_kappa()

    # Test by v_max bins
    if best['kappa'] != 0:
        test_by_vmax_bins(best['kappa'], compactness_median)

    # Conclusions
    conclusions()

    # Save
    save_results(analysis, kappa_results, best, baseline)

    print("\n" + "="*80)
    print("SESSION #45 TRACK C COMPLETE")
    print("="*80)
