"""
Analysis for Compatibility-Synthon Experiment (Phase 2)

Tests predictions:
  1. p_crit ∝ 1/<compatibility>
  2. Compatibility structure affects Hill exponent k
  3. Specialist agents show higher specialization index than generalists
  4. Collective cross-type inference > best individual (synthon signature)
  5. Replacement resilience: coherence recovers after agent swap

2026-03-08
"""

import json
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')


# ---------------------------------------------------------------------------
# Curve fitting (same models as Phase 1)
# ---------------------------------------------------------------------------

def hill(p, k, p_half):
    p = np.clip(p, 1e-10, None)
    return p**k / (p**k + p_half**k)

def logistic(p, k, p_half):
    return 1.0 / (1.0 + np.exp(-k * (p - p_half)))

def tanh_model(p, gamma, p_crit):
    p = np.clip(p, 0, None)
    return np.tanh(gamma * np.log(p / p_crit + 1))


def fit_models(p_vals, c_vals):
    """Fit Hill, logistic, tanh. Return best fit params and AIC."""
    results = {}
    p_arr = np.array(p_vals)
    c_arr = np.array(c_vals)

    models = {
        'hill': (hill, [0.6, 0.05]),
        'logistic': (logistic, [20, 0.05]),
        'tanh': (tanh_model, [5, 0.03]),
    }

    for name, (fn, p0) in models.items():
        try:
            popt, _ = curve_fit(fn, p_arr, c_arr, p0=p0, maxfev=5000,
                                 bounds=([0, 1e-6], [20, 1.0]))
            residuals = c_arr - fn(p_arr, *popt)
            rss = float(np.sum(residuals**2))
            n = len(c_arr)
            k_params = len(popt)
            aic = n * np.log(rss / n + 1e-12) + 2 * k_params
            r2 = float(1 - rss / (np.var(c_arr) * n + 1e-12))
            results[name] = {
                'params': popt.tolist(),
                'rss': round(rss, 6),
                'aic': round(aic, 3),
                'r2': round(r2, 4),
            }
        except Exception as e:
            results[name] = {'error': str(e)}

    return results


def aggregate_by_coupling(runs):
    """Aggregate runs by coupling level, return mean C per p."""
    by_p = {}
    for r in runs:
        p = r['coupling']
        c = r['final'].get('C', 0)
        by_p.setdefault(p, []).append(c)
    p_sorted = sorted(by_p.keys())
    p_vals = p_sorted
    c_vals = [np.mean(by_p[p]) for p in p_sorted]
    return p_vals, c_vals


def find_p_crit_empirical(p_vals, c_vals):
    """Find p* = argmax |d²C/dp²|."""
    if len(p_vals) < 4:
        return None
    p_arr = np.array(p_vals)
    c_arr = np.array(c_vals)
    d2 = np.gradient(np.gradient(c_arr, p_arr), p_arr)
    idx = int(np.argmax(np.abs(d2)))
    return float(p_arr[idx])


# ---------------------------------------------------------------------------
# Experiment A: p_crit vs compatibility
# ---------------------------------------------------------------------------

def analyze_experiment_a(data):
    """Test prediction: p_crit ∝ 1/<compatibility>."""
    print("\n=== Experiment A: p_crit ∝ 1/<compatibility> ===")
    results = data['results']

    compat_values = sorted(set(r['mean_compatibility'] for r in results))
    compat_pcrit = {}
    compat_fits = {}

    for compat in compat_values:
        runs = [r for r in results if r['mean_compatibility'] == compat]
        p_vals, c_vals = aggregate_by_coupling(runs)

        fits = fit_models(p_vals, c_vals)
        p_crit_emp = find_p_crit_empirical(p_vals, c_vals)

        # Best p_crit from Hill fit
        p_crit_hill = None
        if 'hill' in fits and 'params' in fits['hill']:
            p_crit_hill = fits['hill']['params'][1]  # p_half

        c_at_max_p = c_vals[-1] if c_vals else None

        print(f"\n  Compatibility = {compat:.2f}:")
        print(f"    C at max coupling: {c_at_max_p:.3f}")
        print(f"    p_crit (Hill p_half): {p_crit_hill:.4f}" if p_crit_hill else "    p_crit: (fit failed)")
        print(f"    p_crit (empirical): {p_crit_emp:.4f}" if p_crit_emp else "    p_crit empirical: N/A")
        for model, res in fits.items():
            if 'r2' in res:
                print(f"    {model}: R²={res['r2']:.4f}, AIC={res['aic']:.1f}")

        compat_pcrit[compat] = {
            'p_crit_hill': p_crit_hill,
            'p_crit_empirical': p_crit_emp,
            'c_max': c_at_max_p,
            'fits': fits,
        }

    # Test p_crit ∝ 1/<C>
    print("\n  Prediction test: p_crit ∝ 1/<compatibility>")
    valid = [(c, v['p_crit_hill']) for c, v in compat_pcrit.items()
             if v['p_crit_hill'] is not None]
    if len(valid) >= 3:
        c_vals_x = np.array([1.0 / c for c, _ in valid])
        pcrit_y = np.array([p for _, p in valid])
        r, pval = pearsonr(c_vals_x, pcrit_y)
        print(f"    Pearson r(1/compat, p_crit) = {r:.3f}, p = {pval:.4f}")
        if abs(r) > 0.8:
            print("    ✓ Strong support for p_crit ∝ 1/<compatibility>")
        elif abs(r) > 0.5:
            print("    ~ Moderate support")
        else:
            print("    ✗ Weak support — p_crit doesn't scale as predicted")

    return compat_pcrit


# ---------------------------------------------------------------------------
# Experiment B: Compatibility structure
# ---------------------------------------------------------------------------

def analyze_experiment_b(data):
    """Test: does compatibility structure affect Hill exponent k?"""
    print("\n=== Experiment B: Compatibility structure vs Hill exponent k ===")
    results = data['results']

    for structure, runs in results.items():
        p_vals, c_vals = aggregate_by_coupling(runs)
        fits = fit_models(p_vals, c_vals)
        hill_k = None
        if 'hill' in fits and 'params' in fits['hill']:
            hill_k = fits['hill']['params'][0]
        print(f"\n  Structure: {structure}")
        print(f"    C at max coupling: {c_vals[-1]:.3f}" if c_vals else "    No data")
        print(f"    Hill k: {hill_k:.3f}" if hill_k else "    Hill k: (fit failed)")
        for model, res in fits.items():
            if 'r2' in res:
                print(f"    {model}: R²={res['r2']:.4f}, AIC={res['aic']:.1f}")


# ---------------------------------------------------------------------------
# Experiment C: Specialist vs generalist
# ---------------------------------------------------------------------------

def analyze_experiment_c(data):
    """Test: specialists → higher specialization index; comparable or better coherence?"""
    print("\n=== Experiment C: Specialist vs Generalist ===")
    results = data['results']

    for assignment, runs in results.items():
        p_vals, c_vals = aggregate_by_coupling(runs)

        # Mean specialization index at high coupling
        high_p_runs = [r for r in runs if r['coupling'] >= 0.1]
        spec_indices = [r['final'].get('spec_idx', 0) for r in high_p_runs]
        mean_spec = np.mean(spec_indices) if spec_indices else 0

        # Mean emergence ratio from synthon metrics
        emergence_ratios = []
        for r in high_p_runs:
            for key, vals in r.get('synthon', {}).get('cross_type_inference', {}).items():
                emergence_ratios.append(vals.get('emergence_ratio', 1.0))
        mean_emergence = np.mean(emergence_ratios) if emergence_ratios else 1.0

        print(f"\n  Assignment: {assignment}")
        print(f"    C at max coupling: {c_vals[-1]:.3f}" if c_vals else "    No data")
        print(f"    Mean specialization index (p≥0.1): {mean_spec:.4f}")
        print(f"    Mean emergence ratio (collective/individual): {mean_emergence:.3f}")
        if mean_emergence > 1.1:
            print(f"    ✓ Collective exceeds best individual — synthon signature present")
        else:
            print(f"    - No clear synthon emergence above individual capability")


# ---------------------------------------------------------------------------
# Experiment D: Replacement resilience
# ---------------------------------------------------------------------------

def analyze_experiment_d(data):
    """Test: does coherence recover after agent replacement?"""
    print("\n=== Experiment D: Replacement Resilience ===")
    results = data['results']

    ratios = [r['recovery_ratio'] for r in results]
    C_before = [r['C_before_replacement'] for r in results]
    C_after = [r['C_after_replacement'] for r in results]

    print(f"  Mean C before replacement: {np.mean(C_before):.3f} ± {np.std(C_before):.3f}")
    print(f"  Mean C after  replacement: {np.mean(C_after):.3f} ± {np.std(C_after):.3f}")
    print(f"  Mean recovery ratio: {np.mean(ratios):.3f} ± {np.std(ratios):.3f}")
    if np.mean(ratios) > 0.85:
        print("  ✓ High resilience — synthon persists through agent replacement")
    elif np.mean(ratios) > 0.6:
        print("  ~ Moderate resilience — partial recovery")
    else:
        print("  ✗ Low resilience — synthon depends on specific agents")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def find_latest_results(prefix):
    """Find most recent result file matching prefix."""
    files = [f for f in os.listdir(RESULTS_DIR) if f.startswith(prefix)]
    if not files:
        return None
    return os.path.join(RESULTS_DIR, sorted(files)[-1])


if __name__ == '__main__':
    print("=== Compatibility-Synthon Experiment: Analysis ===")

    # Experiment A
    fname_a = find_latest_results('compatibility_uniform_')
    if fname_a:
        with open(fname_a) as f:
            data_a = json.load(f)
        print(f"Loaded: {fname_a} ({len(data_a['results'])} runs)")
        a_summary = analyze_experiment_a(data_a)
    else:
        print("Experiment A results not found.")
        a_summary = {}

    # Experiment B
    fname_b = find_latest_results('compatibility_structure_')
    if fname_b:
        with open(fname_b) as f:
            data_b = json.load(f)
        analyze_experiment_b(data_b)
    else:
        print("Experiment B results not found.")

    # Experiment C
    fname_c = find_latest_results('specialist_vs_generalist_')
    if fname_c:
        with open(fname_c) as f:
            data_c = json.load(f)
        analyze_experiment_c(data_c)
    else:
        print("Experiment C results not found.")

    # Experiment D
    fname_d = find_latest_results('replacement_resilience_')
    if fname_d:
        with open(fname_d) as f:
            data_d = json.load(f)
        analyze_experiment_d(data_d)
    else:
        print("Experiment D results not found.")

    print("\n=== Analysis complete. ===")
