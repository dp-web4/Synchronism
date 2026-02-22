"""
Coupling-Coherence Analysis
Fits 4 sigmoid models, compares via AIC/BIC, tests derived p_crit.

Reads: results/coupling_coherence_results.json
       results/coupling_coherence_variations.json (if present)
Writes: results/coupling_coherence_analysis.json
        results/coupling_coherence_plots/*.png
"""

import numpy as np
import json
import os
import sys
from scipy.optimize import curve_fit
from scipy.special import erf as scipy_erf

# ── Model definitions ──────────────────────────────────────────────

def model_tanh(p, gamma, p_crit):
    """C = tanh(γ · log(p/p_crit + 1))"""
    return np.tanh(gamma * np.log(p / p_crit + 1.0))

def model_logistic(p, k, p_half):
    """C = 1 / (1 + exp(-k·(p - p_half)))"""
    return 1.0 / (1.0 + np.exp(-k * (p - p_half)))

def model_erf(p, k, p_half):
    """C = erf(k · (p - p_half))"""
    return scipy_erf(k * (p - p_half))

def model_hill(p, k, p_half):
    """C = p^k / (p^k + p_half^k)  (Hill function)"""
    # Avoid 0^k issues
    p_safe = np.maximum(p, 1e-10)
    p_half_safe = max(p_half, 1e-10)
    return p_safe**k / (p_safe**k + p_half_safe**k)

MODELS = {
    'tanh': {
        'func': model_tanh,
        'p0': [2.0, 0.02],
        'bounds': ([0.01, 1e-4], [20.0, 1.0]),
        'param_names': ['gamma', 'p_crit'],
    },
    'logistic': {
        'func': model_logistic,
        'p0': [20.0, 0.05],
        'bounds': ([0.1, -0.5], [200.0, 1.0]),
        'param_names': ['k', 'p_half'],
    },
    'erf': {
        'func': model_erf,
        'p0': [5.0, 0.02],
        'bounds': ([0.1, -1.0], [100.0, 1.0]),
        'param_names': ['k', 'p_half'],
    },
    'hill': {
        'func': model_hill,
        'p0': [1.0, 0.02],
        'bounds': ([0.01, 1e-4], [10.0, 1.0]),
        'param_names': ['k', 'p_half'],
    },
}


# ── Data loading / aggregation ─────────────────────────────────────

def load_primary_results(path):
    """Load and aggregate primary experiment results."""
    with open(path) as f:
        data = json.load(f)

    # Group by coupling level, take mean of final C values
    from collections import defaultdict
    by_coupling = defaultdict(list)
    for r in data['results']:
        p = r['coupling']
        by_coupling[p].append(r['final'])

    couplings = sorted(by_coupling.keys())
    mean_C = []
    std_C = []
    mean_conv = []
    mean_corr = []
    for p in couplings:
        finals = by_coupling[p]
        cs = [f['C'] for f in finals]
        convs = [f['C_conv'] for f in finals]
        corrs = [f['C_corr'] for f in finals]
        mean_C.append(np.mean(cs))
        std_C.append(np.std(cs))
        mean_conv.append(np.mean(convs))
        mean_corr.append(np.mean(corrs))

    return {
        'couplings': np.array(couplings),
        'mean_C': np.array(mean_C),
        'std_C': np.array(std_C),
        'mean_conv': np.array(mean_conv),
        'mean_corr': np.array(mean_corr),
        'params': data['params'],
        'raw': data,
    }


# ── Curve fitting ──────────────────────────────────────────────────

def fit_all_models(couplings, mean_C, std_C):
    """Fit all 4 models, return results dict."""
    n = len(couplings)
    results = {}

    for name, spec in MODELS.items():
        try:
            # Use std as sigma weights (avoid zero weights)
            sigma = np.maximum(std_C, 1e-6)
            popt, pcov = curve_fit(
                spec['func'], couplings, mean_C,
                p0=spec['p0'],
                bounds=spec['bounds'],
                sigma=sigma,
                maxfev=10000,
            )
            fitted = spec['func'](couplings, *popt)
            residuals = mean_C - fitted
            rss = float(np.sum(residuals**2))
            ss_tot = float(np.sum((mean_C - np.mean(mean_C))**2))
            r_squared = 1.0 - rss / ss_tot if ss_tot > 0 else 0.0

            # AIC and BIC (2 parameters for all models)
            k_params = 2
            # Log-likelihood assuming Gaussian errors
            log_lik = -n/2 * np.log(2 * np.pi * rss / n) - n/2
            aic = 2 * k_params - 2 * log_lik
            bic = k_params * np.log(n) - 2 * log_lik

            results[name] = {
                'params': {spec['param_names'][i]: round(float(popt[i]), 6) for i in range(len(popt))},
                'rss': round(rss, 6),
                'r_squared': round(r_squared, 6),
                'aic': round(aic, 4),
                'bic': round(bic, 4),
                'fitted_values': fitted.tolist(),
                'success': True,
            }
        except Exception as e:
            results[name] = {
                'error': str(e),
                'success': False,
            }

    return results


# ── Change-point detection ─────────────────────────────────────────

def find_changepoint(couplings, mean_C):
    """Find point of maximum curvature change (numerical 2nd derivative)."""
    if len(couplings) < 5:
        return None

    # Smooth slightly for stable derivatives
    from scipy.ndimage import uniform_filter1d
    smoothed = uniform_filter1d(mean_C, size=3)

    # Numerical second derivative
    dp = np.diff(couplings)
    dC = np.diff(smoothed)
    first_deriv = dC / dp

    dp2 = 0.5 * (dp[:-1] + dp[1:])
    d2C = np.diff(first_deriv) / dp2
    second_deriv = np.abs(d2C)

    # argmax of |d²C/dp²| — the point indices are offset by 1
    idx = np.argmax(second_deriv)
    p_star = float(couplings[idx + 1])  # +1 because second derivative loses one point from each end

    return {
        'p_star': round(p_star, 4),
        'max_curvature': round(float(second_deriv[idx]), 6),
        'second_derivative': [round(float(x), 6) for x in d2C],
    }


# ── Derived p_crit ────────────────────────────────────────────────

def derived_p_crit(eta, H_world, K, m):
    """
    p_crit_derived = η · H(world) / (K · m · (1 - 2η))

    Where:
      η = noise rate
      H(world) = world entropy in bits
      K = number of agents
      m = observations per round per agent
    """
    denominator = K * m * (1.0 - 2.0 * eta)
    if denominator <= 0:
        return float('inf')
    return eta * H_world / denominator


def world_entropy_from_params(n_nodes, n_edges, n_types):
    """Compute world entropy from parameters (matching World.entropy())."""
    total_slots = n_nodes * (n_nodes - 1) * n_types
    p_edge = n_edges / total_slots
    if p_edge <= 0 or p_edge >= 1:
        return 0.0
    return total_slots * (-p_edge * np.log2(p_edge) - (1 - p_edge) * np.log2(1 - p_edge))


# ── Variation analysis ─────────────────────────────────────────────

def analyze_variations(var_path):
    """Analyze variation experiments for derived vs fitted p_crit."""
    if not os.path.exists(var_path):
        return None

    with open(var_path) as f:
        data = json.load(f)

    from collections import defaultdict

    variation_results = []
    for var in data['variations']:
        label = var['label']
        params = var['params']

        # Aggregate by coupling
        by_coupling = defaultdict(list)
        for r in var['results']:
            by_coupling[r['coupling']].append(r['final']['C'])

        couplings = np.array(sorted(by_coupling.keys()))
        mean_C = np.array([np.mean(by_coupling[p]) for p in couplings])
        std_C = np.array([np.std(by_coupling[p]) for p in couplings])

        # Fit tanh model
        try:
            sigma = np.maximum(std_C, 1e-6)
            popt, _ = curve_fit(
                model_tanh, couplings, mean_C,
                p0=[2.0, 0.02],
                bounds=([0.01, 1e-4], [20.0, 1.0]),
                sigma=sigma,
                maxfev=10000,
            )
            fitted_p_crit = float(popt[1])
            fitted_gamma = float(popt[0])
        except Exception:
            fitted_p_crit = None
            fitted_gamma = None

        # Compute derived p_crit
        H = world_entropy_from_params(
            params['n_nodes'], params['n_edges'], params['n_types'])
        d_pcrit = derived_p_crit(
            params['noise_rate'], H,
            params['n_agents'], params['obs_per_round'])

        variation_results.append({
            'label': label,
            'params': params,
            'H_world': round(H, 2),
            'derived_p_crit': round(d_pcrit, 4) if d_pcrit != float('inf') else None,
            'fitted_p_crit': round(fitted_p_crit, 4) if fitted_p_crit else None,
            'fitted_gamma': round(fitted_gamma, 4) if fitted_gamma else None,
            'mean_C_at_couplings': {str(round(p, 2)): round(float(c), 4)
                                     for p, c in zip(couplings, mean_C)},
        })

    # Compute R² of derived vs fitted p_crit
    fitted_vals = [v['fitted_p_crit'] for v in variation_results if v['fitted_p_crit'] is not None and v['derived_p_crit'] is not None]
    derived_vals = [v['derived_p_crit'] for v in variation_results if v['fitted_p_crit'] is not None and v['derived_p_crit'] is not None]

    derivation_r2 = None
    if len(fitted_vals) >= 3:
        fitted_arr = np.array(fitted_vals)
        derived_arr = np.array(derived_vals)
        ss_res = np.sum((fitted_arr - derived_arr)**2)
        ss_tot = np.sum((fitted_arr - np.mean(fitted_arr))**2)
        if ss_tot > 0:
            derivation_r2 = round(float(1.0 - ss_res / ss_tot), 4)

    return {
        'variations': variation_results,
        'derivation_r2': derivation_r2,
        'n_valid_comparisons': len(fitted_vals),
    }


# ── Plot generation ────────────────────────────────────────────────

def generate_plots(couplings, mean_C, std_C, mean_conv, mean_corr,
                   fit_results, changepoint, var_analysis, plot_dir):
    """Generate PNG plots using matplotlib."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot generation")
        return

    os.makedirs(plot_dir, exist_ok=True)

    # ── Plot 1: C(p) with all model fits ──
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(couplings, mean_C, yerr=std_C, fmt='o', markersize=3,
                color='black', alpha=0.5, label='Data (mean ± std)', capsize=2)

    colors = {'tanh': '#e74c3c', 'logistic': '#3498db', 'erf': '#2ecc71', 'hill': '#9b59b6'}
    p_fine = np.linspace(0, 1, 200)
    for name, spec in MODELS.items():
        if name in fit_results and fit_results[name]['success']:
            params = fit_results[name]['params']
            pvals = [params[k] for k in spec['param_names']]
            y_fit = spec['func'](p_fine, *pvals)
            r2 = fit_results[name]['r_squared']
            ax.plot(p_fine, y_fit, color=colors[name], linewidth=2,
                    label=f'{name} (R²={r2:.4f})')

    if changepoint:
        ax.axvline(changepoint['p_star'], color='gray', linestyle='--', alpha=0.5,
                   label=f'Changepoint p*={changepoint["p_star"]:.2f}')

    if 'tanh' in fit_results and fit_results['tanh']['success']:
        p_crit = fit_results['tanh']['params']['p_crit']
        ax.axvline(p_crit, color='#e74c3c', linestyle=':', alpha=0.5,
                   label=f'tanh p_crit={p_crit:.3f}')

    ax.set_xlabel('Coupling (p)', fontsize=12)
    ax.set_ylabel('Coherence C(p)', fontsize=12)
    ax.set_title('Coherence vs Coupling: 4-Model Comparison', fontsize=14)
    ax.legend(fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, 'c_vs_p.png'), dpi=150)
    plt.close(fig)

    # ── Plot 2: Convergence vs Correctness ──
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(mean_conv, mean_corr, c=couplings, cmap='viridis',
                         s=40, edgecolors='black', linewidths=0.5)
    plt.colorbar(scatter, label='Coupling (p)')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='C_conv = C_corr')
    ax.set_xlabel('Convergence (C_conv)', fontsize=12)
    ax.set_ylabel('Correctness (C_corr)', fontsize=12)
    ax.set_title('Convergence vs Correctness\n(detecting "shared wrongness")', fontsize=14)
    ax.legend()
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, 'convergence_vs_correctness.png'), dpi=150)
    plt.close(fig)

    # ── Plot 3: Model comparison (AIC/BIC) ──
    successful = {n: r for n, r in fit_results.items() if r['success']}
    if successful:
        fig, axes = plt.subplots(1, 3, figsize=(14, 5))

        names = list(successful.keys())
        r2_vals = [successful[n]['r_squared'] for n in names]
        aic_vals = [successful[n]['aic'] for n in names]
        bic_vals = [successful[n]['bic'] for n in names]

        bar_colors = [colors.get(n, 'gray') for n in names]

        axes[0].bar(names, r2_vals, color=bar_colors)
        axes[0].set_ylabel('R²')
        axes[0].set_title('Goodness of Fit (R²)')
        axes[0].set_ylim(min(r2_vals) - 0.01, 1.0)

        axes[1].bar(names, aic_vals, color=bar_colors)
        axes[1].set_ylabel('AIC')
        axes[1].set_title('AIC (lower = better)')

        axes[2].bar(names, bic_vals, color=bar_colors)
        axes[2].set_ylabel('BIC')
        axes[2].set_title('BIC (lower = better)')

        fig.suptitle('Model Comparison', fontsize=14)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'model_comparison.png'), dpi=150)
        plt.close(fig)

    # ── Plot 4: Changepoint ──
    if changepoint and changepoint.get('second_derivative'):
        d2 = changepoint['second_derivative']
        p_d2 = couplings[1:-1][:len(d2)]  # Align indices
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax1.plot(couplings, mean_C, 'ko-', markersize=3, label='C(p)')
        ax1.axvline(changepoint['p_star'], color='red', linestyle='--',
                    label=f'p* = {changepoint["p_star"]:.3f}')
        ax1.set_ylabel('C(p)')
        ax1.set_title('Coherence and Change-Point Detection')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2.plot(p_d2, d2, 'b-', linewidth=1.5)
        ax2.axvline(changepoint['p_star'], color='red', linestyle='--')
        ax2.set_xlabel('Coupling (p)')
        ax2.set_ylabel('d²C/dp²')
        ax2.set_title('Second Derivative (curvature)')
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'changepoint.png'), dpi=150)
        plt.close(fig)

    # ── Plot 5: Derived vs Fitted p_crit ──
    if var_analysis and var_analysis['n_valid_comparisons'] >= 2:
        fig, ax = plt.subplots(figsize=(7, 7))
        derived = []
        fitted = []
        labels = []
        for v in var_analysis['variations']:
            if v['fitted_p_crit'] is not None and v['derived_p_crit'] is not None:
                derived.append(v['derived_p_crit'])
                fitted.append(v['fitted_p_crit'])
                labels.append(v['label'])

        ax.scatter(derived, fitted, s=60, c='#e74c3c', edgecolors='black', zorder=5)
        for i, lbl in enumerate(labels):
            ax.annotate(lbl, (derived[i], fitted[i]), textcoords='offset points',
                        xytext=(5, 5), fontsize=8)

        # Identity line
        all_vals = derived + fitted
        lo, hi = min(all_vals) * 0.8, max(all_vals) * 1.2
        ax.plot([lo, hi], [lo, hi], 'k--', alpha=0.3, label='Perfect derivation')

        r2_text = f'R² = {var_analysis["derivation_r2"]:.3f}' if var_analysis['derivation_r2'] is not None else 'R² = N/A'
        ax.set_xlabel('Derived p_crit', fontsize=12)
        ax.set_ylabel('Fitted p_crit', fontsize=12)
        ax.set_title(f'Derived vs Fitted p_crit ({r2_text})', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, 'pcrit_derived_vs_fitted.png'), dpi=150)
        plt.close(fig)

    # ── Plot 6: Convergence + Correctness + Coherence ──
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(couplings, mean_conv, 'b-', linewidth=2, label='C_conv (convergence)')
    ax.plot(couplings, mean_corr, 'g-', linewidth=2, label='C_corr (correctness)')
    ax.plot(couplings, mean_C, 'r-', linewidth=2.5, label='C (coherence = √(conv × corr))')
    ax.fill_between(couplings, mean_C - std_C, mean_C + std_C, color='red', alpha=0.1)
    ax.set_xlabel('Coupling (p)', fontsize=12)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title('Three Metrics: Convergence, Correctness, and Coherence', fontsize=14)
    ax.legend(fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, 'three_metrics.png'), dpi=150)
    plt.close(fig)

    print(f"Plots saved to {plot_dir}/")


# ── Main ───────────────────────────────────────────────────────────

def run_analysis():
    """Run the full analysis pipeline."""
    base_dir = os.path.join(os.path.dirname(__file__), 'results')
    primary_path = os.path.join(base_dir, 'coupling_coherence_results.json')
    var_path = os.path.join(base_dir, 'coupling_coherence_variations.json')
    plot_dir = os.path.join(base_dir, 'coupling_coherence_plots')

    if not os.path.exists(primary_path):
        print(f"ERROR: Primary results not found at {primary_path}")
        print("Run coupling_coherence_experiment.py first.")
        sys.exit(1)

    print("=" * 60)
    print("COUPLING-COHERENCE ANALYSIS")
    print("=" * 60)

    # 1. Load data
    print("\n1. Loading primary results...")
    data = load_primary_results(primary_path)
    print(f"   {len(data['couplings'])} coupling levels, "
          f"C range: [{data['mean_C'].min():.3f}, {data['mean_C'].max():.3f}]")

    # 2. Fit models
    print("\n2. Fitting 4 sigmoid models...")
    fit_results = fit_all_models(data['couplings'], data['mean_C'], data['std_C'])
    for name, res in fit_results.items():
        if res['success']:
            print(f"   {name:10s}: R²={res['r_squared']:.4f}  AIC={res['aic']:.1f}  "
                  f"params={res['params']}")
        else:
            print(f"   {name:10s}: FAILED — {res.get('error', 'unknown')}")

    # Determine winner
    successful = {n: r for n, r in fit_results.items() if r['success']}
    if successful:
        best_aic = min(successful, key=lambda n: successful[n]['aic'])
        best_r2 = max(successful, key=lambda n: successful[n]['r_squared'])
        print(f"\n   Best by AIC: {best_aic}")
        print(f"   Best by R²:  {best_r2}")

        # Check if tanh wins
        if 'tanh' in successful:
            tanh_aic = successful['tanh']['aic']
            for name, res in successful.items():
                if name != 'tanh':
                    delta = tanh_aic - res['aic']
                    verdict = "tanh preferred" if delta < -2 else ("indistinguishable" if abs(delta) <= 2 else "tanh disfavored")
                    print(f"   ΔAIC(tanh - {name}) = {delta:.1f} → {verdict}")

    # 3. Change-point detection
    print("\n3. Change-point detection...")
    changepoint = find_changepoint(data['couplings'], data['mean_C'])
    if changepoint:
        print(f"   p* (max curvature) = {changepoint['p_star']}")
        if 'tanh' in successful:
            p_crit = successful['tanh']['params']['p_crit']
            rel_diff = abs(changepoint['p_star'] - p_crit) / p_crit if p_crit > 0 else float('inf')
            consistent = rel_diff < 0.15
            print(f"   tanh p_crit = {p_crit:.4f}")
            print(f"   |p* - p_crit| / p_crit = {rel_diff:.3f} → "
                  f"{'CONSISTENT' if consistent else 'INCONSISTENT'} (threshold: 0.15)")

    # 4. Derived p_crit for primary experiment
    print("\n4. Derived p_crit (primary experiment)...")
    H = world_entropy_from_params(
        data['params']['n_nodes'], data['params']['n_edges'], data['params']['n_types'])
    d_pcrit = derived_p_crit(
        data['params']['noise_rate'], H,
        data['params']['n_agents'], data['params']['obs_per_round'])
    print(f"   H(world) = {H:.2f} bits")
    print(f"   Derived p_crit = {d_pcrit:.4f}")
    if 'tanh' in successful:
        f_pcrit = successful['tanh']['params']['p_crit']
        print(f"   Fitted  p_crit = {f_pcrit:.4f}")
        print(f"   Ratio (derived/fitted) = {d_pcrit/f_pcrit:.3f}")

    # 5. Variation analysis
    var_analysis = None
    if os.path.exists(var_path):
        print("\n5. Analyzing variations (derived vs fitted p_crit)...")
        var_analysis = analyze_variations(var_path)
        if var_analysis:
            print(f"   Valid comparisons: {var_analysis['n_valid_comparisons']}")
            if var_analysis['derivation_r2'] is not None:
                r2 = var_analysis['derivation_r2']
                verdict = "EXPLANATORY" if r2 > 0.8 else ("MODERATE" if r2 > 0.5 else "WEAK")
                print(f"   Derivation R² = {r2:.4f} → {verdict}")
            for v in var_analysis['variations']:
                print(f"   {v['label']:10s}: derived={v['derived_p_crit']}  "
                      f"fitted={v['fitted_p_crit']}  γ={v['fitted_gamma']}")
    else:
        print("\n5. No variation results found (run experiment with --variations)")

    # 6. Generate plots
    print("\n6. Generating plots...")
    generate_plots(data['couplings'], data['mean_C'], data['std_C'],
                   data['mean_conv'], data['mean_corr'],
                   fit_results, changepoint, var_analysis, plot_dir)

    # 7. Summary and kill criteria assessment
    print("\n" + "=" * 60)
    print("KILL CRITERIA ASSESSMENT")
    print("=" * 60)

    if successful:
        # Kill criterion 1: Does logistic/erf beat tanh?
        if 'tanh' in successful:
            beaten = any(successful[n]['aic'] < successful['tanh']['aic'] - 2
                         for n in successful if n != 'tanh')
            print(f"\n  1. Logistic/erf beats tanh by ΔAIC > 2?  {'YES — tanh NOT preferred' if beaten else 'NO — tanh holds'}")

        # Kill criterion 2: p_crit derivation
        if var_analysis and var_analysis['derivation_r2'] is not None:
            weak = var_analysis['derivation_r2'] < 0.5
            print(f"  2. p_crit derivation R² < 0.5?  {'YES — p_crit is only a fit parameter' if weak else 'NO — derivation has explanatory power'}")
        else:
            print("  2. p_crit derivation: NOT YET TESTED (need variations)")

        # Kill criterion 3: Shared wrongness
        max_gap = max(data['mean_conv'] - data['mean_corr'])
        shared_wrongness = max_gap > 0.3
        print(f"  3. C_conv >> C_corr (shared wrongness)?  {'YES — coherence ≠ correctness' if shared_wrongness else 'NO — convergence tracks correctness'}")
        print(f"     Max gap (C_conv - C_corr) = {max_gap:.3f}")

        # Kill criterion 4: No sigmoid transition
        c_range = data['mean_C'].max() - data['mean_C'].min()
        linear = c_range < 0.3
        print(f"  4. No sigmoid transition (linear)?  {'YES — no phase transition' if linear else 'NO — clear sigmoid observed'}")
        print(f"     C range = {c_range:.3f}")

    # 8. Save analysis
    output = {
        'experiment': 'coupling_coherence_analysis',
        'params': data['params'],
        'fit_results': fit_results,
        'changepoint': changepoint,
        'primary_derived_p_crit': round(d_pcrit, 4),
        'primary_H_world': round(H, 2),
        'variation_analysis': var_analysis,
        'kill_criteria': {
            'tanh_beaten': beaten if 'tanh' in successful else None,
            'pcrit_derivation_r2': var_analysis['derivation_r2'] if var_analysis else None,
            'max_conv_corr_gap': round(float(max_gap), 4),
            'c_range': round(float(c_range), 4),
        } if successful else None,
        'data_summary': {
            'couplings': data['couplings'].tolist(),
            'mean_C': [round(float(x), 4) for x in data['mean_C']],
            'std_C': [round(float(x), 4) for x in data['std_C']],
            'mean_conv': [round(float(x), 4) for x in data['mean_conv']],
            'mean_corr': [round(float(x), 4) for x in data['mean_corr']],
        },
    }

    out_path = os.path.join(base_dir, 'coupling_coherence_analysis.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nAnalysis saved to {out_path}")


if __name__ == '__main__':
    run_analysis()
