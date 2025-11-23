#!/usr/bin/env python3
"""
Session #41: Physical ρ_crit Correlation Analysis

Nova's Recommendation (Session #40 Review):
"The physical interpretation of ρ_crit could be enriched by correlating it
with galaxy properties such as mass, velocity dispersion, and morphology."

Research Question:
What does ρ_crit represent physically?

Hypotheses to Test:
1. Decoherence Scale: ρ_crit ∝ σ_v² (velocity dispersion squared)
2. Mass Scaling: ρ_crit ∝ M_total^α (galaxy mass)
3. Morphology Dependence: Different types have different ρ_crit distributions
4. Density Threshold: ρ_crit ∝ ρ_vis,max (peak visible density)

Goal:
Determine which hypothesis best explains observed ρ_crit values.

Author: CBP Autonomous Synchronism Research
Date: 2025-11-23
Session: #41 - Physical Interpretation Track
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit

sys.path.append(str(Path(__file__).parent))
from synchronism_real_sparc_validation import RealSPARCLoader


def extract_galaxy_properties(galaxy):
    """
    Extract physical properties from SPARC galaxy data.

    Returns dict with:
    - total_luminosity: Total luminosity [L_sun]
    - distance: Distance [Mpc]
    - inclination: Inclination angle [deg]
    - morphology: Morphology type string
    - v_max: Maximum rotation velocity [km/s]
    - v_mean: Mean rotation velocity [km/s]
    - r_max: Maximum radius [kpc]
    - rho_vis_max: Peak visible density [M_sun/pc²]
    - rho_vis_mean: Mean visible density [M_sun/pc²]
    """
    r = galaxy.radius
    v_obs = galaxy.v_obs
    rho_vis = galaxy.total_baryonic_density()

    return {
        'name': galaxy.name,
        'total_luminosity': galaxy.total_luminosity,
        'distance': galaxy.distance,
        'inclination': galaxy.inclination,
        'morphology': galaxy.morphology,
        'v_max': float(np.max(v_obs)),
        'v_mean': float(np.mean(v_obs)),
        'r_max': float(np.max(r)),
        'rho_vis_max': float(np.max(rho_vis)),
        'rho_vis_mean': float(np.mean(rho_vis)),
        # Estimate velocity dispersion (proxy from rotation curve width)
        'sigma_v_estimate': float(np.std(v_obs))
    }


def load_rho_crit_results(results_file='session41_continuous_results.json'):
    """
    Load ρ_crit values from Session #41 continuous optimization.

    If not available, fall back to Session #40 extended grid.
    """
    # Try Session #41 first
    s41_path = Path(__file__).parent / results_file

    if s41_path.exists():
        with open(s41_path) as f:
            data = json.load(f)
        rho_crit_dict = {r['name']: r['rho_crit'] for r in data['results']}
        source = 'Session #41 (continuous)'
    else:
        # Fall back to Session #40
        s40_path = Path(__file__).parent / 'session40_extended_summary.json'
        if s40_path.exists():
            with open(s40_path) as f:
                data = json.load(f)
            # Session #40 only has 84 galaxies - we'll use what's available
            # For full analysis, wait for Session #41 to complete
            print("Warning: Using Session #40 data (84 galaxies only)")
            print("Rerun after Session #41 completes for full 175-galaxy analysis")
            return None, 'Session #40 (incomplete)'
        else:
            raise FileNotFoundError("No ρ_crit results found. Run session41_continuous_rho_crit.py first.")

    return rho_crit_dict, source


def correlate_rho_crit_with_properties(rho_crit_dict):
    """
    Correlate ρ_crit with galaxy physical properties.
    """
    # Load all SPARC galaxies
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()

    # Extract properties and ρ_crit for each galaxy
    data = []
    for galaxy in galaxies:
        if galaxy.name in rho_crit_dict:
            props = extract_galaxy_properties(galaxy)
            props['rho_crit'] = rho_crit_dict[galaxy.name]
            data.append(props)

    print(f"\nLoaded {len(data)} galaxies with both properties and ρ_crit")

    # Convert to arrays for analysis
    rho_crits = np.array([d['rho_crit'] for d in data])
    v_maxs = np.array([d['v_max'] for d in data])
    v_means = np.array([d['v_mean'] for d in data])
    sigma_vs = np.array([d['sigma_v_estimate'] for d in data])
    rho_vis_maxs = np.array([d['rho_vis_max'] for d in data])
    rho_vis_means = np.array([d['rho_vis_mean'] for d in data])
    luminosities = np.array([d['total_luminosity'] for d in data])

    # Test correlations
    correlations = {}

    # Hypothesis 1: Decoherence scale (ρ_crit ∝ σ_v²)
    correlations['sigma_v²'] = test_correlation(
        rho_crits, sigma_vs**2,
        'ρ_crit', 'σ_v²',
        log_x=True, log_y=True
    )

    # Hypothesis 2: Mass/Luminosity scaling
    correlations['luminosity'] = test_correlation(
        rho_crits, luminosities,
        'ρ_crit', 'L_total',
        log_x=True, log_y=True
    )

    # Hypothesis 3: Maximum velocity (proxy for mass)
    correlations['v_max'] = test_correlation(
        rho_crits, v_maxs,
        'ρ_crit', 'v_max',
        log_x=True, log_y=False
    )

    # Hypothesis 4: Peak visible density
    correlations['rho_vis_max'] = test_correlation(
        rho_crits, rho_vis_maxs,
        'ρ_crit', 'ρ_vis,max',
        log_x=True, log_y=True
    )

    # Hypothesis 5: Mean visible density
    correlations['rho_vis_mean'] = test_correlation(
        rho_crits, rho_vis_means,
        'ρ_crit', 'ρ_vis,mean',
        log_x=True, log_y=True
    )

    return data, correlations


def test_correlation(x, y, x_label, y_label, log_x=False, log_y=False):
    """
    Test correlation between two variables.

    Returns:
        dict with Pearson r, Spearman ρ, p-values, and best-fit power law
    """
    # Remove NaN/inf
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 10:
        return {'error': 'Insufficient data'}

    # Pearson correlation (linear)
    r_pearson, p_pearson = pearsonr(x_clean, y_clean)

    # Spearman correlation (rank-based, robust to outliers)
    rho_spearman, p_spearman = spearmanr(x_clean, y_clean)

    # Power law fit: y = A * x^B (in log space: log y = log A + B log x)
    if log_x and log_y:
        log_x_clean = np.log10(x_clean)
        log_y_clean = np.log10(y_clean)

        # Linear fit in log-log space
        coeffs = np.polyfit(log_x_clean, log_y_clean, 1)
        B = coeffs[0]  # Exponent
        log_A = coeffs[1]  # Log of prefactor
        A = 10**log_A

        # R² in log-log space
        y_pred_log = coeffs[0] * log_x_clean + coeffs[1]
        ss_res = np.sum((log_y_clean - y_pred_log)**2)
        ss_tot = np.sum((log_y_clean - np.mean(log_y_clean))**2)
        r_squared = 1 - (ss_res / ss_tot)

        power_law = {
            'A': A,
            'B': B,
            'formula': f'{y_label} = {A:.2e} × {x_label}^{B:.3f}',
            'r_squared': r_squared
        }
    else:
        power_law = None

    return {
        'n': len(x_clean),
        'pearson_r': r_pearson,
        'pearson_p': p_pearson,
        'spearman_rho': rho_spearman,
        'spearman_p': p_spearman,
        'power_law': power_law
    }


def analyze_morphology_dependence(data):
    """
    Analyze ρ_crit distribution by morphology type.
    """
    # Group by morphology
    morphology_groups = {}
    for d in data:
        morph = d['morphology']
        if morph not in morphology_groups:
            morphology_groups[morph] = []
        morphology_groups[morph].append(d['rho_crit'])

    # Statistics for each morphology
    morph_stats = {}
    for morph, rho_crits in morphology_groups.items():
        if len(rho_crits) >= 3:  # Need at least 3 for meaningful stats
            morph_stats[morph] = {
                'n': len(rho_crits),
                'median': np.median(rho_crits),
                'mean': np.mean(rho_crits),
                'std': np.std(rho_crits),
                'min': np.min(rho_crits),
                'max': np.max(rho_crits)
            }

    return morph_stats


def plot_correlations(data, correlations, output_dir='session41_analysis'):
    """
    Visualize ρ_crit correlations with galaxy properties.
    """
    output_path = Path(__file__).parent / output_dir
    output_path.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    # Extract data
    rho_crits = np.array([d['rho_crit'] for d in data])
    sigma_vs = np.array([d['sigma_v_estimate'] for d in data])
    v_maxs = np.array([d['v_max'] for d in data])
    rho_vis_maxs = np.array([d['rho_vis_max'] for d in data])
    rho_vis_means = np.array([d['rho_vis_mean'] for d in data])
    luminosities = np.array([d['total_luminosity'] for d in data])

    # Plot 1: ρ_crit vs σ_v² (Decoherence Hypothesis)
    ax = axes[0]
    if 'sigma_v²' in correlations and correlations['sigma_v²'].get('power_law'):
        mask = (sigma_vs > 0) & (rho_crits > 0)
        ax.scatter(sigma_vs[mask]**2, rho_crits[mask], alpha=0.5, s=30)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('σ_v² [(km/s)²]', fontsize=11)
        ax.set_ylabel('ρ_crit [M☉/pc²]', fontsize=11)

        pl = correlations['sigma_v²']['power_law']
        ax.set_title(f'Decoherence Hypothesis\nρ ∝ {pl["formula"]}\n'
                    f'r = {correlations["sigma_v²"]["spearman_rho"]:.3f}',
                    fontsize=10)
        ax.grid(True, alpha=0.3)

    # Plot 2: ρ_crit vs v_max
    ax = axes[1]
    if 'v_max' in correlations:
        mask = (v_maxs > 0) & (rho_crits > 0)
        ax.scatter(v_maxs[mask], rho_crits[mask], alpha=0.5, s=30)
        ax.set_xscale('linear')
        ax.set_yscale('log')
        ax.set_xlabel('v_max [km/s]', fontsize=11)
        ax.set_ylabel('ρ_crit [M☉/pc²]', fontsize=11)
        ax.set_title(f'Maximum Velocity\n'
                    f'ρ = {correlations["v_max"]["spearman_rho"]:.3f}',
                    fontsize=10)
        ax.grid(True, alpha=0.3)

    # Plot 3: ρ_crit vs ρ_vis,max
    ax = axes[2]
    if 'rho_vis_max' in correlations and correlations['rho_vis_max'].get('power_law'):
        mask = (rho_vis_maxs > 0) & (rho_crits > 0)
        ax.scatter(rho_vis_maxs[mask], rho_crits[mask], alpha=0.5, s=30)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('ρ_vis,max [M☉/pc²]', fontsize=11)
        ax.set_ylabel('ρ_crit [M☉/pc²]', fontsize=11)

        pl = correlations['rho_vis_max']['power_law']
        ax.set_title(f'Peak Density Scaling\n{pl["formula"]}\n'
                    f'ρ = {correlations["rho_vis_max"]["spearman_rho"]:.3f}',
                    fontsize=10)
        ax.grid(True, alpha=0.3)

    # Plot 4: ρ_crit vs luminosity
    ax = axes[3]
    if 'luminosity' in correlations and correlations['luminosity'].get('power_law'):
        mask = (luminosities > 0) & (rho_crits > 0)
        ax.scatter(luminosities[mask], rho_crits[mask], alpha=0.5, s=30)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Total Luminosity [L☉]', fontsize=11)
        ax.set_ylabel('ρ_crit [M☉/pc²]', fontsize=11)

        pl = correlations['luminosity']['power_law']
        ax.set_title(f'Mass/Luminosity Scaling\n{pl["formula"]}\n'
                    f'ρ = {correlations["luminosity"]["spearman_rho"]:.3f}',
                    fontsize=10)
        ax.grid(True, alpha=0.3)

    # Plot 5: ρ_crit vs ρ_vis,mean
    ax = axes[4]
    if 'rho_vis_mean' in correlations and correlations['rho_vis_mean'].get('power_law'):
        mask = (rho_vis_means > 0) & (rho_crits > 0)
        ax.scatter(rho_vis_means[mask], rho_crits[mask], alpha=0.5, s=30)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('ρ_vis,mean [M☉/pc²]', fontsize=11)
        ax.set_ylabel('ρ_crit [M☉/pc²]', fontsize=11)

        pl = correlations['rho_vis_mean']['power_law']
        ax.set_title(f'Mean Density Scaling\n{pl["formula"]}\n'
                    f'ρ = {correlations["rho_vis_mean"]["spearman_rho"]:.3f}',
                    fontsize=10)
        ax.grid(True, alpha=0.3)

    # Plot 6: Summary text
    ax = axes[5]
    ax.axis('off')
    summary_text = "Correlation Summary\n\n"
    for key, corr in correlations.items():
        if 'error' not in corr:
            summary_text += f"{key}:\n"
            summary_text += f"  Spearman ρ = {corr['spearman_rho']:.3f}\n"
            summary_text += f"  p-value = {corr['spearman_p']:.2e}\n\n"
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
           fontsize=9, verticalalignment='top', family='monospace')

    fig.suptitle('Physical ρ_crit Correlations - Session #41',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_file = output_path / 'rho_crit_physical_correlations.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return output_file


def main():
    print("\n" + "="*80)
    print("SESSION #41: PHYSICAL ρ_crit CORRELATION ANALYSIS")
    print("="*80)

    # Load ρ_crit results
    print("\nLoading ρ_crit values...")
    rho_crit_dict, source = load_rho_crit_results()

    if rho_crit_dict is None:
        print("\nWaiting for Session #41 continuous optimization to complete...")
        print("Rerun this script after session41_continuous_rho_crit.py finishes.")
        return

    print(f"Source: {source}")
    print(f"Loaded ρ_crit for {len(rho_crit_dict)} galaxies")

    # Correlate with properties
    print("\nTesting physical correlations...")
    data, correlations = correlate_rho_crit_with_properties(rho_crit_dict)

    # Print results
    print("\n" + "="*80)
    print("CORRELATION RESULTS")
    print("="*80)

    for key, corr in correlations.items():
        if 'error' in corr:
            print(f"\n{key}: {corr['error']}")
            continue

        print(f"\n{key}:")
        print(f"  n = {corr['n']}")
        print(f"  Pearson r = {corr['pearson_r']:.4f} (p = {corr['pearson_p']:.2e})")
        print(f"  Spearman ρ = {corr['spearman_rho']:.4f} (p = {corr['spearman_p']:.2e})")

        if corr.get('power_law'):
            pl = corr['power_law']
            print(f"  Power law: {pl['formula']}")
            print(f"  R² = {pl['r_squared']:.4f}")

            # Interpret strength
            rho = abs(corr['spearman_rho'])
            if rho > 0.7:
                strength = "STRONG"
            elif rho > 0.4:
                strength = "MODERATE"
            else:
                strength = "WEAK"
            print(f"  Correlation strength: {strength}")

    # Morphology analysis
    print("\n" + "="*80)
    print("MORPHOLOGY DEPENDENCE")
    print("="*80)
    morph_stats = analyze_morphology_dependence(data)

    for morph, stats in sorted(morph_stats.items()):
        print(f"\n{morph} (n={stats['n']}):")
        print(f"  Median ρ_crit: {stats['median']:.2f} M☉/pc²")
        print(f"  Mean ± std:    {stats['mean']:.2f} ± {stats['std']:.2f} M☉/pc²")
        print(f"  Range:         [{stats['min']:.2f}, {stats['max']:.2f}]")

    # Generate plots
    print("\nGenerating correlation plots...")
    plot_file = plot_correlations(data, correlations)
    print(f"Saved: {plot_file}")

    print("\n" + "="*80)
    print("PHYSICAL INTERPRETATION ANALYSIS COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
