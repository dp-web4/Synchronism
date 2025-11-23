#!/usr/bin/env python3
"""
Session #40: Simple Rotation Curve Visualizations

Create publication-quality rotation curve plots showing:
1. Observed rotation curves (with errors)
2. Synchronism dark matter predictions
3. Best/typical/worst fits

Uses actual SPARC data structure.

Author: CBP Autonomous Synchronism Research
Date: 2025-11-22
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json

sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import RefinedCoherencePredictor
from synchronism_real_sparc_validation import RealSPARCLoader


def plot_simple_rotation_curve(galaxy, alpha, rho_crit, chi2, output_dir):
    """
    Plot rotation curve using actual SPARC data structure.
    """
    # Get data
    r = galaxy.radius  # kpc
    v_obs = galaxy.v_obs  # km/s
    v_err = galaxy.v_err  # km/s

    # Create predictor and compute predictions using existing infrastructure
    predictor = RefinedCoherencePredictor(rho_crit=rho_crit)

    # Use the predictor's built-in rotation curve method
    v_pred, v_baryon, v_dm = predictor.predict_rotation_curve(galaxy, alpha)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Top: Rotation curve
    ax1.errorbar(r, v_obs, yerr=v_err, fmt='o', color='black',
                 label='Observed', alpha=0.7, markersize=5, capsize=3)
    ax1.plot(r, v_pred, '-', color='red', linewidth=2.5,
             label='Synchronism Prediction', alpha=0.9)
    ax1.plot(r, v_baryon, '--', color='blue', linewidth=1.5,
             label='Baryonic', alpha=0.7)
    ax1.plot(r, v_dm, ':', color='green', linewidth=1.5,
             label='Dark Matter', alpha=0.7)

    ax1.set_xlabel('Radius [kpc]', fontsize=12)
    ax1.set_ylabel('Velocity [km/s]', fontsize=12)
    ax1.set_title(f'{galaxy.name} - Rotation Curve\n' +
                  f'χ² = {chi2:.2f}, α = {alpha:.1f}, ρ_crit = {rho_crit:.1f} M☉/pc²',
                  fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, v_obs.max() * 1.2)

    # Bottom: Residuals
    residuals = v_obs - v_pred
    ax2.errorbar(r, residuals, yerr=v_err, fmt='o', color='black',
                 alpha=0.7, markersize=5, capsize=3)
    ax2.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.fill_between(r, -v_err, v_err, color='gray', alpha=0.2,
                     label='±1σ error')

    ax2.set_xlabel('Radius [kpc]', fontsize=12)
    ax2.set_ylabel('Residual [km/s]', fontsize=12)
    ax2.set_title('Residuals (Observed - Predicted)', fontsize=11)
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save
    output_file = output_dir / f'{galaxy.name}_rotation_curve.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return output_file


def create_comparison_figure(galaxies_dict, fit_params, output_dir):
    """
    Create 3x3 comparison grid: F/UGC/NGC × best/typical/worst
    """
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))

    row_labels = ['F-type Irregular', 'UGC Spiral', 'NGC Massive']
    col_labels = ['Best Fit', 'Typical Fit', 'Worst Fit']

    for i, (cat, gals) in enumerate(galaxies_dict.items()):
        for j, (quality, galaxy) in enumerate(gals):
            ax = axes[i, j]

            # Get fit parameters
            params = fit_params.get(galaxy.name)
            if params is None:
                ax.text(0.5, 0.5, f'{galaxy.name}\nNo data',
                       ha='center', va='center', transform=ax.transAxes)
                continue

            alpha = params['refined_alpha']
            rho_crit = params['refined_rho_crit']
            chi2 = params['refined_chi2_red']

            # Get data
            r = galaxy.radius
            v_obs = galaxy.v_obs
            v_err = galaxy.v_err

            # Compute prediction using built-in method
            predictor = RefinedCoherencePredictor(rho_crit=rho_crit)
            v_pred, _, _ = predictor.predict_rotation_curve(galaxy, alpha)

            # Plot
            ax.errorbar(r, v_obs, yerr=v_err, fmt='o', color='black',
                       alpha=0.5, markersize=3, capsize=2)
            ax.plot(r, v_pred, '-', color='red', linewidth=2, alpha=0.8)

            # Labels
            if j == 0:
                ax.set_ylabel(f'{row_labels[i]}\nV [km/s]', fontsize=10, fontweight='bold')
            else:
                ax.set_ylabel('V [km/s]', fontsize=9)

            if i == 2:
                ax.set_xlabel('R [kpc]', fontsize=9)

            ax.set_title(f'{galaxy.name}\nχ² = {chi2:.2f}', fontsize=9)
            ax.grid(True, alpha=0.2)

    # Super title
    fig.suptitle('Synchronism Dark Matter: Best/Typical/Worst Fits by Galaxy Type',
                fontsize=14, fontweight='bold', y=0.995)

    # Save
    output_file = output_dir / 'comparison_grid_3x3.png'
    fig.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close(fig)

    return output_file


def select_representatives(s38_results_file='session38_refined_coherence_results.json'):
    """
    Select best/typical/worst for each galaxy type.
    """
    # Load results
    with open(Path(__file__).parent / s38_results_file) as f:
        s38_data = json.load(f)

    # Load all galaxies
    loader = RealSPARCLoader()
    all_galaxies = loader.load_all_galaxies()
    gal_dict = {g.name: g for g in all_galaxies}

    # Categorize
    categories = {'F': [], 'UGC': [], 'NGC': []}

    for r in s38_data['results']:
        name = r['name']
        chi2 = r['refined_chi2_red']

        if name.startswith('F'):
            categories['F'].append((name, chi2))
        elif name.startswith('UGC'):
            categories['UGC'].append((name, chi2))
        elif name.startswith('NGC'):
            categories['NGC'].append((name, chi2))

    # Select best/typical/worst
    representatives = {}

    for cat, gals in categories.items():
        if not gals:
            continue

        gals_sorted = sorted(gals, key=lambda x: x[1])

        # Best
        best = gals_sorted[0]

        # Typical (median)
        typical = gals_sorted[len(gals_sorted) // 2]

        # Worst (but reasonable for viz)
        worst_candidates = [g for g in gals_sorted if g[1] < 30]
        worst = worst_candidates[-1] if worst_candidates else gals_sorted[-1]

        representatives[cat] = [
            ('best', gal_dict[best[0]]),
            ('typical', gal_dict[typical[0]]),
            ('worst', gal_dict[worst[0]])
        ]

    return representatives


def main():
    output_dir = Path(__file__).parent / 'session40_rotation_curves'
    output_dir.mkdir(exist_ok=True)

    print("\n" + "="*80)
    print("SESSION #40: ROTATION CURVE VISUALIZATIONS")
    print("="*80)

    # Load fit parameters
    with open(Path(__file__).parent / 'session38_refined_coherence_results.json') as f:
        s38_data = json.load(f)
    fit_params = {r['name']: r for r in s38_data['results']}

    # Select representatives
    print("\nSelecting representative galaxies...")
    representatives = select_representatives()

    # Generate individual plots
    print("\nGenerating individual rotation curves...")
    n_plots = 0
    for cat, gals in representatives.items():
        for quality, galaxy in gals:
            params = fit_params[galaxy.name]
            alpha = params['refined_alpha']
            rho_crit = params['refined_rho_crit']
            chi2 = params['refined_chi2_red']

            output_file = plot_simple_rotation_curve(galaxy, alpha, rho_crit, chi2, output_dir)
            print(f"  {galaxy.name}: χ² = {chi2:.2f} ({quality} {cat})")
            n_plots += 1

    # Create comparison grid
    print("\nGenerating comparison grid...")
    grid_file = create_comparison_figure(representatives, fit_params, output_dir)

    print(f"\n✓ Generated {n_plots} individual rotation curves")
    print(f"✓ Generated comparison grid: {grid_file.name}")
    print(f"\nAll visualizations saved to: {output_dir}/")
    print("="*80)


if __name__ == '__main__':
    main()
