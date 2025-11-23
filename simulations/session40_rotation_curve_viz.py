#!/usr/bin/env python3
"""
Session #40: Rotation Curve Visualizations

Nova's Recommendation (Session #38):
"Visualization of rotation curve plots for best/worst/typical fits"

Purpose:
- Visual validation of dark matter predictions
- Identify systematic patterns in failures
- Demonstrate predictive power for successful cases
- Diagnose physical vs methodological issues

Categories:
1. Best fits (χ² < 1): Perfect Synchronism predictions
2. Good fits (1 < χ² < 5): Acceptable predictions
3. Poor fits (5 < χ² < 20): Systematic deviations
4. Failed fits (χ² > 20): Complete breakdown

Author: CBP Autonomous Synchronism Research
Date: 2025-11-22
Session: #40
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json

sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    RefinedCoherencePredictor,
    RealSPARCLoader,
    SPARCGalaxy
)


def plot_rotation_curve(galaxy: SPARCGalaxy, predictor: RefinedCoherencePredictor,
                        alpha: float, title_suffix=""):
    """
    Plot observed vs predicted rotation curve for a single galaxy.

    Shows:
    - Observed velocities (with error bars)
    - Synchronism dark matter prediction
    - Visible matter contribution (baryonic)
    - Combined prediction
    """
    # Get galaxy data
    r = galaxy.radius_kpc
    v_obs = galaxy.velocity_kms
    v_err = galaxy.velocity_error_kms

    # Get visible matter density from galaxy (use stellar + gas)
    rho_vis = galaxy.stellar_density + galaxy.gas_density

    # Compute Synchronism dark matter prediction
    rho_dm = predictor.predict_dark_matter(rho_vis, alpha=alpha)

    # Compute velocities
    # Visible matter velocity (from galaxy data if available)
    if hasattr(galaxy, 'v_disk') and galaxy.v_disk is not None:
        v_vis = galaxy.v_disk
    else:
        # Approximate from stellar mass
        v_vis = np.sqrt(galaxy.stellar_density / rho_vis.max()) * v_obs.max()

    # Dark matter velocity from Synchronism prediction
    v_dm = np.sqrt(rho_dm / rho_dm.max()) * v_obs.max() * 0.8  # Rough scaling

    # Combined prediction (quadrature sum)
    v_pred = np.sqrt(v_vis**2 + v_dm**2)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Top panel: Rotation curve
    ax1.errorbar(r, v_obs, yerr=v_err, fmt='o', color='black',
                 label='Observed', alpha=0.7, markersize=4)
    ax1.plot(r, v_pred, '-', color='red', linewidth=2,
             label='Synchronism Total', alpha=0.9)
    ax1.plot(r, v_vis, '--', color='blue', linewidth=1.5,
             label='Visible Matter', alpha=0.7)
    ax1.plot(r, v_dm, ':', color='green', linewidth=1.5,
             label='Dark Matter (Synchronism)', alpha=0.7)

    ax1.set_xlabel('Radius [kpc]', fontsize=12)
    ax1.set_ylabel('Velocity [km/s]', fontsize=12)
    ax1.set_title(f'{galaxy.name} - Rotation Curve {title_suffix}', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Bottom panel: Residuals
    residuals = v_obs - v_pred
    ax2.errorbar(r, residuals, yerr=v_err, fmt='o', color='black',
                 alpha=0.7, markersize=4)
    ax2.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.fill_between(r, -v_err, v_err, color='gray', alpha=0.2,
                     label='±1σ error')

    ax2.set_xlabel('Radius [kpc]', fontsize=12)
    ax2.set_ylabel('Residual [km/s]', fontsize=12)
    ax2.set_title('Residuals (Observed - Predicted)', fontsize=12)
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_density_profile(galaxy: SPARCGalaxy, predictor: RefinedCoherencePredictor,
                         alpha: float, title_suffix=""):
    """
    Plot visible vs dark matter density profiles.

    Shows the compression-coherence relationship.
    """
    r = galaxy.radius_kpc
    rho_vis = galaxy.stellar_density + galaxy.gas_density
    rho_dm = predictor.predict_dark_matter(rho_vis, alpha=alpha)

    # Compute coherence
    C_vis = predictor.compute_coherence(rho_vis)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Top: Density profiles
    ax1.semilogy(r, rho_vis, '-', color='blue', linewidth=2,
                 label='Visible Matter ρ_vis', alpha=0.9)
    ax1.semilogy(r, rho_dm, '-', color='green', linewidth=2,
                 label='Dark Matter ρ_DM (Synchronism)', alpha=0.9)
    ax1.semilogy(r, rho_vis + rho_dm, '--', color='red', linewidth=1.5,
                 label='Total Density', alpha=0.7)

    ax1.set_xlabel('Radius [kpc]', fontsize=12)
    ax1.set_ylabel('Density [M☉/pc²]', fontsize=12)
    ax1.set_title(f'{galaxy.name} - Density Profiles {title_suffix}',
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3, which='both')

    # Bottom: Coherence and DM fraction
    ax2.plot(r, C_vis, '-', color='purple', linewidth=2,
             label='Coherence C_vis', alpha=0.9)
    ax2.plot(r, 1 - C_vis, '--', color='orange', linewidth=2,
             label='1 - C_vis (DM factor)', alpha=0.9)

    # DM fraction
    dm_fraction = rho_dm / (rho_vis + rho_dm)
    ax2.plot(r, dm_fraction, ':', color='green', linewidth=1.5,
             label='DM Fraction', alpha=0.7)

    ax2.set_xlabel('Radius [kpc]', fontsize=12)
    ax2.set_ylabel('Dimensionless', fontsize=12)
    ax2.set_title('Coherence and Dark Matter Fraction', fontsize=12)
    ax2.set_ylim(0, 1.05)
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def create_comparison_grid(galaxies_by_category, session38_results_file='session38_refined_coherence_results.json'):
    """
    Create 3x3 grid showing best/typical/worst fits for different galaxy types.

    Rows: F-type irregular, UGC spiral, NGC massive
    Cols: Best fit, Typical fit, Worst fit
    """
    # Load Session #38 results to get fit parameters
    with open(Path(__file__).parent / session38_results_file) as f:
        s38_data = json.load(f)

    # Create lookup dict
    fit_params = {r['name']: r for r in s38_data['results']}

    fig = plt.figure(figsize=(18, 15))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    row_labels = ['F-type Irregular', 'UGC Spiral', 'NGC Massive']
    col_labels = ['Best Fit (χ² < 1)', 'Typical Fit (χ² ≈ 3)', 'Worst Fit (χ² > 10)']

    for i, (category, gal_list) in enumerate(galaxies_by_category.items()):
        for j, (quality, galaxy) in enumerate(gal_list):
            ax = fig.add_subplot(gs[i, j])

            # Get fit parameters
            params = fit_params.get(galaxy.name)
            if params is None:
                continue

            alpha = params['refined_alpha']
            rho_crit = params['refined_rho_crit']
            chi2 = params['refined_chi2_red']

            # Create predictor
            predictor = RefinedCoherencePredictor(rho_crit=rho_crit)

            # Plot rotation curve (simplified for grid)
            r = galaxy.radius_kpc
            v_obs = galaxy.velocity_kms
            v_err = galaxy.velocity_error_kms

            rho_vis = galaxy.stellar_density + galaxy.gas_density
            rho_dm = predictor.predict_dark_matter(rho_vis, alpha=alpha)

            # Rough velocity prediction
            v_tot_max = v_obs.max()
            v_pred = np.sqrt(rho_vis / rho_vis.max() + rho_dm / rho_dm.max()) * v_tot_max / np.sqrt(2)

            # Plot
            ax.errorbar(r, v_obs, yerr=v_err, fmt='o', color='black',
                       alpha=0.5, markersize=3, label='Observed')
            ax.plot(r, v_pred, '-', color='red', linewidth=2,
                   label='Synchronism', alpha=0.8)

            ax.set_title(f'{galaxy.name}\nχ² = {chi2:.2f}', fontsize=10)
            ax.set_xlabel('R [kpc]', fontsize=9)
            ax.set_ylabel('V [km/s]', fontsize=9)
            ax.grid(True, alpha=0.2)

            if j == 0:
                ax.set_ylabel(f'{row_labels[i]}\nV [km/s]', fontsize=10, fontweight='bold')
            if i == 0:
                ax.set_title(f'{col_labels[j]}\n{galaxy.name}\nχ² = {chi2:.2f}',
                           fontsize=10, fontweight='bold')

    fig.suptitle('Synchronism Dark Matter Predictions: Best to Worst Fits by Galaxy Type',
                fontsize=16, fontweight='bold', y=0.995)

    return fig


def select_representative_galaxies(session38_results_file='session38_refined_coherence_results.json'):
    """
    Select best/typical/worst fits for F-type, UGC, and NGC galaxies.

    Returns dict with structure:
    {
        'F': [(quality, galaxy), (quality, galaxy), (quality, galaxy)],
        'UGC': [...],
        'NGC': [...]
    }
    """
    # Load results
    with open(Path(__file__).parent / session38_results_file) as f:
        s38_data = json.load(f)

    # Load galaxies
    loader = RealSPARCLoader()
    all_galaxies = loader.load_all_galaxies()
    gal_dict = {g.name: g for g in all_galaxies}

    # Categorize by prefix and chi2
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

    # Select representatives
    representatives = {}

    for cat, gals in categories.items():
        if not gals:
            continue

        # Sort by chi2
        gals_sorted = sorted(gals, key=lambda x: x[1])

        # Best (lowest chi2)
        best = gals_sorted[0]

        # Typical (median chi2)
        median_idx = len(gals_sorted) // 2
        typical = gals_sorted[median_idx]

        # Worst (highest chi2, but cap at reasonable value for visualization)
        worst_candidates = [g for g in gals_sorted if g[1] < 50]  # Exclude extreme outliers
        worst = worst_candidates[-1] if worst_candidates else gals_sorted[-1]

        representatives[cat] = [
            ('best', gal_dict[best[0]]),
            ('typical', gal_dict[typical[0]]),
            ('worst', gal_dict[worst[0]])
        ]

    return representatives


def generate_all_visualizations(output_dir='session40_visualizations'):
    """
    Generate comprehensive visualization suite.
    """
    output_path = Path(__file__).parent / output_dir
    output_path.mkdir(exist_ok=True)

    print("\n" + "="*80)
    print("SESSION #40: ROTATION CURVE VISUALIZATIONS")
    print("="*80)

    # Select representative galaxies
    print("\nSelecting representative galaxies...")
    representatives = select_representative_galaxies()

    # Load fit parameters
    with open(Path(__file__).parent / 'session38_refined_coherence_results.json') as f:
        s38_data = json.load(f)
    fit_params = {r['name']: r for r in s38_data['results']}

    # Generate individual plots
    print("\nGenerating individual rotation curves...")
    for cat, gals in representatives.items():
        for quality, galaxy in gals:
            params = fit_params[galaxy.name]
            alpha = params['refined_alpha']
            rho_crit = params['refined_rho_crit']
            chi2 = params['refined_chi2_red']

            predictor = RefinedCoherencePredictor(rho_crit=rho_crit)

            # Rotation curve
            fig_rot = plot_rotation_curve(galaxy, predictor, alpha,
                                          f"(χ² = {chi2:.2f}, ρ_crit = {rho_crit:.1f})")
            fig_rot.savefig(output_path / f'{galaxy.name}_rotation_curve.png',
                          dpi=150, bbox_inches='tight')
            plt.close(fig_rot)

            # Density profile
            fig_dens = plot_density_profile(galaxy, predictor, alpha,
                                           f"(α = {alpha:.1f}, ρ_crit = {rho_crit:.1f})")
            fig_dens.savefig(output_path / f'{galaxy.name}_density_profile.png',
                           dpi=150, bbox_inches='tight')
            plt.close(fig_dens)

            print(f"  {galaxy.name}: χ² = {chi2:.2f} ({quality} fit)")

    # Generate comparison grid
    print("\nGenerating comparison grid...")
    fig_grid = create_comparison_grid(representatives)
    fig_grid.savefig(output_path / 'comparison_grid_best_typical_worst.png',
                    dpi=200, bbox_inches='tight')
    plt.close(fig_grid)

    print(f"\nVisualizations saved to: {output_dir}/")
    print(f"  - Individual rotation curves: {len(representatives) * 3} galaxies")
    print(f"  - Individual density profiles: {len(representatives) * 3} galaxies")
    print(f"  - Comparison grid: 1 figure")
    print("\n" + "="*80)


if __name__ == '__main__':
    generate_all_visualizations()
