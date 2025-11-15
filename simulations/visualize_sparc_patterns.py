#!/usr/bin/env python3
"""
Visualization of Systematic Patterns in SPARC Results
=======================================================

Creates publication-quality plots showing where Synchronism works vs fails.

Author: Autonomous Research Agent
Date: 2025-11-14
Session: #17
"""

import numpy as np
import matplotlib.pyplot as plt
from synchronism_real_sparc_validation import RealSPARCLoader, SynchronismPredictor, SPARCValidator


def create_comprehensive_plots():
    """Generate full visualization suite for SPARC results."""

    print("Creating comprehensive visualization plots...")

    # Load and validate
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=None)

    predictor = SynchronismPredictor()
    validator = SPARCValidator(predictor)
    stats = validator.validate_sample(galaxies)

    results = validator.results

    # Extract data
    chi2_reds = np.array([r['chi2_red'] for r in results])
    alphas = np.array([r['alpha_best'] for r in results])
    mass_ratios = np.array([r['M_DM_Mvis'] for r in results])
    names = np.array([r['name'] for r in results])

    # Categorize by catalog
    is_ngc = np.array(['NGC' in name for name in names])
    is_ugc = np.array(['UGC' in name for name in names])
    is_ddo = np.array(['DDO' in name for name in names])
    is_f = np.array([name.startswith('F') for name in names])
    is_other = ~(is_ngc | is_ugc | is_ddo | is_f)

    # Create comprehensive figure
    fig = plt.figure(figsize=(20, 12))

    # Plot 1: χ²_red distribution
    ax1 = plt.subplot(2, 3, 1)
    bins = np.logspace(-1, 3, 50)
    ax1.hist(chi2_reds, bins=bins, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axvline(2, color='green', linestyle='--', linewidth=2, label='Excellent (χ²_red=2)')
    ax1.axvline(5, color='orange', linestyle='--', linewidth=2, label='Acceptable (χ²_red=5)')
    ax1.axvline(np.median(chi2_reds), color='red', linestyle='-', linewidth=2,
                label=f'Median = {np.median(chi2_reds):.1f}')
    ax1.set_xscale('log')
    ax1.set_xlabel('Reduced χ²', fontsize=12)
    ax1.set_ylabel('Number of Galaxies', fontsize=12)
    ax1.set_title('Fit Quality Distribution (175 Galaxies)', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Success rate by catalog
    ax2 = plt.subplot(2, 3, 2)
    catalogs = ['F', 'UGC', 'NGC', 'DDO', 'Other']
    catalog_data = [is_f, is_ugc, is_ngc, is_ddo, is_other]
    excellent_rates = []
    acceptable_rates = []

    for cat_mask in catalog_data:
        cat_chi2 = chi2_reds[cat_mask]
        if len(cat_chi2) > 0:
            excellent_rates.append(100 * np.mean(cat_chi2 < 2))
            acceptable_rates.append(100 * np.mean(cat_chi2 < 5))
        else:
            excellent_rates.append(0)
            acceptable_rates.append(0)

    x = np.arange(len(catalogs))
    width = 0.35
    ax2.bar(x - width/2, excellent_rates, width, label='Excellent (χ²<2)', color='green', alpha=0.7)
    ax2.bar(x + width/2, acceptable_rates, width, label='Acceptable (χ²<5)', color='orange', alpha=0.7)
    ax2.set_xlabel('Galaxy Catalog', fontsize=12)
    ax2.set_ylabel('Success Rate (%)', fontsize=12)
    ax2.set_title('Synchronism Success by Catalog', fontsize=13, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(catalogs)
    ax2.legend(fontsize=10)
    ax2.grid(True, axis='y', alpha=0.3)

    # Plot 3: α distribution
    ax3 = plt.subplot(2, 3, 3)
    ax3.hist(alphas, bins=30, color='purple', alpha=0.7, edgecolor='black')
    ax3.axvline(np.median(alphas), color='red', linestyle='--', linewidth=2,
                label=f'Median = {np.median(alphas):.1f}')
    ax3.set_xlabel('Best-fit α (DM normalization)', fontsize=12)
    ax3.set_ylabel('Number of Galaxies', fontsize=12)
    ax3.set_title('Dark Matter Normalization Distribution', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)

    # Plot 4: M_DM/M_vis distribution
    ax4 = plt.subplot(2, 3, 4)
    ax4.hist(mass_ratios, bins=30, color='brown', alpha=0.7, edgecolor='black')
    ax4.axvline(np.median(mass_ratios), color='red', linestyle='--', linewidth=2,
                label=f'Median = {np.median(mass_ratios):.1f}')
    ax4.axvspan(10, 100, alpha=0.2, color='green', label='Expected range')
    ax4.set_xlabel('M_DM / M_vis', fontsize=12)
    ax4.set_ylabel('Number of Galaxies', fontsize=12)
    ax4.set_title('Dark Matter Fraction', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)

    # Plot 5: χ² vs α correlation
    ax5 = plt.subplot(2, 3, 5)
    colors = []
    for name in names:
        if 'NGC' in name:
            colors.append('blue')
        elif 'UGC' in name:
            colors.append('green')
        elif 'DDO' in name:
            colors.append('red')
        elif name.startswith('F'):
            colors.append('purple')
        else:
            colors.append('gray')

    ax5.scatter(alphas, chi2_reds, c=colors, alpha=0.6, s=30)
    ax5.set_xlabel('Best-fit α', fontsize=12)
    ax5.set_ylabel('Reduced χ²', fontsize=12)
    ax5.set_yscale('log')
    ax5.set_title('Fit Quality vs DM Normalization', fontsize=13, fontweight='bold')
    ax5.axhline(2, color='green', linestyle='--', alpha=0.5)
    ax5.axhline(5, color='orange', linestyle='--', alpha=0.5)
    ax5.grid(True, alpha=0.3)

    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label='NGC'),
                      Patch(facecolor='green', label='UGC'),
                      Patch(facecolor='purple', label='F'),
                      Patch(facecolor='red', label='DDO')]
    ax5.legend(handles=legend_elements, fontsize=9, loc='upper right')

    # Plot 6: Cumulative distribution
    ax6 = plt.subplot(2, 3, 6)
    sorted_chi2 = np.sort(chi2_reds)
    cumulative = np.arange(1, len(sorted_chi2) + 1) / len(sorted_chi2) * 100
    ax6.plot(sorted_chi2, cumulative, 'b-', linewidth=2)
    ax6.axvline(2, color='green', linestyle='--', linewidth=2, label='Excellent threshold')
    ax6.axvline(5, color='orange', linestyle='--', linewidth=2, label='Acceptable threshold')
    ax6.axhline(50, color='red', linestyle=':', linewidth=1.5, label='Median')
    ax6.set_xscale('log')
    ax6.set_xlabel('Reduced χ²', fontsize=12)
    ax6.set_ylabel('Cumulative Percentage', fontsize=12)
    ax6.set_title('Cumulative Fit Quality', fontsize=13, fontweight='bold')
    ax6.legend(fontsize=10)
    ax6.grid(True, alpha=0.3)

    # Add text summary
    text_str = f"""
SPARC Full Sample Results (175 Galaxies)
γ = β = 0.30 (theory-predicted, NO TUNING)

Median χ²_red: {np.median(chi2_reds):.1f}
Excellent (χ²<2): {np.sum(chi2_reds < 2)} ({np.sum(chi2_reds < 2)/175*100:.1f}%)
Acceptable (χ²<5): {np.sum(chi2_reds < 5)} ({np.sum(chi2_reds < 5)/175*100:.1f}%)

Best catalog: F galaxies (75% acceptable)
Median M_DM/M_vis: {np.median(mass_ratios):.1f}
"""
    fig.text(0.02, 0.02, text_str, fontsize=10, family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig('Session17_Full_SPARC_Analysis.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: Session17_Full_SPARC_Analysis.png")

    # Create second figure focusing on best vs worst fits
    fig2, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Get indices of best and worst fits
    sorted_indices = np.argsort(chi2_reds)
    best_6_indices = sorted_indices[:6]
    worst_6_indices = sorted_indices[-6:]

    # Plot best 3 fits
    for i in range(3):
        idx = best_6_indices[i]
        ax = axes[0, i]
        galaxy = galaxies[idx]
        result = results[idx]

        ax.errorbar(galaxy.radius, galaxy.v_obs, yerr=galaxy.v_err,
                   fmt='o', color='black', label='SPARC obs',
                   markersize=4, capsize=2)
        ax.plot(galaxy.radius, result['v_sync'], '-', color='green',
               linewidth=2, label='Synchronism')

        ax.text(0.05, 0.95,
               f"{result['name']}\nχ²_red = {result['chi2_red']:.2f}\n"
               f"α = {result['alpha_best']:.1f}\n"
               f"M_DM/M_vis = {result['M_DM_Mvis']:.1f}",
               transform=ax.transAxes, verticalalignment='top',
               fontsize=9, bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

        ax.set_xlabel('Radius (kpc)', fontsize=10)
        ax.set_ylabel('Velocity (km/s)', fontsize=10)
        ax.set_title(f'BEST FIT #{i+1}', fontsize=11, fontweight='bold', color='green')
        ax.legend(fontsize=8, loc='lower right')
        ax.grid(True, alpha=0.3)

    # Plot worst 3 fits
    for i in range(3):
        idx = worst_6_indices[-(i+1)]
        ax = axes[1, i]
        galaxy = galaxies[idx]
        result = results[idx]

        ax.errorbar(galaxy.radius, galaxy.v_obs, yerr=galaxy.v_err,
                   fmt='o', color='black', label='SPARC obs',
                   markersize=4, capsize=2)
        ax.plot(galaxy.radius, result['v_sync'], '-', color='red',
               linewidth=2, label='Synchronism')

        ax.text(0.05, 0.95,
               f"{result['name']}\nχ²_red = {result['chi2_red']:.1f}\n"
               f"α = {result['alpha_best']:.1f}\n"
               f"M_DM/M_vis = {result['M_DM_Mvis']:.1f}",
               transform=ax.transAxes, verticalalignment='top',
               fontsize=9, bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

        ax.set_xlabel('Radius (kpc)', fontsize=10)
        ax.set_ylabel('Velocity (km/s)', fontsize=10)
        ax.set_title(f'WORST FIT #{i+1}', fontsize=11, fontweight='bold', color='red')
        ax.legend(fontsize=8, loc='lower right')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('Session17_BestWorst_Comparison.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: Session17_BestWorst_Comparison.png")

    print()
    print("Visualization complete!")
    print()


if __name__ == "__main__":
    create_comprehensive_plots()
