#!/usr/bin/env python3
"""
SPARC Correlation Analysis - Session #28
=========================================

Analyze what galaxy properties correlate with good Synchronism fits.

Questions:
1. Do high surface brightness galaxies fit better than LSB?
2. Do massive galaxies fit better than dwarfs?
3. Does gas fraction affect fit quality?
4. Are there morphological patterns?

This will identify WHERE Synchronism works vs where it needs refinement.

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-14
Session: #28 - Cross-domain coherence validation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import json
import os

# Import validator to get results
import sys
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_real_sparc_validation import RealSPARCLoader, SynchronismPredictor, SPARCValidator


def load_sparc_metadata(data_dir="sparc_real_data"):
    """
    Load SPARC galaxy metadata from MassModels table.

    Returns dict: galaxy_name -> {mass, distance, etc.}
    """
    import re

    metadata = {}
    mrt_file = os.path.join(data_dir, "MassModels_Lelli2016c.mrt")

    if not os.path.exists(mrt_file):
        print(f"Warning: Could not find {mrt_file}")
        return metadata

    # Parse MRT table (CDS format)
    with open(mrt_file, 'r') as f:
        in_data = False
        for line in f:
            # Look for start of data
            if line.startswith("--------"):
                in_data = True
                continue

            if not in_data:
                continue

            # Skip empty lines
            if not line.strip():
                continue

            # Parse data line
            parts = line.split()
            if len(parts) < 5:
                continue

            try:
                galaxy = parts[0]
                T_type = int(parts[1])  # Morphological type
                D = float(parts[2])  # Distance (Mpc)
                Inc = float(parts[3])  # Inclination (deg)
                L = float(parts[4])  # Luminosity (10^9 L_sun at 3.6um)

                metadata[galaxy] = {
                    'morphological_type': T_type,
                    'distance_mpc': D,
                    'inclination_deg': Inc,
                    'luminosity_36um': L,
                    'log_luminosity': np.log10(L + 1e-10)  # Avoid log(0)
                }
            except (ValueError, IndexError):
                continue

    return metadata


def analyze_correlations():
    """Analyze correlations between galaxy properties and fit quality."""

    print("=" * 70)
    print("Session #28: SPARC Correlation Analysis")
    print("=" * 70)
    print()
    print("Analyzing what galaxy properties predict good Synchronism fits...")
    print()

    # Load galaxies
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=None)

    print(f"Loaded {len(galaxies)} galaxies")

    # Load metadata
    metadata = load_sparc_metadata()
    print(f"Loaded metadata for {len(metadata)} galaxies")
    print()

    # Run validation
    predictor = SynchronismPredictor(
        gamma=SynchronismPredictor.GAMMA,
        beta=SynchronismPredictor.BETA
    )
    validator = SPARCValidator(predictor)

    print("Running validation to get χ² values...")
    stats = validator.validate_sample(galaxies)

    # Extract validation results for each galaxy
    galaxy_results = []
    for galaxy in galaxies:
        # Re-validate to get individual results
        result = validator.validate_galaxy(galaxy)
        if result is None:
            continue

        # Get metadata
        meta = metadata.get(galaxy.name, {})

        galaxy_results.append({
            'name': galaxy.name,
            'chi2_red': result['chi2_red'],
            'alpha': result['alpha_best'],
            'M_DM': result.get('M_DM', 0),
            'M_vis': result.get('M_vis', 0),
            'M_DM_Mvis': result.get('M_DM', 0) / result.get('M_vis', 1) if result.get('M_vis', 0) > 0 else 0,
            'luminosity': meta.get('luminosity_36um', np.nan),
            'log_luminosity': meta.get('log_luminosity', np.nan),
            'morphological_type': meta.get('morphological_type', np.nan),
            'distance': meta.get('distance_mpc', np.nan),
            'inclination': meta.get('inclination_deg', np.nan),
        })

    print(f"Analyzed {len(galaxy_results)} galaxies with complete data")
    print()

    # Convert to arrays for correlation analysis
    chi2_reds = np.array([g['chi2_red'] for g in galaxy_results])
    log_lums = np.array([g['log_luminosity'] for g in galaxy_results])
    morph_types = np.array([g['morphological_type'] for g in galaxy_results])
    alphas = np.array([g['alpha'] for g in galaxy_results])

    # Remove NaN values for correlation
    valid_lum = ~np.isnan(log_lums) & ~np.isnan(chi2_reds)
    valid_morph = ~np.isnan(morph_types) & ~np.isnan(chi2_reds)

    print("=" * 70)
    print("CORRELATION ANALYSIS")
    print("=" * 70)
    print()

    # Correlation 1: χ²_red vs Luminosity (proxy for mass)
    if np.sum(valid_lum) > 10:
        corr_lum, p_lum = spearmanr(log_lums[valid_lum], chi2_reds[valid_lum])
        print(f"χ²_red vs log(Luminosity):")
        print(f"  Spearman ρ = {corr_lum:.3f} (p = {p_lum:.4f})")
        if p_lum < 0.05:
            if corr_lum < 0:
                print(f"  ✓ SIGNIFICANT: Brighter galaxies fit BETTER")
            else:
                print(f"  ✓ SIGNIFICANT: Brighter galaxies fit WORSE")
        else:
            print(f"  No significant correlation")
        print()

    # Correlation 2: χ²_red vs Morphological Type
    if np.sum(valid_morph) > 10:
        corr_morph, p_morph = spearmanr(morph_types[valid_morph], chi2_reds[valid_morph])
        print(f"χ²_red vs Morphological Type:")
        print(f"  Spearman ρ = {corr_morph:.3f} (p = {p_morph:.4f})")
        print(f"  (T < 0: ellipticals, T ~ 3-7: spirals, T > 8: irregulars)")
        if p_morph < 0.05:
            print(f"  ✓ SIGNIFICANT correlation detected")
        else:
            print(f"  No significant correlation")
        print()

    # Binned analysis: Good fits vs Poor fits
    good_fits = chi2_reds < 2
    poor_fits = chi2_reds > 20

    print("=" * 70)
    print("GOOD FITS vs POOR FITS")
    print("=" * 70)
    print()

    print(f"Good fits (χ²_red < 2): {np.sum(good_fits)} galaxies")
    print(f"Poor fits (χ²_red > 20): {np.sum(poor_fits)} galaxies")
    print()

    # Compare properties
    if np.sum(good_fits & valid_lum) > 0 and np.sum(poor_fits & valid_lum) > 0:
        lum_good = log_lums[good_fits & valid_lum]
        lum_poor = log_lums[poor_fits & valid_lum]

        print("Luminosity comparison:")
        print(f"  Good fits: log(L) = {np.mean(lum_good):.2f} ± {np.std(lum_good):.2f}")
        print(f"  Poor fits: log(L) = {np.mean(lum_poor):.2f} ± {np.std(lum_poor):.2f}")
        print()

    if np.sum(good_fits & valid_morph) > 0 and np.sum(poor_fits & valid_morph) > 0:
        morph_good = morph_types[good_fits & valid_morph]
        morph_poor = morph_types[poor_fits & valid_morph]

        print("Morphological type comparison:")
        print(f"  Good fits: T = {np.median(morph_good):.1f} (median)")
        print(f"  Poor fits: T = {np.median(morph_poor):.1f} (median)")
        print()

    # Create visualization
    create_correlation_plots(galaxy_results)

    # Save results
    save_correlation_results(galaxy_results, {
        'chi2_vs_luminosity': (corr_lum, p_lum) if np.sum(valid_lum) > 10 else (np.nan, np.nan),
        'chi2_vs_morphology': (corr_morph, p_morph) if np.sum(valid_morph) > 10 else (np.nan, np.nan),
    })

    print("=" * 70)
    print("Analysis complete! Results saved to:")
    print("  - Session28_SPARC_Correlations.png")
    print("  - Session28_correlation_results.json")
    print("=" * 70)


def create_correlation_plots(galaxy_results):
    """Create correlation visualization plots."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    chi2_reds = np.array([g['chi2_red'] for g in galaxy_results])
    log_lums = np.array([g['log_luminosity'] for g in galaxy_results])
    morph_types = np.array([g['morphological_type'] for g in galaxy_results])
    alphas = np.array([g['alpha'] for g in galaxy_results])

    # Plot 1: χ²_red vs Luminosity
    ax = axes[0, 0]
    valid = ~np.isnan(log_lums) & ~np.isnan(chi2_reds)
    ax.scatter(log_lums[valid], chi2_reds[valid], alpha=0.6, s=30)
    ax.axhline(y=2, color='g', linestyle='--', label='χ²_red = 2 (good fit)')
    ax.axhline(y=5, color='orange', linestyle='--', label='χ²_red = 5 (acceptable)')
    ax.set_xlabel('log(Luminosity) [log(10⁹ L☉)]')
    ax.set_ylabel('χ²_red')
    ax.set_yscale('log')
    ax.set_title('Fit Quality vs Galaxy Luminosity')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: χ²_red vs Morphological Type
    ax = axes[0, 1]
    valid = ~np.isnan(morph_types) & ~np.isnan(chi2_reds)
    ax.scatter(morph_types[valid], chi2_reds[valid], alpha=0.6, s=30)
    ax.axhline(y=2, color='g', linestyle='--')
    ax.axhline(y=5, color='orange', linestyle='--')
    ax.set_xlabel('Morphological Type T')
    ax.set_ylabel('χ²_red')
    ax.set_yscale('log')
    ax.set_title('Fit Quality vs Morphology')
    ax.text(0.05, 0.95, 'T<0: E/S0\nT~3-7: Spirals\nT>8: Irregulars',
            transform=ax.transAxes, va='top', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.grid(True, alpha=0.3)

    # Plot 3: Histogram of χ²_red
    ax = axes[1, 0]
    chi2_finite = chi2_reds[np.isfinite(chi2_reds)]
    ax.hist(chi2_finite, bins=50, alpha=0.7, edgecolor='black')
    ax.axvline(x=2, color='g', linestyle='--', linewidth=2, label='Good fit threshold')
    ax.axvline(x=np.median(chi2_finite), color='r', linestyle='-', linewidth=2,
               label=f'Median = {np.median(chi2_finite):.1f}')
    ax.set_xlabel('χ²_red')
    ax.set_ylabel('Number of galaxies')
    ax.set_title('Distribution of Fit Quality')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: α vs Luminosity
    ax = axes[1, 1]
    valid = ~np.isnan(log_lums) & ~np.isnan(alphas)
    ax.scatter(log_lums[valid], alphas[valid], alpha=0.6, s=30, c=chi2_reds[valid],
               cmap='RdYlGn_r', norm=plt.matplotlib.colors.LogNorm())
    ax.set_xlabel('log(Luminosity) [log(10⁹ L☉)]')
    ax.set_ylabel('α (DM normalization)')
    ax.set_title('Dark Matter Normalization vs Luminosity')
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('χ²_red')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('Session28_SPARC_Correlations.png', dpi=150, bbox_inches='tight')
    plt.close()


def save_correlation_results(galaxy_results, correlations):
    """Save correlation analysis results to JSON."""

    output = {
        'n_galaxies': len(galaxy_results),
        'correlations': {
            'chi2_vs_luminosity': {
                'rho': correlations['chi2_vs_luminosity'][0],
                'p_value': correlations['chi2_vs_luminosity'][1]
            },
            'chi2_vs_morphology': {
                'rho': correlations['chi2_vs_morphology'][0],
                'p_value': correlations['chi2_vs_morphology'][1]
            }
        },
        'summary_statistics': {
            'median_chi2_red': float(np.median([g['chi2_red'] for g in galaxy_results])),
            'fraction_good_fits': float(np.mean([g['chi2_red'] < 2 for g in galaxy_results])),
            'fraction_acceptable_fits': float(np.mean([g['chi2_red'] < 5 for g in galaxy_results]))
        },
        'best_fits': [
            {k: v for k, v in g.items() if not isinstance(v, float) or not np.isnan(v)}
            for g in sorted(galaxy_results, key=lambda x: x['chi2_red'])[:10]
        ],
        'worst_fits': [
            {k: v for k, v in g.items() if not isinstance(v, float) or not np.isnan(v)}
            for g in sorted(galaxy_results, key=lambda x: x['chi2_red'], reverse=True)[:10]
        ]
    }

    with open('Session28_correlation_results.json', 'w') as f:
        json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)


if __name__ == "__main__":
    analyze_correlations()
