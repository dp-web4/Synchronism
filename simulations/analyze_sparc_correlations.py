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


def load_sparc_metadata(metadata_file="sparc_galaxy_metadata.json"):
    """
    Load SPARC galaxy metadata from JSON file.

    Returns dict: galaxy_name -> {distance, v_max, luminosity, etc.}
    """
    if not os.path.exists(metadata_file):
        print(f"Warning: Could not find {metadata_file}")
        print("Run extract_sparc_metadata.py first to generate metadata file.")
        return {}

    with open(metadata_file, 'r') as f:
        metadata = json.load(f)

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
            'v_max': meta.get('v_max', np.nan),
            'r_last': meta.get('r_last', np.nan),
            'luminosity_estimate': meta.get('luminosity_estimate', np.nan),
            'log_luminosity': meta.get('log_luminosity', np.nan),
            'mass_proxy': meta.get('mass_proxy', np.nan),
            'morphology': meta.get('morphology', 'unknown'),
            'distance': meta.get('distance_mpc', np.nan),
            'n_points': meta.get('n_points', np.nan),
        })

    print(f"Analyzed {len(galaxy_results)} galaxies with complete data")
    print()

    # Convert to arrays for correlation analysis
    chi2_reds = np.array([g['chi2_red'] for g in galaxy_results])
    v_maxes = np.array([g['v_max'] for g in galaxy_results])
    log_lums = np.array([g['log_luminosity'] for g in galaxy_results])
    mass_proxies = np.array([g['mass_proxy'] for g in galaxy_results])
    r_lasts = np.array([g['r_last'] for g in galaxy_results])
    alphas = np.array([g['alpha'] for g in galaxy_results])

    # Remove NaN values for correlation
    valid_lum = ~np.isnan(log_lums) & ~np.isnan(chi2_reds)
    valid_vmax = ~np.isnan(v_maxes) & ~np.isnan(chi2_reds)
    valid_mass = ~np.isnan(mass_proxies) & ~np.isnan(chi2_reds)
    valid_rlast = ~np.isnan(r_lasts) & ~np.isnan(chi2_reds)

    print("=" * 70)
    print("CORRELATION ANALYSIS")
    print("=" * 70)
    print()

    # Correlation 1: χ²_red vs V_max (mass proxy)
    if np.sum(valid_vmax) > 10:
        corr_vmax, p_vmax = spearmanr(v_maxes[valid_vmax], chi2_reds[valid_vmax])
        print(f"χ²_red vs V_max (velocity):")
        print(f"  Spearman ρ = {corr_vmax:.3f} (p = {p_vmax:.4f})")
        if p_vmax < 0.05:
            if corr_vmax < 0:
                print(f"  ✓ SIGNIFICANT: Faster-rotating galaxies fit BETTER")
            else:
                print(f"  ✓ SIGNIFICANT: Faster-rotating galaxies fit WORSE")
        else:
            print(f"  No significant correlation")
        print()

    # Correlation 2: χ²_red vs Mass proxy (V²R)
    if np.sum(valid_mass) > 10:
        corr_mass, p_mass = spearmanr(np.log10(mass_proxies[valid_mass] + 1), chi2_reds[valid_mass])
        print(f"χ²_red vs log(Mass proxy) [V²R]:")
        print(f"  Spearman ρ = {corr_mass:.3f} (p = {p_mass:.4f})")
        if p_mass < 0.05:
            if corr_mass < 0:
                print(f"  ✓ SIGNIFICANT: More massive galaxies fit BETTER")
            else:
                print(f"  ✓ SIGNIFICANT: More massive galaxies fit WORSE")
        else:
            print(f"  No significant correlation")
        print()

    # Correlation 3: χ²_red vs Luminosity (proxy for mass)
    if np.sum(valid_lum) > 10:
        corr_lum, p_lum = spearmanr(log_lums[valid_lum], chi2_reds[valid_lum])
        print(f"χ²_red vs log(Luminosity estimate):")
        print(f"  Spearman ρ = {corr_lum:.3f} (p = {p_lum:.4f})")
        if p_lum < 0.05:
            if corr_lum < 0:
                print(f"  ✓ SIGNIFICANT: Brighter galaxies fit BETTER")
            else:
                print(f"  ✓ SIGNIFICANT: Brighter galaxies fit WORSE")
        else:
            print(f"  No significant correlation")
        print()

    # Correlation 4: χ²_red vs R_last (galaxy size)
    if np.sum(valid_rlast) > 10:
        corr_rlast, p_rlast = spearmanr(r_lasts[valid_rlast], chi2_reds[valid_rlast])
        print(f"χ²_red vs R_last (galaxy size):")
        print(f"  Spearman ρ = {corr_rlast:.3f} (p = {p_rlast:.4f})")
        if p_rlast < 0.05:
            if corr_rlast < 0:
                print(f"  ✓ SIGNIFICANT: Larger galaxies fit BETTER")
            else:
                print(f"  ✓ SIGNIFICANT: Larger galaxies fit WORSE")
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

    # Morphology comparison removed - using string categories now
    # Could add morphology distribution analysis later

    # Create visualization
    create_correlation_plots(galaxy_results)

    # Save results
    save_correlation_results(galaxy_results, {
        'chi2_vs_vmax': (corr_vmax, p_vmax) if np.sum(valid_vmax) > 10 else (np.nan, np.nan),
        'chi2_vs_mass': (corr_mass, p_mass) if np.sum(valid_mass) > 10 else (np.nan, np.nan),
        'chi2_vs_luminosity': (corr_lum, p_lum) if np.sum(valid_lum) > 10 else (np.nan, np.nan),
        'chi2_vs_rlast': (corr_rlast, p_rlast) if np.sum(valid_rlast) > 10 else (np.nan, np.nan),
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
    v_maxes = np.array([g['v_max'] for g in galaxy_results])
    alphas = np.array([g['alpha'] for g in galaxy_results])

    # Plot 1: χ²_red vs Luminosity
    ax = axes[0, 0]
    valid = ~np.isnan(log_lums) & ~np.isnan(chi2_reds)
    ax.scatter(log_lums[valid], chi2_reds[valid], alpha=0.6, s=30)
    ax.axhline(y=2, color='g', linestyle='--', label='χ²_red = 2 (good fit)')
    ax.axhline(y=5, color='orange', linestyle='--', label='χ²_red = 5 (acceptable)')
    ax.set_xlabel('log(Luminosity estimate)')
    ax.set_ylabel('χ²_red')
    ax.set_yscale('log')
    ax.set_title('Fit Quality vs Galaxy Luminosity')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: χ²_red vs V_max
    ax = axes[0, 1]
    valid = ~np.isnan(v_maxes) & ~np.isnan(chi2_reds)
    ax.scatter(v_maxes[valid], chi2_reds[valid], alpha=0.6, s=30)
    ax.axhline(y=2, color='g', linestyle='--')
    ax.axhline(y=5, color='orange', linestyle='--')
    ax.set_xlabel('V_max (km/s)')
    ax.set_ylabel('χ²_red')
    ax.set_yscale('log')
    ax.set_title('Fit Quality vs Maximum Velocity')
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
            key: {
                'rho': value[0] if isinstance(value, tuple) else value,
                'p_value': value[1] if isinstance(value, tuple) and len(value) > 1 else None
            }
            for key, value in correlations.items()
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
