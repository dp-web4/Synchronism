#!/usr/bin/env python3
"""
Session #87: Radial V/V_bar Analysis

The CORRECT test for Synchronism vs MOND:

Synchronism predicts:
- G_eff(r) = G / C(ρ(r))
- C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))
- At each radius, LOCAL density determines C
- Inner disk (high ρ): C → 1, V_obs/V_bar → 1
- Outer disk (low ρ): C < 1, V_obs/V_bar > 1

MOND predicts:
- g = g_bar × μ(g/a₀)
- Transition depends on ACCELERATION, not density
- V_obs/V_bar scales with g/a₀

This analysis tests whether V_obs/V_bar correlates better with:
1. Local surface brightness (Synchronism proxy)
2. Local acceleration (MOND proxy)

Author: CBP Autonomous Synchronism Research
Date: December 5, 2025
"""

import numpy as np
import json
import os
from collections import defaultdict

# Constants
ML_RATIO = 0.5  # M/L ratio at [3.6] micron
G = 4.302e-6  # kpc (km/s)² / M_sun
A0_MOND = 1.2e-10  # m/s² = 3.9e3 (km/s)² / kpc

def parse_sparc_mass_models(filepath):
    """Parse SPARC mass models file to get radial profiles."""
    galaxies = defaultdict(list)

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find data start
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('-----'):
            data_start = i + 1
            break

    for line in lines[data_start:]:
        line = line.strip()
        if len(line) < 50:
            continue

        try:
            parts = line.split()
            if len(parts) < 9:
                continue

            name = parts[0]
            D = float(parts[1])  # Mpc
            R = float(parts[2])  # kpc
            Vobs = float(parts[3])  # km/s
            e_Vobs = float(parts[4])  # km/s
            Vgas = float(parts[5])  # km/s
            Vdisk = float(parts[6])  # km/s (for M/L=1)
            Vbul = float(parts[7])  # km/s (for M/L=1)
            SBdisk = float(parts[8])  # L_sun/pc²
            SBbul = float(parts[9]) if len(parts) > 9 else 0.0

            galaxies[name].append({
                'R': R,
                'Vobs': Vobs,
                'e_Vobs': e_Vobs,
                'Vgas': Vgas,
                'Vdisk': Vdisk,
                'Vbul': Vbul,
                'SBdisk': SBdisk,
                'SBbul': SBbul,
                'D': D
            })
        except (ValueError, IndexError):
            continue

    return dict(galaxies)


def compute_radial_quantities(galaxy_data, ml_ratio=0.5):
    """
    Compute V_bar and derived quantities at each radius.

    Returns list of dicts with:
    - R: radius (kpc)
    - Vobs: observed velocity
    - Vbar: baryonic velocity
    - ratio: Vobs/Vbar (DM enhancement)
    - SB: local surface brightness
    - g_bar: Newtonian acceleration from baryons
    """
    results = []

    for point in galaxy_data:
        # Baryonic velocity: V_bar² = V_gas² + ML×V_disk² + ML×V_bul²
        Vbar_sq = (point['Vgas']**2 +
                   ml_ratio * point['Vdisk']**2 +
                   ml_ratio * point['Vbul']**2)

        if Vbar_sq <= 0:
            continue

        Vbar = np.sqrt(Vbar_sq)

        if Vbar < 1.0 or point['Vobs'] < 1.0:
            continue

        ratio = point['Vobs'] / Vbar

        # Newtonian acceleration: g_bar = V_bar² / R
        R = point['R']
        if R <= 0:
            continue

        g_bar = Vbar_sq / R  # (km/s)² / kpc

        # Convert to m/s² for comparison with a₀
        g_bar_ms2 = g_bar * 1e6 / (3.086e19)  # (km/s)² to m²/s², kpc to m

        results.append({
            'R': R,
            'Vobs': point['Vobs'],
            'Vbar': Vbar,
            'ratio': ratio,
            'SB': point['SBdisk'],
            'g_bar': g_bar,
            'g_bar_ms2': g_bar_ms2,
            'g_ratio': g_bar_ms2 / A0_MOND  # g/a₀
        })

    return results


def analyze_correlation(data_points, x_key, y_key='ratio'):
    """
    Compute correlation between x and y (log-log).
    """
    x_vals = []
    y_vals = []

    for p in data_points:
        x = p[x_key]
        y = p[y_key]
        if x > 0 and y > 0:
            x_vals.append(np.log10(x))
            y_vals.append(np.log10(y))

    if len(x_vals) < 10:
        return None

    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)

    # Pearson correlation
    r = np.corrcoef(x_vals, y_vals)[0, 1]

    # Linear fit
    slope, intercept = np.polyfit(x_vals, y_vals, 1)

    return {
        'n_points': len(x_vals),
        'pearson_r': r,
        'slope': slope,
        'intercept': intercept,
        'x_range': [float(np.min(x_vals)), float(np.max(x_vals))],
        'y_range': [float(np.min(y_vals)), float(np.max(y_vals))]
    }


def bin_by_quantity(data_points, key, n_bins=10):
    """
    Bin data points by a quantity and compute statistics in each bin.
    """
    values = [p[key] for p in data_points if p[key] > 0]
    if len(values) < n_bins * 3:
        n_bins = max(3, len(values) // 3)

    log_values = np.log10(values)
    bin_edges = np.linspace(np.min(log_values), np.max(log_values), n_bins + 1)

    bins = []
    for i in range(n_bins):
        bin_points = [p for p in data_points
                     if p[key] > 0 and
                     bin_edges[i] <= np.log10(p[key]) < bin_edges[i+1]]

        if len(bin_points) < 3:
            continue

        ratios = [p['ratio'] for p in bin_points]

        bins.append({
            'bin_center': 10**((bin_edges[i] + bin_edges[i+1])/2),
            'bin_low': 10**bin_edges[i],
            'bin_high': 10**bin_edges[i+1],
            'n_points': len(bin_points),
            'mean_ratio': np.mean(ratios),
            'std_ratio': np.std(ratios),
            'median_ratio': np.median(ratios)
        })

    return bins


def main():
    print("=" * 70)
    print("SESSION #87: RADIAL V/V_bar ANALYSIS")
    print("=" * 70)
    print()
    print("This tests the CORRECT Synchronism prediction:")
    print("  C depends on LOCAL density at each radius")
    print("  V_obs/V_bar should correlate with local surface brightness")
    print()
    print("MOND prediction:")
    print("  V_obs/V_bar should correlate with g/a₀ (acceleration)")
    print()

    # Load data
    data_file = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/sparc_real_data/MassModels_Lelli2016c.mrt'
    galaxies = parse_sparc_mass_models(data_file)

    print(f"Loaded {len(galaxies)} galaxies with rotation curve data")
    print()

    # Collect all radial points
    all_points = []
    for name, data in galaxies.items():
        points = compute_radial_quantities(data, ML_RATIO)
        for p in points:
            p['galaxy'] = name
        all_points.extend(points)

    print(f"Total radial data points: {len(all_points)}")
    print()

    # Analyze correlations
    print("=" * 70)
    print("CORRELATION ANALYSIS")
    print("=" * 70)
    print()

    # 1. V_obs/V_bar vs Surface Brightness (Synchronism proxy)
    print("1. V_obs/V_bar vs Surface Brightness (SB)")
    print("   Synchronism: NEGATIVE correlation expected")
    print("   (High SB → high ρ → high C → V/V_bar → 1)")
    sb_corr = analyze_correlation(all_points, 'SB', 'ratio')
    if sb_corr:
        print(f"   N = {sb_corr['n_points']}")
        print(f"   Pearson r = {sb_corr['pearson_r']:.3f}")
        print(f"   Slope = {sb_corr['slope']:.3f}")
        if sb_corr['pearson_r'] < -0.1:
            print("   → NEGATIVE correlation: CONSISTENT with Synchronism")
        elif sb_corr['pearson_r'] > 0.1:
            print("   → POSITIVE correlation: CONTRADICTS Synchronism")
        else:
            print("   → No significant correlation")
    print()

    # 2. V_obs/V_bar vs g/a₀ (MOND proxy)
    print("2. V_obs/V_bar vs g/a₀ (acceleration ratio)")
    print("   MOND: NEGATIVE correlation expected")
    print("   (High g → μ → 1 → V/V_bar → 1)")
    g_corr = analyze_correlation(all_points, 'g_ratio', 'ratio')
    if g_corr:
        print(f"   N = {g_corr['n_points']}")
        print(f"   Pearson r = {g_corr['pearson_r']:.3f}")
        print(f"   Slope = {g_corr['slope']:.3f}")
        if g_corr['pearson_r'] < -0.1:
            print("   → NEGATIVE correlation: CONSISTENT with MOND")
        elif g_corr['pearson_r'] > 0.1:
            print("   → POSITIVE correlation: CONTRADICTS MOND")
        else:
            print("   → No significant correlation")
    print()

    # 3. V_obs/V_bar vs Radius (normalized by R_eff?)
    print("3. V_obs/V_bar vs Radius")
    print("   Both theories: ratio increases with R (outer disk more DM-dominated)")
    r_corr = analyze_correlation(all_points, 'R', 'ratio')
    if r_corr:
        print(f"   N = {r_corr['n_points']}")
        print(f"   Pearson r = {r_corr['pearson_r']:.3f}")
        print(f"   Slope = {r_corr['slope']:.3f}")
    print()

    # Compare SB vs g correlation strength
    print("=" * 70)
    print("KEY COMPARISON: Which correlates better with V/V_bar?")
    print("=" * 70)
    print()

    if sb_corr and g_corr:
        r_sb = abs(sb_corr['pearson_r'])
        r_g = abs(g_corr['pearson_r'])

        print(f"  |r| for Surface Brightness: {r_sb:.3f}")
        print(f"  |r| for g/a₀:               {r_g:.3f}")
        print()

        if r_sb > r_g:
            print("  RESULT: Surface brightness is STRONGER correlate")
            print("  → FAVORS Synchronism (density-based coherence)")
        elif r_g > r_sb:
            print("  RESULT: g/a₀ is STRONGER correlate")
            print("  → FAVORS MOND (acceleration-based transition)")
        else:
            print("  RESULT: Similar correlation strength")
            print("  → Cannot discriminate between theories")
    print()

    # Binned analysis
    print("=" * 70)
    print("BINNED ANALYSIS: V/V_bar in SB bins")
    print("=" * 70)
    print()

    sb_bins = bin_by_quantity(all_points, 'SB', n_bins=8)

    print("  SB range (L/pc²)    N    Mean V/V_bar   Std")
    print("-" * 60)
    for b in sb_bins:
        print(f"  {b['bin_low']:8.1f} - {b['bin_high']:8.1f}  {b['n_points']:4d}     {b['mean_ratio']:.3f}    ±{b['std_ratio']:.3f}")
    print()

    # Check trend
    if len(sb_bins) >= 3:
        low_sb_ratio = sb_bins[0]['mean_ratio']
        high_sb_ratio = sb_bins[-1]['mean_ratio']
        trend = high_sb_ratio - low_sb_ratio

        print(f"  Trend: High SB - Low SB = {trend:+.3f}")
        if trend < 0:
            print("  → High SB has LOWER V/V_bar (closer to 1)")
            print("  → CONSISTENT with Synchronism: high ρ → high C → lower enhancement")
        else:
            print("  → High SB has HIGHER V/V_bar")
            print("  → CONTRADICTS Synchronism expectation")
    print()

    # Binned analysis for g/a₀
    print("=" * 70)
    print("BINNED ANALYSIS: V/V_bar in g/a₀ bins")
    print("=" * 70)
    print()

    g_bins = bin_by_quantity(all_points, 'g_ratio', n_bins=8)

    print("  g/a₀ range           N    Mean V/V_bar   Std")
    print("-" * 60)
    for b in g_bins:
        print(f"  {b['bin_low']:8.3f} - {b['bin_high']:8.3f}  {b['n_points']:4d}     {b['mean_ratio']:.3f}    ±{b['std_ratio']:.3f}")
    print()

    if len(g_bins) >= 3:
        low_g_ratio = g_bins[0]['mean_ratio']
        high_g_ratio = g_bins[-1]['mean_ratio']
        trend = high_g_ratio - low_g_ratio

        print(f"  Trend: High g/a₀ - Low g/a₀ = {trend:+.3f}")
        if trend < 0:
            print("  → High g/a₀ has LOWER V/V_bar (closer to 1)")
            print("  → CONSISTENT with MOND: Newtonian regime at high g")
    print()

    # Synchronism-specific analysis: C(ρ) prediction
    print("=" * 70)
    print("SYNCHRONISM PREDICTION: C(ρ) from V/V_bar")
    print("=" * 70)
    print()
    print("If G_eff = G/C, then V_obs² = V_bar² / C")
    print("So C = (V_bar/V_obs)² = 1/(ratio)²")
    print()

    # Compute implied C values
    c_values = []
    sb_values = []
    for p in all_points:
        if p['ratio'] > 0 and p['SB'] > 0:
            C = 1.0 / (p['ratio']**2)
            c_values.append(C)
            sb_values.append(p['SB'])

    c_values = np.array(c_values)
    sb_values = np.array(sb_values)

    print(f"  Implied C range: {np.min(c_values):.3f} to {np.max(c_values):.3f}")
    print(f"  Mean C: {np.mean(c_values):.3f}")
    print(f"  Median C: {np.median(c_values):.3f}")
    print()

    # C vs SB correlation
    log_C = np.log10(c_values[c_values > 0])
    log_SB = np.log10(sb_values[c_values > 0])

    r_C_SB = np.corrcoef(log_SB, log_C)[0, 1]
    print(f"  Correlation log(C) vs log(SB): r = {r_C_SB:.3f}")
    print()

    if r_C_SB > 0.1:
        print("  → POSITIVE: Higher SB → higher C")
        print("  → CONSISTENT with Synchronism: ρ determines C")
    elif r_C_SB < -0.1:
        print("  → NEGATIVE: Higher SB → lower C")
        print("  → CONTRADICTS Synchronism")
    print()

    # Summary
    print("=" * 70)
    print("SESSION #87 SUMMARY")
    print("=" * 70)
    print()

    results = {
        'n_galaxies': len(galaxies),
        'n_data_points': len(all_points),
        'correlations': {
            'ratio_vs_SB': sb_corr,
            'ratio_vs_g': g_corr,
            'ratio_vs_R': r_corr
        },
        'sb_bins': sb_bins,
        'g_bins': g_bins,
        'implied_C': {
            'min': float(np.min(c_values)),
            'max': float(np.max(c_values)),
            'mean': float(np.mean(c_values)),
            'median': float(np.median(c_values)),
            'correlation_with_SB': float(r_C_SB)
        }
    }

    # Key findings
    print("KEY FINDINGS:")
    print()

    if sb_corr and g_corr:
        print(f"1. V/V_bar correlation with SB: r = {sb_corr['pearson_r']:.3f}")
        print(f"2. V/V_bar correlation with g/a₀: r = {g_corr['pearson_r']:.3f}")
        print(f"3. Implied C correlation with SB: r = {r_C_SB:.3f}")
        print()

        # Interpretation
        if sb_corr['pearson_r'] < -0.1 and r_C_SB > 0.1:
            print("INTERPRETATION: CONSISTENT with Synchronism")
            print("  - V/V_bar decreases with SB (negative correlation)")
            print("  - C increases with SB (positive correlation)")
            print("  - This is exactly what C(ρ) predicts!")
            results['interpretation'] = 'consistent_with_synchronism'
        elif sb_corr['pearson_r'] > 0.1 or r_C_SB < -0.1:
            print("INTERPRETATION: CONTRADICTS Synchronism")
            print("  - Wrong correlation direction")
            results['interpretation'] = 'contradicts_synchronism'
        else:
            print("INTERPRETATION: WEAK or AMBIGUOUS")
            print("  - Correlations not strong enough to discriminate")
            results['interpretation'] = 'ambiguous'
    print()

    # Save results
    results_dir = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/results'
    os.makedirs(results_dir, exist_ok=True)

    output_file = os.path.join(results_dir, 'session87_radial_analysis.json')

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Results saved to: {output_file}")


if __name__ == '__main__':
    main()
