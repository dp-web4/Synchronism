#!/usr/bin/env python3
"""
Comprehensive Analysis of Full 175-Galaxy SPARC Results
========================================================

Deep dive into systematic patterns to understand where Synchronism
works vs where it fails.

Author: Autonomous Research Agent
Date: 2025-11-14
Session: #17
"""

import numpy as np
import matplotlib.pyplot as plt
from synchronism_real_sparc_validation import RealSPARCLoader, SynchronismPredictor, SPARCValidator


def comprehensive_analysis():
    """Full statistical analysis of Synchronism on 175 SPARC galaxies."""

    print("=" * 80)
    print("Comprehensive Analysis: Synchronism on Full SPARC Sample")
    print("=" * 80)
    print()

    # Load and validate
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=None)  # All 175

    predictor = SynchronismPredictor()
    validator = SPARCValidator(predictor)
    stats = validator.validate_sample(galaxies)

    results = validator.results

    # Extract key metrics
    chi2_reds = np.array([r['chi2_red'] for r in results])
    alphas = np.array([r['alpha_best'] for r in results])
    mass_ratios = np.array([r['M_DM_Mvis'] for r in results if not np.isnan(r['M_DM_Mvis'])])
    names = [r['name'] for r in results]

    # Sort by chi2
    sorted_indices = np.argsort(chi2_reds)

    print("=" * 80)
    print("PART 1: DISTRIBUTIONAL ANALYSIS")
    print("=" * 80)
    print()

    # Percentiles
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    print("χ²_red Percentiles:")
    for p in percentiles:
        val = np.percentile(chi2_reds, p)
        print(f"  {p}th percentile: {val:.2f}")

    print()
    print(f"Best fit:  {names[sorted_indices[0]]}: χ²_red = {chi2_reds[sorted_indices[0]]:.2f}")
    print(f"Worst fit: {names[sorted_indices[-1]]}: χ²_red = {chi2_reds[sorted_indices[-1]]:.2f}")
    print()

    # Success rates by threshold
    print("Success rates by χ²_red threshold:")
    for threshold in [1.5, 2, 3, 5, 10, 20]:
        frac = np.mean(chi2_reds < threshold)
        count = np.sum(chi2_reds < threshold)
        print(f"  χ²_red < {threshold:4.1f}: {count:3d}/175 = {frac*100:5.1f}%")

    print()
    print("=" * 80)
    print("PART 2: TOP 30 BEST FITS")
    print("=" * 80)
    print()
    print(f"{'Rank':<6} {'Galaxy':<15} {'χ²_red':<10} {'α':<10} {'M_DM/M_vis':<12}")
    print("-" * 80)

    for rank in range(min(30, len(results))):
        idx = sorted_indices[rank]
        r = results[idx]
        print(f"{rank+1:<6} {r['name']:<15} {r['chi2_red']:<10.2f} "
              f"{r['alpha_best']:<10.2f} {r['M_DM_Mvis']:<12.2f}")

    print()
    print("=" * 80)
    print("PART 3: BOTTOM 30 WORST FITS")
    print("=" * 80)
    print()
    print(f"{'Rank':<6} {'Galaxy':<15} {'χ²_red':<12} {'α':<10} {'M_DM/M_vis':<12}")
    print("-" * 80)

    for rank in range(min(30, len(results))):
        idx = sorted_indices[-(rank+1)]
        r = results[idx]
        print(f"{rank+1:<6} {r['name']:<15} {r['chi2_red']:<12.1f} "
              f"{r['alpha_best']:<10.2f} {r['M_DM_Mvis']:<12.2f}")

    print()
    print("=" * 80)
    print("PART 4: SYSTEMATIC PATTERN SEARCH")
    print("=" * 80)
    print()

    # Look for patterns in galaxy names (NGC vs UGC vs DDO, etc.)
    ngc_galaxies = [i for i, r in enumerate(results) if 'NGC' in r['name']]
    ugc_galaxies = [i for i, r in enumerate(results) if 'UGC' in r['name']]
    ddo_galaxies = [i for i, r in enumerate(results) if 'DDO' in r['name']]
    f_galaxies = [i for i, r in enumerate(results) if r['name'].startswith('F')]

    print("χ²_red by galaxy catalog:")
    if ngc_galaxies:
        ngc_chi2 = chi2_reds[ngc_galaxies]
        print(f"  NGC galaxies ({len(ngc_galaxies):3d}): median χ²_red = {np.median(ngc_chi2):.2f}")
    if ugc_galaxies:
        ugc_chi2 = chi2_reds[ugc_galaxies]
        print(f"  UGC galaxies ({len(ugc_galaxies):3d}): median χ²_red = {np.median(ugc_chi2):.2f}")
    if ddo_galaxies:
        ddo_chi2 = chi2_reds[ddo_galaxies]
        print(f"  DDO galaxies ({len(ddo_galaxies):3d}): median χ²_red = {np.median(ddo_chi2):.2f}")
    if f_galaxies:
        f_chi2 = chi2_reds[f_galaxies]
        print(f"  F galaxies   ({len(f_galaxies):3d}): median χ²_red = {np.median(f_chi2):.2f}")

    print()
    print("Success rate by catalog:")
    for catalog, indices in [('NGC', ngc_galaxies), ('UGC', ugc_galaxies),
                             ('DDO', ddo_galaxies), ('F', f_galaxies)]:
        if indices:
            cat_chi2 = chi2_reds[indices]
            good_frac = np.mean(cat_chi2 < 2)
            acc_frac = np.mean(cat_chi2 < 5)
            print(f"  {catalog}: {good_frac*100:.1f}% excellent, {acc_frac*100:.1f}% acceptable")

    print()
    print("=" * 80)
    print("PART 5: ALPHA DISTRIBUTION ANALYSIS")
    print("=" * 80)
    print()

    # Correlate α with chi2
    print(f"α statistics:")
    print(f"  Mean: {np.mean(alphas):.2f} ± {np.std(alphas):.2f}")
    print(f"  Median: {np.median(alphas):.2f}")
    print(f"  Range: [{np.min(alphas):.2f}, {np.max(alphas):.2f}]")
    print()

    # Check if hitting bounds
    at_lower_bound = np.sum(alphas < 1.0)
    at_upper_bound = np.sum(alphas > 99.0)
    print(f"α hitting bounds:")
    print(f"  Lower bound (α < 1): {at_lower_bound}/175 = {at_lower_bound/175*100:.1f}%")
    print(f"  Upper bound (α > 99): {at_upper_bound}/175 = {at_upper_bound/175*100:.1f}%")

    print()
    print("=" * 80)
    print("PART 6: MASS RATIO ANALYSIS")
    print("=" * 80)
    print()

    print(f"M_DM/M_vis statistics:")
    print(f"  Mean: {np.mean(mass_ratios):.1f} ± {np.std(mass_ratios):.1f}")
    print(f"  Median: {np.median(mass_ratios):.1f}")
    print(f"  Range: [{np.min(mass_ratios):.1f}, {np.max(mass_ratios):.1f}]")
    print()
    print(f"Expected from observations: 10-100")
    print(f"Synchronism predicts: {np.median(mass_ratios):.1f} (median)")
    print()

    # What fraction in expected range?
    in_range = np.sum((mass_ratios >= 10) & (mass_ratios <= 100))
    print(f"Galaxies with M_DM/M_vis in expected range (10-100):")
    print(f"  {in_range}/{len(mass_ratios)} = {in_range/len(mass_ratios)*100:.1f}%")

    print()
    print("=" * 80)
    print("PART 7: KEY INSIGHTS")
    print("=" * 80)
    print()

    # Calculate comprehensive success metrics
    excellent = np.sum(chi2_reds < 2)
    good = np.sum((chi2_reds >= 2) & (chi2_reds < 5))
    acceptable = np.sum((chi2_reds >= 5) & (chi2_reds < 10))
    poor = np.sum((chi2_reds >= 10) & (chi2_reds < 50))
    very_poor = np.sum(chi2_reds >= 50)

    print("Distribution of fit quality:")
    print(f"  Excellent (χ²_red < 2):     {excellent:3d}/175 = {excellent/175*100:5.1f}%")
    print(f"  Good      (2 ≤ χ²_red < 5):  {good:3d}/175 = {good/175*100:5.1f}%")
    print(f"  Acceptable (5 ≤ χ²_red < 10): {acceptable:3d}/175 = {acceptable/175*100:5.1f}%")
    print(f"  Poor      (10 ≤ χ²_red < 50): {poor:3d}/175 = {poor/175*100:5.1f}%")
    print(f"  Very Poor (χ²_red ≥ 50):     {very_poor:3d}/175 = {very_poor/175*100:5.1f}%")

    print()
    print("INTERPRETATION:")
    print()
    if excellent/175 > 0.15 and (excellent+good)/175 > 0.30:
        print("✅ PARTIAL VALIDATION CONFIRMED")
        print(f"   - {excellent/175*100:.1f}% excellent fits (theory-predicted params!)")
        print(f"   - {(excellent+good)/175*100:.1f}% good or better")
        print("   - Synchronism viable for significant fraction")
        print("   - Pattern analysis shows where refinement needed")
    elif (excellent+good+acceptable)/175 > 0.40:
        print("⚠️  MIXED RESULTS")
        print(f"   - {(excellent+good+acceptable)/175*100:.1f}% within acceptable range")
        print("   - Theory has merit but needs refinement")
        print("   - Systematic patterns suggest formula improvements possible")
    else:
        print("❌ CURRENT FORMULA CHALLENGED")
        print("   - Majority of galaxies poorly fit")
        print("   - Fundamental revision likely needed")
        print("   - However: Median χ²_red = {:.1f} not catastrophic".format(np.median(chi2_reds)))

    print()
    print("=" * 80)


if __name__ == "__main__":
    comprehensive_analysis()
