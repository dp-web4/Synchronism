#!/usr/bin/env python3
"""
Session #82 Track C: SPARC Environment Proxy Analysis

Explores whether SPARC data contains any proxy for environment that could
be used to test the Synchronism void prediction.

Key insight from Session #80: LSB galaxies may correlate with low-density
environments. If so, we could test for BTFR offset in SPARC's LSB subsample.

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
Session: #82 - Environment Proxy Investigation
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime

# SPARC data directory
SPARC_DIR = Path(__file__).parent / 'sparc_real_data' / 'galaxies'


def load_sparc_galaxies():
    """Load all SPARC galaxies and extract key properties."""
    galaxies = []

    for dat_file in SPARC_DIR.glob('*.dat'):
        try:
            # Read header for distance
            with open(dat_file, 'r') as f:
                lines = f.readlines()

            # Parse header
            name = None
            distance = None
            for line in lines:
                if line.startswith('# Galaxy:'):
                    name = line.split(':')[1].strip()
                elif line.startswith('# Distance:'):
                    distance = float(line.split(':')[1].strip().split()[0])

            if name is None:
                name = dat_file.stem

            # Parse data
            data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
            if not data_lines:
                continue

            radii = []
            v_obs = []
            sb_disk = []

            for line in data_lines:
                parts = line.split()
                if len(parts) >= 7:
                    radii.append(float(parts[0]))
                    v_obs.append(float(parts[1]))
                    sb_disk.append(float(parts[6]))

            if not radii:
                continue

            radii = np.array(radii)
            v_obs = np.array(v_obs)
            sb_disk = np.array(sb_disk)

            # Compute galaxy properties
            v_max = np.max(v_obs)
            r_max = radii[-1]

            # Central surface brightness (first point or extrapolation)
            sb_0 = sb_disk[0] if len(sb_disk) > 0 else np.nan

            # Mean surface brightness
            sb_mean = np.mean(sb_disk[sb_disk > 0]) if np.any(sb_disk > 0) else np.nan

            galaxies.append({
                'name': name,
                'distance': distance,
                'v_max': v_max,
                'r_max': r_max,
                'sb_0': sb_0,
                'sb_mean': sb_mean,
                'n_points': len(radii)
            })

        except Exception as e:
            print(f"Error loading {dat_file}: {e}")
            continue

    return galaxies


def classify_by_surface_brightness(galaxies):
    """
    Classify galaxies as HSB or LSB based on central surface brightness.

    Standard LSB definition: μ_0 > 22.5 mag/arcsec² in B-band
    SPARC uses 3.6μm, so we use L_sun/pc² directly.

    Empirically: LSB galaxies have SB_0 < ~50 L_sun/pc² at 3.6μm
    """

    sb_values = [g['sb_0'] for g in galaxies if not np.isnan(g['sb_0'])]
    sb_median = np.median(sb_values)

    # Use median as dividing line (more robust than fixed threshold)
    lsb = [g for g in galaxies if not np.isnan(g['sb_0']) and g['sb_0'] < sb_median]
    hsb = [g for g in galaxies if not np.isnan(g['sb_0']) and g['sb_0'] >= sb_median]

    return hsb, lsb, sb_median


def compute_btfr_by_type(hsb, lsb):
    """
    Compute BTFR offset between HSB and LSB galaxies.

    If Synchronism is correct and LSB correlates with void environment,
    we expect LSB galaxies to have higher V at fixed mass.
    """

    # For SPARC, we don't have M_bar directly, but v_max is a proxy
    # BTFR: M_bar ∝ V^4, so log(M) ∝ 4 × log(V)
    # At fixed radius, LSB galaxies have lower baryonic density

    hsb_log_v = [np.log10(g['v_max']) for g in hsb if g['v_max'] > 0]
    lsb_log_v = [np.log10(g['v_max']) for g in lsb if g['v_max'] > 0]

    hsb_log_r = [np.log10(g['r_max']) for g in hsb if g['r_max'] > 0]
    lsb_log_r = [np.log10(g['r_max']) for g in lsb if g['r_max'] > 0]

    hsb_sb = [g['sb_0'] for g in hsb]
    lsb_sb = [g['sb_0'] for g in lsb]

    results = {
        'hsb': {
            'count': len(hsb),
            'mean_log_v': np.mean(hsb_log_v),
            'std_log_v': np.std(hsb_log_v),
            'mean_log_r': np.mean(hsb_log_r),
            'mean_sb': np.mean(hsb_sb)
        },
        'lsb': {
            'count': len(lsb),
            'mean_log_v': np.mean(lsb_log_v),
            'std_log_v': np.std(lsb_log_v),
            'mean_log_r': np.mean(lsb_log_r),
            'mean_sb': np.mean(lsb_sb)
        }
    }

    # Delta log(V) - but this is NOT at fixed mass
    # LSB galaxies are typically lower mass, so will naturally have lower V
    results['delta_log_v'] = results['lsb']['mean_log_v'] - results['hsb']['mean_log_v']

    return results


def analyze_environment_proxy():
    """
    Main analysis: can surface brightness serve as environment proxy?
    """
    print("=" * 70)
    print("SESSION #82 TRACK C: SPARC ENVIRONMENT PROXY ANALYSIS")
    print("=" * 70)

    # Load data
    print("\n1. Loading SPARC galaxies...")
    galaxies = load_sparc_galaxies()
    print(f"   Loaded {len(galaxies)} galaxies")

    # Classify by SB
    print("\n2. Classifying by surface brightness...")
    hsb, lsb, sb_median = classify_by_surface_brightness(galaxies)
    print(f"   HSB (SB > {sb_median:.1f} L_sun/pc²): {len(hsb)} galaxies")
    print(f"   LSB (SB < {sb_median:.1f} L_sun/pc²): {len(lsb)} galaxies")

    # Compute BTFR comparison
    print("\n3. Computing BTFR comparison...")
    btfr = compute_btfr_by_type(hsb, lsb)

    print(f"\n   HSB galaxies:")
    print(f"     N = {btfr['hsb']['count']}")
    print(f"     <log V> = {btfr['hsb']['mean_log_v']:.3f} ± {btfr['hsb']['std_log_v']:.3f}")
    print(f"     <SB_0> = {btfr['hsb']['mean_sb']:.1f} L_sun/pc²")

    print(f"\n   LSB galaxies:")
    print(f"     N = {btfr['lsb']['count']}")
    print(f"     <log V> = {btfr['lsb']['mean_log_v']:.3f} ± {btfr['lsb']['std_log_v']:.3f}")
    print(f"     <SB_0> = {btfr['lsb']['mean_sb']:.1f} L_sun/pc²")

    print(f"\n   Δlog(V) = {btfr['delta_log_v']:.3f} dex (LSB - HSB)")

    # Analysis
    print("\n4. Analysis")
    print("=" * 70)

    print("""
    KEY INSIGHT: Surface brightness is NOT a clean environment proxy.

    Problem 1: LSB galaxies have lower mass on average
    - This means they naturally have lower V
    - Cannot directly compare V without matching mass

    Problem 2: LSB ≠ Void
    - LSB galaxies can exist in any environment
    - Environment affects formation, not just current SB

    Problem 3: Need mass-matched samples
    - To test BTFR offset, need galaxies with SAME M_bar
    - Then compare V between environments

    CONCLUSION: SPARC data CANNOT directly test void prediction.

    What we NEED:
    1. Environment classification for each galaxy (δ values)
    2. Cross-match with SDSS or void catalogs
    3. Or wait for BIG-SPARC with 4000 galaxies
    """)

    # Alternative: McGaugh's published BTFR
    print("\n5. Alternative: McGaugh's BTFR Analysis")
    print("=" * 70)
    print("""
    McGaugh et al. have published extensive BTFR analyses:
    - 2000: First systematic BTFR paper
    - 2012: Extended sample with dwarfs
    - 2016: SPARC BTFR (used for our validation)

    Key finding: BTFR scatter is ~0.1 dex

    If void offset were 0.28 dex (extreme voids):
    - Would show as outliers HIGH in V at fixed M
    - These should be visible in BTFR residuals
    - But: only if voids are ~10% of sample

    Prediction: Check BTFR residuals for HIGH-V outliers
    - If these correlate with low-density environment → Support for Synchronism
    - If random distribution → Against Synchronism
    """)

    return {
        'galaxies_loaded': len(galaxies),
        'hsb_count': len(hsb),
        'lsb_count': len(lsb),
        'sb_median': sb_median,
        'btfr_comparison': btfr,
        'conclusion': 'SPARC lacks environment data - cannot directly test void prediction',
        'alternative': 'Cross-match SPARC with void catalogs or use ALFALFA+SDSS'
    }


def main():
    results = analyze_environment_proxy()

    # Save results
    output_path = Path(__file__).parent / 'results' / 'session82_environment_proxy.json'
    output_path.parent.mkdir(exist_ok=True)

    output = {
        'session': 82,
        'track': 'C',
        'title': 'SPARC Environment Proxy Analysis',
        'date': datetime.now().isoformat(),
        'conclusion': results['conclusion'],
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nResults saved to: {output_path}")
    print("\n" + "=" * 70)
    print("SESSION #82 TRACK C COMPLETE")
    print("=" * 70)

    return results


if __name__ == '__main__':
    main()
