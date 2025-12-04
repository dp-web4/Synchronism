#!/usr/bin/env python3
"""
Session #84: ALFALFA α.100 Catalog Download

Downloads the ALFALFA HI source catalog from VizieR for void BTFR analysis.

Data source: J/ApJ/861/49 (Haynes+ 2018)
- 31,502 HI sources
- Key columns: RA, Dec, W50, W20, logMHI, Vhelio

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
Session: #84
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime

try:
    from astroquery.vizier import Vizier
    from astropy.table import Table
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False
    print("WARNING: astroquery not available, will use alternative download method")


def download_alfalfa_catalog():
    """
    Download ALFALFA α.100 catalog from VizieR.
    """
    print("=" * 70)
    print("SESSION #84: ALFALFA α.100 CATALOG DOWNLOAD")
    print("=" * 70)

    if not ASTROQUERY_AVAILABLE:
        print("\nERROR: astroquery not available")
        return None

    # Configure Vizier for full catalog download
    print("\nConfiguring VizieR query...")
    v = Vizier(
        columns=['RAJ2000', 'DEJ2000', 'W50', 'W20', 'logMHI', 'Vhelio', 'SNR', 'Dist', 'AGCNr'],
        row_limit=-1  # No limit - get all rows
    )

    print("Downloading ALFALFA α.100 (J/ApJ/861/49)...")
    print("This may take a few minutes for 31,502 sources...")

    try:
        # Query the catalog
        catalogs = v.get_catalogs('J/ApJ/861/49')

        if catalogs is None or len(catalogs) == 0:
            print("ERROR: No data returned from VizieR")
            return None

        # Get the main table
        alfalfa = catalogs[0]
        print(f"\nDownloaded {len(alfalfa)} sources")

        return alfalfa

    except Exception as e:
        print(f"ERROR downloading catalog: {e}")
        return None


def analyze_catalog(alfalfa):
    """
    Analyze the downloaded ALFALFA catalog.
    """
    print("\n" + "=" * 70)
    print("CATALOG ANALYSIS")
    print("=" * 70)

    print(f"\nTotal sources: {len(alfalfa)}")
    print(f"\nColumns available: {alfalfa.colnames}")

    # Statistics
    print("\n--- Key Statistics ---")

    if 'W50' in alfalfa.colnames:
        w50 = alfalfa['W50']
        w50_valid = w50[~np.isnan(w50)]
        print(f"\nW50 (velocity width at 50%):")
        print(f"  N valid: {len(w50_valid)}")
        print(f"  Min: {np.min(w50_valid):.1f} km/s")
        print(f"  Max: {np.max(w50_valid):.1f} km/s")
        print(f"  Median: {np.median(w50_valid):.1f} km/s")

    if 'logMHI' in alfalfa.colnames:
        logm = alfalfa['logMHI']
        logm_valid = logm[~np.isnan(logm)]
        print(f"\nlog(M_HI) [M_sun]:")
        print(f"  N valid: {len(logm_valid)}")
        print(f"  Min: {np.min(logm_valid):.2f}")
        print(f"  Max: {np.max(logm_valid):.2f}")
        print(f"  Median: {np.median(logm_valid):.2f}")

    if 'Vhelio' in alfalfa.colnames:
        vhel = alfalfa['Vhelio']
        vhel_valid = vhel[~np.isnan(vhel)]
        print(f"\nV_helio (heliocentric velocity):")
        print(f"  N valid: {len(vhel_valid)}")
        print(f"  Min: {np.min(vhel_valid):.0f} km/s")
        print(f"  Max: {np.max(vhel_valid):.0f} km/s")
        print(f"  Median: {np.median(vhel_valid):.0f} km/s")

    if 'SNR' in alfalfa.colnames:
        snr = alfalfa['SNR']
        snr_valid = snr[~np.isnan(snr)]
        print(f"\nSNR (signal-to-noise ratio):")
        print(f"  N valid: {len(snr_valid)}")
        print(f"  Min: {np.min(snr_valid):.1f}")
        print(f"  Max: {np.max(snr_valid):.1f}")
        print(f"  Median: {np.median(snr_valid):.1f}")
        print(f"  N with SNR > 10: {np.sum(snr_valid > 10)}")

    # Sky coverage
    if 'RAJ2000' in alfalfa.colnames and 'DEJ2000' in alfalfa.colnames:
        try:
            ra = np.array([float(x) if x else np.nan for x in alfalfa['RAJ2000']])
            dec = np.array([float(x) if x else np.nan for x in alfalfa['DEJ2000']])
            print(f"\nSky coverage:")
            print(f"  RA range: {np.nanmin(ra):.1f}° to {np.nanmax(ra):.1f}°")
            print(f"  Dec range: {np.nanmin(dec):.1f}° to {np.nanmax(dec):.1f}°")
        except Exception as e:
            print(f"\nSky coverage: Could not parse coordinates ({e})")

    return {
        'n_sources': len(alfalfa),
        'columns': list(alfalfa.colnames)
    }


def apply_quality_cuts(alfalfa):
    """
    Apply quality cuts to ALFALFA catalog for BTFR analysis.
    """
    print("\n" + "=" * 70)
    print("QUALITY CUTS FOR BTFR ANALYSIS")
    print("=" * 70)

    n_initial = len(alfalfa)
    print(f"\nInitial sources: {n_initial}")

    # Create mask for quality cuts
    mask = np.ones(n_initial, dtype=bool)

    # Cut 1: Valid W50
    if 'W50' in alfalfa.colnames:
        w50_valid = ~np.isnan(alfalfa['W50']) & (alfalfa['W50'] > 0)
        mask &= w50_valid
        print(f"After W50 > 0 cut: {np.sum(mask)} ({np.sum(mask)/n_initial*100:.1f}%)")

    # Cut 2: Valid HI mass
    if 'logMHI' in alfalfa.colnames:
        mhi_valid = ~np.isnan(alfalfa['logMHI'])
        mask &= mhi_valid
        print(f"After logMHI valid cut: {np.sum(mask)} ({np.sum(mask)/n_initial*100:.1f}%)")

    # Cut 3: SNR > 10 for reliable measurements
    if 'SNR' in alfalfa.colnames:
        snr_cut = alfalfa['SNR'] > 10
        mask &= snr_cut
        print(f"After SNR > 10 cut: {np.sum(mask)} ({np.sum(mask)/n_initial*100:.1f}%)")

    # Cut 4: W50 > 50 km/s to avoid resolution issues
    if 'W50' in alfalfa.colnames:
        w50_cut = alfalfa['W50'] > 50
        mask &= w50_cut
        print(f"After W50 > 50 km/s cut: {np.sum(mask)} ({np.sum(mask)/n_initial*100:.1f}%)")

    n_final = np.sum(mask)
    print(f"\nFinal sample after all cuts: {n_final} ({n_final/n_initial*100:.1f}%)")

    return alfalfa[mask], mask


def save_catalog(alfalfa, quality_mask=None):
    """
    Save the ALFALFA catalog to local files.
    """
    output_dir = Path(__file__).parent / 'alfalfa_data'
    output_dir.mkdir(exist_ok=True)

    # Save full catalog
    full_path = output_dir / 'alfalfa_full.csv'

    # Convert to simple arrays for saving
    data = {}
    for col in alfalfa.colnames:
        try:
            data[col] = np.array(alfalfa[col])
        except:
            pass

    # Save as CSV manually (more portable)
    with open(full_path, 'w') as f:
        # Header
        cols = list(data.keys())
        f.write(','.join(cols) + '\n')
        # Data
        for i in range(len(alfalfa)):
            row = []
            for col in cols:
                val = data[col][i]
                try:
                    if np.isnan(val):
                        row.append('')
                    else:
                        row.append(str(val))
                except (TypeError, ValueError):
                    # String or other non-numeric type
                    row.append(str(val) if val else '')
            f.write(','.join(row) + '\n')

    print(f"\nFull catalog saved to: {full_path}")

    # Save quality cut version
    if quality_mask is not None:
        quality_alfalfa = alfalfa[quality_mask]
        quality_path = output_dir / 'alfalfa_quality.csv'

        data_quality = {}
        for col in quality_alfalfa.colnames:
            try:
                data_quality[col] = np.array(quality_alfalfa[col])
            except:
                pass

        with open(quality_path, 'w') as f:
            cols = list(data_quality.keys())
            f.write(','.join(cols) + '\n')
            for i in range(len(quality_alfalfa)):
                row = []
                for col in cols:
                    val = data_quality[col][i]
                    try:
                        if np.isnan(val):
                            row.append('')
                        else:
                            row.append(str(val))
                    except (TypeError, ValueError):
                        row.append(str(val) if val else '')
                f.write(','.join(row) + '\n')

        print(f"Quality-cut catalog saved to: {quality_path}")

    return output_dir


def estimate_rotation_velocities(alfalfa):
    """
    Estimate rotation velocities from W50.

    V_rot ≈ W50 / (2 sin i)

    For statistical analysis without individual inclinations,
    use W50 directly or apply average inclination correction.
    """
    print("\n" + "=" * 70)
    print("ROTATION VELOCITY ESTIMATION")
    print("=" * 70)

    if 'W50' not in alfalfa.colnames:
        print("ERROR: W50 column not found")
        return None

    w50 = np.array(alfalfa['W50'])

    # Method 1: Direct W50 (no inclination correction)
    # This adds scatter but no systematic bias for random orientations
    v_rot_uncorrected = w50 / 2

    # Method 2: Statistical correction assuming random inclinations
    # <sin i> ≈ π/4 ≈ 0.785 for randomly oriented disks
    sin_i_mean = np.pi / 4
    v_rot_statistical = w50 / (2 * sin_i_mean)

    print("\nMethod 1: V_rot = W50 / 2 (no inclination correction)")
    print(f"  Adds scatter but no systematic bias")
    print(f"  Median V_rot: {np.nanmedian(v_rot_uncorrected):.1f} km/s")

    print("\nMethod 2: V_rot = W50 / (2 × <sin i>) with <sin i> = π/4")
    print(f"  Statistical correction for random orientations")
    print(f"  Median V_rot: {np.nanmedian(v_rot_statistical):.1f} km/s")

    print("\nFor BTFR analysis, use Method 2 (statistical correction)")
    print("This is standard practice in HI rotation curve literature")

    return {
        'v_rot_uncorrected': v_rot_uncorrected,
        'v_rot_statistical': v_rot_statistical,
        'sin_i_mean': sin_i_mean
    }


def main():
    """Main execution."""
    results = {
        'session': 84,
        'title': 'ALFALFA Catalog Download',
        'date': datetime.now().isoformat(),
        'status': 'STARTED'
    }

    # Download catalog
    alfalfa = download_alfalfa_catalog()

    if alfalfa is None:
        results['status'] = 'FAILED'
        results['error'] = 'Could not download catalog'
        print("\n" + "=" * 70)
        print("SESSION #84 FAILED - Could not download ALFALFA catalog")
        print("=" * 70)
        return results

    # Analyze
    stats = analyze_catalog(alfalfa)
    results['raw_stats'] = stats

    # Quality cuts
    alfalfa_quality, quality_mask = apply_quality_cuts(alfalfa)
    results['quality_n'] = len(alfalfa_quality)

    # Estimate rotation velocities
    v_rot = estimate_rotation_velocities(alfalfa_quality)
    if v_rot:
        results['v_rot_median'] = float(np.nanmedian(v_rot['v_rot_statistical']))

    # Save
    output_dir = save_catalog(alfalfa, quality_mask)
    results['output_dir'] = str(output_dir)
    results['status'] = 'SUCCESS'

    # Summary
    print("\n" + "=" * 70)
    print("SESSION #84 SUMMARY")
    print("=" * 70)
    print(f"""
    ALFALFA α.100 catalog downloaded successfully!

    Raw catalog: {len(alfalfa)} sources
    After quality cuts: {len(alfalfa_quality)} sources

    Quality cuts applied:
    - W50 > 0 (valid velocity width)
    - logMHI valid (HI mass measurement)
    - SNR > 10 (reliable detection)
    - W50 > 50 km/s (avoid resolution issues)

    Files saved to: {output_dir}

    NEXT STEPS:
    1. Download void catalog (Douglass+ 2023)
    2. Cross-match ALFALFA × void catalog by position
    3. Bin by environment δ
    4. Test BTFR offset prediction
    """)

    # Save results
    results_path = Path(__file__).parent / 'results' / 'session84_alfalfa_download.json'
    results_path.parent.mkdir(exist_ok=True)
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"Results saved to: {results_path}")

    return results


if __name__ == '__main__':
    main()
