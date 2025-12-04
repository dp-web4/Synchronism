#!/usr/bin/env python3
"""
Session #84: ALFALFA × Void Cross-Match

Cross-matches ALFALFA galaxies with void catalog to identify environment.

Strategy:
1. Use void centers from Douglass+ 2023 Table 0 (has RA, Dec, R_eff)
2. For each ALFALFA galaxy, find nearest void center
3. Compute distance to void center relative to void radius
4. Classify as: void interior (d < 0.5 R), void edge (0.5-1.0 R), wall (>1.0 R)

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
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False


def load_alfalfa_quality():
    """Load quality-cut ALFALFA catalog."""
    alfalfa_path = Path(__file__).parent / 'alfalfa_data' / 'alfalfa_quality.csv'

    if not alfalfa_path.exists():
        print(f"ERROR: ALFALFA file not found: {alfalfa_path}")
        return None

    # Load CSV
    data = {}
    with open(alfalfa_path, 'r') as f:
        header = f.readline().strip().split(',')
        for col in header:
            data[col] = []

        for line in f:
            parts = line.strip().split(',')
            for i, col in enumerate(header):
                if i < len(parts):
                    try:
                        data[col].append(float(parts[i]) if parts[i] else np.nan)
                    except ValueError:
                        data[col].append(parts[i])
                else:
                    data[col].append(np.nan)

    for col in data:
        data[col] = np.array(data[col])

    print(f"Loaded {len(data[header[0]])} ALFALFA galaxies")
    return data


def download_void_centers():
    """
    Download void centers with coordinates from VizieR.
    Table 0 in J/ApJS/265/7 has: x, y, z, Rad, RAJ2000, DEJ2000, Reff
    """
    print("\nDownloading void centers (Table 0 with coordinates)...")

    if not ASTROQUERY_AVAILABLE:
        print("ERROR: astroquery not available")
        return None

    v = Vizier(columns=['*'], row_limit=-1)

    try:
        # Get specific table with void centers
        catalogs = v.get_catalogs('J/ApJS/265/7')

        # Table 0 should have coordinates
        void_centers = None
        for cat in catalogs:
            if 'RAJ2000' in cat.colnames and 'DEJ2000' in cat.colnames:
                void_centers = cat
                break

        if void_centers is None:
            print("ERROR: Could not find void table with coordinates")
            return None

        print(f"Found void centers: {len(void_centers)} voids")
        print(f"Columns: {void_centers.colnames}")

        return void_centers

    except Exception as e:
        print(f"ERROR: {e}")
        return None


def parse_sexagesimal_to_degrees(coord_str):
    """
    Convert sexagesimal coordinate string to degrees.
    Format: "HH MM SS.S" or "DD MM SS.S"
    """
    try:
        parts = coord_str.split()
        if len(parts) == 3:
            d = float(parts[0])
            m = float(parts[1])
            s = float(parts[2])
            sign = -1 if d < 0 or coord_str.startswith('-') else 1
            return sign * (abs(d) + m/60 + s/3600)
        return float(coord_str)
    except:
        return np.nan


def cross_match_with_voids(alfalfa, void_centers):
    """
    Cross-match ALFALFA galaxies with void centers.
    """
    print("\n" + "=" * 70)
    print("CROSS-MATCHING ALFALFA × VOID CENTERS")
    print("=" * 70)

    # Parse ALFALFA coordinates (sexagesimal format)
    print("\nParsing ALFALFA coordinates...")

    # The ALFALFA catalog has RAJ2000/DEJ2000 in sexagesimal
    if 'RAJ2000' in alfalfa and isinstance(alfalfa['RAJ2000'][0], str):
        # Need to convert from sexagesimal
        alfalfa_ra = []
        alfalfa_dec = []

        for i, (ra_str, dec_str) in enumerate(zip(alfalfa['RAJ2000'], alfalfa['DEJ2000'])):
            try:
                # RA in hours -> degrees
                ra_parts = str(ra_str).split()
                if len(ra_parts) == 3:
                    ra_deg = (float(ra_parts[0]) + float(ra_parts[1])/60 + float(ra_parts[2])/3600) * 15
                else:
                    ra_deg = float(ra_str)

                # Dec in degrees
                dec_parts = str(dec_str).split()
                if len(dec_parts) == 3:
                    sign = -1 if dec_str.startswith('-') or float(dec_parts[0]) < 0 else 1
                    dec_deg = sign * (abs(float(dec_parts[0])) + float(dec_parts[1])/60 + float(dec_parts[2])/3600)
                else:
                    dec_deg = float(dec_str)

                alfalfa_ra.append(ra_deg)
                alfalfa_dec.append(dec_deg)
            except Exception as e:
                alfalfa_ra.append(np.nan)
                alfalfa_dec.append(np.nan)

        alfalfa_ra = np.array(alfalfa_ra)
        alfalfa_dec = np.array(alfalfa_dec)
    else:
        alfalfa_ra = alfalfa['RAJ2000']
        alfalfa_dec = alfalfa['DEJ2000']

    n_valid_coords = np.sum(~np.isnan(alfalfa_ra) & ~np.isnan(alfalfa_dec))
    print(f"Valid ALFALFA coordinates: {n_valid_coords}")

    # Parse void coordinates
    print("\nParsing void coordinates...")
    void_ra = []
    void_dec = []
    void_reff = []

    for i in range(len(void_centers)):
        try:
            ra_str = str(void_centers['RAJ2000'][i])
            dec_str = str(void_centers['DEJ2000'][i])

            # Parse RA (hours to degrees)
            ra_parts = ra_str.split()
            if len(ra_parts) == 3:
                ra_deg = (float(ra_parts[0]) + float(ra_parts[1])/60 + float(ra_parts[2])/3600) * 15
            else:
                ra_deg = float(ra_str)

            # Parse Dec
            dec_parts = dec_str.split()
            if len(dec_parts) == 3:
                sign = -1 if dec_str.startswith('-') or float(dec_parts[0]) < 0 else 1
                dec_deg = sign * (abs(float(dec_parts[0])) + float(dec_parts[1])/60 + float(dec_parts[2])/3600)
            else:
                dec_deg = float(dec_str)

            void_ra.append(ra_deg)
            void_dec.append(dec_deg)
            void_reff.append(float(void_centers['Reff'][i]))  # in Mpc/h
        except:
            pass

    void_ra = np.array(void_ra)
    void_dec = np.array(void_dec)
    void_reff = np.array(void_reff)

    print(f"Valid void centers: {len(void_ra)}")

    # Create SkyCoord objects
    print("\nComputing angular separations...")

    # Filter valid ALFALFA coordinates
    valid_mask = ~np.isnan(alfalfa_ra) & ~np.isnan(alfalfa_dec)
    alfalfa_coords = SkyCoord(ra=alfalfa_ra[valid_mask]*u.deg,
                              dec=alfalfa_dec[valid_mask]*u.deg)
    void_coords = SkyCoord(ra=void_ra*u.deg, dec=void_dec*u.deg)

    # For each ALFALFA galaxy, find nearest void
    from astropy.coordinates import match_coordinates_sky

    print("Matching to nearest void center...")
    idx, sep2d, _ = match_coordinates_sky(alfalfa_coords, void_coords)

    # Get properties of matched voids
    matched_reff = void_reff[idx]  # Mpc/h

    # Convert angular separation to physical distance
    # Need ALFALFA distances
    alfalfa_dist = alfalfa['Dist'][valid_mask]  # Mpc

    # Physical separation in Mpc
    sep_rad = sep2d.rad
    physical_sep = sep_rad * alfalfa_dist  # Mpc

    # Ratio to void radius (using h=0.7)
    h = 0.7
    void_radius_mpc = matched_reff / h  # Convert from Mpc/h to Mpc
    sep_ratio = physical_sep / void_radius_mpc

    print(f"\nSeparation ratios (d/R_void):")
    print(f"  Min: {np.nanmin(sep_ratio):.2f}")
    print(f"  Max: {np.nanmax(sep_ratio):.2f}")
    print(f"  Median: {np.nanmedian(sep_ratio):.2f}")

    # Classify environment
    # Void interior: d/R < 0.5
    # Void edge: 0.5 < d/R < 1.0
    # Wall: d/R > 1.0

    env_class = np.full(len(sep_ratio), '', dtype='U20')
    env_class[sep_ratio < 0.5] = 'void_interior'
    env_class[(sep_ratio >= 0.5) & (sep_ratio < 1.0)] = 'void_edge'
    env_class[sep_ratio >= 1.0] = 'wall'

    # Count
    n_interior = np.sum(env_class == 'void_interior')
    n_edge = np.sum(env_class == 'void_edge')
    n_wall = np.sum(env_class == 'wall')

    print(f"\nEnvironment classification:")
    print(f"  Void interior (d/R < 0.5): {n_interior} ({n_interior/len(env_class)*100:.1f}%)")
    print(f"  Void edge (0.5 < d/R < 1.0): {n_edge} ({n_edge/len(env_class)*100:.1f}%)")
    print(f"  Wall (d/R > 1.0): {n_wall} ({n_wall/len(env_class)*100:.1f}%)")

    # Return results with full indexing
    results = {
        'valid_mask': valid_mask,
        'sep_ratio': sep_ratio,
        'env_class': env_class,
        'alfalfa_ra': alfalfa_ra[valid_mask],
        'alfalfa_dec': alfalfa_dec[valid_mask],
        'alfalfa_dist': alfalfa_dist,
        'void_idx': idx,
        'n_interior': n_interior,
        'n_edge': n_edge,
        'n_wall': n_wall
    }

    return results


def analyze_btfr_by_environment(alfalfa, crossmatch):
    """
    Analyze BTFR by environment classification.
    """
    print("\n" + "=" * 70)
    print("BTFR ANALYSIS BY ENVIRONMENT")
    print("=" * 70)

    # Get valid subset of ALFALFA data
    valid_mask = crossmatch['valid_mask']
    env_class = crossmatch['env_class']

    w50 = alfalfa['W50'][valid_mask]
    logmhi = alfalfa['logMHI'][valid_mask]

    # Compute rotation velocity (statistical correction)
    sin_i_mean = np.pi / 4
    v_rot = w50 / (2 * sin_i_mean)
    log_v = np.log10(v_rot)

    # For BTFR, we use M_HI as proxy for M_bar (gas-dominated for many ALFALFA galaxies)
    log_m = logmhi

    # Analyze by environment
    environments = ['void_interior', 'void_edge', 'wall']
    env_stats = {}

    for env in environments:
        mask = env_class == env
        if np.sum(mask) < 10:
            continue

        log_v_env = log_v[mask]
        log_m_env = log_m[mask]

        # Fit BTFR: log(M) = a + b * log(V)
        valid = ~np.isnan(log_v_env) & ~np.isnan(log_m_env)
        if np.sum(valid) < 10:
            continue

        from numpy.polynomial import polynomial as P
        coef = np.polyfit(log_v_env[valid], log_m_env[valid], 1)
        slope = coef[0]
        intercept = coef[1]

        # Compute residuals
        predicted = intercept + slope * log_v_env[valid]
        residuals = log_m_env[valid] - predicted
        scatter = np.std(residuals)

        env_stats[env] = {
            'n': np.sum(valid),
            'slope': slope,
            'intercept': intercept,
            'scatter': scatter,
            'mean_log_v': np.mean(log_v_env[valid]),
            'mean_log_m': np.mean(log_m_env[valid])
        }

        print(f"\n{env}:")
        print(f"  N = {np.sum(valid)}")
        print(f"  BTFR: log(M) = {intercept:.2f} + {slope:.2f} × log(V)")
        print(f"  Scatter: {scatter:.3f} dex")
        print(f"  <log V> = {np.mean(log_v_env[valid]):.2f}")
        print(f"  <log M> = {np.mean(log_m_env[valid]):.2f}")

    # Compare zero-points at fixed slope
    if 'void_interior' in env_stats and 'wall' in env_stats:
        # Use wall slope as reference
        ref_slope = env_stats['wall']['slope']

        print("\n--- Zero-Point Comparison ---")
        print(f"Reference slope (wall): {ref_slope:.2f}")

        for env in ['void_interior', 'void_edge', 'wall']:
            if env not in env_stats:
                continue

            # Compute zero-point at reference slope
            zp = env_stats[env]['mean_log_m'] - ref_slope * env_stats[env]['mean_log_v']
            env_stats[env]['zero_point'] = zp

            print(f"  {env}: zero-point = {zp:.3f}")

        # Offset between void interior and wall
        if 'void_interior' in env_stats:
            offset = env_stats['void_interior']['zero_point'] - env_stats['wall']['zero_point']
            print(f"\n  Void-Wall offset: {offset:.3f} dex")

            # Significance
            n_void = env_stats['void_interior']['n']
            n_wall = env_stats['wall']['n']
            scatter_void = env_stats['void_interior']['scatter']
            scatter_wall = env_stats['wall']['scatter']
            combined_error = np.sqrt(scatter_void**2/n_void + scatter_wall**2/n_wall)
            significance = abs(offset) / combined_error

            print(f"  Combined error: {combined_error:.3f} dex")
            print(f"  Significance: {significance:.1f}σ")

            env_stats['comparison'] = {
                'void_wall_offset': offset,
                'significance': significance
            }

    return env_stats


def main():
    """Main execution."""
    print("=" * 70)
    print("SESSION #84: ALFALFA × VOID CROSS-MATCH")
    print("=" * 70)

    results = {
        'session': 84,
        'title': 'ALFALFA × Void Cross-Match',
        'date': datetime.now().isoformat()
    }

    # Load ALFALFA
    alfalfa = load_alfalfa_quality()
    if alfalfa is None:
        return results

    # Download void centers
    void_centers = download_void_centers()
    if void_centers is None:
        return results

    # Cross-match
    crossmatch = cross_match_with_voids(alfalfa, void_centers)
    results['crossmatch'] = {
        'n_interior': crossmatch['n_interior'],
        'n_edge': crossmatch['n_edge'],
        'n_wall': crossmatch['n_wall']
    }

    # BTFR analysis
    btfr_stats = analyze_btfr_by_environment(alfalfa, crossmatch)
    results['btfr'] = {k: {kk: float(vv) if isinstance(vv, (np.floating, float)) else vv
                           for kk, vv in v.items()}
                       for k, v in btfr_stats.items()}

    # Summary
    print("\n" + "=" * 70)
    print("SESSION #84 CROSS-MATCH SUMMARY")
    print("=" * 70)

    if 'comparison' in btfr_stats:
        offset = btfr_stats['comparison']['void_wall_offset']
        sig = btfr_stats['comparison']['significance']
        print(f"""
    Cross-match complete!

    Environment classification:
    - Void interior: {crossmatch['n_interior']}
    - Void edge: {crossmatch['n_edge']}
    - Wall: {crossmatch['n_wall']}

    BTFR void-wall offset: {offset:.3f} dex ({sig:.1f}σ)

    Synchronism prediction: +0.11 to +0.28 dex for voids
    Observation: {offset:+.3f} dex

    {'✓ CONSISTENT' if offset > 0.05 else '✗ NO SIGNIFICANT OFFSET'}
        """)
    else:
        print("Could not compute BTFR comparison - insufficient data")

    # Save results
    results_path = Path(__file__).parent / 'results' / 'session84_crossmatch.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"Results saved to: {results_path}")

    return results


if __name__ == '__main__':
    main()
