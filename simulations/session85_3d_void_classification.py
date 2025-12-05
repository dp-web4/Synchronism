#!/usr/bin/env python3
"""
Session #85: Proper 3D Void Membership Classification

Session #84 identified methodology issues with 2D angular classification.
This script implements proper 3D void membership using:
1. Galaxy Cartesian coordinates from RA, Dec, distance
2. Void center Cartesian coordinates from RA, Dec, redshift
3. 3D Euclidean distance to void center
4. Classification based on true 3D distance / void radius

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
Session: #85
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime

try:
    from astropy.cosmology import Planck18
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False


def load_alfalfa_quality():
    """Load quality-cut ALFALFA catalog."""
    alfalfa_path = Path(__file__).parent / 'alfalfa_data' / 'alfalfa_quality.csv'

    if not alfalfa_path.exists():
        print(f"ERROR: ALFALFA file not found: {alfalfa_path}")
        return None

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


def load_void_centers():
    """Load void centers from Session #84 download."""
    void_path = Path(__file__).parent / 'void_data' / 'void_galaxies.csv'

    # Actually we need the void centers, not the galaxies
    # Let me load from VizieR again
    try:
        from astroquery.vizier import Vizier

        v = Vizier(columns=['*'], row_limit=-1)
        catalogs = v.get_catalogs('J/ApJS/265/7')

        # Find table with void centers (has RA, Dec, Reff)
        for cat in catalogs:
            if 'RAJ2000' in cat.colnames and 'Reff' in cat.colnames:
                print(f"Found void centers: {len(cat)} voids")
                return cat

        print("ERROR: Could not find void center table")
        return None

    except Exception as e:
        print(f"ERROR loading void centers: {e}")
        return None


def parse_sexagesimal(coord_str, is_ra=True):
    """
    Parse sexagesimal coordinate string to degrees.
    RA: HH MM SS.S -> degrees (multiply by 15)
    Dec: DD MM SS.S -> degrees
    """
    try:
        parts = str(coord_str).split()
        if len(parts) == 3:
            d = float(parts[0])
            m = float(parts[1])
            s = float(parts[2])

            if is_ra:
                # RA: hours to degrees
                return (abs(d) + m/60 + s/3600) * 15
            else:
                # Dec: degrees
                sign = -1 if d < 0 or str(coord_str).startswith('-') else 1
                return sign * (abs(d) + m/60 + s/3600)
        else:
            return float(coord_str)
    except:
        return np.nan


def compute_cartesian_coordinates(ra_deg, dec_deg, dist_mpc):
    """
    Convert RA, Dec, distance to Cartesian (x, y, z) in Mpc.

    Uses standard astronomical convention:
    x = D * cos(Dec) * cos(RA)
    y = D * cos(Dec) * sin(RA)
    z = D * sin(Dec)
    """
    ra_rad = np.radians(ra_deg)
    dec_rad = np.radians(dec_deg)

    x = dist_mpc * np.cos(dec_rad) * np.cos(ra_rad)
    y = dist_mpc * np.cos(dec_rad) * np.sin(ra_rad)
    z = dist_mpc * np.sin(dec_rad)

    return x, y, z


def classify_void_membership_3d(alfalfa, void_centers):
    """
    Classify ALFALFA galaxies by 3D void membership.
    """
    print("\n" + "=" * 70)
    print("3D VOID MEMBERSHIP CLASSIFICATION")
    print("=" * 70)

    # Parse ALFALFA coordinates
    print("\nParsing ALFALFA coordinates...")

    alfalfa_ra = []
    alfalfa_dec = []

    for ra_str, dec_str in zip(alfalfa['RAJ2000'], alfalfa['DEJ2000']):
        alfalfa_ra.append(parse_sexagesimal(ra_str, is_ra=True))
        alfalfa_dec.append(parse_sexagesimal(dec_str, is_ra=False))

    alfalfa_ra = np.array(alfalfa_ra)
    alfalfa_dec = np.array(alfalfa_dec)
    alfalfa_dist = np.array(alfalfa['Dist'])  # Already in Mpc

    # Compute ALFALFA 3D coordinates
    print("Computing ALFALFA 3D coordinates...")
    alfalfa_x, alfalfa_y, alfalfa_z = compute_cartesian_coordinates(
        alfalfa_ra, alfalfa_dec, alfalfa_dist
    )

    # Valid galaxies mask
    valid_mask = ~np.isnan(alfalfa_x) & ~np.isnan(alfalfa_y) & ~np.isnan(alfalfa_z)
    n_valid = np.sum(valid_mask)
    print(f"Valid ALFALFA galaxies: {n_valid}")

    # Parse void coordinates
    print("\nParsing void coordinates...")

    void_ra = []
    void_dec = []
    void_dist = []  # Need to compute from z or use x,y,z directly
    void_reff = []

    # Check if void catalog has Cartesian coords or just RA/Dec
    if 'x' in void_centers.colnames:
        # Direct Cartesian coordinates available
        void_x = np.array(void_centers['x'])  # Mpc/h
        void_y = np.array(void_centers['y'])
        void_z_cart = np.array(void_centers['z'])
        void_reff = np.array(void_centers['Reff'])  # Mpc/h

        # Convert from Mpc/h to Mpc (h ≈ 0.7)
        h = 0.7
        void_x = void_x / h
        void_y = void_y / h
        void_z_cart = void_z_cart / h
        void_reff = void_reff / h

        print(f"Using direct Cartesian coordinates for {len(void_x)} voids")
    else:
        # Need to convert from RA/Dec
        for i in range(len(void_centers)):
            ra = parse_sexagesimal(void_centers['RAJ2000'][i], is_ra=True)
            dec = parse_sexagesimal(void_centers['DEJ2000'][i], is_ra=False)
            # Need redshift for distance
            # This is more complex - skip for now
            pass

    # Compute 3D distances to all void centers
    print("\nComputing 3D distances to void centers...")

    # For each galaxy, find minimum distance to any void center
    min_dist_ratio = np.full(len(alfalfa_x), np.inf)
    nearest_void = np.full(len(alfalfa_x), -1, dtype=int)

    for i in range(len(void_x)):
        if i % 500 == 0:
            print(f"  Processing void {i}/{len(void_x)}...")

        # 3D distance from all ALFALFA galaxies to this void center
        dx = alfalfa_x - void_x[i]
        dy = alfalfa_y - void_y[i]
        dz = alfalfa_z - void_z_cart[i]

        dist_3d = np.sqrt(dx**2 + dy**2 + dz**2)

        # Distance ratio (d / R_void)
        dist_ratio = dist_3d / void_reff[i]

        # Update minimum
        closer = dist_ratio < min_dist_ratio
        min_dist_ratio[closer] = dist_ratio[closer]
        nearest_void[closer] = i

    # Classify based on 3D distance ratio
    print("\nClassifying by 3D void membership...")

    env_class = np.full(len(min_dist_ratio), 'unknown', dtype='U20')
    env_class[min_dist_ratio < 0.5] = 'void_core'
    env_class[(min_dist_ratio >= 0.5) & (min_dist_ratio < 1.0)] = 'void_interior'
    env_class[(min_dist_ratio >= 1.0) & (min_dist_ratio < 2.0)] = 'void_edge'
    env_class[min_dist_ratio >= 2.0] = 'field'

    # Apply valid mask
    env_class[~valid_mask] = 'invalid'

    # Count
    n_core = np.sum(env_class == 'void_core')
    n_interior = np.sum(env_class == 'void_interior')
    n_edge = np.sum(env_class == 'void_edge')
    n_field = np.sum(env_class == 'field')
    n_invalid = np.sum(env_class == 'invalid')

    total_valid = n_core + n_interior + n_edge + n_field

    print(f"\n3D Environment classification:")
    print(f"  Void core (d/R < 0.5): {n_core} ({n_core/total_valid*100:.1f}%)")
    print(f"  Void interior (0.5 < d/R < 1.0): {n_interior} ({n_interior/total_valid*100:.1f}%)")
    print(f"  Void edge (1.0 < d/R < 2.0): {n_edge} ({n_edge/total_valid*100:.1f}%)")
    print(f"  Field (d/R > 2.0): {n_field} ({n_field/total_valid*100:.1f}%)")
    print(f"  Invalid coordinates: {n_invalid}")

    return {
        'env_class': env_class,
        'min_dist_ratio': min_dist_ratio,
        'nearest_void': nearest_void,
        'valid_mask': valid_mask,
        'n_core': n_core,
        'n_interior': n_interior,
        'n_edge': n_edge,
        'n_field': n_field
    }


def analyze_btfr_3d(alfalfa, classification):
    """
    Analyze BTFR by 3D environment classification.
    """
    print("\n" + "=" * 70)
    print("BTFR ANALYSIS BY 3D ENVIRONMENT")
    print("=" * 70)

    env_class = classification['env_class']
    valid_mask = classification['valid_mask']

    # Get ALFALFA data
    w50 = np.array(alfalfa['W50'])
    logmhi = np.array(alfalfa['logMHI'])

    # Compute rotation velocity
    sin_i_mean = np.pi / 4
    v_rot = w50 / (2 * sin_i_mean)
    log_v = np.log10(v_rot)
    log_m = logmhi

    # Combined void sample (core + interior)
    void_mask = (env_class == 'void_core') | (env_class == 'void_interior')
    field_mask = env_class == 'field'

    # Statistics by environment
    environments = {
        'void_core': env_class == 'void_core',
        'void_interior': env_class == 'void_interior',
        'void_all': void_mask,
        'void_edge': env_class == 'void_edge',
        'field': field_mask
    }

    env_stats = {}

    for env_name, mask in environments.items():
        valid = mask & ~np.isnan(log_v) & ~np.isnan(log_m)
        n = np.sum(valid)

        if n < 10:
            print(f"\n{env_name}: N = {n} (insufficient)")
            continue

        log_v_env = log_v[valid]
        log_m_env = log_m[valid]

        # Fit BTFR
        coef = np.polyfit(log_v_env, log_m_env, 1)
        slope = coef[0]
        intercept = coef[1]

        # Residuals
        predicted = intercept + slope * log_v_env
        residuals = log_m_env - predicted
        scatter = np.std(residuals)

        env_stats[env_name] = {
            'n': n,
            'slope': slope,
            'intercept': intercept,
            'scatter': scatter,
            'mean_log_v': np.mean(log_v_env),
            'mean_log_m': np.mean(log_m_env)
        }

        print(f"\n{env_name}:")
        print(f"  N = {n}")
        print(f"  BTFR: log(M_HI) = {intercept:.2f} + {slope:.2f} × log(V)")
        print(f"  Scatter: {scatter:.3f} dex")

    # Compute void-field offset
    if 'void_all' in env_stats and 'field' in env_stats:
        print("\n--- VOID-FIELD COMPARISON ---")

        # Use field slope as reference
        ref_slope = env_stats['field']['slope']
        print(f"Reference slope (field): {ref_slope:.2f}")

        for env_name in ['void_core', 'void_interior', 'void_all', 'field']:
            if env_name not in env_stats:
                continue
            zp = env_stats[env_name]['mean_log_m'] - ref_slope * env_stats[env_name]['mean_log_v']
            env_stats[env_name]['zero_point'] = zp
            print(f"  {env_name}: zero-point = {zp:.3f}")

        # Offset calculation
        if 'void_all' in env_stats:
            offset = env_stats['void_all']['zero_point'] - env_stats['field']['zero_point']

            # Significance
            n_void = env_stats['void_all']['n']
            n_field = env_stats['field']['n']
            scatter_void = env_stats['void_all']['scatter']
            scatter_field = env_stats['field']['scatter']

            combined_error = np.sqrt(scatter_void**2/n_void + scatter_field**2/n_field)
            significance = abs(offset) / combined_error if combined_error > 0 else 0

            print(f"\n  Void-Field offset: {offset:+.3f} dex")
            print(f"  Combined error: {combined_error:.3f} dex")
            print(f"  Significance: {significance:.1f}σ")

            env_stats['comparison'] = {
                'offset': offset,
                'error': combined_error,
                'significance': significance
            }

    return env_stats


def main():
    """Main execution."""
    print("=" * 70)
    print("SESSION #85: 3D VOID MEMBERSHIP CLASSIFICATION")
    print("=" * 70)

    results = {
        'session': 85,
        'title': '3D Void Membership Classification',
        'date': datetime.now().isoformat()
    }

    # Load data
    alfalfa = load_alfalfa_quality()
    if alfalfa is None:
        return results

    void_centers = load_void_centers()
    if void_centers is None:
        return results

    # 3D classification
    classification = classify_void_membership_3d(alfalfa, void_centers)
    results['classification'] = {
        'n_core': classification['n_core'],
        'n_interior': classification['n_interior'],
        'n_edge': classification['n_edge'],
        'n_field': classification['n_field']
    }

    # BTFR analysis
    btfr_stats = analyze_btfr_3d(alfalfa, classification)

    # Convert to serializable
    results['btfr'] = {}
    for k, v in btfr_stats.items():
        if isinstance(v, dict):
            results['btfr'][k] = {kk: float(vv) if isinstance(vv, (np.floating, float)) else vv
                                  for kk, vv in v.items()}

    # Summary
    print("\n" + "=" * 70)
    print("SESSION #85 3D CLASSIFICATION SUMMARY")
    print("=" * 70)

    total = classification['n_core'] + classification['n_interior'] + classification['n_edge'] + classification['n_field']
    void_frac = (classification['n_core'] + classification['n_interior']) / total * 100

    print(f"""
    3D Classification Results:
    - Void core: {classification['n_core']}
    - Void interior: {classification['n_interior']}
    - Void edge: {classification['n_edge']}
    - Field: {classification['n_field']}

    Void fraction (core + interior): {void_frac:.1f}%
    (Expected: ~5-15% of galaxies in voids)

    BTFR Results:
    """)

    if 'comparison' in btfr_stats:
        offset = btfr_stats['comparison']['offset']
        sig = btfr_stats['comparison']['significance']
        print(f"    Void-Field offset: {offset:+.3f} dex ({sig:.1f}σ)")
        print(f"\n    Synchronism prediction: +0.11 to +0.28 dex")
        print(f"    Observation: {offset:+.3f} dex")

        if offset > 0.05:
            print("\n    ✓ POSITIVE OFFSET - Consistent direction with Synchronism")
        elif offset < -0.05:
            print("\n    ✗ NEGATIVE OFFSET - Opposite to Synchronism prediction")
        else:
            print("\n    ~ NO SIGNIFICANT OFFSET detected")

    # Save
    results_path = Path(__file__).parent / 'results' / 'session85_3d_classification.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults saved to: {results_path}")

    return results


if __name__ == '__main__':
    main()
