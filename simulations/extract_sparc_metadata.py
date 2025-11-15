#!/usr/bin/env python3
"""
SPARC Galaxy Metadata Extractor - Session #29
===============================================

Extract galaxy-level metadata from individual SPARC .dat files for
correlation analysis.

Metadata needed:
- Distance (Mpc)
- Luminosity (derived from surface brightness)
- Mass (derived from rotation curve)
- Max velocity (from rotation curve)
- Galaxy type (from filename patterns)

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-14
Session: #29 - Fix metadata parsing
"""

import numpy as np
import os
import glob
import json
from typing import Dict, Optional


def extract_galaxy_metadata(filepath: str) -> Optional[Dict]:
    """
    Extract metadata from a single SPARC galaxy file.

    Returns dict with:
    - name: Galaxy identifier
    - distance_mpc: Distance in Mpc
    - v_max: Maximum observed velocity (km/s)
    - r_max: Radius at v_max (kpc)
    - r_last: Last measured radius (kpc)
    - n_points: Number of data points
    - luminosity_estimate: Rough luminosity from disk surface brightness
    - morphology: Inferred type (dwarf, spiral, elliptical, etc.)
    """
    try:
        # Extract name from filename
        name = os.path.basename(filepath).replace('.dat', '')

        # Read header for metadata
        distance_mpc = None
        with open(filepath, 'r') as f:
            for line in f:
                if 'Distance:' in line:
                    # Format: "# Distance: 23.70 Mpc"
                    parts = line.split(':')
                    if len(parts) >= 2:
                        distance_str = parts[1].strip().split()[0]
                        try:
                            distance_mpc = float(distance_str)
                        except ValueError:
                            pass
                    break

        if distance_mpc is None:
            print(f"Warning: No distance found for {name}")
            return None

        # Read rotation curve data
        data = np.genfromtxt(
            filepath,
            comments='#',
            names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul',
                   'SBdisk', 'SBbul']
        )

        if len(data) == 0:
            return None

        # Extract key properties
        radii = data['Rad']
        velocities = data['Vobs']
        surface_brightness = data['SBdisk']

        # Find maximum velocity
        v_max = np.max(velocities)
        v_max_idx = np.argmax(velocities)
        r_max = radii[v_max_idx]

        # Last measured radius
        r_last = radii[-1]

        # Number of data points
        n_points = len(data)

        # Rough luminosity estimate from integrated surface brightness
        # L ≈ ∑ SB × 2π r × Δr
        # This is crude but sufficient for correlation analysis
        luminosity_estimate = 0.0
        for i in range(len(radii) - 1):
            dr = radii[i+1] - radii[i]
            r_mid = (radii[i] + radii[i+1]) / 2.0
            sb_mid = (surface_brightness[i] + surface_brightness[i+1]) / 2.0
            luminosity_estimate += sb_mid * 2 * np.pi * r_mid * dr

        # Convert to log scale (avoid log(0))
        log_luminosity = np.log10(luminosity_estimate + 1.0)

        # Infer morphology from name patterns
        morphology = infer_morphology(name, v_max, r_last)

        return {
            'name': name,
            'distance_mpc': distance_mpc,
            'v_max': v_max,
            'r_max': r_max,
            'r_last': r_last,
            'n_points': n_points,
            'luminosity_estimate': luminosity_estimate,
            'log_luminosity': log_luminosity,
            'morphology': morphology,
            'mass_proxy': v_max**2 * r_last,  # Rough mass estimate: M ~ V^2 R / G
        }

    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None


def infer_morphology(name: str, v_max: float, r_last: float) -> str:
    """
    Infer galaxy morphology from name and basic properties.

    Patterns:
    - DDO*, UGC*: Often dwarf irregulars
    - NGC*: Spiral or elliptical
    - ESO*: Various types
    - F*, D*: Specific surveys

    Properties:
    - v_max < 50 km/s: Likely dwarf
    - v_max > 200 km/s: Likely massive spiral/elliptical
    - r_last < 5 kpc: Compact (dwarf or bulge-dominated)
    - r_last > 20 kpc: Extended disk
    """
    name_upper = name.upper()

    # Dwarf indicators
    if name_upper.startswith(('DDO', 'UGC')):
        if v_max < 80:
            return 'dwarf_irregular'
        elif v_max < 150:
            return 'small_spiral'
        else:
            return 'spiral'

    # NGC galaxies
    if name_upper.startswith('NGC'):
        if v_max > 200:
            return 'massive_spiral'
        elif v_max > 100:
            return 'spiral'
        else:
            return 'small_spiral'

    # ESO galaxies
    if name_upper.startswith('ESO'):
        if v_max < 100:
            return 'dwarf_irregular'
        else:
            return 'spiral'

    # Default classification by velocity
    if v_max < 50:
        return 'dwarf'
    elif v_max < 100:
        return 'intermediate'
    elif v_max < 200:
        return 'spiral'
    else:
        return 'massive'


def extract_all_metadata(data_dir: str = "sparc_real_data/galaxies") -> Dict[str, Dict]:
    """
    Extract metadata for all SPARC galaxies.

    Returns dict: galaxy_name -> metadata
    """
    metadata = {}

    # Find all galaxy files
    galaxy_files = glob.glob(os.path.join(data_dir, "*.dat"))

    print(f"Found {len(galaxy_files)} galaxy files")
    print()

    for filepath in sorted(galaxy_files):
        meta = extract_galaxy_metadata(filepath)
        if meta:
            metadata[meta['name']] = meta

    print(f"Successfully extracted metadata for {len(metadata)} galaxies")
    return metadata


def save_metadata(metadata: Dict[str, Dict], output_file: str = "sparc_galaxy_metadata.json"):
    """Save metadata to JSON file."""
    with open(output_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"Metadata saved to {output_file}")


def main():
    """Extract and save SPARC galaxy metadata."""
    print("=" * 70)
    print("SPARC Galaxy Metadata Extractor - Session #29")
    print("=" * 70)
    print()

    # Extract metadata
    metadata = extract_all_metadata()

    # Save to JSON
    save_metadata(metadata)

    # Print summary statistics
    print()
    print("Summary Statistics:")
    print(f"  Total galaxies: {len(metadata)}")

    if len(metadata) > 0:
        v_maxes = [m['v_max'] for m in metadata.values()]
        distances = [m['distance_mpc'] for m in metadata.values()]
        r_lasts = [m['r_last'] for m in metadata.values()]

        print(f"  V_max range: {min(v_maxes):.1f} - {max(v_maxes):.1f} km/s")
        print(f"  Distance range: {min(distances):.1f} - {max(distances):.1f} Mpc")
        print(f"  R_last range: {min(r_lasts):.1f} - {max(r_lasts):.1f} kpc")

        # Morphology distribution
        morphologies = [m['morphology'] for m in metadata.values()]
        from collections import Counter
        morph_counts = Counter(morphologies)
        print()
        print("Morphology Distribution:")
        for morph, count in sorted(morph_counts.items(), key=lambda x: -x[1]):
            print(f"  {morph}: {count}")

    print()
    print("=" * 70)
    print("Metadata extraction complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
