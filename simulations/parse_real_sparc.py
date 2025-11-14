#!/usr/bin/env python3
"""
Parse Real SPARC Mass Models Data
===================================

Converts SPARC Lelli et al. (2016) mass models table into per-galaxy
data files compatible with synchronism_sparc_validation.py

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-14
Session: #16
"""

import numpy as np
import os
from collections import defaultdict


def parse_sparc_mass_models(input_file: str, output_dir: str):
    """
    Parse SPARC MassModels_Lelli2016c.mrt and create per-galaxy files.

    Args:
        input_file: Path to SPARC mass models table
        output_dir: Directory to write individual galaxy files
    """

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Read the data file (skip header)
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Find where data starts (after "---" separator)
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('---'):
            data_start = i + 1
            break

    # Parse data into per-galaxy dictionaries
    galaxies = defaultdict(list)

    for line in lines[data_start:]:
        if not line.strip():
            continue

        # Parse columns (fixed-width format from byte-by-byte description)
        try:
            galaxy_id = line[0:11].strip()
            distance = float(line[12:18])
            radius = float(line[19:25])
            v_obs = float(line[26:32])
            e_vobs = float(line[33:38])
            v_gas = float(line[39:45])
            v_disk = float(line[46:52])
            v_bul = float(line[53:59])
            sb_disk = float(line[60:67])
            sb_bul = float(line[68:76])

            galaxies[galaxy_id].append({
                'distance': distance,
                'radius': radius,
                'v_obs': v_obs,
                'e_vobs': e_vobs,
                'v_gas': v_gas,
                'v_disk': v_disk,
                'v_bul': v_bul,
                'sb_disk': sb_disk,
                'sb_bul': sb_bul
            })

        except (ValueError, IndexError) as e:
            # Skip malformed lines
            continue

    # Write per-galaxy files
    n_galaxies = 0
    for galaxy_id, data_points in galaxies.items():
        # Convert to numpy arrays
        n_points = len(data_points)

        # Skip galaxies with < 3 points
        if n_points < 3:
            continue

        distance = data_points[0]['distance']  # Same for all points

        radius = np.array([p['radius'] for p in data_points])
        v_obs = np.array([p['v_obs'] for p in data_points])
        e_vobs = np.array([p['e_vobs'] for p in data_points])
        v_gas = np.array([p['v_gas'] for p in data_points])
        v_disk = np.array([p['v_disk'] for p in data_points])
        v_bul = np.array([p['v_bul'] for p in data_points])
        sb_disk = np.array([p['sb_disk'] for p in data_points])
        sb_bul = np.array([p['sb_bul'] for p in data_points])

        # Write file in same format as synthetic data
        filename = os.path.join(output_dir, f"{galaxy_id}.dat")
        with open(filename, 'w') as f:
            # Header
            f.write(f"# Galaxy: {galaxy_id}\n")
            f.write(f"# Distance: {distance:.2f} Mpc\n")
            f.write(f"# Source: SPARC (Lelli et al. 2016)\n")
            f.write(f"# Data: Real observational rotation curves\n")
            f.write("# Columns: Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n")
            f.write("# Units: kpc km/s km/s km/s km/s km/s L_☉/pc² L_☉/pc²\n")
            f.write("# Note: Vdisk, Vbul for M/L=1 at 3.6μm\n")

            # Data rows
            for i in range(n_points):
                f.write(f"{radius[i]:8.3f} {v_obs[i]:8.2f} {e_vobs[i]:7.2f} "
                       f"{v_gas[i]:8.2f} {v_disk[i]:8.2f} {v_bul[i]:8.2f} "
                       f"{sb_disk[i]:10.2f} {sb_bul[i]:10.2f}\n")

        n_galaxies += 1

    print(f"Parsed {n_galaxies} galaxies from SPARC mass models")
    print(f"Output directory: {output_dir}")

    return n_galaxies


def main():
    """Parse SPARC data."""
    input_file = "sparc_real_data/MassModels_Lelli2016c.mrt"
    output_dir = "sparc_real_data/galaxies"

    print("=" * 60)
    print("Parsing Real SPARC Mass Models Data")
    print("=" * 60)
    print()
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print()

    n_galaxies = parse_sparc_mass_models(input_file, output_dir)

    print()
    print(f"✓ Successfully parsed {n_galaxies} galaxies")
    print()
    print("Ready for Synchronism validation (Session #16)")
    print()


if __name__ == "__main__":
    main()
