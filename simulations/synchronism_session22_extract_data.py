#!/usr/bin/env python3
"""
Session #22: Extract Session #20 Data for Magnetic Screening Validation
========================================================================

Extracts individual galaxy (ρ_central, ρ_sat_fitted) pairs from Session #20
universality test for validation of Session #21 magnetic screening model.

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-17
Session: #22 - Empirical validation
"""

import numpy as np
import pickle
import os
import sys

# Import Session #20 infrastructure
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_rho_sat_universality import RhoSatUniversalityTester
from synchronism_real_sparc_validation import RealSPARCLoader


def extract_session20_data(save_file: str = "session20_rho_sat_data.pkl"):
    """
    Re-run Session #20 universality test and save individual galaxy results.

    This will take 10-30 minutes but provides real data for Session #22 validation.

    Args:
        save_file: Output filename for pickled results

    Returns:
        Dictionary with rho_centrals, rho_sats_fitted, galaxy_names, galaxy_types
    """
    print("="*80)
    print("SESSION #22: EXTRACTING SESSION #20 DATA")
    print("="*80)

    # Load SPARC galaxies
    print("\nLoading SPARC galaxies...")
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=None)  # All 175
    print(f"✓ Loaded {len(galaxies)} galaxies")

    # Fit with rational formula (Session #20 choice)
    print("\nFitting ρ_sat independently per galaxy (rational formula)...")
    print("This will take 10-30 minutes. Progress will be displayed.")

    tester = RhoSatUniversalityTester(coherence_model='rational')

    # Fit all galaxies
    results_list = []
    for i, galaxy in enumerate(galaxies):
        print(f"[{i+1}/{len(galaxies)}] Fitting {galaxy.name}...")
        result = tester.fit_galaxy_with_rho_sat(galaxy)
        results_list.append(result)

    # Extract successful fits
    successful = [r for r in results_list if r['success']]
    print(f"\n✓ Successfully fit {len(successful)}/{len(galaxies)} galaxies")

    # Extract data arrays
    rho_centrals = np.array([r['rho_central'] for r in successful if not np.isnan(r['rho_central'])])
    rho_sats_fitted = np.array([r['rho_sat'] for r in successful if not np.isnan(r['rho_central'])])
    galaxy_names = [r['name'] for r in successful if not np.isnan(r['rho_central'])]

    # Determine galaxy types
    galaxy_types = []
    for name in galaxy_names:
        if name.startswith('NGC'):
            galaxy_types.append('NGC')
        elif name.startswith('UGC'):
            galaxy_types.append('UGC')
        elif name.startswith('F'):
            galaxy_types.append('F')
        elif name.startswith('DDO'):
            galaxy_types.append('DDO')
        else:
            galaxy_types.append('Other')

    # Package data
    data = {
        'rho_centrals': rho_centrals,
        'rho_sats_fitted': rho_sats_fitted,
        'galaxy_names': galaxy_names,
        'galaxy_types': galaxy_types,
        'all_results': successful,
        'n_galaxies': len(galaxies),
        'n_success': len(successful)
    }

    # Save to disk
    save_path = os.path.join(os.path.dirname(__file__), save_file)
    with open(save_path, 'wb') as f:
        pickle.dump(data, f)

    print(f"\n✓ Data saved to: {save_path}")
    print(f"  {len(rho_centrals)} galaxies with valid (ρ_central, ρ_sat) pairs")

    # Print summary statistics
    print(f"\n{'='*80}")
    print("DATA SUMMARY")
    print(f"{'='*80}")

    print(f"\nρ_central:")
    print(f"  Min:    {np.min(rho_centrals):.2e} M_☉/pc³")
    print(f"  Median: {np.median(rho_centrals):.2e} M_☉/pc³")
    print(f"  Max:    {np.max(rho_centrals):.2e} M_☉/pc³")

    print(f"\nρ_sat_fitted:")
    print(f"  Min:    {np.min(rho_sats_fitted):.2e} M_☉/pc³")
    print(f"  Median: {np.median(rho_sats_fitted):.2e} M_☉/pc³")
    print(f"  Max:    {np.max(rho_sats_fitted):.2e} M_☉/pc³")

    # Correlation
    corr = np.corrcoef(np.log10(rho_centrals), np.log10(rho_sats_fitted))[0, 1]
    print(f"\nCorrelation(log ρ_c, log ρ_sat): r = {corr:.3f}")
    print(f"  Session #20 reported: r = -0.575")

    if abs(corr - (-0.575)) < 0.05:
        print(f"  ✅ Matches Session #20!")
    else:
        print(f"  ⚠️ Differs from Session #20 (expected due to resampling)")

    # Galaxy type breakdown
    print(f"\nGalaxy type distribution:")
    for gtype in ['NGC', 'UGC', 'F', 'DDO', 'Other']:
        count = galaxy_types.count(gtype)
        if count > 0:
            indices = [i for i, t in enumerate(galaxy_types) if t == gtype]
            rho_sat_median = np.median([rho_sats_fitted[i] for i in indices])
            print(f"  {gtype:5s}: {count:3d} galaxies, ρ_sat median = {rho_sat_median:.2e} M_☉/pc³")

    return data


def load_session20_data(save_file: str = "session20_rho_sat_data.pkl"):
    """
    Load previously extracted Session #20 data.

    Args:
        save_file: Filename of pickled results

    Returns:
        Dictionary with extracted data
    """
    save_path = os.path.join(os.path.dirname(__file__), save_file)

    if not os.path.exists(save_path):
        raise FileNotFoundError(
            f"Session #20 data not found: {save_path}\n"
            f"Run extract_session20_data() first to generate it."
        )

    with open(save_path, 'rb') as f:
        data = pickle.load(f)

    print(f"✓ Loaded Session #20 data from: {save_path}")
    print(f"  {len(data['rho_centrals'])} galaxies")

    return data


if __name__ == "__main__":
    # Check if data already exists
    save_file = "session20_rho_sat_data.pkl"
    save_path = os.path.join(os.path.dirname(__file__), save_file)

    if os.path.exists(save_path):
        print(f"Data file already exists: {save_path}")
        print("Loading existing data...")
        data = load_session20_data(save_file)
    else:
        print("Data file not found. Extracting from Session #20...")
        data = extract_session20_data(save_file)

    print("\n✓ Session #22 data extraction complete!")
    print("  Ready for magnetic screening validation")
