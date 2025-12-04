#!/usr/bin/env python3
"""
Session #81 Track A: ALFALFA Data Analysis for Void Galaxy Test

Downloads and analyzes ALFALFA α.100 catalog with SDSS cross-match
to test the void galaxy BTFR prediction.

Data Sources:
- ALFALFA α.100: Haynes et al. (2018), ApJ 861, 49
- ALFALFA-SDSS: Durbala et al. (2020)
- Available at: https://egg.astro.cornell.edu/alfalfa/data/

Author: CBP Autonomous Synchronism Research
Date: December 3, 2025
Session: #81 - Void Galaxy Test Implementation
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime
import urllib.request
import os

# Data directory
DATA_DIR = Path(__file__).parent / 'alfalfa_data'


def print_section(title):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def download_alfalfa_catalog():
    """
    Download ALFALFA α.100 catalog if not present.
    """
    print_section("Step 1: ALFALFA Catalog Access")

    DATA_DIR.mkdir(exist_ok=True)

    # URLs for ALFALFA data
    catalog_url = "https://egg.astro.cornell.edu/alfalfa/data/a100.code12.table2.190808.csv"
    catalog_file = DATA_DIR / "a100_catalog.csv"

    print(f"\n    ALFALFA α.100 Catalog")
    print(f"    URL: {catalog_url}")
    print(f"    Local: {catalog_file}")

    if catalog_file.exists():
        print(f"    Status: Already downloaded ({catalog_file.stat().st_size / 1e6:.1f} MB)")
        return catalog_file
    else:
        print("    Status: Attempting download...")
        try:
            urllib.request.urlretrieve(catalog_url, catalog_file)
            print(f"    Success! Downloaded {catalog_file.stat().st_size / 1e6:.1f} MB")
            return catalog_file
        except Exception as e:
            print(f"    Download failed: {e}")
            print("\n    Manual download instructions:")
            print(f"    1. Visit: https://egg.astro.cornell.edu/alfalfa/data/")
            print(f"    2. Download: a100.code12.table2.190808.csv")
            print(f"    3. Save to: {catalog_file}")
            return None


def analyze_catalog_structure(catalog_file):
    """
    Analyze the structure of the ALFALFA catalog.
    """
    print_section("Step 2: Catalog Structure Analysis")

    if not catalog_file or not catalog_file.exists():
        print("    Catalog file not found. Creating simulated analysis...")
        return analyze_expected_structure()

    # Read first few lines to understand format
    with open(catalog_file, 'r') as f:
        lines = [f.readline() for _ in range(10)]

    print("\n    First 5 lines of catalog:")
    for i, line in enumerate(lines[:5]):
        print(f"    {i}: {line[:100]}...")

    return {'status': 'file_exists', 'path': str(catalog_file)}


def analyze_expected_structure():
    """
    Document expected ALFALFA catalog structure based on Haynes et al. 2018.
    """
    print("""
    Expected ALFALFA α.100 Catalog Columns (Haynes et al. 2018):
    ------------------------------------------------------------

    Key columns for void galaxy analysis:

    1. AGCNr       - Arecibo General Catalog number (ID)
    2. Name        - Object name
    3. RAdeg       - Right Ascension (degrees)
    4. DEdeg       - Declination (degrees)
    5. Vhelio      - Heliocentric velocity (km/s)
    6. W50         - HI line width at 50% peak (km/s)
    7. W20         - HI line width at 20% peak (km/s)
    8. Fc          - Integrated HI flux (Jy km/s)
    9. SNR         - Signal-to-noise ratio
    10. Dist       - Distance (Mpc)
    11. logMHI     - Log HI mass (M_sun)
    12. Code       - Detection code (1 = high quality)

    From ALFALFA-SDSS cross-match (Durbala et al. 2020):
    - Stellar mass (log M_star)
    - Environment metrics
    - Inclination estimates
    """)

    # Expected column info
    columns = {
        'AGCNr': 'Galaxy ID',
        'RAdeg': 'RA (deg)',
        'DEdeg': 'Dec (deg)',
        'Vhelio': 'Heliocentric velocity (km/s)',
        'W50': 'HI width at 50% (km/s)',
        'W20': 'HI width at 20% (km/s)',
        'Fc': 'HI flux (Jy km/s)',
        'SNR': 'Signal-to-noise',
        'Dist': 'Distance (Mpc)',
        'logMHI': 'Log HI mass (M_sun)'
    }

    return {'status': 'expected_structure', 'columns': columns}


def estimate_void_fraction():
    """
    Estimate what fraction of ALFALFA galaxies are in voids.
    """
    print_section("Step 3: Void Fraction Estimation")

    print("""
    Void Fraction Estimates:
    -----------------------

    From cosmic web statistics:
    - Voids occupy ~60% of cosmic volume
    - But galaxies are clustered, so void galaxy fraction is lower
    - Typical estimates: 10-20% of galaxies in voids

    ALFALFA specifics:
    - HI-selected (biased toward gas-rich, late-type)
    - These are more common in low-density environments
    - Expected void fraction: 15-25%

    From Kreckel et al. (2012) void galaxy studies:
    - Used SDSS + ALFALFA
    - Found ~7% of galaxies in "deep voids"
    - ~20% in "void regions" (δ < -0.5)

    For our analysis:
    - Conservative "void" cut: D_neighbor > 5 Mpc, δ < -0.5
    - Expected: ~10-15% of ALFALFA sample
    - With ~31,500 sources: ~3,000-5,000 void galaxies

    After quality cuts (S/N > 6, i > 45°, etc.):
    - Total: ~15,000 galaxies
    - Void: ~2,000 galaxies (expected)
    """)

    estimates = {
        'total_alfalfa': 31500,
        'after_quality_cuts': 15000,
        'void_fraction': 0.13,
        'expected_void_galaxies': 2000,
        'expected_cluster_galaxies': 3000,
        'expected_field_galaxies': 10000
    }

    return estimates


def btfr_analysis_plan():
    """
    Detailed plan for BTFR analysis.
    """
    print_section("Step 4: BTFR Analysis Plan")

    print("""
    BTFR Analysis Pipeline:
    ----------------------

    Step 1: Data Preparation
    ~~~~~~~~~~~~~~~~~~~~~~~~
    - Load ALFALFA α.100 catalog
    - Apply quality cuts:
      * Code = 1 (high quality detections)
      * SNR > 6
      * W50 > 0 (valid line width)
      * Dist > 0 (valid distance)

    Step 2: Compute Physical Properties
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - V_rot = W50 / (2 × sin(i))
      * Need inclination from SDSS (axis ratio b/a)
      * i = arccos(b/a) with disk thickness correction

    - M_HI = 2.356 × 10^5 × D² × F_HI  (in M_sun)
      * D in Mpc, F_HI in Jy km/s

    - M_bar = M_star + 1.4 × M_HI
      * 1.4 factor for helium + metals
      * M_star from SDSS photometry

    Step 3: Environment Classification
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Option A: SDSS-based neighbor density
    - Count L* (M_r < -20) neighbors within 5 Mpc
    - Void: N_neighbor = 0 or D_nearest > 5 Mpc
    - Cluster: N_neighbor > 10 or D_nearest < 2 Mpc

    Option B: Use existing void catalogs
    - Pan et al. (2012) SDSS void catalog
    - Cross-match ALFALFA positions

    Step 4: BTFR Fitting
    ~~~~~~~~~~~~~~~~~~~~
    - Fit: log(M_bar) = α × log(V_rot) + β
    - Separately for void / field / cluster

    - Use orthogonal distance regression (ODR)
    - Include measurement uncertainties

    Step 5: Statistical Comparison
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - Compare β (normalization) between environments
    - Predicted: β_void - β_cluster = -1.44 dex (at fixed V)
    - Or equivalently: Δlog(V) = +0.36 dex at fixed M

    - Compute significance:
      * σ = √(σ_void² + σ_cluster²)
      * z = Δβ / σ
    """)

    plan = {
        'step1': 'Data loading and quality cuts',
        'step2': 'Compute V_rot, M_HI, M_bar',
        'step3': 'Environment classification',
        'step4': 'BTFR fitting by environment',
        'step5': 'Statistical comparison'
    }

    return plan


def simulate_expected_results():
    """
    Simulate what we expect to see if Synchronism is correct.
    """
    print_section("Step 5: Expected Results Simulation")

    np.random.seed(42)

    # Standard BTFR: log(M_bar) = 4 × log(V) + 1.8
    alpha = 4.0
    beta_cluster = 1.8

    # Synchronism prediction: void offset
    delta_logV = 0.36  # At fixed M, void V is 0.36 dex higher
    # Equivalently: at fixed V, void M is 4 × 0.36 = 1.44 dex lower
    delta_logM = -1.44
    beta_void = beta_cluster + delta_logM

    # Generate mock data
    N_cluster = 3000
    N_void = 2000
    scatter = 0.10  # dex

    # Cluster galaxies
    logV_cluster = np.random.uniform(1.5, 2.5, N_cluster)  # 30-300 km/s
    logM_cluster = alpha * logV_cluster + beta_cluster + np.random.normal(0, scatter, N_cluster)

    # Void galaxies (Synchronism prediction)
    logV_void = np.random.uniform(1.5, 2.5, N_void)
    logM_void_sync = alpha * logV_void + beta_void + np.random.normal(0, scatter, N_void)

    # Void galaxies (MOND prediction - no offset)
    logM_void_mond = alpha * logV_void + beta_cluster + np.random.normal(0, scatter, N_void)

    print(f"""
    Mock BTFR Data:
    --------------
    Cluster: N = {N_cluster}
      log(V) range: [{logV_cluster.min():.2f}, {logV_cluster.max():.2f}]
      log(M) range: [{logM_cluster.min():.2f}, {logM_cluster.max():.2f}]
      Scatter: {scatter} dex

    Void (if Synchronism correct): N = {N_void}
      log(V) range: [{logV_void.min():.2f}, {logV_void.max():.2f}]
      log(M) range: [{logM_void_sync.min():.2f}, {logM_void_sync.max():.2f}]
      Offset: Δβ = {delta_logM:.2f} dex

    Void (if MOND correct):
      Same as cluster (no offset)
    """)

    # Compute detection significance
    sigma_mean_cluster = scatter / np.sqrt(N_cluster)
    sigma_mean_void = scatter / np.sqrt(N_void)
    sigma_combined = np.sqrt(sigma_mean_cluster**2 + sigma_mean_void**2)

    # Measure offset in simulated data
    mean_logM_cluster = np.mean(logM_cluster - alpha * logV_cluster)
    mean_logM_void_sync = np.mean(logM_void_sync - alpha * logV_void)
    measured_offset = mean_logM_void_sync - mean_logM_cluster

    significance = abs(measured_offset) / sigma_combined

    print(f"""
    Expected Detection:
    ------------------
    True offset: {delta_logM:.2f} dex
    Measured offset: {measured_offset:.3f} dex
    Combined σ: {sigma_combined:.4f} dex

    Detection significance: {significance:.1f}σ

    If Synchronism is correct:
      We expect ~{significance:.0f}σ detection of the offset

    If MOND is correct:
      We expect NO offset (consistent with 0 at ~{significance:.0f}σ)
    """)

    return {
        'N_cluster': N_cluster,
        'N_void': N_void,
        'predicted_offset': delta_logM,
        'expected_significance': float(significance),
        'scatter': scatter
    }


def implementation_status():
    """
    Current implementation status and next steps.
    """
    print_section("Implementation Status")

    print("""
    Current Status:
    --------------
    ✅ Session #80: Methodology developed
    ✅ Session #81: Data source identified (ALFALFA + SDSS)
    ⏳ Session #81: Download and initial analysis (in progress)

    Data Access:
    -----------
    Primary: https://egg.astro.cornell.edu/alfalfa/data/
    - a100.code12.table2.190808.csv (ALFALFA α.100)
    - durbala2020-table1.21-Sep-2020.fits.gz (SDSS cross-match)
    - durbala2020-table2.21-Sep-2020.fits.gz (Derived properties)

    Alternative: VizieR
    - Catalog: J/ApJ/861/49 (ALFALFA α.100)

    Next Steps:
    ----------
    1. Download full ALFALFA-SDSS cross-matched catalog
    2. Implement quality cuts
    3. Compute V_rot, M_bar for each galaxy
    4. Obtain/compute environment metrics:
       - Option A: Use existing void catalogs
       - Option B: Compute neighbor distances from SDSS
    5. Split by environment
    6. Fit BTFR to each subsample
    7. Measure offset and significance

    Estimated Time:
    --------------
    With data downloaded: 1-2 days to complete analysis
    Full pipeline: ~1 week

    Potential Issues:
    ----------------
    1. Inclination estimates may have systematic errors
    2. Environment classification depends on catalog used
    3. Stellar mass estimates have ~0.1-0.2 dex uncertainty
    4. Distance uncertainties from peculiar velocities

    Robustness Checks:
    -----------------
    - Vary void definition (D > 3, 5, 7 Mpc)
    - Vary inclination cut (i > 30°, 45°, 60°)
    - Compare different stellar mass estimators
    - Split by morphology, gas fraction
    """)

    return {
        'status': 'in_progress',
        'data_identified': True,
        'methodology_ready': True,
        'implementation_started': True,
        'next': 'Download full dataset and implement pipeline'
    }


def main():
    """Run ALFALFA data analysis setup."""
    print("=" * 70)
    print("SESSION #81 TRACK A: ALFALFA DATA ANALYSIS")
    print("Void Galaxy BTFR Test - Data Access")
    print("=" * 70)

    results = {}

    # Step 1: Download/locate catalog
    catalog_file = download_alfalfa_catalog()
    results['catalog'] = analyze_catalog_structure(catalog_file)

    # Step 2: Void fraction estimates
    results['void_fraction'] = estimate_void_fraction()

    # Step 3: Analysis plan
    results['plan'] = btfr_analysis_plan()

    # Step 4: Expected results
    results['expected'] = simulate_expected_results()

    # Step 5: Status
    results['status'] = implementation_status()

    # Save results
    output_path = Path(__file__).parent / 'results' / 'session81_alfalfa_analysis.json'
    output_path.parent.mkdir(exist_ok=True)

    output = {
        'session': 81,
        'track': 'A',
        'title': 'ALFALFA Data Analysis for Void Galaxy Test',
        'date': datetime.now().isoformat(),
        'data_source': 'https://egg.astro.cornell.edu/alfalfa/data/',
        'expected_significance': results['expected']['expected_significance'],
        'next_steps': [
            'Download full ALFALFA-SDSS dataset',
            'Implement quality cuts',
            'Compute V_rot and M_bar',
            'Classify environments',
            'Fit BTFR by environment',
            'Measure offset significance'
        ],
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\n\nResults saved to: {output_path}")

    print("\n" + "=" * 70)
    print("SESSION #81 TRACK A COMPLETE")
    print("=" * 70)

    return results


if __name__ == '__main__':
    main()
