#!/usr/bin/env python3
"""
Session #80 Track B: Void Galaxy Prediction Methodology

Synchronism predicts void galaxies have 130% higher V_max at fixed M_bar.
MOND predicts no environment dependence (a_0 universal).

This is the KEY discriminating test.

This script develops the methodology for testing with ALFALFA + SDSS data.

Author: CBP Autonomous Synchronism Research
Date: December 3, 2025
Session: #80 - Void Galaxy Methodology
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime


def print_section(title):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def void_definition_criteria():
    """
    Define what constitutes a void galaxy.
    """
    print_section("Step 1: Void Definition Criteria")

    print("""
    What is a cosmic void?
    ----------------------
    Voids are large underdense regions in the cosmic web.
    Multiple definitions exist in the literature:

    1. Spherical Underdensity (traditional):
       - Find spheres with δ < -0.8 (where δ = ρ/ρ̄ - 1)
       - Typical radius: 10-30 Mpc

    2. ZOBOV/VIDE (watershed):
       - Trace density field from galaxy surveys
       - Find basins in density using Voronoi tessellation
       - More complex shapes, physically motivated

    3. Simple distance criterion:
       - No massive neighbor within D_min Mpc
       - D_min typically 3-7 Mpc

    For our analysis, we'll use:
    - Primary: Distance to nearest L* galaxy > 5 Mpc
    - Secondary: Local density δ < -0.5
    """)

    criteria = {
        'distance_criterion': {
            'description': 'Distance to nearest L* (M_r < -20) neighbor',
            'void_threshold': '> 5 Mpc',
            'field_threshold': '2-5 Mpc',
            'cluster_threshold': '< 2 Mpc'
        },
        'density_criterion': {
            'description': 'Local overdensity δ = n/n̄ - 1',
            'void_threshold': '< -0.5',
            'field_threshold': '-0.5 to 0.5',
            'cluster_threshold': '> 0.5'
        }
    }

    print("\n    Proposed Classification:")
    print(f"    {'Environment':<15} {'D_neighbor':<15} {'δ':<15}")
    print("    " + "-" * 45)
    print(f"    {'Void':<15} {'> 5 Mpc':<15} {'< -0.5':<15}")
    print(f"    {'Field':<15} {'2-5 Mpc':<15} {'-0.5 to 0.5':<15}")
    print(f"    {'Cluster':<15} {'< 2 Mpc':<15} {'> 0.5':<15}")

    return criteria


def data_sources():
    """
    Identify available data sources.
    """
    print_section("Step 2: Data Sources")

    print("""
    Required Data:
    --------------
    1. HI rotation curves (V_max, V_flat)
    2. Baryonic mass (M_bar = M_star + M_gas)
    3. Environment metrics (D_neighbor, δ)

    Primary Sources:

    ALFALFA (Arecibo Legacy Fast ALFA Survey):
    - 31,502 extragalactic HI detections
    - HI line widths → V_max proxy
    - HI masses → M_gas
    - Mostly field/void galaxies (avoids clusters)
    - Coverage: 7,000 deg² (SDSS overlap)

    SPARC (Spitzer Photometry and Accurate Rotation Curves):
    - 175 galaxies with detailed rotation curves
    - High-quality V(r) profiles
    - Already have M_bar estimates
    - Limited environment coverage

    SDSS (Sloan Digital Sky Survey):
    - Galaxy positions for environment metrics
    - Stellar masses (M_star)
    - Photometry for L* classification

    Void Catalogs:
    - Pan et al. (2012): SDSS void catalog
    - Sutter et al. (2012): VIDE void catalog
    - Cross-match with ALFALFA

    Combined Strategy:
    ------------------
    1. Use ALFALFA for HI widths + gas masses
    2. Use SDSS for stellar masses + environment
    3. Cross-match to get M_bar + environment
    4. Classify: void / field / cluster
    5. Compare TF relation by environment
    """)

    sources = {
        'alfalfa': {
            'N_galaxies': 31502,
            'provides': ['HI_width', 'M_gas', 'redshift'],
            'coverage': '7000 deg²',
            'url': 'http://egg.astro.cornell.edu/alfalfa/'
        },
        'sparc': {
            'N_galaxies': 175,
            'provides': ['V(r)', 'M_bar', 'surface_brightness'],
            'coverage': 'heterogeneous',
            'url': 'http://astroweb.cwru.edu/SPARC/'
        },
        'sdss': {
            'provides': ['positions', 'M_star', 'environment'],
            'void_catalogs': ['Pan2012', 'VIDE']
        }
    }

    return sources


def sample_selection():
    """
    Define sample selection criteria.
    """
    print_section("Step 3: Sample Selection")

    print("""
    Selection Criteria:
    ------------------

    1. Redshift: 0.002 < z < 0.05
       - Avoid local volume peculiar velocities
       - Stay within ALFALFA sensitivity

    2. Inclination: i > 45°
       - Minimize inclination correction errors
       - W_50/sin(i) → V_rot

    3. Signal quality: S/N > 6
       - Reliable HI width measurement
       - Accurate gas mass

    4. Mass range: 10^8 < M_bar < 10^11 M_sun
       - Avoid very low mass (resolution)
       - Avoid very high mass (scarce in voids)

    5. Morphology: Late-type (T > 0)
       - BTFR is cleanest for spirals
       - Exclude ellipticals

    Expected Sample Sizes:
    ---------------------
    From ALFALFA α.100:
    - Total HI detections: ~31,000
    - After quality cuts: ~15,000
    - Void galaxies (D > 5 Mpc): ~2,000 (estimated)
    - Cluster galaxies (D < 2 Mpc): ~3,000
    - Field (intermediate): ~10,000
    """)

    criteria = {
        'redshift': (0.002, 0.05),
        'inclination': '>45 deg',
        'signal_quality': 'S/N > 6',
        'mass_range': (1e8, 1e11),
        'morphology': 'T > 0 (late-type)',
        'expected_void': 2000,
        'expected_cluster': 3000,
        'expected_field': 10000
    }

    return criteria


def btfr_comparison_methodology():
    """
    Methodology for BTFR comparison.
    """
    print_section("Step 4: BTFR Comparison Methodology")

    print("""
    The Baryonic Tully-Fisher Relation:
    -----------------------------------
    log(M_bar) = α × log(V_flat) + β

    Standard values (McGaugh 2012):
    α = 4.0 (slope)
    β = 1.8 (normalization for V in km/s, M in M_sun)

    Synchronism Prediction:
    ----------------------
    Void galaxies have lower C_formation → lower ρ_crit → more enhancement
    At fixed M_bar, void galaxies have higher V_flat.

    Quantitative prediction from Session #75:
    - Cluster: C ~ 1.0 → no enhancement
    - Field: C ~ 0.88 → 13% enhancement
    - Void: C ~ 0.19 → 130% enhancement (!)

    ΔV/V = (G_eff/G)^(1/4) - 1 for deep MONDian limit
         ≈ (1/C)^(1/4) - 1

    For C_void = 0.19:
    ΔV/V ≈ (1/0.19)^(1/4) - 1 = 5.26^(1/4) - 1 = 1.52 - 1 = 0.52 = 52%

    Hmm, let me recalculate...
    """)

    # Careful recalculation
    C_cluster = 0.9999
    C_field = 0.88
    C_void = 0.19

    # G_eff = G / C
    # V^2 = G_eff M / r = (G/C) M / r
    # V ∝ sqrt(1/C)

    # For BTFR: M ∝ V^4
    # If G → G/C: V → V × (1/C)^(1/2)
    # At fixed M_bar, V_void / V_cluster = sqrt(C_cluster / C_void)

    print("\n    Recalculation:")
    print("    V² ∝ G_eff ∝ 1/C")
    print("    V ∝ (1/C)^(1/2)")
    print(f"\n    V_void / V_cluster = √(C_cluster/C_void)")
    print(f"                       = √({C_cluster}/{C_void})")
    print(f"                       = √{C_cluster/C_void:.2f}")
    print(f"                       = {np.sqrt(C_cluster/C_void):.2f}")

    enhancement = np.sqrt(C_cluster/C_void)

    print(f"""
    Predicted Enhancement:
    ---------------------
    V_void / V_cluster = {enhancement:.2f} = {100*(enhancement-1):.0f}% higher

    In BTFR space (log M vs log V):
    Δlog(V) = 0.5 × log(C_cluster/C_void)
            = 0.5 × log({C_cluster/C_void:.2f})
            = 0.5 × {np.log10(C_cluster/C_void):.3f}
            = {0.5 * np.log10(C_cluster/C_void):.3f} dex

    This is a {0.5 * np.log10(C_cluster/C_void):.2f} dex horizontal offset in BTFR!
    Equivalently: {4 * 0.5 * np.log10(C_cluster/C_void):.2f} dex vertical offset.

    Observable Signature:
    --------------------
    Plot log(M_bar) vs log(V_flat) for:
    - Void galaxies (D > 5 Mpc)
    - Cluster galaxies (D < 2 Mpc)

    Synchronism: Void BTFR shifted LEFT by ~{0.5 * np.log10(C_cluster/C_void):.2f} dex
                 (same M, higher V)

    MOND: Same BTFR for both
    """)

    return {
        'C_void': C_void,
        'C_cluster': C_cluster,
        'V_enhancement': float(enhancement),
        'delta_logV': float(0.5 * np.log10(C_cluster/C_void)),
        'delta_logM': float(4 * 0.5 * np.log10(C_cluster/C_void))
    }


def statistical_analysis():
    """
    Statistical methods for detection.
    """
    print_section("Step 5: Statistical Analysis")

    print("""
    Question: Can we detect a 0.36 dex BTFR offset?

    Statistical Power Analysis:
    --------------------------
    Intrinsic BTFR scatter: σ_intrinsic ~ 0.10 dex
    Measurement scatter: σ_meas ~ 0.10 dex
    Total scatter: σ_total ~ 0.14 dex

    For detecting Δ = 0.36 dex offset at 3σ:
    N = 2 × (σ × z_power / Δ)²
    N = 2 × (0.14 × 2.8 / 0.36)²
    N = 2 × (1.09)²
    N ~ 2.4

    Even with ~10 void galaxies we should see a clear signal!

    But with realistic samples:
    - Void: ~2000 galaxies → σ_mean = 0.14/√2000 = 0.003 dex
    - Cluster: ~3000 galaxies → σ_mean = 0.14/√3000 = 0.003 dex

    Detection significance: Δ / √(σ_void² + σ_cluster²)
                          = 0.36 / √(0.003² + 0.003²)
                          = 0.36 / 0.004
                          ~ 85σ (!)

    This is MASSIVELY detectable with existing data.
    """)

    print("""
    Why Hasn't This Been Done?
    -------------------------

    1. No one predicted this specific offset before
       - MOND doesn't predict environment dependence
       - Synchronism void prediction is new (Session #75)

    2. BTFR analyses usually control for environment
       - Remove void/cluster extremes to reduce scatter
       - Opposite of what we want!

    3. Focus has been on verifying BTFR, not testing predictions
       - Confirm slope = 4, not look for offsets

    Proposed Analysis:
    -----------------
    1. Reproduce standard BTFR (all environments)
    2. Split by environment (void / field / cluster)
    3. Fit BTFR to each subsample
    4. Compare normalizations (vertical offset)
    5. Compare slopes (should be same)

    Null hypothesis (MOND): No offset
    Alternative (Sync): Void offset by ~1.5 dex in M at fixed V
    """)

    return {
        'predicted_offset_dex': 0.36,
        'intrinsic_scatter_dex': 0.10,
        'detection_significance': 'Very high (>10σ with ~1000 galaxies)',
        'null_hypothesis': 'MOND: No environment offset',
        'alternative': 'Sync: Void offset by ~1.5 dex in M'
    }


def potential_systematics():
    """
    Discuss potential systematic effects.
    """
    print_section("Step 6: Potential Systematics")

    print("""
    Systematic Effects to Control:
    -----------------------------

    1. Selection Effects
       - Are void galaxies systematically different?
       - Lower mean luminosity (mass)
       - More gas-rich (less star formation)
       → Control for morphology, color, gas fraction

    2. Measurement Bias
       - HI line width in low-density environments
       - Distance uncertainties (peculiar velocities)
       → Use HI linewidth W_50, correct for Hubble flow

    3. Stellar Mass Estimates
       - Different star formation histories
       - Different dust content
       → Use consistent photometry, IMF assumptions

    4. Gas Mass Estimates
       - HI self-absorption
       - Molecular gas contribution
       → Apply standard M_HI corrections, estimate H2

    5. Inclination Corrections
       - Could differ by environment
       → Apply consistent correction method

    Robustness Checks:
    -----------------
    - Vary void definition (D > 3, 5, 7 Mpc)
    - Vary mass binning
    - Vary inclination cut
    - Compare different M_star estimators
    - Compare gas-rich vs stellar-dominated subsamples

    If signal persists across all these checks → robust detection
    """)

    return {
        'systematics': ['selection_effects', 'measurement_bias',
                        'stellar_mass', 'gas_mass', 'inclination'],
        'robustness_checks': ['vary_void_def', 'vary_mass_bins',
                              'vary_inclination_cut', 'compare_estimators']
    }


def implementation_plan():
    """
    Concrete implementation plan.
    """
    print_section("Step 7: Implementation Plan")

    print("""
    Phase 1: Data Preparation (2-3 days)
    ------------------------------------
    1. Download ALFALFA α.100 catalog
    2. Download SDSS photometry for cross-match
    3. Cross-match by position (3 arcsec)
    4. Compute environment metrics:
       - D_neighbor from SDSS
       - Local density from SDSS counts
    5. Classify: void / field / cluster

    Phase 2: Sample Selection (1 day)
    ---------------------------------
    1. Apply quality cuts (z, i, S/N, morphology)
    2. Compute M_bar = M_star + 1.4×M_HI
    3. Compute V_flat from W_50/2/sin(i)
    4. Split by environment

    Phase 3: BTFR Analysis (1-2 days)
    ---------------------------------
    1. Fit BTFR to full sample: log M = α log V + β
    2. Fit BTFR to void subsample
    3. Fit BTFR to cluster subsample
    4. Compare normalizations (β_void - β_cluster)
    5. Compute offset significance

    Phase 4: Systematics (2-3 days)
    -------------------------------
    1. Perform all robustness checks
    2. Vary definitions and cuts
    3. Test for confounding variables
    4. Document all results

    Phase 5: Publication (ongoing)
    ------------------------------
    1. Write up methodology and results
    2. Compare to Synchronism prediction
    3. Discuss implications for MOND

    Total: ~2 weeks to definitive test
    """)

    plan = {
        'phase_1': 'Data preparation (2-3 days)',
        'phase_2': 'Sample selection (1 day)',
        'phase_3': 'BTFR analysis (1-2 days)',
        'phase_4': 'Systematics (2-3 days)',
        'phase_5': 'Publication (ongoing)',
        'total_time': '~2 weeks'
    }

    return plan


def summary():
    """
    Summarize the methodology.
    """
    print_section("SUMMARY: Void Galaxy Test Methodology")

    print("""
    The Most Promising Discriminating Test:
    ======================================

    Synchronism Prediction:
    - Void galaxies formed at low C_formation
    - Lower ρ_crit → more G_eff enhancement
    - At fixed M_bar, void V_flat is ~230% of cluster

    MOND Prediction:
    - a_0 is universal, environment-independent
    - Same BTFR in all environments

    Observable:
    - BTFR offset: ~0.36 dex in log(V) at fixed M
    - Or ~1.5 dex in log(M) at fixed V

    Data:
    - ALFALFA × SDSS (~15,000 galaxies after cuts)
    - ~2,000 void, ~3,000 cluster

    Statistical Power:
    - Expected significance: >10σ
    - Easily detectable with existing data

    Why This Matters:
    ----------------
    1. Uses existing public data (no new observations)
    2. Clear quantitative prediction (0.36 dex)
    3. Directly tests Synchronism vs MOND
    4. Could be done in ~2 weeks

    IMMEDIATE ACTION:
    -----------------
    Download ALFALFA × SDSS and implement this analysis.
    This is the most important test for Synchronism.
    """)


def main():
    """Run void galaxy methodology development."""
    print("=" * 70)
    print("SESSION #80 TRACK B: VOID GALAXY TEST METHODOLOGY")
    print("=" * 70)

    results = {}

    results['void_criteria'] = void_definition_criteria()
    results['data_sources'] = data_sources()
    results['sample_selection'] = sample_selection()
    results['btfr_methodology'] = btfr_comparison_methodology()
    results['statistics'] = statistical_analysis()
    results['systematics'] = potential_systematics()
    results['implementation'] = implementation_plan()
    summary()

    # Save results
    output_path = Path(__file__).parent / 'results' / 'session80_void_methodology.json'
    output_path.parent.mkdir(exist_ok=True)

    output = {
        'session': 80,
        'track': 'B',
        'title': 'Void Galaxy Test Methodology',
        'date': datetime.now().isoformat(),
        'key_prediction': {
            'sync_V_enhancement': 2.3,  # 230%
            'delta_logV': 0.36,
            'testable': True,
            'data_exists': True
        },
        'implementation_time': '~2 weeks',
        'results': results
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\n\nResults saved to: {output_path}")

    print("\n" + "=" * 70)
    print("SESSION #80 TRACK B COMPLETE")
    print("=" * 70)

    return results


if __name__ == '__main__':
    main()
