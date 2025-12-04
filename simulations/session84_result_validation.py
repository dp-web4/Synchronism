#!/usr/bin/env python3
"""
Session #84: Cross-Match Result Validation

The initial cross-match showed:
- Void-wall offset: -0.050 dex (OPPOSITE to Synchronism prediction)
- Significance: 6.7σ

But 60.8% classified as "void interior" seems too high.
Typical void filling factor is ~60% by volume, but galaxies preferentially
reside in denser regions (walls/filaments).

This script validates the methodology and interprets the result.

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
Session: #84
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime


def load_crossmatch_results():
    """Load cross-match results."""
    results_path = Path(__file__).parent / 'results' / 'session84_crossmatch.json'
    with open(results_path, 'r') as f:
        return json.load(f)


def analyze_methodology_issues():
    """
    Analyze potential issues with the cross-match methodology.
    """
    print("=" * 70)
    print("SESSION #84: METHODOLOGY VALIDATION")
    print("=" * 70)

    print("""
    ISSUE 1: Classification based on angular separation only
    ---------------------------------------------------------
    The current method finds the nearest void center and classifies
    based on d/R_void. This doesn't account for:
    - Galaxy redshift (may be behind/in front of void)
    - Void is 3D structure, we used 2D projection
    - A galaxy near void center in angle may not be IN the void

    ISSUE 2: 60.8% "void interior" is too high
    ------------------------------------------
    Observationally, ~5-10% of galaxies reside in voids.
    Our classification must be incorrect.

    ISSUE 3: Using HI mass as proxy for baryonic mass
    -------------------------------------------------
    ALFALFA measures M_HI, not M_bar = M_star + M_HI.
    For massive galaxies, M_star >> M_HI, so HI mass is poor proxy.

    ISSUE 4: Inclination effects
    ----------------------------
    W50 depends on inclination. Without individual corrections,
    this adds scatter but shouldn't create systematic offset.
    """)


def correct_classification_with_redshift():
    """
    Propose corrected classification using 3D positions.
    """
    print("\n" + "=" * 70)
    print("CORRECTED METHODOLOGY PROPOSAL")
    print("=" * 70)

    print("""
    PROPER METHOD:
    1. Convert galaxy RA, Dec, z to Cartesian (x, y, z) in Mpc
    2. For each void, convert center RA, Dec, z to Cartesian
    3. Compute 3D distance between galaxy and void center
    4. Compare to void radius (in Mpc, not angular)
    5. Classify based on 3D distance / R_void

    This requires:
    - Galaxy redshifts (ALFALFA has Vhelio → z)
    - Void center redshifts (available in void catalog)
    - Cosmological distance calculation

    NOTE: The ALFALFA 'Dist' column is distance in Mpc,
    but need to verify if void radii are comoving or proper.
    """)


def interpret_negative_offset():
    """
    Interpret the negative offset finding.
    """
    print("\n" + "=" * 70)
    print("INTERPRETING THE NEGATIVE OFFSET")
    print("=" * 70)

    print("""
    FINDING: Void-wall offset = -0.050 dex at 6.7σ

    This means: At fixed log(V), void galaxies have LOWER M_HI

    Possible interpretations:

    1. SYSTEMATIC EFFECT IN CLASSIFICATION
       - Our "void interior" sample is actually dominated by
         galaxies in moderately underdense regions, not true voids
       - The "wall" sample may include galaxies near void edges

    2. SELECTION BIAS
       - ALFALFA is HI-selected, preferentially finding gas-rich galaxies
       - Void galaxies may be more gas-rich (higher HI at fixed V)
       - But this would predict POSITIVE offset, not negative

    3. ACTUAL PHYSICAL EFFECT
       - Void galaxies have lower M_HI at fixed V
       - This would CONTRADICT Synchronism prediction
       - But also contradicts gas-rich void galaxy observations

    4. STAR FORMATION EFFICIENCY
       - If void galaxies form stars less efficiently
       - They retain more gas → higher M_HI
       - This should give positive offset

    MOST LIKELY: METHODOLOGY ISSUE (#1)
    The 60.8% void classification is clearly wrong.
    Need proper 3D void membership assignment.
    """)


def next_steps():
    """
    Define next steps for proper analysis.
    """
    print("\n" + "=" * 70)
    print("NEXT STEPS FOR PROPER VOID TEST")
    print("=" * 70)

    print("""
    SESSION #85 SHOULD:

    1. PROPER 3D CLASSIFICATION
       - Use ALFALFA Vhelio to compute distance
       - Convert void centers to 3D coordinates
       - Compute true 3D void membership

    2. USE VOID GALAXY CATALOG DIRECTLY
       - Douglass+ 2023 Table 4 has 776,500 galaxies with void membership
       - Cross-match ALFALFA with this table using NSAID
       - This requires NSA cross-match first

    3. COMPUTE M_BAR FROM M_STAR + M_HI
       - Get M_star from SDSS photometry
       - ALFALFA provides M_HI
       - M_bar = M_star + 1.4 × M_HI (helium correction)

    4. USE EXTREME VOIDS ONLY
       - The Synchronism prediction is strongest for δ < -0.8
       - Moderate voids (δ ~ -0.5) have small predicted offset
       - Need to select EXTREME void galaxies

    5. ALTERNATIVE: USE PUBLISHED VOID TF STUDIES
       - Dominguez-Gomez et al. (2019) studied this directly
       - Compare their result to Synchronism prediction
    """)


def summary():
    """
    Session summary.
    """
    print("\n" + "=" * 70)
    print("SESSION #84 VALIDATION SUMMARY")
    print("=" * 70)

    print("""
    PRELIMINARY RESULT: -0.050 dex void-wall offset (6.7σ)

    STATUS: UNRELIABLE due to methodology issues

    ISSUES IDENTIFIED:
    1. 2D angular classification (should be 3D)
    2. 60.8% "void interior" (should be ~5-10%)
    3. Using M_HI instead of M_bar
    4. Not distinguishing extreme vs moderate voids

    CONCLUSION:
    The current analysis CANNOT test Synchronism prediction.
    Need improved methodology for Session #85.

    THIS IS NOT A FALSIFICATION of Synchronism - it's a recognition
    that the test methodology needs refinement.

    KEY INSIGHT:
    Void galaxy BTFR analysis is more complex than anticipated.
    Simple angular cross-match is insufficient.
    """)

    return {
        'preliminary_offset': -0.050,
        'significance': 6.7,
        'status': 'UNRELIABLE',
        'issues': [
            '2D angular classification instead of 3D',
            'Unrealistic void fraction (60.8%)',
            'Using M_HI instead of M_bar',
            'Not selecting extreme voids'
        ],
        'conclusion': 'Methodology insufficient - need refinement',
        'next_priority': 'Proper 3D void membership or use published studies'
    }


def main():
    """Main execution."""
    analyze_methodology_issues()
    correct_classification_with_redshift()
    interpret_negative_offset()
    next_steps()
    result = summary()

    # Save
    output_path = Path(__file__).parent / 'results' / 'session84_validation.json'
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"\nValidation results saved to: {output_path}")

    return result


if __name__ == '__main__':
    main()
