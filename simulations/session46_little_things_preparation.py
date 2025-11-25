#!/usr/bin/env python3
"""
Session #46 Track C: LITTLE THINGS External Validation Preparation

Nova's Session #45 recommendation:
"External validation priority: Focus on LITTLE THINGS (dwarfs)"

This session prepares the data structure and validation framework for testing
Synchronism's dark matter model against the LITTLE THINGS dwarf galaxy survey.

LITTLE THINGS Survey:
- "Local Irregulars That Trace Luminosity Extremes, The H I Nearby Galaxy Survey"
- High-resolution VLA HI survey (~6" angular, <2.6 km/s velocity resolution)
- 26 dwarf galaxies within 11 Mpc
- Rotation curves and mass models published in Oh et al. (2015) AJ 149, 180
- arXiv:1502.01281

Key prediction:
- Synchronism performs best on dwarf galaxies (v_max < 100 km/s)
- Session #44: 67% success for v_max < 100, 81.8% for v_max < 50
- LITTLE THINGS dwarfs should show 60-80% success rate

Author: CBP Autonomous Synchronism Research
Date: 2025-11-25
Session: #46 - LITTLE THINGS Preparation
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime


def document_little_things_survey():
    """
    Document the LITTLE THINGS survey properties.
    """

    print("\n" + "="*80)
    print("LITTLE THINGS SURVEY DOCUMENTATION")
    print("="*80)

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                        LITTLE THINGS SURVEY OVERVIEW                        │
└─────────────────────────────────────────────────────────────────────────────┘

SURVEY NAME:
══════════════════════════════════════════════════════════════════════════════

    LITTLE THINGS = Local Irregulars That Trace Luminosity Extremes,
                    The H I Nearby Galaxy Survey


KEY REFERENCE:
══════════════════════════════════════════════════════════════════════════════

    Oh et al. (2015)
    "High-resolution mass models of dwarf galaxies from LITTLE THINGS"
    Astronomical Journal, 149, 180
    arXiv:1502.01281


SURVEY PROPERTIES:
══════════════════════════════════════════════════════════════════════════════

    Telescope:      Very Large Array (VLA)
    Wavelength:     HI 21 cm
    Resolution:     ~6" angular, <2.6 km/s velocity
    Volume:         Local volume within 11 Mpc
    Sample:         26 dwarf irregular galaxies
    Data Products:  - Rotation curves
                    - Mass models (baryonic + DM decomposition)
                    - Spitzer 3.6μm photometry
                    - Optical U, B, V photometry


KEY SCIENTIFIC FINDINGS:
══════════════════════════════════════════════════════════════════════════════

    1. Most galaxies show LINEAR inner rotation curves
    2. Mean DM density slope: α = -0.32 ± 0.24
    3. Significantly shallower than ΛCDM cusp (α = -1)
    4. Consistent with "cored" profiles (baryonic feedback effects)


WHY LITTLE THINGS FOR SYNCHRONISM VALIDATION:
══════════════════════════════════════════════════════════════════════════════

    1. DWARF GALAXIES: Synchronism excels at v_max < 100 km/s
       - SPARC results: 67% success for v_max < 100
       - SPARC results: 81.8% success for v_max < 50

    2. HIGH RESOLUTION: Accurate rotation curves
       - Better constraint on inner density profiles
       - Less systematic uncertainty

    3. INDEPENDENT SAMPLE: Not from SPARC
       - True external validation
       - Tests generalization of the model

    4. WELL-CHARACTERIZED: Published mass models
       - Can directly compare DM predictions
       - Consistent methodology

""")


def define_validation_framework():
    """
    Define the validation framework for LITTLE THINGS.
    """

    print("\n" + "="*80)
    print("VALIDATION FRAMEWORK FOR LITTLE THINGS")
    print("="*80)

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                    SYNCHRONISM VALIDATION METHODOLOGY                       │
└─────────────────────────────────────────────────────────────────────────────┘

SYNCHRONISM MODEL (0-parameter):
══════════════════════════════════════════════════════════════════════════════

    Step 1: Compute critical density from v_max
            ρ_crit = 0.25 × v_max^1.62   [M_☉/pc³]

    Step 2: Compute coherence profile
            C(r) = tanh(2.0 × log(ρ_vis(r)/ρ_crit + 1))

    Step 3: Predict dark matter density
            ρ_DM(r) = α × (1 - C(r)) × ρ_vis(r)^0.30

    Step 4: Compute rotation curve
            v²(r) = v²_vis(r) + v²_DM(r)

    Step 5: Compare to observed rotation curve


DATA REQUIREMENTS:
══════════════════════════════════════════════════════════════════════════════

    For each LITTLE THINGS galaxy, we need:

    ┌───────────────────────────┬───────────────────────────────────────┐
    │ Parameter                 │ Source                                │
    ├───────────────────────────┼───────────────────────────────────────┤
    │ v_max (km/s)              │ Maximum rotation velocity             │
    │ R (kpc) array             │ Radius sampling                       │
    │ v_obs(R) (km/s)           │ Observed rotation curve               │
    │ σ_v(R) (km/s)             │ Velocity uncertainties                │
    │ Σ_* (M_☉/pc²)             │ Stellar surface density (3.6μm)       │
    │ Σ_gas (M_☉/pc²)           │ Gas surface density (HI + He)         │
    │ Distance (Mpc)            │ Galaxy distance                       │
    │ Inclination (°)           │ Disk inclination                      │
    └───────────────────────────┴───────────────────────────────────────┘


SUCCESS CRITERION:
══════════════════════════════════════════════════════════════════════════════

    Same as SPARC validation: χ² < 5.0

    χ² = Σᵢ [(v_obs,i - v_pred,i) / σ_i]² / N_points

    Success: χ² < 5.0 (acceptable fit)
    Failure: χ² ≥ 5.0 (poor fit)


EXPECTED PERFORMANCE:
══════════════════════════════════════════════════════════════════════════════

    Based on SPARC dwarf results:

    If v_max distribution of LITTLE THINGS is:
    - Mostly v_max < 50 km/s:  Expected 70-80% success
    - Mostly v_max < 100 km/s: Expected 60-70% success
    - Mix up to 100 km/s:      Expected 55-65% success

    PREDICTION: 60-70% success rate on LITTLE THINGS

""")


def list_known_little_things_galaxies():
    """
    List known LITTLE THINGS galaxies from literature.
    """

    print("\n" + "="*80)
    print("LITTLE THINGS GALAXY SAMPLE")
    print("="*80)

    # Known LITTLE THINGS galaxies from Oh et al. (2015)
    # This is the sample from the paper - 26 dwarf irregulars
    little_things_sample = [
        {'name': 'CVnIdwA', 'type': 'dIrr', 'notes': 'Canes Venatici I dwarf A'},
        {'name': 'DDO 43', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 46', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 47', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 50', 'type': 'dIrr', 'notes': 'DDO catalog (Ho II)'},
        {'name': 'DDO 52', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 53', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 63', 'type': 'dIrr', 'notes': 'DDO catalog (UGC 4459)'},
        {'name': 'DDO 70', 'type': 'dIrr', 'notes': 'DDO catalog (Sextans B)'},
        {'name': 'DDO 75', 'type': 'dIrr', 'notes': 'DDO catalog (UGCA 205)'},
        {'name': 'DDO 87', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 101', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 126', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 133', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 154', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 168', 'type': 'dIrr', 'notes': 'DDO catalog'},
        {'name': 'DDO 210', 'type': 'dIrr', 'notes': 'DDO catalog (Aquarius dwarf)'},
        {'name': 'DDO 216', 'type': 'dIrr', 'notes': 'DDO catalog (Pegasus dwarf)'},
        {'name': 'F564-V3', 'type': 'dIrr', 'notes': 'Low surface brightness'},
        {'name': 'IC 1613', 'type': 'dIrr', 'notes': 'Local Group member'},
        {'name': 'NGC 1569', 'type': 'dIrr', 'notes': 'Starburst dwarf'},
        {'name': 'NGC 2366', 'type': 'dIrr', 'notes': 'NGC catalog'},
        {'name': 'UGC 8508', 'type': 'dIrr', 'notes': 'UGC catalog'},
        {'name': 'WLM', 'type': 'dIrr', 'notes': 'Wolf-Lundmark-Melotte'},
        {'name': 'Haro 29', 'type': 'BCG', 'notes': 'Blue compact galaxy'},
        {'name': 'Haro 36', 'type': 'BCG', 'notes': 'Blue compact galaxy'}
    ]

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                    LITTLE THINGS SAMPLE: 26 DWARF GALAXIES                  │
└─────────────────────────────────────────────────────────────────────────────┘

Galaxy sample from Oh et al. (2015) AJ 149, 180:

┌────────────────┬─────────────────────────────────────────────────────────────┐
│ Galaxy         │ Notes                                                       │
├────────────────┼─────────────────────────────────────────────────────────────┤
""")

    for g in little_things_sample:
        print(f"│ {g['name']:<14} │ {g['notes']:<59} │")

    print("""└────────────────┴─────────────────────────────────────────────────────────────┘

SAMPLE CHARACTERISTICS:
══════════════════════════════════════════════════════════════════════════════

    - Primarily dwarf irregular (dIrr) galaxies
    - A few blue compact galaxies (BCG)
    - All within 11 Mpc
    - High-resolution HI rotation curves
    - Multi-wavelength photometry available

EXPECTED v_max RANGE:
══════════════════════════════════════════════════════════════════════════════

    Based on dwarf irregular properties:
    - Typical v_max: 20-80 km/s
    - Range: 10-100 km/s

    This is IDEAL for Synchronism validation!
    (Best performance in v_max < 100 km/s regime)

""")

    return little_things_sample


def create_data_template():
    """
    Create template for LITTLE THINGS data ingestion.
    """

    print("\n" + "="*80)
    print("DATA TEMPLATE FOR LITTLE THINGS")
    print("="*80)

    template = {
        'galaxy_name': 'DDO_154',  # Example
        'survey': 'LITTLE_THINGS',
        'reference': 'Oh et al. (2015) AJ 149, 180',
        'distance_Mpc': 3.7,
        'inclination_deg': 66.0,
        'v_max_km_s': 47.0,
        'rotation_curve': {
            'radius_kpc': [0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0],
            'v_obs_km_s': [10, 20, 28, 35, 42, 45, 47],
            'v_err_km_s': [3, 3, 3, 3, 3, 3, 3]
        },
        'mass_model': {
            'radius_kpc': [0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0],
            'v_star_km_s': [5, 8, 10, 12, 14, 15, 15],
            'v_gas_km_s': [3, 6, 10, 15, 20, 22, 23],
            'v_DM_inferred_km_s': [8, 17, 24, 30, 35, 38, 40]
        },
        'surface_density': {
            'radius_kpc': [0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0],
            'sigma_star_Msun_pc2': [50, 40, 30, 20, 10, 5, 2],
            'sigma_gas_Msun_pc2': [10, 10, 10, 9, 7, 5, 3]
        }
    }

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                    DATA TEMPLATE STRUCTURE                                  │
└─────────────────────────────────────────────────────────────────────────────┘

Each galaxy will be stored as a JSON file with this structure:

{
    "galaxy_name": "DDO_154",
    "survey": "LITTLE_THINGS",
    "reference": "Oh et al. (2015) AJ 149, 180",
    "distance_Mpc": 3.7,
    "inclination_deg": 66.0,
    "v_max_km_s": 47.0,

    "rotation_curve": {
        "radius_kpc": [...],
        "v_obs_km_s": [...],
        "v_err_km_s": [...]
    },

    "mass_model": {
        "radius_kpc": [...],
        "v_star_km_s": [...],
        "v_gas_km_s": [...],
        "v_DM_inferred_km_s": [...]
    },

    "surface_density": {
        "radius_kpc": [...],
        "sigma_star_Msun_pc2": [...],
        "sigma_gas_Msun_pc2": [...]
    }
}


DATA ACQUISITION STRATEGY:
══════════════════════════════════════════════════════════════════════════════

    OPTION 1: Digitize from Oh et al. (2015) tables
              - Published tables in the paper
              - Supplementary data files

    OPTION 2: VizieR catalog search
              - CDS/VizieR hosts many survey catalogs
              - May have machine-readable tables

    OPTION 3: Contact authors
              - Se-Heon Oh (lead author)
              - Request data products directly

    OPTION 4: VLA archive
              - Raw HI data cubes
              - Would require re-reduction

    RECOMMENDED: Start with VizieR/CDS for published tables.

""")

    return template


def outline_validation_pipeline():
    """
    Outline the validation pipeline.
    """

    print("\n" + "="*80)
    print("VALIDATION PIPELINE")
    print("="*80)

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                    SESSION #47+ VALIDATION PIPELINE                         │
└─────────────────────────────────────────────────────────────────────────────┘

STEP 1: DATA ACQUISITION (Session #47)
══════════════════════════════════════════════════════════════════════════════

    [ ] Search VizieR for LITTLE THINGS catalog
    [ ] Download rotation curve tables
    [ ] Download surface density profiles
    [ ] Parse into standard JSON format
    [ ] Validate data quality (NaN checks, unit conversions)


STEP 2: MODEL IMPLEMENTATION (Session #47)
══════════════════════════════════════════════════════════════════════════════

    [ ] Adapt SPARC validation code for LITTLE THINGS format
    [ ] Implement surface density → volume density conversion
    [ ] Test on DDO 154 (well-studied benchmark)


STEP 3: BLIND VALIDATION (Session #48)
══════════════════════════════════════════════════════════════════════════════

    [ ] Run 0-parameter Synchronism model on all 26 galaxies
    [ ] Compute χ² for each galaxy
    [ ] Apply success criterion (χ² < 5.0)
    [ ] Report success rate WITHOUT looking at results first


STEP 4: ANALYSIS (Session #48)
══════════════════════════════════════════════════════════════════════════════

    [ ] Compare to SPARC performance
    [ ] Analyze failure population (if any)
    [ ] Test if v_max correlation holds
    [ ] Compare to ΛCDM predictions from Oh et al. (2015)


STEP 5: PUBLICATION PREPARATION (Session #49+)
══════════════════════════════════════════════════════════════════════════════

    [ ] Write methods section
    [ ] Prepare comparison figures
    [ ] Draft arXiv preprint
    [ ] Address any systematic biases


PREDICTION FOR LITTLE THINGS:
══════════════════════════════════════════════════════════════════════════════

    Based on SPARC dwarf galaxy performance:

    ┌─────────────────────────────────────────────────────────────────┐
    │                                                                 │
    │  PREDICTED SUCCESS RATE: 60-70%                                 │
    │                                                                 │
    │  Confidence: HIGH (dwarf regime is Synchronism's strength)      │
    │                                                                 │
    │  If success < 50%: Model needs revision                         │
    │  If success > 70%: Model validated externally                   │
    │                                                                 │
    └─────────────────────────────────────────────────────────────────┘

""")


def save_results():
    """Save preparation results."""

    little_things_galaxies = [
        'CVnIdwA', 'DDO 43', 'DDO 46', 'DDO 47', 'DDO 50', 'DDO 52', 'DDO 53',
        'DDO 63', 'DDO 70', 'DDO 75', 'DDO 87', 'DDO 101', 'DDO 126', 'DDO 133',
        'DDO 154', 'DDO 168', 'DDO 210', 'DDO 216', 'F564-V3', 'IC 1613',
        'NGC 1569', 'NGC 2366', 'UGC 8508', 'WLM', 'Haro 29', 'Haro 36'
    ]

    output = {
        'session': 46,
        'track': 'C - LITTLE THINGS Preparation',
        'date': datetime.now().isoformat(),
        'survey': {
            'name': 'LITTLE THINGS',
            'full_name': 'Local Irregulars That Trace Luminosity Extremes, The H I Nearby Galaxy Survey',
            'reference': 'Oh et al. (2015) AJ 149, 180',
            'arxiv': '1502.01281',
            'n_galaxies': 26,
            'telescope': 'VLA',
            'wavelength': 'HI 21 cm',
            'volume': '< 11 Mpc'
        },
        'galaxy_sample': little_things_galaxies,
        'validation_framework': {
            'model': 'Synchronism 0-parameter',
            'parameters': {
                'gamma': 2.0,
                'B': 1.62,
                'beta': 0.30,
                'A': 0.25
            },
            'success_criterion': 'chi2 < 5.0'
        },
        'predictions': {
            'expected_success_rate': '60-70%',
            'basis': 'SPARC dwarf results (67% for v_max < 100, 81.8% for v_max < 50)',
            'confidence': 'high',
            'failure_threshold': '50%',
            'validation_threshold': '70%'
        },
        'next_steps': [
            'Session #47: Data acquisition from VizieR/CDS',
            'Session #47: Implement validation pipeline',
            'Session #48: Blind validation test',
            'Session #48: Analysis and comparison',
            'Session #49: Publication preparation'
        ],
        'status': 'Framework prepared, awaiting data acquisition'
    }

    output_path = Path(__file__).parent / 'session46_little_things_preparation_results.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    return output


if __name__ == '__main__':
    print("\n" + "="*80)
    print("SESSION #46 TRACK C: LITTLE THINGS EXTERNAL VALIDATION PREPARATION")
    print("="*80)

    document_little_things_survey()
    define_validation_framework()
    list_known_little_things_galaxies()
    create_data_template()
    outline_validation_pipeline()
    save_results()

    print("\n" + "="*80)
    print("SESSION #46 TRACK C COMPLETE")
    print("="*80)
    print("\nCONCLUSION: LITTLE THINGS validation framework prepared.")
    print("26 dwarf galaxies identified for external validation.")
    print("Predicted success rate: 60-70% (based on SPARC dwarf performance)")
    print("\nNEXT STEP: Data acquisition from VizieR/CDS (Session #47)")
