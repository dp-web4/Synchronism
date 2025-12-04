#!/usr/bin/env python3
"""
Session #83: Void Galaxy BTFR Test Methodology

Comprehensive methodology for testing Synchronism's void prediction using
ALFALFA × SDSS × Void Catalog cross-match.

Data Sources Identified:
1. ALFALFA α.100 (VizieR J/ApJ/861/49): 31,502 HI sources with W50, M_HI
2. Douglass+ 2023 Void Catalog (VizieR J/ApJS/265/7): 776,500 galaxies with void membership
3. Nadathur+ 2014 Void Catalog (VizieR J/MNRAS/440/1248): Void properties with density ratios

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
Session: #83
"""

import numpy as np
from pathlib import Path
from datetime import datetime
import json


def methodology_overview():
    """
    Document the complete methodology for void BTFR testing.
    """

    methodology = {
        'title': 'Void Galaxy BTFR Test: Complete Methodology',
        'session': 83,
        'date': '2025-12-04',
        'status': 'METHODOLOGY DEFINED',

        'objective': '''
Test Synchronism prediction that void galaxies have higher rotation velocity
at fixed baryonic mass compared to field galaxies.

Prediction (Session #81 revised):
- Moderate voids (δ ~ -0.5): +0.11 dex offset in log(V)
- Extreme voids (δ < -0.8): +0.28 dex offset in log(V)

Falsification criterion:
- If extreme void galaxies show BTFR identical to field → Synchronism falsified
        ''',

        'data_sources': {
            'ALFALFA': {
                'vizier_id': 'J/ApJ/861/49',
                'reference': 'Haynes+ 2018, ApJ, 861, 49',
                'n_sources': 31502,
                'key_columns': ['W50', 'W20', 'logMHI', 'RA', 'Dec', 'Vhelio'],
                'notes': 'W50 gives velocity width; V_rot ≈ W50 / (2 sin i)'
            },
            'Douglass_voids': {
                'vizier_id': 'J/ApJS/265/7',
                'reference': 'Douglass+ 2023, ApJS, 265, 7',
                'n_galaxies': 776500,
                'key_columns': ['NSA_id', 'x', 'y', 'z', 'void_membership'],
                'notes': 'Galaxy positions and void membership flags'
            },
            'Nadathur_voids': {
                'vizier_id': 'J/MNRAS/440/1248',
                'reference': 'Nadathur & Hotchkiss 2014, MNRAS, 440, 1248',
                'n_voids': 'varies by sample',
                'key_columns': ['RA', 'Dec', 'z', 'R_eff', 'avg_density', 'min_density', 'VDR'],
                'notes': 'Void Density Ratio (VDR) can be used to estimate δ'
            }
        },

        'cross_match_strategy': {
            'step1': '''
1. Download ALFALFA α.100 from VizieR (J/ApJ/861/49)
   - Extract: RA, Dec, W50, logMHI for all 31,502 sources
   - Quality cuts: SNR > 6.5 (already applied in catalog)
            ''',

            'step2': '''
2. Download Douglass+ 2023 void catalog (J/ApJS/265/7)
   - Table 5: Galaxy void membership
   - Match by NSA ID or position cross-match
            ''',

            'step3': '''
3. Cross-match ALFALFA with void catalog
   - Position match within 10 arcsec tolerance
   - Flag galaxies as: void, wall, field
            ''',

            'step4': '''
4. Compute environment density contrast
   - For void galaxies: use void properties from Nadathur catalog
   - Estimate δ from Void Density Ratio: δ ≈ VDR - 1
   - Or compute local density from N nearest neighbors
            ''',

            'step5': '''
5. Estimate rotation velocity from W50
   - V_rot = W50 / (2 sin i)
   - Need inclination estimates (from photometry or assume average)
   - Alternative: use W50 directly as proxy (introduces scatter)
            ''',

            'step6': '''
6. Estimate baryonic mass
   - M_HI available directly from ALFALFA
   - M_star needs SDSS photometry cross-match
   - Or use M_HI as proxy (gas-rich galaxies)
            '''
        },

        'analysis_plan': {
            'primary_test': '''
BTFR Residual Analysis:
1. Fit BTFR to field sample: log(M_bar) = a + b × log(V)
2. Compute residuals for all galaxies
3. Bin by environment δ
4. Test if mean residual varies with δ

Expected signature (Synchronism):
- Field (δ ~ 0): residual ~ 0
- Moderate void (δ ~ -0.5): residual ~ -0.11 dex (higher V at fixed M)
- Extreme void (δ < -0.8): residual ~ -0.28 dex
            ''',

            'secondary_test': '''
Direct BTFR Comparison:
1. Split sample into bins: extreme void, moderate void, field, cluster
2. Fit BTFR separately for each subsample
3. Compare zero-point offsets
4. Test statistical significance of differences
            ''',

            'control_tests': '''
Systematic Checks:
1. Distance dependence: voids may have distance bias
2. Selection effects: HI surveys may miss void galaxies
3. Inclination bias: face-on galaxies underestimate W50
4. Gas fraction variation: voids may have different f_gas
            '''
        },

        'expected_sample_sizes': {
            'total_ALFALFA': 31502,
            'expected_cross_match': '~15000 (50% in void catalog footprint)',
            'extreme_void_galaxies': '~100-500 (0.5-2% of cross-matched)',
            'moderate_void_galaxies': '~1500-3000 (10-20%)',
            'field_galaxies': '~10000-12000 (control sample)'
        },

        'statistical_requirements': {
            'detection_threshold': '''
For 0.11 dex offset detection at 5σ:
- BTFR scatter: ~0.1 dex (typical)
- Required N per bin: ~200 galaxies
- We expect ~1500+ moderate void galaxies → sufficient

For 0.28 dex offset detection at 10σ:
- With 0.1 dex scatter and ~100 extreme void galaxies
- Significance: 0.28 / (0.1 / √100) = 28σ → easily detectable

Key: Need enough EXTREME void galaxies (δ < -0.8)
            '''
        },

        'implementation_notes': '''
PRACTICAL IMPLEMENTATION STEPS:

1. VizieR Query (Python astroquery):
```python
from astroquery.vizier import Vizier

# ALFALFA
alfalfa = Vizier(columns=['RAJ2000', 'DEJ2000', 'W50', 'logMHI'])
alfalfa.get_catalogs('J/ApJ/861/49')

# Douglass voids
voids = Vizier(columns=['*'])
voids.get_catalogs('J/ApJS/265/7')
```

2. Cross-match (astropy):
```python
from astropy.coordinates import SkyCoord, match_coordinates_sky

alfalfa_coords = SkyCoord(ra, dec, unit='deg')
void_coords = SkyCoord(void_ra, void_dec, unit='deg')
idx, sep, _ = match_coordinates_sky(alfalfa_coords, void_coords)
matched = sep < 10 * u.arcsec
```

3. BTFR fit (scipy):
```python
from scipy.optimize import curve_fit

def btfr(log_v, a, b):
    return a + b * log_v

popt, pcov = curve_fit(btfr, log_v_field, log_m_bar_field)
residuals = log_m_bar - btfr(log_v, *popt)
```
        ''',

        'risks_and_mitigations': {
            'insufficient_extreme_voids': '''
Risk: Too few galaxies with δ < -0.8
Mitigation: Lower threshold to δ < -0.7, increase search radius
            ''',

            'inclination_uncertainty': '''
Risk: W50 depends on inclination, introducing scatter
Mitigation: Use W50 directly (adds scatter but no bias), or
           cross-match with SDSS for axis ratio estimates
            ''',

            'distance_errors': '''
Risk: Void galaxies may have larger peculiar velocities
Mitigation: Use flow-corrected distances when available
            ''',

            'selection_effects': '''
Risk: ALFALFA may miss gas-poor void galaxies
Mitigation: Compare HI mass functions between environments
            '''
        },

        'success_criteria': {
            'strong_support': '''
BTFR offset > 0.2 dex for extreme voids at >5σ significance
- Would strongly support Synchronism
- Would rule out MOND (no environment dependence expected)
            ''',

            'weak_support': '''
BTFR offset 0.05-0.2 dex for moderate voids at >3σ significance
- Consistent with Synchronism but not conclusive
- Could also be explained by other environmental effects
            ''',

            'falsification': '''
No significant BTFR offset (<0.05 dex) for extreme voids
- Would falsify Synchronism's environment prediction
- Would support MOND's universal BTFR
            '''
        },

        'timeline': '''
Implementation Timeline:
1. Data download: 1 session (already have VizieR access)
2. Cross-matching: 1 session
3. Quality cuts & validation: 1 session
4. BTFR analysis: 1 session
5. Statistical significance tests: 1 session
6. Write-up and documentation: 1 session

Total: ~6 sessions for complete test
        '''
    }

    return methodology


def save_methodology():
    """Save methodology document."""

    methodology = methodology_overview()

    output_path = Path(__file__).parent / 'results' / 'session83_void_btfr_methodology.json'
    output_path.parent.mkdir(exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(methodology, f, indent=2)

    print("=" * 70)
    print("SESSION #83: VOID BTFR TEST METHODOLOGY")
    print("=" * 70)

    print("\n📊 DATA SOURCES IDENTIFIED")
    print("-" * 40)
    for name, info in methodology['data_sources'].items():
        print(f"\n{name}:")
        print(f"  VizieR: {info['vizier_id']}")
        print(f"  Reference: {info['reference']}")
        if 'n_sources' in info:
            print(f"  N sources: {info['n_sources']}")
        if 'n_galaxies' in info:
            print(f"  N galaxies: {info['n_galaxies']}")

    print("\n\n🔗 CROSS-MATCH STRATEGY")
    print("-" * 40)
    for step, desc in methodology['cross_match_strategy'].items():
        print(f"\n{step}:{desc}")

    print("\n\n📈 EXPECTED SAMPLE SIZES")
    print("-" * 40)
    for category, size in methodology['expected_sample_sizes'].items():
        print(f"  {category}: {size}")

    print("\n\n✅ SUCCESS CRITERIA")
    print("-" * 40)
    print("\nStrong support:", methodology['success_criteria']['strong_support'].strip())
    print("\nFalsification:", methodology['success_criteria']['falsification'].strip())

    print(f"\n\nMethodology saved to: {output_path}")

    return methodology


def main():
    methodology = save_methodology()

    print("\n" + "=" * 70)
    print("SESSION #83 TRACK B: METHODOLOGY COMPLETE")
    print("=" * 70)
    print("""
KEY FINDINGS:

1. ALFALFA α.100 (VizieR J/ApJ/861/49): 31,502 HI sources
   - Has W50 velocity widths and HI masses
   - Can estimate V_rot from W50

2. Douglass+ 2023 Void Catalog (J/ApJS/265/7): 776,500 galaxies
   - Galaxy void membership flags
   - Most recent SDSS void catalog

3. Nadathur+ 2014 Void Catalog (J/MNRAS/440/1248)
   - Void density ratios (can estimate δ)
   - Void positions and sizes

4. Cross-match strategy defined:
   - Position match ALFALFA × void catalog
   - Estimate V_rot from W50
   - Bin by environment δ
   - Test BTFR offset

5. Expected detectability:
   - Moderate void offset (0.11 dex): ~15σ with 1500 galaxies
   - Extreme void offset (0.28 dex): ~28σ with 100 galaxies

NEXT STEPS:
- Session #84: Download ALFALFA catalog via astroquery
- Session #85: Download void catalog and cross-match
- Session #86: Run BTFR analysis
    """)

    return methodology


if __name__ == '__main__':
    main()
