#!/usr/bin/env python3
"""
Synchronism Research - Session #25
B-Field Literature Analysis for SPARC Galaxies

Following Nova Session #24 recommendation: Test if NGC has systematically
weaker B-fields than F/irregulars.

HYPOTHESIS: If NGC underprediction is physical (not calculation bias), it may
be due to weaker B-fields in NGC spirals → less magnetic screening → lower ρ_sat.

CRITICAL TEST: Do published B-field measurements support this hypothesis?
"""

import numpy as np
import pickle
from pathlib import Path

# ============================================================================
# PART 1: B-FIELD LITERATURE COMPILATION
# ============================================================================

def compile_bfield_literature():
    """
    Compile published B-field measurements from literature review.

    Sources:
    - Niklas 1995: 74 spiral galaxies, <B> = 9 ± 2 µG
    - Fletcher 2010: 21 bright galaxies, <B> = 17 ± 3 µG
    - Chyży et al. 2011: Local Group dwarfs, <B> = 4.2 ± 1.8 µG
    - Beck 2015: Review compilation
    - Van Eck et al. 2015: Table 2 (SFR-B correlation)
    """

    print("="*80)
    print("SESSION #25: B-FIELD LITERATURE COMPILATION")
    print("="*80)
    print()

    # Population-level statistics
    populations = {
        'Spiral galaxies (Niklas 1995)': {
            'n': 74,
            'B_mean': 9.0,  # µG
            'B_std': 2.0,
            'types': ['NGC', 'UGC', 'IC', 'M'],
            'notes': 'Sample average, 1995 data'
        },
        'Bright spirals (Fletcher 2010)': {
            'n': 21,
            'B_mean': 17.0,  # µG (total field)
            'B_std': 3.0,
            'B_ord_mean': 5.0,  # µG (ordered field)
            'B_ord_std': 3.0,
            'types': ['NGC', 'M', 'IC'],
            'notes': '2000-2010 observations, radio-bright'
        },
        'Dwarf irregulars (Chyży+ 2011)': {
            'n': 'LG sample',
            'B_mean': 4.2,  # µG
            'B_std': 1.8,
            'types': ['DDO', 'F', 'IC'],
            'notes': 'Local Group dwarfs, 3× weaker than spirals'
        }
    }

    print("POPULATION-LEVEL STATISTICS")
    print("-" * 80)
    for pop_name, pop_data in populations.items():
        print(f"\n{pop_name}:")
        print(f"  Sample size: {pop_data['n']}")
        print(f"  <B_total> = {pop_data['B_mean']:.1f} ± {pop_data['B_std']:.1f} µG")
        if 'B_ord_mean' in pop_data:
            print(f"  <B_ordered> = {pop_data['B_ord_mean']:.1f} ± {pop_data['B_ord_std']:.1f} µG")
        print(f"  Galaxy types: {', '.join(pop_data['types'])}")
        print(f"  Notes: {pop_data['notes']}")

    print("\n" + "="*80)
    print()

    # Individual galaxy measurements from literature
    # Format: (name, B_total_µG, B_ordered_µG, type, reference, notes)
    galaxies_with_bfield = [
        # Spirals - NGC
        ('NGC 6946', 20.0, 10.0, 'NGC', 'Beck+ 2015', 'High SFR, spiral arms'),
        ('NGC 253', 75.0, None, 'NGC', 'Beck+ 2015', 'Starburst'),
        ('NGC 4254', 15.0, None, 'NGC', 'Chyży+ 2007', 'Virgo cluster spiral'),
        ('NGC 4414', 12.0, None, 'NGC', 'Beck+ 2002', 'Flocculent spiral'),
        ('NGC 2997', 8.0, 4.0, 'NGC', 'Han+ 1999', 'Grand design spiral'),
        ('NGC 4736', 10.0, None, 'NGC', 'Chyży+ 2011', 'Ring galaxy'),
        ('NGC 5775', 15.0, None, 'NGC', 'Soida+ 2011', 'Edge-on spiral'),

        # Messier spirals (often NGC)
        ('M 31', 6.0, 5.0, 'NGC', 'Fletcher+ 2004', 'Andromeda, radio-faint'),
        ('M 33', 6.0, None, 'NGC', 'Tabatabaei+ 2008', 'Triangulum, radio-faint'),
        ('M 51', 25.0, 12.0, 'NGC', 'Fletcher+ 2011', 'High SFR, interacting'),
        ('M 81', 8.0, None, 'NGC', 'Krause+ 1989', 'NGC 3031'),
        ('M 82', 85.0, None, 'NGC', 'Adebahr+ 2013', 'Starburst'),
        ('M 83', 25.0, None, 'NGC', 'Heesen+ 2011', 'High SFR'),

        # Irregulars
        ('IC 342', 7.0, None, 'IC', 'Beck 2015', 'Nearby irregular'),
        ('NGC 4449', 14.0, 7.0, 'NGC', 'Chyży+ 2000', 'Magellanic irregular, starburst'),
        ('NGC 1569', 12.0, None, 'NGC', 'Kepley+ 2010', 'Dwarf starburst'),

        # Dwarf irregulars
        ('IC 10', 10.0, None, 'IC', 'Chyży+ 2003', 'Strongest LG dwarf, starburst'),
        ('DDO 154', 3.0, None, 'DDO', 'Chyży+ 2011', 'Typical dwarf'),
        ('Ho II', 4.0, None, 'F', 'Chyży+ 2011', 'Dwarf irregular'),
        ('IC 2574', 5.0, None, 'IC', 'Chyży+ 2016', 'Dwarf irregular'),
        ('NGC 2366', 8.0, None, 'NGC', 'Chyży+ 2011', 'Dwarf irregular, moderate SFR'),

        # Interacting/peculiar
        ('NGC 4038/9', 75.0, None, 'NGC', 'Chyży+ 2004', 'Antennae, starburst merger'),
    ]

    print("INDIVIDUAL GALAXY B-FIELD MEASUREMENTS")
    print("-" * 80)
    print(f"{'Galaxy':<15} {'B_tot (µG)':<12} {'B_ord (µG)':<12} {'Type':<8} {'Reference':<20} {'Notes'}")
    print("-" * 80)

    for gal_data in galaxies_with_bfield:
        name, B_tot, B_ord, gtype, ref, notes = gal_data
        B_ord_str = f"{B_ord:.1f}" if B_ord is not None else "---"
        print(f"{name:<15} {B_tot:<12.1f} {B_ord_str:<12} {gtype:<8} {ref:<20} {notes}")

    print()
    print("="*80)
    print()

    return populations, galaxies_with_bfield


# ============================================================================
# PART 2: GALAXY TYPE CLASSIFICATION AND ANALYSIS
# ============================================================================

def analyze_by_galaxy_type(galaxies_with_bfield):
    """
    Classify galaxies by SPARC naming convention and compare B-field strengths.

    SPARC galaxy types (from Session #22):
    - NGC: Spirals (n=63)
    - UGC: Spirals (n=79)
    - F: Irregulars (n=16)
    - DDO: Dwarfs (n=5)
    - Other: IC, M, etc. (n=12)
    """

    print("B-FIELD STRENGTH BY GALAXY TYPE")
    print("="*80)
    print()

    # Classify by type
    type_bfields = {
        'NGC_spiral': [],
        'NGC_starburst': [],
        'NGC_dwarf': [],
        'UGC': [],
        'F_irregular': [],
        'DDO': [],
        'IC': [],
        'M': []
    }

    for gal_data in galaxies_with_bfield:
        name, B_tot, B_ord, gtype, ref, notes = gal_data

        # Classify by notes and type
        is_starburst = 'starburst' in notes.lower() or 'merger' in notes.lower()
        is_dwarf = 'dwarf' in notes.lower()
        is_irregular = 'irregular' in notes.lower()

        if gtype == 'NGC':
            if is_starburst:
                type_bfields['NGC_starburst'].append(B_tot)
            elif is_dwarf or is_irregular:
                type_bfields['NGC_dwarf'].append(B_tot)
            else:
                type_bfields['NGC_spiral'].append(B_tot)
        elif gtype == 'F':
            type_bfields['F_irregular'].append(B_tot)
        elif gtype == 'DDO':
            type_bfields['DDO'].append(B_tot)
        elif gtype == 'IC':
            type_bfields['IC'].append(B_tot)
        elif name.startswith('M '):
            type_bfields['M'].append(B_tot)

    # Analyze statistics
    print("STATISTICS BY TYPE")
    print("-" * 80)
    print(f"{'Type':<20} {'n':<6} {'<B> (µG)':<12} {'Median (µG)':<12} {'Range (µG)'}")
    print("-" * 80)

    results = {}
    for type_name, B_values in type_bfields.items():
        if len(B_values) > 0:
            B_mean = np.mean(B_values)
            B_median = np.median(B_values)
            B_min = np.min(B_values)
            B_max = np.max(B_values)

            print(f"{type_name:<20} {len(B_values):<6} {B_mean:<12.1f} {B_median:<12.1f} [{B_min:.1f}, {B_max:.1f}]")

            results[type_name] = {
                'n': len(B_values),
                'mean': B_mean,
                'median': B_median,
                'min': B_min,
                'max': B_max,
                'values': B_values
            }

    print()
    print("="*80)
    print()

    return results


# ============================================================================
# PART 3: CRITICAL HYPOTHESIS TEST
# ============================================================================

def test_ngc_weaker_bfield_hypothesis(type_results):
    """
    Test hypothesis: NGC spirals have WEAKER B-fields than F/DDO irregulars.

    Session #22 found:
    - NGC galaxies: ρ_sat observed ~ 2×10³ M☉/pc³
    - F galaxies: ρ_sat observed ~ 5×10⁶ M☉/pc³
    - F/NGC ratio: 2336×

    Magnetic screening model: ρ_sat = ρ_sat,0 / [1 + (ρ_c/ρ_mag)^δ]

    If B-field weaker → ρ_mag lower → screening MORE effective → ρ_sat LOWER

    HYPOTHESIS: NGC has weaker B-fields than F → explains NGC underprediction

    TEST: Compare <B_NGC> vs <B_F> and <B_DDO>
    """

    print("CRITICAL HYPOTHESIS TEST")
    print("="*80)
    print()

    print("Session #22 NGC Underprediction:")
    print("  NGC observed: ρ_sat ~ 2×10³ M☉/pc³")
    print("  F observed:   ρ_sat ~ 5×10⁶ M☉/pc³")
    print("  F/NGC ratio:  2336×")
    print()

    print("Hypothesis: NGC underprediction due to WEAKER B-fields in NGC spirals")
    print("  If B_NGC < B_F → ρ_mag,NGC < ρ_mag,F → screening stronger → ρ_sat,NGC lower ✓")
    print()

    print("Magnetic screening model:")
    print("  ρ_sat = ρ_sat,0 / [1 + (ρ_c/ρ_mag)^δ]")
    print("  ρ_mag ∝ B² (magnetic energy density)")
    print()

    # Compare NGC spirals vs F/DDO
    B_NGC_spiral = type_results.get('NGC_spiral', {}).get('median', None)
    B_F = type_results.get('F_irregular', {}).get('median', None)
    B_DDO = type_results.get('DDO', {}).get('median', None)

    # Use population means if individual samples too small
    if B_F is None:
        B_F = 4.2  # From Chyży+ 2011 LG dwarf mean
        print("  Using population mean for F/DDO: 4.2 µG (Chyży+ 2011)")

    if B_NGC_spiral is None:
        B_NGC_spiral = 9.0  # From Niklas 1995 spiral mean
        print("  Using population mean for NGC spirals: 9.0 µG (Niklas 1995)")

    print()
    print("RESULTS:")
    print("-" * 80)
    print(f"  <B_NGC_spiral> = {B_NGC_spiral:.1f} µG")
    print(f"  <B_F/DDO>      = {B_F:.1f} µG")
    print()

    ratio_B = B_NGC_spiral / B_F
    print(f"  B_NGC / B_F = {ratio_B:.2f}")
    print()

    # Expected ρ_sat ratio if B-field is the cause
    # ρ_mag ∝ B² → if B_NGC = 2× B_F, then ρ_mag,NGC = 4× ρ_mag,F
    # ρ_sat ∝ 1/[1 + (ρ_c/ρ_mag)^δ]
    # If ρ_mag increases, ρ_sat increases (less screening)

    ratio_rho_mag = ratio_B**2
    print(f"  ρ_mag,NGC / ρ_mag,F ∝ (B_NGC/B_F)² = {ratio_rho_mag:.2f}")
    print()

    # Critical interpretation
    print("INTERPRETATION:")
    print("-" * 80)

    if B_NGC_spiral > B_F:
        print("  ❌ HYPOTHESIS FALSIFIED")
        print()
        print("  Observation: NGC spirals have STRONGER B-fields than F/DDO dwarfs")
        print(f"    B_NGC = {B_NGC_spiral:.1f} µG > B_F = {B_F:.1f} µG")
        print()
        print("  Implication: Stronger B-field → Higher ρ_mag → LESS screening")
        print("    → NGC should have HIGHER ρ_sat, not lower!")
        print()
        print("  Conclusion: Magnetic screening with universal ρ_sat,0 CANNOT")
        print("              explain NGC underprediction via B-field variation.")
        print()
        print("  Physical interpretation:")
        print("    - NGC spirals: Strong B-field (~9 µG) → weak screening → high ρ_sat")
        print("    - F/DDO dwarfs: Weak B-field (~4 µG) → strong screening → low ρ_sat")
        print("    - Observed: OPPOSITE (NGC low, F high)")
        print("    - Paradox: Model predicts inverse of observation!")

    elif B_NGC_spiral < B_F:
        print("  ✓ HYPOTHESIS SUPPORTED")
        print()
        print("  Observation: NGC spirals have WEAKER B-fields than F/DDO")
        print(f"    B_NGC = {B_NGC_spiral:.1f} µG < B_F = {B_F:.1f} µG")
        print()
        print("  Implication: Weaker B-field → Lower ρ_mag → MORE screening")
        print("    → NGC should have LOWER ρ_sat ✓")
        print()
        print("  Conclusion: B-field variation CAN explain NGC underprediction")

    else:
        print("  ⚠️ HYPOTHESIS INCONCLUSIVE")
        print()
        print("  Observation: NGC and F/DDO have similar B-fields")
        print("    → B-field variation cannot explain ρ_sat difference")

    print()
    print("="*80)
    print()

    return {
        'B_NGC': B_NGC_spiral,
        'B_F': B_F,
        'ratio_B': ratio_B,
        'ratio_rho_mag': ratio_rho_mag,
        'hypothesis_supported': B_NGC_spiral < B_F
    }


# ============================================================================
# PART 4: THEORETICAL IMPLICATIONS
# ============================================================================

def analyze_theoretical_implications(test_result):
    """
    Analyze implications for Synchronism magnetic screening model.
    """

    print("THEORETICAL IMPLICATIONS FOR SYNCHRONISM")
    print("="*80)
    print()

    print("Session #22 Magnetic Screening Model:")
    print("  ρ_sat = ρ_sat,0 / [1 + (ρ_central/ρ_mag)^δ]")
    print()
    print("  Parameters (universal fit, R² = 0.406):")
    print("    ρ_sat,0 = (6.93 ± 0.65) × 10⁵ M☉/pc³")
    print("    ρ_mag   = (2.88 ± 0.66) × 10² M☉/pc³")
    print("    δ       = 1.85 ± 0.60")
    print()

    if not test_result['hypothesis_supported']:
        print("FUNDAMENTAL PROBLEM IDENTIFIED:")
        print("-" * 80)
        print()
        print("Literature shows NGC spirals have STRONGER B-fields than F/DDO dwarfs:")
        print(f"  B_NGC ~ {test_result['B_NGC']:.1f} µG")
        print(f"  B_F   ~ {test_result['B_F']:.1f} µG")
        print(f"  Ratio: {test_result['ratio_B']:.2f}")
        print()
        print("Magnetic screening physics:")
        print("  Stronger B-field → Higher ρ_mag → LESS screening → HIGHER ρ_sat")
        print()
        print("Expected from B-field data:")
        print("  NGC should have HIGHER ρ_sat than F (opposite of observation!)")
        print()
        print("Actual observation (Session #20):")
        print("  NGC: ρ_sat ~ 2×10³ M☉/pc³ (LOW)")
        print("  F:   ρ_sat ~ 5×10⁶ M☉/pc³ (HIGH)")
        print()
        print("PARADOX: Model predicts INVERSE of observation!")
        print()
        print("="*80)
        print()

        print("POSSIBLE RESOLUTIONS:")
        print("-" * 80)
        print()
        print("1. MAGNETIC SCREENING MODEL IS INCORRECT")
        print("   - Functional form ρ_sat ∝ 1/[1 + (ρ_c/ρ_mag)^δ] may be wrong")
        print("   - Alternative: ρ_sat ∝ [1 + (ρ_c/ρ_mag)^δ] (opposite dependence)?")
        print("   - Would require rederiving from Synchronism first principles")
        print()

        print("2. ρ_mag IS NOT PROPORTIONAL TO B²")
        print("   - Session #21 assumed ρ_mag ∝ B² (magnetic energy density)")
        print("   - Alternative: ρ_mag ∝ 1/B² (inverse relationship)?")
        print("   - Would contradict standard magnetic screening physics")
        print()

        print("3. ρ_sat,0 IS GALAXY-SPECIFIC, NOT UNIVERSAL")
        print("   - Model assumes universal ρ_sat,0 for all galaxies")
        print("   - If ρ_sat,0,NGC << ρ_sat,0,F, could explain NGC underprediction")
        print("   - But what would cause ρ_sat,0 variation? Not B-field related")
        print("   - Possible: Different coherence physics (Synchronism mechanism)")
        print()

        print("4. NGC GALAXIES HAVE ADDITIONAL SCREENING MECHANISM")
        print("   - Beyond magnetic screening")
        print("   - E.g., AGN feedback, disk thickness, rotation curve shape")
        print("   - Would require new term in model")
        print()

        print("5. ACCEPT NGC AS OUTLIERS")
        print("   - Majority fit (F/DDO/UGC irregulars) is good")
        print("   - NGC spirals are fundamentally different population")
        print("   - Model applies to high-ρ_sat galaxies, not low-ρ_sat spirals")
        print("   - Requires understanding why spirals are different")
        print()

        print("="*80)
        print()

        print("RECOMMENDED NEXT STEPS:")
        print("-" * 80)
        print()
        print("Priority 1: TEST INVERSE MAGNETIC SCREENING MODEL")
        print("  - Refit with ρ_sat ∝ [1 + (ρ_c/ρ_mag)^δ] instead")
        print("  - Check if R² improves and NGC prediction corrects")
        print("  - If yes: Rederive from Synchronism to justify form")
        print()

        print("Priority 2: TEST GALAXY-SPECIFIC ρ_sat,0")
        print("  - Fit ρ_sat,0 separately for NGC vs F/DDO")
        print("  - Check if physical (SFR, morphology, mass-dependent)")
        print("  - Explore Synchronism coherence mechanism differences")
        print()

        print("Priority 3: LITERATURE REVIEW - SFR CORRELATION")
        print("  - B ∝ SFR^0.3 is well-established")
        print("  - NGC spirals: Low SFR → Weaker B? (contradicts compilation)")
        print("  - F irregulars: High SFR → Stronger B? (contradicts compilation)")
        print("  - Resolve apparent contradiction")
        print()

    else:
        print("HYPOTHESIS SUPPORTED:")
        print("-" * 80)
        print()
        print("NGC spirals have weaker B-fields than F/DDO dwarfs")
        print("  → Lower ρ_mag → Stronger screening → Lower ρ_sat ✓")
        print()
        print("Next step: Implement galaxy-specific ρ_mag(B) in Session #22 model")
        print()

    print("="*80)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print()
    print("╔" + "="*78 + "╗")
    print("║" + " "*78 + "║")
    print("║" + "  Synchronism Research - Session #25".center(78) + "║")
    print("║" + "  B-Field Literature Analysis for SPARC Galaxies".center(78) + "║")
    print("║" + " "*78 + "║")
    print("╚" + "="*78 + "╝")
    print()

    # Part 1: Compile literature
    populations, galaxies_with_bfield = compile_bfield_literature()

    # Part 2: Analyze by type
    type_results = analyze_by_galaxy_type(galaxies_with_bfield)

    # Part 3: Critical hypothesis test
    test_result = test_ngc_weaker_bfield_hypothesis(type_results)

    # Part 4: Theoretical implications
    analyze_theoretical_implications(test_result)

    print()
    print("SESSION #25 COMPLETE")
    print("="*80)
    print()
    print("SUMMARY:")
    print(f"  - Compiled {len(galaxies_with_bfield)} individual B-field measurements")
    print(f"  - Analyzed {len(populations)} population-level compilations")
    print(f"  - Classified by SPARC galaxy types")
    print(f"  - Tested NGC weaker B-field hypothesis: {'SUPPORTED ✓' if test_result['hypothesis_supported'] else 'FALSIFIED ❌'}")
    print()

    if not test_result['hypothesis_supported']:
        print("CRITICAL FINDING:")
        print("  NGC spirals have STRONGER B-fields than F/DDO dwarfs")
        print("  → Magnetic screening CANNOT explain NGC underprediction")
        print("  → Model may have incorrect functional form or missing physics")
        print()
        print("Recommended: Test inverse screening model (Session #26)")
        print()

    print("Ready for Session #25 documentation and commit.")
    print()
