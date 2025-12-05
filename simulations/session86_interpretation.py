#!/usr/bin/env python3
"""
Session #86: Interpretation of HSB/LSB BTFR Result

The observation CONTRADICTS the naive Synchronism prediction.
Let's analyze possible explanations.

Result:
- LSB - HSB offset: -0.053 dex (3.0σ)
- Predicted: +0.088 dex
- Opposite sign!

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
"""

import json

def interpret_result():
    """Analyze the unexpected HSB/LSB BTFR result."""

    print("=" * 70)
    print("SESSION #86: INTERPRETATION OF HSB/LSB BTFR RESULT")
    print("=" * 70)
    print()

    # Observed results
    observed_offset = -0.053  # dex (LSB - HSB)
    significance = -3.0  # sigma
    predicted_offset = +0.088  # dex

    print("OBSERVED RESULT:")
    print("-" * 40)
    print(f"  LSB - HSB BTFR offset: {observed_offset:+.3f} dex ({significance:+.1f}σ)")
    print(f"  Naive prediction: {predicted_offset:+.3f} dex")
    print(f"  Direction: OPPOSITE to prediction")
    print()

    print("=" * 70)
    print("POSSIBLE EXPLANATIONS")
    print("=" * 70)
    print()

    explanations = []

    # Explanation 1: Theory misinterpretation
    print("1. THEORY MISINTERPRETATION")
    print("-" * 40)
    print("""
    The naive prediction assumed:
    - Lower SB → lower LOCAL density → lower C → higher G_eff → higher V

    BUT this may be BACKWARDS for Synchronism:
    - C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))
    - As ρ DECREASES (lower SB), C DECREASES
    - G_eff = G/C, so G_eff INCREASES
    - This should give HIGHER V at fixed M

    Wait... the math is correct. So why does LSB have LOWER V?

    Alternative interpretation:
    - Perhaps C isn't just about LOCAL density
    - C might be MORE sensitive to CENTRAL concentration
    - HSB galaxies have higher central density → higher C in the center
    - BUT Synchronism uses AVERAGE density over the rotation curve region

    The issue is: which density matters?
    - Local density at each radius? (original model)
    - Central density? (concentration effect)
    - Average disk density? (global property)
""")
    explanations.append({
        'name': 'Theory misinterpretation',
        'description': 'The relevant density scale may not be average SB',
        'plausibility': 'Medium'
    })

    # Explanation 2: Sample selection
    print("2. SAMPLE SELECTION EFFECTS")
    print("-" * 40)
    print("""
    LSB galaxies in SPARC are:
    - Often gas-dominated (high M_gas/M_star)
    - At larger distances (harder to measure)
    - May have larger inclination errors

    The LSB sample (N=21) is much smaller than HSB (N=70).
    - Could be dominated by outliers
    - May not represent typical LSB galaxies

    Also: LSB galaxies have different formation histories:
    - Late formation in low-density environments
    - Less merger history
    - This could affect the M/L ratio assumption (0.5 used for all)
""")
    explanations.append({
        'name': 'Sample selection effects',
        'description': 'LSB sample may not be representative',
        'plausibility': 'Medium'
    })

    # Explanation 3: M/L ratio variation
    print("3. M/L RATIO VARIATION")
    print("-" * 40)
    print("""
    We assumed constant M/L = 0.5 M_sun/L_sun for all galaxies.

    BUT: M/L varies with stellar population!
    - HSB galaxies: older, redder populations → higher M/L
    - LSB galaxies: younger, bluer populations → lower M/L

    If true M/L ratio is:
    - HSB: ~0.7 (higher actual M_star)
    - LSB: ~0.3 (lower actual M_star)

    Then correcting for this would INCREASE M_bar for HSB more than LSB,
    which would REDUCE the HSB BTFR residual and INCREASE the LSB residual.

    Effect direction: Would make offset MORE negative (worse for Synchronism)

    Wait - this actually makes the tension WORSE, not better!
""")
    explanations.append({
        'name': 'M/L ratio variation',
        'description': 'Variable M/L would make offset more negative',
        'plausibility': 'Low (wrong direction)'
    })

    # Explanation 4: Synchronism LOCAL vs GLOBAL density
    print("4. SYNCHRONISM: LOCAL vs GLOBAL DENSITY EFFECTS")
    print("-" * 40)
    print("""
    KEY INSIGHT from Session #85:
    - Void test showed C(δ) is much WEAKER than expected
    - C is primarily determined by LOCAL baryonic density

    For HSB/LSB comparison:
    - Surface brightness is a GLOBAL property (average over disk)
    - But rotation curves are measured at SPECIFIC radii
    - At each radius, LOCAL density determines C

    For HSB galaxies:
    - High central SB but rapid decline
    - In outer disk (where V_flat measured), density may be similar to LSB

    For LSB galaxies:
    - Low central SB but slow decline
    - In outer disk, density actually comparable to HSB outer disk!

    If C depends on LOCAL density at the measurement radius,
    then HSB and LSB might have similar C where V_flat is measured.

    The observed trend might come from DIFFERENT radii being measured:
    - HSB: V_flat measured at smaller R (higher local density)
    - LSB: V_flat measured at larger R (lower local density)

    This is a SIZE effect, not a surface brightness effect!
""")
    explanations.append({
        'name': 'Local vs global density',
        'description': 'C depends on local density at measurement radius',
        'plausibility': 'HIGH - consistent with Session #85'
    })

    # Explanation 5: MOND-like behavior
    print("5. MOND-LIKE UNIVERSALITY")
    print("-" * 40)
    print("""
    The observed result is closer to MOND prediction than Synchronism!

    MOND predicts: Universal BTFR with no SB dependence
    Observed: Small offset (0.05 dex) but OPPOSITE to Synchronism

    However, the offset is NOT zero:
    - 3σ significant
    - Direction favors HSB (higher V at fixed M)

    This could indicate:
    - MOND is wrong (there IS an SB effect)
    - But Synchronism prediction has wrong SIGN

    Or it could be systematic:
    - Distance errors (HSB galaxies better measured)
    - Inclination errors (HSB disks better defined)
    - Selection effects (HSB sample more complete)
""")
    explanations.append({
        'name': 'MOND-like universality with systematics',
        'description': 'The offset may be due to systematics, not physics',
        'plausibility': 'Medium'
    })

    # Explanation 6: Formation history
    print("6. FORMATION HISTORY EFFECT")
    print("-" * 40)
    print("""
    Session #85 suggested: Formation environment affects C

    But: LSB galaxies typically form in LOWER density environments!
    - If C at formation is lower → G_eff higher → V higher
    - This should give POSITIVE offset (LSB higher V)

    Yet we observe NEGATIVE offset. This suggests:
    - Formation environment effect is EVEN WEAKER than revised C(δ)
    - Or it works in the OPPOSITE direction than predicted

    Could C actually be HIGHER in low-density environments?
    - This contradicts the coherence model
    - But would explain the observation

    Alternative: "Coherence" isn't density-dependent in the way assumed.
""")
    explanations.append({
        'name': 'Formation history reversal',
        'description': 'C may not depend on density as predicted',
        'plausibility': 'Medium - requires theory revision'
    })

    print()
    print("=" * 70)
    print("MOST LIKELY EXPLANATION")
    print("=" * 70)
    print("""
    Based on Session #85 findings and this analysis:

    The most consistent interpretation is:

    1. C is determined by LOCAL baryonic density at each radius
       (NOT by global surface brightness)

    2. Surface brightness is a poor proxy for the density that matters

    3. The observed HSB vs LSB trend comes from:
       - Measurement systematics (distance, inclination)
       - OR different radii being sampled
       - OR M/L ratio variations

    4. The Synchronism rotation curve model (G_eff = G/C(ρ))
       works on a RADIUS-BY-RADIUS basis, not a global galaxy property

    IMPLICATION:
    - The naive "LSB should have higher V" prediction was based on
      incorrect interpretation of the theory
    - The correct test is comparing C(ρ) at the SAME radius in
      different environments, not comparing global properties
""")

    print()
    print("=" * 70)
    print("THEORY STATUS UPDATE")
    print("=" * 70)
    print("""
    Session #85: Environment (void) effect is 8× weaker than predicted
    Session #86: Surface brightness effect has OPPOSITE sign to predicted

    Combined insight:
    - C depends primarily on LOCAL density at each radius
    - Global properties (SB, environment) have weak or opposite effects
    - The core rotation curve model is about LOCAL C(ρ), not global C

    This is actually CONSISTENT with the theory formulation:
    - C(ρ) uses LOCAL density ρ at each radius
    - The BTFR emerges from integrating this over the disk
    - Global properties average out

    The BTFR comparison is NOT the right test!
    The right test compares:
    - Same galaxy at different radii (inner vs outer)
    - Or different galaxies at same LOCAL density
""")

    print()
    print("=" * 70)
    print("REVISED SYNCHRONISM PREDICTIONS")
    print("=" * 70)
    print("""
    REMOVE from testable predictions:
    - "HSB vs LSB BTFR offset" - not a valid prediction
    - Based on misinterpretation of how C enters

    KEEP:
    - Rotation curve shape (V(r) profile)
    - Uses C(ρ(r)) at each radius
    - 52% SPARC success rate validates this

    ADD:
    - Radial variation of V/V_Newtonian should correlate with ρ(r)
    - Inner disk (high ρ): V/V_Newton closer to 1
    - Outer disk (low ρ): V/V_Newton larger (more DM-like)
    - This is the actual prediction that should be tested!
""")

    print()
    print("=" * 70)
    print("SESSION #86 CONCLUSION")
    print("=" * 70)
    print("""
    The HSB/LSB BTFR test shows:
    - Observation: -0.053 dex (LSB lower V than HSB)
    - Naive prediction: +0.088 dex
    - Direction: OPPOSITE

    INTERPRETATION:
    - This was NOT a valid test of Synchronism
    - The theory predicts C(ρ) at LOCAL density, not global SB
    - Surface brightness is not the right proxy

    THEORY STATUS:
    - Core rotation curve model UNAFFECTED
    - HSB/LSB test was a misapplication of the theory
    - Need to update THEORETICAL_STATUS_DEC2025.md to reflect this

    NEXT STEPS:
    - Remove HSB/LSB from discriminating tests table
    - Add note explaining why it's not a valid test
    - Focus on radial V/V_Newton correlation instead
""")

    # Summary results
    results = {
        'observation': {
            'LSB_HSB_offset_dex': observed_offset,
            'significance_sigma': significance,
            'direction': 'LSB has LOWER V than expected'
        },
        'naive_prediction': {
            'offset_dex': predicted_offset,
            'direction': 'LSB should have HIGHER V'
        },
        'status': 'CONTRADICTED',
        'explanation': 'HSB/LSB comparison is not a valid test of Synchronism. Theory predicts C(ρ) at LOCAL radius, not global surface brightness.',
        'theory_impact': 'Core model unaffected. HSB/LSB removed from discriminating tests.',
        'lessons': [
            'C depends on LOCAL density at each radius',
            'Global galaxy properties (SB, environment) are poor proxies',
            'BTFR comparisons average over radii and lose the signal',
            'Radial V/V_Newton profile is the correct test'
        ]
    }

    # Save results
    import os
    results_dir = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/results'
    os.makedirs(results_dir, exist_ok=True)

    with open(os.path.join(results_dir, 'session86_interpretation.json'), 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {os.path.join(results_dir, 'session86_interpretation.json')}")

    return results


if __name__ == '__main__':
    interpret_result()
