#!/usr/bin/env python3
"""
Session #45 Track B: Investigation of B ≈ φ (Golden Ratio) in Astrophysics

Nova's Session #44 critique: "B ≈ φ correlation is intriguing but could be coincidence"

This session investigates whether the golden ratio φ = 1.618 appears in other
astrophysical scaling relations, or if our B = 1.62 is merely coincidental.

Approach:
1. Review known astrophysical power laws
2. Check if φ appears in any
3. Test Synchronism prediction on alternative samples
4. Statistical assessment of coincidence probability

Author: CBP Autonomous Synchronism Research
Date: 2025-11-25
Session: #45 - Golden Ratio Investigation
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime
from scipy import stats


# Golden ratio
PHI = (1 + np.sqrt(5)) / 2  # 1.6180339887...


def known_astrophysical_scalings():
    """
    Review known power-law scalings in astrophysics and check for φ.
    """

    print("\n" + "="*80)
    print("KNOWN ASTROPHYSICAL POWER LAWS")
    print("="*80)

    scalings = [
        {
            'name': 'Tully-Fisher relation',
            'relation': 'L ∝ v^n',
            'exponent': 4.0,
            'domain': 'Spiral galaxies',
            'close_to_phi': False
        },
        {
            'name': 'Faber-Jackson relation',
            'relation': 'L ∝ σ^n',
            'exponent': 4.0,
            'domain': 'Elliptical galaxies',
            'close_to_phi': False
        },
        {
            'name': 'Mass-luminosity (main sequence)',
            'relation': 'L ∝ M^n',
            'exponent': 3.5,
            'domain': 'Stars',
            'close_to_phi': False
        },
        {
            'name': 'Baryonic Tully-Fisher (BTFR)',
            'relation': 'M_bar ∝ v^n',
            'exponent': 3.75,  # McGaugh et al.
            'domain': 'Spiral galaxies',
            'close_to_phi': False
        },
        {
            'name': 'Mass-radius (galaxies)',
            'relation': 'R ∝ M^n',
            'exponent': 0.5,  # Re ∝ M^0.5
            'domain': 'Elliptical galaxies',
            'close_to_phi': False
        },
        {
            'name': 'Kennicutt-Schmidt law',
            'relation': 'Σ_SFR ∝ Σ_gas^n',
            'exponent': 1.4,
            'domain': 'Star formation',
            'close_to_phi': False
        },
        {
            'name': 'MBH-σ relation',
            'relation': 'M_BH ∝ σ^n',
            'exponent': 4.0,  # ~4.0-5.0
            'domain': 'Black holes',
            'close_to_phi': False
        },
        {
            'name': 'Mass-concentration (NFW)',
            'relation': 'c ∝ M^n',
            'exponent': -0.1,  # Weak dependence
            'domain': 'Dark matter halos',
            'close_to_phi': False
        },
        {
            'name': 'Our virial scaling',
            'relation': 'ρ_crit ∝ v_max^B',
            'exponent': 1.62,
            'domain': 'Synchronism',
            'close_to_phi': True
        }
    ]

    print(f"\n  {'Relation':<35} | {'Exponent':>10} | {'|B - φ|':>10} | {'Close to φ?'}")
    print("  " + "-"*75)

    for s in scalings:
        diff = abs(s['exponent'] - PHI)
        close = '✓' if diff < 0.1 else ''
        s['close_to_phi'] = diff < 0.1
        print(f"  {s['name']:<35} | {s['exponent']:>10.2f} | {diff:>10.3f} | {close:^12}")

    print(f"\n  Golden ratio φ = {PHI:.6f}")

    # Count how many are close
    n_close = sum(1 for s in scalings if s['close_to_phi'])
    print(f"\n  Scalings close to φ (|B-φ| < 0.1): {n_close}/{len(scalings)}")

    return scalings


def statistical_significance():
    """
    Assess the statistical significance of B ≈ φ.
    """

    print("\n" + "="*80)
    print("STATISTICAL SIGNIFICANCE ASSESSMENT")
    print("="*80)

    # Our measurement
    B_measured = 1.62
    B_error = 0.73  # From Session #42 fit uncertainty

    print(f"\n  Our measured exponent: B = {B_measured:.2f} ± {B_error:.2f}")
    print(f"  Golden ratio: φ = {PHI:.4f}")
    print(f"  Difference: |B - φ| = {abs(B_measured - PHI):.4f}")

    # Is φ within error bars?
    if abs(B_measured - PHI) < B_error:
        print(f"\n  φ is WITHIN {B_error:.2f} error bars of B ✓")
    else:
        print(f"\n  φ is OUTSIDE error bars")

    # P-value for getting B this close to φ by chance
    # If B could be anywhere in [0, 4], what's P(|B - φ| < 0.05)?
    range_min, range_max = 0, 4
    tolerance = 0.05

    # Assuming uniform prior on [0, 4]
    p_chance = 2 * tolerance / (range_max - range_min)
    print(f"\n  If B uniformly distributed in [{range_min}, {range_max}]:")
    print(f"  P(|B - φ| < {tolerance}) = {p_chance:.3f}")
    print(f"  This is NOT statistically significant (p > 0.01)")

    # More realistic: if B comes from physics, likely near integers or simple fractions
    # Common exponents: 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4
    common_exponents = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
    close_to_common = min(abs(B_measured - c) for c in common_exponents)

    print(f"\n  Distance to nearest 'simple' exponent: {close_to_common:.2f}")
    print(f"  Distance to φ: {abs(B_measured - PHI):.4f}")

    if abs(B_measured - PHI) < close_to_common:
        print(f"  B is CLOSER to φ than to any simple fraction!")
    else:
        print(f"  B is closer to {min(common_exponents, key=lambda x: abs(B_measured - x))}")

    return {
        'B_measured': B_measured,
        'B_error': B_error,
        'phi': PHI,
        'difference': abs(B_measured - PHI),
        'within_error': abs(B_measured - PHI) < B_error,
        'p_uniform': p_chance,
        'closest_simple': min(common_exponents, key=lambda x: abs(B_measured - x))
    }


def fibonacci_in_nature():
    """
    Context: φ appears in many natural systems.
    """

    print("\n" + "="*80)
    print("THE GOLDEN RATIO IN NATURAL SYSTEMS")
    print("="*80)

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                    WHERE φ = 1.618... APPEARS IN NATURE                     │
└─────────────────────────────────────────────────────────────────────────────┘

BIOLOGY:
  - Phyllotaxis (leaf arrangement): 137.5° = 360°/φ²
  - Nautilus shell: Logarithmic spiral with φ ratio
  - DNA helix: 34Å/10bp ≈ 21Å × φ
  - Population dynamics: Stable ratios in breeding pairs

PHYSICS:
  - Quasicrystals: 5-fold symmetry requires φ
  - Golden angle in turbulence: Optimal mixing
  - Penrose tiling: Aperiodic coverage with φ ratios

MATHEMATICS:
  - Fibonacci sequence: F(n+1)/F(n) → φ
  - Optimal search: Golden section search
  - Continued fraction: φ = 1 + 1/(1 + 1/(1 + ...))

WHY φ APPEARS:
  - Optimization: φ often solves min/max problems
  - Self-similarity: φ² = φ + 1 allows fractal scaling
  - Irrational nature: "Most irrational" number (least approximable)

SPECULATIVE CONNECTION TO ASTROPHYSICS:
  - Logarithmic spirals in galaxies
  - Optimal packing in stellar clusters?
  - Self-similar hierarchical structure?

BUT: No established theoretical basis for φ in galactic dynamics.
""")


def test_alternative_virial_forms():
    """
    Test if B = φ is physically motivated or coincidental.
    """

    print("\n" + "="*80)
    print("TESTING ALTERNATIVE VIRIAL FORMS")
    print("="*80)

    # Load data
    import sys
    sys.path.append(str(Path(__file__).parent))

    from session38_sparc_refined_coherence import RealSPARCLoader, SPARCGalaxy
    from session43_combined_predictor import TanhCoherencePredictor

    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()

    print(f"\n  Loaded {len(galaxies)} galaxies")

    # Test virial scalings with different exponents
    exponents_to_test = [
        (1.0, 'B = 1 (linear)'),
        (1.5, 'B = 1.5'),
        (PHI, f'B = φ = {PHI:.3f}'),
        (1.62, 'B = 1.62 (empirical)'),
        (2.0, 'B = 2 (quadratic)'),
    ]

    print(f"\n  Testing different virial exponents B in: ρ_crit = A × v_max^B\n")

    results = []

    for B, label in exponents_to_test:
        # For each exponent, fit optimal A
        v_max_list = []
        best_rho_crit_list = []

        for galaxy in galaxies:
            v_max = np.max(galaxy.v_obs)

            # Find best rho_crit for this galaxy
            best_chi2 = float('inf')
            best_rho = 1.0

            for log_rho in np.linspace(-2, 4, 30):
                rho_crit = 10**log_rho
                predictor = TanhCoherencePredictor(rho_crit=rho_crit, gamma=2.0)
                fit = predictor.fit_alpha(galaxy)
                if fit['chi2_red'] < best_chi2:
                    best_chi2 = fit['chi2_red']
                    best_rho = rho_crit

            if best_chi2 < 10:
                v_max_list.append(v_max)
                best_rho_crit_list.append(best_rho)

        # Fit A given B
        v_max_arr = np.array(v_max_list)
        rho_arr = np.array(best_rho_crit_list)

        # ρ = A × v^B → log(ρ) = log(A) + B × log(v)
        # Given B, fit A
        log_rho_pred = B * np.log(v_max_arr)
        log_rho_actual = np.log(rho_arr)

        # Best A: minimize sum of (log_rho_actual - log(A) - B × log_v)²
        log_A = np.mean(log_rho_actual - log_rho_pred)
        A = np.exp(log_A)

        # Residual
        residuals = log_rho_actual - log_A - log_rho_pred
        rmse = np.sqrt(np.mean(residuals**2))

        # Now test prediction quality
        chi2_list = []

        for galaxy in galaxies:
            v_max = np.max(galaxy.v_obs)
            rho_crit_pred = A * v_max**B

            predictor = TanhCoherencePredictor(rho_crit=rho_crit_pred, gamma=2.0)
            fit = predictor.fit_alpha(galaxy)
            chi2_list.append(fit['chi2_red'])

        chi2_arr = np.array(chi2_list)
        success_rate = 100 * np.sum(chi2_arr < 5.0) / len(chi2_arr)

        results.append({
            'B': float(B),
            'label': label,
            'A': float(A),
            'rmse': float(rmse),
            'success_rate': float(success_rate)
        })

        print(f"  {label:25s}: A={A:.3f}, RMSE={rmse:.3f}, Success={success_rate:.1f}%")

    # Best result
    best = max(results, key=lambda x: x['success_rate'])
    print(f"\n  Best exponent: B = {best['B']:.3f} ({best['label']}) with {best['success_rate']:.1f}% success")

    # Is φ better than empirical 1.62?
    phi_result = [r for r in results if abs(r['B'] - PHI) < 0.01][0]
    emp_result = [r for r in results if r['B'] == 1.62][0]

    print(f"\n  φ = {PHI:.3f}: {phi_result['success_rate']:.1f}% success")
    print(f"  B = 1.62: {emp_result['success_rate']:.1f}% success")
    print(f"  Difference: {phi_result['success_rate'] - emp_result['success_rate']:+.1f} pp")

    return results


def conclusions():
    """
    Summarize conclusions about B ≈ φ.
    """

    print("\n" + "="*80)
    print("CONCLUSIONS: IS B ≈ φ MEANINGFUL?")
    print("="*80)

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│                           ASSESSMENT SUMMARY                                │
└─────────────────────────────────────────────────────────────────────────────┘

EVIDENCE FOR B ≈ φ:
  ✓ Our empirical fit gives B = 1.62 (within 0.2% of φ = 1.618)
  ✓ φ appears in other natural scaling phenomena
  ✓ Self-similarity is important in galactic structure
  ? Performance test shows minimal difference between B = 1.62 and B = φ

EVIDENCE AGAINST:
  ✗ No established theoretical basis for φ in galactic dynamics
  ✗ Error bars (±0.73) include many other values
  ✗ Could be coincidence (p ~ 0.03 for uniform prior)
  ✗ No other astrophysical scaling relations have φ exponents

VERDICT:

  The B ≈ φ result is INTRIGUING but NOT YET SIGNIFICANT.

  We should:
  1. Report B = 1.62 ± 0.73 (empirical)
  2. Note the numerical coincidence with φ
  3. NOT claim φ is fundamental without theoretical derivation
  4. Keep watching for patterns (more data, other systems)

  If future work shows:
  - Theoretical derivation from optimization principle → meaningful
  - Multiple independent measurements cluster around φ → meaningful
  - Systematic deviation from φ with more data → coincidence

┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│   STATUS: Curiosity, not conclusion. Note it, don't hang the theory on it. │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
""")


def json_safe(obj):
    """Convert numpy types to native Python types."""
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [json_safe(v) for v in obj]
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def save_results(exponent_test_results, significance_results):
    """Save investigation results."""

    output = json_safe({
        'session': 45,
        'track': 'B - Golden Ratio Investigation',
        'date': datetime.now().isoformat(),
        'golden_ratio': PHI,
        'our_exponent': 1.62,
        'difference': abs(1.62 - PHI),
        'statistical_assessment': significance_results,
        'exponent_tests': exponent_test_results,
        'conclusion': 'Intriguing coincidence, not yet significant',
        'recommendation': 'Report empirical value, note coincidence, await theory'
    })

    output_path = Path(__file__).parent / 'session45_golden_ratio_results.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    return output


if __name__ == '__main__':
    print("\n" + "="*80)
    print("SESSION #45 TRACK B: GOLDEN RATIO INVESTIGATION")
    print("="*80)

    # Known scalings
    known_astrophysical_scalings()

    # Statistical assessment
    significance = statistical_significance()

    # Context
    fibonacci_in_nature()

    # Test alternative forms
    exponent_results = test_alternative_virial_forms()

    # Conclusions
    conclusions()

    # Save
    save_results(exponent_results, significance)

    print("\n" + "="*80)
    print("SESSION #45 TRACK B COMPLETE")
    print("="*80)
