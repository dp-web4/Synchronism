#!/usr/bin/env python3
"""
Session #79 Track A: Validate BTFR-Derived Parameters on SPARC

From Session #78:
- B = 4 - 3δ where δ ≈ 0.79 → B ≈ 1.63
- A = 3 A_TF / (4π R_0³) where A_TF = 47 M_sun/(km/s)^4

This script:
1. Tests the BTFR-derived B on SPARC rotation curves
2. Finds the optimal A (normalization) for the new B
3. Compares to empirical A=0.25, B=1.62

Author: CBP Autonomous Synchronism Research
Date: December 3, 2025
Session: #79 - BTFR Validation
"""

import numpy as np
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy import stats
from scipy.optimize import minimize_scalar

# Import utilities
sys.path.append(str(Path(__file__).parent))
from synchronism_real_sparc_validation import RealSPARCLoader
from synchronism_sparc_validation import SPARCGalaxy, SynchronismPredictor


class BTFRCoherencePredictor(SynchronismPredictor):
    """
    Synchronism predictor using BTFR-derived ρ_crit formula.

    ρ_crit = A × V^B where B = 4 - 3δ (from BTFR)
    """

    def __init__(self, A: float, B: float, gamma: float = 2.0, beta: float = 0.30):
        super().__init__(gamma=gamma, beta=beta)
        self.A = A
        self.B = B
        self.tanh_gamma = gamma

    def predict_rho_crit(self, v_max: float) -> float:
        """Predict critical density from v_max using BTFR scaling."""
        return self.A * (v_max ** self.B)

    def compute_coherence(self, rho_vis: np.ndarray, rho_crit: float) -> np.ndarray:
        """Compute coherence using tanh formula."""
        arg = self.tanh_gamma * np.log(rho_vis / rho_crit + 1.0)
        C = np.tanh(arg)
        return np.clip(C, 0.0, 0.99999)

    def fit_galaxy(self, galaxy: SPARCGalaxy,
                   alpha_range=np.logspace(-2, 2, 100)) -> dict:
        """
        Fit alpha parameter for galaxy using BTFR-predicted rho_crit.
        """
        v_max = np.max(galaxy.v_obs)
        rho_crit = self.predict_rho_crit(v_max)

        best_chi2 = np.inf
        best_alpha = 1.0

        for alpha in alpha_range:
            rho_vis = galaxy.total_baryonic_density()
            C_vis = self.compute_coherence(rho_vis, rho_crit)
            rho_DM = self.compute_dark_matter_density(rho_vis, C_vis, alpha)
            rho_total = rho_vis + rho_DM
            M_enc = self.compute_enclosed_mass(galaxy.radius, rho_total)
            v_pred = np.sqrt(self.G * M_enc / np.maximum(galaxy.radius, 1e-6))

            residuals = (galaxy.v_obs - v_pred) ** 2
            weights = 1.0 / (galaxy.v_err ** 2 + 1e-6)
            chi2 = np.sum(residuals * weights)
            dof = len(galaxy.v_obs) - 1
            chi2_red = chi2 / max(dof, 1)

            if chi2_red < best_chi2:
                best_chi2 = chi2_red
                best_alpha = alpha

        return {
            'name': galaxy.name,
            'v_max': float(v_max),
            'rho_crit': float(rho_crit),
            'alpha_best': float(best_alpha),
            'chi2_red': float(best_chi2)
        }


def get_sparc_data_dir():
    """Get absolute path to SPARC data directory."""
    script_dir = Path(__file__).parent.absolute()
    return str(script_dir / "sparc_real_data" / "galaxies")


def test_btfr_derived_parameters():
    """
    Test the BTFR-derived B parameter on SPARC galaxies.
    """
    print("=" * 70)
    print("Session #79 Track A: BTFR-Derived Parameter Validation")
    print("=" * 70)

    # Load SPARC galaxies with absolute path
    data_dir = get_sparc_data_dir()
    loader = RealSPARCLoader(data_dir=data_dir)
    galaxies = loader.load_all_galaxies()

    print(f"\nLoaded {len(galaxies)} galaxies from SPARC")

    # BTFR-derived B
    delta = 0.79  # From size-velocity scaling
    B_btfr = 4 - 3 * delta

    print(f"\nBTFR-Derived Parameters:")
    print(f"  δ = {delta:.2f} (from R ∝ V^δ)")
    print(f"  B = 4 - 3δ = {B_btfr:.3f}")
    print(f"  (Empirical B = 1.62 for comparison)")

    # Test range of A values with BTFR-derived B
    A_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

    results = {}

    for A in A_values:
        print(f"\n{'='*60}")
        print(f"Testing A = {A:.2f}, B = {B_btfr:.3f}")
        print(f"{'='*60}")

        predictor = BTFRCoherencePredictor(A=A, B=B_btfr, gamma=2.0, beta=0.30)

        galaxy_results = []
        for i, galaxy in enumerate(galaxies):
            try:
                result = predictor.fit_galaxy(galaxy)
                galaxy_results.append(result)

                if (i + 1) % 50 == 0:
                    print(f"  Processed {i+1}/{len(galaxies)} galaxies...")
            except Exception as e:
                print(f"  Error on {galaxy.name}: {e}")

        # Analyze results
        chi2_values = np.array([r['chi2_red'] for r in galaxy_results])
        success_rate = 100.0 * np.sum(chi2_values < 5.0) / len(chi2_values)

        results[A] = {
            'A': A,
            'B': B_btfr,
            'success_rate': float(success_rate),
            'median_chi2': float(np.median(chi2_values)),
            'mean_chi2': float(np.mean(chi2_values)),
            'n_success': int(np.sum(chi2_values < 5.0)),
            'n_total': len(galaxy_results)
        }

        print(f"\n  Success rate (χ² < 5): {success_rate:.1f}%")
        print(f"  Median χ²: {np.median(chi2_values):.2f}")

    # Find best A
    best_A = max(results.keys(), key=lambda a: results[a]['success_rate'])

    print("\n" + "=" * 70)
    print("SUMMARY: BTFR-Derived B with Varying A")
    print("=" * 70)

    print(f"\n{'A':>8} {'B':>8} {'Success %':>12} {'Median χ²':>12}")
    print("-" * 44)
    for A in sorted(results.keys()):
        r = results[A]
        marker = " ← BEST" if A == best_A else ""
        print(f"{r['A']:>8.2f} {r['B']:>8.3f} {r['success_rate']:>12.1f} {r['median_chi2']:>12.2f}{marker}")

    print(f"\nBest configuration: A = {best_A:.2f}, B = {B_btfr:.3f}")
    print(f"Empirical reference: A = 0.25, B = 1.62")

    return results, best_A, B_btfr


def compare_to_empirical():
    """
    Compare BTFR-derived to empirical parameters.
    """
    print("\n" + "=" * 70)
    print("COMPARISON: BTFR-Derived vs Empirical")
    print("=" * 70)

    data_dir = get_sparc_data_dir()
    loader = RealSPARCLoader(data_dir=data_dir)
    galaxies = loader.load_all_galaxies()

    # BTFR-derived
    delta = 0.79
    B_btfr = 4 - 3 * delta
    A_btfr = 0.25  # Use same A for fair comparison

    # Empirical
    A_emp = 0.25
    B_emp = 1.62

    print(f"\nBTFR-Derived: A = {A_btfr}, B = {B_btfr:.3f}")
    print(f"Empirical:    A = {A_emp}, B = {B_emp}")

    # Test both
    predictor_btfr = BTFRCoherencePredictor(A=A_btfr, B=B_btfr, gamma=2.0, beta=0.30)
    predictor_emp = BTFRCoherencePredictor(A=A_emp, B=B_emp, gamma=2.0, beta=0.30)

    results_btfr = []
    results_emp = []

    print("\nProcessing galaxies...")
    for i, galaxy in enumerate(galaxies):
        try:
            r_btfr = predictor_btfr.fit_galaxy(galaxy)
            r_emp = predictor_emp.fit_galaxy(galaxy)
            results_btfr.append(r_btfr)
            results_emp.append(r_emp)
        except Exception as e:
            pass

        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(galaxies)}...")

    chi2_btfr = np.array([r['chi2_red'] for r in results_btfr])
    chi2_emp = np.array([r['chi2_red'] for r in results_emp])

    success_btfr = 100.0 * np.sum(chi2_btfr < 5.0) / len(chi2_btfr)
    success_emp = 100.0 * np.sum(chi2_emp < 5.0) / len(chi2_emp)

    print("\n" + "-" * 50)
    print("RESULTS:")
    print("-" * 50)
    print(f"\n{'Model':<20} {'Success %':>12} {'Median χ²':>12}")
    print("-" * 44)
    print(f"{'BTFR-Derived':<20} {success_btfr:>12.1f} {np.median(chi2_btfr):>12.2f}")
    print(f"{'Empirical':<20} {success_emp:>12.1f} {np.median(chi2_emp):>12.2f}")

    print(f"\nDifference: {success_btfr - success_emp:+.1f} percentage points")

    if abs(success_btfr - success_emp) < 2.0:
        print("\n✅ BTFR-derived B performs EQUIVALENTLY to empirical B!")
        print("   This validates the theoretical derivation.")
    elif success_btfr > success_emp:
        print("\n✅ BTFR-derived B performs BETTER than empirical B!")
    else:
        print(f"\n⚠️ BTFR-derived B is {success_emp - success_btfr:.1f}pp worse than empirical.")

    return {
        'btfr': {
            'A': A_btfr, 'B': B_btfr,
            'success_rate': float(success_btfr),
            'median_chi2': float(np.median(chi2_btfr))
        },
        'empirical': {
            'A': A_emp, 'B': B_emp,
            'success_rate': float(success_emp),
            'median_chi2': float(np.median(chi2_emp))
        }
    }


def explore_optimal_delta():
    """
    Find the optimal δ that maximizes SPARC success.
    """
    print("\n" + "=" * 70)
    print("EXPLORATION: Finding Optimal δ (R ∝ V^δ)")
    print("=" * 70)

    data_dir = get_sparc_data_dir()
    loader = RealSPARCLoader(data_dir=data_dir)
    galaxies = loader.load_all_galaxies()

    # Test range of delta values
    delta_values = np.linspace(0.5, 1.0, 11)
    A_fixed = 0.25  # Keep A fixed to isolate δ effect

    results = {}

    for delta in delta_values:
        B = 4 - 3 * delta

        predictor = BTFRCoherencePredictor(A=A_fixed, B=B, gamma=2.0, beta=0.30)

        galaxy_results = []
        for galaxy in galaxies:
            try:
                result = predictor.fit_galaxy(galaxy)
                galaxy_results.append(result)
            except:
                pass

        chi2_values = np.array([r['chi2_red'] for r in galaxy_results])
        success_rate = 100.0 * np.sum(chi2_values < 5.0) / len(chi2_values)

        results[delta] = {
            'delta': float(delta),
            'B': float(B),
            'success_rate': float(success_rate),
            'median_chi2': float(np.median(chi2_values))
        }

    # Find optimal
    best_delta = max(results.keys(), key=lambda d: results[d]['success_rate'])

    print(f"\n{'δ':>8} {'B = 4-3δ':>12} {'Success %':>12} {'Median χ²':>12}")
    print("-" * 48)
    for delta in sorted(results.keys()):
        r = results[delta]
        marker = " ← BEST" if delta == best_delta else ""
        marker = " ← THEORY" if abs(delta - 0.79) < 0.02 else marker
        marker = " ← EMP" if abs(r['B'] - 1.62) < 0.02 else marker
        print(f"{r['delta']:>8.2f} {r['B']:>12.3f} {r['success_rate']:>12.1f} {r['median_chi2']:>12.2f}{marker}")

    print(f"\nOptimal δ = {best_delta:.2f} → B = {4 - 3*best_delta:.3f}")
    print(f"Theoretical δ = 0.79 → B = {4 - 3*0.79:.3f}")
    print(f"For B = 1.62 (empirical), need δ = {(4 - 1.62)/3:.3f}")

    return results, best_delta


if __name__ == '__main__':
    print("=" * 70)
    print("SESSION #79: BTFR-DERIVED PARAMETER VALIDATION")
    print("=" * 70)

    # Run validation tests
    results_A, best_A, B_btfr = test_btfr_derived_parameters()

    # Compare to empirical
    comparison = compare_to_empirical()

    # Explore optimal delta
    delta_results, best_delta = explore_optimal_delta()

    # Save all results
    output = {
        'session': 79,
        'track': 'A',
        'title': 'BTFR-Derived Parameter Validation',
        'date': datetime.now().isoformat(),

        'btfr_derivation': {
            'delta_theory': 0.79,
            'B_theory': 4 - 3*0.79,
            'B_empirical': 1.62
        },

        'A_sweep_results': {str(k): v for k, v in results_A.items()},
        'best_A': best_A,

        'comparison': comparison,

        'delta_exploration': {str(k): v for k, v in delta_results.items()},
        'best_delta': float(best_delta)
    }

    output_path = Path(__file__).parent / 'results' / 'session79_btfr_validation.json'
    output_path.parent.mkdir(exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n\nResults saved to: {output_path}")

    print("\n" + "=" * 70)
    print("SESSION #79 TRACK A COMPLETE")
    print("=" * 70)
