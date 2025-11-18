#!/usr/bin/env python3
"""
Session #24: ρ_central Calculation Bias Test
==============================================

Tests whether NGC underprediction (Session #22) is due to systematic bias
in how we calculate ρ_central from galaxy data.

Hypothesis:
  NGC spirals have bulges that dominate central density.
  Current method (first data point) may systematically OVERESTIMATE ρ_central.
  This would make model predict LOW ρ_sat (high ρ_c → low ρ_sat in screening model).

Test:
  Compare ρ_central calculation methods:
  1. First point: galaxy.total_baryonic_density()[0] (current)
  2. Peak density: max(galaxy.total_baryonic_density())
  3. Average central 1 kpc
  4. Median central 3 kpc

  If NGC systematic difference disappears, bias confirmed!

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-17
Session: #24 - Quick diagnostic following Session #23 null result
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle
import os
import sys

# Import Session #20/22 infrastructure
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_session22_extract_data import load_session20_data
from synchronism_real_sparc_validation import RealSPARCLoader


# ============================================================================
# Part 1: Alternative ρ_central Calculation Methods
# ============================================================================

def compute_rho_central_methods(galaxy) -> dict:
    """
    Compute ρ_central using multiple methods.

    Args:
        galaxy: SPARCGalaxy object

    Returns:
        Dictionary with different ρ_central estimates
    """
    rho_baryon = galaxy.total_baryonic_density()  # Surface density [M_☉/pc²]
    radius = galaxy.radius  # [kpc]

    # Convert to volume density (approximate)
    # For thin disk: ρ_vol ~ Σ / h where h ~ 1 kpc (scale height)
    # This is crude but consistent across methods
    h_scale = 0.3  # kpc (disk scale height estimate)
    rho_vol = rho_baryon / h_scale  # [M_☉/pc³]

    methods = {}

    # Method 1: First point (current Session #20 method)
    methods['first_point'] = rho_vol[0] if len(rho_vol) > 0 else np.nan

    # Method 2: Peak density (maximum value)
    methods['peak'] = np.max(rho_vol) if len(rho_vol) > 0 else np.nan

    # Method 3: Average over central 1 kpc
    central_1kpc = radius <= 1.0
    if np.sum(central_1kpc) > 0:
        methods['avg_1kpc'] = np.mean(rho_vol[central_1kpc])
    else:
        methods['avg_1kpc'] = rho_vol[0] if len(rho_vol) > 0 else np.nan

    # Method 4: Median over central 3 kpc
    central_3kpc = radius <= 3.0
    if np.sum(central_3kpc) > 0:
        methods['median_3kpc'] = np.median(rho_vol[central_3kpc])
    else:
        methods['median_3kpc'] = rho_vol[0] if len(rho_vol) > 0 else np.nan

    # Method 5: Weighted average (weight by 1/r to emphasize center)
    if len(radius) > 0:
        weights = 1.0 / (radius + 0.1)  # Add offset to avoid div by zero
        methods['weighted'] = np.average(rho_vol, weights=weights)
    else:
        methods['weighted'] = np.nan

    return methods


# ============================================================================
# Part 2: Recompute Session #20 Data with Alternative Methods
# ============================================================================

class RhoCentralBiasTester:
    """
    Tests whether ρ_central calculation method affects NGC underprediction.
    """

    def __init__(self):
        """Initialize tester."""
        self.loader = RealSPARCLoader()

    def extract_rho_central_all_methods(self) -> dict:
        """
        Extract ρ_central using all methods for 175 SPARC galaxies.

        Returns:
            Dictionary with arrays for each method
        """
        print("="*80)
        print("SESSION #24: ρ_central CALCULATION BIAS TEST")
        print("="*80)

        print("\nLoading 175 SPARC galaxies...")
        galaxies = self.loader.load_all_galaxies(limit=None)
        print(f"✓ Loaded {len(galaxies)} galaxies")

        print("\nComputing ρ_central with multiple methods...")

        results = {
            'galaxy_names': [],
            'galaxy_types': [],
            'first_point': [],
            'peak': [],
            'avg_1kpc': [],
            'median_3kpc': [],
            'weighted': []
        }

        for i, galaxy in enumerate(galaxies):
            if (i+1) % 50 == 0:
                print(f"  [{i+1}/{len(galaxies)}] Processed...")

            # Get galaxy name and type
            name = galaxy.name
            gtype = 'NGC' if name.startswith('NGC') else \
                    'UGC' if name.startswith('UGC') else \
                    'F' if name.startswith('F') else \
                    'DDO' if name.startswith('DDO') else 'Other'

            # Compute ρ_central with all methods
            methods = compute_rho_central_methods(galaxy)

            results['galaxy_names'].append(name)
            results['galaxy_types'].append(gtype)
            results['first_point'].append(methods['first_point'])
            results['peak'].append(methods['peak'])
            results['avg_1kpc'].append(methods['avg_1kpc'])
            results['median_3kpc'].append(methods['median_3kpc'])
            results['weighted'].append(methods['weighted'])

        # Convert to arrays
        for key in ['first_point', 'peak', 'avg_1kpc', 'median_3kpc', 'weighted']:
            results[key] = np.array(results[key])

        print(f"\n✓ Computed ρ_central with 5 methods for {len(results['galaxy_names'])} galaxies")

        return results

    def compare_methods(self, rho_central_data: dict):
        """
        Compare ρ_central distributions by galaxy type.

        Args:
            rho_central_data: Output from extract_rho_central_all_methods()
        """
        print(f"\n{'='*80}")
        print("ρ_central DISTRIBUTION BY METHOD AND GALAXY TYPE")
        print(f"{'='*80}")

        methods = ['first_point', 'peak', 'avg_1kpc', 'median_3kpc', 'weighted']
        types = ['NGC', 'UGC', 'F', 'DDO', 'Other']

        for gtype in types:
            indices = [i for i, t in enumerate(rho_central_data['galaxy_types']) if t == gtype]
            if len(indices) == 0:
                continue

            print(f"\n{gtype} galaxies (n={len(indices)}):")
            print(f"  {'Method':<15} {'Median ρ_c':<20} {'Ratio to First':<15}")
            print(f"  {'-'*50}")

            first_median = np.median([rho_central_data['first_point'][i] for i in indices])

            for method in methods:
                values = [rho_central_data[method][i] for i in indices]
                median_val = np.median(values)
                ratio = median_val / first_median if first_median > 0 else np.nan

                print(f"  {method:<15} {median_val:<20.2e} {ratio:<15.3f}")

        # Test NGC vs F difference
        print(f"\n{'='*80}")
        print("NGC vs F RATIO COMPARISON")
        print(f"{'='*80}")

        ngc_indices = [i for i, t in enumerate(rho_central_data['galaxy_types']) if t == 'NGC']
        f_indices = [i for i, t in enumerate(rho_central_data['galaxy_types']) if t == 'F']

        if len(ngc_indices) > 0 and len(f_indices) > 0:
            print(f"\n{'Method':<15} {'NGC median':<15} {'F median':<15} {'NGC/F ratio':<15}")
            print(f"{'-'*60}")

            for method in methods:
                ngc_median = np.median([rho_central_data[method][i] for i in ngc_indices])
                f_median = np.median([rho_central_data[method][i] for i in f_indices])
                ratio = ngc_median / f_median if f_median > 0 else np.nan

                print(f"{method:<15} {ngc_median:<15.2e} {f_median:<15.2e} {ratio:<15.2f}")

            print(f"\nInterpretation:")
            print(f"  Session #20 used 'first_point' method")
            print(f"  If NGC/F ratio changes significantly with method → BIAS present")
            print(f"  If NGC/F ratio stable across methods → bias NOT the issue")

    def test_magnetic_screening_with_method(self, method_name: str, rho_central_data: dict):
        """
        Refit Session #22 magnetic screening model with alternative ρ_central.

        Args:
            method_name: Which ρ_central method to use
            rho_central_data: Output from extract_rho_central_all_methods()

        Returns:
            Fit results
        """
        print(f"\n{'='*80}")
        print(f"MAGNETIC SCREENING FIT WITH '{method_name.upper()}' METHOD")
        print(f"{'='*80}")

        # Load Session #20 fitted ρ_sat values
        session20_data = load_session20_data()
        rho_sats_fitted = session20_data['rho_sats_fitted']

        # Use alternative ρ_central
        rho_centrals_alt = rho_central_data[method_name]

        # Remove NaN values
        valid = ~np.isnan(rho_centrals_alt) & ~np.isnan(rho_sats_fitted)
        rho_c = rho_centrals_alt[valid]
        rho_sat = rho_sats_fitted[valid]

        print(f"\nValid galaxies: {np.sum(valid)}/{len(valid)}")

        # Magnetic screening model (Session #22)
        def rho_sat_model(rho_c, rho_sat_0, rho_mag, delta):
            return rho_sat_0 / (1.0 + (rho_c / rho_mag)**delta)

        # Fit
        try:
            p0 = [7e5, 90, 1.85]  # Session #22 values
            bounds = ([1e4, 1e1, 0.5], [1e7, 1e6, 3.0])

            popt, pcov = curve_fit(
                rho_sat_model, rho_c, rho_sat,
                p0=p0, bounds=bounds, maxfev=10000
            )

            rho_sat_0, rho_mag, delta = popt
            errors = np.sqrt(np.diag(pcov))

            # Predict
            rho_sat_pred = rho_sat_model(rho_c, rho_sat_0, rho_mag, delta)

            # Statistics
            residuals = rho_sat - rho_sat_pred
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((rho_sat - np.mean(rho_sat))**2)
            r_squared = 1 - ss_res / ss_tot

            log_residuals = np.log10(rho_sat) - np.log10(rho_sat_pred)
            log_rms = np.sqrt(np.mean(log_residuals**2))

            corr = np.corrcoef(np.log10(rho_c), np.log10(rho_sat))[0, 1]

            print(f"\n✅ FIT SUCCESSFUL")
            print(f"\nBest-fit parameters:")
            print(f"  ρ_sat,0 = {rho_sat_0:.2e} ± {errors[0]:.2e} M_☉/pc³")
            print(f"  ρ_mag   = {rho_mag:.2e} ± {errors[1]:.2e} M_☉/pc³")
            print(f"  δ       = {delta:.2f} ± {errors[2]:.2f}")

            print(f"\nFit quality:")
            print(f"  R² = {r_squared:.3f}")
            print(f"  RMS(log) = {log_rms:.3f}")
            print(f"  Correlation: r = {corr:.3f}")

            return {
                'success': True,
                'method': method_name,
                'rho_sat_0': rho_sat_0,
                'rho_mag': rho_mag,
                'delta': delta,
                'r_squared': r_squared,
                'log_rms': log_rms,
                'correlation': corr
            }

        except Exception as e:
            print(f"\n❌ FIT FAILED: {e}")
            return {'success': False, 'method': method_name}

    def comprehensive_comparison(self, rho_central_data: dict):
        """
        Test all ρ_central methods and compare results.

        Args:
            rho_central_data: Output from extract_rho_central_all_methods()
        """
        methods = ['first_point', 'peak', 'avg_1kpc', 'median_3kpc', 'weighted']

        results = []
        for method in methods:
            result = self.test_magnetic_screening_with_method(method, rho_central_data)
            results.append(result)

        # Summary table
        print(f"\n{'='*80}")
        print("SUMMARY: ρ_central METHOD COMPARISON")
        print(f"{'='*80}")

        print(f"\n{'Method':<15} {'R²':<10} {'RMS(log)':<10} {'Correlation':<12} {'Δ R² from S#22':<20}")
        print(f"{'-'*80}")

        baseline_r2 = 0.406  # Session #22 (first_point method)

        for res in results:
            if res['success']:
                delta_r2 = res['r_squared'] - baseline_r2
                print(f"{res['method']:<15} {res['r_squared']:<10.3f} {res['log_rms']:<10.3f} {res['correlation']:<12.3f} {delta_r2:+20.3f}")

        print(f"\nSession #22 baseline (first_point): R² = {baseline_r2:.3f}")

        # Interpretation
        print(f"\n{'='*80}")
        print("INTERPRETATION")
        print(f"{'='*80}")

        r2_values = [r['r_squared'] for r in results if r['success']]
        r2_range = max(r2_values) - min(r2_values)

        print(f"\nR² range across methods: {r2_range:.3f}")

        if r2_range < 0.05:
            print(f"  → MINIMAL variation (< 0.05)")
            print(f"  → ρ_central calculation method NOT the issue")
            print(f"  → NGC underprediction has PHYSICAL cause (not measurement artifact)")
        elif r2_range < 0.15:
            print(f"  → MODERATE variation (< 0.15)")
            print(f"  → Some sensitivity to ρ_central method")
            print(f"  → But unlikely to fully explain NGC underprediction")
        else:
            print(f"  → LARGE variation (> 0.15)")
            print(f"  → ρ_central calculation method MATTERS")
            print(f"  → NGC underprediction may be systematic bias!")

        # Best method
        best_idx = np.argmax([r['r_squared'] for r in results if r['success']])
        best_method = [r for r in results if r['success']][best_idx]

        print(f"\nBest method: {best_method['method']}")
        print(f"  R² = {best_method['r_squared']:.3f} (Session #22: {baseline_r2:.3f})")

        if best_method['r_squared'] > baseline_r2 + 0.05:
            print(f"  ✅ SIGNIFICANT IMPROVEMENT - Consider using this method!")
        else:
            print(f"  ⚠️ No significant improvement - Stick with original method")

        return results


# ============================================================================
# Part 3: Main Execution
# ============================================================================

def test_rho_central_bias():
    """
    Main test: Check if ρ_central calculation bias explains NGC underprediction.
    """
    print("="*80)
    print("SESSION #24: ρ_central CALCULATION BIAS TEST")
    print("Testing Nova's Session #23 recommendation")
    print("="*80)

    tester = RhoCentralBiasTester()

    # Extract ρ_central with all methods
    rho_central_data = tester.extract_rho_central_all_methods()

    # Compare methods
    tester.compare_methods(rho_central_data)

    # Test magnetic screening with each method
    results = tester.comprehensive_comparison(rho_central_data)

    # Final conclusion
    print(f"\n{'='*80}")
    print("SESSION #24 CONCLUSION")
    print(f"{'='*80}")

    r2_values = [r['r_squared'] for r in results if r['success']]
    r2_range = max(r2_values) - min(r2_values)

    print(f"\nHypothesis: NGC underprediction due to ρ_central calculation bias")

    if r2_range < 0.05:
        print(f"  Result: ❌ REJECTED")
        print(f"  Evidence: R² variation < 0.05 across methods")
        print(f"  Conclusion: NGC underprediction has PHYSICAL cause, not bias")
        print(f"\nNext steps:")
        print(f"  1. Literature B-field compilation (Session #25)")
        print(f"  2. SFR-based B-field proxy")
        print(f"  3. Accept NGC as outliers, focus on majority fit")
    else:
        print(f"  Result: ✅ SUPPORTED")
        print(f"  Evidence: R² variation {r2_range:.3f} across methods")
        print(f"  Conclusion: ρ_central method affects fit quality")
        print(f"\nNext steps:")
        print(f"  1. Adopt best method: {[r for r in results if r['success']][np.argmax(r2_values)]['method']}")
        print(f"  2. Rerun Session #20-22 analysis with corrected ρ_central")
        print(f"  3. Check if NGC underprediction resolves")

    return results


if __name__ == "__main__":
    # Run ρ_central bias test
    results = test_rho_central_bias()

    print("\n✓ Session #24 ρ_central bias test complete!")
    print("  Ready for Session #24 documentation")
