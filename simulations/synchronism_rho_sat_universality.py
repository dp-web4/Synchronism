#!/usr/bin/env python3
"""
Synchronism ρ_sat Universality Test - Session #20 Priority 2
=============================================================

Tests whether saturation density ρ_sat is universal across galaxies by:
1. Fitting ρ_sat independently for each galaxy
2. Measuring scatter σ(ρ_sat)/⟨ρ_sat⟩
3. Testing universality hypothesis: scatter < factor of 2-3

Expected outcome (Session #18 prediction):
- Universal ρ_sat ≈ 2×10^4 M_☉/pc³ ± 30% scatter
- If scatter < 50%: Universal screening mechanism ✅
- If scatter > factor of 10: Galaxy-dependent, mechanism wrong ❌

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-16
Session: #20 Priority 2 - ρ_sat universality validation
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple
import os
import sys

# Import Session #19-20 infrastructure
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_saturation_test import (
    SynchronismSaturationPredictor,
    SaturationValidator
)
from synchronism_real_sparc_validation import (
    RealSPARCLoader,
    SPARCGalaxy
)
from synchronism_galaxy_type_analysis import GalaxyTypeCategorizer


# ============================================================================
# Part 1: Per-Galaxy ρ_sat Fitting
# ============================================================================

class RhoSatUniversalityTester:
    """
    Tests ρ_sat universality by fitting independently per galaxy.

    Universality criteria:
    - Strong: σ(ρ_sat)/⟨ρ_sat⟩ < 0.3 (30% scatter)
    - Moderate: σ(ρ_sat)/⟨ρ_sat⟩ < 0.5 (factor of 2)
    - Weak: σ(ρ_sat)/⟨ρ_sat⟩ < 1.0 (factor of 3)
    - Failed: σ(ρ_sat)/⟨ρ_sat⟩ > 2.0 (factor of 10+, galaxy-dependent)
    """

    def __init__(self, coherence_model: str = 'rational'):
        """
        Initialize tester.

        Args:
            coherence_model: 'rational', 'exponential', or 'logarithmic'
        """
        self.coherence_model = coherence_model
        self.categorizer = GalaxyTypeCategorizer()

    def fit_galaxy_with_rho_sat(self, galaxy: SPARCGalaxy) -> Dict:
        """
        Fit galaxy including ρ_sat as free parameter.

        Args:
            galaxy: SPARC galaxy data

        Returns:
            Fit results with ρ_sat, C_max, α, χ²_red
        """
        from scipy.optimize import differential_evolution

        # Create predictor
        predictor = SynchronismSaturationPredictor(coherence_model=self.coherence_model)

        # Define objective function
        def objective(params):
            """Chi-squared for ρ_sat fitting."""
            if self.coherence_model == 'rational':
                # Fit: (alpha, rho_0, rho_sat)
                # C_max derived from rho_sat/rho_0
                alpha, rho_0, rho_sat = params
                C_max = None  # Will be derived in predictor
            else:
                # Fit: (alpha, rho_0, rho_sat, C_max)
                alpha, rho_0, rho_sat, C_max = params

            # Predict rotation curve
            try:
                v_pred, _, _ = predictor.predict_rotation_curve(
                    galaxy, alpha, rho_0, rho_sat, C_max
                )

                # Chi-squared
                residuals = (galaxy.v_obs - v_pred) / galaxy.v_err
                chi2 = np.sum(residuals**2)

                return chi2
            except:
                return np.inf

        # Set bounds
        if self.coherence_model == 'rational':
            # bounds: alpha [1, 100], rho_0 [0.1, 100], rho_sat [1e2, 1e6]
            bounds = [(1, 100), (0.1, 100), (1e2, 1e6)]
        else:
            # bounds: alpha, rho_0, rho_sat, C_max
            bounds = [(1, 100), (0.1, 100), (1e2, 1e6), (0.5, 0.99)]

        # Optimize
        try:
            result = differential_evolution(
                objective,
                bounds,
                maxiter=500,
                atol=1e-7,
                seed=42,
                workers=1
            )

            # Extract parameters
            if self.coherence_model == 'rational':
                alpha_best, rho_0_best, rho_sat_best = result.x
                C_max_best = min(0.99, (rho_sat_best / rho_0_best) ** predictor.gamma)
            else:
                alpha_best, rho_0_best, rho_sat_best, C_max_best = result.x

            # Compute reduced chi-squared
            n_data = len(galaxy.radius)
            n_params = len(result.x)
            chi2_red = result.fun / (n_data - n_params)

            # Get final prediction
            v_pred, v_dm, C_vis = predictor.predict_rotation_curve(
                galaxy, alpha_best, rho_0_best, rho_sat_best, C_max_best
            )

            # Compute central density for analysis
            rho_central = galaxy.total_baryonic_density()[0] if len(galaxy.radius) > 0 else np.nan

            return {
                'success': True,
                'name': galaxy.name,
                'alpha': alpha_best,
                'rho_0': rho_0_best,
                'rho_sat': rho_sat_best,
                'C_max': C_max_best,
                'chi2_red': chi2_red,
                'rho_central': rho_central,
                'v_pred': v_pred,
                'v_dm': v_dm,
                'C_vis': C_vis
            }

        except Exception as e:
            print(f"  Fit failed for {galaxy.name}: {e}")
            return {
                'success': False,
                'name': galaxy.name,
                'alpha': np.nan,
                'rho_0': np.nan,
                'rho_sat': np.nan,
                'C_max': np.nan,
                'chi2_red': np.inf,
                'rho_central': np.nan
            }

    def analyze_universality(self, galaxies: List[SPARCGalaxy]) -> Dict:
        """
        Fit all galaxies and analyze ρ_sat distribution.

        Args:
            galaxies: List of SPARC galaxies

        Returns:
            Universality statistics
        """
        print(f"\n{'='*80}")
        print(f"ρ_sat UNIVERSALITY TEST - {self.coherence_model.upper()}")
        print(f"{'='*80}")
        print(f"\nFitting ρ_sat independently for {len(galaxies)} galaxies...")
        print("This may take 10-30 minutes depending on sample size.")

        # Fit all galaxies
        results = []
        for i, galaxy in enumerate(galaxies):
            print(f"[{i+1}/{len(galaxies)}] Fitting {galaxy.name}...")
            result = self.fit_galaxy_with_rho_sat(galaxy)
            results.append(result)

        # Extract successful fits
        successful = [r for r in results if r['success']]
        n_success = len(successful)

        print(f"\n✓ Successfully fit {n_success}/{len(galaxies)} galaxies")

        if n_success < 3:
            print("ERROR: Insufficient successful fits for statistical analysis")
            return {'success': False, 'n_galaxies': n_success}

        # Extract distributions
        rho_sats = np.array([r['rho_sat'] for r in successful])
        C_maxs = np.array([r['C_max'] for r in successful])
        alphas = np.array([r['alpha'] for r in successful])
        chi2_reds = np.array([r['chi2_red'] for r in successful])
        rho_centrals = np.array([r['rho_central'] for r in successful if not np.isnan(r['rho_central'])])

        # Compute statistics
        rho_sat_median = np.median(rho_sats)
        rho_sat_mean = np.mean(rho_sats)
        rho_sat_std = np.std(rho_sats)
        rho_sat_scatter = rho_sat_std / rho_sat_mean  # Fractional scatter

        # Log-scale statistics (more appropriate for multiplicative scatter)
        log_rho_sats = np.log10(rho_sats)
        log_median = np.median(log_rho_sats)
        log_std = np.std(log_rho_sats)
        factor_scatter = 10**log_std  # Typical factor variation

        # Percentiles
        rho_sat_p16 = np.percentile(rho_sats, 16)
        rho_sat_p84 = np.percentile(rho_sats, 84)

        # Success rate
        n_excellent = np.sum(chi2_reds < 2.0)
        n_good = np.sum((chi2_reds >= 2.0) & (chi2_reds < 5.0))
        n_acceptable = np.sum((chi2_reds >= 5.0) & (chi2_reds < 10.0))
        pct_total_good = 100 * (n_excellent + n_good + n_acceptable) / len(galaxies)

        # Universality assessment
        if rho_sat_scatter < 0.3:
            universality = "STRONG"
            verdict = "✅ Universal screening mechanism confirmed"
        elif rho_sat_scatter < 0.5:
            universality = "MODERATE"
            verdict = "✅ Universal with ~factor of 2 variation"
        elif rho_sat_scatter < 1.0:
            universality = "WEAK"
            verdict = "⚠️ Loosely universal (factor of 3)"
        else:
            universality = "FAILED"
            verdict = "❌ Galaxy-dependent, not universal"

        stats = {
            'success': True,
            'model': self.coherence_model,
            'n_galaxies': len(galaxies),
            'n_success': n_success,
            'rho_sat_median': rho_sat_median,
            'rho_sat_mean': rho_sat_mean,
            'rho_sat_std': rho_sat_std,
            'rho_sat_scatter': rho_sat_scatter,
            'rho_sat_p16': rho_sat_p16,
            'rho_sat_p84': rho_sat_p84,
            'log_median': 10**log_median,
            'log_std': log_std,
            'factor_scatter': factor_scatter,
            'C_max_median': np.median(C_maxs),
            'C_max_std': np.std(C_maxs),
            'alpha_median': np.median(alphas),
            'chi2_red_median': np.median(chi2_reds),
            'pct_total_good': pct_total_good,
            'universality': universality,
            'verdict': verdict,
            'results': results
        }

        # Print summary
        self._print_summary(stats)

        return stats

    def _print_summary(self, stats: Dict):
        """Print formatted summary of universality test."""
        print(f"\n{'='*80}")
        print("UNIVERSALITY TEST RESULTS")
        print(f"{'='*80}")

        print(f"\nSample:")
        print(f"  Total galaxies: {stats['n_galaxies']}")
        print(f"  Successful fits: {stats['n_success']}")
        print(f"  Success rate: {100*stats['n_success']/stats['n_galaxies']:.1f}%")

        print(f"\nρ_sat Distribution:")
        print(f"  Median: {stats['rho_sat_median']:.2e} M_☉/pc³")
        print(f"  Mean:   {stats['rho_sat_mean']:.2e} M_☉/pc³")
        print(f"  Std:    {stats['rho_sat_std']:.2e} M_☉/pc³")
        print(f"  16-84 percentile: [{stats['rho_sat_p16']:.2e}, {stats['rho_sat_p84']:.2e}]")

        print(f"\nScatter Analysis:")
        print(f"  Linear scatter: σ/⟨ρ_sat⟩ = {stats['rho_sat_scatter']:.3f} ({100*stats['rho_sat_scatter']:.1f}%)")
        print(f"  Log scatter: σ(log ρ_sat) = {stats['log_std']:.3f}")
        print(f"  Factor variation: ×{stats['factor_scatter']:.2f}")

        print(f"\nUniversality Assessment:")
        print(f"  Classification: {stats['universality']}")
        print(f"  Verdict: {stats['verdict']}")

        print(f"\nSession #18 Prediction:")
        print(f"  Predicted: ρ_sat ≈ 2×10^4 M_☉/pc³")
        print(f"  Measured:  ρ_sat ≈ {stats['rho_sat_median']:.2e} M_☉/pc³")

        if abs(np.log10(stats['rho_sat_median']) - np.log10(2e4)) < 0.5:
            print(f"  STATUS: ✅ Within factor of 2-3 of prediction")
        else:
            print(f"  STATUS: ⚠️ Differs from prediction")

        print(f"\nFit Quality:")
        print(f"  Median χ²_red: {stats['chi2_red_median']:.2f}")
        print(f"  Total success (χ² < 10): {stats['pct_total_good']:.1f}%")

    def analyze_by_galaxy_type(self, galaxies: List[SPARCGalaxy]) -> Dict:
        """
        Analyze ρ_sat universality broken down by galaxy type.

        Args:
            galaxies: List of SPARC galaxies

        Returns:
            Type-specific universality statistics
        """
        print(f"\n{'='*80}")
        print(f"ρ_sat UNIVERSALITY BY GALAXY TYPE - {self.coherence_model.upper()}")
        print(f"{'='*80}")

        # Categorize galaxies
        categorized = self.categorizer.categorize_sample(galaxies)

        # Analyze each type
        type_stats = {}

        for galaxy_type, type_galaxies in categorized.items():
            if len(type_galaxies) < 5:  # Need at least 5 for statistics
                print(f"\nSkipping {galaxy_type} (only {len(type_galaxies)} galaxies)")
                continue

            print(f"\n{'-'*80}")
            print(f"Analyzing {galaxy_type} galaxies ({len(type_galaxies)} total)")
            print(f"{'-'*80}")

            stats = self.analyze_universality(type_galaxies)
            type_stats[galaxy_type] = stats

        # Overall statistics
        print(f"\n{'='*80}")
        print("OVERALL UNIVERSALITY (All Galaxies)")
        print(f"{'='*80}")

        overall_stats = self.analyze_universality(galaxies)

        return {
            'model': self.coherence_model,
            'overall': overall_stats,
            'by_type': type_stats,
            'categorized': categorized
        }


# ============================================================================
# Part 2: Correlation Analysis
# ============================================================================

def analyze_rho_sat_correlations(results: List[Dict]) -> Dict:
    """
    Analyze correlations between ρ_sat and galaxy properties.

    Tests:
    - ρ_sat vs ρ_central (density dependence?)
    - ρ_sat vs galaxy type (morphology dependence?)
    - C_max vs ρ_sat (consistency with rational formula?)

    Args:
        results: List of fit results from universality test

    Returns:
        Correlation statistics
    """
    print(f"\n{'='*80}")
    print("ρ_sat CORRELATION ANALYSIS")
    print(f"{'='*80}")

    # Extract successful fits
    successful = [r for r in results if r['success']]

    if len(successful) < 10:
        print("ERROR: Insufficient data for correlation analysis")
        return {'success': False}

    # Extract arrays
    rho_sats = np.array([r['rho_sat'] for r in successful])
    rho_centrals = np.array([r['rho_central'] for r in successful if not np.isnan(r['rho_central'])])
    C_maxs = np.array([r['C_max'] for r in successful])

    # Correlation: ρ_sat vs ρ_central
    if len(rho_centrals) > 0:
        corr_rho = np.corrcoef(np.log10(rho_sats[:len(rho_centrals)]),
                               np.log10(rho_centrals))[0, 1]
        print(f"\nρ_sat vs ρ_central correlation: r = {corr_rho:.3f}")

        if abs(corr_rho) < 0.3:
            print("  → Weak correlation: ρ_sat largely independent of central density ✅")
        elif abs(corr_rho) < 0.6:
            print("  → Moderate correlation: Some density dependence")
        else:
            print("  → Strong correlation: ρ_sat varies with galaxy density ❌")

    # Correlation: C_max vs ρ_sat (rational formula check)
    corr_Cmax = np.corrcoef(C_maxs, np.log10(rho_sats))[0, 1]
    print(f"\nC_max vs log(ρ_sat) correlation: r = {corr_Cmax:.3f}")

    if abs(corr_Cmax) < 0.3:
        print("  → Weak correlation: C_max and ρ_sat independent")
    else:
        print("  → Correlation present: May indicate rational formula structure")

    return {
        'success': True,
        'corr_rho_central': corr_rho if len(rho_centrals) > 0 else np.nan,
        'corr_C_max': corr_Cmax
    }


# ============================================================================
# Part 3: Main Execution
# ============================================================================

def test_rho_sat_universality():
    """
    Main test: Fit ρ_sat per galaxy and test universality.

    Expected (Session #18):
    - Universal ρ_sat ≈ 2×10^4 M_☉/pc³
    - Scatter σ/⟨ρ_sat⟩ < 0.5 (factor of 2-3)
    """
    print("=" * 80)
    print("SESSION #20 PRIORITY 2: ρ_sat UNIVERSALITY TEST")
    print("Fitting saturation density independently per galaxy")
    print("=" * 80)

    # Load SPARC data
    print("\nLoading SPARC galaxies...")
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=None)  # All 175
    print(f"Loaded {len(galaxies)} galaxies from SPARC database")

    # Test rational formula (recommended from Session #18)
    print("\n" + "="*80)
    print("Testing RATIONAL coherence formula (Session #18 recommended)")
    print("="*80)

    tester = RhoSatUniversalityTester(coherence_model='rational')

    # Overall analysis
    overall_stats = tester.analyze_universality(galaxies)

    # Correlation analysis
    if overall_stats['success']:
        corr_stats = analyze_rho_sat_correlations(overall_stats['results'])

    # Galaxy-type breakdown
    print(f"\n{'='*80}")
    print("GALAXY-TYPE BREAKDOWN")
    print(f"{'='*80}")

    type_stats = tester.analyze_by_galaxy_type(galaxies)

    # Final summary
    print(f"\n{'='*80}")
    print("SESSION #20 PRIORITY 2 COMPLETE")
    print(f"{'='*80}")

    if overall_stats['success']:
        print(f"\nUniversality: {overall_stats['universality']}")
        print(f"Verdict: {overall_stats['verdict']}")
        print(f"\nMedian ρ_sat: {overall_stats['rho_sat_median']:.2e} M_☉/pc³")
        print(f"Scatter: {100*overall_stats['rho_sat_scatter']:.1f}%")
        print(f"Factor variation: ×{overall_stats['factor_scatter']:.2f}")

        print(f"\nSession #18 Prediction Check:")
        print(f"  Predicted: ρ_sat ≈ 2×10^4 M_☉/pc³, scatter < 50%")

        if overall_stats['rho_sat_scatter'] < 0.5:
            print(f"  Result: ✅ VALIDATED - Universal screening mechanism")
        else:
            print(f"  Result: ❌ FAILED - Galaxy-dependent saturation")

    return overall_stats


if __name__ == "__main__":
    # Run ρ_sat universality test
    results = test_rho_sat_universality()

    print("\n✓ Session #20 Priority 2 complete!")
    print("  ρ_sat universality tested across 175 SPARC galaxies")
    print("  Ready for Session #20 final documentation")
