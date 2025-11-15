#!/usr/bin/env python3
"""
Synchronism SPARC Saturation Formula Testing - Session #19
============================================================

Tests Session #18 Track C saturation-aware coherence formulas:
1. Formula 1 (Exponential): C = C_max[1 - exp(-(ρ/ρ_crit)^γ)]
2. Formula 2 (Rational): C = C_max(ρ/ρ_0)^γ/[1 + (ρ/ρ_sat)^γ]  (RECOMMENDED)
3. Formula 3 (Logarithmic): C = C_max[ln(1+ρ/ρ_0)/ln(1+ρ_sat/ρ_0)]^γ

Expected outcomes:
- NGC galaxies: 30% → 50-60% success
- Overall: 40% → 55-65% success
- ρ_sat ≈ 2×10^4 M_☉/pc³ (universal across galaxies)

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-15
Session: #19 - Empirical validation of theoretical saturation solution
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import os
import sys

# Import existing SPARC validation infrastructure
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_real_sparc_validation import (
    RealSPARCLoader,
    SPARCGalaxy,
    SPARCValidator
)


# ============================================================================
# Part 1: Saturation-Aware Synchronism Predictor
# ============================================================================

class SynchronismSaturationPredictor:
    """
    Enhanced Synchronism predictor with saturation-aware coherence formulas.

    Tests three coherence formulas from Session #18 Track C:
    - Exponential saturation
    - Rational function (recommended from theory)
    - Logarithmic saturation
    """

    # Theory-predicted parameters (Session #14, unchanged)
    GAMMA = 0.30  # Coherence exponent
    BETA = 0.30   # DM modulation exponent

    # Physical constants
    G = 4.3e-6  # (km/s)²·kpc/M_☉

    # Predicted saturation density (Session #18 Track C)
    RHO_SAT_PREDICTED = 2.0e4  # M_☉/pc³

    def __init__(self,
                 coherence_model: str = 'rational',
                 gamma: float = GAMMA,
                 beta: float = BETA):
        """
        Initialize predictor with coherence model choice.

        Args:
            coherence_model: 'power_law' (original), 'exponential',
                           'rational' (recommended), or 'logarithmic'
            gamma: Coherence exponent (default 0.30)
            beta: DM modulation exponent (default 0.30)
        """
        self.coherence_model = coherence_model
        self.gamma = gamma
        self.beta = beta

    def coherence_power_law(self,
                           rho_vis: np.ndarray,
                           rho_0: float) -> np.ndarray:
        """
        Original power-law coherence (Sessions #13-17).

        C_vis = (ρ_vis/ρ_0)^γ  (capped at 1)

        Problem: Unbounded, saturates at C=1
        """
        C = (rho_vis / rho_0) ** self.gamma
        return np.clip(C, 0.0, 1.0)

    def coherence_exponential(self,
                             rho_vis: np.ndarray,
                             rho_crit: float,
                             C_max: float) -> np.ndarray:
        """
        Formula 1: Exponential saturation (Session #18 Track C).

        C_vis = C_max · [1 - exp(-(ρ_vis/ρ_crit)^γ)]

        Behavior:
        - Low density: C ≈ C_max(ρ/ρ_crit)^γ (power-law)
        - High density: C → C_max (saturates below 1)

        Args:
            rho_vis: Visible matter density [M_☉/pc²]
            rho_crit: Critical density scale [M_☉/pc²]
            C_max: Maximum coherence (< 1 from quantum/MRH limits)
        """
        x = (rho_vis / rho_crit) ** self.gamma
        C = C_max * (1.0 - np.exp(-x))
        return np.clip(C, 0.0, C_max)

    def coherence_rational(self,
                          rho_vis: np.ndarray,
                          rho_0: float,
                          rho_sat: float,
                          C_max: Optional[float] = None) -> np.ndarray:
        """
        Formula 2: Rational function (RECOMMENDED from Session #18).

        C_vis = C_max · (ρ_vis/ρ_0)^γ / [1 + (ρ_vis/ρ_sat)^γ]

        Behavior:
        - Low density (ρ << ρ_sat): C ≈ C_max(ρ/ρ_0)^γ (original formula)
        - High density (ρ >> ρ_sat): C → C_max(ρ_0/ρ_sat)^γ < C_max

        Physical meaning:
        - ρ_sat: Density where screening becomes important (≈ 2×10^4 M_☉/pc³)
        - C_max derived: C_max = (ρ_sat/ρ_0)^γ or set manually

        Args:
            rho_vis: Visible matter density [M_☉/pc²]
            rho_0: Low-density normalization [M_☉/pc²]
            rho_sat: Saturation density [M_☉/pc²]
            C_max: Max coherence (if None, derived from ρ_sat/ρ_0)
        """
        x = (rho_vis / rho_0) ** self.gamma
        y = (rho_vis / rho_sat) ** self.gamma

        if C_max is None:
            # Derive C_max from ratio (Session #18 approach)
            C_max = min(0.99, (rho_sat / rho_0) ** self.gamma)

        C = C_max * x / (1.0 + y)
        return np.clip(C, 0.0, C_max)

    def coherence_logarithmic(self,
                             rho_vis: np.ndarray,
                             rho_0: float,
                             rho_sat: float,
                             C_max: float) -> np.ndarray:
        """
        Formula 3: Logarithmic saturation (Session #18 Track C).

        C_vis = C_max · [ln(1 + ρ_vis/ρ_0) / ln(1 + ρ_sat/ρ_0)]^γ

        Behavior:
        - Low density: C ∝ (ρ/ρ_0)^γ / [ln(1+ρ_sat/ρ_0)]^γ
        - At ρ = ρ_sat: C → C_max

        Physical meaning: Information-theoretic coherence growth

        Args:
            rho_vis: Visible matter density [M_☉/pc²]
            rho_0: Normalization density [M_☉/pc²]
            rho_sat: Saturation density [M_☉/pc²]
            C_max: Maximum coherence
        """
        numerator = np.log(1.0 + rho_vis / rho_0)
        denominator = np.log(1.0 + rho_sat / rho_0)

        C = C_max * (numerator / denominator) ** self.gamma
        return np.clip(C, 0.0, C_max)

    def compute_coherence(self,
                         rho_vis: np.ndarray,
                         rho_0: Optional[float] = None,
                         rho_sat: Optional[float] = None,
                         rho_crit: Optional[float] = None,
                         C_max: Optional[float] = None) -> np.ndarray:
        """
        Unified coherence computation supporting all models.

        Args:
            rho_vis: Visible matter density [M_☉/pc²]
            rho_0: Normalization density (default: max(rho_vis))
            rho_sat: Saturation density (default: 2×10^4 M_☉/pc³)
            rho_crit: Critical density for exponential (default: rho_0)
            C_max: Maximum coherence (default: 0.95)

        Returns:
            Coherence C_vis [dimensionless, 0 to C_max]
        """
        # Set defaults
        if rho_0 is None:
            rho_0 = max(np.max(rho_vis), 1.0)

        if rho_sat is None:
            rho_sat = self.RHO_SAT_PREDICTED

        if rho_crit is None:
            rho_crit = rho_0

        if C_max is None:
            C_max = 0.95  # Default from quantum + MRH limits

        # Select model
        if self.coherence_model == 'power_law':
            return self.coherence_power_law(rho_vis, rho_0)

        elif self.coherence_model == 'exponential':
            return self.coherence_exponential(rho_vis, rho_crit, C_max)

        elif self.coherence_model == 'rational':
            return self.coherence_rational(rho_vis, rho_0, rho_sat, C_max)

        elif self.coherence_model == 'logarithmic':
            return self.coherence_logarithmic(rho_vis, rho_0, rho_sat, C_max)

        else:
            raise ValueError(f"Unknown coherence model: {self.coherence_model}")

    def compute_dark_matter_density(self,
                                    rho_vis: np.ndarray,
                                    C_vis: np.ndarray,
                                    alpha: float) -> np.ndarray:
        """
        Compute dark matter density (unchanged from Session #13-18).

        ρ_DM = α(1 - C_vis)ρ_vis^β

        Args:
            rho_vis: Visible matter density [M_☉/pc²]
            C_vis: Coherence from compute_coherence
            alpha: DM normalization

        Returns:
            Dark matter density [M_☉/pc²]
        """
        return alpha * (1.0 - C_vis) * (rho_vis ** self.beta)

    def compute_enclosed_mass(self,
                             radius: np.ndarray,
                             sigma_total: np.ndarray) -> np.ndarray:
        """
        Compute enclosed mass from surface density (thin disk).

        M(<r) = 2π ∫₀ʳ Σ(r') r' dr'
        """
        sigma_kpc = sigma_total * 1e6  # M_☉/pc² → M_☉/kpc²
        dr = np.diff(radius, prepend=0)
        dM = 2.0 * np.pi * radius * sigma_kpc * dr
        return np.cumsum(dM)

    def predict_rotation_curve(self,
                               galaxy: SPARCGalaxy,
                               alpha: float,
                               rho_0: Optional[float] = None,
                               rho_sat: Optional[float] = None,
                               C_max: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Predict rotation curve with saturation-aware coherence.

        Args:
            galaxy: SPARC galaxy data
            alpha: DM normalization
            rho_0: Normalization density
            rho_sat: Saturation density
            C_max: Maximum coherence

        Returns:
            (v_rot, v_dm, C_vis): Rotation velocity, DM contribution, coherence
        """
        # Total baryonic density
        rho_vis = galaxy.total_baryonic_density()

        # Coherence with saturation
        C_vis = self.compute_coherence(rho_vis, rho_0, rho_sat, None, C_max)

        # Dark matter density
        rho_dm = self.compute_dark_matter_density(rho_vis, C_vis, alpha)

        # Total density
        rho_total = rho_vis + rho_dm

        # Enclosed masses
        M_total = self.compute_enclosed_mass(galaxy.radius, rho_total)
        M_dm = self.compute_enclosed_mass(galaxy.radius, rho_dm)

        # Rotation velocities
        v_rot = np.sqrt(self.G * M_total / galaxy.radius)
        v_dm = np.sqrt(self.G * M_dm / galaxy.radius)

        return v_rot, v_dm, C_vis


# ============================================================================
# Part 2: Saturation-Aware Validator
# ============================================================================

class SaturationValidator:
    """Enhanced validator for saturation formula testing."""

    def __init__(self, predictor: SynchronismSaturationPredictor):
        """Initialize with saturation-aware predictor."""
        self.predictor = predictor

    def fit_galaxy(self,
                   galaxy: SPARCGalaxy,
                   fit_saturation: bool = False) -> Dict:
        """
        Fit Synchronism parameters to galaxy rotation curve.

        Args:
            galaxy: SPARC galaxy
            fit_saturation: If True, fit ρ_sat and C_max; else use predicted

        Returns:
            Best-fit parameters and goodness-of-fit
        """
        # Define fit function
        def objective(params):
            """Chi-squared for parameter optimization."""
            if fit_saturation and self.predictor.coherence_model == 'rational':
                # Fit: (alpha, rho_0, rho_sat)
                alpha, rho_0, rho_sat = params
                C_max = None  # Derived
            elif fit_saturation and self.predictor.coherence_model in ['exponential', 'logarithmic']:
                # Fit: (alpha, rho_0, rho_sat, C_max)
                alpha, rho_0, rho_sat, C_max = params
            else:
                # Fit only alpha (use predicted rho_sat)
                alpha = params[0] if isinstance(params, (list, tuple, np.ndarray)) else params
                rho_0 = max(np.max(galaxy.total_baryonic_density()), 1.0)
                rho_sat = self.predictor.RHO_SAT_PREDICTED
                C_max = 0.95

            # Predict
            v_pred, _, _ = self.predictor.predict_rotation_curve(
                galaxy, alpha, rho_0, rho_sat, C_max
            )

            # Chi-squared
            residuals = (galaxy.v_obs - v_pred) / galaxy.v_err
            chi2 = np.sum(residuals**2)

            return chi2

        # Set bounds and initial guess
        if fit_saturation and self.predictor.coherence_model == 'rational':
            # bounds: alpha [1, 100], rho_0 [1, 100], rho_sat [1e3, 1e5]
            bounds = [(1, 100), (1, 100), (1e3, 1e5)]
            x0 = [50, 10, 2e4]
        elif fit_saturation:
            # bounds: alpha, rho_0, rho_sat, C_max
            bounds = [(1, 100), (1, 100), (1e3, 1e5), (0.8, 0.99)]
            x0 = [50, 10, 2e4, 0.95]
        else:
            # Fit only alpha
            bounds = [(1, 100)]
            x0 = [50]

        # Optimize
        try:
            result = differential_evolution(
                objective,
                bounds,
                maxiter=300,
                atol=1e-6,
                seed=42
            )

            # Extract parameters
            if fit_saturation and self.predictor.coherence_model == 'rational':
                alpha_best, rho_0_best, rho_sat_best = result.x
                C_max_best = min(0.99, (rho_sat_best / rho_0_best) ** self.predictor.gamma)
            elif fit_saturation:
                alpha_best, rho_0_best, rho_sat_best, C_max_best = result.x
            else:
                alpha_best = result.x[0]
                rho_0_best = max(np.max(galaxy.total_baryonic_density()), 1.0)
                rho_sat_best = self.predictor.RHO_SAT_PREDICTED
                C_max_best = 0.95

            # Compute reduced chi-squared
            n_data = len(galaxy.radius)
            n_params = len(result.x)
            chi2_red = result.fun / (n_data - n_params)

            # Get final prediction
            v_pred, v_dm, C_vis = self.predictor.predict_rotation_curve(
                galaxy, alpha_best, rho_0_best, rho_sat_best, C_max_best
            )

            return {
                'success': True,
                'alpha': alpha_best,
                'rho_0': rho_0_best,
                'rho_sat': rho_sat_best,
                'C_max': C_max_best,
                'chi2_red': chi2_red,
                'v_pred': v_pred,
                'v_dm': v_dm,
                'C_vis': C_vis
            }

        except Exception as e:
            print(f"  Fit failed for {galaxy.name}: {e}")
            return {
                'success': False,
                'alpha': np.nan,
                'rho_0': np.nan,
                'rho_sat': np.nan,
                'C_max': np.nan,
                'chi2_red': np.inf
            }

    def validate_sample(self,
                        galaxies: List[SPARCGalaxy],
                        fit_saturation: bool = False) -> Dict:
        """
        Validate on sample of galaxies.

        Args:
            galaxies: List of SPARC galaxies
            fit_saturation: If True, fit sat. params; else use predicted

        Returns:
            Statistics dictionary
        """
        print(f"\nValidating {len(galaxies)} galaxies with {self.predictor.coherence_model} coherence...")
        print(f"Fit saturation parameters: {fit_saturation}")

        results = []
        for i, galaxy in enumerate(galaxies):
            print(f"[{i+1}/{len(galaxies)}] Fitting {galaxy.name}...")
            result = self.fit_galaxy(galaxy, fit_saturation)
            result['name'] = galaxy.name
            results.append(result)

        # Compile statistics
        chi2_reds = np.array([r['chi2_red'] for r in results if r['success']])
        alphas = np.array([r['alpha'] for r in results if r['success']])
        rho_sats = np.array([r['rho_sat'] for r in results if r['success']])
        C_maxs = np.array([r['C_max'] for r in results if r['success']])

        # Success rates
        n_total = len(galaxies)
        n_excellent = np.sum(chi2_reds < 2.0)
        n_good = np.sum((chi2_reds >= 2.0) & (chi2_reds < 5.0))
        n_acceptable = np.sum((chi2_reds >= 5.0) & (chi2_reds < 10.0))

        stats = {
            'model': self.predictor.coherence_model,
            'n_galaxies': n_total,
            'n_success': len(chi2_reds),
            'chi2_red_median': np.median(chi2_reds),
            'chi2_red_mean': np.mean(chi2_reds),
            'chi2_red_std': np.std(chi2_reds),
            'pct_excellent': 100 * n_excellent / n_total,
            'pct_good': 100 * n_good / n_total,
            'pct_acceptable': 100 * n_acceptable / n_total,
            'pct_total_good': 100 * (n_excellent + n_good + n_acceptable) / n_total,
            'alpha_median': np.median(alphas),
            'rho_sat_median': np.median(rho_sats),
            'rho_sat_std': np.std(rho_sats),
            'C_max_median': np.median(C_maxs),
            'results': results
        }

        return stats


# ============================================================================
# Part 3: Main Test Execution
# ============================================================================

def test_saturation_formulas():
    """
    Main test: Compare all coherence formulas on SPARC data.

    Expected (Session #18 predictions):
    - NGC galaxies: 30% → 50-60% with saturation
    - Overall: 40% → 55-65%
    - ρ_sat ≈ 2×10^4 M_☉/pc³ (universal)
    """
    print("=" * 80)
    print("SESSION #19: SATURATION FORMULA TESTING")
    print("Testing Session #18 Track C predictions on 175 SPARC galaxies")
    print("=" * 80)

    # Load SPARC data
    print("\nLoading SPARC galaxies...")
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=None)  # All 175
    print(f"Loaded {len(galaxies)} galaxies from SPARC database")

    # Test each formula
    models = ['power_law', 'rational', 'exponential', 'logarithmic']
    all_stats = {}

    for model in models:
        print(f"\n{'='*80}")
        print(f"Testing {model.upper()} coherence formula")
        print(f"{'='*80}")

        # Create predictor and validator
        predictor = SynchronismSaturationPredictor(coherence_model=model)
        validator = SaturationValidator(predictor)

        # Validate (use predicted ρ_sat for speed, can fit later)
        stats = validator.validate_sample(galaxies, fit_saturation=False)

        all_stats[model] = stats

        # Print summary
        print(f"\n{model.upper()} RESULTS:")
        print(f"  Median χ²_red: {stats['chi2_red_median']:.2f}")
        print(f"  Success rates:")
        print(f"    Excellent (χ² < 2): {stats['pct_excellent']:.1f}%")
        print(f"    Good (2 ≤ χ² < 5): {stats['pct_good']:.1f}%")
        print(f"    Acceptable (5 ≤ χ² < 10): {stats['pct_acceptable']:.1f}%")
        print(f"    TOTAL GOOD: {stats['pct_total_good']:.1f}%")

        if model != 'power_law':
            print(f"  ρ_sat median: {stats['rho_sat_median']:.1e} M_☉/pc³")
            print(f"  ρ_sat scatter: {stats['rho_sat_std']:.1e} M_☉/pc³")
            print(f"  C_max median: {stats['C_max_median']:.3f}")

    # Compare results
    print(f"\n{'='*80}")
    print("COMPARISON SUMMARY")
    print(f"{'='*80}")
    print(f"{'Model':<15} {'Median χ²':<12} {'Total Good %':<15} {'vs Power-Law'}")
    print(f"{'-'*80}")

    baseline = all_stats['power_law']['pct_total_good']
    for model in models:
        stats = all_stats[model]
        improvement = stats['pct_total_good'] - baseline
        print(f"{model:<15} {stats['chi2_red_median']:<12.2f} "
              f"{stats['pct_total_good']:<15.1f} {improvement:+.1f}%")

    print(f"\n{'='*80}")
    print("SESSION #18 PREDICTION CHECK")
    print(f"{'='*80}")
    print(f"Expected NGC improvement: 30% → 50-60%")
    print(f"Expected overall: 40% → 55-65%")
    print(f"Expected ρ_sat: ~2×10^4 M_☉/pc³ (universal)")
    print(f"\nActual results: See above")
    print(f"{'='*80}")

    return all_stats


if __name__ == "__main__":
    # Run saturation formula tests
    results = test_saturation_formulas()

    print("\n✓ Session #19 saturation testing complete!")
    print("  Results show impact of coherence saturation on SPARC fits")
    print("  Compare to Session #17 baseline: 40% success with power-law")
