#!/usr/bin/env python3
"""
Session #23: Galaxy-Specific ρ_mag Model
=========================================

Tests whether galaxy-specific ρ_mag (based on B-field strength) can explain
Session #22's NGC underprediction.

Hypothesis:
  ρ_mag = ρ_mag,0 × (B_ref / B)^α

  NGC galaxies (high B ~ 10 μG) → low ρ_mag → low ρ_sat ✓
  F galaxies (low B ~ 3 μG) → high ρ_mag → high ρ_sat ✓

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-17
Session: #23 - Refinement following Session #22 NGC underprediction
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle
import os
import sys

# Import Session #22 data
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_session22_extract_data import load_session20_data


# ============================================================================
# Part 1: B-field Estimation from Morphology
# ============================================================================

def estimate_b_field_morphology(galaxy_name: str) -> float:
    """
    Estimate magnetic field strength from galaxy morphology.

    Based on literature (Beck & Wielebinski 2013):
    - NGC/UGC (spirals): B ~ 8-12 μG
    - F (irregulars): B ~ 2-4 μG
    - DDO (dwarfs): B ~ 1-2 μG
    - IC/ESO (mixed): B ~ 5-8 μG

    Args:
        galaxy_name: Galaxy identifier (e.g., "NGC2403")

    Returns:
        B-field strength [μG]
    """
    # Extract prefix
    if galaxy_name.startswith('NGC'):
        return 10.0  # Spirals
    elif galaxy_name.startswith('UGC'):
        return 10.0  # Spirals
    elif galaxy_name.startswith('F'):
        return 3.0   # Irregulars
    elif galaxy_name.startswith('DDO'):
        return 1.5   # Dwarfs
    elif galaxy_name.startswith('IC'):
        return 6.0   # Mixed
    elif galaxy_name.startswith('ESO'):
        return 6.0   # Mixed
    else:
        return 5.0   # Default


# ============================================================================
# Part 2: Galaxy-Specific ρ_mag Model
# ============================================================================

def rho_mag_galaxy_specific(B: np.ndarray,
                            rho_mag_0: float,
                            B_ref: float,
                            alpha: float) -> np.ndarray:
    """
    Galaxy-specific ρ_mag from B-field strength.

    ρ_mag = ρ_mag,0 × (B_ref / B)^α

    Physical interpretation:
    - Strong B-field → low ρ_mag → early saturation → low ρ_sat
    - Weak B-field → high ρ_mag → late saturation → high ρ_sat

    Args:
        B: Magnetic field strength [μG]
        rho_mag_0: Baseline ρ_mag for B = B_ref [M_☉/pc³]
        B_ref: Reference B-field [μG]
        alpha: Power-law index [dimensionless]

    Returns:
        Galaxy-specific ρ_mag [M_☉/pc³]
    """
    return rho_mag_0 * (B_ref / B)**alpha


def rho_sat_galaxy_specific_model(rho_central: np.ndarray,
                                   B_fields: np.ndarray,
                                   rho_sat_0: float,
                                   rho_mag_0: float,
                                   B_ref: float,
                                   alpha: float,
                                   delta: float) -> np.ndarray:
    """
    Magnetic screening model with galaxy-specific ρ_mag.

    ρ_sat = ρ_sat,0 / [1 + (ρ_c / ρ_mag(B))^δ]

    where ρ_mag(B) = ρ_mag,0 × (B_ref / B)^α

    Args:
        rho_central: Central density [M_☉/pc³]
        B_fields: Magnetic field strengths [μG]
        rho_sat_0: Baseline saturation density [M_☉/pc³]
        rho_mag_0: Baseline ρ_mag [M_☉/pc³]
        B_ref: Reference B-field [μG]
        alpha: B-field power-law index
        delta: Density power-law index

    Returns:
        ρ_sat [M_☉/pc³]
    """
    # Compute galaxy-specific ρ_mag
    rho_mag = rho_mag_galaxy_specific(B_fields, rho_mag_0, B_ref, alpha)

    # Magnetic screening formula
    rho_sat = rho_sat_0 / (1.0 + (rho_central / rho_mag)**delta)

    return rho_sat


# ============================================================================
# Part 3: Model Fitting and Validation
# ============================================================================

class GalaxySpecificValidator:
    """
    Tests galaxy-specific ρ_mag model against Session #22 real data.
    """

    def __init__(self, data: dict):
        """
        Initialize with Session #20 data.

        Args:
            data: Dictionary from load_session20_data()
        """
        self.rho_centrals = data['rho_centrals']
        self.rho_sats_fitted = data['rho_sats_fitted']
        self.galaxy_names = data['galaxy_names']
        self.galaxy_types = data['galaxy_types']

        # Estimate B-fields
        self.B_fields = np.array([estimate_b_field_morphology(name) for name in self.galaxy_names])

        print(f"✓ Loaded {len(self.rho_centrals)} galaxies from Session #20")
        print(f"  B-field range: {self.B_fields.min():.1f} - {self.B_fields.max():.1f} μG")

    def fit_galaxy_specific_model(self) -> dict:
        """
        Fit galaxy-specific ρ_mag model to real data.

        Free parameters: (ρ_sat,0, ρ_mag,0, B_ref, α, δ)

        Returns:
            Best-fit parameters and statistics
        """
        print("\n" + "="*80)
        print("FITTING GALAXY-SPECIFIC ρ_mag MODEL")
        print("="*80)

        print(f"\nModel: ρ_sat = ρ_sat,0 / [1 + (ρ_c / ρ_mag(B))^δ]")
        print(f"       ρ_mag(B) = ρ_mag,0 × (B_ref / B)^α")
        print(f"Free parameters: (ρ_sat,0, ρ_mag,0, B_ref, α, δ)")

        # Wrapper for curve_fit
        def model_wrapper(rho_c, rho_sat_0, rho_mag_0, B_ref, alpha, delta):
            return rho_sat_galaxy_specific_model(
                rho_c, self.B_fields, rho_sat_0, rho_mag_0, B_ref, alpha, delta
            )

        # Initial guess (informed by Session #22)
        # Session #22: ρ_sat,0 = 6.93e5, ρ_mag = 8.65e1, δ = 1.85
        # New: ρ_mag,0 ~ 1e3, B_ref ~ 5 μG, α ~ 1.5
        p0 = [7e5, 1e3, 5.0, 1.5, 1.85]

        # Bounds
        bounds = (
            [1e4, 1e2, 1.0, 0.5, 0.5],   # lower
            [1e7, 1e5, 20.0, 3.0, 3.0]   # upper
        )

        # Fit
        try:
            popt, pcov = curve_fit(
                model_wrapper,
                self.rho_centrals,
                self.rho_sats_fitted,
                p0=p0,
                bounds=bounds,
                maxfev=10000
            )

            rho_sat_0, rho_mag_0, B_ref, alpha, delta = popt
            errors = np.sqrt(np.diag(pcov))

            print(f"\n✅ FIT SUCCESSFUL!")

            print(f"\n{'='*80}")
            print("BEST-FIT PARAMETERS (Galaxy-Specific ρ_mag)")
            print(f"{'='*80}")
            print(f"  ρ_sat,0  = {rho_sat_0:.2e} ± {errors[0]:.2e} M_☉/pc³")
            print(f"  ρ_mag,0  = {rho_mag_0:.2e} ± {errors[1]:.2e} M_☉/pc³  (at B = B_ref)")
            print(f"  B_ref    = {B_ref:.2f} ± {errors[2]:.2f} μG")
            print(f"  α        = {alpha:.2f} ± {errors[3]:.2f}  (B-field index)")
            print(f"  δ        = {delta:.2f} ± {errors[4]:.2f}  (density index)")

            # Predict
            rho_sats_predicted = rho_sat_galaxy_specific_model(
                self.rho_centrals, self.B_fields,
                rho_sat_0, rho_mag_0, B_ref, alpha, delta
            )

            # Statistics
            residuals = self.rho_sats_fitted - rho_sats_predicted
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.rho_sats_fitted - np.mean(self.rho_sats_fitted))**2)
            r_squared = 1 - ss_res / ss_tot

            # Log-scale residuals
            log_residuals = np.log10(self.rho_sats_fitted) - np.log10(rho_sats_predicted)
            log_rms = np.sqrt(np.mean(log_residuals**2))

            # Correlation
            corr = np.corrcoef(np.log10(self.rho_centrals), np.log10(self.rho_sats_fitted))[0, 1]

            print(f"\n{'='*80}")
            print("FIT QUALITY")
            print(f"{'='*80}")
            print(f"  R²                              = {r_squared:.3f}")
            print(f"  RMS(log₁₀ residuals)            = {log_rms:.3f}")
            print(f"  Correlation(log ρ_c, log ρ_sat) = {corr:.3f}")

            # Compare with Session #22 universal model
            print(f"\n{'='*80}")
            print("COMPARISON WITH SESSION #22 (Universal ρ_mag)")
            print(f"{'='*80}")
            print(f"  Metric           Session #22     Session #23")
            print(f"  R²               0.406           {r_squared:.3f}   {'+' if r_squared > 0.406 else '-'}{abs(r_squared - 0.406):.3f}")
            print(f"  RMS(log)         1.833           {log_rms:.3f}   {'-' if log_rms < 1.833 else '+'}{abs(log_rms - 1.833):.3f}")

            if r_squared > 0.406:
                improvement = ((r_squared - 0.406) / 0.406) * 100
                print(f"\n  ✅ R² IMPROVED by {improvement:.1f}%")
            else:
                print(f"\n  ❌ No R² improvement")

            return {
                'success': True,
                'rho_sat_0': rho_sat_0,
                'rho_mag_0': rho_mag_0,
                'B_ref': B_ref,
                'alpha': alpha,
                'delta': delta,
                'errors': errors,
                'r_squared': r_squared,
                'log_rms': log_rms,
                'correlation': corr,
                'predicted': rho_sats_predicted,
                'residuals': residuals,
                'log_residuals': log_residuals
            }

        except Exception as e:
            print(f"\n❌ FIT FAILED: {e}")
            return {'success': False}

    def analyze_by_galaxy_type(self, fit_results: dict):
        """
        Analyze galaxy-specific ρ_mag values by type.

        Args:
            fit_results: Output from fit_galaxy_specific_model()
        """
        print(f"\n{'='*80}")
        print("GALAXY-TYPE ρ_mag DISTRIBUTION")
        print(f"{'='*80}")

        # Compute galaxy-specific ρ_mag for each galaxy
        rho_mags = rho_mag_galaxy_specific(
            self.B_fields,
            fit_results['rho_mag_0'],
            fit_results['B_ref'],
            fit_results['alpha']
        )

        # Breakdown by type
        unique_types = sorted(set(self.galaxy_types))

        print(f"\n{'Type':<10} {'N':<5} {'B_avg':<10} {'ρ_mag (median)':<20} {'ρ_sat (obs)':<20} {'ρ_sat (pred)':<20}")
        print("-"*100)

        for gtype in unique_types:
            indices = [i for i, t in enumerate(self.galaxy_types) if t == gtype]
            if len(indices) == 0:
                continue

            B_avg = np.mean([self.B_fields[i] for i in indices])
            rho_mag_median = np.median([rho_mags[i] for i in indices])
            rho_sat_obs_median = np.median([self.rho_sats_fitted[i] for i in indices])
            rho_sat_pred_median = np.median([fit_results['predicted'][i] for i in indices])

            print(f"{gtype:<10} {len(indices):<5} {B_avg:<10.1f} {rho_mag_median:<20.2e} {rho_sat_obs_median:<20.2e} {rho_sat_pred_median:<20.2e}")

        # Test F/NGC ratio improvement
        print(f"\n{'='*80}")
        print("F/NGC RATIO TEST")
        print(f"{'='*80}")

        ngc_indices = [i for i, t in enumerate(self.galaxy_types) if t == 'NGC']
        f_indices = [i for i, t in enumerate(self.galaxy_types) if t == 'F']

        if len(ngc_indices) > 0 and len(f_indices) > 0:
            # Observed
            rho_sat_ngc_obs = np.median([self.rho_sats_fitted[i] for i in ngc_indices])
            rho_sat_f_obs = np.median([self.rho_sats_fitted[i] for i in f_indices])
            ratio_obs = rho_sat_f_obs / rho_sat_ngc_obs

            # Predicted (galaxy-specific)
            rho_sat_ngc_pred = np.median([fit_results['predicted'][i] for i in ngc_indices])
            rho_sat_f_pred = np.median([fit_results['predicted'][i] for i in f_indices])
            ratio_pred = rho_sat_f_pred / rho_sat_ngc_pred

            # Session #22 (universal)
            ratio_pred_s22 = 11.4  # From Session #22

            print(f"\nObserved F/NGC ratio:           {ratio_obs:>8.1f}×")
            print(f"Predicted (Session #22):        {ratio_pred_s22:>8.1f}×   (error: {ratio_obs/ratio_pred_s22:>6.1f}×)")
            print(f"Predicted (Session #23):        {ratio_pred:>8.1f}×   (error: {ratio_obs/ratio_pred:>6.1f}×)")

            if ratio_obs / ratio_pred < ratio_obs / ratio_pred_s22:
                improvement = ratio_pred_s22 / ratio_pred
                print(f"\n  ✅ F/NGC RATIO IMPROVED by factor of {improvement:.1f}×")
            else:
                print(f"\n  ❌ No F/NGC ratio improvement")

    def plot_results(self, fit_results: dict, output_file: str = "session23_galaxy_specific_rho_mag.png"):
        """
        Create comprehensive validation plots.

        Args:
            fit_results: Output from fit_galaxy_specific_model()
            output_file: Output filename
        """
        if not fit_results['success']:
            print("Cannot plot - fit failed")
            return

        print(f"\nCreating validation plots: {output_file}")

        fig = plt.figure(figsize=(18, 12))
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

        # Color by galaxy type
        type_colors = {'NGC': 'red', 'UGC': 'blue', 'F': 'green', 'DDO': 'orange', 'Other': 'gray'}

        # Plot 1: Data and model fit (colored by type)
        ax1 = fig.add_subplot(gs[0, 0])

        for gtype, color in type_colors.items():
            indices = [i for i, t in enumerate(self.galaxy_types) if t == gtype]
            if len(indices) > 0:
                ax1.loglog(
                    [self.rho_centrals[i] for i in indices],
                    [self.rho_sats_fitted[i] for i in indices],
                    'o', color=color, alpha=0.6, ms=6, label=f"{gtype} (n={len(indices)})"
                )

        # Model prediction (show curves for different B-fields)
        rho_c_range = np.logspace(0, 5, 200)
        for B_val, label in [(1.5, 'DDO (B=1.5μG)'), (3.0, 'F (B=3μG)'), (10.0, 'NGC (B=10μG)')]:
            B_array = np.full_like(rho_c_range, B_val)
            rho_sat_model = rho_sat_galaxy_specific_model(
                rho_c_range, B_array,
                fit_results['rho_sat_0'],
                fit_results['rho_mag_0'],
                fit_results['B_ref'],
                fit_results['alpha'],
                fit_results['delta']
            )
            ax1.loglog(rho_c_range, rho_sat_model, '--', lw=2, label=label)

        ax1.set_xlabel("Central Density ρ_central [M_☉/pc³]", fontsize=13)
        ax1.set_ylabel("Saturation Density ρ_sat [M_☉/pc³]", fontsize=13)
        ax1.set_title(f"Session #23: Galaxy-Specific ρ_mag (R² = {fit_results['r_squared']:.3f})", fontsize=14, fontweight='bold')
        ax1.legend(fontsize=9, loc='best', ncol=2)
        ax1.grid(True, alpha=0.3, which='both')

        # Plot 2: Residuals vs B-field
        ax2 = fig.add_subplot(gs[0, 1])

        for gtype, color in type_colors.items():
            indices = [i for i, t in enumerate(self.galaxy_types) if t == gtype]
            if len(indices) > 0:
                ax2.scatter(
                    [self.B_fields[i] for i in indices],
                    [fit_results['log_residuals'][i] for i in indices],
                    c=color, alpha=0.6, s=40
                )

        ax2.axhline(0, ls='-', c='k', lw=2)
        ax2.axhline(fit_results['log_rms'], ls='--', c='r', lw=1.5, label=f"RMS = {fit_results['log_rms']:.3f}")
        ax2.axhline(-fit_results['log_rms'], ls='--', c='r', lw=1.5)

        ax2.set_xlabel("Magnetic Field B [μG]", fontsize=13)
        ax2.set_ylabel("log₁₀(ρ_sat,obs / ρ_sat,model)", fontsize=13)
        ax2.set_title("Residuals vs B-field", fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        # Plot 3: Observed vs Predicted
        ax3 = fig.add_subplot(gs[0, 2])

        for gtype, color in type_colors.items():
            indices = [i for i, t in enumerate(self.galaxy_types) if t == gtype]
            if len(indices) > 0:
                ax3.loglog(
                    [fit_results['predicted'][i] for i in indices],
                    [self.rho_sats_fitted[i] for i in indices],
                    'o', color=color, alpha=0.6, ms=6
                )

        # 1:1 line
        lims = [
            min(fit_results['predicted'].min(), self.rho_sats_fitted.min()),
            max(fit_results['predicted'].max(), self.rho_sats_fitted.max())
        ]
        ax3.loglog(lims, lims, 'k--', lw=2, alpha=0.7, label='1:1 line')
        ax3.loglog(lims, [10*l for l in lims], 'gray', ls=':', lw=1, alpha=0.5)
        ax3.loglog(lims, [0.1*l for l in lims], 'gray', ls=':', lw=1, alpha=0.5)

        ax3.set_xlabel("ρ_sat (model) [M_☉/pc³]", fontsize=13)
        ax3.set_ylabel("ρ_sat (observed) [M_☉/pc³]", fontsize=13)
        ax3.set_title("Observed vs Predicted", fontsize=14, fontweight='bold')
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3, which='both')

        # Plot 4: ρ_mag distribution by galaxy type
        ax4 = fig.add_subplot(gs[1, 0])

        rho_mags = rho_mag_galaxy_specific(
            self.B_fields,
            fit_results['rho_mag_0'],
            fit_results['B_ref'],
            fit_results['alpha']
        )

        for gtype, color in type_colors.items():
            indices = [i for i, t in enumerate(self.galaxy_types) if t == gtype]
            if len(indices) > 0:
                rho_mag_vals = [rho_mags[i] for i in indices]
                ax4.hist(np.log10(rho_mag_vals), bins=10, alpha=0.6, color=color, label=gtype)

        ax4.set_xlabel("log₁₀(ρ_mag) [M_☉/pc³]", fontsize=13)
        ax4.set_ylabel("Count", fontsize=13)
        ax4.set_title("ρ_mag Distribution by Type", fontsize=14, fontweight='bold')
        ax4.legend(fontsize=10)
        ax4.grid(True, alpha=0.3, axis='y')

        # Plot 5: ρ_mag vs B-field
        ax5 = fig.add_subplot(gs[1, 1])

        B_range = np.linspace(1, 15, 100)
        rho_mag_curve = rho_mag_galaxy_specific(
            B_range,
            fit_results['rho_mag_0'],
            fit_results['B_ref'],
            fit_results['alpha']
        )

        ax5.loglog(B_range, rho_mag_curve, 'k-', lw=2, label=f"ρ_mag ∝ B^(-{fit_results['alpha']:.2f})")

        for gtype, color in type_colors.items():
            indices = [i for i, t in enumerate(self.galaxy_types) if t == gtype]
            if len(indices) > 0:
                ax5.scatter(
                    [self.B_fields[i] for i in indices],
                    [rho_mags[i] for i in indices],
                    c=color, alpha=0.6, s=40, label=gtype
                )

        ax5.set_xlabel("Magnetic Field B [μG]", fontsize=13)
        ax5.set_ylabel("ρ_mag [M_☉/pc³]", fontsize=13)
        ax5.set_title(f"ρ_mag vs B-field (α = {fit_results['alpha']:.2f})", fontsize=14, fontweight='bold')
        ax5.legend(fontsize=10, ncol=2)
        ax5.grid(True, alpha=0.3, which='both')

        # Plot 6: R² comparison
        ax6 = fig.add_subplot(gs[1, 2])

        models = ['Session #22\n(Universal)', 'Session #23\n(Galaxy-Specific)']
        r_squared_values = [0.406, fit_results['r_squared']]
        colors_bar = ['lightblue', 'lightgreen' if fit_results['r_squared'] > 0.406 else 'lightcoral']

        bars = ax6.bar(models, r_squared_values, color=colors_bar, edgecolor='black', linewidth=2)
        ax6.set_ylabel("R²", fontsize=13)
        ax6.set_title("Model Comparison", fontsize=14, fontweight='bold')
        ax6.set_ylim(0, max(r_squared_values) * 1.2)
        ax6.grid(True, alpha=0.3, axis='y')

        # Add value labels on bars
        for bar, val in zip(bars, r_squared_values):
            height = bar.get_height()
            ax6.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ Validation plots saved: {output_file}")

        plt.close()


# ============================================================================
# Part 4: Main Execution
# ============================================================================

def test_galaxy_specific_rho_mag():
    """
    Main test: Validate galaxy-specific ρ_mag model.
    """
    print("="*80)
    print("SESSION #23: GALAXY-SPECIFIC ρ_mag MODEL")
    print("Testing whether ρ_mag ∝ 1/B^α explains NGC underprediction")
    print("="*80)

    # Load Session #20 data
    data = load_session20_data()

    # Initialize validator
    validator = GalaxySpecificValidator(data)

    # Fit galaxy-specific model
    fit_results = validator.fit_galaxy_specific_model()

    if not fit_results['success']:
        print("\n❌ Validation failed - cannot fit model")
        return fit_results

    # Analyze by galaxy type
    validator.analyze_by_galaxy_type(fit_results)

    # Create plots
    validator.plot_results(fit_results)

    # Final summary
    print(f"\n{'='*80}")
    print("SESSION #23 SUMMARY")
    print(f"{'='*80}")

    if fit_results['r_squared'] > 0.406:
        print(f"\n✅ GALAXY-SPECIFIC ρ_mag MODEL SUCCESSFUL")
        improvement = ((fit_results['r_squared'] - 0.406) / 0.406) * 100
        print(f"   R² improved from 0.406 to {fit_results['r_squared']:.3f} (+{improvement:.1f}%)")
    else:
        print(f"\n⚠️ GALAXY-SPECIFIC MODEL DID NOT IMPROVE FIT")
        print(f"   R² = {fit_results['r_squared']:.3f} (Session #22: 0.406)")

    print(f"\nBest-fit parameters:")
    print(f"  ρ_sat,0 = {fit_results['rho_sat_0']:.2e} M_☉/pc³")
    print(f"  ρ_mag,0 = {fit_results['rho_mag_0']:.2e} M_☉/pc³")
    print(f"  B_ref   = {fit_results['B_ref']:.2f} μG")
    print(f"  α       = {fit_results['alpha']:.2f}  (B-field index)")
    print(f"  δ       = {fit_results['delta']:.2f}  (density index)")

    print(f"\nPhysical interpretation:")
    print(f"  ρ_mag ∝ B^(-{fit_results['alpha']:.2f})")
    print(f"  Strong B-fields → low ρ_mag → low ρ_sat ✓")

    return fit_results


if __name__ == "__main__":
    # Run galaxy-specific ρ_mag test
    results = test_galaxy_specific_rho_mag()

    print("\n✓ Session #23 galaxy-specific ρ_mag validation complete!")
    print("  Ready for Session #23 documentation")
