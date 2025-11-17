#!/usr/bin/env python3
"""
Session #22: Dark Matter β(r) Radial Profile Test
==================================================

Tests Session #21 prediction that dark matter power-law index β varies
with radius:
- Inner regions (high ρ, C → C_max): β → 1.0
- Outer regions (low ρ, low C): β → 0.3

This tests whether ρ_DM = α(1-C)ρ_vis^β has β = β(r), not β = constant.

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-17
Session: #22 - Testing dark matter axiom predictions
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import os

# Import Session #20/22 infrastructure
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_session22_extract_data import load_session20_data
from synchronism_real_sparc_validation import RealSPARCLoader


# ============================================================================
# Part 1: Radial β(r) Measurement
# ============================================================================

class DarkMatterBetaTester:
    """
    Tests whether dark matter power-law index β varies with radius.

    Session #21 prediction:
    ρ_DM = α(1-C) ρ_vis^β(r)

    where β(r) transitions from ~0.3 (outer) to ~1.0 (inner)
    """

    def __init__(self):
        """Initialize tester."""
        self.loader = RealSPARCLoader()

    def measure_beta_profile(self, galaxy, n_bins: int = 5):
        """
        Measure β(r) by fitting ρ_DM vs ρ_vis in radial bins.

        Args:
            galaxy: SPARC galaxy object
            n_bins: Number of radial bins

        Returns:
            (r_bins, beta_bins, beta_errors)
        """
        # Get rotation curve components
        r = galaxy.radius  # kpc
        v_obs = galaxy.v_obs
        v_disk = galaxy.v_disk
        v_gas = galaxy.v_gas if hasattr(galaxy, 'v_gas') else np.zeros_like(r)
        v_bulge = galaxy.v_bulge if hasattr(galaxy, 'v_bulge') else np.zeros_like(r)

        # Visible matter velocity
        v_vis = np.sqrt(v_disk**2 + v_gas**2 + v_bulge**2)

        # Dark matter velocity (from residual)
        v_dm_squared = np.maximum(0, v_obs**2 - v_vis**2)
        v_dm = np.sqrt(v_dm_squared)

        # Convert to densities (approximate)
        # ρ ∝ v² / r² (for NFW-like profiles)
        rho_vis = v_vis**2 / (r**2 + 1e-3)  # Add small offset to avoid r=0
        rho_dm = v_dm**2 / (r**2 + 1e-3)

        # Only keep points with significant dark matter
        mask = (v_dm > 0.1 * v_obs) & (rho_vis > 0)
        if np.sum(mask) < n_bins:
            return None, None, None

        r = r[mask]
        rho_vis = rho_vis[mask]
        rho_dm = rho_dm[mask]

        # Bin by radius
        r_edges = np.percentile(r, np.linspace(0, 100, n_bins + 1))
        r_bins = []
        beta_bins = []
        beta_errors = []

        for i in range(n_bins):
            # Select points in this radial bin
            bin_mask = (r >= r_edges[i]) & (r < r_edges[i+1])
            if np.sum(bin_mask) < 3:
                continue

            r_bin = r[bin_mask]
            rho_vis_bin = rho_vis[bin_mask]
            rho_dm_bin = rho_dm[bin_mask]

            # Fit ρ_DM = A * ρ_vis^β in this bin
            def power_law(rho_v, A, beta):
                return A * rho_v**beta

            try:
                popt, pcov = curve_fit(
                    power_law,
                    rho_vis_bin,
                    rho_dm_bin,
                    p0=[1.0, 0.5],
                    bounds=([1e-3, 0.01], [1e3, 2.0]),
                    maxfev=5000
                )

                A_best, beta_best = popt
                beta_err = np.sqrt(pcov[1, 1])

                r_bins.append(np.median(r_bin))
                beta_bins.append(beta_best)
                beta_errors.append(beta_err)

            except:
                continue

        if len(r_bins) == 0:
            return None, None, None

        return np.array(r_bins), np.array(beta_bins), np.array(beta_errors)

    def test_beta_variation(self, n_galaxies: int = 20, n_bins: int = 4):
        """
        Test β(r) variation across sample of galaxies.

        Args:
            n_galaxies: Number of galaxies to analyze
            n_bins: Number of radial bins per galaxy

        Returns:
            Results dictionary
        """
        print("="*80)
        print("SESSION #22: DARK MATTER β(r) RADIAL PROFILE TEST")
        print("="*80)

        print(f"\nSession #21 prediction:")
        print(f"  Inner regions (high ρ): β → 1.0")
        print(f"  Outer regions (low ρ): β → 0.3")

        print(f"\nTesting {n_galaxies} galaxies with {n_bins} radial bins each...")

        # Load galaxies
        galaxies = self.loader.load_all_galaxies(limit=None)

        # Select diverse sample (different types and sizes)
        np.random.seed(42)
        sample_indices = np.random.choice(len(galaxies), min(n_galaxies, len(galaxies)), replace=False)
        sample = [galaxies[i] for i in sample_indices]

        # Measure β(r) for each galaxy
        all_r_bins = []
        all_beta_bins = []
        all_beta_errors = []
        galaxy_names = []

        for i, galaxy in enumerate(sample):
            print(f"[{i+1}/{len(sample)}] Measuring β(r) for {galaxy.name}...")

            r_bins, beta_bins, beta_errors = self.measure_beta_profile(galaxy, n_bins=n_bins)

            if r_bins is None:
                print(f"  Skipped (insufficient data)")
                continue

            print(f"  ✓ {len(r_bins)} radial bins")

            all_r_bins.extend(r_bins)
            all_beta_bins.extend(beta_bins)
            all_beta_errors.extend(beta_errors)
            galaxy_names.extend([galaxy.name] * len(r_bins))

        print(f"\n✓ Measured {len(all_r_bins)} radial bins across {len(set(galaxy_names))} galaxies")

        # Convert to arrays
        all_r_bins = np.array(all_r_bins)
        all_beta_bins = np.array(all_beta_bins)
        all_beta_errors = np.array(all_beta_errors)

        # Analyze radial trend
        print(f"\n{'='*80}")
        print("RADIAL TREND ANALYSIS")
        print(f"{'='*80}")

        # Split into inner and outer regions
        r_median = np.median(all_r_bins)
        inner_mask = all_r_bins < r_median
        outer_mask = all_r_bins >= r_median

        beta_inner = all_beta_bins[inner_mask]
        beta_outer = all_beta_bins[outer_mask]

        print(f"\nInner regions (r < {r_median:.2f} kpc):")
        print(f"  N = {len(beta_inner)}")
        print(f"  β_median = {np.median(beta_inner):.2f}")
        print(f"  β_mean = {np.mean(beta_inner):.2f} ± {np.std(beta_inner):.2f}")
        print(f"  Session #21 prediction: β → 1.0")

        print(f"\nOuter regions (r >= {r_median:.2f} kpc):")
        print(f"  N = {len(beta_outer)}")
        print(f"  β_median = {np.median(beta_outer):.2f}")
        print(f"  β_mean = {np.mean(beta_outer):.2f} ± {np.std(beta_outer):.2f}")
        print(f"  Session #21 prediction: β → 0.3")

        # Test if β increases inward
        if np.median(beta_inner) > np.median(beta_outer):
            delta_beta = np.median(beta_inner) - np.median(beta_outer)
            print(f"\n✅ TREND CONFIRMED: β increases inward")
            print(f"   Δβ = {delta_beta:+.2f} (inner - outer)")
        else:
            print(f"\n❌ Trend NOT confirmed: β does not increase inward")

        # Correlation test
        if len(all_r_bins) > 10:
            corr = np.corrcoef(all_r_bins, all_beta_bins)[0, 1]
            print(f"\nCorrelation(r, β): {corr:.3f}")

            if corr < -0.2:
                print(f"  ✅ Negative correlation: β decreases with radius")
            elif abs(corr) < 0.2:
                print(f"  ⚠️ Weak correlation: No clear radial trend")
            else:
                print(f"  ❌ Positive correlation: β increases with radius (unexpected)")

        return {
            'r_bins': all_r_bins,
            'beta_bins': all_beta_bins,
            'beta_errors': all_beta_errors,
            'galaxy_names': galaxy_names,
            'beta_inner_median': np.median(beta_inner),
            'beta_outer_median': np.median(beta_outer),
            'r_median': r_median,
            'correlation': corr if len(all_r_bins) > 10 else np.nan
        }

    def plot_beta_profile(self, results: dict, output_file: str = "session22_beta_profile.png"):
        """
        Plot β(r) results.

        Args:
            results: Output from test_beta_variation()
            output_file: Output filename
        """
        print(f"\nCreating β(r) profile plot: {output_file}")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: β vs r (scatter)
        ax1.errorbar(
            results['r_bins'],
            results['beta_bins'],
            yerr=results['beta_errors'],
            fmt='o', alpha=0.5, ms=6, elinewidth=1, capsize=2
        )

        # Session #21 predictions
        ax1.axhline(1.0, ls='--', c='red', lw=2, label='Inner prediction (β → 1.0)')
        ax1.axhline(0.3, ls='--', c='blue', lw=2, label='Outer prediction (β → 0.3)')

        # Median trend (binned)
        r_edges = np.percentile(results['r_bins'], [0, 33, 67, 100])
        for i in range(len(r_edges)-1):
            mask = (results['r_bins'] >= r_edges[i]) & (results['r_bins'] < r_edges[i+1])
            if np.sum(mask) > 0:
                r_mean = np.mean(results['r_bins'][mask])
                beta_mean = np.mean(results['beta_bins'][mask])
                beta_std = np.std(results['beta_bins'][mask])
                ax1.errorbar([r_mean], [beta_mean], yerr=[beta_std],
                           fmt='s', color='red', ms=10, elinewidth=2, capsize=4)

        ax1.set_xlabel("Radius r [kpc]", fontsize=13)
        ax1.set_ylabel("Power-law Index β", fontsize=13)
        ax1.set_title(f"Dark Matter β(r) Profile (n={len(results['r_bins'])} bins)", fontsize=14, fontweight='bold')
        ax1.legend(fontsize=11)
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-0.2, 2.0)

        # Plot 2: β distribution (inner vs outer)
        r_median = results['r_median']
        inner_mask = results['r_bins'] < r_median
        outer_mask = results['r_bins'] >= r_median

        ax2.hist(results['beta_bins'][inner_mask], bins=15, alpha=0.6, color='red',
                label=f'Inner (r < {r_median:.1f} kpc)', density=True)
        ax2.hist(results['beta_bins'][outer_mask], bins=15, alpha=0.6, color='blue',
                label=f'Outer (r ≥ {r_median:.1f} kpc)', density=True)

        # Session #21 predictions
        ax2.axvline(1.0, ls='--', c='red', lw=2)
        ax2.axvline(0.3, ls='--', c='blue', lw=2)

        ax2.set_xlabel("Power-law Index β", fontsize=13)
        ax2.set_ylabel("Probability Density", fontsize=13)
        ax2.set_title("β Distribution: Inner vs Outer", fontsize=14, fontweight='bold')
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ β(r) profile plot saved: {output_file}")

        plt.close()


# ============================================================================
# Part 2: Main Execution
# ============================================================================

def test_dark_matter_beta_profile():
    """
    Main test: Measure β(r) and test Session #21 prediction.
    """
    print("="*80)
    print("SESSION #22: DARK MATTER β(r) PREDICTION TEST")
    print("Testing Session #21 axiom-derived dark matter formula")
    print("="*80)

    # Initialize tester
    tester = DarkMatterBetaTester()

    # Test β(r) variation
    results = tester.test_beta_variation(n_galaxies=30, n_bins=5)

    # Plot results
    tester.plot_beta_profile(results)

    # Summary
    print(f"\n{'='*80}")
    print("SESSION #22 DARK MATTER β(r) TEST SUMMARY")
    print(f"{'='*80}")

    print(f"\nSession #21 prediction: β varies from ~0.3 (outer) to ~1.0 (inner)")

    print(f"\nMeasured:")
    print(f"  Inner β_median = {results['beta_inner_median']:.2f}")
    print(f"  Outer β_median = {results['beta_outer_median']:.2f}")
    print(f"  Δβ = {results['beta_inner_median'] - results['beta_outer_median']:+.2f}")
    print(f"  Correlation(r, β) = {results['correlation']:.3f}")

    if results['beta_inner_median'] > results['beta_outer_median']:
        print(f"\n✅ PREDICTION PARTIALLY CONFIRMED")
        print(f"   β increases inward as predicted")
    else:
        print(f"\n⚠️ Prediction requires further investigation")

    print(f"\nNote: This is a preliminary test. Full validation requires:")
    print(f"  - Better dark matter extraction (deconvolution)")
    print(f"  - Coherence C(r) measurements")
    print(f"  - Individual galaxy analysis")

    return results


if __name__ == "__main__":
    # Run dark matter β(r) test
    results = test_dark_matter_beta_profile()

    print("\n✓ Session #22 dark matter β(r) test complete!")
    print("  Ready for final documentation")
