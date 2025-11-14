#!/usr/bin/env python3
"""
Synchronism REAL SPARC Database Validation - Session #16
==========================================================

THE CRITICAL TEST: Does Synchronism or NFW better match real observations?

Tests Synchronism dark matter predictions (γ=β=0.30, theory-predicted) against
175 REAL galaxy rotation curves from SPARC (Lelli et al. 2016).

This is THE decisive test between Synchronism and ΛCDM/NFW cosmology.

Key Difference from Session #15:
- Session #15: Synthetic data with NFW halos (circular test)
- Session #16: REAL observational data (no model assumptions!)

Author: Autonomous Research Agent (Claude Code)
Date: 2025-11-14
Session: #16 - Where theory meets reality
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import os
import glob


# Import from Session #15 validation code
import sys
sys.path.insert(0, os.path.dirname(__file__))
from synchronism_sparc_validation import (
    SPARCGalaxy,
    SynchronismPredictor,
    SPARCValidator
)


# ============================================================================
# Part 1: Real SPARC Data Loader
# ============================================================================

class RealSPARCLoader:
    """Load real SPARC galaxy data from Lelli et al. (2016)."""

    def __init__(self, data_dir: str = "sparc_real_data/galaxies"):
        """Initialize loader with SPARC data directory."""
        self.data_dir = data_dir

    def load_galaxy(self, filepath: str) -> Optional[SPARCGalaxy]:
        """
        Load one real SPARC galaxy from data file.

        SPARC format (columns):
        Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul

        Key differences from synthetic:
        - Vdisk, Vbul are velocities (not densities) for M/L=1
        - SBdisk, SBbul are surface brightness (not mass density)
        - Need to convert to mass densities for Synchronism
        """
        try:
            # Read data
            data = np.genfromtxt(
                filepath,
                comments='#',
                names=['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul',
                       'SBdisk', 'SBbul']
            )

            # Extract galaxy name from filename
            name = os.path.basename(filepath).replace('.dat', '')

            # Read distance from header (if available)
            distance = 10.0  # Default Mpc
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Distance:' in line:
                        try:
                            distance = float(line.split(':')[1].split()[0])
                        except:
                            pass
                        break

            # Convert surface brightness to surface density
            # SPARC provides SB in L_☉/pc² at 3.6μm
            # Assume M/L = 0.5 M_☉/L_☉ (typical for 3.6μm, stellar populations)
            M_L = 0.5

            sigma_disk = data['SBdisk'] * M_L  # M_☉/pc²
            sigma_bulge = data['SBbul'] * M_L if 'SBbul' in data.dtype.names else None

            # Convert gas velocity to surface density
            # V_gas² = G M_gas(<r) / r, for thin disk: M(<r) = 2π Σ r²
            # So: Σ_gas = V_gas² / (2π G r)
            G = 4.3e-6  # (km/s)²·kpc/M_☉
            sigma_gas = data['Vgas']**2 / (2.0 * np.pi * G * data['Rad']) / 1e6  # M_☉/pc²

            # Handle division by zero at r=0
            sigma_gas = np.where(data['Rad'] > 0.01, sigma_gas, 0.0)

            # Create galaxy object
            galaxy = SPARCGalaxy(
                name=name,
                distance=distance,
                inclination=45.0,  # Not critical for our analysis
                radius=data['Rad'],
                v_obs=data['Vobs'],
                v_err=data['errV'],
                sigma_disk=sigma_disk,
                sigma_bulge=sigma_bulge if sigma_bulge is not None else None,
                sigma_gas=sigma_gas,
                morphology="SPARC",
                total_luminosity=0.0  # Not needed for analysis
            )

            return galaxy

        except Exception as e:
            print(f"Error loading {filepath}: {e}")
            return None

    def load_all_galaxies(self, limit: Optional[int] = None) -> List[SPARCGalaxy]:
        """
        Load all SPARC galaxies from data directory.

        Args:
            limit: Optional limit on number of galaxies (for testing)

        Returns:
            List of SPARCGalaxy objects
        """
        # Find all .dat files
        pattern = os.path.join(self.data_dir, "*.dat")
        files = glob.glob(pattern)

        if limit:
            files = files[:limit]

        print(f"Loading {len(files)} galaxies from SPARC...")

        galaxies = []
        for filepath in files:
            galaxy = self.load_galaxy(filepath)
            if galaxy is not None:
                galaxies.append(galaxy)

        print(f"✓ Loaded {len(galaxies)} galaxies successfully")

        return galaxies


# ============================================================================
# Part 2: Main Validation
# ============================================================================

def main():
    """Session #16: The critical test - Synchronism vs reality!"""

    print("=" * 70)
    print("Session #16: Synchronism REAL SPARC Validation")
    print("=" * 70)
    print()
    print("THE CRITICAL TEST:")
    print("  Do real galaxies follow Synchronism or ΛCDM/NFW?")
    print()
    print("Theory-predicted parameters (Session #14, NO TUNING!):")
    print(f"  γ (coherence exponent) = {SynchronismPredictor.GAMMA:.2f}")
    print(f"  β (DM modulation exponent) = {SynchronismPredictor.BETA:.2f}")
    print()
    print("Free parameter:")
    print("  α (DM normalization) - fitted to each galaxy")
    print()
    print("Data source:")
    print("  SPARC: 175 real galaxies (Lelli et al. 2016)")
    print("  REAL observational data (no NFW assumptions!)")
    print()
    print("=" * 70)
    print()

    # Load real SPARC data
    loader = RealSPARCLoader()

    # Start with subset for initial testing
    print("Phase 1: Testing on 20 representative galaxies...")
    print()
    galaxies = loader.load_all_galaxies(limit=20)

    if len(galaxies) == 0:
        print("ERROR: No galaxies loaded!")
        return

    # Initialize Synchronism predictor with theory-predicted parameters
    predictor = SynchronismPredictor(
        gamma=SynchronismPredictor.GAMMA,  # 0.30 from Session #14
        beta=SynchronismPredictor.BETA      # 0.30 from Session #14
    )

    # Initialize validator
    validator = SPARCValidator(predictor)

    # Validate sample
    print("=" * 70)
    print("Running Synchronism validation on REAL galaxies...")
    print("=" * 70)
    print()

    stats = validator.validate_sample(galaxies)

    # Print results
    print()
    print("=" * 70)
    print("SESSION #16 RESULTS: Synchronism vs REAL Observational Data")
    print("=" * 70)
    print()
    print(f"Sample size: {stats['n_galaxies']} real galaxies")
    print()
    print("Goodness-of-fit:")
    print(f"  Mean χ²_red = {stats['mean_chi2_red']:.3f} ± {stats['std_chi2_red']:.3f}")
    print(f"  Median χ²_red = {stats['median_chi2_red']:.3f}")
    print(f"  Fraction with χ²_red < 2 (good fit): {stats['fraction_good_fit']*100:.1f}%")
    print(f"  Fraction with χ²_red < 3 (acceptable): {stats['fraction_acceptable_fit']*100:.1f}%")
    print()
    print("Dark matter normalization:")
    print(f"  Mean α = {stats['alpha_mean']:.3f} ± {stats['alpha_std']:.3f}")
    print(f"  Median α = {stats['alpha_median']:.3f}")
    print()
    print("Dark matter fraction:")
    print(f"  Mean M_DM/M_vis = {stats['M_DM_Mvis_mean']:.1f} ± {stats['M_DM_Mvis_std']:.1f}")
    print(f"  Expected range: 10-100 (observed from other methods)")
    print()

    # Compare to Session #15 (NFW synthetic)
    print("=" * 70)
    print("COMPARISON TO SESSION #15 (NFW Synthetic Data)")
    print("=" * 70)
    print()
    print("Session #15 (NFW synthetic):")
    print("  χ²_red = 277.5 ± 14.3")
    print("  Fraction good fit: 0.0%")
    print("  Result: Synchronism ≠ NFW (as expected)")
    print()
    print("Session #16 (REAL data):")
    print(f"  χ²_red = {stats['mean_chi2_red']:.1f} ± {stats['std_chi2_red']:.1f}")
    print(f"  Fraction good fit: {stats['fraction_good_fit']*100:.1f}%")
    print()

    # Interpret results
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print()

    if stats['mean_chi2_red'] < 2.0 and stats['fraction_good_fit'] > 0.6:
        print("✅ STRONG VALIDATION OF SYNCHRONISM!")
        print()
        print("Synchronism's theory-predicted γ = β = 0.30 produces")
        print("EXCELLENT fits to real galaxy rotation curves.")
        print()
        print("Implication: Dark matter = coherence (no exotic particles!)")
        print()
        print("This is a REVOLUTIONARY result if validated across full sample!")
        print()
        print("Next: Test on all 175 SPARC galaxies for statistical robustness.")

    elif stats['mean_chi2_red'] < 5.0 and stats['fraction_acceptable_fit'] > 0.4:
        print("⚠️  PARTIAL VALIDATION")
        print()
        print("Synchronism shows promise but not perfect fits.")
        print()
        print(f"Possible reasons (χ²_red = {stats['mean_chi2_red']:.1f}):")
        print("  1. γ, β need galaxy-type refinement")
        print("  2. Formula needs additional terms (∇ρ_vis?)")
        print("  3. M/L ratio assumptions (0.5 may vary)")
        print("  4. Observational systematics")
        print()
        print("Next: Analyze per-galaxy residuals, identify systematic patterns.")

    else:
        print("❌ SYNCHRONISM CHALLENGED")
        print()
        print("Current formula does not match real observations well.")
        print()
        print(f"With χ²_red = {stats['mean_chi2_red']:.1f}, Synchronism needs:")
        print("  1. Major formula revision, OR")
        print("  2. Different interpretation of Ξ^DM, OR")
        print("  3. Acknowledgment that DM ≠ coherence")
        print()
        print("However: Synchronism's EM and gravity predictions (Sessions #8-12)")
        print("remain valid independent of dark matter interpretation.")

    print()
    print("=" * 70)

    # Create visualizations
    print()
    print("Creating visualization plots...")
    validator.plot_sample_results(
        galaxies,
        filename="Session16_RealSPARC_Validation.png"
    )

    print()
    print("=" * 70)
    print("Session #16 Phase 1 complete!")
    print("=" * 70)
    print()
    print("Critical finding will be documented and pushed to repository.")
    print()


if __name__ == "__main__":
    main()
