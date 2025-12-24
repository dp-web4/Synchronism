#!/usr/bin/env python3
"""
SESSION #176: CLUSTER VELOCITY DISPERSION PROFILES - THEORETICAL FRAMEWORK
===========================================================================
Date: December 24, 2025

FOLLOW-UP TO SESSION #175:
--------------------------
Session #175 concluded the Real Data Application arc with the finding that
peculiar velocities are the WRONG observable for testing the coherence function
due to MRH mismatch (kpc vs Mpc scales).

However, Session #173 found a promising result: cluster velocity dispersion
ratio (outer/inner) = 3.28, in the Synchronism direction.

This session develops the THEORETICAL FRAMEWORK for cluster velocity dispersion
profiles under Synchronism vs ΛCDM.

KEY INSIGHT:
------------
Cluster internal dynamics DOES probe the correct MRH:
- Cluster virial radius: R_vir ~ 1-3 Mpc
- Velocity dispersion profile: r ~ 0.1 - 10 Mpc
- This overlaps with coherence function transition scale

PREDICTION:
-----------
Synchronism: σ_v(r) should INCREASE with radius (enhanced G_eff in outskirts)
ΛCDM: σ_v(r) should DECREASE with radius (standard virial equilibrium)

This is a DISCRIMINATING prediction.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq

print("=" * 70)
print("SESSION #176: CLUSTER VELOCITY DISPERSION PROFILES")
print("Theoretical Framework for Synchronism vs ΛCDM")
print("=" * 70)

# =============================================================================
# 1. SYNCHRONISM COHERENCE FUNCTION
# =============================================================================

print("\n" + "=" * 70)
print("1. SYNCHRONISM COHERENCE FUNCTION")
print("=" * 70)

phi = (1 + np.sqrt(5)) / 2  # Golden ratio
Omega_m = 0.3

def coherence(rho_ratio):
    """
    Coherence function C(ρ) from Synchronism whitepaper.

    C(ρ) = Ω_m + (1 - Ω_m) × (ρ/ρ_t)^(1/φ) / [1 + (ρ/ρ_t)^(1/φ)]

    where:
    - ρ_t = 1 (transition density in units of cosmic mean)
    - φ = golden ratio
    - Ω_m = 0.3

    Returns: C between Ω_m (low density) and 1 (high density)
    """
    if np.isscalar(rho_ratio):
        if rho_ratio <= 0:
            return Omega_m
        x = rho_ratio ** (1/phi)
        return Omega_m + (1 - Omega_m) * x / (1 + x)
    else:
        result = np.full_like(rho_ratio, Omega_m)
        pos = rho_ratio > 0
        x = np.zeros_like(rho_ratio)
        x[pos] = rho_ratio[pos] ** (1/phi)
        result[pos] = Omega_m + (1 - Omega_m) * x[pos] / (1 + x[pos])
        return result

def G_eff_ratio(rho_ratio):
    """
    Effective gravitational constant ratio.

    G_eff/G = 1/C(ρ)

    In low-density regions: G_eff > G (enhanced gravity)
    In high-density regions: G_eff ≈ G (standard gravity)
    """
    return 1.0 / coherence(rho_ratio)

# Demonstrate coherence function
print("\nCoherence function behavior:")
print("-" * 50)
densities = [0.01, 0.1, 1.0, 10, 100, 1000]
for rho in densities:
    C = coherence(rho)
    G_ratio = G_eff_ratio(rho)
    print(f"  ρ/ρ̄ = {rho:>6.2f}: C = {C:.4f}, G_eff/G = {G_ratio:.4f}")

# =============================================================================
# 2. CLUSTER DENSITY PROFILE
# =============================================================================

print("\n" + "=" * 70)
print("2. CLUSTER DENSITY PROFILE (NFW)")
print("=" * 70)

def NFW_profile(r, r_s, rho_s):
    """
    Navarro-Frenk-White (NFW) density profile.

    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]

    where:
    - r_s = scale radius
    - ρ_s = characteristic density
    """
    x = r / r_s
    return rho_s / (x * (1 + x)**2)

def cluster_density_profile(r, M_200, c_200=4.0):
    """
    NFW density profile for a cluster with given virial mass and concentration.

    Parameters:
    - M_200: virial mass in M_sun
    - c_200: concentration parameter (typical: 4-6 for clusters)

    Returns density in M_sun/Mpc³
    """
    # Critical density at z=0
    H0 = 70  # km/s/Mpc
    G = 4.302e-9  # Mpc³/(M_sun·Gyr²)
    rho_crit = 3 * (H0/977.8)**2 / (8 * np.pi * G)  # M_sun/Mpc³

    # R_200 (virial radius) from M_200
    R_200 = (3 * M_200 / (4 * np.pi * 200 * rho_crit))**(1/3)  # Mpc

    # Scale radius
    r_s = R_200 / c_200

    # Characteristic density
    delta_c = (200/3) * c_200**3 / (np.log(1 + c_200) - c_200/(1 + c_200))
    rho_s = delta_c * rho_crit

    return NFW_profile(r, r_s, rho_s), R_200, r_s

# Example: Coma-like cluster (M ~ 10^15 M_sun)
M_200 = 1e15  # M_sun
c_200 = 4.0

radii = np.logspace(-2, 1, 100)  # 0.01 to 10 Mpc

density, R_200, r_s = cluster_density_profile(radii, M_200, c_200)

print(f"\nExample cluster (M_200 = {M_200:.0e} M_sun):")
print(f"  Virial radius R_200 = {R_200:.2f} Mpc")
print(f"  Scale radius r_s = {r_s:.2f} Mpc")
print(f"  Concentration c = {c_200}")

# Cosmic mean density
rho_cosmic = 2.77e11 * Omega_m  # M_sun/Mpc³

print(f"\nDensity profile (relative to cosmic mean):")
sample_r = [0.1, 0.3, 1.0, R_200, 3*R_200, 10*R_200]
for r in sample_r:
    rho, _, _ = cluster_density_profile(r, M_200, c_200)
    rho_ratio = rho / rho_cosmic
    print(f"  r = {r:.2f} Mpc: ρ/ρ̄ = {rho_ratio:.1f}")

# =============================================================================
# 3. VELOCITY DISPERSION IN ΛCDM
# =============================================================================

print("\n" + "=" * 70)
print("3. ΛCDM VELOCITY DISPERSION PROFILE")
print("=" * 70)

def jeans_equation_sigma_LCDM(r, M_200, c_200=4.0, beta=0.0):
    """
    Velocity dispersion from spherical Jeans equation in ΛCDM.

    For an NFW profile with constant anisotropy β:

    σ_r²(r) = (1/ρ(r)) ∫_r^∞ ρ(s) × (G M(<s)/s²) × (s/r)^(2β) ds

    Simplified for β=0 (isotropic):

    σ²(r) = (G/ρ(r)) ∫_r^∞ ρ(s) × M(<s)/s² ds
    """
    # Get profile parameters
    _, R_200, r_s = cluster_density_profile(1.0, M_200, c_200)

    # Critical density
    H0 = 70
    G = 4.302e-9
    rho_crit = 3 * (H0/977.8)**2 / (8 * np.pi * G)

    # NFW characteristic density
    delta_c = (200/3) * c_200**3 / (np.log(1 + c_200) - c_200/(1 + c_200))
    rho_s = delta_c * rho_crit

    # Enclosed mass for NFW: M(<r) = 4π ρ_s r_s³ × [ln(1 + r/r_s) - (r/r_s)/(1 + r/r_s)]
    def M_enclosed(s):
        x = s / r_s
        return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x/(1 + x))

    def rho_NFW(s):
        x = s / r_s
        return rho_s / (x * (1 + x)**2)

    # Integrand for Jeans equation (isotropic case)
    def integrand(s):
        return rho_NFW(s) * G * M_enclosed(s) / s**2

    # Numerical integration
    rho_r = rho_NFW(r)

    # Integrate from r to large radius (100*R_200 as proxy for infinity)
    r_max = 100 * R_200
    integral, _ = quad(integrand, r, r_max, limit=100)

    sigma_r_squared = integral / rho_r

    # Convert to 1D velocity dispersion (km/s)
    # σ_1D ≈ σ_r for isotropic case
    sigma_1D = np.sqrt(sigma_r_squared) * 977.8  # Convert Mpc/Gyr to km/s

    return sigma_1D

# Calculate ΛCDM profile
print(f"\nΛCDM velocity dispersion profile (Jeans equation):")
print("-" * 50)

sigma_LCDM = []
radii_sample = np.array([0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0])
for r in radii_sample:
    sigma = jeans_equation_sigma_LCDM(r, M_200, c_200)
    sigma_LCDM.append(sigma)
    print(f"  r = {r:.1f} Mpc (r/R_200 = {r/R_200:.2f}): σ_v = {sigma:.0f} km/s")

sigma_LCDM = np.array(sigma_LCDM)

# =============================================================================
# 4. VELOCITY DISPERSION IN SYNCHRONISM
# =============================================================================

print("\n" + "=" * 70)
print("4. SYNCHRONISM VELOCITY DISPERSION PROFILE")
print("=" * 70)

def jeans_equation_sigma_Synchronism(r, M_200, c_200=4.0, beta=0.0):
    """
    Velocity dispersion from spherical Jeans equation in Synchronism.

    Key modification: Replace G with G_eff(ρ) = G/C(ρ)

    σ_r²(r) = (1/ρ(r)) ∫_r^∞ ρ(s) × (G_eff(s) M(<s)/s²) × (s/r)^(2β) ds
    """
    # Get profile parameters
    _, R_200, r_s = cluster_density_profile(1.0, M_200, c_200)

    # Critical density
    H0 = 70
    G = 4.302e-9
    rho_crit = 3 * (H0/977.8)**2 / (8 * np.pi * G)

    # Cosmic mean density
    rho_cosmic = rho_crit * Omega_m

    # NFW characteristic density
    delta_c = (200/3) * c_200**3 / (np.log(1 + c_200) - c_200/(1 + c_200))
    rho_s = delta_c * rho_crit

    # Enclosed mass for NFW
    def M_enclosed(s):
        x = s / r_s
        return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x/(1 + x))

    def rho_NFW(s):
        x = s / r_s
        return rho_s / (x * (1 + x)**2)

    # Integrand with G_eff
    def integrand(s):
        rho = rho_NFW(s)
        rho_ratio = rho / rho_cosmic
        G_eff = G * G_eff_ratio(rho_ratio)  # G_eff = G/C(ρ)
        return rho * G_eff * M_enclosed(s) / s**2

    # Numerical integration
    rho_r = rho_NFW(r)

    r_max = 100 * R_200
    integral, _ = quad(integrand, r, r_max, limit=100)

    sigma_r_squared = integral / rho_r
    sigma_1D = np.sqrt(sigma_r_squared) * 977.8

    return sigma_1D

# Calculate Synchronism profile
print(f"\nSynchronism velocity dispersion profile (modified Jeans equation):")
print("-" * 50)

sigma_Sync = []
for r in radii_sample:
    sigma = jeans_equation_sigma_Synchronism(r, M_200, c_200)
    sigma_Sync.append(sigma)

    # Also calculate G_eff at this radius
    rho, _, _ = cluster_density_profile(r, M_200, c_200)
    rho_ratio = rho / rho_cosmic
    G_ratio = G_eff_ratio(rho_ratio)

    print(f"  r = {r:.1f} Mpc: σ_v = {sigma:.0f} km/s (G_eff/G = {G_ratio:.3f})")

sigma_Sync = np.array(sigma_Sync)

# =============================================================================
# 5. COMPARISON AND PREDICTIONS
# =============================================================================

print("\n" + "=" * 70)
print("5. ΛCDM vs SYNCHRONISM COMPARISON")
print("=" * 70)

print(f"\nVelocity dispersion ratio (Synchronism/ΛCDM):")
print("-" * 50)
for i, r in enumerate(radii_sample):
    ratio = sigma_Sync[i] / sigma_LCDM[i]
    print(f"  r = {r:.1f} Mpc: σ_Sync/σ_ΛCDM = {ratio:.3f}")

# Key prediction: Outer/Inner ratio
inner_r = 0.2  # Mpc
outer_r = 3.0  # Mpc

sigma_inner_LCDM = jeans_equation_sigma_LCDM(inner_r, M_200, c_200)
sigma_outer_LCDM = jeans_equation_sigma_LCDM(outer_r, M_200, c_200)
sigma_inner_Sync = jeans_equation_sigma_Synchronism(inner_r, M_200, c_200)
sigma_outer_Sync = jeans_equation_sigma_Synchronism(outer_r, M_200, c_200)

ratio_LCDM = sigma_outer_LCDM / sigma_inner_LCDM
ratio_Sync = sigma_outer_Sync / sigma_inner_Sync

print(f"\n*** KEY DISCRIMINATING PREDICTION ***")
print(f"\nOuter (r={outer_r} Mpc) / Inner (r={inner_r} Mpc) dispersion ratio:")
print(f"  ΛCDM:       σ_outer/σ_inner = {ratio_LCDM:.3f}")
print(f"  Synchronism: σ_outer/σ_inner = {ratio_Sync:.3f}")
print(f"  Difference: {(ratio_Sync - ratio_LCDM)/ratio_LCDM * 100:.1f}%")

# Session #173 observed ratio
observed_ratio = 3.28
print(f"\n  Session #173 observed (CF4 clusters): {observed_ratio:.2f}")
print(f"  Closer to: {'Synchronism' if abs(observed_ratio - ratio_Sync) < abs(observed_ratio - ratio_LCDM) else 'ΛCDM'}")

# =============================================================================
# 6. RADIAL PROFILE BEHAVIOR
# =============================================================================

print("\n" + "=" * 70)
print("6. RADIAL PROFILE BEHAVIOR")
print("=" * 70)

# Compute gradient
inner_slope_LCDM = (sigma_LCDM[1] - sigma_LCDM[0]) / (radii_sample[1] - radii_sample[0])
inner_slope_Sync = (sigma_Sync[1] - sigma_Sync[0]) / (radii_sample[1] - radii_sample[0])

print(f"\nInner gradient (dσ/dr at r ~ 0.1-0.2 Mpc):")
print(f"  ΛCDM:       {inner_slope_LCDM:.1f} km/s/Mpc")
print(f"  Synchronism: {inner_slope_Sync:.1f} km/s/Mpc")

outer_slope_LCDM = (sigma_LCDM[-1] - sigma_LCDM[-2]) / (radii_sample[-1] - radii_sample[-2])
outer_slope_Sync = (sigma_Sync[-1] - sigma_Sync[-2]) / (radii_sample[-1] - radii_sample[-2])

print(f"\nOuter gradient (dσ/dr at r ~ 3-5 Mpc):")
print(f"  ΛCDM:       {outer_slope_LCDM:.1f} km/s/Mpc")
print(f"  Synchronism: {outer_slope_Sync:.1f} km/s/Mpc")

# Key qualitative difference
print(f"\n*** QUALITATIVE PREDICTION ***")
print(f"  ΛCDM: σ(r) decreases monotonically outward (standard virial)")
print(f"  Synchronism: σ(r) enhanced in outskirts (may show upturn or slower decline)")

# =============================================================================
# 7. VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("7. GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Density profile with coherence regions
ax1 = axes[0, 0]
r_plot = np.logspace(-2, 1.3, 200)
rho_plot, _, _ = cluster_density_profile(r_plot, M_200, c_200)
rho_ratio_plot = rho_plot / rho_cosmic

ax1.loglog(r_plot, rho_ratio_plot, 'k-', linewidth=2, label='NFW profile')
ax1.axhline(1.0, color='orange', linestyle='--', label='ρ = ρ_cosmic (C transition)')
ax1.axhline(10.0, color='green', linestyle=':', alpha=0.7, label='ρ = 10 ρ_cosmic')
ax1.axvline(R_200, color='red', linestyle='--', alpha=0.7, label=f'R_200 = {R_200:.1f} Mpc')
ax1.fill_between(r_plot, 0.01, rho_ratio_plot, where=(rho_ratio_plot > 1),
                  alpha=0.2, color='blue', label='C → 1 (ΛCDM regime)')
ax1.fill_between(r_plot, 0.01, rho_ratio_plot, where=(rho_ratio_plot < 1),
                  alpha=0.2, color='orange', label='C < 1 (enhanced G_eff)')
ax1.set_xlabel('Radius (Mpc)')
ax1.set_ylabel('ρ / ρ_cosmic')
ax1.set_title(f'Cluster Density Profile (M_200 = {M_200:.0e} M☉)')
ax1.legend(loc='upper right', fontsize=8)
ax1.set_ylim(0.01, 1e6)
ax1.set_xlim(0.01, 20)
ax1.grid(True, alpha=0.3)

# Panel 2: G_eff profile
ax2 = axes[0, 1]
G_eff_plot = G_eff_ratio(rho_ratio_plot)
C_plot = coherence(rho_ratio_plot)

ax2.semilogx(r_plot, G_eff_plot, 'b-', linewidth=2, label='G_eff/G')
ax2.semilogx(r_plot, C_plot, 'r--', linewidth=2, label='C(ρ)')
ax2.axhline(1.0, color='gray', linestyle=':')
ax2.axhline(1/Omega_m, color='orange', linestyle='--', alpha=0.5,
            label=f'G_eff/G max = {1/Omega_m:.2f}')
ax2.axvline(R_200, color='red', linestyle='--', alpha=0.7)
ax2.set_xlabel('Radius (Mpc)')
ax2.set_ylabel('Ratio')
ax2.set_title('Coherence and Effective Gravity')
ax2.legend(loc='upper right')
ax2.set_ylim(0, 4)
ax2.grid(True, alpha=0.3)

# Panel 3: Velocity dispersion profiles
ax3 = axes[1, 0]
ax3.plot(radii_sample, sigma_LCDM, 'b-o', linewidth=2, markersize=8, label='ΛCDM')
ax3.plot(radii_sample, sigma_Sync, 'r-s', linewidth=2, markersize=8, label='Synchronism')
ax3.axvline(R_200, color='gray', linestyle='--', alpha=0.7, label=f'R_200')
ax3.set_xlabel('Radius (Mpc)')
ax3.set_ylabel('Velocity dispersion σ_v (km/s)')
ax3.set_title('Velocity Dispersion Profiles')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 6)

# Panel 4: Ratio of Synchronism to ΛCDM
ax4 = axes[1, 1]
ratio_plot = sigma_Sync / sigma_LCDM
ax4.plot(radii_sample, ratio_plot, 'g-o', linewidth=2, markersize=8)
ax4.axhline(1.0, color='gray', linestyle='--')
ax4.axvline(R_200, color='red', linestyle='--', alpha=0.7, label=f'R_200')
ax4.set_xlabel('Radius (Mpc)')
ax4.set_ylabel('σ_Sync / σ_ΛCDM')
ax4.set_title('Synchronism Enhancement Factor')
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_ylim(0.8, 2.0)
ax4.set_xlim(0, 6)

plt.suptitle('Session #176: Cluster Velocity Dispersion Profiles\nSynchronism vs ΛCDM Theoretical Predictions',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session176_theory.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session176_theory.png")

# =============================================================================
# 8. OBSERVABLE PREDICTIONS
# =============================================================================

print("\n" + "=" * 70)
print("8. OBSERVATIONAL PREDICTIONS")
print("=" * 70)

print("""
TESTABLE PREDICTIONS FOR CLUSTER VELOCITY DISPERSIONS
======================================================

1. RADIAL PROFILE SHAPE:
   - ΛCDM: σ(r) decreases monotonically ∝ r^(-0.1 to -0.3)
   - Synchronism: σ(r) shows shallower decline or slight upturn in outskirts

2. OUTER/INNER RATIO:
   - ΛCDM predicts: σ(3 Mpc) / σ(0.2 Mpc) ≈ {:.2f}
   - Synchronism predicts: σ(3 Mpc) / σ(0.2 Mpc) ≈ {:.2f}
   - Session #173 observed: 3.28 (from CF4 cluster candidates)

3. MASS-DEPENDENT EFFECT:
   - Lower mass clusters have steeper density gradients
   - Synchronism effect should be STRONGER for low-mass systems

4. REDSHIFT EVOLUTION:
   - Cosmic mean density higher at high z
   - Transition to enhanced G_eff occurs at higher cluster density
   - Synchronism effect should be WEAKER at high redshift

5. DATA SOURCES FOR TESTING:
   - Planck PSZ2 catalog (388 clusters with σ_v and M_dyn)
   - DES redMaPPer clusters
   - SDSS spectroscopic cluster members
   - Individual massive clusters (Coma, A2029, A426)

6. EXPECTED SIGNATURE:
   - In stacked cluster profiles: excess dispersion at r > R_200
   - In dynamical mass estimates: M_dyn > M_lensing at large radii
   - In cluster infall regions: higher peculiar velocities than ΛCDM predicts
""".format(ratio_LCDM, ratio_Sync))

# =============================================================================
# 9. COMPARISON WITH SESSION #173 RESULT
# =============================================================================

print("\n" + "=" * 70)
print("9. COMPARISON WITH SESSION #173 RESULT")
print("=" * 70)

print(f"""
SESSION #173 CLUSTER INFALL TEST:
=================================

Methodology:
- Identified cluster candidates from CF4 overdensities (top 1%)
- Computed radial velocities relative to cluster systemic velocity
- Compared inner (<5 Mpc) vs outer (20-60 Mpc) regions

Result:
- Inner dispersion: σ_inner ≈ 900 km/s
- Outer dispersion: σ_outer ≈ 2950 km/s
- Ratio (outer/inner): 3.28

This session's theoretical prediction:
- ΛCDM ratio: {ratio_LCDM:.2f}
- Synchronism ratio: {ratio_Sync:.2f}

INTERPRETATION:
---------------
The Session #173 observed ratio (3.28) is MUCH higher than either prediction.

However, Session #173 used DIFFERENT radii than this theoretical calculation:
- Session #173: 5 Mpc (inner) vs 20-60 Mpc (outer)
- This session: 0.2 Mpc (inner) vs 3 Mpc (outer)

At 20-60 Mpc, we're in the COSMIC FIELD, not cluster outskirts.
The density there is ~1× cosmic mean, so G_eff/G → {1/Omega_m:.2f}

This could explain the high observed ratio!

Let me recalculate with Session #173's radii...
""")

# Recalculate with Session #173 radii
inner_173 = 5.0  # Mpc
outer_173 = 40.0  # Mpc (middle of 20-60 range)

sigma_inner_173_LCDM = jeans_equation_sigma_LCDM(inner_173, M_200, c_200)
sigma_inner_173_Sync = jeans_equation_sigma_Synchronism(inner_173, M_200, c_200)

# For outer region at 40 Mpc, cluster potential is negligible
# Use field velocity dispersion from cosmic structure
# σ_field ≈ 300-600 km/s from LSS
sigma_field = 400  # km/s (approximate field velocity)

print(f"\nRecalculation with Session #173 radii:")
print(f"  Inner (5 Mpc):")
print(f"    ΛCDM: σ = {sigma_inner_173_LCDM:.0f} km/s")
print(f"    Synchronism: σ = {sigma_inner_173_Sync:.0f} km/s")
print(f"  Outer (20-60 Mpc):")
print(f"    This is FIELD, not cluster!")
print(f"    Typical field σ ≈ {sigma_field} km/s")

print(f"""

CRITICAL INSIGHT:
-----------------
Session #173's "outer region" (20-60 Mpc) is NOT cluster outskirts.
It's the COSMIC FIELD where different physics applies:
- No longer bound to cluster
- Velocity dispersion from large-scale structure
- Hubble flow + peculiar velocities

The ratio of 3.28 is comparing:
- Cluster-bound dynamics (inner)
- Field peculiar velocities (outer)

These are DIFFERENT dynamical regimes!

RECOMMENDATION:
---------------
For proper Synchronism test, need velocity dispersion profiles
WITHIN cluster virial radius (r < 3 Mpc):
- Use spectroscopic cluster surveys
- Stack similar-mass clusters
- Measure σ(r) from 0.1 to 3 R_200
""")

# =============================================================================
# 10. SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #176: SUMMARY")
print("=" * 70)

print(f"""
THEORETICAL FRAMEWORK FOR CLUSTER VELOCITY DISPERSIONS
======================================================

1. COHERENCE FUNCTION APPLICATION:
   - Applied C(ρ) = Ω_m + (1-Ω_m) × (ρ/ρ_t)^(1/φ) / [1 + (ρ/ρ_t)^(1/φ)]
   - G_eff = G/C(ρ) → enhanced in low-density regions

2. CLUSTER DENSITY PROFILE:
   - Used NFW profile for Coma-like cluster (M = 10^15 M_sun)
   - ρ/ρ_cosmic ranges from 10^5 (core) to 0.1 (outskirts)

3. JEANS EQUATION ANALYSIS:
   - Derived σ(r) from spherical Jeans equation
   - Modified with G_eff(ρ) for Synchronism

4. KEY PREDICTIONS:
   - ΛCDM: σ(3 Mpc) / σ(0.2 Mpc) = {ratio_LCDM:.3f}
   - Synchronism: σ(3 Mpc) / σ(0.2 Mpc) = {ratio_Sync:.3f}
   - Enhancement at r > R_200: {(ratio_Sync/ratio_LCDM - 1)*100:.1f}%

5. SESSION #173 INTERPRETATION:
   - Observed ratio = 3.28 used different radii (5 vs 40 Mpc)
   - Outer region is FIELD, not cluster outskirts
   - Need intra-cluster σ(r) profile for proper test

6. RECOMMENDED OBSERVATIONAL TESTS:
   - Planck PSZ2 clusters with spectroscopic members
   - σ(r) profiles from SDSS/DESI spectroscopy
   - Stacked profiles for statistical significance

7. NEXT STEPS:
   - Download PSZ2 velocity dispersion data
   - Analyze σ(r) profiles for individual clusters
   - Compare dynamical vs lensing masses

THEORETICAL OUTPUT:
- session176_theory.png: Visualization of predictions
- Framework for interpreting cluster velocity data

MRH CHECK:
----------
✓ Cluster dynamics probes r ~ 0.1-10 Mpc
✓ This overlaps with coherence function transition scale
✓ This IS the correct MRH for testing Synchronism
  (Unlike peculiar velocities which probe wrong scale)
""")

print("=" * 70)
print("SESSION #176 THEORETICAL FRAMEWORK COMPLETE")
print("=" * 70)
