#!/usr/bin/env python3
"""
SESSION #137: ULTRA-DIFFUSE GALAXIES (UDGs) IN SYNCHRONISM
===========================================================

Date: December 17, 2025
Focus: Comprehensive analysis of UDGs as tests of Synchronism

Background:
- UDGs are galaxies with stellar masses of dwarfs but sizes of giants
- Effective radii: R_e ~ 1.5-5 kpc (like Milky Way)
- Stellar masses: M* ~ 10^7-10^8 M_sun (like dwarfs)
- Surface brightness: very low (μ_0 > 24 mag/arcsec²)

Previous Synchronism sessions on UDGs:
- Session #90: TDG analysis (not discriminating due to similar Σ)
- Session #93: UDG test formulated (V/V_bar 30% higher)
- Session #97: DF2/DF4 analyzed (consistent with tidal stripping)

Challenge:
UDGs present puzzles for both ΛCDM and MOND:
1. Some UDGs appear "dark matter free" (DF2, DF4)
2. Some UDGs have extreme dark matter content
3. MOND predicts specific velocity dispersions based on mass
4. How does Synchronism explain this diversity?

This session will:
1. Model C(ρ) profiles for typical UDGs
2. Predict velocity dispersions
3. Analyze the "dark matter free" galaxies
4. Compare to ΛCDM and MOND predictions
5. Identify discriminating tests
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

print("=" * 70)
print("SESSION #137: ULTRA-DIFFUSE GALAXIES IN SYNCHRONISM")
print("=" * 70)
print("Date: December 17, 2025")
print("Focus: UDGs as stress tests for Synchronism")
print("=" * 70)

# Physical constants
G = 6.674e-11  # m³/kg/s²
c = 3e8  # m/s
M_sun = 2e30  # kg
pc = 3.086e16  # m
kpc = 1000 * pc
rho_crit = 9.2e-27  # kg/m³

# Golden ratio
phi = (1 + np.sqrt(5)) / 2

# MOND acceleration
a0 = 1.2e-10  # m/s²

print("\n" + "=" * 70)
print("PART 1: UDG PROPERTIES AND SAMPLE")
print("=" * 70)

print("""
ULTRA-DIFFUSE GALAXIES:
=======================

Definition (van Dokkum et al. 2015):
- Central surface brightness: μ_0 > 24 mag/arcsec² (g-band)
- Effective radius: R_e > 1.5 kpc

Physical properties:
- Stellar mass: M* ~ 10^7 - 10^8 M_sun
- Very low surface density: Σ* ~ 1-10 M_sun/pc²
- Gas-poor (typically)
- Found in clusters and field

Key examples:
1. Dragonfly 44 (DF44): Extreme dark matter content
2. NGC 1052-DF2: Apparently "dark matter free"
3. NGC 1052-DF4: Also "dark matter free"
4. VCC 1287: Cluster UDG with normal DM
""")

# UDG sample (representative cases)
udg_sample = {
    'DF44': {
        'M_star': 3e8 * M_sun,      # Stellar mass
        'R_e': 4.6 * kpc,            # Effective radius
        'sigma_obs': 47,             # Observed velocity dispersion (km/s)
        'sigma_err': 8,              # Error (km/s)
        'environment': 'Coma cluster',
        'notes': 'Extreme DM content claimed'
    },
    'DF2': {
        'M_star': 2e8 * M_sun,
        'R_e': 2.2 * kpc,
        'sigma_obs': 8.5,            # Very low!
        'sigma_err': 2.3,
        'environment': 'NGC 1052 group',
        'notes': 'Dark matter free?'
    },
    'DF4': {
        'M_star': 1.5e8 * M_sun,
        'R_e': 1.6 * kpc,
        'sigma_obs': 4.2,            # Extremely low
        'sigma_err': 2.2,
        'environment': 'NGC 1052 group',
        'notes': 'Dark matter free?'
    },
    'VCC1287': {
        'M_star': 4.3e8 * M_sun,
        'R_e': 3.3 * kpc,
        'sigma_obs': 33,
        'sigma_err': 10,
        'environment': 'Virgo cluster',
        'notes': 'Typical cluster UDG'
    },
    'DGSAT_I': {
        'M_star': 2.3e7 * M_sun,
        'R_e': 4.7 * kpc,
        'sigma_obs': 56,             # High for its mass
        'sigma_err': 10,
        'environment': 'Field',
        'notes': 'Isolated UDG'
    }
}

print(f"\nUDG Sample Properties:")
print(f"{'Name':<12} {'M*':<12} {'R_e':<10} {'σ_obs':<12} {'Environment':<15}")
print("-" * 65)
for name, props in udg_sample.items():
    print(f"{name:<12} {props['M_star']/M_sun:.1e}   {props['R_e']/kpc:.1f} kpc    "
          f"{props['sigma_obs']:.1f} ± {props['sigma_err']:.1f}  {props['environment']:<15}")

print("\n" + "=" * 70)
print("PART 2: SYNCHRONISM COHERENCE PROFILES FOR UDGs")
print("=" * 70)

def coherence(rho, Omega_m=0.315, B=phi, rho_t=1e-21):
    """
    Derived coherence function (Session #131):
    C(ρ) = Ω_m + (1 - Ω_m) × (ρ/ρ_t)^(1/B) / [1 + (ρ/ρ_t)^(1/B)]
    """
    if np.any(rho <= 0):
        return np.where(rho > 0, coherence(np.maximum(rho, 1e-30)), 0.315)
    x = (rho / rho_t) ** (1/B)
    C = Omega_m + (1 - Omega_m) * x / (1 + x)
    return np.clip(C, Omega_m, 1.0)

def sersic_density_3d(r, M_star, R_e, n=1.0):
    """
    3D density from deprojected Sérsic profile.

    For UDGs, typically n ~ 0.7-1.5 (nearly exponential)

    Using approximation:
    ρ(r) ≈ (M_star / (4π R_e³)) × b_n^(3n) × exp(-b_n × (r/R_e)^(1/n)) / (r/R_e)^(1-1/n)

    Simplified: exponential (n=1) → ρ(r) = ρ_0 × exp(-r/R_d)
    where R_d = R_e / 1.678 for n=1
    """
    b_n = 1.999 * n - 0.327  # Approximation for b_n
    R_d = R_e / b_n  # Scale length

    # Central density (approximate)
    rho_0 = M_star / (8 * np.pi * R_d**3)  # For exponential disk

    # 3D density profile
    rho = rho_0 * np.exp(-r / R_d)

    return rho

def compute_udg_profile(name, props, r_range=None):
    """
    Compute coherence and G_eff profiles for a UDG.
    """
    M_star = props['M_star']
    R_e = props['R_e']

    if r_range is None:
        r_range = np.logspace(-1, 2, 200) * kpc  # 0.1 to 100 kpc

    # Stellar density profile
    rho_star = sersic_density_3d(r_range, M_star, R_e, n=1.0)

    # Convert to kg/m³
    rho_star_kg = rho_star  # Already in kg/m³

    # Coherence profile
    C_profile = coherence(rho_star_kg)

    # G_eff/G profile
    G_eff_ratio = 1 / C_profile

    return r_range, rho_star_kg, C_profile, G_eff_ratio

print("""
DENSITY AND COHERENCE PROFILES:
===============================

UDGs have very LOW surface densities, which means:
- Lower average 3D densities than normal galaxies
- C will be LOWER throughout
- G_eff will be HIGHER → enhanced gravity

This could explain the diversity:
- Isolated UDGs: Low environmental density → Low C → High G_eff
- Cluster UDGs: Higher environmental density → Higher C → Lower G_eff
""")

# Compute profiles for all UDGs
r_plot = np.logspace(-1, 2, 200) * kpc

print(f"\nCoherence at R_e for each UDG:")
print(f"{'Name':<12} {'ρ(R_e) kg/m³':<15} {'C(R_e)':<10} {'G_eff/G':<10}")
print("-" * 50)

udg_profiles = {}
for name, props in udg_sample.items():
    r_range, rho, C, G_ratio = compute_udg_profile(name, props, r_plot)
    udg_profiles[name] = {'r': r_range, 'rho': rho, 'C': C, 'G_ratio': G_ratio}

    # Find values at R_e
    R_e = props['R_e']
    idx = np.argmin(np.abs(r_range - R_e))
    rho_Re = rho[idx]
    C_Re = C[idx]
    G_Re = G_ratio[idx]

    print(f"{name:<12} {rho_Re:<15.2e} {C_Re:<10.4f} {G_Re:<10.3f}")

print("\n" + "=" * 70)
print("PART 3: VELOCITY DISPERSION PREDICTIONS")
print("=" * 70)

print("""
VELOCITY DISPERSION IN SYNCHRONISM:
===================================

For a pressure-supported system:
σ² ≈ G_eff × M / R

In Synchronism:
σ²_sync = (G / C_avg) × M / R = σ²_Newton / C_avg

For UDGs with low C (low density):
→ σ_sync > σ_Newton

PREDICTIONS:
============
1. Standard Newtonian: σ² = G × M / R
2. MOND: σ = (G × M × a0)^(1/4)  (deep MOND)
3. Synchronism: σ² = (G / C_avg) × M / R
""")

def sigma_newtonian(M, R):
    """Newtonian velocity dispersion (virial theorem)."""
    return np.sqrt(G * M / R) / 1000  # km/s

def sigma_mond(M):
    """MOND velocity dispersion (deep MOND regime)."""
    return (G * M * a0) ** 0.25 / 1000  # km/s

def sigma_synchronism(M, R, C_avg):
    """Synchronism velocity dispersion."""
    G_eff = G / C_avg
    return np.sqrt(G_eff * M / R) / 1000  # km/s

def compute_average_C(name, props):
    """
    Compute mass-weighted average coherence within R_e.
    """
    M_star = props['M_star']
    R_e = props['R_e']

    # Integrate C × ρ × 4πr² dr from 0 to R_e, divide by M(<R_e)
    def integrand_C(r):
        rho = sersic_density_3d(r, M_star, R_e, n=1.0)
        C = coherence(rho)
        return C * rho * 4 * np.pi * r**2

    def integrand_M(r):
        rho = sersic_density_3d(r, M_star, R_e, n=1.0)
        return rho * 4 * np.pi * r**2

    numerator, _ = quad(integrand_C, 1e3, R_e, limit=100)  # Start from 1 km to avoid singularity
    denominator, _ = quad(integrand_M, 1e3, R_e, limit=100)

    if denominator > 0:
        return numerator / denominator
    return 0.315  # Fallback to minimum

print(f"\nVelocity Dispersion Predictions:")
print(f"{'Name':<12} {'σ_obs':<12} {'σ_Newton':<12} {'σ_MOND':<12} {'C_avg':<10} {'σ_Sync':<12}")
print("-" * 75)

predictions = []
for name, props in udg_sample.items():
    M_star = props['M_star']
    R_e = props['R_e']
    sigma_obs = props['sigma_obs']
    sigma_err = props['sigma_err']

    # Predictions
    sig_N = sigma_newtonian(M_star, R_e)
    sig_M = sigma_mond(M_star)
    C_avg = compute_average_C(name, props)
    sig_S = sigma_synchronism(M_star, R_e, C_avg)

    predictions.append({
        'name': name,
        'sigma_obs': sigma_obs,
        'sigma_err': sigma_err,
        'sigma_N': sig_N,
        'sigma_M': sig_M,
        'sigma_S': sig_S,
        'C_avg': C_avg
    })

    print(f"{name:<12} {sigma_obs:>5.1f} ± {sigma_err:<4.1f}  {sig_N:>8.1f}      {sig_M:>8.1f}      "
          f"{C_avg:<10.4f} {sig_S:>8.1f}")

print("\n" + "=" * 70)
print("PART 4: ANALYSIS OF 'DARK MATTER FREE' GALAXIES")
print("=" * 70)

print("""
THE DF2/DF4 PUZZLE:
===================

NGC 1052-DF2 and DF4 show unusually LOW velocity dispersions:
- DF2: σ = 8.5 ± 2.3 km/s (expected ~20 km/s for M*)
- DF4: σ = 4.2 ± 2.2 km/s (expected ~15 km/s for M*)

Interpretations:
1. ΛCDM: These galaxies lack dark matter halos
2. MOND: Problematic - should show MOND effects
3. External Field Effect (EFE): Environmental suppression

SYNCHRONISM INTERPRETATION:
===========================

These galaxies are in the NGC 1052 GROUP, which has:
- Massive elliptical NGC 1052 as host
- Environmental density higher than field
- UDGs may be tidally influenced

Key insight: C depends on LOCAL density, which includes:
1. Galaxy's own stellar density
2. Environmental contribution (group medium)
3. Tidal heating effects

If these UDGs are TIDALLY STRIPPED:
- Lower stellar mass → Lower ρ → Lower C → BUT...
- Tidal disruption increases velocity dispersion
- If C is elevated by environment → σ could be LOWER

Let's model this:
""")

def sigma_sync_with_environment(M_star, R_e, C_stellar, rho_env):
    """
    Synchronism prediction including environmental density contribution.

    The effective C is a combination of:
    - Stellar contribution (from galaxy itself)
    - Environmental contribution (from group/cluster medium)
    """
    # Environmental coherence
    C_env = coherence(rho_env)

    # Total coherence (dominant term wins due to logistic function)
    # At any point, C = f(ρ_total) where ρ_total = ρ_star + ρ_env
    # For simplicity, use max(C_stellar, C_env) as approximation
    C_total = max(C_stellar, C_env)

    G_eff = G / C_total
    sigma = np.sqrt(G_eff * M_star / R_e) / 1000
    return sigma, C_total

# NGC 1052 group environmental density
# NGC 1052 is at ~20 Mpc, group has hot gas
rho_NGC1052_group = 1e-24  # kg/m³ (intragroup medium, rough estimate)

# For comparison, Coma cluster core
rho_Coma_core = 1e-23  # kg/m³

# Virgo cluster
rho_Virgo = 5e-24  # kg/m³

print(f"\nEnvironmental densities:")
print(f"  NGC 1052 group: ρ ~ {rho_NGC1052_group:.0e} kg/m³ → C_env = {coherence(rho_NGC1052_group):.4f}")
print(f"  Virgo cluster: ρ ~ {rho_Virgo:.0e} kg/m³ → C_env = {coherence(rho_Virgo):.4f}")
print(f"  Coma cluster: ρ ~ {rho_Coma_core:.0e} kg/m³ → C_env = {coherence(rho_Coma_core):.4f}")

print(f"\n\nDF2/DF4 with Environmental Correction:")
print(f"{'Galaxy':<10} {'C_stellar':<12} {'C_env':<10} {'C_total':<10} {'σ_pred':<10} {'σ_obs':<12}")
print("-" * 70)

for name in ['DF2', 'DF4']:
    props = udg_sample[name]
    M_star = props['M_star']
    R_e = props['R_e']

    # Get stellar C_avg
    C_stellar = compute_average_C(name, props)

    # Environmental contribution
    C_env = coherence(rho_NGC1052_group)

    # With environment
    sigma_pred, C_total = sigma_sync_with_environment(M_star, R_e, C_stellar, rho_NGC1052_group)

    print(f"{name:<10} {C_stellar:<12.4f} {C_env:<10.4f} {C_total:<10.4f} "
          f"{sigma_pred:<10.1f} {props['sigma_obs']:.1f} ± {props['sigma_err']:.1f}")

print("""
INSIGHT:
========

The environmental density in the NGC 1052 group (~10⁻²⁴ kg/m³) gives
C_env ~ 0.32-0.35, which is LOWER than what we need to explain the
low velocity dispersions.

However, if these galaxies are:
1. TIDALLY STRIPPED (lost outer mass)
2. NOT in virial equilibrium
3. Have complex velocity fields

Then the simple σ² = GM/R relation doesn't apply directly.

Alternative explanation:
- DF2/DF4 may have formed from tidal debris
- Their kinematics reflect formation history, not virial equilibrium
- Synchronism + tidal history could explain the diversity
""")

print("\n" + "=" * 70)
print("PART 5: CLUSTER VS FIELD UDGs")
print("=" * 70)

print("""
ENVIRONMENTAL DEPENDENCE:
=========================

Synchronism predicts:
- Field UDGs: Lower environmental C → Higher G_eff → Higher σ
- Cluster UDGs: Higher environmental C → Lower G_eff → Lower σ

This is OPPOSITE to tidal heating expectations in ΛCDM!
In ΛCDM, cluster UDGs should be tidally heated → higher σ.

Let's test this prediction with our sample:
""")

# Classify by environment
field_udgs = ['DGSAT_I']
cluster_udgs = ['DF44', 'VCC1287']
group_udgs = ['DF2', 'DF4']

print(f"\nEnvironmental Classification:")
print(f"\n  FIELD UDGs (isolated): {field_udgs}")
print(f"  CLUSTER UDGs (Coma, Virgo): {cluster_udgs}")
print(f"  GROUP UDGs (NGC 1052): {group_udgs}")

# Compute predictions with environment
env_densities = {
    'field': 1e-26,           # Cosmic void
    'group': 1e-24,           # Group medium
    'cluster': 5e-24          # Cluster outskirts
}

print(f"\n\nPredictions vs Observations by Environment:")
print(f"{'Category':<12} {'UDG':<12} {'ρ_env':<12} {'C_env':<10} {'σ_pred':<10} {'σ_obs':<12} {'Match?':<8}")
print("-" * 85)

results_by_env = []
for category, udg_list in [('Field', field_udgs), ('Cluster', cluster_udgs), ('Group', group_udgs)]:
    for name in udg_list:
        props = udg_sample[name]
        M_star = props['M_star']
        R_e = props['R_e']
        sigma_obs = props['sigma_obs']
        sigma_err = props['sigma_err']

        # Environmental density
        if category == 'Field':
            rho_env = env_densities['field']
        elif category == 'Cluster':
            rho_env = env_densities['cluster']
        else:
            rho_env = env_densities['group']

        C_env = coherence(rho_env)
        C_stellar = compute_average_C(name, props)

        # Prediction
        sigma_pred, C_total = sigma_sync_with_environment(M_star, R_e, C_stellar, rho_env)

        # Check if matches within 2σ
        match = abs(sigma_pred - sigma_obs) < 2 * sigma_err
        match_str = "✓" if match else "✗"

        results_by_env.append({
            'category': category,
            'name': name,
            'sigma_pred': sigma_pred,
            'sigma_obs': sigma_obs,
            'sigma_err': sigma_err,
            'match': match
        })

        print(f"{category:<12} {name:<12} {rho_env:<12.0e} {C_env:<10.4f} "
              f"{sigma_pred:<10.1f} {sigma_obs:.1f} ± {sigma_err:<5.1f} {match_str:<8}")

# Summary
n_match = sum(1 for r in results_by_env if r['match'])
n_total = len(results_by_env)
print(f"\nMatches: {n_match}/{n_total} within 2σ")

print("\n" + "=" * 70)
print("PART 6: DISCRIMINATING PREDICTIONS")
print("=" * 70)

print("""
SYNCHRONISM vs MOND vs ΛCDM:
============================

For UDGs, the three frameworks make DIFFERENT predictions:

1. ΛCDM (with DM halos):
   - σ determined by DM halo mass
   - Scatter from halo-to-halo variation
   - No systematic trend with environment (except tidal heating)

2. MOND:
   - σ = (G M a0)^(1/4) in deep MOND
   - External Field Effect in groups/clusters
   - σ_MOND should be LOWER in dense environments (EFE)

3. SYNCHRONISM:
   - σ² = G M / (C × R)
   - C increases with density
   - σ_Sync should be LOWER in dense environments (higher C)

KEY DIFFERENCE:
===============
Both MOND and Synchronism predict LOWER σ in dense environments,
but for different reasons:
- MOND: External field breaks self-gravity enhancement
- Synchronism: Higher C reduces G_eff

The QUANTITATIVE predictions differ:
""")

# Compare predictions for all UDGs
print(f"\nQuantitative Comparison:")
print(f"{'UDG':<12} {'σ_ΛCDM':<10} {'σ_MOND':<10} {'σ_Sync':<10} {'σ_obs':<12} {'Best Fit':<12}")
print("-" * 70)

for p in predictions:
    name = p['name']
    sigma_obs = p['sigma_obs']
    sigma_err = p['sigma_err']

    # ΛCDM: Use M*/L ratio to estimate halo mass, then predict σ
    # Simplified: σ_ΛCDM ≈ σ_Newton × (M_halo/M*)^0.3
    # Assume typical M_halo/M* ~ 100 for UDGs
    sigma_LCDM = p['sigma_N'] * (100)**0.3  # Rough estimate

    sigma_MOND = p['sigma_M']
    sigma_Sync = p['sigma_S']

    # Find best fit
    residuals = {
        'ΛCDM': abs(sigma_LCDM - sigma_obs) / sigma_err,
        'MOND': abs(sigma_MOND - sigma_obs) / sigma_err,
        'Sync': abs(sigma_Sync - sigma_obs) / sigma_err
    }
    best = min(residuals, key=residuals.get)

    print(f"{name:<12} {sigma_LCDM:<10.1f} {sigma_MOND:<10.1f} {sigma_Sync:<10.1f} "
          f"{sigma_obs:.1f} ± {sigma_err:<5.1f} {best:<12}")

print("""
ANALYSIS:
=========

1. DF44 (high σ): ΛCDM can explain with massive DM halo
   Synchronism predicts lower (better match requires lower C)

2. DF2/DF4 (low σ): Challenge for all theories
   - ΛCDM: Requires DM-free formation
   - MOND: Violates unless EFE is strong
   - Synchronism: Requires high C (environmental?)

3. VCC1287, DGSAT_I: Intermediate cases
   All theories can fit with appropriate parameters

CONCLUSION:
===========
The current UDG sample does NOT clearly discriminate between theories.
However, systematic studies of:
- σ vs environment for many UDGs
- σ vs surface brightness
- σ vs distance from cluster center

Could provide discriminating tests.
""")

print("\n" + "=" * 70)
print("PART 7: TESTABLE PREDICTIONS FOR FUTURE OBSERVATIONS")
print("=" * 70)

print("""
SYNCHRONISM-SPECIFIC PREDICTIONS FOR UDGs:
==========================================

1. SURFACE BRIGHTNESS - VELOCITY DISPERSION RELATION:
   Lower Σ → Lower ρ → Lower C → Higher G_eff → Higher σ

   Prediction: σ should INCREASE with decreasing surface brightness
   (at fixed stellar mass)

   This is OPPOSITE to what you'd expect from virial theorem alone!

2. RADIAL COHERENCE GRADIENT:
   C should increase toward center (higher ρ)

   Prediction: Inner regions have LOWER G_eff than outer regions
   → Velocity dispersion profile should be STEEPER than Newtonian

3. ENVIRONMENTAL TREND:
   σ_field > σ_cluster (at fixed M*, R_e)

   Prediction: Field UDGs should have systematically higher σ
   than cluster UDGs of same mass and size

4. DISTANCE FROM CLUSTER CENTER:
   UDGs near cluster centers (higher ρ_env, higher C) should have
   LOWER σ than UDGs in cluster outskirts

   Prediction: Anti-correlation of σ with cluster-centric distance
""")

# Quantitative predictions
print(f"\nQuantitative Predictions:")

# 1. Surface brightness - sigma relation
print(f"\n1. Surface Brightness Effect:")
print(f"   For UDG with M* = 10⁸ M_☉, R_e = 3 kpc:")

M_test = 1e8 * M_sun
R_e_test = 3 * kpc

for Sigma_label, rho_central in [('High Σ', 1e-20), ('Medium Σ', 1e-21), ('Low Σ', 1e-22)]:
    C_avg = coherence(rho_central)
    sigma = sigma_synchronism(M_test, R_e_test, C_avg)
    print(f"   {Sigma_label}: ρ_central ~ {rho_central:.0e} → C ~ {C_avg:.3f} → σ ~ {sigma:.1f} km/s")

# 2. Environmental trend
print(f"\n2. Environmental Effect:")
print(f"   For UDG with M* = 10⁸ M_☉, R_e = 3 kpc:")

for env_label, rho_env in [('Field', 1e-26), ('Group', 1e-24), ('Cluster', 1e-23)]:
    C_env = coherence(rho_env)
    # Assume stellar gives C ~ 0.4
    C_total = max(0.4, C_env)
    sigma = sigma_synchronism(M_test, R_e_test, C_total)
    print(f"   {env_label}: ρ_env ~ {rho_env:.0e} → C_eff ~ {C_total:.3f} → σ ~ {sigma:.1f} km/s")

print("\n" + "=" * 70)
print("PART 8: VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# 1. Coherence profiles
ax1 = axes[0, 0]
colors = plt.cm.viridis(np.linspace(0, 1, len(udg_sample)))
for (name, data), color in zip(udg_profiles.items(), colors):
    ax1.plot(data['r']/kpc, data['C'], color=color, label=name, lw=2)
ax1.axhline(0.315, color='gray', ls='--', alpha=0.5, label='C_min = Ω_m')
ax1.set_xlabel('Radius (kpc)')
ax1.set_ylabel('Coherence C(r)')
ax1.set_title('Coherence Profiles of UDGs')
ax1.legend(fontsize=8)
ax1.set_xlim(0.1, 20)
ax1.set_xscale('log')
ax1.grid(True, alpha=0.3)

# 2. G_eff/G profiles
ax2 = axes[0, 1]
for (name, data), color in zip(udg_profiles.items(), colors):
    ax2.plot(data['r']/kpc, data['G_ratio'], color=color, label=name, lw=2)
ax2.axhline(1.0, color='gray', ls='--', alpha=0.5, label='Newtonian')
ax2.set_xlabel('Radius (kpc)')
ax2.set_ylabel('G_eff / G')
ax2.set_title('Effective Gravity Enhancement')
ax2.legend(fontsize=8)
ax2.set_xlim(0.1, 20)
ax2.set_xscale('log')
ax2.grid(True, alpha=0.3)

# 3. Predicted vs observed sigma
ax3 = axes[0, 2]
obs = [p['sigma_obs'] for p in predictions]
pred_N = [p['sigma_N'] for p in predictions]
pred_M = [p['sigma_M'] for p in predictions]
pred_S = [p['sigma_S'] for p in predictions]
err = [p['sigma_err'] for p in predictions]
names = [p['name'] for p in predictions]

x = np.arange(len(predictions))
width = 0.2

ax3.bar(x - 1.5*width, pred_N, width, label='Newtonian', alpha=0.7)
ax3.bar(x - 0.5*width, pred_M, width, label='MOND', alpha=0.7)
ax3.bar(x + 0.5*width, pred_S, width, label='Synchronism', alpha=0.7)
ax3.errorbar(x + 1.5*width, obs, yerr=err, fmt='ko', capsize=3, label='Observed')
ax3.set_xticks(x)
ax3.set_xticklabels(names, rotation=45, ha='right')
ax3.set_ylabel('Velocity Dispersion (km/s)')
ax3.set_title('Predictions vs Observations')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)

# 4. C vs environment
ax4 = axes[1, 0]
rho_env_range = np.logspace(-28, -20, 100)
C_env_range = coherence(rho_env_range)
ax4.plot(rho_env_range, C_env_range, 'b-', lw=2)
ax4.axvline(1e-26, color='green', ls=':', label='Field')
ax4.axvline(1e-24, color='orange', ls=':', label='Group')
ax4.axvline(1e-23, color='red', ls=':', label='Cluster')
ax4.set_xscale('log')
ax4.set_xlabel('Environmental Density (kg/m³)')
ax4.set_ylabel('Coherence C')
ax4.set_title('C vs Environmental Density')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. Predicted environmental trend
ax5 = axes[1, 1]
envs = ['Field', 'Group', 'Cluster']
C_envs = [coherence(1e-26), coherence(1e-24), coherence(1e-23)]
sigma_envs = [sigma_synchronism(M_test, R_e_test, max(0.4, c)) for c in C_envs]

ax5.plot(envs, sigma_envs, 'bo-', lw=2, markersize=10)
ax5.set_xlabel('Environment')
ax5.set_ylabel('Predicted σ (km/s)')
ax5.set_title('Synchronism: σ vs Environment\n(M* = 10⁸ M_☉, R_e = 3 kpc)')
ax5.grid(True, alpha=0.3)

# 6. Summary
ax6 = axes[1, 2]
ax6.axis('off')
summary_text = """
SESSION #137 KEY FINDINGS:
==========================

ULTRA-DIFFUSE GALAXIES IN SYNCHRONISM:

1. UDG COHERENCE PROFILES:
   • Very low densities → C ~ 0.3-0.5
   • G_eff enhanced by 2-3× throughout

2. VELOCITY DISPERSION:
   • Synchronism predicts higher σ for lower Σ
   • Environmental C can reduce predicted σ

3. DF2/DF4 PUZZLE:
   • Low σ requires high C or non-equilibrium
   • Tidal origin may explain kinematics

4. DISCRIMINATING PREDICTIONS:
   • σ increases with decreasing Σ
   • Field UDGs: higher σ than cluster UDGs
   • Anti-correlation with cluster-centric distance

5. CURRENT STATUS:
   • Sample too small to discriminate theories
   • Systematic surveys needed

NEXT STEPS:
• Analyze larger UDG samples
• Test environmental trends
• Model tidal formation scenarios
"""
ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=9,
         verticalalignment='top', fontfamily='monospace')

plt.suptitle('Session #137: Ultra-Diffuse Galaxies in Synchronism', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('session137_udgs.png', dpi=150, bbox_inches='tight')
print("\nVisualization saved to session137_udgs.png")

print("\n" + "=" * 70)
print("SESSION #137 SUMMARY")
print("=" * 70)

summary = """
ULTRA-DIFFUSE GALAXY ANALYSIS COMPLETE:
=======================================

KEY FINDINGS:
=============

1. UDG COHERENCE:
   • UDGs have very low densities → C ~ 0.35-0.5
   • G_eff enhanced by factor 2-3
   • Should show enhanced velocity dispersions

2. PREDICTIONS VS OBSERVATIONS:
   • Pure stellar model underpredicts σ for some UDGs
   • Environmental density affects predictions
   • DF2/DF4 low σ requires special explanation

3. DF2/DF4 INTERPRETATION:
   • May be tidally formed/disrupted
   • Not in virial equilibrium
   • Environmental C doesn't fully explain

4. ENVIRONMENTAL TRENDS:
   • Synchronism predicts: σ_field > σ_cluster
   • Same direction as MOND (EFE)
   • Quantitative differences could discriminate

5. TESTABLE PREDICTIONS:
   • σ vs surface brightness (inverse correlation)
   • σ vs cluster-centric distance (anti-correlation)
   • Field vs cluster UDG systematics

CURRENT STATUS:
===============
⚠️ Sample too small for definitive test
⚠️ Environmental modeling needs refinement
✓ Framework makes clear predictions
✓ Testable with larger surveys

NEXT STEPS:
===========
1. Analyze SDSS/HSC UDG catalogs
2. Model tidal formation in Synchronism
3. Detailed DF2/DF4 N-body simulations
"""
print(summary)

results = {
    'n_udgs_analyzed': 5,
    'predictions_within_2sigma': f"{n_match}/{n_total}",
    'key_prediction': 'σ vs Σ inverse correlation',
    'environmental_trend': 'σ_field > σ_cluster',
    'DF2_DF4_status': 'Requires non-equilibrium or tidal model',
    'status': 'UDG analysis complete - larger samples needed'
}

print(f"\nFinal results: {results}")
