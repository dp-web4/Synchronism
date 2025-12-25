#!/usr/bin/env python3
"""
SESSION #177b: VOID vs CLUSTER GALAXY PREDICTION
=================================================
Date: December 24, 2025

Following Session #177's scale-dependent formalism, this script derives
a CONCRETE testable prediction:

PREDICTION:
Void galaxies should show MORE apparent "dark matter" than cluster galaxies
of the same stellar mass.

Quantification:
- How much stronger should the effect be?
- What observational signature would we see?
- How does this compare to ΛCDM predictions?
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("SESSION #177b: VOID vs CLUSTER GALAXY PREDICTION")
print("=" * 70)

# =============================================================================
# 1. SCALE-DEPENDENT COHERENCE FUNCTION
# =============================================================================

phi = (1 + np.sqrt(5)) / 2
Omega_m = 0.3

# From Session #177 fit
A = 124.84  # M_sun/pc³
alpha = -3.033

def rho_transition(L_kpc):
    """Scale-dependent transition density"""
    return A * L_kpc**alpha

def coherence_scaled(rho_local, L_kpc):
    """Scale-dependent coherence function"""
    rho_t = rho_transition(L_kpc)
    rho_ratio = rho_local / rho_t

    if np.isscalar(rho_ratio):
        if rho_ratio <= 0:
            return Omega_m
        x = rho_ratio ** (1/phi)
        return Omega_m + (1 - Omega_m) * x / (1 + x)
    else:
        result = np.full_like(rho_ratio, Omega_m, dtype=float)
        pos = rho_ratio > 0
        x = np.zeros_like(rho_ratio, dtype=float)
        x[pos] = rho_ratio[pos] ** (1/phi)
        result[pos] = Omega_m + (1 - Omega_m) * x[pos] / (1 + x[pos])
        return result

def G_eff_scaled(rho_local, L_kpc):
    """Scale-dependent effective gravity ratio"""
    return 1.0 / coherence_scaled(rho_local, L_kpc)

# =============================================================================
# 2. MODEL GALAXY PROFILES
# =============================================================================

print("\n" + "=" * 70)
print("1. MODEL GALAXY (MW-like, M* = 5×10^10 M☉)")
print("=" * 70)

# Milky Way-like galaxy parameters
M_star = 5e10  # M_sun
R_disk = 3.0   # kpc (scale length)

def galaxy_density_profile(r_kpc, M_star, R_disk):
    """
    Simplified exponential disk + bulge + gas profile.
    Returns density in M_sun/pc³
    """
    # Exponential disk
    Sigma_disk = (M_star / (2 * np.pi * R_disk**2)) * np.exp(-r_kpc / R_disk)  # M_sun/kpc²
    h_disk = 0.3  # kpc (disk thickness)
    rho_disk = Sigma_disk / (2 * h_disk * 1e9)  # Convert to M_sun/pc³

    # Bulge (simplified)
    r_bulge = 1.0  # kpc
    M_bulge = 0.2 * M_star
    rho_bulge = (M_bulge / (4/3 * np.pi * r_bulge**3)) / 1e9 * np.exp(-r_kpc/r_bulge)

    return rho_disk + rho_bulge

# Test at various radii
radii = np.array([1, 2, 5, 10, 15, 20, 30, 50, 100])  # kpc
densities = np.array([galaxy_density_profile(r, M_star, R_disk) for r in radii])

print("\nBaryonic density profile (before environment effects):")
print("-" * 50)
for r, rho in zip(radii, densities):
    print(f"  r = {r:>5} kpc: ρ_baryon = {rho:.2e} M☉/pc³")

# =============================================================================
# 3. VOID vs CLUSTER ENVIRONMENT
# =============================================================================

print("\n" + "=" * 70)
print("2. ENVIRONMENT EFFECT")
print("=" * 70)

print("""
In Synchronism, the ENVIRONMENT affects the transition density interpretation:

CLUSTER GALAXY (ρ_env ~ 100 ρ_cosmic):
- Galaxy sits in dense environment
- Environment "fills in" the density deficit at large radii
- Effective ρ at large r is HIGHER → less enhancement

VOID GALAXY (ρ_env ~ 0.2 ρ_cosmic):
- Galaxy sits in underdense environment
- No environmental contribution at large radii
- Effective ρ at large r is LOWER → more enhancement

KEY INSIGHT:
The environment modifies the density at large radii where
the galaxy profile has dropped below the baryonic contribution.
""")

# Cosmic mean density for reference
rho_cosmic = 4e-8  # M_sun/pc³ (approximately)

# Environment densities
rho_env_cluster = 100 * rho_cosmic   # In cluster environment
rho_env_field = 1 * rho_cosmic       # Field galaxy
rho_env_void = 0.2 * rho_cosmic      # In void

print(f"\nEnvironment densities:")
print(f"  Cluster: {rho_env_cluster:.2e} M☉/pc³ (ρ/ρ_cosmic = 100)")
print(f"  Field: {rho_env_field:.2e} M☉/pc³ (ρ/ρ_cosmic = 1)")
print(f"  Void: {rho_env_void:.2e} M☉/pc³ (ρ/ρ_cosmic = 0.2)")

# =============================================================================
# 4. EFFECTIVE DENSITY WITH ENVIRONMENT
# =============================================================================

print("\n" + "=" * 70)
print("3. EFFECTIVE DENSITY PROFILES")
print("=" * 70)

def effective_density(r_kpc, M_star, R_disk, rho_env):
    """
    Effective density = max(baryonic, environment)

    At small r: baryonic dominates
    At large r: environment floor matters
    """
    rho_baryon = galaxy_density_profile(r_kpc, M_star, R_disk)
    return np.maximum(rho_baryon, rho_env)

print("\nEffective densities at r = 50 kpc:")
for name, rho_env in [('Cluster', rho_env_cluster), ('Field', rho_env_field), ('Void', rho_env_void)]:
    rho_eff = effective_density(50, M_star, R_disk, rho_env)
    rho_baryon = galaxy_density_profile(50, M_star, R_disk)
    print(f"  {name}: ρ_eff = {rho_eff:.2e} M☉/pc³ (ρ_baryon = {rho_baryon:.2e})")

# =============================================================================
# 5. G_eff PROFILES BY ENVIRONMENT
# =============================================================================

print("\n" + "=" * 70)
print("4. G_eff PROFILES BY ENVIRONMENT")
print("=" * 70)

r_fine = np.logspace(0, 2, 100)  # 1 to 100 kpc

G_eff_cluster = np.zeros_like(r_fine)
G_eff_field = np.zeros_like(r_fine)
G_eff_void = np.zeros_like(r_fine)

for i, r in enumerate(r_fine):
    rho_cl = effective_density(r, M_star, R_disk, rho_env_cluster)
    rho_fl = effective_density(r, M_star, R_disk, rho_env_field)
    rho_vo = effective_density(r, M_star, R_disk, rho_env_void)

    G_eff_cluster[i] = G_eff_scaled(rho_cl, r)
    G_eff_field[i] = G_eff_scaled(rho_fl, r)
    G_eff_void[i] = G_eff_scaled(rho_vo, r)

print("\nG_eff/G at selected radii:")
print("-" * 70)
print(f"{'r (kpc)':>10} {'Cluster':>15} {'Field':>15} {'Void':>15} {'Void/Cluster':>15}")
print("-" * 70)

for r in [5, 10, 20, 30, 50, 100]:
    idx = np.argmin(np.abs(r_fine - r))
    ratio = G_eff_void[idx] / G_eff_cluster[idx]
    print(f"{r:>10} {G_eff_cluster[idx]:>15.3f} {G_eff_field[idx]:>15.3f} {G_eff_void[idx]:>15.3f} {ratio:>15.3f}")

# =============================================================================
# 6. ROTATION CURVE PREDICTION
# =============================================================================

print("\n" + "=" * 70)
print("5. ROTATION CURVE PREDICTION")
print("=" * 70)

# Circular velocity: v² = G_eff × G × M(<r) / r
# For simplicity, use enclosed mass from exponential disk

def enclosed_mass(r_kpc, M_star, R_disk):
    """Approximate enclosed mass for exponential disk"""
    x = r_kpc / R_disk
    # Approximate: M(<r) = M_star × [1 - (1 + x)exp(-x)]
    return M_star * (1 - (1 + x) * np.exp(-x))

def rotation_velocity(r_kpc, M_star, R_disk, rho_env):
    """Rotation velocity including G_eff enhancement"""
    M_enc = enclosed_mass(r_kpc, M_star, R_disk)
    rho_eff = effective_density(r_kpc, M_star, R_disk, rho_env)
    G_ratio = G_eff_scaled(rho_eff, r_kpc)

    # G in (km/s)² kpc / M_sun = 4.302e-6
    G = 4.302e-6  # kpc (km/s)² / M_sun

    v_squared = G_ratio * G * M_enc / r_kpc
    return np.sqrt(v_squared)

# Compute rotation curves
v_cluster = np.array([rotation_velocity(r, M_star, R_disk, rho_env_cluster) for r in r_fine])
v_field = np.array([rotation_velocity(r, M_star, R_disk, rho_env_field) for r in r_fine])
v_void = np.array([rotation_velocity(r, M_star, R_disk, rho_env_void) for r in r_fine])

# Also compute Newtonian (no enhancement)
v_newton = np.array([rotation_velocity(r, M_star, R_disk, 1e10) for r in r_fine])  # High rho → G_eff ≈ G

print("\nRotation velocities at selected radii:")
print("-" * 80)
print(f"{'r (kpc)':>10} {'Newtonian':>12} {'Cluster':>12} {'Field':>12} {'Void':>12} {'Void/Newton':>12}")
print("-" * 80)

for r in [5, 10, 20, 30, 50, 100]:
    idx = np.argmin(np.abs(r_fine - r))
    ratio = v_void[idx] / v_newton[idx]
    print(f"{r:>10} {v_newton[idx]:>12.1f} {v_cluster[idx]:>12.1f} {v_field[idx]:>12.1f} {v_void[idx]:>12.1f} {ratio:>12.3f}")

# =============================================================================
# 7. DARK MATTER FRACTION PREDICTION
# =============================================================================

print("\n" + "=" * 70)
print("6. INFERRED 'DARK MATTER' FRACTION")
print("=" * 70)

print("""
If an observer uses Newtonian gravity to infer mass:

M_dyn = v² × r / G

And compares to baryonic mass:

f_DM = (M_dyn - M_baryon) / M_dyn = 1 - M_baryon/M_dyn

In Synchronism, the "excess" mass is actually G_eff enhancement.
""")

def dark_matter_fraction(r_kpc, M_star, R_disk, rho_env):
    """Inferred dark matter fraction assuming Newtonian gravity"""
    v = rotation_velocity(r_kpc, M_star, R_disk, rho_env)
    M_baryon = enclosed_mass(r_kpc, M_star, R_disk)

    G = 4.302e-6
    M_dyn = v**2 * r_kpc / G

    f_DM = 1 - M_baryon / M_dyn
    return max(0, f_DM)

print("\nInferred dark matter fraction at r = 50 kpc:")
print("-" * 50)
for name, rho_env in [('Cluster', rho_env_cluster), ('Field', rho_env_field), ('Void', rho_env_void)]:
    f_DM = dark_matter_fraction(50, M_star, R_disk, rho_env)
    print(f"  {name}: f_DM = {f_DM*100:.1f}%")

print("\nInferred dark matter fraction at r = 100 kpc:")
print("-" * 50)
for name, rho_env in [('Cluster', rho_env_cluster), ('Field', rho_env_field), ('Void', rho_env_void)]:
    f_DM = dark_matter_fraction(100, M_star, R_disk, rho_env)
    print(f"  {name}: f_DM = {f_DM*100:.1f}%")

# =============================================================================
# 8. KEY TESTABLE PREDICTION
# =============================================================================

print("\n" + "=" * 70)
print("7. KEY TESTABLE PREDICTION")
print("=" * 70)

# Calculate the difference at r = 50 kpc
f_DM_cluster = dark_matter_fraction(50, M_star, R_disk, rho_env_cluster)
f_DM_void = dark_matter_fraction(50, M_star, R_disk, rho_env_void)
v_ratio_50 = np.interp(50, r_fine, v_void) / np.interp(50, r_fine, v_cluster)

print(f"""
*** QUANTITATIVE PREDICTION ***

For identical galaxies (same M*, same R_disk) in different environments:

At r = 50 kpc:
  Cluster galaxy: f_DM = {f_DM_cluster*100:.1f}%, v_rot = {np.interp(50, r_fine, v_cluster):.1f} km/s
  Void galaxy: f_DM = {f_DM_void*100:.1f}%, v_rot = {np.interp(50, r_fine, v_void):.1f} km/s

Predicted differences:
  Δf_DM = {(f_DM_void - f_DM_cluster)*100:.1f}% (void has {(f_DM_void - f_DM_cluster)*100:.1f}% MORE "dark matter")
  v_void / v_cluster = {v_ratio_50:.3f} (void has {(v_ratio_50-1)*100:.1f}% HIGHER rotation velocity)

OBSERVATIONAL TEST:
==================

1. Select galaxies with similar:
   - Stellar mass (from photometry)
   - Morphology (disk galaxies)
   - Surface brightness profile

2. Classify by environment:
   - Cluster (ρ_env > 10 ρ_cosmic)
   - Field (ρ_env ~ 1 ρ_cosmic)
   - Void (ρ_env < 0.5 ρ_cosmic)

3. Measure rotation curves (HI or Hα)

4. Compare:
   - Outer rotation velocity
   - Inferred dark matter fraction

EXPECTED SIGNATURE:
- Void galaxies: HIGHER v_rot at large r
- Void galaxies: HIGHER inferred f_DM
- The effect should scale with environment density

ΛCDM PREDICTION:
In ΛCDM, environment affects dark matter halo properties:
- Cluster galaxies: May have truncated halos (tidal stripping)
- Field/void galaxies: Should have full halos

But ΛCDM predicts LESS dark matter in void galaxies
(lower NFW concentration, shallower profiles)

SYNCHRONISM predicts the OPPOSITE:
MORE apparent dark matter in void galaxies.

This is a DISCRIMINATING test!
""")

# =============================================================================
# 9. VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("8. GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: G_eff profiles
ax1 = axes[0, 0]
ax1.semilogx(r_fine, G_eff_cluster, 'b-', linewidth=2, label='Cluster (ρ_env = 100 ρ_cosmic)')
ax1.semilogx(r_fine, G_eff_field, 'g-', linewidth=2, label='Field (ρ_env = ρ_cosmic)')
ax1.semilogx(r_fine, G_eff_void, 'r-', linewidth=2, label='Void (ρ_env = 0.2 ρ_cosmic)')
ax1.axhline(1.0, color='gray', linestyle=':', linewidth=1)
ax1.set_xlabel('Radius (kpc)')
ax1.set_ylabel('G_eff / G')
ax1.set_title('Effective Gravity Enhancement by Environment')
ax1.legend(loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_ylim(0.8, 2.5)

# Panel 2: Rotation curves
ax2 = axes[0, 1]
ax2.semilogx(r_fine, v_newton, 'k--', linewidth=2, label='Newtonian (no enhancement)')
ax2.semilogx(r_fine, v_cluster, 'b-', linewidth=2, label='Cluster')
ax2.semilogx(r_fine, v_field, 'g-', linewidth=2, label='Field')
ax2.semilogx(r_fine, v_void, 'r-', linewidth=2, label='Void')
ax2.set_xlabel('Radius (kpc)')
ax2.set_ylabel('Rotation velocity (km/s)')
ax2.set_title(f'Rotation Curves (M* = {M_star:.0e} M☉)')
ax2.legend(loc='lower right')
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, 300)

# Panel 3: Dark matter fraction profiles
ax3 = axes[1, 0]
r_sample = np.logspace(0.5, 2, 50)
f_DM_c = np.array([dark_matter_fraction(r, M_star, R_disk, rho_env_cluster) for r in r_sample])
f_DM_f = np.array([dark_matter_fraction(r, M_star, R_disk, rho_env_field) for r in r_sample])
f_DM_v = np.array([dark_matter_fraction(r, M_star, R_disk, rho_env_void) for r in r_sample])

ax3.semilogx(r_sample, f_DM_c * 100, 'b-', linewidth=2, label='Cluster')
ax3.semilogx(r_sample, f_DM_f * 100, 'g-', linewidth=2, label='Field')
ax3.semilogx(r_sample, f_DM_v * 100, 'r-', linewidth=2, label='Void')
ax3.set_xlabel('Radius (kpc)')
ax3.set_ylabel('Inferred dark matter fraction (%)')
ax3.set_title('Environment-Dependent "Dark Matter" Fraction')
ax3.legend(loc='lower right')
ax3.grid(True, alpha=0.3)
ax3.set_ylim(0, 100)

# Panel 4: Environment density effect
ax4 = axes[1, 1]
rho_env_range = np.logspace(-9, -5, 50)  # M_sun/pc³
f_DM_50 = np.array([dark_matter_fraction(50, M_star, R_disk, rho) for rho in rho_env_range])
rho_cosmic_ratio = rho_env_range / rho_cosmic

ax4.semilogx(rho_cosmic_ratio, f_DM_50 * 100, 'b-', linewidth=2)
ax4.axvline(0.2, color='red', linestyle='--', alpha=0.7, label='Void')
ax4.axvline(1.0, color='green', linestyle='--', alpha=0.7, label='Field')
ax4.axvline(100, color='blue', linestyle='--', alpha=0.7, label='Cluster')
ax4.set_xlabel('Environment density (ρ_env / ρ_cosmic)')
ax4.set_ylabel('Inferred f_DM at r = 50 kpc (%)')
ax4.set_title('Dark Matter Fraction vs Environment Density')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.suptitle('Session #177b: Void vs Cluster Galaxy Prediction\nSynchronism: More "Dark Matter" in Voids',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session177b_void_cluster.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session177b_void_cluster.png")

# =============================================================================
# 10. SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #177b: SUMMARY")
print("=" * 70)

print(f"""
VOID vs CLUSTER GALAXY PREDICTION
=================================

1. THEORETICAL BASIS:
   Scale-dependent coherence function with environment effect:
   - Lower environment density → lower effective density at large r
   - Lower density → stronger G_eff enhancement
   - Stronger enhancement → higher rotation velocity

2. QUANTITATIVE PREDICTION (MW-like galaxy, r = 50 kpc):
   Cluster: f_DM = {f_DM_cluster*100:.1f}%, v_rot = {np.interp(50, r_fine, v_cluster):.1f} km/s
   Void: f_DM = {f_DM_void*100:.1f}%, v_rot = {np.interp(50, r_fine, v_void):.1f} km/s
   Difference: {(f_DM_void-f_DM_cluster)*100:.1f}% more "dark matter" in voids

3. OBSERVATIONAL TEST:
   - Compare rotation curves of matched galaxies in different environments
   - Void galaxies should show HIGHER v_rot at fixed r
   - This is OPPOSITE to ΛCDM prediction (truncated halos in clusters)

4. DISCRIMINATING POWER:
   ΛCDM: Cluster galaxies have LESS dark matter (tidal stripping)
   Synchronism: Void galaxies have MORE apparent dark matter (G_eff enhancement)

5. DATA SOURCES:
   - SPARC rotation curve database + environment classification
   - THINGS/LITTLE THINGS HI surveys
   - SDSS void catalogs

FILES CREATED:
- session177b_void_cluster.png

SIGNIFICANCE:
This is a CLEAR DISCRIMINATING TEST between Synchronism and ΛCDM.
The predictions are opposite in sign, not just magnitude.
""")

print("=" * 70)
print("SESSION #177b COMPLETE")
print("=" * 70)
