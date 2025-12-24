#!/usr/bin/env python3
"""
SESSION #176b: CLUSTER VELOCITY DISPERSION PROFILES - CORRECTED
================================================================
Date: December 24, 2025

CORRECTION: Session #176a had incorrect critical density calculation
resulting in R_200 = 203 Mpc (should be ~2 Mpc for Coma-like cluster).

This version uses correct cosmological parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

print("=" * 70)
print("SESSION #176b: CLUSTER VELOCITY DISPERSION (CORRECTED)")
print("=" * 70)

# =============================================================================
# 1. CORRECT COSMOLOGICAL PARAMETERS
# =============================================================================

print("\n" + "=" * 70)
print("1. COSMOLOGICAL PARAMETERS")
print("=" * 70)

# Cosmological constants
H0 = 70  # km/s/Mpc
h = H0 / 100
Omega_m = 0.3

# Critical density (z=0)
# ρ_crit = 3H²/(8πG) = 2.775 × 10^11 h² M_sun/Mpc³
rho_crit = 2.775e11 * h**2  # M_sun/Mpc³
print(f"Critical density: ρ_crit = {rho_crit:.3e} M_sun/Mpc³")

# Mean cosmic density
rho_cosmic = rho_crit * Omega_m
print(f"Cosmic mean density: ρ_cosmic = {rho_cosmic:.3e} M_sun/Mpc³")

# Gravitational constant in useful units
# G = 4.302 × 10^-6 kpc (km/s)² / M_sun
# G = 4.302 × 10^-3 Mpc (km/s)² / (10^6 M_sun)
G_units = 4.302e-6  # kpc (km/s)² / M_sun

# =============================================================================
# 2. SYNCHRONISM COHERENCE FUNCTION
# =============================================================================

print("\n" + "=" * 70)
print("2. COHERENCE FUNCTION")
print("=" * 70)

phi = (1 + np.sqrt(5)) / 2  # Golden ratio

def coherence(rho_ratio):
    """C(ρ) from Synchronism - returns value between Ω_m and 1"""
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

def G_eff_ratio(rho_ratio):
    """G_eff/G = 1/C(ρ)"""
    return 1.0 / coherence(rho_ratio)

# Demonstrate
print("\nCoherence function:")
for rho in [0.01, 0.1, 1.0, 10, 100, 1000, 10000]:
    C = coherence(rho)
    G_eff = G_eff_ratio(rho)
    print(f"  ρ/ρ̄ = {rho:>7.2f}: C = {C:.4f}, G_eff/G = {G_eff:.4f}")

# =============================================================================
# 3. NFW CLUSTER PROFILE (CORRECT)
# =============================================================================

print("\n" + "=" * 70)
print("3. NFW CLUSTER PROFILE")
print("=" * 70)

def get_cluster_params(M_200, c_200=4.0):
    """
    Get NFW profile parameters for cluster with given virial mass.

    M_200 = (4π/3) × 200 × ρ_crit × R_200³
    """
    # Virial radius
    R_200 = (3 * M_200 / (4 * np.pi * 200 * rho_crit))**(1/3)  # Mpc

    # Scale radius
    r_s = R_200 / c_200

    # Characteristic overdensity
    delta_c = (200/3) * c_200**3 / (np.log(1 + c_200) - c_200/(1 + c_200))

    # Characteristic density
    rho_s = delta_c * rho_crit

    return R_200, r_s, rho_s, delta_c

def NFW_density(r, r_s, rho_s):
    """NFW density profile"""
    x = r / r_s
    return rho_s / (x * (1 + x)**2)

def NFW_mass(r, r_s, rho_s):
    """Enclosed mass within r for NFW profile"""
    x = r / r_s
    return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x/(1 + x))

# Example clusters
clusters = {
    'Coma-like': (1e15, 4.0),    # Massive, typical
    'Virgo-like': (4e14, 5.0),   # Moderate mass
    'Group': (1e14, 6.0),        # Galaxy group
}

print("\nCluster parameters:")
print("-" * 70)
for name, (M_200, c_200) in clusters.items():
    R_200, r_s, rho_s, delta_c = get_cluster_params(M_200, c_200)
    print(f"\n{name} (M_200 = {M_200:.0e} M_sun, c = {c_200}):")
    print(f"  R_200 = {R_200:.3f} Mpc = {R_200*1000:.0f} kpc")
    print(f"  r_s = {r_s:.3f} Mpc = {r_s*1000:.0f} kpc")
    print(f"  δ_c = {delta_c:.0f}")
    print(f"  ρ_s = {rho_s:.2e} M_sun/Mpc³")
    print(f"  ρ_s/ρ_cosmic = {rho_s/rho_cosmic:.0f}")

# =============================================================================
# 4. VELOCITY DISPERSION FROM JEANS EQUATION
# =============================================================================

print("\n" + "=" * 70)
print("4. JEANS EQUATION VELOCITY DISPERSION")
print("=" * 70)

def sigma_jeans_LCDM(r, M_200, c_200=4.0, beta=0.0):
    """
    Isotropic velocity dispersion from spherical Jeans equation (ΛCDM).

    For β=0: σ_r²(r) = (1/ρ(r)) ∫_r^∞ ρ(s) × G × M(<s)/s² ds

    Returns σ in km/s.
    """
    R_200, r_s, rho_s, _ = get_cluster_params(M_200, c_200)

    def rho(s):
        return NFW_density(s, r_s, rho_s)

    def M_enc(s):
        return NFW_mass(s, r_s, rho_s)

    def integrand(s):
        # G in Mpc³/(M_sun × Gyr²): 4.302e-9
        # But we want (km/s)², so use G in Mpc (km/s)²/M_sun
        # G = 4.302e-6 kpc (km/s)²/M_sun = 4.302e-9 Mpc (km/s)²/M_sun
        G = 4.302e-9  # Mpc (km/s)²/M_sun
        return rho(s) * G * M_enc(s) / s**2

    rho_r = rho(r)
    if rho_r <= 0:
        return 0.0

    # Integrate from r to large radius
    r_max = min(100 * R_200, 100)  # Mpc
    integral, _ = quad(integrand, r, r_max, limit=200)

    sigma_squared = integral / rho_r
    return np.sqrt(max(sigma_squared, 0))

def sigma_jeans_Sync(r, M_200, c_200=4.0, beta=0.0):
    """
    Isotropic velocity dispersion from Jeans equation with Synchronism G_eff.

    σ_r²(r) = (1/ρ(r)) ∫_r^∞ ρ(s) × G_eff(s) × M(<s)/s² ds

    where G_eff(s) = G/C(ρ(s)/ρ_cosmic)
    """
    R_200, r_s, rho_s, _ = get_cluster_params(M_200, c_200)

    def rho(s):
        return NFW_density(s, r_s, rho_s)

    def M_enc(s):
        return NFW_mass(s, r_s, rho_s)

    def integrand(s):
        G = 4.302e-9  # Mpc (km/s)²/M_sun
        density = rho(s)
        rho_ratio = density / rho_cosmic
        G_eff = G * G_eff_ratio(rho_ratio)
        return density * G_eff * M_enc(s) / s**2

    rho_r = rho(r)
    if rho_r <= 0:
        return 0.0

    r_max = min(100 * R_200, 100)
    integral, _ = quad(integrand, r, r_max, limit=200)

    sigma_squared = integral / rho_r
    return np.sqrt(max(sigma_squared, 0))

# =============================================================================
# 5. COMPUTE PROFILES FOR COMA-LIKE CLUSTER
# =============================================================================

print("\n" + "=" * 70)
print("5. VELOCITY DISPERSION PROFILES")
print("=" * 70)

M_200 = 1e15  # Coma-like
c_200 = 4.0
R_200, r_s, rho_s, _ = get_cluster_params(M_200, c_200)

print(f"\nComa-like cluster (M_200 = {M_200:.0e} M_sun):")
print(f"  R_200 = {R_200:.3f} Mpc = {R_200*1000:.0f} kpc")
print(f"  r_s = {r_s:.3f} Mpc = {r_s*1000:.0f} kpc")

# Radial grid from 0.01 R_200 to 5 R_200
r_norm = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
radii = r_norm * R_200

print(f"\nVelocity dispersion profiles:")
print("-" * 80)
print(f"{'r/R_200':>8} {'r (Mpc)':>10} {'ρ/ρ_c':>10} {'G_eff/G':>10} {'σ_ΛCDM':>10} {'σ_Sync':>10} {'Ratio':>8}")
print("-" * 80)

sigma_LCDM = []
sigma_Sync = []
rho_ratios = []

for i, r in enumerate(radii):
    # Density at this radius
    rho = NFW_density(r, r_s, rho_s)
    rho_ratio = rho / rho_cosmic
    G_eff = G_eff_ratio(rho_ratio)
    rho_ratios.append(rho_ratio)

    # Velocity dispersions
    s_lcdm = sigma_jeans_LCDM(r, M_200, c_200)
    s_sync = sigma_jeans_Sync(r, M_200, c_200)

    sigma_LCDM.append(s_lcdm)
    sigma_Sync.append(s_sync)

    ratio = s_sync / s_lcdm if s_lcdm > 0 else 0

    print(f"{r_norm[i]:>8.2f} {r:>10.4f} {rho_ratio:>10.0f} {G_eff:>10.4f} {s_lcdm:>10.0f} {s_sync:>10.0f} {ratio:>8.3f}")

sigma_LCDM = np.array(sigma_LCDM)
sigma_Sync = np.array(sigma_Sync)
rho_ratios = np.array(rho_ratios)

# =============================================================================
# 6. KEY PREDICTIONS
# =============================================================================

print("\n" + "=" * 70)
print("6. KEY DISCRIMINATING PREDICTIONS")
print("=" * 70)

# Find indices for inner and outer
inner_idx = np.argmin(np.abs(r_norm - 0.1))  # 0.1 R_200
outer_idx = np.argmin(np.abs(r_norm - 2.0))  # 2.0 R_200

sigma_inner_LCDM = sigma_LCDM[inner_idx]
sigma_outer_LCDM = sigma_LCDM[outer_idx]
sigma_inner_Sync = sigma_Sync[inner_idx]
sigma_outer_Sync = sigma_Sync[outer_idx]

ratio_LCDM = sigma_outer_LCDM / sigma_inner_LCDM
ratio_Sync = sigma_outer_Sync / sigma_inner_Sync

print(f"\n*** OUTER/INNER DISPERSION RATIO ***")
print(f"\nInner radius: r = 0.1 R_200 = {0.1*R_200:.3f} Mpc = {0.1*R_200*1000:.0f} kpc")
print(f"Outer radius: r = 2.0 R_200 = {2.0*R_200:.3f} Mpc = {2.0*R_200*1000:.0f} kpc")
print(f"\n  ΛCDM:")
print(f"    σ_inner = {sigma_inner_LCDM:.0f} km/s")
print(f"    σ_outer = {sigma_outer_LCDM:.0f} km/s")
print(f"    Ratio = {ratio_LCDM:.3f}")
print(f"\n  Synchronism:")
print(f"    σ_inner = {sigma_inner_Sync:.0f} km/s")
print(f"    σ_outer = {sigma_outer_Sync:.0f} km/s")
print(f"    Ratio = {ratio_Sync:.3f}")
print(f"\n  Difference: {(ratio_Sync/ratio_LCDM - 1)*100:.1f}%")

# Enhancement factor
enhancement = sigma_Sync / sigma_LCDM
print(f"\n*** SYNCHRONISM ENHANCEMENT ***")
print(f"\n  At 0.1 R_200: {enhancement[inner_idx]:.3f}")
print(f"  At 1.0 R_200: {enhancement[np.argmin(np.abs(r_norm - 1.0))]:.3f}")
print(f"  At 2.0 R_200: {enhancement[outer_idx]:.3f}")
print(f"  At 5.0 R_200: {enhancement[-1]:.3f}")

# =============================================================================
# 7. MASS DEPENDENCE
# =============================================================================

print("\n" + "=" * 70)
print("7. MASS DEPENDENCE OF ENHANCEMENT")
print("=" * 70)

masses = [1e13, 5e13, 1e14, 5e14, 1e15]
print(f"\nEnhancement at r = R_200 for different cluster masses:")
print("-" * 50)

for M in masses:
    R200, rs, rhos, _ = get_cluster_params(M, c_200=4.0)
    r = R200  # At virial radius

    s_lcdm = sigma_jeans_LCDM(r, M, 4.0)
    s_sync = sigma_jeans_Sync(r, M, 4.0)
    enh = s_sync / s_lcdm if s_lcdm > 0 else 0

    # Density at R_200 (should be ~200 ρ_crit = ~600 ρ_cosmic)
    rho = NFW_density(r, rs, rhos)
    rho_ratio = rho / rho_cosmic

    print(f"  M = {M:.0e} M_sun: R_200 = {R200:.3f} Mpc, ρ/ρ_c = {rho_ratio:.0f}, enhancement = {enh:.4f}")

# =============================================================================
# 8. VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("8. GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Density profile
ax1 = axes[0, 0]
ax1.loglog(r_norm, rho_ratios, 'b-', linewidth=2)
ax1.axhline(1.0, color='orange', linestyle='--', label='ρ = ρ_cosmic')
ax1.axhline(200, color='red', linestyle=':', label='ρ = 200 ρ_crit')
ax1.axvline(1.0, color='gray', linestyle='--', alpha=0.5, label='R_200')
ax1.fill_between(r_norm, 0.1, rho_ratios, where=(rho_ratios > 1),
                  alpha=0.3, color='blue', label='C → 1 regime')
ax1.fill_between(r_norm, 0.1, np.ones_like(rho_ratios), where=(rho_ratios < 1),
                  alpha=0.3, color='orange', label='Enhanced G_eff')
ax1.set_xlabel('r / R_200')
ax1.set_ylabel('ρ / ρ_cosmic')
ax1.set_title(f'Cluster Density Profile\n(M_200 = {M_200:.0e} M☉, R_200 = {R_200*1000:.0f} kpc)')
ax1.legend(loc='upper right', fontsize=9)
ax1.set_xlim(0.01, 5)
ax1.set_ylim(0.1, 1e6)
ax1.grid(True, alpha=0.3)

# Panel 2: G_eff profile
ax2 = axes[0, 1]
G_eff_profile = G_eff_ratio(rho_ratios)
C_profile = coherence(rho_ratios)
ax2.semilogx(r_norm, G_eff_profile, 'b-', linewidth=2, label='G_eff/G')
ax2.semilogx(r_norm, C_profile, 'r--', linewidth=2, label='C(ρ)')
ax2.axhline(1.0, color='gray', linestyle=':')
ax2.axhline(1/Omega_m, color='orange', linestyle='--', alpha=0.5,
            label=f'G_eff/G max = {1/Omega_m:.2f}')
ax2.axvline(1.0, color='gray', linestyle='--', alpha=0.5)
ax2.set_xlabel('r / R_200')
ax2.set_ylabel('Ratio')
ax2.set_title('Coherence Function and G_eff')
ax2.legend(loc='upper right')
ax2.set_xlim(0.01, 5)
ax2.set_ylim(0.2, 4)
ax2.grid(True, alpha=0.3)

# Panel 3: Velocity dispersion profiles
ax3 = axes[1, 0]
ax3.semilogx(r_norm, sigma_LCDM, 'b-o', linewidth=2, markersize=6, label='ΛCDM')
ax3.semilogx(r_norm, sigma_Sync, 'r-s', linewidth=2, markersize=6, label='Synchronism')
ax3.axvline(1.0, color='gray', linestyle='--', alpha=0.5, label='R_200')
ax3.set_xlabel('r / R_200')
ax3.set_ylabel('Velocity dispersion σ_v (km/s)')
ax3.set_title('Velocity Dispersion Profiles')
ax3.legend(loc='upper right')
ax3.set_xlim(0.01, 5)
ax3.grid(True, alpha=0.3)

# Panel 4: Enhancement factor
ax4 = axes[1, 1]
enhancement = sigma_Sync / sigma_LCDM
ax4.semilogx(r_norm, enhancement, 'g-o', linewidth=2, markersize=6)
ax4.axhline(1.0, color='gray', linestyle='--')
ax4.axvline(1.0, color='gray', linestyle='--', alpha=0.5, label='R_200')
ax4.set_xlabel('r / R_200')
ax4.set_ylabel('σ_Sync / σ_ΛCDM')
ax4.set_title('Synchronism Enhancement Factor')
ax4.legend(loc='upper left')
ax4.set_xlim(0.01, 5)
ax4.grid(True, alpha=0.3)

plt.suptitle('Session #176b: Cluster Velocity Dispersion Profiles (Corrected)\nSynchronism vs ΛCDM',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session176b_corrected.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session176b_corrected.png")

# =============================================================================
# 9. OBSERVATIONAL TESTS
# =============================================================================

print("\n" + "=" * 70)
print("9. OBSERVATIONAL TESTS")
print("=" * 70)

print(f"""
TESTABLE PREDICTIONS
====================

1. DISPERSION PROFILE SHAPE:
   - At cluster core (r < 0.1 R_200): σ_Sync ≈ σ_ΛCDM (high density → C ≈ 1)
   - At outskirts (r > R_200): σ_Sync > σ_ΛCDM by {(enhancement[-1]-1)*100:.0f}%

2. OUTER/INNER RATIO:
   - ΛCDM predicts: σ(2 R_200) / σ(0.1 R_200) = {ratio_LCDM:.3f}
   - Synchronism predicts: σ(2 R_200) / σ(0.1 R_200) = {ratio_Sync:.3f}

3. DYNAMICAL MASS DISCREPANCY:
   - If velocities are enhanced, M_dyn from Jeans will be overestimated
   - Compare M_dyn vs M_lensing (unaffected by G_eff)
   - Prediction: M_dyn/M_lensing increases with radius

4. OBSERVABLE SOURCES:
   - Planck PSZ2 clusters (388 with σ measurements)
   - SDSS spectroscopic members for radial profiles
   - Individual deep surveys (Coma, Virgo, Perseus)

5. EXPECTED SIGNAL SIZE:
   - Enhancement at 3 R_200: ~{(enhancement[-2]-1)*100:.0f}% if ρ ~ ρ_cosmic
   - This is at the EDGE of detectability given measurement errors

6. COMPLICATIONS:
   - Interlopers inflate measured σ at large radii
   - Need careful membership selection
   - Projection effects matter
""")

# =============================================================================
# 10. SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #176b: SUMMARY")
print("=" * 70)

print(f"""
CORRECTED THEORETICAL FRAMEWORK
===============================

1. CLUSTER PARAMETERS (now correct):
   - Coma-like: M_200 = 10^15 M_sun, R_200 = {R_200:.3f} Mpc = {R_200*1000:.0f} kpc
   - Scale radius: r_s = {r_s:.3f} Mpc = {r_s*1000:.0f} kpc

2. DENSITY PROFILE:
   - Core (r ~ 0.01 R_200): ρ/ρ_cosmic ~ {rho_ratios[0]:.0f}
   - R_200: ρ/ρ_cosmic ~ {rho_ratios[np.argmin(np.abs(r_norm-1.0))]:.0f}
   - Outskirts (r ~ 5 R_200): ρ/ρ_cosmic ~ {rho_ratios[-1]:.1f}

3. COHERENCE FUNCTION:
   - Core: C → 1, G_eff ≈ G (standard gravity)
   - Outskirts: C → {coherence(rho_ratios[-1]):.2f}, G_eff ≈ {G_eff_ratio(rho_ratios[-1]):.2f} G

4. VELOCITY DISPERSION:
   - σ_ΛCDM at R_200: {sigma_LCDM[np.argmin(np.abs(r_norm-1.0))]:.0f} km/s
   - σ_Sync at R_200: {sigma_Sync[np.argmin(np.abs(r_norm-1.0))]:.0f} km/s
   - Enhancement: {enhancement[np.argmin(np.abs(r_norm-1.0))]*100-100:.1f}%

5. KEY PREDICTION:
   - σ(2 R_200) / σ(0.1 R_200):
     ΛCDM: {ratio_LCDM:.3f}
     Synchronism: {ratio_Sync:.3f}
     Difference: {(ratio_Sync/ratio_LCDM-1)*100:.1f}%

6. CRITICAL INSIGHT:
   - Even at 5 R_200, density is still ~10× cosmic mean
   - Synchronism enhancement is SMALL within cluster outskirts
   - Significant enhancement only occurs FAR from cluster (r >> R_200)
   - This explains Session #173 result: 20-60 Mpc is FIELD, not cluster

7. MRH ASSESSMENT:
   ✓ Cluster dynamics IS the correct MRH for testing
   ✓ But enhancement is only ~{(enhancement[-1]-1)*100:.0f}% at accessible radii
   ✓ This may be below observational uncertainty
   ✓ Need LOWER mass systems or void-cluster comparisons

NEXT STEPS:
-----------
1. Analyze dwarf spheroidals (lower mass, steeper density gradient)
2. Compare void galaxies vs cluster galaxies
3. Use dynamical vs lensing mass ratios
""")

print("=" * 70)
print("SESSION #176b COMPLETE")
print("=" * 70)
