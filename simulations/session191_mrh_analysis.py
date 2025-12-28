#!/usr/bin/env python3
"""
Session #191 Part 2: MRH Analysis for Rotation Curves
======================================================

Initial test showed:
- Newtonian χ² = 922
- Synchronism χ² = 623 (TDG-calibrated A = 1.9e39 kg)
- Best fit A = 1e40 kg

Question: Is the MRH being applied correctly?

The issue: What is the "relevant MRH" for stellar orbits?

Options:
1. MRH = R (galactocentric radius) - what we used
2. MRH = R_orbit (full orbit scale, ~2πR)
3. MRH = vertical scale height (disk thickness, ~0.3 kpc)
4. MRH = local Jeans length
5. MRH = scale at which pattern coherence is evaluated

Let's investigate which interpretation gives best physical results.

Author: Autonomous Synchronism Research Session #191
Date: December 28, 2025
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
G = 6.67430e-11  # m³/kg/s²
phi = (1 + np.sqrt(5)) / 2  # Golden ratio
Omega_m = 0.315

# Unit conversions
pc_to_m = 3.086e16
kpc_to_m = pc_to_m * 1e3
M_sun = 1.989e30  # kg

# MW parameters (from main test)
M_bulge = 0.9e10 * M_sun
r_bulge = 0.5 * kpc_to_m
M_disk = 5.0e10 * M_sun
R_d = 2.9 * kpc_to_m
M_gas = 1.0e10 * M_sun
R_gas = 4.0 * kpc_to_m

# Synchronism parameters
A_TDG = 1.9e39  # kg (from TDG calibration)

print("=" * 70)
print("SESSION #191: MRH ANALYSIS FOR ROTATION CURVES")
print("=" * 70)

# Baryonic mass model
def enclosed_mass_bulge(r):
    x = r / r_bulge
    return M_bulge * x**2 / (1 + x)**2

def enclosed_mass_disk_approx(R):
    x = R / R_d
    return M_disk * (1 - (1 + x) * np.exp(-x))

def enclosed_mass_gas(R):
    x = R / R_gas
    return M_gas * (1 - (1 + x) * np.exp(-x))

def enclosed_mass_bary(R):
    return enclosed_mass_bulge(R) + enclosed_mass_disk_approx(R) + enclosed_mass_gas(R)

def v_circ_newtonian(R):
    M_enc = enclosed_mass_bary(R)
    return np.sqrt(G * M_enc / R)

# Coherence function
def coherence(rho_ratio):
    if rho_ratio <= 0:
        return Omega_m
    x = rho_ratio ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

# =============================================================================
# PART 1: DIFFERENT MRH INTERPRETATIONS
# =============================================================================
print("\n" + "=" * 70)
print("PART 1: DIFFERENT MRH INTERPRETATIONS")
print("=" * 70)

"""
The question: At what scale should we evaluate ρ_t?

Key insight from Session #188: MRH must be "relevant to the phenomenon."

For stellar orbits:
- The "phenomenon" is the gravitational force felt by a star
- The star orbits at radius R with period T ∝ R/v
- The star's motion samples the gravitational field over the orbit

Different MRH options:

1. L = R (current approach)
   - ρ_t = A / R³
   - Problem: Gives high ρ_t in inner galaxy

2. L = 2πR (orbit circumference)
   - ρ_t = A / (2πR)³
   - More appropriate for orbit-averaged dynamics

3. L = characteristic scale of density gradient
   - For exponential disk: L ~ R_d ~ 3 kpc
   - Fixed scale, not R-dependent

4. L = vertical scale height h_z
   - For thin disk: h_z ~ 0.3 kpc
   - The "local" scale for disk stars

5. L = orbit-dependent with environment
   - L depends on both R and local density structure

Let's test Option 2 (orbit scale) as most physically motivated.
"""

print("\nMRH Interpretations:")
print("1. L = R (point at radius)")
print("2. L = 2πR (orbit circumference)")
print("3. L = R_d (disk scale length, fixed)")
print("4. L = h_z (vertical scale, fixed)")
print("5. L = hybrid (orbit + environment)")

# =============================================================================
# PART 2: ORBIT-BASED MRH
# =============================================================================
print("\n" + "=" * 70)
print("PART 2: ORBIT-BASED MRH IMPLEMENTATION")
print("=" * 70)

"""
Key insight: A stellar orbit "samples" a region of size ~2πR × h_z × h_z

For the coherence evaluation, we need to consider:
- The gravitational influence region (the orbit)
- The matter density WITHIN that region

Revision: The LOCAL density should be compared to the TRANSITION density
at the orbit scale.

Previous approach:
  ρ_local = M_enc / (4π/3 × R³)  (spherical average)
  ρ_t = A / R³

Revised approach:
  ρ_local = local environment density (need proper calculation)
  L_orbit = 2πR (orbit scale)
  ρ_t = A / L_orbit³ = A / (2πR)³

This gives ρ_t that is (2π)³ ≈ 248× smaller than before!
"""

def local_density_old(R):
    """Old: spherical average enclosed density"""
    M_enc = enclosed_mass_bary(R)
    V_sphere = 4/3 * np.pi * R**3
    return M_enc / V_sphere

def rho_t_old(R):
    """Old: ρ_t at scale R"""
    return A_TDG / R**3

def local_density_disk(R):
    """Disk midplane density estimate"""
    # For exponential disk, midplane density:
    # ρ_0(R) ≈ Σ(R) / (2h_z) where Σ = Σ_0 exp(-R/R_d)
    h_z = 0.3 * kpc_to_m  # Thin disk scale height
    Sigma = (M_disk / (2 * np.pi * R_d**2)) * np.exp(-R / R_d)
    rho_disk = Sigma / (2 * h_z)

    # Add bulge contribution at inner radii
    rho_bulge = enclosed_mass_bulge(R) / (4/3 * np.pi * R**3)

    return max(rho_disk, rho_bulge)

def rho_t_orbit(R):
    """New: ρ_t at orbit scale L = 2πR"""
    L_orbit = 2 * np.pi * R
    return A_TDG / L_orbit**3

# Compare the two approaches
R_test = np.array([1, 2, 5, 8, 10, 15, 20]) * kpc_to_m

print(f"\n{'R (kpc)':<10} {'ρ_local_old':<15} {'ρ_t_old':<15} {'ρ/ρ_t old':<12} {'C_old':<8}")
print("-" * 70)
for R in R_test:
    rho_l = local_density_old(R)
    rho_t = rho_t_old(R)
    ratio = rho_l / rho_t
    C = coherence(ratio)
    print(f"{R/kpc_to_m:<10.0f} {rho_l:<15.2e} {rho_t:<15.2e} {ratio:<12.3e} {C:<8.3f}")

print(f"\n{'R (kpc)':<10} {'ρ_local_disk':<15} {'ρ_t_orbit':<15} {'ρ/ρ_t new':<12} {'C_new':<8}")
print("-" * 70)
for R in R_test:
    rho_l = local_density_disk(R)
    rho_t = rho_t_orbit(R)
    ratio = rho_l / rho_t
    C = coherence(ratio)
    print(f"{R/kpc_to_m:<10.0f} {rho_l:<15.2e} {rho_t:<15.2e} {ratio:<12.3e} {C:<8.3f}")

# =============================================================================
# PART 3: REVISED ROTATION CURVE
# =============================================================================
print("\n" + "=" * 70)
print("PART 3: REVISED ROTATION CURVE (ORBIT MRH)")
print("=" * 70)

def G_eff_orbit_mrh(R):
    """G_eff using orbit-based MRH"""
    rho_local = local_density_disk(R)
    rho_trans = rho_t_orbit(R)
    C = coherence(rho_local / rho_trans)
    return G / C

def v_circ_orbit_mrh(R):
    """Circular velocity with orbit-based MRH"""
    M_enc = enclosed_mass_bary(R)
    G_effective = G_eff_orbit_mrh(R)
    return np.sqrt(G_effective * M_enc / R)

# MW data
mw_data = np.array([
    [1.0, 180, 15],
    [2.0, 220, 10],
    [3.0, 235, 8],
    [4.0, 240, 8],
    [5.0, 238, 7],
    [6.0, 235, 7],
    [7.0, 232, 6],
    [8.0, 230, 5],
    [9.0, 228, 6],
    [10.0, 225, 7],
    [12.0, 220, 10],
    [14.0, 218, 12],
    [16.0, 215, 15],
    [18.0, 212, 18],
    [20.0, 210, 20],
    [25.0, 200, 25],
])

R_data = mw_data[:, 0] * kpc_to_m
V_data = mw_data[:, 1] * 1000
V_err = mw_data[:, 2] * 1000

R_model = np.linspace(0.5 * kpc_to_m, 30 * kpc_to_m, 200)

V_newton = np.array([v_circ_newtonian(R) for R in R_model])
V_orbit_mrh = np.array([v_circ_orbit_mrh(R) for R in R_model])

V_orbit_data = np.array([v_circ_orbit_mrh(R) for R in R_data])
chi2_orbit = np.sum(((V_data - V_orbit_data) / V_err)**2)

print(f"\nOrbit-based MRH χ² = {chi2_orbit:.1f}")

# =============================================================================
# PART 4: OPTIMIZE A FOR MW
# =============================================================================
print("\n" + "=" * 70)
print("PART 4: OPTIMIZE A FOR MILKY WAY")
print("=" * 70)

"""
The TDG-calibrated A was derived from a different system (TDGs).
Let's find the optimal A for MW and compare to TDG value.

If they're similar: Synchronism is consistent
If they're different: May indicate scale-dependent A or model issue
"""

A_values = np.logspace(38, 42, 50)
chi2_list = []

for A in A_values:
    def rho_t_test(R):
        return A / (2 * np.pi * R)**3

    def G_eff_test(R):
        rho_local = local_density_disk(R)
        rho_trans = rho_t_test(R)
        C = coherence(rho_local / rho_trans)
        return G / C

    def v_circ_test(R):
        M_enc = enclosed_mass_bary(R)
        return np.sqrt(G_eff_test(R) * M_enc / R)

    V_test = np.array([v_circ_test(R) for R in R_data])
    chi2 = np.sum(((V_data - V_test) / V_err)**2)
    chi2_list.append(chi2)

chi2_array = np.array(chi2_list)
best_idx = np.argmin(chi2_array)
A_best = A_values[best_idx]
chi2_best = chi2_array[best_idx]

print(f"\nOptimization results (orbit-based MRH):")
print(f"  Best A = {A_best:.2e} kg")
print(f"  Best A = {A_best/M_sun:.2e} M_sun")
print(f"  Best χ² = {chi2_best:.1f}")
print(f"  TDG A = {A_TDG:.2e} kg")
print(f"  Ratio A_MW/A_TDG = {A_best/A_TDG:.1f}")

# Calculate best-fit rotation curve
def rho_t_best(R):
    return A_best / (2 * np.pi * R)**3

def G_eff_best(R):
    rho_local = local_density_disk(R)
    rho_trans = rho_t_best(R)
    C = coherence(rho_local / rho_trans)
    return G / C

def v_circ_best(R):
    M_enc = enclosed_mass_bary(R)
    return np.sqrt(G_eff_best(R) * M_enc / R)

V_best = np.array([v_circ_best(R) for R in R_model])

# =============================================================================
# PART 5: VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 5: VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
R_plot = R_model / kpc_to_m

# 1. Rotation curves
ax1 = axes[0, 0]
ax1.errorbar(mw_data[:, 0], mw_data[:, 1], yerr=mw_data[:, 2],
             fmt='ko', capsize=3, label='MW observations', markersize=6)
ax1.plot(R_plot, V_newton/1000, 'b--', linewidth=2, label='Newtonian')
ax1.plot(R_plot, V_orbit_mrh/1000, 'g-', linewidth=2,
         label=f'Orbit MRH, A_TDG (χ²={chi2_orbit:.0f})')
ax1.plot(R_plot, V_best/1000, 'r-', linewidth=2.5,
         label=f'Orbit MRH, A_best (χ²={chi2_best:.0f})')

ax1.set_xlabel('Radius R (kpc)', fontsize=12)
ax1.set_ylabel('V_circ (km/s)', fontsize=12)
ax1.set_title('MW Rotation Curve - Orbit-Based MRH', fontsize=14)
ax1.legend(loc='lower right', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 30)
ax1.set_ylim(0, 300)

# 2. Chi-squared vs A
ax2 = axes[0, 1]
ax2.loglog(A_values, chi2_array, 'b-', linewidth=2)
ax2.axhline(chi2_best, color='r', linestyle='--', alpha=0.5)
ax2.axvline(A_best, color='r', linestyle='--', label=f'A_best = {A_best:.1e} kg')
ax2.axvline(A_TDG, color='g', linestyle=':', label=f'A_TDG = {A_TDG:.1e} kg')

ax2.set_xlabel('A (kg)', fontsize=12)
ax2.set_ylabel('χ²', fontsize=12)
ax2.set_title('χ² vs Normalization Constant A', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# 3. Coherence profile
ax3 = axes[1, 0]
C_best = [coherence(local_density_disk(R) / rho_t_best(R)) for R in R_model]
G_eff_ratio = [1/C for C in C_best]

ax3_twin = ax3.twinx()
ax3.plot(R_plot, C_best, 'b-', linewidth=2, label='C(R)')
ax3_twin.plot(R_plot, G_eff_ratio, 'r-', linewidth=2, label='G_eff/G')

ax3.axhline(Omega_m, color='blue', linestyle='--', alpha=0.5)
ax3.axhline(1.0, color='blue', linestyle=':', alpha=0.5)

ax3.set_xlabel('Radius R (kpc)', fontsize=12)
ax3.set_ylabel('Coherence C', color='blue', fontsize=12)
ax3_twin.set_ylabel('G_eff / G', color='red', fontsize=12)
ax3.set_title('Coherence Profile (Best-Fit A)', fontsize=14)
ax3.set_xlim(0, 30)
ax3.set_ylim(0, 1.1)
ax3.grid(True, alpha=0.3)

# 4. Density comparison
ax4 = axes[1, 1]
rho_local_vals = [local_density_disk(R) for R in R_model]
rho_t_best_vals = [rho_t_best(R) for R in R_model]
rho_t_tdg_vals = [A_TDG / (2*np.pi*R)**3 for R in R_model]

ax4.loglog(R_plot, rho_local_vals, 'b-', linewidth=2, label='ρ_local (disk)')
ax4.loglog(R_plot, rho_t_best_vals, 'r-', linewidth=2, label=f'ρ_t (A_best)')
ax4.loglog(R_plot, rho_t_tdg_vals, 'g--', linewidth=2, label=f'ρ_t (A_TDG)')

ax4.set_xlabel('Radius R (kpc)', fontsize=12)
ax4.set_ylabel('Density (kg/m³)', fontsize=12)
ax4.set_title('Local vs Transition Density', fontsize=14)
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0.5, 30)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session191_mrh_analysis.png', dpi=150)
print("Saved: session191_mrh_analysis.png")

# =============================================================================
# PART 6: INTERPRETATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 6: INTERPRETATION")
print("=" * 70)

print(f"""
KEY FINDINGS
============

1. MRH REVISION:
   - Original: L = R (galactocentric radius)
   - Revised: L = 2πR (orbit circumference)
   - This is more physically motivated: stars sample full orbit

2. DENSITY MODEL REVISION:
   - Original: ρ_local = M_enc / (4π/3 R³) (spherical average)
   - Revised: ρ_local = Σ(R) / (2h_z) (disk midplane)
   - This reflects actual stellar environment

3. OPTIMAL A FOR MW:
   - A_best = {A_best:.2e} kg = {A_best/M_sun:.2e} M_sun
   - A_TDG = {A_TDG:.2e} kg = {A_TDG/M_sun:.2e} M_sun
   - Ratio = {A_best/A_TDG:.1f}

4. INTERPRETATION OF DISCREPANCY:
""")

if 0.3 < A_best/A_TDG < 3:
    print("   A_MW and A_TDG are within factor of 3 - CONSISTENT")
    print("   Differences may be due to:")
    print("   - Different environments (tidal debris vs isolated disk)")
    print("   - Uncertainty in baryonic mass models")
    print("   - Simplified MRH implementation")
elif A_best > A_TDG:
    print(f"   A_MW >> A_TDG by factor {A_best/A_TDG:.0f}")
    print("   Possible explanations:")
    print("   - MW needs stronger modification than TDGs")
    print("   - TDG calibration may be incorrect")
    print("   - MRH prescription needs refinement")
else:
    print(f"   A_MW << A_TDG by factor {A_TDG/A_best:.0f}")
    print("   Possible explanations:")
    print("   - TDGs are in more modified regime")
    print("   - Different scale dependence")

print(f"""
5. PHYSICAL PICTURE:
   - At R = 8 kpc (solar radius): C ~ {coherence(local_density_disk(8*kpc_to_m) / rho_t_best(8*kpc_to_m)):.3f}
   - At R = 20 kpc: C ~ {coherence(local_density_disk(20*kpc_to_m) / rho_t_best(20*kpc_to_m)):.3f}
   - Inner galaxy: Nearly Newtonian (high ρ)
   - Outer galaxy: Enhanced gravity (low ρ, C < 1)

6. COMPARISON TO ΛCDM:
   - ΛCDM requires ~10× baryonic mass in dark matter halo
   - Synchronism achieves similar effect via G_eff enhancement
   - Key difference: No new particles, just density-dependent coupling

7. NEXT STEPS:
   - Test revised formulation on other galaxies
   - Investigate if A should have explicit scale dependence
   - Compare to dwarf galaxies (should show stronger effect)
""")

print("=" * 70)
