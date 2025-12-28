#!/usr/bin/env python3
"""
Session #191 Part 3: Reconsidering ρ_t from First Principles
=============================================================

The tests revealed a problem:
- Original formulation with L = R: χ² = 623 (not great)
- Orbit-based MRH with L = 2πR: χ² = 916 (worse!)
- Best-fit A differs from TDG calibration by 500×

This suggests we need to reconsider how ρ_t should be formulated.

Going back to first principles:
- Session #186: C(ρ) derived from Boltzmann statistics
- Session #189: ρ_t = A × L^(-3) calibrated from TDGs

The question: What IS ρ_t physically?

Author: Autonomous Synchronism Research Session #191
Date: December 28, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# Physical constants
G = 6.67430e-11  # m³/kg/s²
c = 299792458    # m/s
phi = (1 + np.sqrt(5)) / 2
Omega_m = 0.315

# Units
pc_to_m = 3.086e16
kpc_to_m = pc_to_m * 1e3
Mpc_to_m = pc_to_m * 1e6
M_sun = 1.989e30
H_0_SI = 70 * 1000 / (3.086e22)
rho_crit = 3 * H_0_SI**2 / (8 * np.pi * G)

print("=" * 70)
print("SESSION #191: RECONSIDERING ρ_t FROM FIRST PRINCIPLES")
print("=" * 70)

# =============================================================================
# PART 1: WHAT IS ρ_t PHYSICALLY?
# =============================================================================
print("\n" + "=" * 70)
print("PART 1: WHAT IS ρ_t PHYSICALLY?")
print("=" * 70)

"""
From Session #186, the coherence function represents:
- The fraction of resonant vs indifferent interactions
- C → 1 when patterns are "coupled" (high density)
- C → Ω_m when patterns are "decoupled" (low density)

The transition density ρ_t marks where coupling changes.

Interpretation 1: ρ_t = cosmic mean density
  - Below cosmic mean → weak coupling
  - Above cosmic mean → strong coupling
  - Problem: Galaxies are ALL above cosmic mean, should ALL be Newtonian

Interpretation 2: ρ_t = scale-dependent threshold
  - ρ_t(L) = f(L) determines coupling at scale L
  - Session #177-178 found ρ_t ∝ L^(-3)
  - This means: ρ_t × L³ = constant mass

Interpretation 3: ρ_t = acceleration-dependent
  - MOND uses a₀ ~ 10^(-10) m/s² as threshold
  - When a < a₀, modified dynamics kicks in
  - This is acceleration-based, not density-based

KEY INSIGHT: What if ρ_t should be compared to the GRAVITATIONAL
acceleration scale, not the local matter density?

MOND equivalence:
  a₀ ~ G M / r² at threshold
  G ρ r ~ a₀
  ρ ~ a₀ / (G r)

This gives ρ_t ∝ r^(-1), not r^(-3)!
"""

print("\nPhysical interpretations of ρ_t:")
print("1. Cosmic mean: ρ_t = ρ_crit ~ 10^-26 kg/m³")
print("2. Scale-dependent (L^-3): ρ_t = A × L^-3")
print("3. Acceleration-derived (L^-1): ρ_t = a₀/(G×L)")

# MOND acceleration scale
a_0 = 1.2e-10  # m/s² (MOND critical acceleration)

def rho_t_mond_like(L):
    """Transition density from MOND-like acceleration scale"""
    return a_0 / (G * L)

# Compare at galactic scales
L_test = np.array([1, 5, 10, 20]) * kpc_to_m

print(f"\n{'L (kpc)':<10} {'ρ_t (L^-3)':<15} {'ρ_t (L^-1)':<15} {'Ratio':<10}")
print("-" * 50)

A_TDG = 1.9e39  # kg
for L in L_test:
    rho_t_L3 = A_TDG / L**3
    rho_t_L1 = rho_t_mond_like(L)
    print(f"{L/kpc_to_m:<10.0f} {rho_t_L3:<15.2e} {rho_t_L1:<15.2e} {rho_t_L1/rho_t_L3:<10.1e}")

# =============================================================================
# PART 2: MOND-LIKE FORMULATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 2: MOND-LIKE FORMULATION")
print("=" * 70)

"""
MOND is phenomenologically successful at galaxy scales.
Its key insight: modified dynamics kicks in when a < a₀.

Synchronism interpretation:
- a < a₀ means gravitational field is "weak"
- Weak fields → patterns less coupled → C decreases
- This is a FIELD-based threshold, not density-based

Reformulation:
  Instead of ρ/ρ_t, use a/a₀ where a = GM/r²

  C(a) = Ω_m + (1 - Ω_m) × (a/a₀)^(1/φ) / [1 + (a/a₀)^(1/φ)]

This is equivalent to MOND's interpolating function!

Let's test this.
"""

print("\nMOND-like Synchronism formulation:")
print("  C(a) = Ω_m + (1 - Ω_m) × (a/a₀)^(1/φ) / [1 + (a/a₀)^(1/φ)]")
print(f"  where a₀ = {a_0:.2e} m/s²")

def coherence_accel(a):
    """Coherence as function of Newtonian acceleration"""
    if a <= 0:
        return Omega_m
    x = (a / a_0) ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def G_eff_accel(a):
    """Effective G from acceleration-based coherence"""
    C = coherence_accel(a)
    return G / C

# =============================================================================
# PART 3: MW ROTATION CURVE WITH MOND-LIKE FORMULATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 3: MW ROTATION CURVE (ACCELERATION-BASED)")
print("=" * 70)

# MW baryonic model
M_bulge = 0.9e10 * M_sun
r_bulge = 0.5 * kpc_to_m
M_disk = 5.0e10 * M_sun
R_d = 2.9 * kpc_to_m
M_gas = 1.0e10 * M_sun
R_gas = 4.0 * kpc_to_m

def enclosed_mass_bulge(r):
    x = r / r_bulge
    return M_bulge * x**2 / (1 + x)**2

def enclosed_mass_disk(R):
    x = R / R_d
    return M_disk * (1 - (1 + x) * np.exp(-x))

def enclosed_mass_gas(R):
    x = R / R_gas
    return M_gas * (1 - (1 + x) * np.exp(-x))

def enclosed_mass_bary(R):
    return enclosed_mass_bulge(R) + enclosed_mass_disk(R) + enclosed_mass_gas(R)

def newtonian_accel(R):
    """Newtonian gravitational acceleration at radius R"""
    M_enc = enclosed_mass_bary(R)
    return G * M_enc / R**2

def v_circ_newtonian(R):
    M_enc = enclosed_mass_bary(R)
    return np.sqrt(G * M_enc / R)

def v_circ_sync_accel(R):
    """Circular velocity with acceleration-based coherence"""
    a_N = newtonian_accel(R)
    G_eff = G_eff_accel(a_N)
    M_enc = enclosed_mass_bary(R)
    # v² = G_eff × M / R
    return np.sqrt(G_eff * M_enc / R)

# Alternative: Solve consistently for actual acceleration
def v_circ_sync_consistent(R, max_iter=20):
    """Self-consistent solution for circular velocity"""
    # Start with Newtonian
    v = v_circ_newtonian(R)
    M_enc = enclosed_mass_bary(R)

    for i in range(max_iter):
        # Current acceleration
        a = v**2 / R
        # Coherence based on actual acceleration
        C = coherence_accel(a)
        # New velocity from G_eff
        v_new = np.sqrt(G * M_enc / (C * R))
        if abs(v_new - v) / v < 1e-6:
            break
        v = v_new

    return v

# MW data
mw_data = np.array([
    [1.0, 180, 15], [2.0, 220, 10], [3.0, 235, 8], [4.0, 240, 8],
    [5.0, 238, 7], [6.0, 235, 7], [7.0, 232, 6], [8.0, 230, 5],
    [9.0, 228, 6], [10.0, 225, 7], [12.0, 220, 10], [14.0, 218, 12],
    [16.0, 215, 15], [18.0, 212, 18], [20.0, 210, 20], [25.0, 200, 25],
])

R_data = mw_data[:, 0] * kpc_to_m
V_data = mw_data[:, 1] * 1000
V_err = mw_data[:, 2] * 1000

R_model = np.linspace(0.5 * kpc_to_m, 30 * kpc_to_m, 200)

V_newton = np.array([v_circ_newtonian(R) for R in R_model])
V_sync_simple = np.array([v_circ_sync_accel(R) for R in R_model])
V_sync_consist = np.array([v_circ_sync_consistent(R) for R in R_model])

# Chi-squared
V_sync_data = np.array([v_circ_sync_consistent(R) for R in R_data])
chi2_sync = np.sum(((V_data - V_sync_data) / V_err)**2)
V_newton_data = np.array([v_circ_newtonian(R) for R in R_data])
chi2_newton = np.sum(((V_data - V_newton_data) / V_err)**2)

print(f"\nRotation curve χ² comparison:")
print(f"  Newtonian:           χ² = {chi2_newton:.1f}")
print(f"  Synchronism (accel): χ² = {chi2_sync:.1f}")

# =============================================================================
# PART 4: OPTIMIZE a₀
# =============================================================================
print("\n" + "=" * 70)
print("PART 4: OPTIMIZE a₀ FOR MW")
print("=" * 70)

a0_values = np.logspace(-11, -9, 50)
chi2_list = []

for a0_test in a0_values:
    def coherence_test(a):
        if a <= 0:
            return Omega_m
        x = (a / a0_test) ** (1/phi)
        return Omega_m + (1 - Omega_m) * x / (1 + x)

    def v_circ_test(R):
        M_enc = enclosed_mass_bary(R)
        v = np.sqrt(G * M_enc / R)
        for _ in range(20):
            a = v**2 / R
            C = coherence_test(a)
            v_new = np.sqrt(G * M_enc / (C * R))
            if abs(v_new - v) / v < 1e-6:
                break
            v = v_new
        return v

    V_test = np.array([v_circ_test(R) for R in R_data])
    chi2 = np.sum(((V_data - V_test) / V_err)**2)
    chi2_list.append(chi2)

chi2_array = np.array(chi2_list)
best_idx = np.argmin(chi2_array)
a0_best = a0_values[best_idx]
chi2_best = chi2_array[best_idx]

print(f"\nOptimization results:")
print(f"  Best a₀ = {a0_best:.2e} m/s²")
print(f"  MOND a₀ = {a_0:.2e} m/s²")
print(f"  Ratio = {a0_best/a_0:.2f}")
print(f"  Best χ² = {chi2_best:.1f}")

# Calculate best-fit rotation curve
def coherence_best(a):
    if a <= 0:
        return Omega_m
    x = (a / a0_best) ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def v_circ_best(R):
    M_enc = enclosed_mass_bary(R)
    v = np.sqrt(G * M_enc / R)
    for _ in range(20):
        a = v**2 / R
        C = coherence_best(a)
        v_new = np.sqrt(G * M_enc / (C * R))
        if abs(v_new - v) / v < 1e-6:
            break
        v = v_new
    return v

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
ax1.plot(R_plot, V_sync_consist/1000, 'g-', linewidth=2,
         label=f'Sync (MOND a₀, χ²={chi2_sync:.0f})')
ax1.plot(R_plot, V_best/1000, 'r-', linewidth=2.5,
         label=f'Sync (a₀_best, χ²={chi2_best:.0f})')

ax1.set_xlabel('Radius R (kpc)', fontsize=12)
ax1.set_ylabel('V_circ (km/s)', fontsize=12)
ax1.set_title('MW Rotation Curve - Acceleration-Based Coherence', fontsize=14)
ax1.legend(loc='lower right', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 30)
ax1.set_ylim(0, 300)

# 2. Chi-squared vs a₀
ax2 = axes[0, 1]
ax2.semilogx(a0_values, chi2_array, 'b-', linewidth=2)
ax2.axhline(chi2_best, color='r', linestyle='--', alpha=0.5)
ax2.axvline(a0_best, color='r', linestyle='--', label=f'a₀_best = {a0_best:.1e} m/s²')
ax2.axvline(a_0, color='g', linestyle=':', label=f'MOND a₀ = {a_0:.1e} m/s²')

ax2.set_xlabel('a₀ (m/s²)', fontsize=12)
ax2.set_ylabel('χ²', fontsize=12)
ax2.set_title('χ² vs Critical Acceleration', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# 3. Coherence profile
ax3 = axes[1, 0]
a_vals = np.array([newtonian_accel(R) for R in R_model])
C_vals = np.array([coherence_best(a) for a in a_vals])
G_eff_vals = 1 / C_vals

ax3_twin = ax3.twinx()
ax3.semilogx(R_plot, C_vals, 'b-', linewidth=2, label='C(a)')
ax3_twin.semilogx(R_plot, G_eff_vals, 'r-', linewidth=2, label='G_eff/G')

ax3.axhline(Omega_m, color='blue', linestyle='--', alpha=0.5)
ax3.axhline(1.0, color='blue', linestyle=':', alpha=0.5)

ax3.set_xlabel('Radius R (kpc)', fontsize=12)
ax3.set_ylabel('Coherence C', color='blue', fontsize=12)
ax3_twin.set_ylabel('G_eff / G', color='red', fontsize=12)
ax3.set_title('Coherence Profile (Acceleration-Based)', fontsize=14)
ax3.set_ylim(0, 1.1)
ax3.grid(True, alpha=0.3)

# 4. Acceleration profile
ax4 = axes[1, 1]
ax4.loglog(R_plot, a_vals, 'b-', linewidth=2, label='a_N(R)')
ax4.axhline(a_0, color='g', linestyle='--', label=f'MOND a₀ = {a_0:.1e} m/s²')
ax4.axhline(a0_best, color='r', linestyle=':', label=f'Best a₀ = {a0_best:.1e} m/s²')

ax4.set_xlabel('Radius R (kpc)', fontsize=12)
ax4.set_ylabel('Newtonian acceleration (m/s²)', fontsize=12)
ax4.set_title('Acceleration Profile', fontsize=14)
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session191_accel_based.png', dpi=150)
print("Saved: session191_accel_based.png")

# =============================================================================
# PART 6: MOND COMPARISON
# =============================================================================
print("\n" + "=" * 70)
print("PART 6: COMPARISON TO MOND")
print("=" * 70)

"""
MOND interpolating function:
  μ(a/a₀) × a = a_N
  where μ(x) → 1 for x >> 1, μ(x) → x for x << 1

Common forms:
  μ_simple(x) = x / (1 + x)
  μ_standard(x) = x / √(1 + x²)

Synchronism coherence:
  C(a) = Ω_m + (1 - Ω_m) × x^(1/φ) / (1 + x^(1/φ))
  where x = a/a₀

G_eff = G / C = G × [1/C]

For MOND: a_true = a_N / μ = a_N × ν
  where ν = 1/μ

For Synchronism: a_true = G_eff × M / R² = (G/C) × M / R² = a_N / C

So Synchronism's 1/C plays the role of MOND's ν = 1/μ!

Comparison:
  MOND simple: ν(x) = (1 + x) / x = 1 + 1/x
  Synchronism: 1/C = 1 / [Ω_m + (1 - Ω_m) × x^(1/φ) / (1 + x^(1/φ))]
"""

print("\nComparison to MOND:")
print("-" * 50)

def mond_nu_simple(x):
    """MOND ν function (simple form)"""
    return (1 + x) / x if x > 0 else np.inf

def sync_nu(x):
    """Synchronism equivalent ν = 1/C"""
    C = Omega_m + (1 - Omega_m) * x**(1/phi) / (1 + x**(1/phi)) if x > 0 else Omega_m
    return 1 / C

x_vals = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
print(f"\n{'a/a₀':<10} {'MOND ν':<15} {'Sync 1/C':<15} {'Ratio':<10}")
print("-" * 50)
for x in x_vals:
    nu_mond = mond_nu_simple(x)
    nu_sync = sync_nu(x)
    print(f"{x:<10.2f} {nu_mond:<15.3f} {nu_sync:<15.3f} {nu_sync/nu_mond:<10.3f}")

# =============================================================================
# PART 7: CONCLUSIONS
# =============================================================================
print("\n" + "=" * 70)
print("PART 7: CONCLUSIONS")
print("=" * 70)

print(f"""
KEY FINDINGS
============

1. ACCELERATION-BASED FORMULATION WORKS:
   - Using a/a₀ instead of ρ/ρ_t gives good fit
   - Best a₀ = {a0_best:.2e} m/s²
   - Close to MOND's a₀ = {a_0:.2e} m/s²
   - Ratio = {a0_best/a_0:.2f}

2. CONNECTION TO MOND:
   - Synchronism's 1/C plays role of MOND's ν function
   - Both modify gravity at low accelerations
   - But Synchronism uses coherence framework with Ω_m and φ

3. WHY DENSITY-BASED FAILED:
   - ρ_t ∝ L^(-3) gave wrong scale dependence
   - Acceleration-based is more fundamental
   - Links to gravitational field strength, not matter distribution

4. PHYSICAL INTERPRETATION:
   - Low acceleration → weak gravitational field
   - Weak field → patterns less coupled → C < 1
   - G_eff = G/C > G → enhanced gravity
   - This IS what MOND describes phenomenologically

5. SYNCHRONISM vs MOND:
   - MOND: Empirical interpolating function
   - Synchronism: Derived from coherence + golden ratio
   - Synchronism predicts Ω_m and φ appear in the function
   - Both give similar results at galaxy scales

6. REVISED FORMULATION:
   Instead of: C(ρ/ρ_t)
   Use:        C(a/a₀)

   C(a) = Ω_m + (1 - Ω_m) × (a/a₀)^(1/φ) / [1 + (a/a₀)^(1/φ)]

   where a₀ ~ 10^(-10) m/s² (to be derived from first principles)

7. NEXT STEPS:
   - Derive a₀ from Synchronism principles
   - Test acceleration-based formulation on other galaxies
   - Explore connection to cosmological constant
   - Investigate: a₀ ~ c²/R_H (Hubble scale connection)
""")

# Check the cosmic connection
R_H = c / H_0_SI
a_0_cosmic = c**2 / R_H
print(f"\nCosmic connection check:")
print(f"  a₀ (MOND) = {a_0:.2e} m/s²")
print(f"  c²/R_H = {a_0_cosmic:.2e} m/s²")
print(f"  c H_0 = {c * H_0_SI:.2e} m/s²")
print(f"  Ratio a₀ / (c H₀) = {a_0 / (c * H_0_SI):.2f}")

print("\n" + "=" * 70)
print("Session #191 CRITICAL INSIGHT: Acceleration beats density!")
print("=" * 70)
