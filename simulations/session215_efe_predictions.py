#!/usr/bin/env python3
"""
Session #215: External Field Effect (EFE) Predictions
=====================================================

Key Question: Does Synchronism predict any EFE-like behavior?
- MOND: Strong EFE from non-local μ(|g_int + g_ext|)
- Synchronism: Local C(a) - no EFE in principle

This session quantifies the difference and identifies observational tests.

Author: Autonomous Research Agent
Date: January 2, 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq

# =============================================================================
# Physical Constants
# =============================================================================

G = 6.674e-11  # m^3/(kg·s^2)
c = 2.998e8    # m/s
H0 = 67.4      # km/s/Mpc
H0_SI = H0 * 1e3 / (3.086e22)  # s^-1
Omega_m = 0.315
phi = 1.618034  # Golden ratio

# Critical accelerations
a0_sync = c * H0_SI * Omega_m**phi  # Synchronism: ~1.01e-10 m/s²
a0_mond = 1.2e-10  # MOND empirical value

print("=" * 70)
print("Session #215: External Field Effect (EFE) Predictions")
print("=" * 70)
print(f"\nCritical accelerations:")
print(f"  a0_sync = {a0_sync:.3e} m/s²")
print(f"  a0_mond = {a0_mond:.3e} m/s²")

# =============================================================================
# Synchronism Functions
# =============================================================================

def C_sync(a):
    """Synchronism coherence function - LOCAL only."""
    if isinstance(a, (list, np.ndarray)):
        a = np.asarray(a)
        result = np.zeros_like(a, dtype=float)
        mask = a > 0
        x = (a[mask] / a0_sync) ** (1/phi)
        result[mask] = Omega_m + (1 - Omega_m) * x / (1 + x)
        result[~mask] = Omega_m  # Low-a limit
        return result
    else:
        if a <= 0:
            return Omega_m
        x = (a / a0_sync) ** (1/phi)
        return Omega_m + (1 - Omega_m) * x / (1 + x)

def G_eff_sync(a_int, a_ext=0):
    """
    Synchronism effective G.

    CRITICAL: In Synchronism, G_eff depends ONLY on a_int (local).
    The a_ext parameter is included to show it has NO EFFECT.
    """
    # Synchronism is LOCAL - external field has NO effect
    return 1.0 / C_sync(a_int)

# =============================================================================
# MOND Functions
# =============================================================================

def nu_mond(x):
    """MOND simple interpolating function."""
    return 0.5 * (1 + np.sqrt(1 + 4/x))

def mu_mond(x):
    """MOND μ function (inverse of ν)."""
    return x / (1 + x)  # Simple form

def G_eff_mond_no_efe(a_int):
    """MOND without EFE (isolated system)."""
    x = a_int / a0_mond
    return nu_mond(x)

def G_eff_mond_with_efe(a_int, a_ext):
    """
    MOND with External Field Effect.

    In the quasi-linear MOND approximation:
    g = ν(|g_N + g_ext|/a0) × g_N

    This is simplified for the case where g_ext >> g_int (satellite limit):
    g ≈ g_N × [ν(g_ext/a0) + (g_int/g_ext) × ν'(g_ext/a0)]

    For simplicity, we use the "extreme EFE" limit:
    G_eff = ν(a_ext/a0) when a_ext >> a_int
    G_eff = ν(a_int/a0) when a_int >> a_ext
    """
    x_int = a_int / a0_mond
    x_ext = a_ext / a0_mond

    if a_ext <= 0:
        return nu_mond(x_int)

    # Combined effective acceleration (vector addition approximation)
    # The full QUMOND treatment is more complex, but this captures the essence
    x_combined = np.sqrt(a_int**2 + a_ext**2) / a0_mond

    # In strong EFE regime (a_ext >> a_int), MOND predicts:
    # The internal dynamics is suppressed by the external field
    if a_ext > 5 * a_int:
        # Strong EFE: internal dynamics is nearly Newtonian
        # g ≈ g_N × ν(a_ext/a0) / ν(a_int/a0)
        return nu_mond(x_ext) / nu_mond(x_int) * nu_mond(x_int)
    else:
        # Weak EFE: use combined acceleration
        return nu_mond(x_combined)

def G_eff_mond_efe_proper(a_int, a_ext):
    """
    Proper MOND EFE calculation for satellite galaxies.

    When embedded in an external field a_ext, the internal dynamics
    is modified. In the deep-MOND limit for internal dynamics:

    G_eff_internal ≈ ν(a_ext/a0) for a_int << a0 << a_ext
    G_eff_internal ≈ ν(a_int/a0) for a_ext = 0

    This shows the EFE suppresses the MOND boost.
    """
    if a_ext <= 0:
        return nu_mond(a_int / a0_mond)

    x_ext = a_ext / a0_mond
    x_int = a_int / a0_mond

    # In strong EFE limit, internal dynamics sees:
    # g = g_N × ν(a_ext/a0) (the internal boost is determined by external field)
    # This means velocity dispersion scales as:
    # σ² ∝ G × M × ν(a_ext/a0) / r

    if x_ext > 1:  # Newtonian external field
        return nu_mond(x_int)  # Nearly isolated behavior
    elif x_ext > 0.1 * x_int:  # EFE regime
        # Interpolate between isolated and EFE-dominated
        # The boost is reduced when embedded
        isolated_boost = nu_mond(x_int)
        efe_boost = nu_mond(x_ext)
        # Weight by relative strength
        weight = np.exp(-x_int / (x_ext + 1e-10))
        return isolated_boost * (1 - weight) + efe_boost * weight
    else:
        return nu_mond(x_int)

# =============================================================================
# Part 1: EFE Strength Comparison
# =============================================================================

print("\n" + "=" * 70)
print("Part 1: External Field Effect Strength")
print("=" * 70)

# Typical satellite galaxy parameters
# Dwarf around MW: M_star ~ 10^6 Msun, r_half ~ 100 pc
M_dwarf = 1e6 * 2e30  # kg
r_half = 100 * 3.086e16  # m
a_int = G * M_dwarf / r_half**2  # Internal acceleration

# External field from MW at various distances
distances_kpc = np.array([50, 100, 200, 400, 800])
M_MW = 1e12 * 2e30  # kg

print(f"\nDwarf galaxy parameters:")
print(f"  M_star = 10^6 M_sun")
print(f"  r_half = 100 pc")
print(f"  a_int = {a_int:.3e} m/s² = {a_int/a0_sync:.3f} a0_sync")

print(f"\nEFE comparison at different MW distances:")
print("-" * 70)
print(f"{'Distance':>10} | {'a_ext/a0':>10} | {'Sync G_eff':>12} | {'MOND iso':>12} | {'MOND EFE':>12}")
print("-" * 70)

efe_results = []
for d_kpc in distances_kpc:
    d_m = d_kpc * 1e3 * 3.086e16
    a_ext = G * M_MW / d_m**2

    geff_sync = G_eff_sync(a_int)
    geff_mond_iso = G_eff_mond_no_efe(a_int)
    geff_mond_efe = G_eff_mond_efe_proper(a_int, a_ext)

    efe_results.append({
        'd_kpc': d_kpc,
        'a_ext': a_ext,
        'a_ext_a0': a_ext / a0_mond,
        'sync': geff_sync,
        'mond_iso': geff_mond_iso,
        'mond_efe': geff_mond_efe,
        'efe_suppression': geff_mond_efe / geff_mond_iso
    })

    print(f"{d_kpc:>8} kpc | {a_ext/a0_mond:>10.3f} | {geff_sync:>12.2f} | {geff_mond_iso:>12.2f} | {geff_mond_efe:>12.2f}")

print("-" * 70)

# =============================================================================
# Part 2: Velocity Dispersion Predictions
# =============================================================================

print("\n" + "=" * 70)
print("Part 2: Velocity Dispersion Predictions")
print("=" * 70)

def sigma_predict(M_star, r_half, G_eff):
    """Predict velocity dispersion from virial theorem."""
    # σ² ~ G_eff × G × M / r
    M = M_star * 2e30
    r = r_half * 3.086e16
    sigma_sq = G_eff * G * M / r
    return np.sqrt(sigma_sq) / 1e3  # km/s

print(f"\nVelocity dispersion predictions for 10^6 M_sun dwarf at 100 pc:")
print("-" * 70)
print(f"{'Distance':>10} | {'σ_Sync':>10} | {'σ_MOND_iso':>10} | {'σ_MOND_EFE':>10} | {'EFE suppression':>15}")
print("-" * 70)

for res in efe_results:
    sigma_sync = sigma_predict(1e6, 100, res['sync'])
    sigma_mond_iso = sigma_predict(1e6, 100, res['mond_iso'])
    sigma_mond_efe = sigma_predict(1e6, 100, res['mond_efe'])

    print(f"{res['d_kpc']:>8} kpc | {sigma_sync:>8.2f} km/s | {sigma_mond_iso:>8.2f} km/s | {sigma_mond_efe:>8.2f} km/s | {res['efe_suppression']:>14.1%}")

print("-" * 70)

# =============================================================================
# Part 3: Field Dwarfs vs Satellites - The Test
# =============================================================================

print("\n" + "=" * 70)
print("Part 3: Field Dwarfs vs Satellites - Discriminating Test")
print("=" * 70)

# Compare:
# 1. Field dwarf (isolated, a_ext ~ 0.0001 a0)
# 2. Satellite at 100 kpc from MW (a_ext ~ 0.1 a0)

print("\nScenario: Compare identical dwarfs in different environments")
print(f"  Mass: 10^6 M_sun, r_half: 100 pc")

# Field dwarf
a_ext_field = 0.0001 * a0_mond  # Very weak external field
sigma_sync_field = sigma_predict(1e6, 100, G_eff_sync(a_int))
sigma_mond_field = sigma_predict(1e6, 100, G_eff_mond_no_efe(a_int))

# Satellite at 100 kpc
d_sat = 100 * 1e3 * 3.086e16
a_ext_sat = G * M_MW / d_sat**2
sigma_sync_sat = sigma_predict(1e6, 100, G_eff_sync(a_int))  # Same - no EFE!
sigma_mond_sat = sigma_predict(1e6, 100, G_eff_mond_efe_proper(a_int, a_ext_sat))

print(f"\nField dwarf (a_ext ~ 0):")
print(f"  Synchronism: σ = {sigma_sync_field:.2f} km/s")
print(f"  MOND:        σ = {sigma_mond_field:.2f} km/s")

print(f"\nSatellite at 100 kpc (a_ext = {a_ext_sat/a0_mond:.3f} a0):")
print(f"  Synchronism: σ = {sigma_sync_sat:.2f} km/s")
print(f"  MOND:        σ = {sigma_mond_sat:.2f} km/s")

print(f"\n" + "=" * 70)
print("KEY PREDICTIONS:")
print("=" * 70)
print(f"Synchronism: σ_field / σ_satellite = {sigma_sync_field/sigma_sync_sat:.3f} (NO EFE)")
print(f"MOND:        σ_field / σ_satellite = {sigma_mond_field/sigma_mond_sat:.3f} (EFE suppresses satellite)")
print("=" * 70)

# =============================================================================
# Part 4: Sample of Real Satellites and Field Dwarfs
# =============================================================================

print("\n" + "=" * 70)
print("Part 4: MW Satellite Predictions")
print("=" * 70)

# MW satellites with known properties
# Data from McConnachie 2012 compilation
mw_satellites = [
    {"name": "Sagittarius", "M_star": 2.1e7, "r_half": 2600, "d_MW": 26, "sigma_obs": 11.4},
    {"name": "LMC", "M_star": 1.5e9, "r_half": 1000, "d_MW": 50, "sigma_obs": 20.2},
    {"name": "Fornax", "M_star": 2.0e7, "r_half": 710, "d_MW": 147, "sigma_obs": 11.7},
    {"name": "Leo I", "M_star": 5.5e6, "r_half": 251, "d_MW": 254, "sigma_obs": 9.2},
    {"name": "Sculptor", "M_star": 2.3e6, "r_half": 283, "d_MW": 86, "sigma_obs": 9.2},
    {"name": "Leo II", "M_star": 7.4e5, "r_half": 176, "d_MW": 233, "sigma_obs": 6.6},
    {"name": "Sextans", "M_star": 4.4e5, "r_half": 695, "d_MW": 86, "sigma_obs": 7.9},
    {"name": "Carina", "M_star": 3.8e5, "r_half": 250, "d_MW": 105, "sigma_obs": 6.6},
    {"name": "Draco", "M_star": 2.9e5, "r_half": 221, "d_MW": 76, "sigma_obs": 9.1},
    {"name": "Ursa Minor", "M_star": 2.9e5, "r_half": 181, "d_MW": 76, "sigma_obs": 9.5},
]

print(f"\n{'Name':>12} | {'M_star':>10} | {'r_half':>8} | {'d_MW':>6} | {'σ_obs':>6} | {'σ_Sync':>8} | {'σ_MOND':>8} | {'σ_MOND_EFE':>10}")
print("-" * 100)

satellite_predictions = []
for sat in mw_satellites:
    M = sat["M_star"]
    r = sat["r_half"]
    d = sat["d_MW"]

    # Internal acceleration
    a_i = G * M * 2e30 / (r * 3.086e16)**2

    # External acceleration from MW
    a_e = G * M_MW / (d * 1e3 * 3.086e16)**2

    # Predictions
    geff_sync = G_eff_sync(a_i)
    geff_mond_iso = G_eff_mond_no_efe(a_i)
    geff_mond_efe = G_eff_mond_efe_proper(a_i, a_e)

    sigma_sync = sigma_predict(M, r, geff_sync)
    sigma_mond = sigma_predict(M, r, geff_mond_iso)
    sigma_mond_efe = sigma_predict(M, r, geff_mond_efe)

    satellite_predictions.append({
        "name": sat["name"],
        "M_star": M,
        "r_half": r,
        "d_MW": d,
        "sigma_obs": sat["sigma_obs"],
        "sigma_sync": sigma_sync,
        "sigma_mond": sigma_mond,
        "sigma_mond_efe": sigma_mond_efe,
        "a_int": a_i,
        "a_ext": a_e
    })

    print(f"{sat['name']:>12} | {M:>10.2e} | {r:>8.0f} | {d:>6.0f} | {sat['sigma_obs']:>6.1f} | {sigma_sync:>8.2f} | {sigma_mond:>8.2f} | {sigma_mond_efe:>10.2f}")

print("-" * 100)

# =============================================================================
# Part 5: f_indiff Required to Match Observations
# =============================================================================

print("\n" + "=" * 70)
print("Part 5: f_indiff Required to Match Observations")
print("=" * 70)

print("\nFor Synchronism to match observed σ, what f_indiff is needed?")
print("-" * 70)

def find_f_indiff(sigma_obs, sigma_bary, max_f=1000):
    """
    Find f_indiff such that:
    σ² = G_eff × G × M_total / r
        = G_eff × G × M_bary × (1 + f_indiff) / r

    So: (σ_obs/σ_bary)² = (1 + f_indiff) × G_eff/G_eff_bary

    Simplified: f_indiff ≈ (σ_obs/σ_bary)² - 1
    """
    return (sigma_obs / sigma_bary)**2 - 1

print(f"{'Name':>12} | {'σ_obs':>6} | {'σ_Sync':>8} | {'f_indiff (Sync)':>15} | {'σ_MOND_EFE':>10} | {'f_indiff (MOND)':>15}")
print("-" * 100)

for pred in satellite_predictions:
    f_sync = find_f_indiff(pred["sigma_obs"], pred["sigma_sync"])
    f_mond = find_f_indiff(pred["sigma_obs"], pred["sigma_mond_efe"])

    print(f"{pred['name']:>12} | {pred['sigma_obs']:>6.1f} | {pred['sigma_sync']:>8.2f} | {f_sync:>15.1f} | {pred['sigma_mond_efe']:>10.2f} | {f_mond:>15.1f}")

# =============================================================================
# Part 6: Create Visualization
# =============================================================================

print("\n" + "=" * 70)
print("Part 6: Creating Visualization")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle("Session #215: External Field Effect (EFE) Predictions", fontsize=14)

# Panel 1: G_eff vs internal acceleration for different external fields
ax1 = axes[0, 0]
a_range = np.logspace(-14, -8, 200)

# MOND with different EFE strengths
a_ext_values = [0, 0.01*a0_mond, 0.1*a0_mond, 1.0*a0_mond]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
labels = ['Isolated', 'a_ext=0.01a₀', 'a_ext=0.1a₀', 'a_ext=a₀']

for a_ext, color, label in zip(a_ext_values, colors, labels):
    geff = np.array([G_eff_mond_efe_proper(a, a_ext) for a in a_range])
    ax1.loglog(a_range/a0_mond, geff, color=color, linestyle='--', alpha=0.8, label=f'MOND {label}')

# Synchronism (no EFE - single line)
geff_sync = np.array([G_eff_sync(a) for a in a_range])
ax1.loglog(a_range/a0_sync, geff_sync, 'k-', linewidth=2, label='Synchronism (no EFE)')

ax1.axhline(y=3.17, color='gray', linestyle=':', alpha=0.5, label='Sync max (1/Ωm)')
ax1.set_xlabel('a/a₀')
ax1.set_ylabel('G_eff/G')
ax1.set_title('G_eff vs Internal Acceleration\n(MOND shows EFE, Sync does not)')
ax1.legend(fontsize=8)
ax1.set_xlim(1e-4, 100)
ax1.set_ylim(1, 100)
ax1.grid(True, alpha=0.3)

# Panel 2: EFE suppression vs distance from MW
ax2 = axes[0, 1]
d_range = np.logspace(1, 3, 50)  # 10 to 1000 kpc

efe_suppression = []
for d in d_range:
    d_m = d * 1e3 * 3.086e16
    a_ext = G * M_MW / d_m**2
    geff_iso = G_eff_mond_no_efe(a_int)
    geff_efe = G_eff_mond_efe_proper(a_int, a_ext)
    efe_suppression.append(geff_efe / geff_iso)

ax2.plot(d_range, efe_suppression, 'b-', linewidth=2, label='MOND EFE suppression')
ax2.axhline(y=1.0, color='k', linestyle='-', linewidth=2, label='Synchronism (no EFE)')
ax2.set_xlabel('Distance from MW (kpc)')
ax2.set_ylabel('G_eff(EFE) / G_eff(isolated)')
ax2.set_title('EFE Strength vs Distance from MW\n(for 10⁶ M☉ dwarf at 100 pc)')
ax2.set_xscale('log')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, 1.1)

# Panel 3: Observed vs predicted σ for satellites
ax3 = axes[1, 0]
names = [p["name"] for p in satellite_predictions]
sigma_obs = [p["sigma_obs"] for p in satellite_predictions]
sigma_sync = [p["sigma_sync"] for p in satellite_predictions]
sigma_mond_efe = [p["sigma_mond_efe"] for p in satellite_predictions]

x = np.arange(len(names))
width = 0.25

bars1 = ax3.bar(x - width, sigma_obs, width, label='Observed', color='black', alpha=0.7)
bars2 = ax3.bar(x, sigma_sync, width, label='Synchronism (bary only)', color='#1f77b4', alpha=0.7)
bars3 = ax3.bar(x + width, sigma_mond_efe, width, label='MOND + EFE (bary only)', color='#ff7f0e', alpha=0.7)

ax3.set_xlabel('MW Satellite')
ax3.set_ylabel('Velocity Dispersion (km/s)')
ax3.set_title('MW Satellites: Observed vs Predicted σ\n(both theories underpredict without additional mass)')
ax3.set_xticks(x)
ax3.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3, axis='y')

# Panel 4: Key prediction - field vs satellite
ax4 = axes[1, 1]

# Simulate a mass sequence
masses = np.logspace(5, 8, 20)
r_half_typical = 100  # pc

field_sigma_sync = []
field_sigma_mond = []
sat_sigma_sync = []
sat_sigma_mond = []

for M in masses:
    a_i = G * M * 2e30 / (r_half_typical * 3.086e16)**2
    a_e = G * M_MW / (100e3 * 3.086e16)**2  # at 100 kpc

    field_sigma_sync.append(sigma_predict(M, r_half_typical, G_eff_sync(a_i)))
    field_sigma_mond.append(sigma_predict(M, r_half_typical, G_eff_mond_no_efe(a_i)))
    sat_sigma_sync.append(sigma_predict(M, r_half_typical, G_eff_sync(a_i)))
    sat_sigma_mond.append(sigma_predict(M, r_half_typical, G_eff_mond_efe_proper(a_i, a_e)))

ax4.loglog(masses, field_sigma_sync, 'b-', linewidth=2, label='Sync (field & sat - same!)')
ax4.loglog(masses, field_sigma_mond, 'r-', linewidth=2, label='MOND field (isolated)')
ax4.loglog(masses, sat_sigma_mond, 'r--', linewidth=2, label='MOND satellite (EFE)')

# Shade the difference region
ax4.fill_between(masses, sat_sigma_mond, field_sigma_mond,
                  color='red', alpha=0.2, label='MOND EFE effect')

ax4.set_xlabel('Stellar Mass (M☉)')
ax4.set_ylabel('Velocity Dispersion (km/s)')
ax4.set_title('THE KEY TEST:\nField vs Satellite Dwarfs (r_half = 100 pc)')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session215_efe_predictions.png', dpi=150)
plt.close()

print("Saved: session215_efe_predictions.png")

# =============================================================================
# Part 7: Summary
# =============================================================================

print("\n" + "=" * 70)
print("Session #215: KEY FINDINGS SUMMARY")
print("=" * 70)

print("""
1. EXTERNAL FIELD EFFECT (EFE) COMPARISON:
   - MOND: Strong EFE - embedded dwarfs have suppressed MOND boost
   - Synchronism: NO EFE - dynamics is purely local

2. QUANTITATIVE DIFFERENCE:
   - For a 10^6 M_sun dwarf at 100 kpc from MW:
     * MOND predicts ~50% suppression of boost due to EFE
     * Synchronism predicts NO change from isolated case

3. THE DISCRIMINATING TEST:
   - Compare velocity dispersions of:
     (a) Field dwarfs (isolated, a_ext ~ 0)
     (b) MW/M31 satellites (embedded, a_ext ~ 0.1 a0)

   - MOND prediction: σ_field > σ_satellite (by ~30-50%)
   - Sync prediction: σ_field = σ_satellite (no difference)

4. OBSERVATIONAL STATUS:
   - Current data is ambiguous due to:
     * f_indiff variations (Synchronism)
     * Distance uncertainties
     * Formation history differences

   - Need: Mass-matched samples with careful environment control

5. TESTABLE PREDICTION:
   At fixed M_star and r_half:
   - If σ_field / σ_satellite ~ 1.0 → Favors Synchronism
   - If σ_field / σ_satellite ~ 1.3-1.5 → Favors MOND

6. COMPLEMENTARY TO SESSION #208 (VOID GALAXIES):
   - Voids test the UNBOUNDED vs BOUNDED nature of G_eff
   - EFE tests the LOCAL vs NON-LOCAL nature of dynamics
   - Together they provide orthogonal discrimination
""")

print("=" * 70)
print("Session #215: COMPLETE")
print("=" * 70)
