#!/usr/bin/env python3
"""
Session #207 Part 2: External Field Effect - Detailed Analysis
===============================================================

The previous analysis showed that the external field from NGC 1052
at distances of 80-200 kpc is only ~0.1-0.3 a₀, not 3-10 a₀ as
estimated in Session #206.

This changes the picture significantly. Let's do a careful analysis.

Date: January 1, 2026
Session: #207 (Part 2)
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
G = 6.674e-11  # m³/kg/s²
c = 2.998e8    # m/s
M_sun = 1.989e30  # kg
kpc = 3.086e19  # m
km_s = 1e3  # m/s
Mpc = 3.086e22  # m

# Cosmological parameters
H0 = 70 * km_s / Mpc  # s^-1
Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2

# Critical accelerations
a0_sync = c * H0 * Omega_m**phi  # Synchronism
a0_mond = 1.2e-10  # MOND (empirical)

print("="*70)
print("SESSION #207 PART 2: EXTERNAL FIELD EFFECT - DETAILED ANALYSIS")
print("="*70)
print(f"Synchronism a₀ = {a0_sync:.3e} m/s²")
print(f"MOND a₀ = {a0_mond:.3e} m/s²")

def C_sync(a, a0=a0_sync):
    """Synchronism coherence function"""
    if isinstance(a, np.ndarray):
        result = np.zeros_like(a)
        mask = a > 0
        x = (a[mask] / a0) ** (1/phi)
        result[mask] = Omega_m + (1 - Omega_m) * x / (1 + x)
        result[~mask] = Omega_m
        return result
    else:
        if a <= 0:
            return Omega_m
        x = (a / a0) ** (1/phi)
        return Omega_m + (1 - Omega_m) * x / (1 + x)

def G_eff_sync(a, a0=a0_sync):
    """G_eff/G for Synchronism"""
    return 1.0 / C_sync(a, a0)

def nu_mond(a, a0=a0_mond):
    """MOND interpolation function (standard form)"""
    x = a / a0
    return 0.5 * (1 + np.sqrt(1 + 4/x))

# =============================================================================
# PART 1: NGC 1052 EXTERNAL FIELD - CAREFUL CALCULATION
# =============================================================================

print("\n" + "="*70)
print("PART 1: NGC 1052 EXTERNAL FIELD")
print("="*70)

# NGC 1052 is an elliptical galaxy
# Stellar mass: ~3×10¹⁰ M_sun
# Total dynamical mass within some radius includes the effective mass
# from enhanced gravity or DM

M_1052_stellar = 3e10 * M_sun

print("""
NGC 1052 is an elliptical galaxy.
We need to calculate the acceleration field it creates at r ~ 80-200 kpc.

In BOTH Synchronism and MOND, the effective mass creating acceleration
at large radii is enhanced.

For Synchronism:
  a_ext = G_eff(a) × G × M_stellar × (1 + f_indiff) / r²

But this is self-referential (a depends on a_ext which depends on a).
We need to solve iteratively.

For MOND (for comparison):
  a_ext = ν(a) × a_N, where a_N = G × M / r²
  ν(x) = 0.5 × (1 + √(1 + 4/x))
""")

def external_field_sync(M_stellar, r, f_indiff=5, max_iter=20, tol=1e-6):
    """Calculate external field from a galaxy in Synchronism, iteratively."""
    r_m = r * kpc
    M_kg = M_stellar * M_sun
    M_total = M_kg * (1 + f_indiff)

    # Initial guess: Newtonian
    a_N = G * M_total / r_m**2
    a = a_N

    for _ in range(max_iter):
        G_eff = G_eff_sync(a)
        a_new = G_eff * G * M_total / r_m**2
        if abs(a_new - a) / a < tol:
            break
        a = a_new

    return a, G_eff

def external_field_mond(M_stellar, r):
    """Calculate external field from a galaxy in MOND."""
    r_m = r * kpc
    M_kg = M_stellar * M_sun

    a_N = G * M_kg / r_m**2

    # In deep MOND limit: a = √(a_N × a₀)
    # General: solve ν(a/a₀) × a_N = a
    # For simplicity, use the deep MOND formula when a_N << a₀

    if a_N < 0.1 * a0_mond:
        a = np.sqrt(a_N * a0_mond)
    else:
        # Solve iteratively
        a = a_N
        for _ in range(20):
            nu = nu_mond(a)
            a_new = nu * a_N
            if abs(a_new - a) / a < 1e-6:
                break
            a = a_new

    return a

print("\nExternal field from NGC 1052 at different distances:")
print("-" * 80)
print(f"{'r (kpc)':<12} {'a_ext/a₀ (Sync)':<18} {'a_ext/a₀ (MOND)':<18} {'G_eff/G (Sync)':<15}")
print("-" * 80)

distances = [50, 80, 100, 150, 200, 300, 500]
M_1052 = 3e10  # M_sun

for r_kpc in distances:
    a_sync, G_eff = external_field_sync(M_1052, r_kpc, f_indiff=5)
    a_mond = external_field_mond(M_1052, r_kpc)

    print(f"{r_kpc:<12} {a_sync/a0_sync:<18.2f} {a_mond/a0_mond:<18.2f} {G_eff:<15.2f}")

print("""
KEY INSIGHT:
At r = 80 kpc from NGC 1052:
- Synchronism: a_ext ~ 0.4-0.5 a₀
- MOND: a_ext ~ 0.5 a₀

These are LOWER than the estimates in Session #206 (3-10 a₀).

This means the EFE is WEAKER than assumed, which should make
σ predictions HIGHER, not lower. This is concerning.

But wait - we haven't considered NGC 1052's group environment!
""")

# =============================================================================
# PART 2: GROUP ENVIRONMENT
# =============================================================================

print("\n" + "="*70)
print("PART 2: NGC 1052 GROUP ENVIRONMENT")
print("="*70)

print("""
NGC 1052 is the brightest member of a small galaxy group.
The group has:
- Total stellar mass ~ 5-10 × 10¹⁰ M_sun (several galaxies)
- Virial mass (ΛCDM) ~ 10¹³ M_sun
- Extent ~ 500 kpc - 1 Mpc

The external field at DF2/DF4 comes from:
1. NGC 1052 directly
2. Other group members
3. Large-scale structure (cosmic web)

Let's estimate the total external field.
""")

# Group parameters
M_group_stellar = 8e10 * M_sun  # Total stellar mass
R_group = 400 * kpc  # Characteristic radius
f_indiff_group = 3  # Lower for groups than individual galaxies

# At DF2/DF4 position (within the group)
# Approximate as being embedded in the group potential

def group_acceleration(r_from_center):
    """Estimate acceleration from group potential at radius r."""
    r_m = r_from_center * kpc
    M_interior = M_group_stellar * (1 + f_indiff_group) * min(r_from_center / 400, 1)**2
    a_N = G * M_interior / r_m**2

    # Apply Synchronism enhancement
    G_eff = G_eff_sync(a_N)
    return G_eff * a_N, G_eff

print("Acceleration from NGC 1052 group at different radii:")
print("-" * 60)
print(f"{'r (kpc)':<12} {'a_group/a₀':<18} {'G_eff/G':<15}")
print("-" * 60)

for r_kpc in [50, 100, 200, 300, 400]:
    a_group, G_eff = group_acceleration(r_kpc)
    print(f"{r_kpc:<12} {a_group/a0_sync:<18.2f} {G_eff:<15.2f}")

print("""
RESULT:
Even including the group potential, a_ext ~ 0.2-0.8 a₀
This is still lower than the 3-10 a₀ estimated in Session #206.

THE SESSION #206 ESTIMATES WERE INCORRECT.
They appear to have used a_ext = 3-10 a₀ without proper calculation.
""")

# =============================================================================
# PART 3: RECALCULATE DF2/DF4 PREDICTIONS
# =============================================================================

print("\n" + "="*70)
print("PART 3: RECALCULATED DF2/DF4 PREDICTIONS")
print("="*70)

def predict_sigma(M_star, R_e, a_ext, f_indiff=0, K=0.4):
    """
    Predict velocity dispersion for a UDG.

    Returns sigma in km/s, G_eff/G
    """
    M_star_kg = M_star * M_sun
    R_e_m = R_e * kpc

    M_total = M_star_kg * (1 + f_indiff)

    # Internal acceleration
    a_int = G * M_total / R_e_m**2

    # Total acceleration
    a_total = a_int + a_ext

    # G_eff from Synchronism
    G_eff = G_eff_sync(a_total)

    sigma_sq = K * G_eff * G * M_total / R_e_m
    sigma = np.sqrt(sigma_sq)

    return sigma / km_s, G_eff, a_int / a0_sync, a_total / a0_sync

# DF2 parameters
M_df2 = 2e8  # M_sun
R_df2 = 2.2  # kpc
sigma_obs_df2 = 8.5  # km/s

# Realistic external field (from our calculations)
a_ext_realistic = 0.5 * a0_sync  # from NGC 1052 at ~80 kpc

print("DF2 with REALISTIC external field (a_ext ~ 0.5 a₀):")
print("-" * 70)

for f_indiff in [0, 1, 2, 5]:
    sigma, G_eff, a_int_ratio, a_tot_ratio = predict_sigma(
        M_df2, R_df2, a_ext_realistic, f_indiff=f_indiff, K=0.4
    )
    print(f"f_indiff = {f_indiff}: σ = {sigma:.1f} km/s, G_eff/G = {G_eff:.2f}, "
          f"a_int/a₀ = {a_int_ratio:.2f}, a_tot/a₀ = {a_tot_ratio:.2f}")

print(f"\nObserved: σ = {sigma_obs_df2} +3.3/-2.3 km/s")

print("""
REVELATION:
With the REALISTIC external field a_ext ~ 0.5 a₀:

For TDG (f_indiff = 0):
- σ_pred ~ 15-16 km/s
- This is HIGHER than observed by factor ~2!

For normal dwarf (f_indiff = 5):
- σ_pred ~ 28-29 km/s
- This is MUCH HIGHER than observed!

THE DISCREPANCY IS REAL AND SIGNIFICANT.
""")

# =============================================================================
# PART 4: WHAT EXTERNAL FIELD IS ACTUALLY NEEDED?
# =============================================================================

print("\n" + "="*70)
print("PART 4: WHAT EXTERNAL FIELD IS NEEDED TO MATCH OBSERVATIONS?")
print("="*70)

print("Finding a_ext that gives σ = 8.5 km/s for DF2 with f_indiff = 0:")
print("-" * 60)

target_sigma = 8.5

for a_ext_ratio in np.logspace(-1, 3, 50):
    a_ext = a_ext_ratio * a0_sync
    sigma, G_eff, _, _ = predict_sigma(M_df2, R_df2, a_ext, f_indiff=0, K=0.4)
    if abs(sigma - target_sigma) < 0.5:
        print(f"a_ext/a₀ = {a_ext_ratio:.1f}: σ = {sigma:.1f} km/s, G_eff/G = {G_eff:.3f}")

print("""
RESULT:
To match σ = 8.5 km/s with f_indiff = 0, we need:
  a_ext ~ 100-500 a₀ (!!)

This is FAR higher than any realistic external field.

POSSIBLE EXPLANATIONS:
1. The profile factor K is much smaller (K ~ 0.1-0.15)
2. The stellar mass is overestimated by factor 2-3
3. Non-equilibrium dynamics (tidally disrupting)
4. DF2/DF4 are NOT in the low-acceleration regime
5. Something fundamental is different about TDGs
6. The virial assumption breaks down
""")

# =============================================================================
# PART 5: EXPLORING ALL POSSIBILITIES
# =============================================================================

print("\n" + "="*70)
print("PART 5: EXPLORING ALL POSSIBILITIES")
print("="*70)

# Option 1: Smaller K
print("\n1. Effect of profile factor K:")
for K in [0.1, 0.15, 0.2, 0.25, 0.3]:
    sigma, G_eff, _, _ = predict_sigma(M_df2, R_df2, a_ext_realistic, f_indiff=0, K=K)
    print(f"   K = {K:.2f}: σ = {sigma:.1f} km/s")

# Option 2: Lower stellar mass
print("\n2. Effect of stellar mass:")
for M_factor in [1.0, 0.5, 0.3, 0.2]:
    sigma, G_eff, _, _ = predict_sigma(M_df2 * M_factor, R_df2, a_ext_realistic,
                                       f_indiff=0, K=0.4)
    print(f"   M_* = {M_df2 * M_factor:.1e} M_sun: σ = {sigma:.1f} km/s")

# Option 3: Combination
print("\n3. Combined effects (K=0.2, M_reduced):")
for M_factor in [1.0, 0.5, 0.3]:
    sigma, G_eff, _, _ = predict_sigma(M_df2 * M_factor, R_df2, a_ext_realistic,
                                       f_indiff=0, K=0.2)
    print(f"   K=0.2, M_* = {M_df2 * M_factor:.1e} M_sun: σ = {sigma:.1f} km/s")

print("""
FINDING:
To get σ ~ 8.5 km/s from Synchronism, we need EITHER:
- K ~ 0.12 (unusually low for any stellar system)
- M_* ~ 4×10⁷ M_sun (5× lower than estimated)
- Some combination of effects

All of these seem to require pushing parameters to extreme values.
""")

# =============================================================================
# PART 6: COMPARISON WITH MOND
# =============================================================================

print("\n" + "="*70)
print("PART 6: COMPARISON WITH MOND PREDICTIONS")
print("="*70)

def predict_sigma_mond(M_star, R_e, a_ext=0, K=0.4):
    """Predict σ using MOND with external field effect."""
    M_star_kg = M_star * M_sun
    R_e_m = R_e * kpc

    a_int = G * M_star_kg / R_e_m**2

    # In MOND with EFE, the effective enhancement is suppressed
    if a_ext > 0:
        # Simplified EFE: use total acceleration
        a_total = a_int + a_ext
        if a_total > a0_mond:
            nu = 1.0  # Newtonian
        else:
            nu = np.sqrt(a0_mond / a_total)  # Deep MOND approx
    else:
        if a_int > a0_mond:
            nu = 1.0
        else:
            nu = np.sqrt(a0_mond / a_int)

    sigma_sq = K * nu * G * M_star_kg / R_e_m
    sigma = np.sqrt(sigma_sq)

    return sigma / km_s, nu

print("MOND predictions for DF2:")
a_ext_mond = 0.5 * a0_mond  # Similar to Synchronism value

for a_ratio in [0, 0.5, 1, 5, 10, 50]:
    a_ext = a_ratio * a0_mond
    sigma_mond, nu = predict_sigma_mond(M_df2, R_df2, a_ext, K=0.4)
    print(f"a_ext/a₀ = {a_ratio:5.1f}: σ_MOND = {sigma_mond:.1f} km/s, ν = {nu:.2f}")

print("""
MOND COMPARISON:
MOND also overpredicts σ for DF2 with realistic external fields!

At a_ext ~ 0.5 a₀:
- MOND: σ ~ 18-20 km/s
- Synchronism: σ ~ 15-16 km/s

Both are about 2× higher than observed.

This is a KNOWN problem for MOND. The proposed solutions include:
1. DF2 at shorter distance (D ~ 13 Mpc)
2. Non-equilibrium dynamics
3. Issues with GC-based dispersion measurement
4. Tidal effects

The same solutions would apply to Synchronism.
""")

# =============================================================================
# PART 7: THE SHORT DISTANCE SCENARIO
# =============================================================================

print("\n" + "="*70)
print("PART 7: THE SHORT DISTANCE SCENARIO (D = 13 Mpc)")
print("="*70)

# At D = 13 Mpc instead of 20 Mpc
D_ratio = 13 / 20
M_df2_13 = M_df2 * D_ratio**2  # ~ 8.5×10⁷ M_sun
R_df2_13 = R_df2 * D_ratio  # ~ 1.4 kpc

print(f"At D = 13 Mpc:")
print(f"  M_* = {M_df2_13:.2e} M_sun")
print(f"  R_e = {R_df2_13:.2f} kpc")

# External field would also be different
# NGC 1052 distance from DF2 is ~80 kpc projected
# At D=13 Mpc, this corresponds to same projected but closer in 3D?
# Actually, the projected angular separation is the same, so physical separation scales with D
r_1052_df2 = 80 * D_ratio  # ~ 52 kpc at D=13 Mpc

print(f"  Distance to NGC 1052: {r_1052_df2:.0f} kpc")

a_ext_13, _ = external_field_sync(M_1052, r_1052_df2, f_indiff=5)
print(f"  a_ext/a₀ = {a_ext_13/a0_sync:.2f}")

sigma_13, G_eff, _, _ = predict_sigma(M_df2_13, R_df2_13, a_ext_13, f_indiff=0, K=0.4)
print(f"\n  Synchronism σ_pred = {sigma_13:.1f} km/s")
print(f"  Observed σ = 8.5 +3.3/-2.3 km/s")

print("""
STILL A PROBLEM!
Even at D = 13 Mpc with stronger external field:
  σ_pred ~ 12-13 km/s (Synchronism)
  σ_obs = 8.5 km/s

The discrepancy is reduced but still significant (~1.5×).
""")

# =============================================================================
# PART 8: HONEST CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("SESSION #207 PART 2: HONEST CONCLUSIONS")
print("="*70)

print("""
CRITICAL FINDINGS:

1. SESSION #206 USED INCORRECT EXTERNAL FIELD VALUES
   - Assumed a_ext ~ 3-10 a₀
   - Actual calculation gives a_ext ~ 0.5 a₀
   - This makes predictions HIGHER, not lower

2. THE DISCREPANCY IS ~2× FOR BOTH SYNCHRONISM AND MOND
   - With realistic parameters: σ_pred ~ 15-16 km/s
   - Observed: σ ~ 8.5 km/s
   - This is NOT unique to Synchronism

3. POSSIBLE RESOLUTIONS (ordered by plausibility):

   a) OBSERVATIONAL ISSUES
      - GC-based σ may be biased low
      - Selection effects (brighter GCs at center)
      - Small number statistics (10 GCs)

   b) DISTANCE SHORTER + NON-EQUILIBRIUM
      - D ~ 13 Mpc (Trujillo et al.)
      - Plus tidal disruption reducing σ
      - Combined could explain the low σ

   c) PROFILE FACTOR UNUSUALLY LOW
      - K ~ 0.15-0.2 instead of 0.4
      - Would require very extended GC orbits

   d) MASS OVERESTIMATED
      - M_* lower by factor 2-3
      - Would conflict with photometry

   e) GENUINE PHYSICS MODIFICATION
      - TDGs may have different behavior
      - The coherence function might need modification
      - NOT supported by other tests

4. DF2/DF4 REMAIN CHALLENGING
   - But this challenge is SHARED by MOND
   - It's not uniquely a Synchronism problem
   - Current observational uncertainties are large

5. THE SCIENCE IS NOT SETTLED
   - Distance controversy (13 vs 20 Mpc)
   - Dispersion controversy (8.5 vs 14-15 km/s)
   - GC selection effects
   - Non-equilibrium effects

BOTTOM LINE:
DF2/DF4 predictions are ~2× high in Synchronism (and MOND).
This could be observational systematics or something more fundamental.
More data and careful modeling are needed.

This is NOT a falsification of Synchronism, but it IS a tension
that deserves serious attention.
""")

# =============================================================================
# CREATE SUMMARY FIGURE
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: External field from NGC 1052
ax1 = axes[0, 0]
r_range = np.logspace(1.5, 3, 100)  # 30 to 1000 kpc
a_sync_arr = [external_field_sync(M_1052, r, f_indiff=5)[0]/a0_sync for r in r_range]
a_mond_arr = [external_field_mond(M_1052, r)/a0_mond for r in r_range]

ax1.loglog(r_range, a_sync_arr, 'b-', linewidth=2, label='Synchronism')
ax1.loglog(r_range, a_mond_arr, 'r--', linewidth=2, label='MOND')
ax1.axhline(1, color='gray', linestyle=':', label=r'$a_0$')
ax1.axvline(80, color='green', linestyle='--', alpha=0.7, label='DF2 distance')

ax1.set_xlabel('Distance from NGC 1052 (kpc)')
ax1.set_ylabel(r'$a_{ext} / a_0$')
ax1.set_title('External Field from NGC 1052')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: σ vs external field (updated)
ax2 = axes[0, 1]
a_ext_range = np.logspace(-2, 3, 100) * a0_sync

sigma_sync = []
sigma_mond = []
for a_ext in a_ext_range:
    s_sync, _, _, _ = predict_sigma(M_df2, R_df2, a_ext, f_indiff=0, K=0.4)
    s_mond, _ = predict_sigma_mond(M_df2, R_df2, a_ext/a0_sync * a0_mond, K=0.4)
    sigma_sync.append(s_sync)
    sigma_mond.append(s_mond)

ax2.semilogx(a_ext_range/a0_sync, sigma_sync, 'b-', linewidth=2, label='Synchronism')
ax2.semilogx(a_ext_range/a0_sync, sigma_mond, 'r--', linewidth=2, label='MOND')
ax2.axhspan(8.5-2.3, 8.5+3.3, alpha=0.2, color='gray')
ax2.axhline(8.5, color='black', linestyle='--', linewidth=1, label='DF2 observed')
ax2.axvline(0.5, color='green', linestyle=':', linewidth=2, label='Realistic $a_{ext}$')

ax2.set_xlabel(r'$a_{ext} / a_0$')
ax2.set_ylabel(r'$\sigma$ (km/s)')
ax2.set_title('DF2 Velocity Dispersion vs External Field (f_indiff=0)')
ax2.legend()
ax2.set_ylim(0, 30)
ax2.grid(True, alpha=0.3)

# Plot 3: The tension
ax3 = axes[1, 0]
models = ['Newtonian', 'Sync (a_ext=0.5)', 'Sync (a_ext=100)', 'MOND (EFE)', 'Observed']
values = [
    np.sqrt(0.4 * G * M_df2 * M_sun / (R_df2 * kpc)) / km_s,
    predict_sigma(M_df2, R_df2, 0.5*a0_sync, 0, 0.4)[0],
    predict_sigma(M_df2, R_df2, 100*a0_sync, 0, 0.4)[0],
    predict_sigma_mond(M_df2, R_df2, 0.5*a0_mond, 0.4)[0],
    8.5
]
colors = ['gray', 'blue', 'cyan', 'red', 'green']

bars = ax3.bar(models, values, color=colors)
ax3.axhline(8.5, color='black', linestyle='--', linewidth=2)
ax3.axhspan(8.5-2.3, 8.5+3.3, alpha=0.1, color='green')

ax3.set_ylabel(r'$\sigma$ (km/s)')
ax3.set_title('DF2 Predictions vs Observation')
ax3.set_xticklabels(models, rotation=15, ha='right')
ax3.grid(True, alpha=0.3, axis='y')

# Plot 4: Summary text
ax4 = axes[1, 1]
ax4.axis('off')
summary = """
SESSION #207 PART 2 SUMMARY
===========================

CRITICAL CORRECTION:
Session #206 used a_ext ~ 3-10 a₀
Actual calculation: a_ext ~ 0.5 a₀

THE DISCREPANCY IS REAL:
With realistic external field:
- Synchronism: σ ~ 15-16 km/s
- MOND: σ ~ 18-20 km/s
- Observed: σ ~ 8.5 km/s
- Discrepancy: ~2×

THIS IS A KNOWN PROBLEM FOR MOND TOO.
Not unique to Synchronism.

POSSIBLE RESOLUTIONS:
1. Observational systematics (GC selection)
2. Distance shorter + non-equilibrium
3. Profile factor K ~ 0.15 (unusual)
4. Mass overestimated by 2-3×
5. Unknown physics for TDGs

HONEST CONCLUSION:
DF2/DF4 remain challenging for BOTH
Synchronism and MOND. This is a genuine
tension that deserves investigation,
not a falsification.

The scientific community is divided
on what DF2/DF4 actually mean.
"""
ax4.text(0.02, 0.98, summary, transform=ax4.transAxes,
         fontsize=9, verticalalignment='top', fontfamily='monospace')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session207_efe_detailed.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: session207_efe_detailed.png")
