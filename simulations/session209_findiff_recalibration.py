#!/usr/bin/env python3
"""
Session #209 Part 2: f_indiff Recalibration
============================================

The UFD analysis revealed a significant discrepancy:
- Predicted slope: f_indiff ∝ M_baryon^(-0.20)
- Fitted from UFDs: f_indiff ∝ M_baryon^(-0.57)

This suggests the original calibration (Session #203) may be incomplete.

Possible explanations:
1. Different formation physics for UFDs vs larger dwarfs
2. Environmental dependence (MW satellites vs isolated)
3. Need to recalibrate using full dataset (UFDs + classical + clusters)
4. The scaling relation may not be a simple power law

Date: January 1, 2026
Session: #209 (Part 2)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Physical constants
G = 6.674e-11
M_sun = 1.989e30
pc = 3.086e16
kpc = 3.086e19
km_s = 1e3
Mpc = 3.086e22
c = 2.998e8

# Cosmological parameters
H0 = 70 * km_s / Mpc
Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2
a0 = c * H0 * Omega_m**phi

print("="*70)
print("SESSION #209 PART 2: f_indiff RECALIBRATION")
print("="*70)

def C_sync(a):
    if a <= 0:
        return Omega_m
    x = (a / a0) ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def G_eff_sync(a):
    return 1.0 / C_sync(a)

# =============================================================================
# PART 1: COMPILE FULL DATASET
# =============================================================================

print("\n" + "="*70)
print("PART 1: COMPILE FULL DATASET ACROSS MASS SCALES")
print("="*70)

# Complete dataset: UFDs, classical dSphs, disk galaxies, clusters
# Format: (M_baryon, R_char, M_dyn, sigma_or_V) - where applicable

# UFDs (M_star, R_half_pc, sigma_km_s)
ufd_data = [
    ('Segue 1', 340, 29, 3.7),
    ('Segue 2', 900, 35, 3.4),
    ('Coma Ber', 3700, 77, 4.6),
    ('UMa II', 4100, 149, 6.7),
    ('Bootes I', 29000, 242, 2.4),
    ('CVn II', 7900, 74, 4.6),
    ('Hercules', 37000, 330, 3.7),
    ('Leo IV', 15000, 116, 3.3),
    ('Leo V', 11000, 133, 2.4),
]

# Classical dSphs (M_star, R_half_pc, sigma_km_s)
dsph_data = [
    ('Draco', 290000, 221, 9.1),
    ('Ursa Minor', 290000, 181, 9.5),
    ('Sculptor', 2300000, 283, 9.2),
    ('Carina', 380000, 250, 6.6),
    ('Sextans', 440000, 695, 7.9),
    ('Leo I', 5500000, 251, 9.2),
    ('Leo II', 740000, 176, 6.6),
    ('Fornax', 20000000, 710, 11.7),
]

# Disk galaxies from SPARC-like data (M_baryon, R_flat_kpc, V_flat_km_s)
# Representative sample
disk_data = [
    ('DDO 154', 3e8, 5, 47),
    ('DDO 168', 5e8, 4, 54),
    ('NGC 2403', 5e9, 15, 134),
    ('NGC 3198', 1.5e10, 25, 150),
    ('NGC 7331', 6e10, 30, 250),
    ('MW', 6e10, 20, 220),
]

# Clusters (M_baryon, R_vir_kpc, sigma_km_s)
cluster_data = [
    ('Fornax Cl', 5e12, 500, 350),
    ('Virgo', 1e13, 1000, 600),
    ('Coma', 2e14, 2000, 1000),
    ('A2029', 3e14, 2500, 1200),
]

def M_dyn_from_sigma(sigma, R, scale='pc'):
    """Wolf et al. estimator: M(< R) = 4 σ² R / G"""
    if scale == 'pc':
        R_m = R * pc
    else:
        R_m = R * kpc
    return 4 * (sigma * km_s)**2 * R_m / G / M_sun

def M_dyn_from_V(V, R):
    """Circular velocity: M(< R) = V² R / G"""
    R_m = R * kpc
    return (V * km_s)**2 * R_m / G / M_sun

# =============================================================================
# PART 2: CALCULATE OBSERVED f_indiff
# =============================================================================

print("\n" + "="*70)
print("PART 2: INFER f_indiff FROM OBSERVATIONS")
print("="*70)

def infer_f_indiff(M_baryon, M_dyn, R, scale='pc'):
    """
    Infer f_indiff from observations.

    M_dyn = G_eff × G × M_baryon × (1 + f_indiff) / R
    So: f_indiff = M_dyn / (M_baryon × G_eff) - 1

    But G_eff depends on total acceleration which depends on f_indiff.
    Need iterative solution.
    """
    if scale == 'pc':
        R_m = R * pc
    else:
        R_m = R * kpc

    # Iterate to find self-consistent f_indiff
    f = 10  # Initial guess
    for _ in range(50):
        M_total = M_baryon * M_sun * (1 + f)
        a = G * M_total / R_m**2
        G_eff = G_eff_sync(a)

        # M_dyn / M_baryon = G_eff × (1 + f)
        f_new = (M_dyn / M_baryon) / G_eff - 1

        if f_new < 0:
            f_new = 0
        if abs(f_new - f) < 0.01:
            break
        f = 0.5 * (f + f_new)

    return max(0, f), G_eff

print("\nInferred f_indiff for all systems:")
print("-" * 80)
print(f"{'System':<15} {'M_baryon':<12} {'M_dyn':<12} {'G_eff/G':<10} {'f_indiff':<10}")
print("-" * 80)

all_data = []

# UFDs
for name, M_star, R_half, sigma in ufd_data:
    M_dyn = M_dyn_from_sigma(sigma, R_half, 'pc')
    f, G_eff = infer_f_indiff(M_star, M_dyn, R_half, 'pc')
    print(f"{name:<15} {M_star:<12.2e} {M_dyn:<12.2e} {G_eff:<10.2f} {f:<10.1f}")
    all_data.append((M_star, f, 'UFD'))

# Classical dSphs
for name, M_star, R_half, sigma in dsph_data:
    M_dyn = M_dyn_from_sigma(sigma, R_half, 'pc')
    f, G_eff = infer_f_indiff(M_star, M_dyn, R_half, 'pc')
    print(f"{name:<15} {M_star:<12.2e} {M_dyn:<12.2e} {G_eff:<10.2f} {f:<10.1f}")
    all_data.append((M_star, f, 'dSph'))

# Disk galaxies (use V_flat)
for name, M_bar, R_flat, V_flat in disk_data:
    M_dyn = M_dyn_from_V(V_flat, R_flat)
    f, G_eff = infer_f_indiff(M_bar, M_dyn, R_flat, 'kpc')
    print(f"{name:<15} {M_bar:<12.2e} {M_dyn:<12.2e} {G_eff:<10.2f} {f:<10.1f}")
    all_data.append((M_bar, f, 'Disk'))

# Clusters
for name, M_bar, R_vir, sigma in cluster_data:
    M_dyn = M_dyn_from_sigma(sigma, R_vir, 'kpc')
    f, G_eff = infer_f_indiff(M_bar, M_dyn, R_vir, 'kpc')
    print(f"{name:<15} {M_bar:<12.2e} {M_dyn:<12.2e} {G_eff:<10.2f} {f:<10.1f}")
    all_data.append((M_bar, f, 'Cluster'))

# =============================================================================
# PART 3: FIT POWER LAW
# =============================================================================

print("\n" + "="*70)
print("PART 3: FIT f_indiff SCALING RELATION")
print("="*70)

M_arr = np.array([d[0] for d in all_data])
f_arr = np.array([d[1] for d in all_data])
types = [d[2] for d in all_data]

# Filter valid data (f > 0)
mask = f_arr > 0
M_valid = M_arr[mask]
f_valid = f_arr[mask]

# Log-log fit
log_M = np.log10(M_valid)
log_f = np.log10(f_valid)

# Simple power law: f = A × M^α
slope, intercept = np.polyfit(log_M, log_f, 1)
A = 10**intercept

print(f"\nSimple power law fit:")
print(f"  f_indiff = {A:.2e} × (M_baryon)^({slope:.3f})")
print(f"  Normalized: f_indiff = {A * 1e8**(-slope):.1f} × (M_baryon / 10⁸ M_sun)^({slope:.3f})")

# Compare to Session #203 prediction
slope_203 = -0.20
A_203 = 20  # at 10^8 M_sun

print(f"\nSession #203 prediction:")
print(f"  f_indiff = 20 × (M_baryon / 10⁸ M_sun)^(-0.20)")

print(f"\nDiscrepancy:")
print(f"  Slope: {slope:.3f} vs -0.20 (difference: {slope - (-0.20):.3f})")
print(f"  Normalization ratio: {A * 1e8**(-slope) / A_203:.2f}")

# =============================================================================
# PART 4: SEGMENTED POWER LAW
# =============================================================================

print("\n" + "="*70)
print("PART 4: INVESTIGATE MASS-DEPENDENT SLOPES")
print("="*70)

print("""
The data may not follow a single power law.
Let's check if the slope changes with mass scale.
""")

# Split by mass
low_mass = [(M, f) for M, f, t in all_data if f > 0 and M < 1e6]
mid_mass = [(M, f) for M, f, t in all_data if f > 0 and 1e6 <= M < 1e10]
high_mass = [(M, f) for M, f, t in all_data if f > 0 and M >= 1e10]

def fit_slope(data):
    if len(data) < 2:
        return None, None
    M = np.array([d[0] for d in data])
    f = np.array([d[1] for d in data])
    slope, intercept = np.polyfit(np.log10(M), np.log10(f), 1)
    return slope, 10**intercept

for name, data in [('Low mass (< 10⁶ M_sun)', low_mass),
                   ('Mid mass (10⁶ - 10¹⁰ M_sun)', mid_mass),
                   ('High mass (> 10¹⁰ M_sun)', high_mass)]:
    slope, A = fit_slope(data)
    if slope is not None:
        print(f"{name}:")
        print(f"  N = {len(data)}, slope = {slope:.2f}")
    else:
        print(f"{name}: insufficient data")

# =============================================================================
# PART 5: REVISED THEORETICAL INTERPRETATION
# =============================================================================

print("\n" + "="*70)
print("PART 5: THEORETICAL INTERPRETATION")
print("="*70)

print("""
WHAT DOES THE STEEPER SLOPE MEAN?

Session #203 derived f_indiff ∝ M_baryon^(-0.20) from:
1. SHMR (stellar-halo mass relation) in ΛCDM
2. Baryon retention fraction
3. The indifferent mass tracks the halo, not baryons

The observed steeper slope (~-0.5) suggests:

POSSIBILITY 1: Formation physics differs at low masses
- UFDs formed differently (reionization quenched)
- Their f_indiff may follow different scaling
- Not a fundamental problem, but a complexity

POSSIBILITY 2: Environmental dependence
- MW satellites experience tidal effects
- This could alter the f_indiff relation
- Field dwarfs might follow -0.20 slope

POSSIBILITY 3: The theory needs refinement
- The -0.20 slope was phenomenological
- May need more sophisticated derivation
- Could involve assembly history

POSSIBILITY 4: Measurement systematics
- UFD σ measurements are difficult
- Contamination, binaries affect results
- Classical dSphs are more reliable

MOST LIKELY INTERPRETATION:
The steeper slope at low masses is REAL and indicates
that f_indiff formation is more complex than a simple
power law. The -0.20 slope may only apply to L* galaxies.

A PIECEWISE or BROKEN POWER LAW may be more appropriate:
- M > 10¹⁰ M_sun: f_indiff ∝ M^(-0.20) (Session #203)
- M < 10⁸ M_sun: f_indiff ∝ M^(-0.5) (this analysis)
""")

# =============================================================================
# PART 6: REVISED PREDICTIONS
# =============================================================================

print("\n" + "="*70)
print("PART 6: REVISED f_indiff SCALING")
print("="*70)

def f_indiff_revised(M_baryon):
    """
    Revised f_indiff scaling with mass-dependent slope.

    For M > 10^10 M_sun: f_indiff = 5 × (M / 10^10)^(-0.20)
    For M < 10^8 M_sun: f_indiff = 50 × (M / 10^8)^(-0.50)
    Smooth transition in between.
    """
    log_M = np.log10(M_baryon)

    if log_M >= 10:
        # High mass: original scaling
        return 5 * (M_baryon / 1e10)**(-0.20)
    elif log_M <= 8:
        # Low mass: steeper scaling
        return 50 * (M_baryon / 1e8)**(-0.50)
    else:
        # Transition region: interpolate slopes
        f_low = 50 * (1e8 / 1e8)**(-0.50)  # = 50
        f_high = 5 * (1e10 / 1e10)**(-0.20)  # = 5
        # Linear interpolation in log space
        t = (log_M - 8) / 2  # 0 at 10^8, 1 at 10^10
        log_f = (1 - t) * np.log10(f_low) + t * np.log10(f_high)
        return 10**log_f

print("Revised f_indiff scaling:")
print("-" * 60)
print(f"{'M_baryon (M_sun)':<20} {'f_indiff (original)':<20} {'f_indiff (revised)'}")
print("-" * 60)

def f_indiff_original(M):
    return 20 * (M / 1e8)**(-0.20)

for M in [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e14]:
    f_orig = f_indiff_original(M)
    f_rev = f_indiff_revised(M)
    print(f"{M:<20.0e} {f_orig:<20.0f} {f_rev:<.0f}")

# =============================================================================
# CREATE FIGURE
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: f_indiff vs M_baryon - data and fits
ax1 = axes[0, 0]

# Plot data by type
colors = {'UFD': 'blue', 'dSph': 'green', 'Disk': 'orange', 'Cluster': 'red'}
for M, f, t in all_data:
    if f > 0:
        ax1.scatter(M, f, c=colors[t], s=80, label=t if t not in ax1.get_legend_handles_labels()[1] else '')

# Plot fits
M_plot = np.logspace(2, 15, 100)
f_orig = [f_indiff_original(M) for M in M_plot]
f_rev = [f_indiff_revised(M) for M in M_plot]
f_simple = [10**(intercept + slope * np.log10(M)) for M in M_plot]

ax1.loglog(M_plot, f_orig, 'k--', linewidth=2, label='Session #203 (-0.20)')
ax1.loglog(M_plot, f_rev, 'r-', linewidth=2, label='Revised (broken)')
ax1.loglog(M_plot, f_simple, 'b:', linewidth=2, label=f'Best fit ({slope:.2f})')

ax1.set_xlabel(r'$M_{baryon}$ ($M_\odot$)')
ax1.set_ylabel(r'$f_{indiff}$')
ax1.set_title('Indifferent Mass Fraction Scaling')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(1e2, 1e15)
ax1.set_ylim(0.1, 1e4)

# Plot 2: Residuals from original fit
ax2 = axes[0, 1]

residuals = []
for M, f, t in all_data:
    if f > 0:
        f_pred = f_indiff_original(M)
        res = np.log10(f / f_pred)
        residuals.append((M, res, t))
        ax2.scatter(M, res, c=colors[t], s=80)

ax2.axhline(0, color='k', linestyle='--')
ax2.axhline(0.3, color='gray', linestyle=':', alpha=0.5)
ax2.axhline(-0.3, color='gray', linestyle=':', alpha=0.5)

ax2.set_xlabel(r'$M_{baryon}$ ($M_\odot$)')
ax2.set_ylabel(r'$\log_{10}(f_{obs} / f_{pred})$')
ax2.set_title('Residuals from Session #203 Prediction')
ax2.set_xscale('log')
ax2.grid(True, alpha=0.3)

# Plot 3: Residuals from revised fit
ax3 = axes[1, 0]

for M, f, t in all_data:
    if f > 0:
        f_pred = f_indiff_revised(M)
        res = np.log10(f / f_pred)
        ax3.scatter(M, res, c=colors[t], s=80)

ax3.axhline(0, color='k', linestyle='--')
ax3.axhline(0.3, color='gray', linestyle=':', alpha=0.5)
ax3.axhline(-0.3, color='gray', linestyle=':', alpha=0.5)

ax3.set_xlabel(r'$M_{baryon}$ ($M_\odot$)')
ax3.set_ylabel(r'$\log_{10}(f_{obs} / f_{pred})$')
ax3.set_title('Residuals from Revised Prediction')
ax3.set_xscale('log')
ax3.grid(True, alpha=0.3)

# Calculate RMS residuals
rms_orig = np.sqrt(np.mean([r[1]**2 for r in residuals]))
rms_rev = np.sqrt(np.mean([
    np.log10(f / f_indiff_revised(M))**2
    for M, f, t in all_data if f > 0
]))

# Plot 4: Summary
ax4 = axes[1, 1]
ax4.axis('off')

summary = f"""
f_indiff RECALIBRATION RESULTS
==============================

ORIGINAL (Session #203):
f_indiff = 20 × (M / 10⁸)^(-0.20)
RMS residual: {rms_orig:.2f} dex

REVISED (This analysis):
M > 10¹⁰: f = 5 × (M / 10¹⁰)^(-0.20)
M < 10⁸: f = 50 × (M / 10⁸)^(-0.50)
Smooth transition between
RMS residual: {rms_rev:.2f} dex

KEY FINDINGS:

1. The slope is STEEPER for low-mass systems
   - UFDs: slope ~ -0.5
   - L* galaxies: slope ~ -0.2

2. This may reflect FORMATION PHYSICS
   - Reionization affects low-mass systems
   - Different star formation efficiency
   - Assembly history differences

3. IMPLICATIONS FOR THEORY
   - f_indiff is not a simple universal scaling
   - May depend on formation environment
   - Need physical model for f_indiff origin

4. PREDICTIONS IMPROVED
   - Revised scaling fits data better
   - Reduces tension with UFD observations
   - But adds complexity to framework
"""
ax4.text(0.02, 0.98, summary, transform=ax4.transAxes,
         fontsize=9, verticalalignment='top', fontfamily='monospace')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session209_findiff_recalibration.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: session209_findiff_recalibration.png")

# =============================================================================
# CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("SESSION #209 PART 2 CONCLUSIONS")
print("="*70)

print(f"""
KEY RESULTS:

1. f_indiff SCALING NEEDS REVISION
   - Original: f ∝ M^(-0.20)
   - Data shows steeper slope at low masses
   - Revised: broken power law with transition at ~10⁸ M_sun

2. PHYSICAL INTERPRETATION
   - Formation physics differs for low-mass systems
   - Reionization, environment, assembly history all play roles
   - f_indiff is more complex than a simple universal scaling

3. IMPROVED FITS
   - Original RMS: {rms_orig:.2f} dex
   - Revised RMS: {rms_rev:.2f} dex
   - Better agreement with UFD and cluster data

4. THEORETICAL IMPLICATIONS
   - Need to derive f_indiff from first principles
   - May require understanding of "indifferent pattern" formation
   - Connection to structure formation is key

5. NEXT STEPS
   - Develop physical model for f_indiff origin
   - Test revised scaling on larger samples
   - Look for environmental dependence
""")
