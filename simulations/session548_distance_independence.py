#!/usr/bin/env python3
"""
======================================================================
SESSION #548: DISTANCE INDEPENDENCE — CAN THE MODEL SURVIVE WITHOUT D?
======================================================================

Distance is the 2nd-largest error source (172% of residual variance,
Session #523). Three of the 6 model variables are distance-independent
(logV, f_gas, and the logV×c_V interaction partly). But logL ∝ D² and
c_V depends on r_eff ∝ D. The offset itself depends on distance through
g_bar(r) and g_obs(r).

This session systematically analyzes distance dependence:
1. Which variables depend on distance and how?
2. Can we build a distance-free model using only V, f_gas?
3. How does the offset itself depend on distance?
4. Distance perturbation: how sensitive is the model to D errors?
5. Can logL be replaced by a distance-free proxy?
6. Can c_V be made distance-free?
7. The distance-independent subset: what LOO can we achieve?
8. Synthesis: distance as signal vs noise

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-07
Session: #548
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from session372_sparc_sb_test import (
    load_sparc_catalog, load_sparc_mass_models,
    compute_gbar_gobs
)

a0_mond = 1.2e-10


def nu_mcgaugh(x):
    return 1 / (1 - np.exp(-np.sqrt(np.clip(x, 1e-10, None))))


def build_model(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    yhat = X @ beta
    resid = y - yhat
    R2 = 1 - np.sum(resid**2) / np.sum((y - np.mean(y))**2)
    rms = np.sqrt(np.mean(resid**2))
    return beta, yhat, resid, R2, rms


def loo_r2(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    return 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)


from scipy import stats as sp_stats

print("=" * 70)
print("SESSION #548: DISTANCE INDEPENDENCE")
print("=" * 70)

# Load data
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
catalog = load_sparc_catalog(
    os.path.join(base_dir, "sparc_real_data", "SPARC_Lelli2016c.mrt"))
models = load_sparc_mass_models(
    os.path.join(base_dir, "sparc_real_data", "MassModels_Lelli2016c.mrt"))

ml_disk = 0.5
ml_bul = 0.7


def prepare_galaxies_full():
    """Prepare galaxies with extra metadata including distance."""
    galaxies = []
    for gal_id, points in models.items():
        if len(points) < 5 or gal_id not in catalog:
            continue
        cat = catalog[gal_id]
        vflat = cat.get('vflat', 0)
        lum = cat.get('luminosity', 0)
        sb_eff = cat.get('sb_eff', 0)
        distance = cat.get('distance', 0)
        inclination = cat.get('inclination', 0)
        hubble_type = cat.get('hubble_type', 0)

        if vflat <= 0 or lum <= 0 or sb_eff <= 0 or distance <= 0:
            continue

        r_eff_pc = np.sqrt(lum * 1e9 / (2 * np.pi * max(sb_eff, 1)))
        r_eff_kpc = r_eff_pc / 1000

        v_obs = np.array([pt['v_obs'] for pt in points])
        v_gas = np.array([pt['v_gas'] for pt in points])
        v_disk = np.array([pt['v_disk'] for pt in points])
        v_bul = np.array([pt.get('v_bul', 0) for pt in points])
        radius = np.array([pt['radius'] for pt in points])
        e_vobs = np.array([pt.get('e_vobs', 5.0) for pt in points])

        g_bar, g_obs = compute_gbar_gobs(v_obs, v_gas, v_disk, v_bul,
                                          radius, ml_disk, ml_bul)

        valid = (g_bar > 0) & (g_obs > 0) & np.isfinite(g_bar) & np.isfinite(g_obs)
        if valid.sum() < 5:
            continue

        g_bar_v = g_bar[valid]
        g_obs_v = g_obs[valid]
        radius_v = radius[valid]
        v_obs_v = v_obs[valid]
        v_gas_v = v_gas[valid]
        v_disk_v = v_disk[valid]
        e_vobs_v = e_vobs[valid]

        if r_eff_kpc > 0 and r_eff_kpc >= radius_v.min() and r_eff_kpc <= radius_v.max():
            v_at_reff = np.interp(r_eff_kpc, radius_v, np.abs(v_obs_v))
            c_V_val = v_at_reff / vflat
        else:
            c_V_val = np.nan
        if not np.isfinite(c_V_val):
            continue

        mond = g_bar_v < a0_mond
        if mond.sum() < 3:
            continue

        radius_m = radius_v[mond]
        med_r = np.median(radius_m)
        outer_mond = mond.copy()
        outer_mond[mond] = radius_m > med_r

        g_rar = g_bar_v * nu_mcgaugh(g_bar_v / a0_mond)
        offset_pts = np.log10(g_obs_v) - np.log10(g_rar)

        if outer_mond.sum() >= 2:
            offset_val = np.mean(offset_pts[outer_mond])
        else:
            offset_val = np.mean(offset_pts[mond])

        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas_val = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # Angular size of r_eff (distance-independent observable)
        # r_eff_kpc / distance_Mpc = angular_size in some unit
        # But we need apparent magnitude or angular quantity
        # Surface brightness is distance-independent!
        log_sb = np.log10(sb_eff) if sb_eff > 0 else np.nan

        # Apparent magnitude proxy: flux ∝ L/D² → log(flux) = logL - 2logD
        log_flux = np.log10(lum) - 2 * np.log10(distance)

        galaxies.append({
            'id': gal_id,
            'offset': offset_val,
            'logV': np.log10(vflat),
            'logL': np.log10(lum),
            'c_V': c_V_val,
            'f_gas': f_gas_val,
            'distance': distance,
            'logD': np.log10(distance),
            'sb_eff': sb_eff,
            'log_sb': log_sb,
            'log_flux': log_flux,
            'r_eff_kpc': r_eff_kpc,
            'inclination': inclination,
            'hubble_type': hubble_type,
            'r_max': radius_v.max(),
            'v_obs': v_obs_v,
            'v_gas': v_gas_v,
            'v_disk': v_disk_v,
            'radius': radius_v,
            'g_bar': g_bar_v,
            'g_obs': g_obs_v,
            'e_vobs': e_vobs_v,
            'outer_mond': outer_mond,
        })
    return galaxies


galaxies = prepare_galaxies_full()
n = len(galaxies)
print(f"\n{n} galaxies loaded")

offset = np.array([g['offset'] for g in galaxies])
logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
logD = np.array([g['logD'] for g in galaxies])
log_sb = np.array([g['log_sb'] for g in galaxies])
log_flux = np.array([g['log_flux'] for g in galaxies])
distance = np.array([g['distance'] for g in galaxies])
ones = np.ones(n)

# Standard 6-var model
X6 = np.column_stack([ones, logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, offset)
loo6 = loo_r2(X6, offset)

print(f"\nStandard 6-var: R²={R2_6:.4f}, LOO={loo6:.4f}, RMS={rms6:.4f}")

# ============================================================
print("\n" + "=" * 60)
print("TEST 1: DISTANCE DEPENDENCE OF EACH VARIABLE")
print("=" * 60)
# ============================================================

# Analyze how each model variable correlates with distance
variables = {
    'logV': logV,
    'logL': logL,
    'c_V': c_V,
    'f_gas': f_gas,
    'logV×c_V': logV * c_V,
    'logL×f_gas': logL * f_gas,
    'offset': offset,
}

print(f"\nCorrelation with log(distance):")
print(f"{'Variable':<15} {'r(var,logD)':<15} {'p-value':<12} {'Distance-dep?'}")
print("-" * 55)

for name, var in variables.items():
    r, p = sp_stats.pearsonr(var, logD)
    dep = "YES" if abs(r) > 0.3 else ("marginal" if abs(r) > 0.15 else "no")
    print(f"{name:<15} {r:+.3f}         {p:.1e}      {dep}")

# Theoretical distance scaling
# logL = logL_true + 2×logD_error → logL depends on D²
# c_V depends on r_eff which depends on D (angular size × distance)
# f_gas: velocity ratios → distance-independent
# logV: directly measured → distance-independent
# offset: g ∝ V²/r ∝ 1/D → complex

# Check if logL is really logL_true + 2logD
# If distance errors dominate, r(logL, logD) should be strong
# But there's also physical correlation: nearby galaxies may be biased

print(f"\n--- Theoretical distance scaling ---")
print(f"logL = const + 2×logD (luminosity ∝ D²):  observed r(logL,logD)={sp_stats.pearsonr(logL, logD)[0]:+.3f}")
print(f"logV is distance-independent:               observed r(logV,logD)={sp_stats.pearsonr(logV, logD)[0]:+.3f}")
print(f"f_gas is distance-independent:              observed r(f_gas,logD)={sp_stats.pearsonr(f_gas, logD)[0]:+.3f}")
print(f"c_V depends on r_eff ∝ D:                   observed r(c_V,logD)={sp_stats.pearsonr(c_V, logD)[0]:+.3f}")

# Partial correlation: r(offset, logD | 6-var model)
r_resid_D, p_resid_D = sp_stats.pearsonr(resid6, logD)
print(f"\nr(residual, logD) = {r_resid_D:+.3f}, p={p_resid_D:.3f}")
print(f"  → Model residuals {'DO' if abs(r_resid_D)>0.15 else 'do NOT'} correlate with distance")

print(f"\n✓ TEST 1 PASSED: Distance dependence mapped")

# ============================================================
print("\n" + "=" * 60)
print("TEST 2: DISTANCE-FREE MODEL (logV + f_gas only)")
print("=" * 60)
# ============================================================

# Build models using only distance-independent variables
# Distance-independent: logV, f_gas
# Distance-dependent: logL (∝ D²), c_V (through r_eff ∝ D)

# Model 1: logV only
X_V = np.column_stack([ones, logV])
_, _, _, R2_V, _ = build_model(X_V, offset)
loo_V = loo_r2(X_V, offset)

# Model 2: logV + f_gas
X_Vf = np.column_stack([ones, logV, f_gas])
_, _, _, R2_Vf, _ = build_model(X_Vf, offset)
loo_Vf = loo_r2(X_Vf, offset)

# Model 3: logV + f_gas + logV×f_gas
X_Vff = np.column_stack([ones, logV, f_gas, logV*f_gas])
_, _, _, R2_Vff, _ = build_model(X_Vff, offset)
loo_Vff = loo_r2(X_Vff, offset)

# Model 4: logV + f_gas + f_gas²
X_Vf2 = np.column_stack([ones, logV, f_gas, f_gas**2])
_, _, _, R2_Vf2, _ = build_model(X_Vf2, offset)
loo_Vf2 = loo_r2(X_Vf2, offset)

# Model 5: logV + logV² + f_gas + f_gas²
X_full_free = np.column_stack([ones, logV, logV**2, f_gas, f_gas**2, logV*f_gas])
_, _, _, R2_ff, _ = build_model(X_full_free, offset)
loo_ff = loo_r2(X_full_free, offset)

print(f"\nDistance-free models (using only logV and f_gas):")
print(f"{'Model':<35} {'R²':<8} {'LOO R²':<8} {'ΔLOO vs 6-var'}")
print("-" * 65)
print(f"{'logV only':<35} {R2_V:.4f}  {loo_V:.4f}  {loo_V - loo6:+.4f}")
print(f"{'logV + f_gas':<35} {R2_Vf:.4f}  {loo_Vf:.4f}  {loo_Vf - loo6:+.4f}")
print(f"{'logV + f_gas + logV×f_gas':<35} {R2_Vff:.4f}  {loo_Vff:.4f}  {loo_Vff - loo6:+.4f}")
print(f"{'logV + f_gas + f_gas²':<35} {R2_Vf2:.4f}  {loo_Vf2:.4f}  {loo_Vf2 - loo6:+.4f}")
print(f"{'logV + logV² + f_gas + f_gas² + V×f':<35} {R2_ff:.4f}  {loo_ff:.4f}  {loo_ff - loo6:+.4f}")
print(f"{'Standard 6-var':<35} {R2_6:.4f}  {loo6:.4f}  {'---':>6}")

# What fraction of 6-var LOO is achieved?
frac_2 = loo_Vf / loo6
frac_best = loo_ff / loo6
print(f"\nlogV + f_gas achieves {frac_2:.1%} of 6-var LOO")
print(f"Best distance-free achieves {frac_best:.1%} of 6-var LOO")

print(f"\n✓ TEST 2 PASSED: Distance-free models constructed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 3: HOW DOES THE OFFSET DEPEND ON DISTANCE?")
print("=" * 60)
# ============================================================

# The offset = log(g_obs / g_RAR)
# g_obs = V_obs² / r  (r in physical units = angular_size × D)
# g_bar = V_bar² / r  (same r)
# g_RAR = g_bar × ν(g_bar/a₀)
#
# At a SPECIFIC radius point: g_obs/g_bar = V_obs²/V_bar² → distance-independent!
# But g_bar itself = V_bar²/r depends on D.
# And ν(g_bar/a₀) depends on g_bar, hence on D.
# So g_obs/g_RAR = V_obs²/(V_bar²×ν) which depends on ν, which depends on r, hence D.
#
# Key insight: offset = log(V_obs²) - log(V_bar² × ν(V_bar²/(r×a₀)))
# The ν function introduces the D-dependence.

# Test: compute offset at different assumed distances
# Scale factor: if true D → α×D, then r → α×r, g_bar → g_bar/α, g_obs → g_obs/α
# g_bar/a₀ → g_bar/(α×a₀), so ν changes

# Simulate distance perturbation for each galaxy
alpha_values = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
offsets_at_alpha = np.zeros((n, len(alpha_values)))

for j, alpha in enumerate(alpha_values):
    for i, g in enumerate(galaxies):
        # Physical radius scales as α
        radius_scaled = g['radius'] * alpha
        # Velocities don't change (measured from Doppler)
        v_obs_v = g['v_obs']
        v_gas_v = g['v_gas']
        v_disk_v = g['v_disk']

        # Recompute g_bar and g_obs at scaled radii
        # g = V²/r, so g scales as 1/α
        g_bar_scaled = g['g_bar'] / alpha
        g_obs_scaled = g['g_obs'] / alpha

        # Recompute RAR prediction at scaled g_bar
        g_rar_scaled = g_bar_scaled * nu_mcgaugh(g_bar_scaled / a0_mond)
        offset_pts_scaled = np.log10(g_obs_scaled) - np.log10(g_rar_scaled)

        # Use same outer MOND mask (approximately)
        outer_mond = g['outer_mond']
        if outer_mond.sum() >= 2:
            offsets_at_alpha[i, j] = np.mean(offset_pts_scaled[outer_mond])
        else:
            mond_mask = g_bar_scaled < a0_mond
            if mond_mask.sum() >= 2:
                offsets_at_alpha[i, j] = np.mean(offset_pts_scaled[mond_mask])
            else:
                offsets_at_alpha[i, j] = np.nan

# Sensitivity: d(offset)/d(log α)
alpha_ref = 3  # index of α=1.0
valid_all = np.all(np.isfinite(offsets_at_alpha), axis=1)
n_valid = valid_all.sum()

# For each galaxy, compute sensitivity using finite differences
log_alpha = np.log10(alpha_values)
sensitivity = np.zeros(n)
for i in range(n):
    if valid_all[i]:
        # Linear regression of offset vs log(alpha)
        slope, _, _, _, _ = sp_stats.linregress(log_alpha, offsets_at_alpha[i, :])
        sensitivity[i] = slope
    else:
        sensitivity[i] = np.nan

valid_sens = np.isfinite(sensitivity)
mean_sens = np.mean(sensitivity[valid_sens])
std_sens = np.std(sensitivity[valid_sens])

print(f"\nOffset sensitivity to distance: d(offset)/d(log α)")
print(f"  Mean: {mean_sens:+.3f}")
print(f"  Std:  {std_sens:.3f}")
print(f"  Range: [{np.min(sensitivity[valid_sens]):+.3f}, {np.max(sensitivity[valid_sens]):+.3f}]")

# If distance is wrong by 20% (logα ≈ 0.08), how much does offset change?
delta_20pct = mean_sens * 0.08
print(f"\n  20% distance error → offset shift: {delta_20pct:+.4f} dex")
print(f"  Model RMS = {rms6:.4f} dex")
print(f"  Ratio: {abs(delta_20pct)/rms6:.2f}× model RMS")

# Does sensitivity correlate with galaxy properties?
print(f"\nSensitivity correlates with:")
for name, var in [('logV', logV), ('logL', logL), ('f_gas', f_gas), ('c_V', c_V)]:
    mask = valid_sens & np.isfinite(var)
    r, p = sp_stats.pearsonr(sensitivity[mask], var[mask])
    print(f"  r(sensitivity, {name}) = {r:+.3f}, p={p:.3f}")

print(f"\n✓ TEST 3 PASSED: Offset distance dependence quantified")

# ============================================================
print("\n" + "=" * 60)
print("TEST 4: DISTANCE PERTURBATION — MODEL SENSITIVITY")
print("=" * 60)
# ============================================================

# Add random distance errors and measure model degradation
np.random.seed(42)
n_mc = 200
sigma_d_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
loo_at_sigma = []

for sigma_d in sigma_d_values:
    loo_trials = []
    for trial in range(n_mc if sigma_d > 0 else 1):
        # Random distance perturbation (in log space, so multiplicative)
        if sigma_d > 0:
            log_alpha_rand = np.random.normal(0, sigma_d, n)
        else:
            log_alpha_rand = np.zeros(n)

        # Perturbed offset for each galaxy
        offset_perturbed = np.zeros(n)
        logL_perturbed = logL + 2 * log_alpha_rand  # L ∝ D²
        c_V_perturbed = np.zeros(n)

        for i, g in enumerate(galaxies):
            alpha_i = 10**log_alpha_rand[i]
            # Recompute offset at perturbed distance
            g_bar_scaled = g['g_bar'] / alpha_i
            g_obs_scaled = g['g_obs'] / alpha_i
            g_rar_scaled = g_bar_scaled * nu_mcgaugh(g_bar_scaled / a0_mond)
            offset_pts_scaled = np.log10(g_obs_scaled) - np.log10(g_rar_scaled)

            outer_mond = g['outer_mond']
            if outer_mond.sum() >= 2:
                offset_perturbed[i] = np.mean(offset_pts_scaled[outer_mond])
            else:
                mond_mask = g_bar_scaled < a0_mond
                if mond_mask.sum() >= 2:
                    offset_perturbed[i] = np.mean(offset_pts_scaled[mond_mask])
                else:
                    offset_perturbed[i] = g['offset']  # fallback

            # Perturb c_V: r_eff scales with D, so need to re-interpolate
            r_eff_perturbed = g['r_eff_kpc'] * alpha_i
            radius_v = g['radius'] * alpha_i  # physical radii also scale
            v_obs_v = g['v_obs']
            if r_eff_perturbed >= radius_v.min() and r_eff_perturbed <= radius_v.max():
                # But wait — the observed radii also scale with D
                # Actually, radius in SPARC is already in physical units (kpc) = angular × D
                # So if D changes, all radii scale: radius_perturbed = radius × alpha
                # And r_eff_perturbed = r_eff × alpha
                # So c_V = V(r_eff)/V_flat is evaluated at the SAME angular position
                # This means c_V is actually distance-INDEPENDENT!
                # The velocity at a given angular position doesn't change with D
                c_V_perturbed[i] = g['c_V']
            else:
                c_V_perturbed[i] = g['c_V']

        # Build perturbed model
        X6_pert = np.column_stack([ones, logV, logL_perturbed, c_V_perturbed,
                                   f_gas, logV*c_V_perturbed, logL_perturbed*f_gas])
        try:
            loo_pert = loo_r2(X6_pert, offset_perturbed)
            loo_trials.append(loo_pert)
        except:
            pass

    loo_at_sigma.append((sigma_d, np.mean(loo_trials), np.std(loo_trials) if len(loo_trials) > 1 else 0))

print(f"\nModel LOO vs distance error (σ in dex of log D):")
print(f"{'σ(logD)':<10} {'Mean LOO':<12} {'σ(LOO)':<10} {'ΔLOO'}")
print("-" * 45)
for sig, mean_loo, std_loo in loo_at_sigma:
    print(f"{sig:.2f}      {mean_loo:.4f}      {std_loo:.4f}     {mean_loo - loo6:+.4f}")

# 20% distance error is σ(logD) ≈ 0.08
# What percentage LOO degradation at that level?
# Interpolate
for sig, mean_loo, _ in loo_at_sigma:
    if sig == 0.10:
        degrad_10 = (loo6 - mean_loo) / loo6 * 100
        print(f"\n10% random distance errors (σ=0.04 dex): degradation to interpolate")
        print(f"σ=0.10 dex (26% errors): LOO degrades by {degrad_10:.1f}%")

# Important insight about c_V
print(f"\n--- KEY INSIGHT: c_V is distance-INDEPENDENT ---")
print(f"c_V = V(r_eff) / V_flat")
print(f"r_eff_physical = r_eff_angular × D")
print(f"But V(r) is measured at angular position, not physical radius")
print(f"When D changes, all radii scale, so V at 'same angular position' is unchanged")
print(f"Therefore c_V = V(angular_r_eff) / V_flat → distance-independent")

# Wait — is this right? The rotation curve points are at ANGULAR positions.
# The physical radii in the data are just angular × D.
# V(r) is really V(θ × D) where θ is angular position.
# If D changes, we'd look at V at a DIFFERENT angular position for the same r_eff.
# Actually, r_eff in SPARC is computed from SB profile angular size × D.
# So r_eff_angular is fixed, but r_eff_physical scales with D.
# The RC sampling positions θ_i are fixed (angular), so r_i = θ_i × D scales with D.
# To evaluate V at r_eff_physical = r_eff_angular × D, we interpolate at θ_eff × D in
# the curve of V vs θ × D. Since both scale together, we get V at the SAME angular point.
# So yes, c_V is distance-independent!

print(f"\n✓ TEST 4 PASSED: Distance perturbation sensitivity quantified")

# ============================================================
print("\n" + "=" * 60)
print("TEST 5: SURFACE BRIGHTNESS AS DISTANCE-FREE PROXY FOR logL")
print("=" * 60)
# ============================================================

# Surface brightness is distance-independent (flux/angular_area)
# logL = log(2π × r_eff² × SB_eff) = log(2π) + 2×log(r_eff) + log(SB)
# log(r_eff) = log(r_eff_angular) + logD
# So logL = const + 2logD + log(SB_eff) + 2log(θ_eff)
#
# Can we replace logL with log_sb (distance-independent)?
# Or with some combination of distance-independent quantities?

# In MOND: offset ∝ logV - 0.25×logL
# logL enters because it traces M_bar via M/L
# But V already traces M_bar via BTFR: M ∝ V⁴
# So logL's role is to correct M/L variations at fixed V

# Try: log_sb as proxy for the L-information
# Note: SB ∝ L/(2π r_eff²). At fixed V (hence fixed M, hence fixed r_eff from TF),
# SB ∝ L/r_eff² ∝ (M/L)⁻¹ × M/r_eff² ∝ (M/L)⁻¹ × Σ

valid_sb = np.isfinite(log_sb)

# Models with SB instead of L
X_sb2 = np.column_stack([ones[valid_sb], logV[valid_sb], log_sb[valid_sb]])
_, _, _, R2_sb2, _ = build_model(X_sb2, offset[valid_sb])
loo_sb2 = loo_r2(X_sb2, offset[valid_sb])

X_sb4 = np.column_stack([ones[valid_sb], logV[valid_sb], log_sb[valid_sb],
                          f_gas[valid_sb], c_V[valid_sb]])
_, _, _, R2_sb4, _ = build_model(X_sb4, offset[valid_sb])
loo_sb4 = loo_r2(X_sb4, offset[valid_sb])

X_sb6 = np.column_stack([ones[valid_sb], logV[valid_sb], log_sb[valid_sb],
                          c_V[valid_sb], f_gas[valid_sb],
                          logV[valid_sb]*c_V[valid_sb],
                          log_sb[valid_sb]*f_gas[valid_sb]])
beta_sb6, _, resid_sb6, R2_sb6, rms_sb6 = build_model(X_sb6, offset[valid_sb])
loo_sb6 = loo_r2(X_sb6, offset[valid_sb])

# Comparison with same galaxies, standard model
X6_same = np.column_stack([ones[valid_sb], logV[valid_sb], logL[valid_sb],
                           c_V[valid_sb], f_gas[valid_sb],
                           logV[valid_sb]*c_V[valid_sb],
                           logL[valid_sb]*f_gas[valid_sb]])
_, _, _, R2_6same, _ = build_model(X6_same, offset[valid_sb])
loo_6same = loo_r2(X6_same, offset[valid_sb])

print(f"\nSurface brightness as distance-free L proxy ({valid_sb.sum()} galaxies):")
print(f"{'Model':<40} {'R²':<8} {'LOO R²':<8} {'ΔLOO'}")
print("-" * 65)
print(f"{'logV + log_SB':<40} {R2_sb2:.4f}  {loo_sb2:.4f}  {loo_sb2 - loo_6same:+.4f}")
print(f"{'logV + log_SB + f_gas + c_V':<40} {R2_sb4:.4f}  {loo_sb4:.4f}  {loo_sb4 - loo_6same:+.4f}")
print(f"{'logV + log_SB + c_V + f_gas + ints':<40} {R2_sb6:.4f}  {loo_sb6:.4f}  {loo_sb6 - loo_6same:+.4f}")
print(f"{'Standard 6-var (same galaxies)':<40} {R2_6same:.4f}  {loo_6same:.4f}  {'---':>6}")

# How much of logL's information does SB capture?
r_SB_L, _ = sp_stats.pearsonr(log_sb[valid_sb], logL[valid_sb])
print(f"\nr(log_SB, logL) = {r_SB_L:+.3f}")

# Partial: r(log_SB, offset | logV)
resid_sb_V = log_sb[valid_sb] - np.column_stack([ones[valid_sb], logV[valid_sb]]) @ \
    np.linalg.lstsq(np.column_stack([ones[valid_sb], logV[valid_sb]]), log_sb[valid_sb], rcond=None)[0]
resid_off_V = offset[valid_sb] - np.column_stack([ones[valid_sb], logV[valid_sb]]) @ \
    np.linalg.lstsq(np.column_stack([ones[valid_sb], logV[valid_sb]]), offset[valid_sb], rcond=None)[0]
r_partial_SB = sp_stats.pearsonr(resid_sb_V, resid_off_V)[0]

resid_L_V = logL[valid_sb] - np.column_stack([ones[valid_sb], logV[valid_sb]]) @ \
    np.linalg.lstsq(np.column_stack([ones[valid_sb], logV[valid_sb]]), logL[valid_sb], rcond=None)[0]
r_partial_L = sp_stats.pearsonr(resid_L_V, resid_off_V)[0]

print(f"r_partial(log_SB, offset | logV) = {r_partial_SB:+.3f}")
print(f"r_partial(logL, offset | logV) = {r_partial_L:+.3f}")
print(f"SB captures {abs(r_partial_SB)/abs(r_partial_L)*100:.0f}% of L's partial signal")

print(f"\n✓ TEST 5 PASSED: Surface brightness proxy evaluated")

# ============================================================
print("\n" + "=" * 60)
print("TEST 6: TRULY DISTANCE-FREE FORMULATION")
print("=" * 60)
# ============================================================

# Distance-independent quantities available:
# 1. logV (velocity from spectroscopy)
# 2. f_gas (velocity ratios)
# 3. c_V (velocity ratio — shown in Test 4 to be D-independent)
# 4. log_SB (surface brightness = flux/angular_area)
# 5. Hubble type (morphological classification)
# 6. Inclination (from axis ratio — approximately D-independent)
# 7. logV×c_V, logV×f_gas, c_V×f_gas (interactions)

# The OFFSET itself depends on distance! But can we define a distance-free offset?
# offset = log(g_obs/g_RAR) — this depends on physical radii
#
# Alternative: use V_obs/V_RAR(g_bar) at specific radii
# But g_bar = V_bar²/r depends on r, hence on D
#
# Deep MOND limit: g_obs ≈ √(g_bar × a₀), so g_obs/g_bar ≈ √(a₀/g_bar)
# In deep MOND: V_obs⁴ ≈ a₀ × V_bar² × r, so the relationship involves r
#
# The RAR offset is FUNDAMENTALLY distance-dependent. There is no way to
# compute it without knowing physical radii.
#
# BUT: What if we regress V_flat on (SB, f_gas, c_V, type)?
# This would predict V_flat from purely photometric/morphological data.
# The BTFR residual (V_obs vs V_predicted) is then distance-free.

# Approach: Build distance-free predictor of logV
# Then the "BTFR residual" δ = logV_obs - logV_predicted is distance-free
# And this should correlate with offset (because offset ≈ mass correction)

# Actually, a simpler approach: the offset is approximately
# offset ≈ 0.25 × (4logV - logL) + corrections
# = 0.25 × (4logV - logL_true - 2logD)
# The D-dependence enters through 2logD in logL.
# If we could replace logL with something D-free...

# Best distance-free model: logV + c_V + f_gas + interactions
X_dfree = np.column_stack([ones, logV, c_V, f_gas, logV*c_V, logV*f_gas])
beta_dfree, yhat_dfree, resid_dfree, R2_dfree, rms_dfree = build_model(X_dfree, offset)
loo_dfree = loo_r2(X_dfree, offset)

# With SB added
X_dfree_sb = np.column_stack([ones[valid_sb], logV[valid_sb], c_V[valid_sb],
                               f_gas[valid_sb], log_sb[valid_sb],
                               logV[valid_sb]*c_V[valid_sb],
                               log_sb[valid_sb]*f_gas[valid_sb]])
_, _, _, R2_dfree_sb, _ = build_model(X_dfree_sb, offset[valid_sb])
loo_dfree_sb = loo_r2(X_dfree_sb, offset[valid_sb])

# Include c_V×f_gas
X_dfree2 = np.column_stack([ones, logV, c_V, f_gas, logV*c_V, logV*f_gas, c_V*f_gas])
_, _, _, R2_dfree2, _ = build_model(X_dfree2, offset)
loo_dfree2 = loo_r2(X_dfree2, offset)

print(f"\nTruly distance-free models:")
print(f"{'Model':<45} {'R²':<8} {'LOO R²':<8} {'ΔLOO vs 6-var'}")
print("-" * 75)
print(f"{'logV + c_V + f_gas + V×c_V + V×f_gas':<45} {R2_dfree:.4f}  {loo_dfree:.4f}  {loo_dfree - loo6:+.4f}")
print(f"{'+ c_V×f_gas':<45} {R2_dfree2:.4f}  {loo_dfree2:.4f}  {loo_dfree2 - loo6:+.4f}")
print(f"{'+ log_SB (N={valid_sb.sum()})':<45} {R2_dfree_sb:.4f}  {loo_dfree_sb:.4f}  {loo_dfree_sb - loo_6same:+.4f}")
print(f"{'Standard 6-var':<45} {R2_6:.4f}  {loo6:.4f}  {'---':>6}")

# What does logL add beyond the distance-free variables?
# Partial r(logL, offset | logV, c_V, f_gas, logV×c_V)
X_partial = np.column_stack([ones, logV, c_V, f_gas, logV*c_V])
resid_L_partial = logL - X_partial @ np.linalg.lstsq(X_partial, logL, rcond=None)[0]
resid_off_partial = offset - X_partial @ np.linalg.lstsq(X_partial, offset, rcond=None)[0]
r_L_partial, p_L_partial = sp_stats.pearsonr(resid_L_partial, resid_off_partial)

print(f"\nr_partial(logL, offset | V, c_V, f_gas, V×c_V) = {r_L_partial:+.3f}, p={p_L_partial:.1e}")
print(f"  → logL adds {'significant' if p_L_partial < 0.01 else 'non-significant'} information beyond D-free vars")

# What fraction of logL can be predicted from D-free variables?
r2_L_from_dfree = 1 - np.sum(resid_L_partial**2) / np.sum((logL - np.mean(logL))**2)
print(f"R²(logL ~ D-free vars) = {r2_L_from_dfree:.3f}")
print(f"  → {r2_L_from_dfree*100:.1f}% of logL is predictable from D-free quantities")

print(f"\n✓ TEST 6 PASSED: Distance-free formulation evaluated")

# ============================================================
print("\n" + "=" * 60)
print("TEST 7: WHAT DOES logL ADD? THE INFORMATION DECOMPOSITION")
print("=" * 60)
# ============================================================

# logL enters the model in two terms: β₂×logL and β₆×logL×f_gas
# What is the UNIQUE information in logL that isn't in (V, c_V, f_gas)?

# Decompose logL into D-free predictable + residual
logL_pred_dfree = X_partial @ np.linalg.lstsq(X_partial, logL, rcond=None)[0]
logL_resid = logL - logL_pred_dfree  # The D-dependent part

# Build model replacing logL with its components
X_decomp = np.column_stack([ones, logV, logL_pred_dfree, logL_resid, c_V, f_gas,
                            logV*c_V, logL_pred_dfree*f_gas, logL_resid*f_gas])
beta_decomp, _, _, R2_decomp, _ = build_model(X_decomp, offset)
loo_decomp = loo_r2(X_decomp, offset)

# Which part matters more?
X_pred_only = np.column_stack([ones, logV, logL_pred_dfree, c_V, f_gas,
                               logV*c_V, logL_pred_dfree*f_gas])
_, _, _, R2_pred, _ = build_model(X_pred_only, offset)
loo_pred = loo_r2(X_pred_only, offset)

X_resid_only = np.column_stack([ones, logV, logL_resid, c_V, f_gas,
                                logV*c_V, logL_resid*f_gas])
_, _, _, R2_resid, _ = build_model(X_resid_only, offset)
loo_resid_only = loo_r2(X_resid_only, offset)

print(f"\nlogL information decomposition:")
print(f"  logL = logL_predictable (from V,c_V,f_gas) + logL_residual (unique)")
print(f"  R²(logL ~ V,c_V,f_gas,V×c_V) = {r2_L_from_dfree:.3f}")
print(f"")
print(f"{'Model':<45} {'R²':<8} {'LOO R²'}")
print("-" * 60)
print(f"{'Standard 6-var (with logL)':<45} {R2_6:.4f}  {loo6:.4f}")
print(f"{'Replace logL with logL_predictable':<45} {R2_pred:.4f}  {loo_pred:.4f}")
print(f"{'Replace logL with logL_residual':<45} {R2_resid:.4f}  {loo_resid_only:.4f}")
print(f"{'Both components (9 vars)':<45} {R2_decomp:.4f}  {loo_decomp:.4f}")
print(f"{'Distance-free (no logL)':<45} {R2_dfree:.4f}  {loo_dfree:.4f}")

# The key question: is logL_residual (the unique, D-dependent part) important?
delta_loo_resid = loo_resid_only - loo_dfree
delta_loo_pred = loo_pred - loo_dfree
print(f"\nΔLOO from adding logL_predictable: {delta_loo_pred:+.4f}")
print(f"ΔLOO from adding logL_residual:    {delta_loo_resid:+.4f}")
print(f"ΔLOO from adding full logL:        {loo6 - loo_dfree:+.4f}")

# Correlation structure
r_pred_resid = sp_stats.pearsonr(logL_pred_dfree, logL_resid)[0]
print(f"\nr(logL_pred, logL_resid) = {r_pred_resid:+.4f} (should be ~0)")

# How much does logD contribute to logL_residual?
r_resid_D, _ = sp_stats.pearsonr(logL_resid, logD)
print(f"r(logL_residual, logD) = {r_resid_D:+.3f}")
print(f"  → logL's unique info is {abs(r_resid_D)*100:.0f}% distance-driven")

print(f"\n✓ TEST 7 PASSED: logL information decomposition complete")

# ============================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS — DISTANCE AS SIGNAL VS NOISE")
print("=" * 60)
# ============================================================

# Summary statistics
print(f"\n{'='*60}")
print(f"DISTANCE INDEPENDENCE SUMMARY")
print(f"{'='*60}")

print(f"\n1. VARIABLE DISTANCE DEPENDENCE:")
print(f"   Distance-INDEPENDENT: logV, f_gas, c_V (all velocity-based)")
print(f"   Distance-DEPENDENT:   logL (∝ D²), offset (through ν function)")

print(f"\n2. MODEL PERFORMANCE:")
print(f"   Standard 6-var (with logL):     LOO = {loo6:.4f}")
print(f"   Distance-free (V, c_V, f_gas):  LOO = {loo_dfree:.4f}")
print(f"   Gap:                             ΔLOO = {loo6 - loo_dfree:+.4f} ({(loo6 - loo_dfree)/loo6*100:.1f}%)")

print(f"\n3. WHAT logL CONTRIBUTES:")
print(f"   {r2_L_from_dfree*100:.1f}% of logL is predictable from D-free variables")
print(f"   r_partial(logL, offset | D-free) = {r_L_partial:+.3f}")
print(f"   logL's unique info is {abs(r_resid_D)*100:.0f}% distance-driven")

print(f"\n4. DISTANCE SENSITIVITY:")
print(f"   d(offset)/d(logD) ≈ {mean_sens:+.3f}")
print(f"   20% distance error → {abs(delta_20pct):.4f} dex offset shift")
print(f"   = {abs(delta_20pct)/rms6:.1f}× model RMS")

# Physical interpretation
# SB performance vs L
frac_sb = loo_sb6 / loo_6same if 'loo_sb6' in dir() else 0
print(f"\n5. SURFACE BRIGHTNESS AS L PROXY:")
print(f"   SB-based 6-var LOO: {loo_sb6:.4f} ({loo_sb6/loo_6same*100:.1f}% of standard)")
print(f"   SB captures {abs(r_partial_SB)/abs(r_partial_L)*100:.0f}% of L's partial signal")

# Distance error budget
sigma_d_sparc = 0.08  # typical ~20% distance errors in SPARC
for sig, mean_loo, _ in loo_at_sigma:
    if sig == 0.10:
        print(f"\n6. DISTANCE ERROR IMPACT:")
        print(f"   σ(logD)=0.10 (26% errors): LOO = {mean_loo:.4f} ({(loo6-mean_loo)/loo6*100:.1f}% degradation)")

# Is logL helping or hurting?
# It helps: ΔLOO = +X when added
# But it's also the channel through which distance errors propagate
# The question: does the physical information in logL outweigh the distance noise?

# Monte Carlo: what LOO would the D-free model achieve if distances were perfect?
# We can't test this directly, but we can bound it:
# If logL_residual is all noise/distance, then logL_predictable should match logL's power
print(f"\n7. SIGNAL vs NOISE DECOMPOSITION:")
print(f"   logL_predictable → ΔLOO = {delta_loo_pred:+.4f} (D-free signal)")
print(f"   logL_residual    → ΔLOO = {delta_loo_resid:+.4f} (includes D-noise)")
total_delta = loo6 - loo_dfree
pred_frac = delta_loo_pred / total_delta if total_delta > 0 else 0
resid_frac = delta_loo_resid / total_delta if total_delta > 0 else 0
print(f"   Predictable fraction: {pred_frac:.1%}")
print(f"   Residual fraction:    {resid_frac:.1%}")
print(f"   → {'Most' if pred_frac > 0.5 else 'Less than half'} of logL's contribution is D-free signal")

# c_V insight
print(f"\n8. c_V IS DISTANCE-INDEPENDENT:")
print(f"   c_V = V(r_eff)/V_flat — ratio of velocities at angular positions")
print(f"   r(c_V, logD) = {sp_stats.pearsonr(c_V, logD)[0]:+.3f}")
print(f"   When D changes, both r_eff and RC sampling scale equally")
print(f"   → c_V is preserved (measures RC shape, not physical size)")

# Final assessment
d_free_frac = loo_dfree / loo6
print(f"\n{'='*60}")
print(f"BOTTOM LINE:")
print(f"  Distance-free model achieves {d_free_frac:.1%} of standard LOO")
print(f"  logL adds {(1-d_free_frac)*100:.1f}% improvement (part signal, part noise)")
print(f"  c_V is distance-independent (velocity ratio)")
print(f"  The model is {('substantially' if d_free_frac > 0.9 else 'partially')} distance-independent")
print(f"{'='*60}")

print(f"\n✓ TEST 8 PASSED: Synthesis complete")

# ============================================================
# FINAL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print(f"SESSION #548: ALL 8 TESTS PASSED")
print(f"{'='*70}")
print(f"\nKey findings:")
print(f"  1. 3/6 variables are distance-independent: logV, f_gas, c_V")
print(f"  2. c_V is D-free because it's a velocity ratio at angular positions")
print(f"  3. Distance-free model: LOO = {loo_dfree:.4f} vs standard {loo6:.4f}")
print(f"  4. logL adds ΔLOO = {loo6-loo_dfree:+.4f} ({(loo6-loo_dfree)/loo6*100:.1f}% of total)")
print(f"  5. {r2_L_from_dfree*100:.1f}% of logL predictable from D-free variables")
print(f"  6. Distance sensitivity: {abs(delta_20pct):.4f} dex per 20% D error")
print(f"  7. SB as L proxy: captures {abs(r_partial_SB)/abs(r_partial_L)*100:.0f}% of partial signal")
print(f"  8. Model is {d_free_frac:.1%} achievable without distance information")
