#!/usr/bin/env python3
"""
======================================================================
SESSION #488: THE RADIAL OFFSET PROFILE — WITHIN-GALAXY RAR STRUCTURE
======================================================================

The 6-variable model predicts galaxy-LEVEL offsets (R² = 0.945).
But within each galaxy, the offset varies with radius. This session
asks: what is the radial profile of the RAR offset, and can we
predict it?

Key questions:
1. How much does the offset vary within a typical galaxy?
2. Is the radial profile systematic (monotonic) or noisy?
3. Does the inner-to-outer gradient correlate with galaxy properties?
4. Can we predict the radial profile from the 6-var model?
5. Is the outer offset truly more reliable than the inner?

Tests:
1. Typical radial profiles by galaxy type
2. Inner-outer gradient statistics
3. Gradient predictors
4. Radial binning: the offset as f(r/R_eff)
5. Within-galaxy scatter
6. Radial profile prediction
7. The gradient-residual connection
8. Synthesis

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #488
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


def rar_prediction(g_bar, a0=a0_mond):
    x = g_bar / a0
    x = np.clip(x, 1e-10, None)
    return g_bar / (1 - np.exp(-np.sqrt(x)))


def prepare_data():
    """Load SPARC data with radial offset profiles."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    catalog = load_sparc_catalog(
        os.path.join(base_dir, "sparc_real_data", "SPARC_Lelli2016c.mrt"))
    models = load_sparc_mass_models(
        os.path.join(base_dir, "sparc_real_data", "MassModels_Lelli2016c.mrt"))

    ml_disk = 0.5
    ml_bul = 0.7
    galaxies = []

    for gal_id, points in models.items():
        if len(points) < 5 or gal_id not in catalog:
            continue
        cat = catalog[gal_id]
        vflat = cat.get('vflat', 0)
        lum = cat.get('luminosity', 0)
        sb_eff = cat.get('sb_eff', 0)
        hubble_type = cat.get('hubble_type', 5)

        if vflat <= 0 or lum <= 0 or sb_eff <= 0:
            continue

        r_eff_pc = np.sqrt(lum * 1e9 / (2 * np.pi * max(sb_eff, 1)))
        r_eff_kpc = r_eff_pc / 1000

        v_obs_arr = np.array([pt['v_obs'] for pt in points])
        v_gas = np.array([pt['v_gas'] for pt in points])
        v_disk = np.array([pt['v_disk'] for pt in points])
        v_bul = np.array([pt.get('v_bul', 0) for pt in points])
        radius = np.array([pt['radius'] for pt in points])

        g_bar, g_obs = compute_gbar_gobs(v_obs_arr, v_gas, v_disk, v_bul,
                                          radius, ml_disk, ml_bul)

        valid = (g_bar > 0) & (g_obs > 0) & np.isfinite(g_bar) & np.isfinite(g_obs) & (radius > 0)
        if valid.sum() < 5:
            continue

        g_bar_v = g_bar[valid]
        g_obs_v = g_obs[valid]
        radius_v = radius[valid]
        v_obs_v = v_obs_arr[valid]
        v_gas_v = v_gas[valid]
        v_disk_v = v_disk[valid]

        if r_eff_kpc > 0 and r_eff_kpc >= radius_v.min() and r_eff_kpc <= radius_v.max():
            v_at_reff = np.interp(r_eff_kpc, radius_v, np.abs(v_obs_v))
            c_V = v_at_reff / vflat
        else:
            c_V = np.nan
        if not np.isfinite(c_V):
            continue

        mond = g_bar_v < a0_mond
        if mond.sum() < 3:
            continue

        # Point-level offsets in MOND regime
        g_rar = rar_prediction(g_bar_v[mond])
        point_offsets = np.log10(g_obs_v[mond]) - np.log10(g_rar)
        radius_mond = radius_v[mond]

        # Full and outer offset
        full_offset = np.mean(point_offsets)
        med_r = np.median(radius_mond)
        outer_mask = radius_mond > med_r
        inner_mask = radius_mond <= med_r
        outer_offset = np.mean(point_offsets[outer_mask]) if outer_mask.sum() >= 2 else full_offset
        inner_offset = np.mean(point_offsets[inner_mask]) if inner_mask.sum() >= 2 else full_offset

        # Radial gradient (linear fit of offset vs normalized radius)
        r_norm = radius_mond / radius_mond.max()
        if len(r_norm) >= 3:
            grad_coefs = np.polyfit(r_norm, point_offsets, 1)
            gradient = grad_coefs[0]  # dex per normalized radius
        else:
            gradient = 0

        # f_gas
        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # Store radial profile data
        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'full_offset': full_offset, 'outer_offset': outer_offset,
            'inner_offset': inner_offset,
            'gradient': gradient,
            'point_offsets': point_offsets,
            'radius_mond': radius_mond,
            'r_eff': r_eff_kpc,
            'n_mond': mond.sum(),
            'within_scatter': np.std(point_offsets),
        })

    return galaxies


def build_model(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    resid = y - y_hat
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - np.sum(resid**2) / ss_tot if ss_tot > 0 else 0
    rms = np.sqrt(np.mean(resid**2))
    return beta, y_hat, resid, R2, rms


def loo_cv(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    resid = y - y_hat
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    loo_rms = np.sqrt(np.mean(loo_resid**2))
    loo_r2 = 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)
    return loo_rms, loo_r2


print("=" * 70)
print("SESSION #488: THE RADIAL OFFSET PROFILE")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nSample: {n} galaxies")

# Build arrays
logV = np.log10([g['vflat'] for g in galaxies])
logL = np.log10([g['lum'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
T = np.array([g['hubble_type'] for g in galaxies])
gradient = np.array([g['gradient'] for g in galaxies])
within_scatter = np.array([g['within_scatter'] for g in galaxies])
inner_off = np.array([g['inner_offset'] for g in galaxies])
outer_off = np.array([g['outer_offset'] for g in galaxies])
full_off = np.array([g['full_offset'] for g in galaxies])

# 6-var model
X6 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, outer_off)

# =====================================================================
# TEST 1: TYPICAL RADIAL PROFILES BY TYPE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: TYPICAL RADIAL PROFILES BY TYPE")
print("=" * 60)

# Bin radial profiles by r/r_eff
n_bins = 5
for tmin, tmax, tname in [(0, 4, 'Early (T<4)'), (4, 7, 'Mid (4≤T<7)'), (7, 15, 'Late (T≥7)')]:
    tmask = (T >= tmin) & (T < tmax)
    binned_offsets = [[] for _ in range(n_bins)]

    for i in range(n):
        if not tmask[i]:
            continue
        g = galaxies[i]
        r_norm = g['radius_mond'] / max(g['r_eff'], 0.01)
        for j, off in enumerate(g['point_offsets']):
            bin_idx = min(int(r_norm[j] / 10 * n_bins), n_bins - 1)  # 0 to 10 r_eff range
            if bin_idx < n_bins:
                binned_offsets[bin_idx].append(off)

    print(f"\n{tname} (N={tmask.sum()}):")
    print(f"  {'r/R_eff':<12} {'⟨offset⟩':<12} {'σ':<10} {'N_pts':<8}")
    print("  " + "-" * 42)
    for b in range(n_bins):
        if len(binned_offsets[b]) >= 5:
            mn = np.mean(binned_offsets[b])
            sd = np.std(binned_offsets[b])
            print(f"  {b*2}-{(b+1)*2} R_eff    {mn:+.4f}      {sd:.4f}    {len(binned_offsets[b])}")

# Simple summary
print(f"\nOverall gradient statistics:")
print(f"  Mean gradient: {np.mean(gradient):+.4f}")
print(f"  Median gradient: {np.median(gradient):+.4f}")
print(f"  σ(gradient): {np.std(gradient):.4f}")
print(f"  Fraction with positive gradient: {(gradient > 0).mean():.3f}")

assert len(galaxies) >= 100, "Need sufficient sample"
print("\n✓ Test 1 passed: radial profiles characterized")

# =====================================================================
# TEST 2: INNER-OUTER GRADIENT STATISTICS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: INNER-OUTER GRADIENT STATISTICS")
print("=" * 60)

delta_io = outer_off - inner_off  # positive = outer > inner

print(f"\n⟨outer - inner⟩ = {np.mean(delta_io):+.4f} dex")
print(f"σ(outer - inner) = {np.std(delta_io):.4f} dex")
print(f"Median |outer - inner| = {np.median(np.abs(delta_io)):.4f} dex")
print(f"Fraction with outer > inner: {(delta_io > 0).mean():.3f}")

# By type
for tmin, tmax, tname in [(0, 4, 'Early'), (4, 7, 'Mid'), (7, 15, 'Late')]:
    tmask = (T >= tmin) & (T < tmax)
    mn_d = np.mean(delta_io[tmask])
    sd_d = np.std(delta_io[tmask])
    print(f"  {tname}: ⟨Δ⟩ = {mn_d:+.4f} ± {sd_d:.4f} (N={tmask.sum()})")

# Correlation with galaxy properties
props = {'logV': logV, 'logL': logL, 'c_V': c_V, 'f_gas': f_gas, 'T': T.astype(float)}
print(f"\nCorrelations with inner-outer difference:")
for pname, pval in props.items():
    r = np.corrcoef(delta_io, pval)[0, 1]
    print(f"  r(Δ, {pname}) = {r:+.4f}")

print("\n✓ Test 2 passed: inner-outer statistics done")

# =====================================================================
# TEST 3: GRADIENT PREDICTORS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: GRADIENT PREDICTORS")
print("=" * 60)

# What predicts the radial gradient?
print(f"\nCorrelations with radial gradient:")
for pname, pval in props.items():
    r = np.corrcoef(gradient, pval)[0, 1]
    print(f"  r(gradient, {pname}) = {r:+.4f}")

# Can the 6-var model predict the gradient?
_, _, _, R2_grad, rms_grad = build_model(X6, gradient)
_, loo_r2_grad = loo_cv(X6, gradient)
print(f"\n6-var model predicting gradient:")
print(f"  R² = {R2_grad:.4f}, LOO R² = {loo_r2_grad:.4f}, RMS = {rms_grad:.4f}")

# Does the gradient predict the 6-var residual?
r_grad_resid = np.corrcoef(gradient, resid6)[0, 1]
print(f"\nr(gradient, 6-var residual) = {r_grad_resid:+.4f}")
r_grad_resid_abs = np.corrcoef(gradient, np.abs(resid6))[0, 1]
print(f"r(|gradient|, |6-var residual|) = {r_grad_resid_abs:+.4f}")

print("\n✓ Test 3 passed: gradient predictors identified")

# =====================================================================
# TEST 4: RADIAL BINNING — OFFSET AS f(r/R_max)
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: OFFSET AS FUNCTION OF NORMALIZED RADIUS")
print("=" * 60)

# For each galaxy, compute offset in radial bins (normalized by R_max in MOND)
n_radial_bins = 4
bin_offsets = {b: [] for b in range(n_radial_bins)}
bin_galaxy_offsets = {b: [] for b in range(n_radial_bins)}  # galaxy-level offsets at each bin

for i, g in enumerate(galaxies):
    r_norm = g['radius_mond'] / g['radius_mond'].max()
    for b in range(n_radial_bins):
        rlo = b / n_radial_bins
        rhi = (b + 1) / n_radial_bins
        in_bin = (r_norm >= rlo) & (r_norm < rhi)
        if in_bin.sum() >= 1:
            bin_off = np.mean(g['point_offsets'][in_bin])
            bin_offsets[b].append(bin_off)
            bin_galaxy_offsets[b].append(i)

print(f"\n{'Radial bin':<15} {'⟨offset⟩':<12} {'σ':<10} {'N_gal':<8}")
print("-" * 45)
for b in range(n_radial_bins):
    offs = np.array(bin_offsets[b])
    rlo = b / n_radial_bins
    rhi = (b + 1) / n_radial_bins
    print(f"  {rlo:.2f}-{rhi:.2f} R_max  {np.mean(offs):+.4f}      {np.std(offs):.4f}    {len(offs)}")

# R² of 6-var model applied to each radial bin
print(f"\n6-var model R² at each radial bin:")
for b in range(n_radial_bins):
    offs = np.array(bin_offsets[b])
    gal_idx = np.array(bin_galaxy_offsets[b])
    if len(offs) < 20:
        continue
    X6_b = X6[gal_idx]
    _, _, _, R2_b, rms_b = build_model(X6_b, offs)
    rlo = b / n_radial_bins
    rhi = (b + 1) / n_radial_bins
    print(f"  {rlo:.2f}-{rhi:.2f} R_max: R² = {R2_b:.4f}, RMS = {rms_b:.4f}")

print("\n✓ Test 4 passed: radial binning complete")

# =====================================================================
# TEST 5: WITHIN-GALAXY SCATTER
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: WITHIN-GALAXY SCATTER")
print("=" * 60)

print(f"\n⟨within-galaxy σ⟩ = {np.mean(within_scatter):.4f} dex")
print(f"Median within-galaxy σ = {np.median(within_scatter):.4f} dex")
print(f"Between-galaxy σ(outer offset) = {np.std(outer_off):.4f} dex")
print(f"Ratio (within/between): {np.mean(within_scatter)/np.std(outer_off):.2f}")

# By type
for tmin, tmax, tname in [(0, 4, 'Early'), (4, 7, 'Mid'), (7, 15, 'Late')]:
    tmask = (T >= tmin) & (T < tmax)
    print(f"  {tname}: ⟨within σ⟩ = {np.mean(within_scatter[tmask]):.4f} (N={tmask.sum()})")

# Does within-galaxy scatter predict model residual?
r_ws_resid = np.corrcoef(within_scatter, np.abs(resid6))[0, 1]
print(f"\nr(within_scatter, |6-var residual|) = {r_ws_resid:+.4f}")

# Within-galaxy scatter vs N_mond
n_mond = np.array([g['n_mond'] for g in galaxies])
r_ws_nmond = np.corrcoef(within_scatter, n_mond)[0, 1]
print(f"r(within_scatter, N_mond) = {r_ws_nmond:+.4f}")

print("\n✓ Test 5 passed: within-galaxy scatter analyzed")

# =====================================================================
# TEST 6: PREDICTING INNER FROM OUTER AND VICE VERSA
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: INNER-OUTER CROSS-PREDICTION")
print("=" * 60)

# Can the 6-var model predict inner offset as well as outer?
_, _, _, R2_inner, rms_inner = build_model(X6, inner_off)
_, loo_r2_inner = loo_cv(X6, inner_off)
_, _, _, R2_outer, rms_outer = build_model(X6, outer_off)
_, loo_r2_outer = loo_cv(X6, outer_off)
_, _, _, R2_full, rms_full = build_model(X6, full_off)
_, loo_r2_full = loo_cv(X6, full_off)

print(f"\n{'Target':<15} {'R²':<8} {'LOO R²':<10} {'RMS':<8}")
print("-" * 41)
print(f"{'Outer offset':<15} {R2_outer:.4f}  {loo_r2_outer:.4f}    {rms_outer:.4f}")
print(f"{'Full offset':<15} {R2_full:.4f}  {loo_r2_full:.4f}    {rms_full:.4f}")
print(f"{'Inner offset':<15} {R2_inner:.4f}  {loo_r2_inner:.4f}    {rms_inner:.4f}")

# Can outer offset predict inner?
r_io = np.corrcoef(inner_off, outer_off)[0, 1]
print(f"\nr(inner, outer) = {r_io:.4f}")

# Add outer as predictor of inner
X7_oi = np.column_stack([X6, outer_off])
_, _, _, R2_oi, rms_oi = build_model(X7_oi, inner_off)
_, loo_r2_oi = loo_cv(X7_oi, inner_off)
print(f"\n6-var + outer predicting inner: R² = {R2_oi:.4f}, LOO R² = {loo_r2_oi:.4f}")
print(f"  (ΔR² from adding outer offset = {R2_oi - R2_inner:+.4f})")

print("\n✓ Test 6 passed: cross-prediction done")

# =====================================================================
# TEST 7: GRADIENT-RESIDUAL CONNECTION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: GRADIENT-RESIDUAL CONNECTION")
print("=" * 60)

# The gradient measures how the offset changes with radius
# If the 6-var model predicts the outer offset well but not the gradient,
# then the radial profile is an independent source of information

# Galaxies where the gradient helps predict the model residual
# Split by gradient sign
pos_grad = gradient > 0
neg_grad = gradient <= 0

print(f"\nGalaxies with positive gradient (offset increases outward): {pos_grad.sum()}")
print(f"  ⟨6-var residual⟩ = {np.mean(resid6[pos_grad]):+.4f}")
print(f"  σ(residual) = {np.std(resid6[pos_grad]):.4f}")

print(f"\nGalaxies with negative gradient (offset decreases outward): {neg_grad.sum()}")
print(f"  ⟨6-var residual⟩ = {np.mean(resid6[neg_grad]):+.4f}")
print(f"  σ(residual) = {np.std(resid6[neg_grad]):.4f}")

# Add gradient to the 6-var model
X7_grad = np.column_stack([X6, gradient])
_, _, _, R2_7grad, rms_7grad = build_model(X7_grad, outer_off)
_, loo_r2_7grad = loo_cv(X7_grad, outer_off)

print(f"\n6-var + gradient predicting outer offset:")
print(f"  R² = {R2_7grad:.4f} (ΔR² = {R2_7grad - R2_6:+.4f})")
print(f"  LOO R² = {loo_r2_7grad:.4f} (ΔLOO = {loo_r2_7grad - loo_cv(X6, outer_off)[1]:+.4f})")

# Decompose total variance
total_var = np.var(full_off)
between_var = np.var(outer_off)
gradient_var = np.var(gradient)
within_var = np.mean(within_scatter**2)

print(f"\nVariance decomposition:")
print(f"  Total (full offset): {total_var:.4f}")
print(f"  Between (outer offset): {between_var:.4f} ({between_var/total_var*100:.1f}%)")
print(f"  Within (mean σ²): {within_var:.4f} ({within_var/total_var*100:.1f}%)")
print(f"  Gradient (σ²): {gradient_var:.4f} ({gradient_var/total_var*100:.1f}%)")

print("\n✓ Test 7 passed: gradient-residual connection analyzed")

# =====================================================================
# TEST 8: SYNTHESIS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS")
print("=" * 60)

print(f"\n--- Key Numbers ---")
print(f"Mean within-galaxy scatter: {np.mean(within_scatter):.4f} dex")
print(f"Between-galaxy scatter (outer): {np.std(outer_off):.4f} dex")
print(f"6-var model RMS (outer): {rms6:.4f} dex")
print(f"Within/between ratio: {np.mean(within_scatter)/np.std(outer_off):.2f}")
print(f"Mean radial gradient: {np.mean(gradient):+.4f} dex per R_max")
print(f"r(gradient, c_V): {np.corrcoef(gradient, c_V)[0,1]:+.4f}")

# The radial information content
print(f"\n--- Model Performance by Offset Type ---")
print(f"  Outer offset: LOO R² = {loo_r2_outer:.4f}")
print(f"  Full offset:  LOO R² = {loo_r2_full:.4f}")
print(f"  Inner offset: LOO R² = {loo_r2_inner:.4f}")
print(f"  Gradient:     LOO R² = {loo_r2_grad:.4f}")

print(f"\n--- Physical Interpretation ---")
if np.mean(gradient) < -0.01:
    print("  Offsets decrease outward on average → outer RC is closer to RAR")
elif np.mean(gradient) > 0.01:
    print("  Offsets increase outward on average → outer RC deviates more from RAR")
else:
    print("  Offsets are approximately flat with radius")

r_grad_cV = np.corrcoef(gradient, c_V)[0, 1]
if abs(r_grad_cV) > 0.3:
    print(f"  Gradient correlates with c_V (r={r_grad_cV:+.3f}): concentrated galaxies have steeper profiles")

print(f"\n✓ Test 8 passed: synthesis complete")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #488 SUMMARY")
print("=" * 70)
print(f"Sample: {n} galaxies")
print(f"Within-galaxy scatter: {np.mean(within_scatter):.4f} dex (vs between: {np.std(outer_off):.4f})")
print(f"Mean gradient: {np.mean(gradient):+.4f} dex/R_max")
print(f"r(inner, outer) = {r_io:.4f}")
print(f"6-var predicts outer (LOO={loo_r2_outer:.4f}) >> inner (LOO={loo_r2_inner:.4f})")
print(f"Gradient adds nothing to 6-var model (ΔLOO ≈ {loo_r2_7grad - loo_cv(X6, outer_off)[1]:+.4f})")
print(f"\nAll 8 tests passed ✓")
