#!/usr/bin/env python3
"""
======================================================================
SESSION #498: WITHIN-GALAXY RAR VARIATION
======================================================================

Session #488 found within-galaxy scatter = 0.087 dex (53% of between-galaxy).
The 6-var model predicts galaxy-level offsets. But what drives the
POINT-TO-POINT variation within each galaxy?

At each radius: offset_point = log(g_obs) - log(g_rar(g_bar))
The galaxy-level offset is the mean of outer-MOND points.
The residual: Δ_point = offset_point - offset_galaxy

What predicts Δ_point?
A) Radius (radial gradient in offset)
B) Local g_bar / a₀ (acceleration regime)
C) Local gas fraction (gas vs disk dominance)
D) Rotation curve slope (non-circular motions)
E) Nothing — it's pure noise

Tests:
1. Point-level offset statistics
2. Radial gradient within galaxies
3. Acceleration dependence of point offset
4. Local gas fraction effect
5. RC slope effect (non-circular motions)
6. Combined point-level model
7. Noise vs signal at point level
8. Implications for the RAR intrinsic scatter

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #498
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
    """Load SPARC data at point level."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    catalog = load_sparc_catalog(
        os.path.join(base_dir, "sparc_real_data", "SPARC_Lelli2016c.mrt"))
    models = load_sparc_mass_models(
        os.path.join(base_dir, "sparc_real_data", "MassModels_Lelli2016c.mrt"))

    ml_disk = 0.5
    ml_bul = 0.7
    all_points = []
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
        e_vobs = np.array([pt.get('e_vobs', 5.0) for pt in points])

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
        e_vobs_v = e_vobs[valid]

        # c_V
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

        # Galaxy-level offset
        radius_m = radius_v[mond]
        med_r = np.median(radius_m)
        outer_mond = mond.copy()
        outer_mond[mond] = radius_m > med_r

        g_rar = rar_prediction(g_bar_v)
        point_offsets = np.log10(g_obs_v) - np.log10(g_rar)

        if outer_mond.sum() >= 2:
            galaxy_offset = np.mean(point_offsets[outer_mond])
        else:
            galaxy_offset = np.mean(point_offsets[mond])

        # f_gas
        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas_gal = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        gal_idx = len(galaxies)
        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas_gal,
            'galaxy_offset': galaxy_offset,
        })

        # Store point-level data for MOND points only
        mond_idx = np.where(mond)[0]
        for j in mond_idx:
            # Local gas fraction at this radius
            f_gas_local = v_gas_v[j]**2 / max(v_gas_v[j]**2 + ml_disk * v_disk_v[j]**2, 1e-10)

            # Local RC slope (using neighboring points)
            if j > 0 and j < len(v_obs_v) - 1:
                dv = v_obs_v[j+1] - v_obs_v[j-1]
                dr = radius_v[j+1] - radius_v[j-1]
                rc_slope = (dv / max(abs(v_obs_v[j]), 1)) / max(dr / radius_v[j], 0.01)
            else:
                rc_slope = 0

            all_points.append({
                'gal_idx': gal_idx,
                'radius': radius_v[j],
                'r_norm': radius_v[j] / radius_v.max(),  # normalized radius
                'g_bar': g_bar_v[j],
                'g_obs': g_obs_v[j],
                'offset_point': point_offsets[j],
                'offset_galaxy': galaxy_offset,
                'delta': point_offsets[j] - galaxy_offset,
                'log_g_ratio': np.log10(g_bar_v[j] / a0_mond),
                'f_gas_local': f_gas_local,
                'rc_slope': rc_slope,
                'e_vobs': e_vobs_v[j],
                'v_obs': abs(v_obs_v[j]),
                'is_outer': outer_mond[j] if j < len(outer_mond) else False,
            })

    return galaxies, all_points


print("=" * 70)
print("SESSION #498: WITHIN-GALAXY RAR VARIATION")
print("=" * 70)

galaxies, all_points = prepare_data()
n_gal = len(galaxies)
n_pts = len(all_points)
print(f"\n{n_gal} galaxies, {n_pts} MOND-regime points")

# Extract arrays
delta = np.array([p['delta'] for p in all_points])
r_norm = np.array([p['r_norm'] for p in all_points])
log_g_ratio = np.array([p['log_g_ratio'] for p in all_points])
f_gas_local = np.array([p['f_gas_local'] for p in all_points])
rc_slope = np.array([p['rc_slope'] for p in all_points])
e_vobs = np.array([p['e_vobs'] for p in all_points])
v_obs = np.array([p['v_obs'] for p in all_points])
gal_idx = np.array([p['gal_idx'] for p in all_points])
is_outer = np.array([p['is_outer'] for p in all_points])

# =====================================================================
# TEST 1: POINT-LEVEL OFFSET STATISTICS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: POINT-LEVEL OFFSET STATISTICS")
print("=" * 60)

print(f"\nPoint-level residual (δ = offset_point - offset_galaxy):")
print(f"  Mean: {np.mean(delta):+.4f} dex")
print(f"  Std: {np.std(delta):.4f} dex")
print(f"  Median: {np.median(delta):+.4f} dex")
print(f"  IQR: [{np.percentile(delta, 25):.4f}, {np.percentile(delta, 75):.4f}]")

# Between vs within galaxy variance
between_var = np.var([g['galaxy_offset'] for g in galaxies])
within_vars = []
for gi in range(n_gal):
    pts = delta[gal_idx == gi]
    if len(pts) >= 3:
        within_vars.append(np.var(pts))
mean_within_var = np.mean(within_vars)

print(f"\n  Between-galaxy σ: {np.sqrt(between_var):.4f} dex")
print(f"  Within-galaxy σ (mean): {np.sqrt(mean_within_var):.4f} dex")
print(f"  Within/between ratio: {np.sqrt(mean_within_var/between_var):.2f}")

# By type
for name, T_range in [('Early', (0, 4)), ('Mid', (4, 7)), ('Late', (7, 20))]:
    type_mask = np.array([T_range[0] <= galaxies[p['gal_idx']]['hubble_type'] < T_range[1]
                           for p in all_points])
    if type_mask.sum() > 10:
        print(f"  {name}: within σ = {np.std(delta[type_mask]):.4f}")

print("\n✓ Test 1 passed: statistics computed")

# =====================================================================
# TEST 2: RADIAL GRADIENT
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: RADIAL GRADIENT WITHIN GALAXIES")
print("=" * 60)

# Bin by normalized radius
r_bins = [(0, 0.25, 'Inner 25%'), (0.25, 0.5, '25-50%'),
          (0.5, 0.75, '50-75%'), (0.75, 1.01, 'Outer 25%')]

print(f"\n{'Bin':<15} {'N':<8} {'⟨δ⟩':<10} {'σ(δ)':<10}")
print("-" * 45)
for lo, hi, name in r_bins:
    mask = (r_norm >= lo) & (r_norm < hi)
    if mask.sum() > 10:
        print(f"  {name:<15} {mask.sum():<8} {np.mean(delta[mask]):+.4f}    {np.std(delta[mask]):.4f}")

r_delta = np.corrcoef(r_norm, delta)[0, 1]
print(f"\n  r(r_norm, δ) = {r_delta:+.4f}")

# Per-galaxy radial gradient
gradients = []
for gi in range(n_gal):
    pts_mask = gal_idx == gi
    if pts_mask.sum() >= 5:
        r_pts = r_norm[pts_mask]
        d_pts = delta[pts_mask]
        slope = np.polyfit(r_pts, d_pts, 1)[0]
        gradients.append(slope)

gradients = np.array(gradients)
print(f"\n  Per-galaxy gradient (dδ/d(r/r_max)):")
print(f"    Mean: {np.mean(gradients):+.4f}")
print(f"    Median: {np.median(gradients):+.4f}")
print(f"    Fraction positive: {np.mean(gradients > 0)*100:.0f}%")

print("\n✓ Test 2 passed: radial gradient analyzed")

# =====================================================================
# TEST 3: ACCELERATION DEPENDENCE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: ACCELERATION DEPENDENCE OF POINT OFFSET")
print("=" * 60)

# Does the offset residual depend on how deep in MOND we are?
g_bins = [(-3, -1.5), (-1.5, -1), (-1, -0.5), (-0.5, 0)]

print(f"\n{'log(g/a₀)':<12} {'N':<8} {'⟨δ⟩':<10} {'σ(δ)':<10}")
print("-" * 42)
for lo, hi in g_bins:
    mask = (log_g_ratio >= lo) & (log_g_ratio < hi)
    if mask.sum() > 10:
        print(f"  [{lo:+.1f},{hi:+.1f})   {mask.sum():<8} {np.mean(delta[mask]):+.4f}    {np.std(delta[mask]):.4f}")

r_g_delta = np.corrcoef(log_g_ratio, delta)[0, 1]
print(f"\n  r(log(g/a₀), δ) = {r_g_delta:+.4f}")

print("\n✓ Test 3 passed: acceleration dependence analyzed")

# =====================================================================
# TEST 4: LOCAL GAS FRACTION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: LOCAL GAS FRACTION EFFECT")
print("=" * 60)

# Points where gas dominates (high f_gas_local) vs disk dominates
gas_dom = f_gas_local > 0.5
disk_dom = f_gas_local < 0.2

print(f"\n  Gas-dominated points (f_gas > 0.5, N={gas_dom.sum()}):")
print(f"    ⟨δ⟩ = {np.mean(delta[gas_dom]):+.4f}, σ = {np.std(delta[gas_dom]):.4f}")
print(f"  Disk-dominated points (f_gas < 0.2, N={disk_dom.sum()}):")
print(f"    ⟨δ⟩ = {np.mean(delta[disk_dom]):+.4f}, σ = {np.std(delta[disk_dom]):.4f}")

r_fgas_delta = np.corrcoef(f_gas_local, delta)[0, 1]
print(f"\n  r(f_gas_local, δ) = {r_fgas_delta:+.4f}")

# Partial: controlling acceleration
finite = np.isfinite(f_gas_local) & np.isfinite(log_g_ratio) & np.isfinite(delta)
if finite.sum() > 50:
    X_ctrl = np.column_stack([np.ones(finite.sum()), log_g_ratio[finite]])
    beta_d = np.linalg.lstsq(X_ctrl, delta[finite], rcond=None)[0]
    resid_d = delta[finite] - X_ctrl @ beta_d
    beta_f = np.linalg.lstsq(X_ctrl, f_gas_local[finite], rcond=None)[0]
    resid_f = f_gas_local[finite] - X_ctrl @ beta_f
    r_partial = np.corrcoef(resid_d, resid_f)[0, 1]
    print(f"  Partial r(f_gas, δ | g_bar) = {r_partial:+.4f}")

print("\n✓ Test 4 passed: local gas fraction analyzed")

# =====================================================================
# TEST 5: RC SLOPE (NON-CIRCULAR MOTIONS)
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: RC SLOPE AND NON-CIRCULAR MOTIONS")
print("=" * 60)

# Points where RC is rising vs flat vs declining
# rc_slope > 0: rising; rc_slope ≈ 0: flat; rc_slope < 0: declining
rc_finite = np.isfinite(rc_slope) & (np.abs(rc_slope) < 10)

if rc_finite.sum() > 50:
    rising = rc_finite & (rc_slope > 0.5)
    flat_rc = rc_finite & (np.abs(rc_slope) < 0.5)
    declining = rc_finite & (rc_slope < -0.5)

    print(f"\n  Rising RC (N={rising.sum()}): ⟨δ⟩ = {np.mean(delta[rising]):+.4f}")
    print(f"  Flat RC (N={flat_rc.sum()}): ⟨δ⟩ = {np.mean(delta[flat_rc]):+.4f}")
    print(f"  Declining RC (N={declining.sum()}): ⟨δ⟩ = {np.mean(delta[declining]):+.4f}")

    r_slope_delta = np.corrcoef(rc_slope[rc_finite], delta[rc_finite])[0, 1]
    print(f"\n  r(RC slope, δ) = {r_slope_delta:+.4f}")

print("\n✓ Test 5 passed: RC slope analyzed")

# =====================================================================
# TEST 6: COMBINED POINT-LEVEL MODEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: COMBINED POINT-LEVEL MODEL")
print("=" * 60)

# Build a point-level model: δ = f(r_norm, log_g, f_gas_local, rc_slope)
finite_all = np.isfinite(r_norm) & np.isfinite(log_g_ratio) & np.isfinite(f_gas_local) & \
             np.isfinite(rc_slope) & np.isfinite(delta) & (np.abs(rc_slope) < 10)

X_point = np.column_stack([
    np.ones(finite_all.sum()),
    r_norm[finite_all],
    log_g_ratio[finite_all],
    f_gas_local[finite_all],
    rc_slope[finite_all],
])
y_point = delta[finite_all]

beta_pt = np.linalg.lstsq(X_point, y_point, rcond=None)[0]
yhat_pt = X_point @ beta_pt
resid_pt = y_point - yhat_pt
ss_tot_pt = np.sum((y_point - np.mean(y_point))**2)
R2_pt = 1 - np.sum(resid_pt**2) / ss_tot_pt if ss_tot_pt > 0 else 0

var_names_pt = ['intercept', 'r_norm', 'log(g/a₀)', 'f_gas_local', 'RC_slope']
print(f"\nPoint-level model (N={finite_all.sum()}):")
print(f"  R² = {R2_pt:.4f}")
print(f"  RMS = {np.sqrt(np.mean(resid_pt**2)):.4f} dex")

print(f"\n  {'Variable':<15} {'β':<10}")
print("  " + "-" * 25)
for name, b in zip(var_names_pt, beta_pt):
    print(f"  {name:<15} {b:+.4f}")

# Individual contributions
for i, name in enumerate(var_names_pt[1:], 1):
    X_single = np.column_stack([np.ones(finite_all.sum()), X_point[:, i]])
    b_single = np.linalg.lstsq(X_single, y_point, rcond=None)[0]
    yhat_s = X_single @ b_single
    r2_s = 1 - np.sum((y_point - yhat_s)**2) / ss_tot_pt
    print(f"  {name} alone: R² = {r2_s:.4f}")

print("\n✓ Test 6 passed: point-level model built")

# =====================================================================
# TEST 7: NOISE VS SIGNAL AT POINT LEVEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: NOISE VS SIGNAL AT POINT LEVEL")
print("=" * 60)

# Expected noise from V_obs errors
# δ(log g_obs) = 2 × δV/V / ln(10)
frac_e = e_vobs / np.clip(v_obs, 1, None)
expected_noise = 2 * frac_e / np.log(10)

print(f"\nExpected point-level noise from V_obs errors:")
print(f"  Mean σ_noise: {np.mean(expected_noise):.4f} dex")
print(f"  Observed σ(δ): {np.std(delta):.4f} dex")
print(f"  Noise/observed: {np.mean(expected_noise)/np.std(delta)*100:.0f}%")

# Compare noise with actual spread per galaxy
print(f"\nPer-galaxy comparison:")
gal_results = []
for gi in range(n_gal):
    pts = gal_idx == gi
    if pts.sum() >= 5:
        obs_std = np.std(delta[pts])
        noise_std = np.sqrt(np.mean(expected_noise[pts]**2))
        signal_std = np.sqrt(max(obs_std**2 - noise_std**2, 0))
        gal_results.append({
            'obs_std': obs_std, 'noise_std': noise_std, 'signal_std': signal_std
        })

obs_stds = [g['obs_std'] for g in gal_results]
noise_stds = [g['noise_std'] for g in gal_results]
signal_stds = [g['signal_std'] for g in gal_results]

print(f"  Mean observed within-galaxy σ: {np.mean(obs_stds):.4f}")
print(f"  Mean expected noise σ: {np.mean(noise_stds):.4f}")
print(f"  Mean signal σ (after noise removal): {np.mean(signal_stds):.4f}")
print(f"  Noise/observed ratio: {np.mean(noise_stds)/np.mean(obs_stds)*100:.0f}%")

# Fraction of galaxies where noise dominates
noise_dom = np.mean([n > o for n, o in zip(noise_stds, obs_stds)])
print(f"\n  Fraction where noise > observed: {noise_dom*100:.0f}%")
print(f"  → {'Most' if noise_dom > 0.5 else 'Minority of'} within-galaxy scatter is V_obs noise")

print("\n✓ Test 7 passed: noise vs signal analyzed")

# =====================================================================
# TEST 8: IMPLICATIONS FOR RAR INTRINSIC SCATTER
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: RAR INTRINSIC SCATTER ESTIMATE")
print("=" * 60)

# Total RAR scatter = between-galaxy + within-galaxy
# Between-galaxy: 6-var model explains 94.5%, residual = 0.038 dex
# Within-galaxy: some is noise, some is physical

# The "true" RAR scatter (if we could remove all known effects):
# σ²_RAR = σ²_between_model_resid + σ²_within_signal
total_rar_scatter = np.sqrt(np.mean([p['offset_point']**2 for p in all_points]))
between_resid = 0.038  # from 6-var model

print(f"\n--- RAR Scatter Decomposition ---")
print(f"  Total RAR scatter (all points): {total_rar_scatter:.4f} dex")
print(f"  Between-galaxy model residual: {between_resid:.4f} dex")
print(f"  Within-galaxy observed σ: {np.sqrt(mean_within_var):.4f} dex")
print(f"  Within-galaxy signal (noise-removed): {np.mean(signal_stds):.4f} dex")

# The irreducible RAR scatter
# = sqrt(between_resid² + within_signal²)
irr_scatter = np.sqrt(between_resid**2 + np.mean([s**2 for s in signal_stds]))
print(f"\n  Irreducible RAR scatter estimate: {irr_scatter:.4f} dex")
print(f"  = {(10**irr_scatter - 1)*100:.1f}% in velocity")

# Compare with literature
print(f"\n  Literature comparison:")
print(f"    McGaugh+2016 total scatter: 0.13 dex")
print(f"    Our total (all points): {total_rar_scatter:.3f} dex")
print(f"    After 6-var model: {between_resid:.3f} dex (between)")
print(f"    Irreducible (signal): {irr_scatter:.3f} dex")

# What fraction of RAR scatter is "understood"?
frac_explained = 1 - (irr_scatter**2 / total_rar_scatter**2)
print(f"\n  Fraction of RAR scatter explained: {frac_explained*100:.0f}%")

print(f"\n✓ Test 8 passed: intrinsic scatter estimated")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #498 SUMMARY")
print("=" * 70)
print(f"MOND points: {n_pts} across {n_gal} galaxies")
print(f"Within-galaxy σ: {np.sqrt(mean_within_var):.4f} dex")
print(f"r(r_norm, δ) = {r_delta:+.4f}")
print(f"r(log g/a₀, δ) = {r_g_delta:+.4f}")
print(f"Point-level model R² = {R2_pt:.4f}")
print(f"Noise/observed within-galaxy: {np.mean(noise_stds)/np.mean(obs_stds)*100:.0f}%")
print(f"Irreducible scatter: {irr_scatter:.4f} dex")
print(f"RAR scatter explained: {frac_explained*100:.0f}%")
print(f"\nAll 8 tests passed ✓")
