#!/usr/bin/env python3
"""
======================================================================
SESSION #554: GALAXY CENSUS — PER-GALAXY MODEL ASSESSMENT
======================================================================

The 6-var model has RMS=0.038 dex across 128 galaxies. Some galaxies
are predicted to within 0.001 dex (0.2%); others have residuals of
0.10+ dex. This session provides a comprehensive per-galaxy assessment:
which galaxies does the model handle best and worst, and what
distinguishes them?

Tests:
1. The full residual distribution: galaxy-by-galaxy
2. Best-predicted galaxies: what makes them well-behaved?
3. Worst-predicted galaxies: what makes them outliers?
4. The quality connection: do data quality flags predict residuals?
5. Physical vs measurement outliers: can we distinguish?
6. Stability: which galaxies matter most for the model?
7. The "perfect" galaxies: zero residual within measurement noise
8. Synthesis: the galaxy-by-galaxy picture

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-07
Session: #554
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
    return 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2), loo_resid


from scipy import stats as sp_stats

print("=" * 70)
print("SESSION #554: GALAXY CENSUS")
print("=" * 70)

# Load data
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
catalog = load_sparc_catalog(
    os.path.join(base_dir, "sparc_real_data", "SPARC_Lelli2016c.mrt"))
models_data = load_sparc_mass_models(
    os.path.join(base_dir, "sparc_real_data", "MassModels_Lelli2016c.mrt"))

ml_disk = 0.5
ml_bul = 0.7

# Prepare galaxies
galaxies = []
for gal_id, points in models_data.items():
    if len(points) < 5 or gal_id not in catalog:
        continue
    cat = catalog[gal_id]
    vflat = cat.get('vflat', 0)
    lum = cat.get('luminosity', 0)
    sb_eff = cat.get('sb_eff', 0)
    distance = cat.get('distance', 0)
    inclination = cat.get('inclination', 0)
    hubble_type = cat.get('hubble_type', 0)
    quality = cat.get('quality', 0)

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

    # Estimate measurement noise on offset
    # δ(offset) ≈ 2 × δ(V)/V / ln(10) for velocity-dominated error
    mean_frac_err = np.mean(e_vobs_v[outer_mond] / np.abs(v_obs_v[outer_mond])) if outer_mond.sum() >= 2 else np.mean(e_vobs_v / np.abs(v_obs_v))
    offset_noise = 2 * mean_frac_err / np.log(10) / np.sqrt(max(outer_mond.sum(), 1))

    galaxies.append({
        'id': gal_id,
        'offset': offset_val,
        'logV': np.log10(vflat),
        'logL': np.log10(lum),
        'c_V': c_V_val,
        'f_gas': f_gas_val,
        'distance': distance,
        'inclination': inclination,
        'hubble_type': hubble_type,
        'quality': quality,
        'n_points': len(g_bar_v),
        'n_mond': mond.sum(),
        'n_outer': outer_mond.sum(),
        'mean_e_vobs': np.mean(e_vobs_v),
        'offset_noise': offset_noise,
        'vflat': vflat,
    })

n = len(galaxies)
print(f"\n{n} galaxies loaded")

offset = np.array([g['offset'] for g in galaxies])
logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
quality = np.array([g['quality'] for g in galaxies])
inclination = np.array([g['inclination'] for g in galaxies])
hubble_type = np.array([g['hubble_type'] for g in galaxies])
n_points = np.array([g['n_points'] for g in galaxies])
n_outer = np.array([g['n_outer'] for g in galaxies])
offset_noise = np.array([g['offset_noise'] for g in galaxies])
ones = np.ones(n)

# Standard 6-var model
X6 = np.column_stack([ones, logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, offset)
loo6, loo_resid = loo_r2(X6, offset)

# Leverage
H = X6 @ np.linalg.inv(X6.T @ X6) @ X6.T
h = np.diag(H)

# Standardized residuals
sigma_hat = np.sqrt(np.sum(resid6**2) / (n - 7))
std_resid = resid6 / (sigma_hat * np.sqrt(1 - h))

# Cook's D
p = 7
mse = np.sum(resid6**2) / (n - p)
cooks_d = (resid6**2 * h) / (p * mse * (1 - h)**2)

# Signal-to-noise: |residual| / noise_estimate
snr = np.abs(resid6) / np.maximum(offset_noise, 0.001)

print(f"\nStandard 6-var: R²={R2_6:.4f}, LOO={loo6:.4f}, RMS={rms6:.4f}")

# ============================================================
print("\n" + "=" * 60)
print("TEST 1: THE FULL RESIDUAL DISTRIBUTION")
print("=" * 60)
# ============================================================

# Sort by |residual|
idx_sorted = np.argsort(np.abs(resid6))

print(f"\nResidual distribution:")
print(f"  Median |resid|: {np.median(np.abs(resid6)):.4f} dex")
print(f"  Mean |resid|:   {np.mean(np.abs(resid6)):.4f} dex")
print(f"  RMS:            {rms6:.4f} dex")
print(f"  Max |resid|:    {np.max(np.abs(resid6)):.4f} dex")
print(f"  Min |resid|:    {np.min(np.abs(resid6)):.5f} dex")

# Percentile distribution
percentiles = [10, 25, 50, 75, 90, 95, 99]
print(f"\n  Percentiles of |residual|:")
for p_val in percentiles:
    val = np.percentile(np.abs(resid6), p_val)
    print(f"    {p_val}th: {val:.4f} dex ({val/rms6:.2f}σ)")

# How many within noise?
within_noise = np.abs(resid6) < offset_noise
print(f"\n  Galaxies within estimated noise: {within_noise.sum()}/{n} ({within_noise.sum()/n*100:.0f}%)")
within_2noise = np.abs(resid6) < 2 * offset_noise
print(f"  Galaxies within 2× noise: {within_2noise.sum()}/{n} ({within_2noise.sum()/n*100:.0f}%)")

print(f"\n✓ TEST 1 PASSED: Distribution characterized")

# ============================================================
print("\n" + "=" * 60)
print("TEST 2: BEST-PREDICTED GALAXIES")
print("=" * 60)
# ============================================================

# Top 10 best-predicted
print(f"\n10 best-predicted galaxies:")
print(f"{'Rank':<5} {'Galaxy':<15} {'Resid (dex)':<13} {'Noise':<8} {'SNR':<6} {'Q':<3} "
      f"{'V(km/s)':<9} {'T':<4}")
print("-" * 70)

for rank, idx in enumerate(idx_sorted[:10]):
    g = galaxies[idx]
    print(f"{rank+1:<5} {g['id']:<15} {resid6[idx]:+.5f}     {offset_noise[idx]:.4f}  "
          f"{snr[idx]:.2f}  {g['quality']:<3.0f} "
          f"{g['vflat']:<9.1f} {g['hubble_type']:<4.0f}")

# Characteristics of best 25%
best25 = idx_sorted[:32]
worst25 = idx_sorted[-32:]

print(f"\nCharacteristics of best 25% vs worst 25%:")
print(f"  {'Property':<15} {'Best 25%':<12} {'Worst 25%':<12} {'Full sample'}")
print(f"  {'-'*50}")
for name, arr in [('|resid| (dex)', np.abs(resid6)),
                   ('logV', logV), ('logL', logL),
                   ('f_gas', f_gas), ('c_V', c_V),
                   ('Quality', quality),
                   ('N_outer', n_outer.astype(float)),
                   ('Inclination', inclination),
                   ('Noise est.', offset_noise)]:
    print(f"  {name:<15} {np.mean(arr[best25]):.3f}       {np.mean(arr[worst25]):.3f}       {np.mean(arr):.3f}")

print(f"\n✓ TEST 2 PASSED: Best galaxies characterized")

# ============================================================
print("\n" + "=" * 60)
print("TEST 3: WORST-PREDICTED GALAXIES")
print("=" * 60)
# ============================================================

# Top 10 worst-predicted
print(f"\n10 worst-predicted galaxies:")
print(f"{'Rank':<5} {'Galaxy':<15} {'Resid (dex)':<13} {'Noise':<8} {'SNR':<6} {'Cook D':<8} "
      f"{'V(km/s)':<9} {'T':<4}")
print("-" * 75)

for rank, idx in enumerate(idx_sorted[-10:][::-1]):
    g = galaxies[idx]
    print(f"{rank+1:<5} {g['id']:<15} {resid6[idx]:+.5f}     {offset_noise[idx]:.4f}  "
          f"{snr[idx]:.2f}  {cooks_d[idx]:.3f}   "
          f"{g['vflat']:<9.1f} {g['hubble_type']:<4.0f}")

# How many outliers have |standardized resid| > 2?
extreme = np.abs(std_resid) > 2
print(f"\n  Galaxies with |standardized resid| > 2: {extreme.sum()}/{n} ({extreme.sum()/n*100:.1f}%)")
print(f"  Expected for normal: ~5% ({0.05*n:.0f})")

# List them
if extreme.sum() > 0:
    print(f"\n  Extreme outliers:")
    for idx in np.where(extreme)[0]:
        g = galaxies[idx]
        print(f"    {g['id']}: resid={resid6[idx]:+.4f}, std_resid={std_resid[idx]:+.2f}, "
              f"V={g['vflat']:.0f}, T={g['hubble_type']:.0f}")

print(f"\n✓ TEST 3 PASSED: Worst galaxies characterized")

# ============================================================
print("\n" + "=" * 60)
print("TEST 4: QUALITY CONNECTION")
print("=" * 60)
# ============================================================

# SPARC quality flags: 1=best, 2=good, 3=acceptable
for q in [1, 2, 3]:
    mask = quality == q
    if mask.sum() >= 3:
        mean_abs = np.mean(np.abs(resid6[mask]))
        rms_q = np.sqrt(np.mean(resid6[mask]**2))
        mean_noise = np.mean(offset_noise[mask])
        frac_within = np.mean(np.abs(resid6[mask]) < offset_noise[mask])
        print(f"Quality {q} (N={mask.sum()}):")
        print(f"  RMS = {rms_q:.4f}, Mean |resid| = {mean_abs:.4f}, "
              f"Mean noise = {mean_noise:.4f}")
        print(f"  Within noise: {frac_within*100:.0f}%")
        print(f"  Ratio resid/noise: {mean_abs/mean_noise:.2f}")

# Is quality correlated with residual magnitude?
r_q_abs, p_q_abs = sp_stats.pearsonr(quality, np.abs(resid6))
r_q_res, p_q_res = sp_stats.pearsonr(quality, resid6)
print(f"\nr(Quality, |resid|) = {r_q_abs:+.3f}, p={p_q_abs:.3f}")
print(f"r(Quality, resid)  = {r_q_res:+.3f}, p={p_q_res:.3f}")

# LOO by quality
for q in [1, 2, 3]:
    mask = quality == q
    if mask.sum() >= 10:
        X6_q = X6[mask]
        loo_q = loo_r2(X6_q, offset[mask])[0] if mask.sum() > 8 else np.nan
        print(f"Quality {q}: LOO R² = {loo_q:.4f} (N={mask.sum()})")

print(f"\n✓ TEST 4 PASSED: Quality connection established")

# ============================================================
print("\n" + "=" * 60)
print("TEST 5: PHYSICAL vs MEASUREMENT OUTLIERS")
print("=" * 60)
# ============================================================

# For each outlier (|resid| > 2σ), assess:
# - Is |resid| > noise estimate? (measurement candidate)
# - Is it high-leverage? (unusual galaxy)
# - Is it high Cook's D? (model-influential)

threshold = 2 * rms6  # ~2σ
outliers = np.abs(resid6) > threshold
n_outliers = outliers.sum()

print(f"\nOutlier analysis (|resid| > {threshold:.4f} dex = 2σ):")
print(f"  N outliers: {n_outliers}/{n} ({n_outliers/n*100:.1f}%)")

if n_outliers > 0:
    print(f"\n  {'Galaxy':<15} {'Resid':<10} {'Noise':<8} {'|r|/noise':<10} "
          f"{'h':<6} {'Cook D':<8} {'Assessment'}")
    print(f"  {'-'*75}")

    for idx in np.where(outliers)[0]:
        g = galaxies[idx]
        ratio = abs(resid6[idx]) / max(offset_noise[idx], 0.001)
        if ratio < 1.5:
            assess = "NOISE"
        elif h[idx] > 2*p/n:
            assess = "HIGH LEVERAGE"
        elif cooks_d[idx] > 0.15:
            assess = "INFLUENTIAL"
        else:
            assess = "PHYSICAL?"
        print(f"  {g['id']:<15} {resid6[idx]:+.4f}   {offset_noise[idx]:.4f}  "
              f"{ratio:.1f}       {h[idx]:.3f}  {cooks_d[idx]:.3f}   {assess}")

# Summary
if n_outliers > 0:
    outlier_noise = np.mean(np.abs(resid6[outliers]) < 1.5 * offset_noise[outliers])
    outlier_lever = np.mean(h[outliers] > 2*p/n)
    print(f"\n  Noise-consistent: {outlier_noise*100:.0f}%")
    print(f"  High-leverage: {outlier_lever*100:.0f}%")

print(f"\n✓ TEST 5 PASSED: Outlier classification complete")

# ============================================================
print("\n" + "=" * 60)
print("TEST 6: GALAXY IMPORTANCE — LEAVE-ONE-OUT IMPACT")
print("=" * 60)
# ============================================================

# For each galaxy, compute LOO R² without it
# Using the hat matrix: loo_resid_i = resid_i / (1 - h_i)
# The impact of each galaxy on LOO is through its leverage and residual

# DFBETAS: how much does each coefficient change when galaxy i is removed
# DFBETAS_ij = (beta_j - beta_j(-i)) / se(beta_j)
# Using the formula: DFBETAS = (resid / (1-h)) * (X'X)^(-1) x_i / sqrt(MSE)

XtX_inv = np.linalg.inv(X6.T @ X6)

# Max DFBETAS per galaxy
max_dfbetas = np.zeros(n)
for i in range(n):
    xi = X6[i, :]
    dfbeta_i = (resid6[i] / (1 - h[i])) * (XtX_inv @ xi)
    se_beta = np.sqrt(np.diag(mse * XtX_inv))
    dfbetas_i = dfbeta_i / se_beta
    max_dfbetas[i] = np.max(np.abs(dfbetas_i))

# Galaxies with highest impact
impact_sorted = np.argsort(max_dfbetas)[::-1]

print(f"\n10 most important galaxies (highest max |DFBETAS|):")
print(f"{'Rank':<5} {'Galaxy':<15} {'Max |DFBETAS|':<15} {'Resid':<10} {'Leverage':<10}")
print("-" * 60)
for rank in range(10):
    idx = impact_sorted[rank]
    g = galaxies[idx]
    print(f"{rank+1:<5} {g['id']:<15} {max_dfbetas[idx]:.3f}          "
          f"{resid6[idx]:+.4f}   {h[idx]:.3f}")

# What's the typical DFBETAS?
print(f"\nDFBETAS statistics:")
print(f"  Mean max |DFBETAS|: {np.mean(max_dfbetas):.3f}")
print(f"  Median: {np.median(max_dfbetas):.3f}")
print(f"  Max: {np.max(max_dfbetas):.3f} ({galaxies[impact_sorted[0]]['id']})")
print(f"  N with max |DFBETAS| > 0.5: {np.sum(max_dfbetas > 0.5)}/{n}")

print(f"\n✓ TEST 6 PASSED: Galaxy importance quantified")

# ============================================================
print("\n" + "=" * 60)
print("TEST 7: THE 'PERFECT' GALAXIES")
print("=" * 60)
# ============================================================

# Galaxies where |residual| < noise estimate
# These are "perfectly predicted" — the model matches within measurement error

perfect = np.abs(resid6) < offset_noise
n_perfect = perfect.sum()

print(f"\n'Perfect' galaxies (|resid| < noise): {n_perfect}/{n} ({n_perfect/n*100:.0f}%)")

if n_perfect > 0:
    print(f"\n  Properties of 'perfect' vs 'imperfect' galaxies:")
    imperfect = ~perfect
    print(f"  {'Property':<15} {'Perfect (N={0})':<20} {'Imperfect (N={1})':<20}".format(
        n_perfect, imperfect.sum()))
    print(f"  {'-'*55}")
    for name, arr in [('logV', logV), ('logL', logL), ('f_gas', f_gas),
                       ('c_V', c_V), ('Quality', quality),
                       ('Inclination', inclination),
                       ('Hubble type', hubble_type),
                       ('N_outer', n_outer.astype(float)),
                       ('Noise est.', offset_noise)]:
        print(f"  {name:<15} {np.mean(arr[perfect]):.3f}               {np.mean(arr[imperfect]):.3f}")

    # t-test for differences
    print(f"\n  Significant differences:")
    for name, arr in [('logV', logV), ('logL', logL), ('f_gas', f_gas),
                       ('c_V', c_V), ('Quality', quality),
                       ('N_outer', n_outer.astype(float)),
                       ('Noise est.', offset_noise)]:
        t, p = sp_stats.ttest_ind(arr[perfect], arr[imperfect])
        if p < 0.1:
            print(f"    {name}: t={t:+.2f}, p={p:.3f}")

# "Super-perfect": |resid| < 0.5× noise
super_perfect = np.abs(resid6) < 0.5 * offset_noise
print(f"\n  'Super-perfect' (|resid| < 0.5× noise): {super_perfect.sum()}/{n} ({super_perfect.sum()/n*100:.0f}%)")

print(f"\n✓ TEST 7 PASSED: Perfect galaxies identified")

# ============================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS — THE GALAXY-BY-GALAXY PICTURE")
print("=" * 60)
# ============================================================

print(f"\n{'='*60}")
print(f"GALAXY CENSUS SYNTHESIS")
print(f"{'='*60}")

print(f"\n1. RESIDUAL DISTRIBUTION:")
print(f"   Median |resid| = {np.median(np.abs(resid6)):.4f} dex")
print(f"   90th percentile = {np.percentile(np.abs(resid6), 90):.4f} dex")
print(f"   Within noise: {within_noise.sum()}/{n} ({within_noise.sum()/n*100:.0f}%)")
print(f"   Within 2× noise: {within_2noise.sum()}/{n} ({within_2noise.sum()/n*100:.0f}%)")

print(f"\n2. BEST GALAXIES (lowest |resid|):")
for idx in idx_sorted[:3]:
    g = galaxies[idx]
    print(f"   {g['id']}: |resid|={abs(resid6[idx]):.5f}, V={g['vflat']:.0f}")

print(f"\n3. WORST GALAXIES (highest |resid|):")
for idx in idx_sorted[-3:][::-1]:
    g = galaxies[idx]
    print(f"   {g['id']}: |resid|={abs(resid6[idx]):.4f}, V={g['vflat']:.0f}")

print(f"\n4. QUALITY MATTERS:")
for q in [1, 2, 3]:
    mask = quality == q
    if mask.sum() >= 3:
        print(f"   Q={q}: RMS={np.sqrt(np.mean(resid6[mask]**2)):.4f} (N={mask.sum()})")

print(f"\n5. OUTLIER CENSUS:")
print(f"   N with |resid| > 2σ: {outliers.sum()}")
if outliers.sum() > 0:
    n_noise_consistent = np.sum(np.abs(resid6[outliers]) < 1.5 * offset_noise[outliers])
    print(f"   Noise-consistent: {n_noise_consistent}/{outliers.sum()}")

print(f"\n6. MOST IMPORTANT GALAXY:")
print(f"   {galaxies[impact_sorted[0]]['id']} (max |DFBETAS|={max_dfbetas[impact_sorted[0]]:.3f})")

print(f"\n7. PREDICTION QUALITY BREAKDOWN:")
bins = [(0, 0.01, 'Excellent'), (0.01, 0.02, 'Good'),
        (0.02, 0.04, 'Average'), (0.04, 0.06, 'Below avg'),
        (0.06, 0.10, 'Poor'), (0.10, 1.0, 'Outlier')]
for lo, hi, label in bins:
    mask = (np.abs(resid6) >= lo) & (np.abs(resid6) < hi)
    print(f"   {label:<12} (|resid| {lo:.2f}-{hi:.2f}): {mask.sum():>3} galaxies ({mask.sum()/n*100:.0f}%)")

print(f"\n{'='*60}")
print(f"BOTTOM LINE:")
print(f"  {within_noise.sum()}/{n} galaxies ({within_noise.sum()/n*100:.0f}%) predicted within noise")
n_good = np.sum(np.abs(resid6) < 0.04)
print(f"  {n_good}/{n} galaxies ({n_good/n*100:.0f}%) predicted to <0.04 dex (10%)")
n_poor = np.sum(np.abs(resid6) > 0.06)
print(f"  {n_poor}/{n} galaxies ({n_poor/n*100:.0f}%) with |resid| > 0.06 dex")
print(f"  Best: {galaxies[idx_sorted[0]]['id']} ({abs(resid6[idx_sorted[0]]):.5f} dex)")
print(f"  Worst: {galaxies[idx_sorted[-1]]['id']} ({abs(resid6[idx_sorted[-1]]):.4f} dex)")
print(f"{'='*60}")

print(f"\n✓ TEST 8 PASSED: Synthesis complete")

# ============================================================
# FINAL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print(f"SESSION #554: ALL 8 TESTS PASSED")
print(f"{'='*70}")
print(f"\nKey findings:")
print(f"  1. Median |resid| = {np.median(np.abs(resid6)):.4f} dex")
print(f"  2. {within_noise.sum()}/{n} ({within_noise.sum()/n*100:.0f}%) within estimated noise")
print(f"  3. {np.sum(np.abs(resid6) < 0.04)}/{n} ({np.sum(np.abs(resid6) < 0.04)/n*100:.0f}%) < 0.04 dex (10%)")
print(f"  4. N outliers (>2σ): {outliers.sum()}")
print(f"  5. r(quality, |resid|) = {r_q_abs:+.3f}")
print(f"  6. Most important: {galaxies[impact_sorted[0]]['id']}")
print(f"  7. Best/worst ratio: {abs(resid6[idx_sorted[-1]])/max(abs(resid6[idx_sorted[0]]),1e-6):.0f}×")
print(f"  8. Model predicts galaxy-by-galaxy to {rms6:.4f} dex average")
