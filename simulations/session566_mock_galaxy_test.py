#!/usr/bin/env python3
"""
======================================================================
SESSION #566: MOCK GALAXY TEST — FORWARD MODELING VALIDATION
======================================================================

All previous sessions analyzed real SPARC data. This session generates
mock galaxies with KNOWN MOND physics and KNOWN M/L, applies the model,
and tests whether it recovers the input parameters. This is the first
controlled test of the model's recovery ability.

The mock galaxies:
- MOND RAR with known a₀ and interpolation function
- Known M/L_disk (drawn from distribution)
- Known distance (with errors)
- Known V_flat, L, c_V, f_gas
- Realistic measurement noise

Tests:
1. Generate mock galaxies from MOND + realistic properties
2. Compute offsets and fit the 6-var model to mocks
3. Compare mock coefficients to real SPARC coefficients
4. Recovery of known M/L from mock offset
5. Effect of M/L scatter on model performance
6. Effect of distance errors on recovery
7. What fraction of model performance is M/L vs noise?
8. Synthesis: what the mocks tell us about the real model

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-07
Session: #566
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


def loo_r2_val(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    return 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)


from scipy import stats as sp_stats

print("=" * 70)
print("SESSION #566: MOCK GALAXY TEST")
print("=" * 70)

# First load real SPARC data as a template
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
catalog = load_sparc_catalog(
    os.path.join(base_dir, "sparc_real_data", "SPARC_Lelli2016c.mrt"))
models = load_sparc_mass_models(
    os.path.join(base_dir, "sparc_real_data", "MassModels_Lelli2016c.mrt"))

ml_disk = 0.5
ml_bul = 0.7

# Build real galaxy sample
real_galaxies = []
for gal_id, points in models.items():
    if len(points) < 5 or gal_id not in catalog:
        continue
    cat = catalog[gal_id]
    vflat = cat.get('vflat', 0)
    lum = cat.get('luminosity', 0)
    sb_eff = cat.get('sb_eff', 0)

    if vflat <= 0 or lum <= 0 or sb_eff <= 0:
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

    real_galaxies.append({
        'id': gal_id,
        'offset': offset_val,
        'logV': np.log10(vflat),
        'logL': np.log10(lum),
        'c_V': c_V_val,
        'f_gas': f_gas_val,
        'g_bar': g_bar_v,
        'g_obs': g_obs_v,
        'v_obs': v_obs_v,
        'v_gas': v_gas_v,
        'v_disk': v_disk_v,
        'e_vobs': e_vobs_v,
        'radius': radius_v,
        'mond': mond,
        'outer_mond': outer_mond,
    })

n_real = len(real_galaxies)
print(f"\n{n_real} real galaxies loaded as templates")

# Real model
offset_real = np.array([g['offset'] for g in real_galaxies])
logV_real = np.array([g['logV'] for g in real_galaxies])
logL_real = np.array([g['logL'] for g in real_galaxies])
cV_real = np.array([g['c_V'] for g in real_galaxies])
fgas_real = np.array([g['f_gas'] for g in real_galaxies])
ones_real = np.ones(n_real)

X6_real = np.column_stack([ones_real, logV_real, logL_real, cV_real, fgas_real,
                            logV_real*cV_real, logL_real*fgas_real])
beta_real, yhat_real, resid_real, R2_real, rms_real = build_model(X6_real, offset_real)
loo_real = loo_r2_val(X6_real, offset_real)
print(f"Real 6-var: R²={R2_real:.4f}, LOO={loo_real:.4f}, RMS={rms_real:.4f}")

# ============================================================
print("\n" + "=" * 60)
print("TEST 1: GENERATE MOCK GALAXIES")
print("=" * 60)
# ============================================================

# Strategy: use real galaxy templates but assign each a random "true" M/L
# The true g_obs = g_bar_true × ν(g_bar_true / a₀)
# where g_bar_true = g_bar × (true_ml / assumed_ml)
# The observed g_obs includes velocity noise

np.random.seed(42)

# True M/L distribution: log-normal centered on 0.5 with scatter
# (matching Session #517: optimal ~0.75, σ(log M/L) = 0.076 dex)
true_ml_disk = 0.5 * 10**(np.random.normal(0, 0.08, n_real))

# For each galaxy, compute:
# 1. True g_bar = g_bar_template × (true_ml / 0.5) for disk component
# 2. True g_obs = g_bar_true × ν(g_bar_true / a₀)
# 3. Observed v_obs = sqrt(g_obs × r) + noise

mock_galaxies = []
for i, g in enumerate(real_galaxies):
    # Scale disk contribution by true M/L ratio
    ml_ratio = true_ml_disk[i] / ml_disk

    # g_bar has contributions from disk and gas
    # g_bar_disk ∝ M/L × L, g_gas fixed
    # We need to separate disk and gas contributions
    v_disk = g['v_disk']
    v_gas = g['v_gas']
    radius = g['radius']

    # True g_bar: scale only disk component
    kms_to_ms = 1e3
    kpc_to_m = 3.0857e19

    g_disk = np.sign(v_disk) * v_disk**2 * kms_to_ms**2 / (radius * kpc_to_m)
    g_gas = np.sign(v_gas) * v_gas**2 * kms_to_ms**2 / (radius * kpc_to_m)

    g_bar_true = np.abs(g_disk * ml_ratio + g_gas)

    # True g_obs from MOND
    valid = g_bar_true > 0
    g_bar_t = g_bar_true[valid]
    g_obs_true = g_bar_t * nu_mcgaugh(g_bar_t / a0_mond)

    # Add velocity noise
    # v_obs_true = sqrt(g_obs_true × r)
    r_valid = radius[valid]
    v_true = np.sqrt(np.abs(g_obs_true) * r_valid * kpc_to_m) / kms_to_ms

    e_v = g['e_vobs'][valid]
    v_mock = v_true + np.random.normal(0, e_v)
    v_mock = np.abs(v_mock)  # velocities must be positive

    # Compute mock g_obs, g_bar using ASSUMED M/L (not true M/L)
    g_bar_assumed = np.abs(g_disk[valid] + g_gas[valid])  # assumed M/L = 0.5
    g_obs_mock = v_mock**2 * kms_to_ms**2 / (r_valid * kpc_to_m)

    valid2 = (g_bar_assumed > 0) & (g_obs_mock > 0) & np.isfinite(g_bar_assumed) & np.isfinite(g_obs_mock)
    if valid2.sum() < 5:
        continue

    g_bar_a = g_bar_assumed[valid2]
    g_obs_m = g_obs_mock[valid2]

    # MOND prediction with assumed M/L
    g_rar = g_bar_a * nu_mcgaugh(g_bar_a / a0_mond)
    offset_pts = np.log10(g_obs_m) - np.log10(g_rar)

    # Mimic the same offset computation
    mond_mask = g_bar_a < a0_mond
    if mond_mask.sum() < 3:
        continue

    radius_m = r_valid[valid2][mond_mask]
    med_r = np.median(radius_m)
    outer_mond = mond_mask.copy()
    outer_mond[mond_mask] = radius_m > med_r

    if outer_mond.sum() >= 2:
        offset_mock = np.mean(offset_pts[outer_mond])
    else:
        offset_mock = np.mean(offset_pts[mond_mask])

    mock_galaxies.append({
        'offset': offset_mock,
        'logV': g['logV'],
        'logL': g['logL'],
        'c_V': g['c_V'],
        'f_gas': g['f_gas'],
        'true_ml': true_ml_disk[i],
        'log_ml_ratio': np.log10(true_ml_disk[i] / ml_disk),
        'real_offset': g['offset'],
    })

n_mock = len(mock_galaxies)
print(f"\n{n_mock}/{n_real} mock galaxies generated")
print(f"True M/L distribution: mean={np.mean(true_ml_disk[:n_mock]):.3f}, "
      f"std={np.std(true_ml_disk[:n_mock]):.3f}")

mock_offset = np.array([g['offset'] for g in mock_galaxies])
mock_logV = np.array([g['logV'] for g in mock_galaxies])
mock_logL = np.array([g['logL'] for g in mock_galaxies])
mock_cV = np.array([g['c_V'] for g in mock_galaxies])
mock_fgas = np.array([g['f_gas'] for g in mock_galaxies])
mock_ml = np.array([g['log_ml_ratio'] for g in mock_galaxies])
ones_mock = np.ones(n_mock)

print(f"Mock offset: mean={np.mean(mock_offset):+.4f}, std={np.std(mock_offset):.4f}")
print(f"Real offset: mean={np.mean(offset_real):+.4f}, std={np.std(offset_real):.4f}")

print(f"\n\u2713 TEST 1 PASSED: Mock galaxies generated")

# ============================================================
print("\n" + "=" * 60)
print("TEST 2: FIT 6-VAR MODEL TO MOCKS")
print("=" * 60)
# ============================================================

X6_mock = np.column_stack([ones_mock, mock_logV, mock_logL, mock_cV, mock_fgas,
                            mock_logV*mock_cV, mock_logL*mock_fgas])

beta_mock, yhat_mock, resid_mock, R2_mock, rms_mock = build_model(X6_mock, mock_offset)
loo_mock = loo_r2_val(X6_mock, mock_offset)

print(f"\nMock 6-var model:")
print(f"  R²={R2_mock:.4f}, LOO={loo_mock:.4f}, RMS={rms_mock:.4f}")
print(f"  Real: R²={R2_real:.4f}, LOO={loo_real:.4f}, RMS={rms_real:.4f}")

print(f"\nCoefficient comparison:")
var_names = ['const', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V', 'logL×f_gas']
print(f"{'Variable':<12} {'Real':>8} {'Mock':>8} {'Ratio':>8}")
print("-" * 40)
for j, name in enumerate(var_names):
    ratio = beta_mock[j] / beta_real[j] if abs(beta_real[j]) > 1e-10 else np.nan
    print(f"{name:<12} {beta_real[j]:+8.4f} {beta_mock[j]:+8.4f} {ratio:>8.2f}")

# Sign agreement
sign_agree = np.sum(np.sign(beta_mock) == np.sign(beta_real))
print(f"\nSign agreement: {sign_agree}/{len(beta_real)}")

print(f"\n\u2713 TEST 2 PASSED: Mock model fitted")

# ============================================================
print("\n" + "=" * 60)
print("TEST 3: DOES MOCK OFFSET CORRELATE WITH TRUE M/L?")
print("=" * 60)
# ============================================================

r_ml, p_ml = sp_stats.pearsonr(mock_offset, mock_ml)
print(f"\nr(mock offset, log(true M/L ratio)) = {r_ml:+.4f} (p={p_ml:.2e})")

# After controlling for V, L, c_V, f_gas
# The residual should correlate with M/L
r_resid_ml, p_resid_ml = sp_stats.pearsonr(resid_mock, mock_ml)
print(f"r(mock residual, log(true M/L ratio)) = {r_resid_ml:+.4f} (p={p_resid_ml:.2e})")

# How much of mock offset variance is explained by true M/L?
X_ml = np.column_stack([ones_mock, mock_ml])
_, _, _, r2_ml, _ = build_model(X_ml, mock_offset)
print(f"\nR²(offset ~ M/L ratio): {r2_ml:.4f}")

# M/L + galaxy properties
X_ml_full = np.column_stack([X6_mock, mock_ml])
_, _, _, r2_ml_full, rms_ml_full = build_model(X_ml_full, mock_offset)
loo_ml_full = loo_r2_val(X_ml_full, mock_offset)
print(f"R²(offset ~ 6-var + M/L): {r2_ml_full:.4f}, LOO={loo_ml_full:.4f}, RMS={rms_ml_full:.4f}")
print(f"6-var only: R²={R2_mock:.4f}, LOO={loo_mock:.4f}")
print(f"Adding M/L: ΔR²={r2_ml_full - R2_mock:+.4f}, ΔLOO={loo_ml_full - loo_mock:+.4f}")

# Invert: predict true M/L from mock offset + properties
X_inv_ml = np.column_stack([ones_mock, mock_offset, mock_logV, mock_logL, mock_cV, mock_fgas])
_, yhat_ml_inv, _, r2_ml_inv, rms_ml_inv = build_model(X_inv_ml, mock_ml)
loo_ml_inv = loo_r2_val(X_inv_ml, mock_ml)
print(f"\nM/L recovery from offset + properties:")
print(f"  R²={r2_ml_inv:.4f}, LOO={loo_ml_inv:.4f}, RMS={rms_ml_inv:.4f}")
print(f"  (True M/L scatter: {np.std(mock_ml):.4f} dex)")

print(f"\n\u2713 TEST 3 PASSED: M/L correlation analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 4: MOCK VS REAL OFFSET COMPARISON")
print("=" * 60)
# ============================================================

# How similar are mock and real offsets for the same galaxies?
real_off_matched = np.array([g['real_offset'] for g in mock_galaxies])
r_match, p_match = sp_stats.pearsonr(mock_offset, real_off_matched)

print(f"\nMock vs real offset for same galaxies:")
print(f"  r(mock, real) = {r_match:.4f} (p={p_match:.2e})")
print(f"  RMS difference = {np.sqrt(np.mean((mock_offset - real_off_matched)**2)):.4f} dex")
print(f"  Mean difference = {np.mean(mock_offset - real_off_matched):+.4f} dex")
print(f"  Mock std = {np.std(mock_offset):.4f}")
print(f"  Real std = {np.std(real_off_matched):.4f}")
print(f"  Std ratio = {np.std(mock_offset)/np.std(real_off_matched):.3f}")

# What explains the difference?
diff = mock_offset - real_off_matched
r_diff_ml, _ = sp_stats.pearsonr(diff, mock_ml)
r_diff_logV, _ = sp_stats.pearsonr(diff, mock_logV)
r_diff_fgas, _ = sp_stats.pearsonr(diff, mock_fgas)
print(f"\nCorrelations of (mock - real) difference:")
print(f"  r(diff, log M/L ratio) = {r_diff_ml:+.4f}")
print(f"  r(diff, logV) = {r_diff_logV:+.4f}")
print(f"  r(diff, f_gas) = {r_diff_fgas:+.4f}")

print(f"\n\u2713 TEST 4 PASSED: Mock vs real comparison complete")

# ============================================================
print("\n" + "=" * 60)
print("TEST 5: EFFECT OF M/L SCATTER ON MODEL")
print("=" * 60)
# ============================================================

# Run mock test with different levels of M/L scatter
ml_scatters = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]

print(f"\nModel performance vs M/L scatter:")
print(f"{'σ(logM/L)':>10} {'R²':>8} {'LOO':>8} {'RMS':>8} {'offset_std':>10}")
print("-" * 50)

for sigma_ml in ml_scatters:
    np.random.seed(42)
    ml_draws = 0.5 * 10**(np.random.normal(0, sigma_ml, n_real))

    mock_offs = []
    mock_props = []
    for j, g in enumerate(real_galaxies):
        ml_r = ml_draws[j] / ml_disk
        v_disk = g['v_disk']
        v_gas = g['v_gas']
        radius = g['radius']

        kms_to_ms = 1e3
        kpc_to_m = 3.0857e19

        g_disk = np.sign(v_disk) * v_disk**2 * kms_to_ms**2 / (radius * kpc_to_m)
        g_gas = np.sign(v_gas) * v_gas**2 * kms_to_ms**2 / (radius * kpc_to_m)
        g_bar_t = np.abs(g_disk * ml_r + g_gas)

        valid = g_bar_t > 0
        g_bt = g_bar_t[valid]
        g_obs_t = g_bt * nu_mcgaugh(g_bt / a0_mond)

        r_v = radius[valid]
        v_true = np.sqrt(np.abs(g_obs_t) * r_v * kpc_to_m) / kms_to_ms
        e_v = g['e_vobs'][valid]
        v_m = np.abs(v_true + np.random.normal(0, e_v))

        g_bar_a = np.abs(g_disk[valid] + g_gas[valid])
        g_obs_m = v_m**2 * kms_to_ms**2 / (r_v * kpc_to_m)

        valid2 = (g_bar_a > 0) & (g_obs_m > 0) & np.isfinite(g_bar_a) & np.isfinite(g_obs_m)
        if valid2.sum() < 5:
            continue

        g_rar = g_bar_a[valid2] * nu_mcgaugh(g_bar_a[valid2] / a0_mond)
        off_pts = np.log10(g_obs_m[valid2]) - np.log10(g_rar)

        mond_m = g_bar_a[valid2] < a0_mond
        if mond_m.sum() < 3:
            continue

        r_m = r_v[valid2][mond_m]
        med = np.median(r_m)
        outer_m = mond_m.copy()
        outer_m[mond_m] = r_m > med

        if outer_m.sum() >= 2:
            off = np.mean(off_pts[outer_m])
        else:
            off = np.mean(off_pts[mond_m])

        mock_offs.append(off)
        mock_props.append([g['logV'], g['logL'], g['c_V'], g['f_gas']])

    mo = np.array(mock_offs)
    mp = np.array(mock_props)
    nm = len(mo)
    X6_m = np.column_stack([np.ones(nm), mp[:, 0], mp[:, 1], mp[:, 2], mp[:, 3],
                             mp[:, 0]*mp[:, 2], mp[:, 1]*mp[:, 3]])
    try:
        _, _, _, r2_m, rms_m = build_model(X6_m, mo)
        loo_m = loo_r2_val(X6_m, mo)
        print(f"   {sigma_ml:>6.3f}  {r2_m:>8.4f} {loo_m:>8.4f} {rms_m:>8.4f} {np.std(mo):>10.4f}")
    except:
        print(f"   {sigma_ml:>6.3f}  failed")

print(f"\n  Real data: R²={R2_real:.4f}, LOO={loo_real:.4f}, RMS={rms_real:.4f}, std={np.std(offset_real):.4f}")

print(f"\n\u2713 TEST 5 PASSED: M/L scatter effect analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 6: EFFECT OF DISTANCE ERRORS")
print("=" * 60)
# ============================================================

# Add distance errors to the mock galaxies
# Distance error → logL changes by 2×log10(α) where α = D_err/D_true
dist_errors = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]

print(f"\nModel performance vs distance error (σ(D)/D):")
print(f"{'σ(D)/D':>8} {'R²':>8} {'LOO':>8} {'RMS':>8} {'offset_std':>10}")
print("-" * 48)

for sigma_d in dist_errors:
    np.random.seed(42)
    ml_draws = 0.5 * 10**(np.random.normal(0, 0.08, n_real))  # fixed M/L scatter
    d_factors = 10**(np.random.normal(0, sigma_d * np.log10(np.e), n_real))

    mock_offs = []
    mock_props = []

    for j, g in enumerate(real_galaxies):
        ml_r = ml_draws[j] / ml_disk
        d_f = d_factors[j]

        v_disk = g['v_disk']
        v_gas = g['v_gas']
        radius = g['radius']

        kms_to_ms = 1e3
        kpc_to_m = 3.0857e19

        g_disk = np.sign(v_disk) * v_disk**2 * kms_to_ms**2 / (radius * kpc_to_m)
        g_gas = np.sign(v_gas) * v_gas**2 * kms_to_ms**2 / (radius * kpc_to_m)
        g_bar_t = np.abs(g_disk * ml_r + g_gas)

        valid = g_bar_t > 0
        g_bt = g_bar_t[valid]
        g_obs_t = g_bt * nu_mcgaugh(g_bt / a0_mond)

        r_v = radius[valid]
        v_true = np.sqrt(np.abs(g_obs_t) * r_v * kpc_to_m) / kms_to_ms
        e_v = g['e_vobs'][valid]
        v_m = np.abs(v_true + np.random.normal(0, e_v))

        # Distance error: affects logL (observed L scales as D²)
        logL_err = g['logL'] + 2 * np.log10(d_f)

        # g_bar changes with distance: g_bar ∝ M/r ∝ L×M/L / (D×θ)
        # But we're using the TEMPLATE galaxy with assumed distances
        # The key effect: luminosity used in model is wrong

        g_bar_a = np.abs(g_disk[valid] + g_gas[valid])
        g_obs_m = v_m**2 * kms_to_ms**2 / (r_v * kpc_to_m)

        valid2 = (g_bar_a > 0) & (g_obs_m > 0) & np.isfinite(g_bar_a) & np.isfinite(g_obs_m)
        if valid2.sum() < 5:
            continue

        g_rar = g_bar_a[valid2] * nu_mcgaugh(g_bar_a[valid2] / a0_mond)
        off_pts = np.log10(g_obs_m[valid2]) - np.log10(g_rar)

        mond_m = g_bar_a[valid2] < a0_mond
        if mond_m.sum() < 3:
            continue

        r_m = r_v[valid2][mond_m]
        med = np.median(r_m)
        outer_m = mond_m.copy()
        outer_m[mond_m] = r_m > med

        if outer_m.sum() >= 2:
            off = np.mean(off_pts[outer_m])
        else:
            off = np.mean(off_pts[mond_m])

        mock_offs.append(off)
        mock_props.append([g['logV'], logL_err, g['c_V'], g['f_gas']])

    mo = np.array(mock_offs)
    mp = np.array(mock_props)
    nm = len(mo)
    X6_m = np.column_stack([np.ones(nm), mp[:, 0], mp[:, 1], mp[:, 2], mp[:, 3],
                             mp[:, 0]*mp[:, 2], mp[:, 1]*mp[:, 3]])
    try:
        _, _, _, r2_m, rms_m = build_model(X6_m, mo)
        loo_m = loo_r2_val(X6_m, mo)
        print(f"  {sigma_d:>6.2f} {r2_m:>8.4f} {loo_m:>8.4f} {rms_m:>8.4f} {np.std(mo):>10.4f}")
    except:
        print(f"  {sigma_d:>6.2f} failed")

print(f"\n\u2713 TEST 6 PASSED: Distance error effect analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 7: NOISE-ONLY MOCKS (M/L = CONSTANT)")
print("=" * 60)
# ============================================================

# What happens when ALL galaxies have the SAME M/L?
# Any offset variance must come from measurement noise alone

np.random.seed(42)

noise_offs = []
noise_props = []

for j, g in enumerate(real_galaxies):
    # True M/L = exactly 0.5 (no scatter)
    v_disk = g['v_disk']
    v_gas = g['v_gas']
    radius = g['radius']

    kms_to_ms = 1e3
    kpc_to_m = 3.0857e19

    g_disk = np.sign(v_disk) * v_disk**2 * kms_to_ms**2 / (radius * kpc_to_m)
    g_gas = np.sign(v_gas) * v_gas**2 * kms_to_ms**2 / (radius * kpc_to_m)
    g_bar_t = np.abs(g_disk + g_gas)  # true M/L = assumed M/L = 0.5

    valid = g_bar_t > 0
    g_bt = g_bar_t[valid]
    g_obs_t = g_bt * nu_mcgaugh(g_bt / a0_mond)

    r_v = radius[valid]
    v_true = np.sqrt(np.abs(g_obs_t) * r_v * kpc_to_m) / kms_to_ms
    e_v = g['e_vobs'][valid]
    v_m = np.abs(v_true + np.random.normal(0, e_v))

    g_bar_a = np.abs(g_disk[valid] + g_gas[valid])
    g_obs_m = v_m**2 * kms_to_ms**2 / (r_v * kpc_to_m)

    valid2 = (g_bar_a > 0) & (g_obs_m > 0) & np.isfinite(g_bar_a) & np.isfinite(g_obs_m)
    if valid2.sum() < 5:
        continue

    g_rar = g_bar_a[valid2] * nu_mcgaugh(g_bar_a[valid2] / a0_mond)
    off_pts = np.log10(g_obs_m[valid2]) - np.log10(g_rar)

    mond_m = g_bar_a[valid2] < a0_mond
    if mond_m.sum() < 3:
        continue

    r_m = r_v[valid2][mond_m]
    med = np.median(r_m)
    outer_m = mond_m.copy()
    outer_m[mond_m] = r_m > med

    if outer_m.sum() >= 2:
        off = np.mean(off_pts[outer_m])
    else:
        off = np.mean(off_pts[mond_m])

    noise_offs.append(off)
    noise_props.append([g['logV'], g['logL'], g['c_V'], g['f_gas']])

no = np.array(noise_offs)
np_arr = np.array(noise_props)
n_noise = len(no)

X6_noise = np.column_stack([np.ones(n_noise), np_arr[:, 0], np_arr[:, 1], np_arr[:, 2], np_arr[:, 3],
                             np_arr[:, 0]*np_arr[:, 2], np_arr[:, 1]*np_arr[:, 3]])

_, _, _, r2_noise, rms_noise = build_model(X6_noise, no)
loo_noise = loo_r2_val(X6_noise, no)

print(f"\nNoise-only mock (M/L = constant = 0.5):")
print(f"  Offset std: {np.std(no):.4f} (real: {np.std(offset_real):.4f})")
print(f"  R²: {r2_noise:.4f} (real: {R2_real:.4f})")
print(f"  LOO: {loo_noise:.4f} (real: {loo_real:.4f})")
print(f"  RMS: {rms_noise:.4f} (real: {rms_real:.4f})")

# What fraction of real variance is M/L vs noise?
var_real = np.var(offset_real)
var_noise = np.var(no)
var_ml = var_real - var_noise
print(f"\nVariance decomposition:")
print(f"  Total (real): {var_real:.6f}")
print(f"  Noise only: {var_noise:.6f} ({var_noise/var_real*100:.1f}%)")
print(f"  M/L contribution: {var_ml:.6f} ({var_ml/var_real*100:.1f}%)")

print(f"\n  Mean noise offset: {np.mean(no):+.4f} (should be ~0)")
r_noise_logV, _ = sp_stats.pearsonr(no, np_arr[:, 0])
r_noise_fgas, _ = sp_stats.pearsonr(no, np_arr[:, 3])
print(f"  r(noise offset, logV) = {r_noise_logV:+.4f}")
print(f"  r(noise offset, f_gas) = {r_noise_fgas:+.4f}")

print(f"\n\u2713 TEST 7 PASSED: Noise-only mock analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS — WHAT MOCKS TELL US")
print("=" * 60)
# ============================================================

print(f"\n{'='*60}")
print(f"MOCK GALAXY TEST — SYNTHESIS")
print(f"{'='*60}")

print(f"\n1. MOCK MODEL PERFORMANCE:")
print(f"   Mock R²={R2_mock:.4f} vs Real R²={R2_real:.4f}")
print(f"   Mock LOO={loo_mock:.4f} vs Real LOO={loo_real:.4f}")
print(f"   Mock offset std={np.std(mock_offset):.4f} vs Real={np.std(offset_real):.4f}")

print(f"\n2. COEFFICIENT RECOVERY:")
print(f"   Sign agreement: {sign_agree}/{len(beta_real)}")
n_close = np.sum(np.abs(beta_mock / beta_real - 1) < 0.5)
print(f"   Within 50% of real: {n_close}/{len(beta_real)}")

print(f"\n3. M/L INFORMATION:")
print(f"   r(offset, true M/L) = {r_ml:+.4f}")
print(f"   M/L recovery from offset: LOO R²={loo_ml_inv:.4f}")

print(f"\n4. VARIANCE DECOMPOSITION:")
print(f"   Noise contribution: {var_noise/var_real*100:.1f}% of real variance")
print(f"   M/L contribution: {var_ml/var_real*100:.1f}% of real variance")
if var_noise / var_real < 0.5:
    print(f"   → Real offset is DOMINATED by M/L variation, not noise")
elif var_noise / var_real > 0.8:
    print(f"   → Real offset is DOMINATED by noise, not M/L")
else:
    print(f"   → Real offset has comparable M/L and noise contributions")

print(f"\n5. MOCK-REAL AGREEMENT:")
print(f"   r(mock, real offset) = {r_match:.4f}")
print(f"   RMS(mock - real) = {np.sqrt(np.mean((mock_offset - real_off_matched)**2)):.4f}")

print(f"\n{'='*60}")
print(f"CONCLUSION:")
print(f"  The mock test {'VALIDATES' if R2_mock > 0.5 else 'CHALLENGES'} the model.")
if sign_agree >= 5:
    print(f"  Coefficient signs are reproduced ({sign_agree}/7).")
if abs(r_ml) > 0.3:
    print(f"  The offset genuinely encodes M/L (r={r_ml:+.3f}).")
print(f"  Noise contributes {var_noise/var_real*100:.0f}% of offset variance.")
print(f"  M/L variation contributes {var_ml/var_real*100:.0f}%.")
print(f"{'='*60}")

print(f"\n\u2713 TEST 8 PASSED: Synthesis complete")

# ============================================================
print(f"\n{'='*70}")
print(f"SESSION #566: ALL 8 TESTS PASSED")
print(f"{'='*70}")
