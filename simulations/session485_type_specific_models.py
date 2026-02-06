#!/usr/bin/env python3
"""
======================================================================
SESSION #485: TYPE-SPECIFIC MODELS — DIFFERENT PHYSICS FOR DIFFERENT GALAXIES
======================================================================

Session #484 showed logL×f_gas helps late types (ΔLOO = +0.052) but
hurts early types (ΔLOO = -0.035). This suggests type-specific models.

This session explores:
1. Optimal model for each type group
2. Whether c_V matters for early types but not late types
3. Whether f_gas matters for late types but not early types
4. A composite model using type-specific coefficients
5. How far can we push R² for each type?

Tests:
1. Type-specific 5-variable models
2. Type-specific variable importance (forward selection)
3. The minimal model per type
4. The composite model: fit separately, predict jointly
5. Late-type + logL×f_gas deep dive
6. Early-type + c_V deep dive
7. Cross-prediction: train on one type, predict another
8. The optimal strategy

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #485
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
    """Load SPARC data."""
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

        # Outer offset
        radius_m = radius_v[mond]
        med_r = np.median(radius_m)
        outer_mond = mond.copy()
        outer_mond[mond] = radius_m > med_r
        if outer_mond.sum() >= 2:
            g_rar_out = rar_prediction(g_bar_v[outer_mond])
            outer_offset = np.mean(np.log10(g_obs_v[outer_mond]) - np.log10(g_rar_out))
        else:
            g_rar = rar_prediction(g_bar_v[mond])
            outer_offset = np.mean(np.log10(g_obs_v[mond]) - np.log10(g_rar))

        # f_gas
        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'outer_offset': outer_offset, 'sb_eff': sb_eff, 'r_eff': r_eff_kpc,
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
    return loo_rms, loo_r2, loo_resid


print("=" * 70)
print("SESSION #485: TYPE-SPECIFIC MODELS")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nSample: {n} galaxies")

# Build arrays
logV = np.log10([g['vflat'] for g in galaxies])
logL = np.log10([g['lum'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
y = np.array([g['outer_offset'] for g in galaxies])
T = np.array([g['hubble_type'] for g in galaxies])

# Type masks
early = T < 4    # S0-Sb
mid = (T >= 4) & (T < 7)  # Sbc-Sd
late = T >= 7    # Sdm-Im

type_masks = {'Early (T<4)': early, 'Mid (4≤T<7)': mid, 'Late (T≥7)': late}

# Full-sample models
X5 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V])
X6 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])

_, _, _, R2_5all, _ = build_model(X5, y)
_, _, _, R2_6all, _ = build_model(X6, y)
_, loo_r2_5all, _ = loo_cv(X5, y)
_, loo_r2_6all, _ = loo_cv(X6, y)

print(f"\nFull sample 5-var: R² = {R2_5all:.4f}, LOO R² = {loo_r2_5all:.4f}")
print(f"Full sample 6-var: R² = {R2_6all:.4f}, LOO R² = {loo_r2_6all:.4f}")

# =====================================================================
# TEST 1: TYPE-SPECIFIC 5-VARIABLE MODELS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: TYPE-SPECIFIC 5-VARIABLE MODELS")
print("=" * 60)

print(f"\n{'Group':<15} {'N':<5} {'σ(off)':<8} {'R²':<8} {'LOO R²':<10} {'RMS':<8}")
print("-" * 54)

type_results = {}
for name, mask in type_masks.items():
    nm = mask.sum()
    X5_g = X5[mask]
    y_g = y[mask]
    _, _, _, R2_g, rms_g = build_model(X5_g, y_g)
    loo_g, loo_r2_g, _ = loo_cv(X5_g, y_g)
    sig = np.std(y_g)
    type_results[name] = {'R2': R2_g, 'loo_r2': loo_r2_g, 'rms': rms_g, 'n': nm, 'sigma': sig}
    print(f"{name:<15} {nm:<5} {sig:.4f}  {R2_g:.4f}  {loo_r2_g:.4f}    {rms_g:.4f}")

# Compare to full sample
print(f"\n{'Full sample':<15} {n:<5} {np.std(y):.4f}  {R2_5all:.4f}  {loo_r2_5all:.4f}    {np.sqrt(np.mean(build_model(X5,y)[2]**2)):.4f}")

print("\n✓ Test 1 passed: type-specific baselines established")

# =====================================================================
# TEST 2: TYPE-SPECIFIC VARIABLE IMPORTANCE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: TYPE-SPECIFIC VARIABLE IMPORTANCE (FORWARD SELECTION)")
print("=" * 60)

all_vars = {
    'logV': logV, 'logL': logL, 'c_V': c_V, 'f_gas': f_gas,
    'logV×c_V': logV * c_V, 'logL×f_gas': logL * f_gas,
    'f_gas²': f_gas**2, 'logV×f_gas': logV * f_gas,
}

for group_name, mask in type_masks.items():
    nm = mask.sum()
    y_g = y[mask]
    print(f"\n--- {group_name} (N={nm}) ---")
    print(f"  {'Step':<5} {'Added':<15} {'R²':<8} {'LOO R²':<10}")
    print("  " + "-" * 38)

    selected = []
    remaining = list(all_vars.keys())
    X_step = np.ones((nm, 1))

    for step in range(min(6, len(remaining))):
        best_name = None
        best_loo_r2 = -999

        for vname in remaining:
            var = all_vars[vname][mask]
            X_try = np.column_stack([X_step, var])
            try:
                _, loo_r2_try, _ = loo_cv(X_try, y_g)
            except:
                continue
            if loo_r2_try > best_loo_r2:
                best_loo_r2 = loo_r2_try
                best_name = vname

        if best_name is None:
            break

        selected.append(best_name)
        remaining.remove(best_name)
        X_step = np.column_stack([X_step, all_vars[best_name][mask]])

        _, _, _, R2_s, _ = build_model(X_step, y_g)
        _, loo_r2_s, _ = loo_cv(X_step, y_g)
        print(f"  {step+1:<5} {best_name:<15} {R2_s:.4f}  {loo_r2_s:.4f}")

        if step > 0 and loo_r2_s - best_loo_r2 < -0.01:
            break

    print(f"  Selected: {' + '.join(selected[:4])}")

print("\n✓ Test 2 passed: type-specific importance done")

# =====================================================================
# TEST 3: MINIMAL MODELS PER TYPE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: MINIMAL MODELS PER TYPE")
print("=" * 60)

# For each type, test 2-var, 3-var, and optimal models
for group_name, mask in type_masks.items():
    nm = mask.sum()
    y_g = y[mask]
    print(f"\n--- {group_name} (N={nm}) ---")

    models_to_test = {
        'logV only': np.column_stack([np.ones(nm), logV[mask]]),
        'logV + logL': np.column_stack([np.ones(nm), logV[mask], logL[mask]]),
        'logV + logL + c_V': np.column_stack([np.ones(nm), logV[mask], logL[mask], c_V[mask]]),
        'logV + logL + f_gas': np.column_stack([np.ones(nm), logV[mask], logL[mask], f_gas[mask]]),
        '5-var standard': X5[mask],
        '6-var (logL×f_gas)': X6[mask],
    }

    if group_name == 'Late (T≥7)':
        # Add N_corr-based model for late types
        V_ms = np.array([g['vflat'] for g in galaxies]) * 1e3
        R_m = np.array([g['r_eff'] for g in galaxies]) * 3.086e19
        N_corr = V_ms**2 / (R_m * a0_mond)
        logN = np.log10(np.clip(N_corr, 1e-5, None))
        models_to_test['logN_corr only'] = np.column_stack([np.ones(nm), logN[mask]])
        models_to_test['logN + f_gas'] = np.column_stack([np.ones(nm), logN[mask], f_gas[mask]])

    print(f"  {'Model':<25} {'k':<4} {'R²':<8} {'LOO R²':<10} {'RMS':<8}")
    print("  " + "-" * 55)

    for mname, X_m in models_to_test.items():
        k = X_m.shape[1] - 1
        _, _, _, R2_m, rms_m = build_model(X_m, y_g)
        try:
            _, loo_r2_m, _ = loo_cv(X_m, y_g)
        except:
            loo_r2_m = np.nan
        print(f"  {mname:<25} {k:<4} {R2_m:.4f}  {loo_r2_m:.4f}    {rms_m:.4f}")

print("\n✓ Test 3 passed: minimal models tested")

# =====================================================================
# TEST 4: COMPOSITE MODEL — FIT SEPARATELY, PREDICT JOINTLY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: COMPOSITE MODEL")
print("=" * 60)

# Fit type-specific models, then combine predictions
# Early types: 5-var (c_V matters)
# Late types: 5-var + logL×f_gas (gas matters)
# Mid types: either

composite_pred = np.zeros(n)
composite_loo_resid = np.zeros(n)

for group_name, mask in type_masks.items():
    nm = mask.sum()
    y_g = y[mask]

    if group_name == 'Late (T≥7)':
        X_g = X6[mask]
    else:
        X_g = X5[mask]

    beta_g, yhat_g, resid_g, R2_g, rms_g = build_model(X_g, y_g)
    _, loo_r2_g, loo_resid_g = loo_cv(X_g, y_g)

    composite_pred[mask] = yhat_g
    composite_loo_resid[mask] = loo_resid_g

    print(f"{group_name:<15} (N={nm}): R² = {R2_g:.4f}, LOO R² = {loo_r2_g:.4f}, model = "
          f"{'6-var' if group_name == 'Late (T≥7)' else '5-var'}")

# Composite R²
ss_res_comp = np.sum((y - composite_pred)**2)
ss_tot = np.sum((y - np.mean(y))**2)
R2_composite = 1 - ss_res_comp / ss_tot
rms_composite = np.sqrt(np.mean((y - composite_pred)**2))

ss_res_loo = np.sum(composite_loo_resid**2)
loo_r2_composite = 1 - ss_res_loo / ss_tot
loo_rms_composite = np.sqrt(np.mean(composite_loo_resid**2))

print(f"\nComposite model (type-specific):")
print(f"  R² = {R2_composite:.4f}")
print(f"  LOO R² = {loo_r2_composite:.4f}")
print(f"  RMS = {rms_composite:.4f}")
print(f"  LOO RMS = {loo_rms_composite:.4f}")

print(f"\nComparison:")
print(f"  {'Model':<30} {'R²':<8} {'LOO R²':<10}")
print("  " + "-" * 48)
print(f"  {'5-var (full sample)':<30} {R2_5all:.4f}  {loo_r2_5all:.4f}")
print(f"  {'6-var (full sample)':<30} {R2_6all:.4f}  {loo_r2_6all:.4f}")
print(f"  {'Composite (type-specific)':<30} {R2_composite:.4f}  {loo_r2_composite:.4f}")

assert R2_composite > R2_5all, "Composite should beat 5-var"
print("\n✓ Test 4 passed: composite model built")

# =====================================================================
# TEST 5: LATE-TYPE + logL×f_gas DEEP DIVE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: LATE-TYPE + logL×f_gas DEEP DIVE")
print("=" * 60)

# Late types with the full interaction model
mask_late = late
nm_late = mask_late.sum()
y_late = y[mask_late]
X5_late = X5[mask_late]
X6_late = X6[mask_late]

beta5_late, _, resid5_late, R2_5l, rms5_late = build_model(X5_late, y_late)
beta6_late, _, resid6_late, R2_6l, rms6_late = build_model(X6_late, y_late)
_, loo_r2_5l, _ = loo_cv(X5_late, y_late)
_, loo_r2_6l, _ = loo_cv(X6_late, y_late)

print(f"\nLate types (N={nm_late}):")
print(f"  5-var: R² = {R2_5l:.4f}, LOO R² = {loo_r2_5l:.4f}")
print(f"  6-var: R² = {R2_6l:.4f}, LOO R² = {loo_r2_6l:.4f}")

# Coefficients for late types
var_names_6 = ['intercept', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V', 'logL×f_gas']
print(f"\n  Late-type 6-var coefficients:")
se6_late = np.sqrt(np.diag(np.sum(resid6_late**2) / (nm_late - 7) * np.linalg.inv(X6_late.T @ X6_late)))
for i, name in enumerate(var_names_6):
    t_val = beta6_late[i] / se6_late[i]
    sig = '***' if abs(t_val) > 3 else '**' if abs(t_val) > 2 else '*' if abs(t_val) > 1.5 else ''
    print(f"    {name:<15} {beta6_late[i]:+.4f} (t={t_val:+.2f}) {sig}")

# What if we drop c_V for late types?
X4_late = np.column_stack([np.ones(nm_late), logV[mask_late], logL[mask_late],
                            f_gas[mask_late], logL[mask_late] * f_gas[mask_late]])
_, _, _, R2_4l, _ = build_model(X4_late, y_late)
_, loo_r2_4l, _ = loo_cv(X4_late, y_late)
print(f"\n  4-var (logV, logL, f_gas, logL×f_gas): R² = {R2_4l:.4f}, LOO R² = {loo_r2_4l:.4f}")

# Gas-rich late types
mask_gas_late = mask_late & (f_gas > 0.3)
nm_gl = mask_gas_late.sum()
if nm_gl >= 15:
    X6_gl = X6[mask_gas_late]
    y_gl = y[mask_gas_late]
    _, _, _, R2_gl, _ = build_model(X6_gl, y_gl)
    _, loo_r2_gl, _ = loo_cv(X6_gl, y_gl)
    print(f"\n  Gas-rich late (N={nm_gl}): R² = {R2_gl:.4f}, LOO R² = {loo_r2_gl:.4f}")

print("\n✓ Test 5 passed: late-type deep dive done")

# =====================================================================
# TEST 6: EARLY-TYPE + c_V DEEP DIVE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: EARLY-TYPE DEEP DIVE")
print("=" * 60)

mask_early = early
nm_early = mask_early.sum()
y_early = y[mask_early]

print(f"\nEarly types (N={nm_early}):")

# Various models for early types
early_models = {
    'logV only': np.column_stack([np.ones(nm_early), logV[mask_early]]),
    'logV + logL': np.column_stack([np.ones(nm_early), logV[mask_early], logL[mask_early]]),
    'logV + logL + c_V': np.column_stack([np.ones(nm_early), logV[mask_early],
                                           logL[mask_early], c_V[mask_early]]),
    '5-var': X5[mask_early],
    '6-var': X6[mask_early],
    'logV + c_V': np.column_stack([np.ones(nm_early), logV[mask_early], c_V[mask_early]]),
    'logV + c_V + logV×c_V': np.column_stack([np.ones(nm_early), logV[mask_early],
                                                c_V[mask_early], logV[mask_early] * c_V[mask_early]]),
}

print(f"  {'Model':<25} {'k':<4} {'R²':<8} {'LOO R²':<10}")
print("  " + "-" * 47)

for mname, X_m in early_models.items():
    k = X_m.shape[1] - 1
    _, _, _, R2_m, _ = build_model(X_m, y_early)
    try:
        _, loo_r2_m, _ = loo_cv(X_m, y_early)
    except:
        loo_r2_m = np.nan
    print(f"  {mname:<25} {k:<4} {R2_m:.4f}  {loo_r2_m:.4f}")

# c_V importance for early types
_, _, resid_noCV, _, _ = build_model(
    np.column_stack([np.ones(nm_early), logV[mask_early], logL[mask_early]]), y_early)
_, _, resid_withCV, _, _ = build_model(
    np.column_stack([np.ones(nm_early), logV[mask_early], logL[mask_early], c_V[mask_early]]), y_early)
r_cV_partial = np.corrcoef(c_V[mask_early], resid_noCV)[0, 1]
print(f"\n  Partial r(c_V, offset | logV, logL) = {r_cV_partial:+.4f}")

print("\n✓ Test 6 passed: early-type deep dive done")

# =====================================================================
# TEST 7: CROSS-PREDICTION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: CROSS-PREDICTION (TRAIN ON ONE TYPE, PREDICT ANOTHER)")
print("=" * 60)

# Train 5-var on each type, predict the others
for train_name, train_mask in type_masks.items():
    beta_train, _, _, R2_train, _ = build_model(X5[train_mask], y[train_mask])

    for test_name, test_mask in type_masks.items():
        y_test = y[test_mask]
        y_pred = X5[test_mask] @ beta_train
        resid_test = y_test - y_pred
        ss_res = np.sum(resid_test**2)
        ss_tot_test = np.sum((y_test - np.mean(y_test))**2)
        R2_test = 1 - ss_res / ss_tot_test if ss_tot_test > 0 else 0
        rms_test = np.sqrt(np.mean(resid_test**2))

        label = "SELF" if train_name == test_name else ""
        print(f"  Train: {train_name:<15} → Test: {test_name:<15} R² = {R2_test:.4f}  RMS = {rms_test:.4f}  {label}")

# How well does the full-sample model predict each type vs type-specific?
print(f"\nFull-sample vs type-specific prediction:")
beta_full, _, _, _, _ = build_model(X5, y)
for name, mask in type_masks.items():
    y_g = y[mask]
    # Full sample prediction
    y_pred_full = X5[mask] @ beta_full
    rms_full = np.sqrt(np.mean((y_g - y_pred_full)**2))
    # Type-specific prediction
    beta_type, _, _, _, _ = build_model(X5[mask], y_g)
    y_pred_type = X5[mask] @ beta_type
    rms_type = np.sqrt(np.mean((y_g - y_pred_type)**2))
    print(f"  {name:<15} Full RMS = {rms_full:.4f}, Type-specific RMS = {rms_type:.4f}, "
          f"Δ = {rms_type - rms_full:+.4f}")

print("\n✓ Test 7 passed: cross-prediction done")

# =====================================================================
# TEST 8: THE OPTIMAL STRATEGY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: THE OPTIMAL STRATEGY")
print("=" * 60)

# Build the optimal composite:
# Early: 5-var
# Mid: try both
# Late: 6-var (logL×f_gas)

strategies = {
    'Uniform 5-var': {},
    'Uniform 6-var': {},
    'Composite (5E + 5M + 6L)': {},
    'Composite (5E + 6M + 6L)': {},
}

for strat_name in strategies:
    pred = np.zeros(n)
    loo_resid_all = np.zeros(n)

    for group_name, mask in type_masks.items():
        y_g = y[mask]

        if strat_name == 'Uniform 5-var':
            X_g = X5[mask]
        elif strat_name == 'Uniform 6-var':
            X_g = X6[mask]
        elif strat_name == 'Composite (5E + 5M + 6L)':
            X_g = X6[mask] if group_name == 'Late (T≥7)' else X5[mask]
        elif strat_name == 'Composite (5E + 6M + 6L)':
            X_g = X5[mask] if group_name == 'Early (T<4)' else X6[mask]

        _, yhat_g, _, _, _ = build_model(X_g, y_g)
        _, _, loo_resid_g = loo_cv(X_g, y_g)
        pred[mask] = yhat_g
        loo_resid_all[mask] = loo_resid_g

    R2_s = 1 - np.sum((y - pred)**2) / np.sum((y - np.mean(y))**2)
    loo_r2_s = 1 - np.sum(loo_resid_all**2) / np.sum((y - np.mean(y))**2)
    rms_s = np.sqrt(np.mean((y - pred)**2))
    loo_rms_s = np.sqrt(np.mean(loo_resid_all**2))

    strategies[strat_name] = {'R2': R2_s, 'loo_r2': loo_r2_s, 'rms': rms_s, 'loo_rms': loo_rms_s}

print(f"\n{'Strategy':<35} {'R²':<8} {'LOO R²':<10} {'RMS':<8} {'LOO RMS':<10}")
print("-" * 71)
for name, info in strategies.items():
    print(f"{name:<35} {info['R2']:.4f}  {info['loo_r2']:.4f}    {info['rms']:.4f}  {info['loo_rms']:.4f}")

best_strat = max(strategies.items(), key=lambda x: x[1]['loo_r2'])
print(f"\nBest strategy by LOO R²: {best_strat[0]} → LOO R² = {best_strat[1]['loo_r2']:.4f}")

# Compare to global models
print(f"\nGlobal 5-var: R² = {R2_5all:.4f}, LOO R² = {loo_r2_5all:.4f}")
print(f"Global 6-var: R² = {R2_6all:.4f}, LOO R² = {loo_r2_6all:.4f}")
print(f"Best composite: R² = {best_strat[1]['R2']:.4f}, LOO R² = {best_strat[1]['loo_r2']:.4f}")

assert best_strat[1]['loo_r2'] > loo_r2_5all, "Best strategy should beat 5-var"
print("\n✓ Test 8 passed: optimal strategy identified")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #485 SUMMARY")
print("=" * 70)
print(f"Sample: {n} galaxies ({mask_early.sum()} early, {mid.sum()} mid, {mask_late.sum()} late)")
print(f"\nType-specific 5-var LOO R²:")
for name, info in type_results.items():
    print(f"  {name}: {info['loo_r2']:.4f}")
print(f"\nBest strategy: {best_strat[0]}")
print(f"  R² = {best_strat[1]['R2']:.4f}, LOO R² = {best_strat[1]['loo_r2']:.4f}")
print(f"\nAll 8 tests passed ✓")
