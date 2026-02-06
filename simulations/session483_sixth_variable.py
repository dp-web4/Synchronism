#!/usr/bin/env python3
"""
======================================================================
SESSION #483: THE SIXTH VARIABLE — CAN WE BREAK R² = 0.92?
======================================================================

The outer-only 5-variable model has R² = 0.913 (LOO R² = 0.898).
Session #482 found NN autocorrelation r = +0.46 in residuals, meaning
similar galaxies have similar residuals — structure remains.

This session systematically tests candidate 6th variables:
- Surface brightness (SB_eff)
- log(N_corr) = log(V²/(R×a₀))
- Effective radius (R_eff)
- Maximum radius (R_max)
- Hubble type (T)
- RC roughness
- Disk surface brightness (SB_disk)

And then tests 2-variable additions, non-linear terms, and the LOO
performance to avoid overfitting.

Tests:
1. Single variable screening: all candidates vs outer residual
2. 6-variable models: add each candidate to the 5-var model
3. LOO validation of top candidates
4. Two-variable additions (7-variable models)
5. Non-linear terms: squared and interaction effects
6. By galaxy type: which types benefit from 6th variable?
7. Stepwise selection: let the data choose
8. Synthesis: the optimal model

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #483
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
    """Load SPARC data with all needed quantities."""
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
        sb_disk = cat.get('sb_disk', 0)
        hubble_type = cat.get('hubble_type', 5)
        distance = cat.get('distance', 0)
        inclination = cat.get('inclination', 0)
        quality = cat.get('quality', 1)

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

        # Full + outer offset
        g_rar = rar_prediction(g_bar_v[mond])
        full_offset = np.mean(np.log10(g_obs_v[mond]) - np.log10(g_rar))

        radius_m = radius_v[mond]
        med_r = np.median(radius_m)
        outer_mond = mond.copy()
        outer_mond[mond] = radius_m > med_r
        if outer_mond.sum() >= 2:
            g_rar_out = rar_prediction(g_bar_v[outer_mond])
            outer_offset = np.mean(np.log10(g_obs_v[outer_mond]) - np.log10(g_rar_out))
        else:
            outer_offset = full_offset

        # f_gas
        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # RC roughness
        point_resids = np.log10(g_obs_v[mond]) - np.log10(g_rar)
        rc_roughness = np.std(point_resids)

        # N_corr
        V_ms = vflat * 1e3
        R_m = r_eff_kpc * 3.086e19
        N_corr = V_ms**2 / (R_m * a0_mond) if R_m > 0 else np.nan

        # R_max from RC
        R_max = radius_v.max()

        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'full_offset': full_offset, 'outer_offset': outer_offset,
            'distance': distance, 'inclination': inclination,
            'quality': quality, 'sb_eff': sb_eff, 'sb_disk': sb_disk,
            'r_eff': r_eff_kpc, 'R_max': R_max,
            'n_mond': mond.sum(), 'rc_roughness': rc_roughness,
            'N_corr': N_corr,
        })

    return galaxies


def build_model(X, y):
    """OLS regression."""
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    resid = y - y_hat
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - np.sum(resid**2) / ss_tot if ss_tot > 0 else 0
    rms = np.sqrt(np.mean(resid**2))
    return beta, y_hat, resid, R2, rms


def loo_rms(X, y):
    """Leave-one-out RMS via hat matrix."""
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    resid = y - y_hat
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    loo_r2 = 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)
    return np.sqrt(np.mean(loo_resid**2)), loo_r2


def adjusted_r2(R2, n, k):
    """Adjusted R² penalizing for extra parameters."""
    if n - k - 1 <= 0:
        return R2
    return 1 - (1 - R2) * (n - 1) / (n - k - 1)


print("=" * 70)
print("SESSION #483: THE SIXTH VARIABLE — CAN WE BREAK R² = 0.92?")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nSample: {n} galaxies")

# Build base 5-variable model (outer offset)
logV = np.log10([g['vflat'] for g in galaxies])
logL = np.log10([g['lum'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
y_outer = np.array([g['outer_offset'] for g in galaxies])

X5 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V])
beta5, yhat5, resid5, R2_5, rms5 = build_model(X5, y_outer)
loo5, loo_r2_5 = loo_rms(X5, y_outer)
adj_r2_5 = adjusted_r2(R2_5, n, 5)

print(f"\n5-variable baseline: R² = {R2_5:.4f}, adj R² = {adj_r2_5:.4f}")
print(f"  LOO RMS = {loo5:.4f}, LOO R² = {loo_r2_5:.4f}")

# Prepare candidate 6th variables
candidates = {}

# SB
sb_eff = np.array([g['sb_eff'] for g in galaxies])
log_sb = np.log10(np.clip(sb_eff, 1, None))
candidates['log_SB_eff'] = log_sb

# SB_disk
sb_disk = np.array([g['sb_disk'] for g in galaxies])
valid_sb_disk = sb_disk > 0
if valid_sb_disk.sum() > n * 0.5:
    log_sb_disk = np.log10(np.clip(sb_disk, 1, None))
    candidates['log_SB_disk'] = log_sb_disk

# N_corr
N_corr = np.array([g['N_corr'] for g in galaxies])
log_N = np.log10(np.clip(N_corr, 1e-5, None))
candidates['log_N_corr'] = log_N

# R_eff
r_eff = np.array([g['r_eff'] for g in galaxies])
log_reff = np.log10(np.clip(r_eff, 0.001, None))
candidates['log_R_eff'] = log_reff

# R_max
R_max = np.array([g['R_max'] for g in galaxies])
log_rmax = np.log10(np.clip(R_max, 0.01, None))
candidates['log_R_max'] = log_rmax

# Hubble type
T = np.array([g['hubble_type'] for g in galaxies], dtype=float)
candidates['T'] = T

# RC roughness
roughness = np.array([g['rc_roughness'] for g in galaxies])
candidates['roughness'] = roughness

# N_mond
n_mond = np.array([g['n_mond'] for g in galaxies], dtype=float)
candidates['N_mond'] = n_mond

# =====================================================================
# TEST 1: Single variable screening vs outer model residual
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: CANDIDATE SCREENING — CORRELATION WITH RESIDUAL")
print("=" * 60)

print(f"\n{'Variable':<15} {'r(X,resid)':<12} {'r(X,|resid|)':<14} {'r²(partial)':<12}")
print("-" * 53)

screening = {}
for name, var in candidates.items():
    finite = np.isfinite(var) & np.isfinite(resid5)
    if finite.sum() < 20:
        continue
    r_signed = np.corrcoef(var[finite], resid5[finite])[0, 1]
    r_abs = np.corrcoef(var[finite], np.abs(resid5[finite]))[0, 1]

    # Partial correlation: residual of var ~ X5, then correlate with resid5
    _, _, var_resid, _, _ = build_model(X5[finite], var[finite])
    r_partial = np.corrcoef(var_resid, resid5[finite])[0, 1]

    screening[name] = {
        'r_signed': r_signed, 'r_abs': r_abs,
        'r_partial': r_partial, 'var': var
    }
    print(f"{name:<15} {r_signed:+.4f}      {r_abs:+.4f}        {r_partial:+.4f}")

# Rank by |partial correlation|
ranked = sorted(screening.items(), key=lambda x: abs(x[1]['r_partial']), reverse=True)
print(f"\nTop candidates by |partial r|:")
for i, (name, info) in enumerate(ranked[:5]):
    print(f"  {i+1}. {name}: partial r = {info['r_partial']:+.4f}")

assert len(screening) >= 5, "Need at least 5 candidates"
print("\n✓ Test 1 passed: all candidates screened")

# =====================================================================
# TEST 2: 6-variable models — add each candidate
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: 6-VARIABLE MODELS")
print("=" * 60)

print(f"\n{'6th variable':<15} {'R²':<8} {'adj R²':<8} {'ΔR²':<8} {'RMS':<8} {'ΔRMS':<8}")
print("-" * 55)

six_var_results = {}
for name, info in screening.items():
    var = info['var']
    X6 = np.column_stack([X5, var])
    _, _, _, R2_6, rms6 = build_model(X6, y_outer)
    adj6 = adjusted_r2(R2_6, n, 6)
    dR2 = R2_6 - R2_5
    drms = rms6 - rms5
    six_var_results[name] = {'R2': R2_6, 'adj_R2': adj6, 'dR2': dR2, 'rms': rms6}
    print(f"{name:<15} {R2_6:.4f}  {adj6:.4f}  {dR2:+.4f}  {rms6:.4f}  {drms:+.4f}")

best_6th = max(six_var_results.items(), key=lambda x: x[1]['adj_R2'])
print(f"\nBest 6th variable (by adj R²): {best_6th[0]} → R² = {best_6th[1]['R2']:.4f} (ΔR² = {best_6th[1]['dR2']:+.4f})")

assert best_6th[1]['R2'] >= R2_5, "6th variable should not decrease R²"
print("✓ Test 2 passed: 6-variable models built")

# =====================================================================
# TEST 3: LOO validation of top candidates
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: LOO VALIDATION OF TOP CANDIDATES")
print("=" * 60)

# Take top 5 by ΔR²
top5 = sorted(six_var_results.items(), key=lambda x: x[1]['dR2'], reverse=True)[:5]

print(f"\n{'6th variable':<15} {'R²':<8} {'LOO RMS':<10} {'LOO R²':<10} {'ΔLOO R²':<10}")
print("-" * 53)

loo_results = {}
for name, info in top5:
    var = screening[name]['var']
    X6 = np.column_stack([X5, var])
    loo6, loo_r2_6 = loo_rms(X6, y_outer)
    d_loo_r2 = loo_r2_6 - loo_r2_5
    loo_results[name] = {'loo_rms': loo6, 'loo_r2': loo_r2_6, 'd_loo_r2': d_loo_r2}
    print(f"{name:<15} {info['R2']:.4f}  {loo6:.4f}    {loo_r2_6:.4f}    {d_loo_r2:+.4f}")

best_loo = max(loo_results.items(), key=lambda x: x[1]['loo_r2'])
print(f"\nBest by LOO R²: {best_loo[0]} → LOO R² = {best_loo[1]['loo_r2']:.4f} (ΔLOO R² = {best_loo[1]['d_loo_r2']:+.4f})")
print(f"Baseline LOO R²: {loo_r2_5:.4f}")

# Check if any candidate genuinely improves LOO
any_improve = any(v['d_loo_r2'] > 0.001 for v in loo_results.values())
print(f"\nAny candidate improves LOO R² by > 0.001? {'YES' if any_improve else 'NO'}")

assert len(loo_results) == 5, "Need LOO for top 5"
print("✓ Test 3 passed: LOO validation complete")

# =====================================================================
# TEST 4: Two-variable additions (7-variable models)
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: TWO-VARIABLE ADDITIONS (7-VARIABLE MODELS)")
print("=" * 60)

# Test all pairs of top candidates
top_names = [name for name, _ in top5[:4]]
pair_results = {}

print(f"\n{'Pair':<30} {'R²':<8} {'adj R²':<8} {'ΔR² vs 5':<10} {'LOO RMS':<10}")
print("-" * 66)

for i in range(len(top_names)):
    for j in range(i+1, len(top_names)):
        name_i = top_names[i]
        name_j = top_names[j]
        var_i = screening[name_i]['var']
        var_j = screening[name_j]['var']
        X7 = np.column_stack([X5, var_i, var_j])
        _, _, _, R2_7, rms7 = build_model(X7, y_outer)
        adj7 = adjusted_r2(R2_7, n, 7)
        loo7, loo_r2_7 = loo_rms(X7, y_outer)
        pair_name = f"{name_i} + {name_j}"
        dR2 = R2_7 - R2_5
        pair_results[pair_name] = {'R2': R2_7, 'adj_R2': adj7, 'loo_rms': loo7, 'loo_r2': loo_r2_7}
        print(f"{pair_name:<30} {R2_7:.4f}  {adj7:.4f}  {dR2:+.4f}    {loo7:.4f}")

best_pair = max(pair_results.items(), key=lambda x: x[1]['loo_r2'])
print(f"\nBest pair by LOO R²: {best_pair[0]} → LOO R² = {best_pair[1]['loo_r2']:.4f}")

assert len(pair_results) >= 3, "Need at least 3 pairs tested"
print("✓ Test 4 passed: 7-variable models tested")

# =====================================================================
# TEST 5: Non-linear terms — squared and interaction effects
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: NON-LINEAR TERMS")
print("=" * 60)

# Test quadratic and interaction terms of existing variables
nonlinear_terms = {
    'logV²': logV**2,
    'logL²': logL**2,
    'c_V²': c_V**2,
    'f_gas²': f_gas**2,
    'logV×logL': logV * logL,
    'logV×f_gas': logV * f_gas,
    'logL×c_V': logL * c_V,
    'logL×f_gas': logL * f_gas,
    'c_V×f_gas': c_V * f_gas,
}

# Also add the best external candidate as an interaction
if ranked:
    best_ext_name = ranked[0][0]
    best_ext_var = ranked[0][1]['var']
    nonlinear_terms[f'{best_ext_name}×logV'] = best_ext_var * logV
    nonlinear_terms[f'{best_ext_name}²'] = best_ext_var**2

print(f"\n{'Term':<20} {'R²':<8} {'adj R²':<8} {'ΔR²':<8} {'LOO RMS':<10}")
print("-" * 54)

nl_results = {}
for name, term in nonlinear_terms.items():
    if not np.all(np.isfinite(term)):
        continue
    X6nl = np.column_stack([X5, term])
    _, _, _, R2_nl, rms_nl = build_model(X6nl, y_outer)
    adj_nl = adjusted_r2(R2_nl, n, 6)
    loo_nl, loo_r2_nl = loo_rms(X6nl, y_outer)
    dR2 = R2_nl - R2_5
    nl_results[name] = {'R2': R2_nl, 'adj_R2': adj_nl, 'loo_rms': loo_nl, 'loo_r2': loo_r2_nl}
    print(f"{name:<20} {R2_nl:.4f}  {adj_nl:.4f}  {dR2:+.4f}  {loo_nl:.4f}")

best_nl = max(nl_results.items(), key=lambda x: x[1]['loo_r2'])
print(f"\nBest non-linear term by LOO R²: {best_nl[0]} → LOO R² = {best_nl[1]['loo_r2']:.4f}")

assert len(nl_results) >= 5, "Need at least 5 non-linear terms"
print("✓ Test 5 passed: non-linear terms tested")

# =====================================================================
# TEST 6: By galaxy type — which types benefit?
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: TYPE-DEPENDENT BENEFIT")
print("=" * 60)

# Use the best 6th variable from Test 2
best_var_name = best_6th[0]
best_var = screening[best_var_name]['var']

T_arr = np.array([g['hubble_type'] for g in galaxies])

type_groups = {
    'Early (T<4)': T_arr < 4,
    'Mid (4≤T<7)': (T_arr >= 4) & (T_arr < 7),
    'Late (T≥7)': T_arr >= 7,
}

print(f"\nAdding {best_var_name} as 6th variable:")
print(f"\n{'Type':<15} {'N':<5} {'R²_5var':<10} {'R²_6var':<10} {'ΔR²':<10} {'LOO_5':<10} {'LOO_6':<10}")
print("-" * 70)

for group_name, mask in type_groups.items():
    nm = mask.sum()
    if nm < 10:
        print(f"{group_name:<15} {nm:<5} (too few)")
        continue
    X5_g = X5[mask]
    X6_g = np.column_stack([X5_g, best_var[mask]])
    y_g = y_outer[mask]

    _, _, _, R2_5g, _ = build_model(X5_g, y_g)
    _, _, _, R2_6g, _ = build_model(X6_g, y_g)
    loo_5g, _ = loo_rms(X5_g, y_g)
    loo_6g, _ = loo_rms(X6_g, y_g)
    dR2 = R2_6g - R2_5g
    print(f"{group_name:<15} {nm:<5} {R2_5g:.4f}    {R2_6g:.4f}    {dR2:+.4f}    {loo_5g:.4f}    {loo_6g:.4f}")

print("✓ Test 6 passed: type-dependent analysis done")

# =====================================================================
# TEST 7: Forward stepwise selection
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: FORWARD STEPWISE SELECTION")
print("=" * 60)

# Start with intercept, add the variable that gives best LOO at each step
all_vars = {
    'logV': logV, 'logL': logL, 'c_V': c_V, 'f_gas': f_gas, 'logV×c_V': logV * c_V,
}
# Add external candidates
for name, info in screening.items():
    if name not in all_vars:
        all_vars[name] = info['var']

# Forward selection
selected = []
remaining = list(all_vars.keys())
X_step = np.ones((n, 1))  # intercept only

print(f"\n{'Step':<5} {'Added':<20} {'R²':<8} {'LOO R²':<10} {'ΔLOO R²':<10}")
print("-" * 53)

step_results = []
for step in range(min(10, len(remaining))):
    best_name = None
    best_loo_r2 = -999

    for name in remaining:
        var = all_vars[name]
        if not np.all(np.isfinite(var)):
            continue
        X_try = np.column_stack([X_step, var])
        try:
            _, loo_r2_try = loo_rms(X_try, y_outer)
        except:
            continue
        if loo_r2_try > best_loo_r2:
            best_loo_r2 = loo_r2_try
            best_name = name

    if best_name is None:
        break

    selected.append(best_name)
    remaining.remove(best_name)
    X_step = np.column_stack([X_step, all_vars[best_name]])

    _, _, _, R2_step, _ = build_model(X_step, y_outer)
    _, loo_r2_step = loo_rms(X_step, y_outer)

    prev_loo = step_results[-1]['loo_r2'] if step_results else 0
    d_loo = loo_r2_step - prev_loo

    step_results.append({'name': best_name, 'R2': R2_step, 'loo_r2': loo_r2_step})
    print(f"{step+1:<5} {best_name:<20} {R2_step:.4f}  {loo_r2_step:.4f}    {d_loo:+.4f}")

    # Stop if marginal improvement < 0.001
    if step > 0 and d_loo < 0.001:
        print(f"  → Stopped: marginal LOO R² improvement < 0.001")
        break

print(f"\nStepwise selected model ({len(selected)} variables): {' + '.join(selected)}")

# Verify it matches or exceeds the 5-var model
if len(step_results) >= 5:
    assert step_results[-1]['R2'] >= 0.85, f"Stepwise model should have R² ≥ 0.85"
print("✓ Test 7 passed: stepwise selection complete")

# =====================================================================
# TEST 8: SYNTHESIS — THE OPTIMAL MODEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS — THE OPTIMAL MODEL")
print("=" * 60)

# Compare all models
print("\n--- Model Comparison ---")
print(f"\n{'Model':<35} {'k':<4} {'R²':<8} {'adj R²':<8} {'LOO R²':<10} {'LOO RMS':<10}")
print("-" * 75)

# Baseline 5-var
print(f"{'5-var (logV,logL,c_V,f_gas,V×c_V)':<35} {5:<4} {R2_5:.4f}  {adj_r2_5:.4f}  {loo_r2_5:.4f}    {loo5:.4f}")

# Best single addition
if best_loo[0] in screening:
    var_b = screening[best_loo[0]]['var']
    X6b = np.column_stack([X5, var_b])
    _, _, _, R2_b, _ = build_model(X6b, y_outer)
    adj_b = adjusted_r2(R2_b, n, 6)
    loo_b, loo_r2_b = loo_rms(X6b, y_outer)
    print(f"{'5-var + ' + best_loo[0]:<35} {6:<4} {R2_b:.4f}  {adj_b:.4f}  {loo_r2_b:.4f}    {loo_b:.4f}")

# Best non-linear
best_nl_name = best_nl[0]
best_nl_info = best_nl[1]
print(f"{'5-var + ' + best_nl_name:<35} {6:<4} {best_nl_info['R2']:.4f}  {best_nl_info['adj_R2']:.4f}  {best_nl_info['loo_r2']:.4f}    {best_nl_info['loo_rms']:.4f}")

# Best pair
bp_name = best_pair[0]
bp_info = best_pair[1]
print(f"{'5-var + ' + bp_name:<35} {7:<4} {bp_info['R2']:.4f}  {bp_info['adj_R2']:.4f}  {bp_info['loo_r2']:.4f}    {bp_info['loo_rms']:.4f}")

# Stepwise
if step_results:
    last = step_results[-1]
    sw_name = f"Stepwise ({len(selected)} vars)"
    # Rebuild for adj R² and LOO
    X_sw = np.column_stack([np.ones(n)] + [all_vars[s] for s in selected])
    _, _, _, R2_sw, _ = build_model(X_sw, y_outer)
    adj_sw = adjusted_r2(R2_sw, n, len(selected))
    loo_sw, loo_r2_sw = loo_rms(X_sw, y_outer)
    print(f"{sw_name:<35} {len(selected):<4} {R2_sw:.4f}  {adj_sw:.4f}  {loo_r2_sw:.4f}    {loo_sw:.4f}")

# Summary analysis
print("\n--- Key Conclusions ---")

# Is the best 6th variable genuine?
if best_loo[1]['d_loo_r2'] > 0.005:
    print(f"\n✓ GENUINE IMPROVEMENT: {best_loo[0]} improves LOO R² by {best_loo[1]['d_loo_r2']:+.4f}")
    print(f"  New LOO R² = {best_loo[1]['loo_r2']:.4f} vs baseline {loo_r2_5:.4f}")
elif best_loo[1]['d_loo_r2'] > 0:
    print(f"\n~ MARGINAL IMPROVEMENT: {best_loo[0]} improves LOO R² by only {best_loo[1]['d_loo_r2']:+.4f}")
    print(f"  Not convincingly better than 5-variable model")
else:
    print(f"\n✗ NO IMPROVEMENT: Best candidate ({best_loo[0]}) does not improve LOO R²")
    print(f"  5-variable model is already optimal for this sample")

# NN autocorrelation test on best 6-var model
if best_loo[0] in screening:
    var_best = screening[best_loo[0]]['var']
    X6_best = np.column_stack([X5, var_best])
    _, _, resid_6best, _, _ = build_model(X6_best, y_outer)

    # NN autocorrelation
    from scipy.spatial.distance import cdist
    feat = np.column_stack([logV, logL, c_V, f_gas])
    feat_std = (feat - feat.mean(axis=0)) / np.clip(feat.std(axis=0), 1e-10, None)
    dist_matrix = cdist(feat_std, feat_std, 'euclidean')
    np.fill_diagonal(dist_matrix, np.inf)
    nn_idx = np.argmin(dist_matrix, axis=1)

    r_nn_5 = np.corrcoef(resid5, resid5[nn_idx])[0, 1]
    r_nn_6 = np.corrcoef(resid_6best, resid_6best[nn_idx])[0, 1]

    print(f"\nNN autocorrelation:")
    print(f"  5-var model: r = {r_nn_5:+.3f}")
    print(f"  6-var model ({best_loo[0]}): r = {r_nn_6:+.3f}")

    if r_nn_6 < r_nn_5 - 0.05:
        print(f"  → 6th variable REDUCES autocorrelation by {r_nn_5 - r_nn_6:.3f}")
    else:
        print(f"  → 6th variable does NOT reduce autocorrelation")

# Check if the 5-var model is already at the information ceiling
noise_est = 0.025  # from Session 482
irreducible = np.sqrt(max(rms5**2 - noise_est**2, 0))
print(f"\nInformation ceiling analysis:")
print(f"  5-var RMS: {rms5:.4f} dex")
print(f"  Estimated noise: {noise_est:.4f} dex")
print(f"  Irreducible: {irreducible:.4f} dex")
print(f"  Ceiling R² (if all reducible variance captured): "
      f"{1 - noise_est**2 / np.var(y_outer):.4f}")

print("\n✓ Test 8 passed: synthesis complete")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #483 SUMMARY")
print("=" * 70)
print(f"Sample: {n} galaxies")
print(f"5-var baseline: R² = {R2_5:.4f}, LOO R² = {loo_r2_5:.4f}")
print(f"Best 6th variable (in-sample): {best_6th[0]} (ΔR² = {best_6th[1]['dR2']:+.4f})")
print(f"Best 6th variable (LOO): {best_loo[0]} (ΔLOO R² = {best_loo[1]['d_loo_r2']:+.4f})")
print(f"Best non-linear term: {best_nl_name} (LOO R² = {best_nl_info['loo_r2']:.4f})")
print(f"Best pair: {bp_name} (LOO R² = {bp_info['loo_r2']:.4f})")
print(f"Stepwise: {len(selected)} variables, LOO R² = {step_results[-1]['loo_r2']:.4f}")
print(f"\nAll 8 tests passed ✓")
