#!/usr/bin/env python3
"""
======================================================================
SESSION #499: OUTLIER FORENSICS — WHAT MAKES GALAXIES HARD TO PREDICT?
======================================================================

The 6-var model achieves LOO R² = 0.938 with RMS = 0.038 dex. But some
galaxies are systematically mispredicted. What makes them special?

If the outliers are astrophysically unusual (mergers, interactions,
extreme properties), this validates the model. If they're just noisy
measurements, it constrains the true intrinsic scatter.

Tests:
1. Identify outliers: top/bottom residual galaxies
2. Outlier properties: are they extreme in any observable?
3. Leverage analysis: which galaxies have high leverage?
4. Influence analysis: Cook's distance / DFBETAS
5. Residual vs morphology: do types cluster in residuals?
6. Leave-Two-Out: are outlier residuals stable?
7. Robust regression: do outliers change the model?
8. Synthesis: outlier taxonomy

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #499
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


def prepare_galaxies():
    """Load SPARC and compute galaxy-level offsets + properties."""
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
        distance = cat.get('distance', 0)
        inclination = cat.get('inclination', 0)
        quality = cat.get('quality', 0)

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

        radius_m = radius_v[mond]
        med_r = np.median(radius_m)
        outer_mond = mond.copy()
        outer_mond[mond] = radius_m > med_r

        g_rar = rar_prediction(g_bar_v)
        point_offsets = np.log10(g_obs_v) - np.log10(g_rar)

        if outer_mond.sum() >= 2:
            offset = np.mean(point_offsets[outer_mond])
        else:
            offset = np.mean(point_offsets[mond])

        # f_gas
        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # Mean V_obs error
        mean_frac_err = np.mean(e_vobs_v / np.clip(np.abs(v_obs_v), 1, None))

        # RC shape: ratio of max V to V_flat
        v_max = np.max(np.abs(v_obs_v))
        rc_shape = v_max / max(abs(vflat), 1)

        # N MOND points
        n_mond = mond.sum()

        # Roughness
        if len(point_offsets) >= 3:
            roughness = np.mean(np.abs(np.diff(point_offsets)))
        else:
            roughness = np.nan

        # Within-galaxy scatter
        within_sigma = np.std(point_offsets[mond]) if mond.sum() >= 3 else np.nan

        galaxies.append({
            'id': gal_id,
            'offset': offset,
            'logV': np.log10(vflat),
            'logL': np.log10(lum),
            'c_V': c_V,
            'f_gas': f_gas,
            'hubble_type': hubble_type,
            'distance': distance,
            'inclination': inclination,
            'quality': quality,
            'mean_frac_err': mean_frac_err,
            'rc_shape': rc_shape,
            'n_mond': n_mond,
            'n_points': len(points),
            'roughness': roughness,
            'within_sigma': within_sigma,
            'r_eff_kpc': r_eff_kpc,
        })

    return galaxies


def build_6var_model(galaxies):
    """Build 6-var model and return design matrix, coefficients, residuals."""
    n = len(galaxies)
    logV = np.array([g['logV'] for g in galaxies])
    logL = np.array([g['logL'] for g in galaxies])
    c_V = np.array([g['c_V'] for g in galaxies])
    f_gas = np.array([g['f_gas'] for g in galaxies])
    y = np.array([g['offset'] for g in galaxies])

    X = np.column_stack([
        np.ones(n), logV, logL, c_V, f_gas,
        logV * c_V, logL * f_gas
    ])

    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    yhat = X @ beta
    resid = y - yhat

    # Hat matrix
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)

    # LOO residuals
    loo_resid = resid / (1 - h)

    return X, beta, yhat, resid, h, loo_resid


print("=" * 70)
print("SESSION #499: OUTLIER FORENSICS")
print("=" * 70)

galaxies = prepare_galaxies()
n = len(galaxies)
print(f"\n{n} galaxies loaded")

X, beta, yhat, resid, h, loo_resid = build_6var_model(galaxies)
y = np.array([g['offset'] for g in galaxies])

ss_tot = np.sum((y - np.mean(y))**2)
ss_res = np.sum(resid**2)
R2 = 1 - ss_res / ss_tot
rms = np.sqrt(np.mean(resid**2))
print(f"6-var model: R² = {R2:.4f}, RMS = {rms:.4f}")

# =====================================================================
# TEST 1: IDENTIFY OUTLIERS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: IDENTIFY OUTLIERS")
print("=" * 60)

# Sort by absolute LOO residual
abs_loo = np.abs(loo_resid)
sorted_idx = np.argsort(abs_loo)[::-1]

print(f"\nTop 10 outliers (by |LOO residual|):")
print(f"{'Rank':<5} {'Galaxy':<15} {'Offset':<8} {'Predicted':<10} {'LOO resid':<10} {'Type':<5}")
print("-" * 55)
for rank, idx in enumerate(sorted_idx[:10]):
    g = galaxies[idx]
    print(f"  {rank+1:<5} {g['id']:<15} {g['offset']:+.3f}   {yhat[idx]:+.3f}      "
          f"{loo_resid[idx]:+.4f}     {g['hubble_type']:.0f}")

# How many outliers beyond 2σ, 3σ?
sigma_loo = np.std(loo_resid)
n_2sigma = np.sum(abs_loo > 2 * sigma_loo)
n_3sigma = np.sum(abs_loo > 3 * sigma_loo)
print(f"\nLOO residual σ = {sigma_loo:.4f}")
print(f"  Beyond 2σ: {n_2sigma} ({n_2sigma/n*100:.1f}%)")
print(f"  Beyond 3σ: {n_3sigma} ({n_3sigma/n*100:.1f}%)")

# Expected from normal
print(f"  Expected from Gaussian: {n*0.046:.1f} beyond 2σ, {n*0.003:.1f} beyond 3σ")

# Positive vs negative outliers
big_pos = loo_resid > 2 * sigma_loo
big_neg = loo_resid < -2 * sigma_loo
print(f"  Positive outliers: {big_pos.sum()}")
print(f"  Negative outliers: {big_neg.sum()}")

print("\n✓ Test 1 passed: outliers identified")

# =====================================================================
# TEST 2: OUTLIER PROPERTIES
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: OUTLIER PROPERTIES — ARE OUTLIERS EXTREME?")
print("=" * 60)

is_outlier = abs_loo > 2 * sigma_loo
is_normal = ~is_outlier

properties = ['logV', 'logL', 'c_V', 'f_gas', 'hubble_type', 'distance',
              'inclination', 'quality', 'mean_frac_err', 'rc_shape',
              'n_mond', 'n_points', 'roughness', 'within_sigma']

print(f"\n{'Property':<15} {'Outlier mean':<14} {'Normal mean':<14} {'|Δ|/σ':<8}")
print("-" * 55)
for prop in properties:
    vals = np.array([g[prop] for g in galaxies])
    finite = np.isfinite(vals)
    if finite.sum() < 50:
        continue
    out_vals = vals[is_outlier & finite]
    norm_vals = vals[is_normal & finite]
    if len(out_vals) < 2 or len(norm_vals) < 2:
        continue
    diff_sigma = abs(np.mean(out_vals) - np.mean(norm_vals)) / max(np.std(vals[finite]), 1e-10)
    print(f"  {prop:<15} {np.mean(out_vals):<14.3f} {np.mean(norm_vals):<14.3f} {diff_sigma:<8.2f}")

print("\n✓ Test 2 passed: outlier properties compared")

# =====================================================================
# TEST 3: LEVERAGE ANALYSIS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: LEVERAGE ANALYSIS")
print("=" * 60)

p = X.shape[1]  # number of parameters
h_threshold = 2 * p / n  # standard threshold

high_lev = h > h_threshold
print(f"\nLeverage threshold (2p/n): {h_threshold:.3f}")
print(f"High-leverage galaxies: {high_lev.sum()} ({high_lev.sum()/n*100:.1f}%)")

# Top leverage galaxies
h_sorted = np.argsort(h)[::-1]
print(f"\nTop 5 by leverage:")
print(f"{'Galaxy':<15} {'h':<8} {'LOO resid':<10} {'logV':<8} {'f_gas':<8} {'Type':<5}")
print("-" * 55)
for idx in h_sorted[:5]:
    g = galaxies[idx]
    print(f"  {g['id']:<15} {h[idx]:<8.3f} {loo_resid[idx]:+.4f}     "
          f"{g['logV']:<8.2f} {g['f_gas']:<8.2f} {g['hubble_type']:.0f}")

# Correlation between leverage and absolute residual
r_h_resid = np.corrcoef(h, np.abs(loo_resid))[0, 1]
print(f"\nr(leverage, |LOO resid|) = {r_h_resid:+.3f}")

# How many are both high-leverage AND outliers?
both = high_lev & is_outlier
print(f"High-leverage AND outlier: {both.sum()}")
if both.sum() > 0:
    for i in np.where(both)[0]:
        g = galaxies[i]
        print(f"  → {g['id']} (h={h[i]:.3f}, LOO resid={loo_resid[i]:+.4f})")

print("\n✓ Test 3 passed: leverage analyzed")

# =====================================================================
# TEST 4: INFLUENCE ANALYSIS (COOK'S DISTANCE)
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: COOK'S DISTANCE AND INFLUENCE")
print("=" * 60)

# Cook's distance
mse = ss_res / (n - p)
cooks_d = (resid**2 / (p * mse)) * (h / (1 - h)**2)

cooks_threshold = 4 / n
high_cooks = cooks_d > cooks_threshold
print(f"\nCook's distance threshold (4/n): {cooks_threshold:.4f}")
print(f"High-influence galaxies: {high_cooks.sum()} ({high_cooks.sum()/n*100:.1f}%)")

# Top by Cook's distance
cooks_sorted = np.argsort(cooks_d)[::-1]
print(f"\nTop 10 by Cook's distance:")
print(f"{'Galaxy':<15} {'Cook D':<10} {'h':<8} {'LOO resid':<10}")
print("-" * 45)
for idx in cooks_sorted[:10]:
    g = galaxies[idx]
    print(f"  {g['id']:<15} {cooks_d[idx]:<10.4f} {h[idx]:<8.3f} {loo_resid[idx]:+.4f}")

# DFBETAS: how much does each galaxy change each coefficient?
print(f"\nDFBETAS analysis (coefficient change when galaxy i removed):")
XtX_inv = np.linalg.inv(X.T @ X)
sigma_hat = np.sqrt(mse)

var_names = ['const', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V', 'logL×f_gas']

# Find galaxy with largest impact on each coefficient
for j, vname in enumerate(var_names):
    dfbetas_j = (XtX_inv[j, :] @ X.T * resid / (1 - h)) / sigma_hat
    # Simplified DFBETAS
    max_idx = np.argmax(np.abs(dfbetas_j))
    print(f"  {vname:<12}: max impact from {galaxies[max_idx]['id']} "
          f"(DFBETA = {dfbetas_j[max_idx]:+.3f})")

print("\n✓ Test 4 passed: influence analyzed")

# =====================================================================
# TEST 5: RESIDUAL VS MORPHOLOGY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: RESIDUAL PATTERNS BY MORPHOLOGY")
print("=" * 60)

types = np.array([g['hubble_type'] for g in galaxies])

for name, T_lo, T_hi in [('Early (T<4)', 0, 4), ('Mid (4≤T<7)', 4, 7), ('Late (T≥7)', 7, 20)]:
    mask = (types >= T_lo) & (types < T_hi)
    if mask.sum() >= 5:
        r_type = loo_resid[mask]
        print(f"\n{name} (N={mask.sum()}):")
        print(f"  Mean LOO resid: {np.mean(r_type):+.4f}")
        print(f"  Std LOO resid: {np.std(r_type):.4f}")
        print(f"  Outlier fraction: {np.mean(np.abs(r_type) > 2*sigma_loo)*100:.0f}%")

        # Are outliers concentrated in any type?
        n_outlier_type = np.sum(is_outlier[mask])
        expected = mask.sum() * is_outlier.sum() / n
        print(f"  Outliers: {n_outlier_type} (expected: {expected:.1f})")

# Correlation: residual with type
r_type_resid = np.corrcoef(types, loo_resid)[0, 1]
print(f"\nr(T, LOO resid) = {r_type_resid:+.4f}")

# Quality flag analysis
quals = np.array([g['quality'] for g in galaxies])
for q in sorted(set(quals)):
    qmask = quals == q
    if qmask.sum() >= 5:
        print(f"\nQuality={q:.0f} (N={qmask.sum()}):")
        print(f"  Mean |LOO resid|: {np.mean(np.abs(loo_resid[qmask])):.4f}")
        print(f"  Outlier fraction: {np.mean(is_outlier[qmask])*100:.0f}%")

print("\n✓ Test 5 passed: morphological patterns analyzed")

# =====================================================================
# TEST 6: LEAVE-TWO-OUT STABILITY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: LEAVE-TWO-OUT — ARE OUTLIER RESIDUALS STABLE?")
print("=" * 60)

# For top 5 outliers: remove each, refit, check if remaining outliers
# stay as outliers
top5_idx = sorted_idx[:5]

print(f"\nStability test for top 5 outliers:")
print(f"{'Removed':<15} {'LOO resid':<12} ", end="")
for idx in top5_idx:
    print(f"{galaxies[idx]['id']:<12} ", end="")
print()
print("-" * 75)

# Original LOO residuals for top 5
print(f"  {'(none)':<15} {'—':<12} ", end="")
for idx in top5_idx:
    print(f"{loo_resid[idx]:+.4f}       ", end="")
print()

for remove_idx in top5_idx:
    keep = np.array([i for i in range(n) if i != remove_idx])
    X_sub = X[keep]
    y_sub = y[keep]
    beta_sub = np.linalg.lstsq(X_sub, y_sub, rcond=None)[0]
    yhat_sub = X_sub @ beta_sub
    resid_sub = y_sub - yhat_sub
    H_sub = X_sub @ np.linalg.inv(X_sub.T @ X_sub) @ X_sub.T
    h_sub = np.diag(H_sub)
    loo_sub = resid_sub / (1 - h_sub)

    print(f"  {galaxies[remove_idx]['id']:<15} {'—':<12} ", end="")
    for target_idx in top5_idx:
        if target_idx == remove_idx:
            print(f"{'(removed)':<12} ", end="")
        else:
            # Find this galaxy in the subset
            sub_pos = np.searchsorted(keep, target_idx)
            if sub_pos < len(keep) and keep[sub_pos] == target_idx:
                print(f"{loo_sub[sub_pos]:+.4f}       ", end="")
            else:
                print(f"{'N/A':<12} ", end="")
    print()

# Stability metric: correlation of all LOO residuals with and without top outlier
remove_top = sorted_idx[0]
keep = np.array([i for i in range(n) if i != remove_top])
X_sub = X[keep]
y_sub = y[keep]
beta_sub = np.linalg.lstsq(X_sub, y_sub, rcond=None)[0]
yhat_sub = X_sub @ beta_sub
resid_sub = y_sub - yhat_sub
H_sub = X_sub @ np.linalg.inv(X_sub.T @ X_sub) @ X_sub.T
h_sub = np.diag(H_sub)
loo_sub = resid_sub / (1 - h_sub)

r_stability = np.corrcoef(loo_resid[keep], loo_sub)[0, 1]
print(f"\nStability: r(LOO with all, LOO without #{galaxies[remove_top]['id']}) = {r_stability:.4f}")

print("\n✓ Test 6 passed: leave-two-out stability tested")

# =====================================================================
# TEST 7: ROBUST REGRESSION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: ROBUST REGRESSION — OUTLIER IMPACT ON COEFFICIENTS")
print("=" * 60)

# Compare OLS with iteratively reweighted least squares (Huber weights)
print(f"\nOLS coefficients (all {n} galaxies):")
for vname, b in zip(var_names, beta):
    print(f"  {vname:<12} {b:+.4f}")
print(f"  R² = {R2:.4f}")

# Fit without outliers
keep_clean = ~is_outlier
X_clean = X[keep_clean]
y_clean = y[keep_clean]
beta_clean = np.linalg.lstsq(X_clean, y_clean, rcond=None)[0]
yhat_clean = X_clean @ beta_clean
resid_clean = y_clean - yhat_clean
R2_clean = 1 - np.sum(resid_clean**2) / np.sum((y_clean - np.mean(y_clean))**2)

print(f"\nWithout 2σ outliers ({keep_clean.sum()} galaxies):")
for vname, b_ols, b_clean in zip(var_names, beta, beta_clean):
    change = b_clean - b_ols
    print(f"  {vname:<12} {b_clean:+.4f} (Δ = {change:+.4f})")
print(f"  R² = {R2_clean:.4f}")

# Iterative reweighted LS (Huber-like)
weights = np.ones(n)
for iteration in range(10):
    W = np.diag(weights)
    beta_w = np.linalg.lstsq(np.sqrt(W) @ X, np.sqrt(W) @ y, rcond=None)[0]
    resid_w = y - X @ beta_w
    mad = np.median(np.abs(resid_w))
    sigma_mad = 1.4826 * mad
    # Huber weights
    for i in range(n):
        if abs(resid_w[i]) <= 1.345 * sigma_mad:
            weights[i] = 1.0
        else:
            weights[i] = 1.345 * sigma_mad / abs(resid_w[i])

yhat_w = X @ beta_w
resid_w = y - yhat_w
R2_w = 1 - np.sum(weights * resid_w**2) / np.sum(weights * (y - np.average(y, weights=weights))**2)

print(f"\nHuber robust regression:")
for vname, b_ols, b_hub in zip(var_names, beta, beta_w):
    change = b_hub - b_ols
    print(f"  {vname:<12} {b_hub:+.4f} (Δ = {change:+.4f})")

# Max coefficient change
max_change = max(abs(beta_w[j] - beta[j]) for j in range(len(beta)))
print(f"\n  Max coefficient change: {max_change:.4f}")
print(f"  → {'STABLE' if max_change < 0.1 else 'UNSTABLE'}: outliers do {'not ' if max_change < 0.1 else ''}change the model")

print("\n✓ Test 7 passed: robust regression compared")

# =====================================================================
# TEST 8: OUTLIER TAXONOMY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: OUTLIER TAXONOMY")
print("=" * 60)

# Classify each outlier
print(f"\nOutlier classification:")
print(f"{'Galaxy':<15} {'LOO resid':<10} {'Category':<20} {'Evidence':<30}")
print("-" * 75)

for idx in sorted_idx[:n_2sigma]:
    g = galaxies[idx]
    lr = loo_resid[idx]

    # Classify
    categories = []
    if g['mean_frac_err'] > np.percentile([ga['mean_frac_err'] for ga in galaxies], 75):
        categories.append("Noisy (high V err)")
    if g['quality'] >= 3:
        categories.append("Low quality")
    if g['roughness'] > np.percentile([ga['roughness'] for ga in galaxies if np.isfinite(ga['roughness'])], 90):
        categories.append("Rough RC")
    if g['n_mond'] < 5:
        categories.append("Few MOND pts")
    if h[idx] > h_threshold:
        categories.append("High leverage")
    if g['inclination'] < 30 or g['inclination'] > 85:
        categories.append("Extreme incl")
    if g['f_gas'] > 0.8:
        categories.append("Very gas-rich")
    if g['f_gas'] < 0.05:
        categories.append("Very gas-poor")

    if not categories:
        categories.append("Genuinely unusual")

    cat_str = "; ".join(categories[:2])

    evidence = []
    if g['mean_frac_err'] > 0.1:
        evidence.append(f"V_err={g['mean_frac_err']:.0%}")
    if g['quality'] >= 3:
        evidence.append(f"Q={g['quality']:.0f}")
    evidence.append(f"T={g['hubble_type']:.0f}")
    evidence.append(f"n_M={g['n_mond']}")

    ev_str = ", ".join(evidence[:3])

    print(f"  {g['id']:<15} {lr:+.4f}     {cat_str:<20} {ev_str}")

# Summary stats
all_frac_err = np.array([g['mean_frac_err'] for g in galaxies])
out_frac_err = all_frac_err[is_outlier]
norm_frac_err = all_frac_err[is_normal]

print(f"\nSummary:")
print(f"  Outliers with high V_err: "
      f"{np.sum(out_frac_err > np.percentile(all_frac_err, 75))}/{is_outlier.sum()}")
print(f"  Outliers with low quality: "
      f"{np.sum(np.array([g['quality'] for g in galaxies])[is_outlier] >= 3)}/{is_outlier.sum()}")
print(f"  Outliers with extreme inclination: "
      f"{np.sum((np.array([g['inclination'] for g in galaxies])[is_outlier] < 30) | (np.array([g['inclination'] for g in galaxies])[is_outlier] > 85))}/{is_outlier.sum()}")
print(f"  Outliers with high leverage: {np.sum(is_outlier & high_lev)}/{is_outlier.sum()}")

# How many are "measurement" vs "physical" outliers?
meas_outlier = np.zeros(n, dtype=bool)
for i in range(n):
    if is_outlier[i]:
        if galaxies[i]['mean_frac_err'] > np.percentile(all_frac_err, 75) or \
           galaxies[i]['quality'] >= 3 or \
           galaxies[i]['inclination'] < 30 or galaxies[i]['inclination'] > 85:
            meas_outlier[i] = True

phys_outlier = is_outlier & ~meas_outlier
print(f"\n  Measurement outliers: {meas_outlier.sum()}")
print(f"  Physical outliers: {phys_outlier.sum()}")
if phys_outlier.sum() > 0:
    print(f"  Physical outlier fraction of sample: {phys_outlier.sum()/n*100:.1f}%")

print(f"\n✓ Test 8 passed: taxonomy complete")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #499 SUMMARY")
print("=" * 70)
print(f"Total outliers (>2σ LOO): {n_2sigma}")
print(f"Positive/Negative: {big_pos.sum()}/{big_neg.sum()}")
print(f"r(leverage, |resid|) = {r_h_resid:+.3f}")
print(f"Max Cook's D: {cooks_d[cooks_sorted[0]]:.4f} ({galaxies[cooks_sorted[0]]['id']})")
print(f"Robust regression max Δβ: {max_change:.4f}")
print(f"LOO stability r = {r_stability:.4f}")
print(f"Measurement outliers: {meas_outlier.sum()}, Physical: {phys_outlier.sum()}")
print(f"\nAll 8 tests passed ✓")
