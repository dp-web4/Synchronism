#!/usr/bin/env python3
"""
======================================================================
SESSION #509: RESIDUAL FORENSICS — WHAT'S LEFT IN THE 0.038 DEX?
======================================================================

The 6-var model leaves RMS = 0.038 dex. Session #491 showed noise = 28%
of total variance. The remaining ~72% of residual should be:
- M/L errors (estimated RMS ≈ 0.03 dex from session #491)
- Distance errors
- Inclination errors
- MOND interpolation function imperfection

This session forensically examines the residual to identify what
information it still contains, focusing on observational properties
NOT in the model.

Tests:
1. Residual vs all available catalog properties
2. Residual spatial distribution: distance, environment
3. Residual vs data quality indicators
4. Residual vs inclination (systematic bias test)
5. Nearest-neighbor residual correlation (spatial clustering)
6. Residual normality and outlier structure
7. Residual vs MOND regime (interpolation function test)
8. The irreducible floor: what fraction is truly random?

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #509
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


def build_model(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    yhat = X @ beta
    resid = y - yhat
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot
    rms = np.sqrt(np.mean(resid**2))
    return beta, yhat, resid, R2, rms


def loo_r2(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    return 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)


def prepare_galaxies():
    """Load SPARC with full catalog properties."""
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
        e_vobs = np.array([pt.get('e_vobs', 5) for pt in points])

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
            mean_g_bar = np.mean(g_bar_v[outer_mond])
            n_mond_pts = outer_mond.sum()
        else:
            offset = np.mean(point_offsets[mond])
            mean_g_bar = np.mean(g_bar_v[mond])
            n_mond_pts = mond.sum()

        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # Number of RC points
        n_points = valid.sum()

        # Mean velocity error fraction
        mean_e_frac = np.mean(e_vobs_v / np.abs(v_obs_v + 1e-10))

        # Radial extent
        r_max_kpc = radius_v.max()

        # Has bulge?
        has_bulge = 1 if np.any(v_bul != 0) and np.max(np.abs(v_bul)) > 0.1 else 0

        # RC roughness: RMS of v_obs variations around smooth trend
        if n_points >= 5:
            from numpy.polynomial import polynomial as P
            coeffs = P.polyfit(radius_v, v_obs_v, min(3, n_points-1))
            smooth = P.polyval(radius_v, coeffs)
            roughness = np.sqrt(np.mean((v_obs_v - smooth)**2)) / vflat
        else:
            roughness = np.nan

        # Within-galaxy scatter of point offsets (outer MOND)
        if n_mond_pts >= 3:
            if outer_mond.sum() >= 3:
                within_scatter = np.std(point_offsets[outer_mond])
            else:
                within_scatter = np.std(point_offsets[mond])
        else:
            within_scatter = np.nan

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
            'log_g_ratio': np.log10(mean_g_bar / a0_mond),
            'n_points': n_points,
            'n_mond_pts': n_mond_pts,
            'mean_e_frac': mean_e_frac,
            'r_max_kpc': r_max_kpc,
            'has_bulge': has_bulge,
            'sb_eff': sb_eff,
            'sb_disk': sb_disk if sb_disk > 0 else np.nan,
            'roughness': roughness,
            'within_scatter': within_scatter,
            'vflat': vflat,
            'lum': lum,
        })

    return galaxies


print("=" * 70)
print("SESSION #509: RESIDUAL FORENSICS — WHAT'S LEFT IN THE 0.038 DEX?")
print("=" * 70)

galaxies = prepare_galaxies()
n = len(galaxies)
print(f"\n{n} galaxies loaded")

# Extract arrays
offset = np.array([g['offset'] for g in galaxies])
logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])

# Build reference model and get residuals
X6 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])
beta6, yhat6, resid, R2_6, rms_6 = build_model(X6, offset)

print(f"\n6-var model: R² = {R2_6:.4f}, RMS = {rms_6:.4f}")

# =====================================================================
# TEST 1: RESIDUAL VS ALL CATALOG PROPERTIES
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: RESIDUAL vs ALL AVAILABLE PROPERTIES")
print("=" * 60)

# Properties NOT in the model
properties = {
    'hubble_type': np.array([g['hubble_type'] for g in galaxies]),
    'distance': np.array([g['distance'] for g in galaxies]),
    'inclination': np.array([g['inclination'] for g in galaxies]),
    'quality': np.array([g['quality'] for g in galaxies]),
    'n_points': np.array([g['n_points'] for g in galaxies]),
    'n_mond_pts': np.array([g['n_mond_pts'] for g in galaxies]),
    'mean_e_frac': np.array([g['mean_e_frac'] for g in galaxies]),
    'r_max_kpc': np.array([g['r_max_kpc'] for g in galaxies]),
    'has_bulge': np.array([g['has_bulge'] for g in galaxies]),
    'sb_eff': np.array([g['sb_eff'] for g in galaxies]),
    'within_scatter': np.array([g['within_scatter'] for g in galaxies]),
    'log_g_ratio': np.array([g['log_g_ratio'] for g in galaxies]),
}

# Also include derived properties
properties['log_distance'] = np.log10(np.maximum(properties['distance'], 0.1))
properties['sin_incl'] = np.sin(np.radians(properties['inclination']))
properties['log_rmax'] = np.log10(np.maximum(properties['r_max_kpc'], 0.1))

print(f"\n{'Property':<20} {'r(resid)':>10} {'p-value':>10} {'Note'}")
print("-" * 55)

significant = []
for name, vals in sorted(properties.items()):
    valid_mask = np.isfinite(vals)
    if valid_mask.sum() < 20:
        continue
    r = np.corrcoef(resid[valid_mask], vals[valid_mask])[0, 1]
    # Approximate p-value for correlation
    n_valid = valid_mask.sum()
    t_stat = r * np.sqrt(n_valid - 2) / np.sqrt(1 - r**2 + 1e-10)
    # Two-sided p from t-distribution approximation
    from scipy import stats as sp_stats
    p_val = 2 * sp_stats.t.sf(abs(t_stat), n_valid - 2)
    sig = "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
    print(f"  {name:<20} {r:>+10.4f} {p_val:>10.4f} {sig}")
    if p_val < 0.05:
        significant.append((name, r, p_val))

print(f"\n  Significant at p<0.05: {len(significant)} / {len(properties)}")
for name, r, p in significant:
    print(f"    {name}: r = {r:+.4f} (p = {p:.4f})")

print("\n✓ Test 1 passed: residual correlations computed")

# =====================================================================
# TEST 2: DISTANCE AND ENVIRONMENT
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: DISTANCE AND ENVIRONMENTAL EFFECTS")
print("=" * 60)

distance = properties['distance']
log_dist = properties['log_distance']

r_dist = np.corrcoef(resid, distance)[0, 1]
r_logdist = np.corrcoef(resid, log_dist)[0, 1]

print(f"\n  r(resid, distance) = {r_dist:+.4f}")
print(f"  r(resid, log distance) = {r_logdist:+.4f}")

# Split near vs far
med_dist = np.median(distance)
near = distance <= med_dist
far = distance > med_dist

print(f"\n  Median distance: {med_dist:.1f} Mpc")
print(f"  Near ({near.sum()}) mean resid: {np.mean(resid[near]):+.4f} ± {np.std(resid[near])/np.sqrt(near.sum()):.4f}")
print(f"  Far ({far.sum()}) mean resid: {np.mean(resid[far]):+.4f} ± {np.std(resid[far])/np.sqrt(far.sum()):.4f}")

# t-test
t_dist = (np.mean(resid[near]) - np.mean(resid[far])) / np.sqrt(
    np.var(resid[near])/near.sum() + np.var(resid[far])/far.sum())
p_dist = 2 * sp_stats.t.sf(abs(t_dist), n - 2)
print(f"  t = {t_dist:+.3f}, p = {p_dist:.4f}")

# Also check distance×V interaction (nearby massive vs distant massive)
# This could reveal resolution effects
near_massive = near & (logV > np.median(logV))
far_massive = far & (logV > np.median(logV))
if near_massive.sum() >= 5 and far_massive.sum() >= 5:
    print(f"\n  Massive galaxies (logV > median):")
    print(f"  Near ({near_massive.sum()}) mean resid: {np.mean(resid[near_massive]):+.4f}")
    print(f"  Far ({far_massive.sum()}) mean resid: {np.mean(resid[far_massive]):+.4f}")

print("\n✓ Test 2 passed: distance effects tested")

# =====================================================================
# TEST 3: DATA QUALITY INDICATORS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: DATA QUALITY INDICATORS")
print("=" * 60)

quality = properties['quality']
mean_e_frac = properties['mean_e_frac']
n_points = properties['n_points']

# Quality flags
for q in sorted(np.unique(quality)):
    mask = quality == q
    if mask.sum() >= 3:
        print(f"  Quality={int(q)}: N={mask.sum():>3}, mean resid={np.mean(resid[mask]):+.5f}, std={np.std(resid[mask]):.4f}")

# Error fraction
r_efrac = np.corrcoef(resid, mean_e_frac)[0, 1]
print(f"\n  r(resid, mean error fraction) = {r_efrac:+.4f}")

# Number of points
r_npts = np.corrcoef(resid, n_points)[0, 1]
print(f"  r(resid, N_points) = {r_npts:+.4f}")

# Within-galaxy scatter
ws = properties['within_scatter']
ws_valid = np.isfinite(ws)
if ws_valid.sum() > 20:
    r_ws = np.corrcoef(resid[ws_valid], ws[ws_valid])[0, 1]
    print(f"  r(resid, within-galaxy scatter) = {r_ws:+.4f}")

    # Do noisy galaxies have larger |residuals|?
    r_abs_ws = np.corrcoef(np.abs(resid[ws_valid]), ws[ws_valid])[0, 1]
    print(f"  r(|resid|, within-galaxy scatter) = {r_abs_ws:+.4f}")

# RC roughness
roughness = np.array([g['roughness'] for g in galaxies])
rough_valid = np.isfinite(roughness)
if rough_valid.sum() > 20:
    r_rough = np.corrcoef(resid[rough_valid], roughness[rough_valid])[0, 1]
    r_abs_rough = np.corrcoef(np.abs(resid[rough_valid]), roughness[rough_valid])[0, 1]
    print(f"  r(resid, RC roughness) = {r_rough:+.4f}")
    print(f"  r(|resid|, RC roughness) = {r_abs_rough:+.4f}")

print("\n✓ Test 3 passed: quality indicators tested")

# =====================================================================
# TEST 4: INCLINATION SYSTEMATICS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: INCLINATION SYSTEMATICS")
print("=" * 60)

inclination = properties['inclination']
sin_incl = properties['sin_incl']

r_incl = np.corrcoef(resid, inclination)[0, 1]
r_sin = np.corrcoef(resid, sin_incl)[0, 1]

print(f"\n  r(resid, inclination) = {r_incl:+.4f}")
print(f"  r(resid, sin(incl)) = {r_sin:+.4f}")

# Inclination-dependent systematic: v_corr = v_obs / sin(i)
# If inclination is wrong by δi, offset shifts by ~2δi/tan(i)
# This should show up as a correlation at low inclination

# Split by inclination
low_incl = inclination < 45
mid_incl = (inclination >= 45) & (inclination <= 70)
high_incl = inclination > 70

for label, mask in [('Low (<45°)', low_incl), ('Mid (45-70°)', mid_incl), ('High (>70°)', high_incl)]:
    if mask.sum() >= 5:
        print(f"  {label}: N={mask.sum():>3}, mean resid={np.mean(resid[mask]):+.5f}, std={np.std(resid[mask]):.4f}")

# Test: do face-on galaxies have more scatter?
if low_incl.sum() >= 5 and high_incl.sum() >= 5:
    F_var = np.var(resid[low_incl]) / np.var(resid[high_incl])
    print(f"\n  Variance ratio (low/high incl): {F_var:.3f}")
    print(f"  Face-on galaxies {'MORE' if F_var > 1 else 'LESS'} scattered")

# Partial correlation with inclination controlling V and L
X_vl = np.column_stack([np.ones(n), logV, logL])
beta_r_vl = np.linalg.lstsq(X_vl, resid, rcond=None)[0]
resid_r_vl = resid - X_vl @ beta_r_vl
beta_i_vl = np.linalg.lstsq(X_vl, inclination, rcond=None)[0]
resid_i_vl = inclination - X_vl @ beta_i_vl
partial_r_incl = np.corrcoef(resid_r_vl, resid_i_vl)[0, 1]
print(f"\n  Partial r(resid, inclination | V, L) = {partial_r_incl:+.4f}")

print("\n✓ Test 4 passed: inclination systematics tested")

# =====================================================================
# TEST 5: NEAREST-NEIGHBOR RESIDUAL CORRELATION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: NEAREST-NEIGHBOR RESIDUAL CORRELATION")
print("=" * 60)

# Are galaxies with similar properties residual-correlated?
# Use the 6-var predictor space for nearest neighbors

X_pred = X6[:, 1:]  # 6 predictors without intercept
X_std = (X_pred - X_pred.mean(axis=0)) / X_pred.std(axis=0)

# Find nearest neighbor for each galaxy
nn_resids = np.zeros(n)
nn_dists = np.zeros(n)
for i in range(n):
    dists = np.sqrt(np.sum((X_std - X_std[i])**2, axis=1))
    dists[i] = np.inf  # exclude self
    j = np.argmin(dists)
    nn_resids[i] = resid[j]
    nn_dists[i] = dists[j]

r_nn = np.corrcoef(resid, nn_resids)[0, 1]
print(f"\n  r(resid, nearest-neighbor resid) = {r_nn:+.4f}")

# This should be ~0 if the model is complete (no missing variables)
# Positive = systematic patterns remain; Negative = noise-driven

# 5 nearest neighbors
nn5_resids = np.zeros(n)
for i in range(n):
    dists = np.sqrt(np.sum((X_std - X_std[i])**2, axis=1))
    dists[i] = np.inf
    nn5_idx = np.argsort(dists)[:5]
    nn5_resids[i] = np.mean(resid[nn5_idx])

r_nn5 = np.corrcoef(resid, nn5_resids)[0, 1]
print(f"  r(resid, 5-NN mean resid) = {r_nn5:+.4f}")

# Type-based clustering: galaxies of the same Hubble type
types = properties['hubble_type']
type_mean_resid = np.zeros(n)
for i in range(n):
    same_type = (types == types[i]) & (np.arange(n) != i)
    if same_type.sum() > 0:
        type_mean_resid[i] = np.mean(resid[same_type])
    else:
        type_mean_resid[i] = 0

r_type_nn = np.corrcoef(resid, type_mean_resid)[0, 1]
print(f"  r(resid, same-type mean resid) = {r_type_nn:+.4f}")

# Moran's I for spatial autocorrelation (using predictor-space distance)
W = 1.0 / (np.sqrt(np.sum((X_std[:, None, :] - X_std[None, :, :])**2, axis=2)) + 1e-10)
np.fill_diagonal(W, 0)
W_sum = W.sum()
resid_c = resid - np.mean(resid)
morans_I = n / W_sum * np.sum(W * np.outer(resid_c, resid_c)) / np.sum(resid_c**2)
# Expected under null: -1/(n-1)
E_I = -1 / (n - 1)
print(f"\n  Moran's I (predictor space): {morans_I:+.4f}")
print(f"  Expected (null): {E_I:+.4f}")
print(f"  Excess: {morans_I - E_I:+.4f}")

print("\n✓ Test 5 passed: nearest-neighbor correlation tested")

# =====================================================================
# TEST 6: RESIDUAL DISTRIBUTION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: RESIDUAL NORMALITY AND STRUCTURE")
print("=" * 60)

# Shapiro-Wilk test
w_stat, w_p = sp_stats.shapiro(resid)
print(f"\n  Shapiro-Wilk test: W = {w_stat:.4f}, p = {w_p:.4f}")
print(f"  {'Normal' if w_p > 0.05 else 'Non-normal'} at 5% level")

# D'Agostino-Pearson test
k2_stat, k2_p = sp_stats.normaltest(resid)
print(f"  D'Agostino-Pearson: K² = {k2_stat:.4f}, p = {k2_p:.4f}")

# Skewness and kurtosis
skew = sp_stats.skew(resid)
kurt = sp_stats.kurtosis(resid)
print(f"\n  Skewness: {skew:+.4f} (normal: 0)")
print(f"  Kurtosis: {kurt:+.4f} (normal: 0)")

# Heavy tails?
frac_1sigma = np.mean(np.abs(resid) < np.std(resid))
frac_2sigma = np.mean(np.abs(resid) < 2 * np.std(resid))
frac_3sigma = np.mean(np.abs(resid) < 3 * np.std(resid))
print(f"\n  Fraction within:")
print(f"    1σ: {frac_1sigma:.3f} (normal: 0.683)")
print(f"    2σ: {frac_2sigma:.3f} (normal: 0.954)")
print(f"    3σ: {frac_3sigma:.3f} (normal: 0.997)")

# Outliers: galaxies with |resid| > 2σ
sigma = np.std(resid)
outliers_2s = np.where(np.abs(resid) > 2 * sigma)[0]
print(f"\n  Outliers beyond 2σ ({2*sigma:.4f} dex): {len(outliers_2s)}")
for idx in outliers_2s:
    g = galaxies[idx]
    print(f"    {g['id']}: resid={resid[idx]:+.4f}, logV={g['logV']:.3f}, type={g['hubble_type']}")

# QQ plot statistics: deviation from normal at tails
sorted_resid = np.sort(resid)
normal_quantiles = sp_stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n) * np.std(resid) + np.mean(resid)
qq_correlation = np.corrcoef(sorted_resid, normal_quantiles)[0, 1]
print(f"\n  QQ plot correlation: {qq_correlation:.6f} (1.0 = perfect normal)")

print("\n✓ Test 6 passed: residual distribution analyzed")

# =====================================================================
# TEST 7: MOND REGIME DEPENDENCE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: RESIDUAL vs MOND REGIME")
print("=" * 60)

log_g = properties['log_g_ratio']

r_g_resid = np.corrcoef(resid, log_g)[0, 1]
print(f"\n  r(resid, log(g_bar/a₀)) = {r_g_resid:+.4f}")
print(f"  (Should be ~0 if interpolation function is adequate)")

# Partial correlation controlling all model variables
# The 6-var model residual should already be orthogonal to the predictors
# But log_g_ratio is NOT a predictor, so it can still correlate
# This tests the interpolation function

# Split by MOND depth
g_terciles = np.percentile(log_g, [33, 67])
deep = log_g <= g_terciles[0]
moderate = (log_g > g_terciles[0]) & (log_g <= g_terciles[1])
shallow = log_g > g_terciles[1]

print(f"\n  By MOND depth:")
for label, mask in [('Deep', deep), ('Moderate', moderate), ('Shallow', shallow)]:
    print(f"    {label}: N={mask.sum()}, mean resid={np.mean(resid[mask]):+.5f}, std={np.std(resid[mask]):.4f}")

# Does the residual have a quadratic dependence on log_g?
X_g = np.column_stack([np.ones(n), log_g])
X_g2 = np.column_stack([np.ones(n), log_g, log_g**2])
_, _, _, R2_g1, _ = build_model(X_g, resid)
_, _, _, R2_g2, _ = build_model(X_g2, resid)
F_quad = (R2_g2 - R2_g1) / (1 - R2_g2) * (n - 3)
p_quad = 1 - sp_stats.f.cdf(F_quad, 1, n - 3)

print(f"\n  Linear fit of resid on log(g/a₀): R² = {R2_g1:.4f}")
print(f"  Adding quadratic: R² = {R2_g2:.4f}")
print(f"  F-test for quadratic term: F = {F_quad:.3f}, p = {p_quad:.4f}")

# If significant, the interpolation function has a regime-dependent bias
if p_quad < 0.05:
    print(f"  → Quadratic term SIGNIFICANT: interpolation function shows regime dependence")
else:
    print(f"  → Quadratic term not significant: interpolation function adequate")

print("\n✓ Test 7 passed: MOND regime dependence tested")

# =====================================================================
# TEST 8: THE IRREDUCIBLE FLOOR
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: THE IRREDUCIBLE RESIDUAL FLOOR")
print("=" * 60)

# Method 1: Split-half stability
# If residuals are pure noise, split-half correlation = 0
# If residuals contain signal, split-half correlation > 0
rng = np.random.RandomState(509)

# For each galaxy, if we had multiple independent measurements,
# we could compute reliability. Instead, use the within-galaxy scatter
# as a proxy for measurement noise

ws = properties['within_scatter']
ws_valid = np.isfinite(ws)

# Expected noise contribution to residual:
# σ_noise ≈ within_scatter / sqrt(n_mond_pts) (error of mean)
n_mond = properties['n_mond_pts']
sigma_noise_est = np.zeros(n)
for i in range(n):
    if np.isfinite(ws[i]) and n_mond[i] > 0:
        sigma_noise_est[i] = ws[i] / np.sqrt(n_mond[i])
    else:
        sigma_noise_est[i] = np.nan

valid_noise = np.isfinite(sigma_noise_est)
mean_noise_est = np.mean(sigma_noise_est[valid_noise])
rms_noise_est = np.sqrt(np.mean(sigma_noise_est[valid_noise]**2))

print(f"\n  Estimated measurement noise in offset:")
print(f"  Mean σ_noise: {mean_noise_est:.4f} dex")
print(f"  RMS σ_noise: {rms_noise_est:.4f} dex")
print(f"  Model residual RMS: {rms_6:.4f} dex")

noise_fraction = rms_noise_est**2 / rms_6**2
print(f"\n  Noise variance fraction: {noise_fraction*100:.1f}%")
print(f"  Signal variance fraction: {(1-noise_fraction)*100:.1f}%")

# Method 2: Bootstrap residual reproducibility
# Resample galaxies, refit model, correlate residuals of held-out galaxies
B_boot = 1000
resid_corrs = []
for b in range(B_boot):
    idx_train = rng.choice(n, size=n, replace=True)
    in_train = np.zeros(n, dtype=bool)
    in_train[idx_train] = True
    idx_test = np.where(~in_train)[0]

    if len(idx_test) < 10:
        continue

    X_train = X6[idx_train]
    y_train = offset[idx_train]
    beta_b = np.linalg.lstsq(X_train, y_train, rcond=None)[0]

    pred_test = X6[idx_test] @ beta_b
    resid_test = offset[idx_test] - pred_test

    # Compare with original residuals at test points
    r_boot = np.corrcoef(resid_test, resid[idx_test])[0, 1]
    resid_corrs.append(r_boot)

resid_corrs = np.array(resid_corrs)
print(f"\n  Bootstrap OOB residual correlation (B={B_boot}):")
print(f"  Mean r(resid_OOB, resid_full): {np.mean(resid_corrs):+.4f}")
print(f"  Std: {np.std(resid_corrs):.4f}")
print(f"  If close to 1: residuals are stable (signal)")
print(f"  If close to 0: residuals are noise")

# Method 3: Augmented model — add ALL catalog properties
# If R² doesn't improve, the residual is truly uninformative
extra_vars = []
extra_names = []
for name in ['distance', 'inclination', 'hubble_type', 'n_points', 'log_g_ratio']:
    vals = properties[name]
    if np.all(np.isfinite(vals)):
        extra_vars.append(vals)
        extra_names.append(name)

if extra_vars:
    X_aug = np.column_stack([X6] + [v.reshape(-1, 1) for v in extra_vars])
    _, _, _, R2_aug, rms_aug = build_model(X_aug, offset)
    loo_aug = loo_r2(X_aug, offset)

    print(f"\n  Augmented model (+{len(extra_names)} catalog vars):")
    print(f"  R² = {R2_aug:.4f} (was {R2_6:.4f}, Δ = {R2_aug - R2_6:+.4f})")
    print(f"  LOO = {loo_aug:.4f} (was {loo_r2(X6, offset):.4f})")
    print(f"  Added: {', '.join(extra_names)}")

    # F-test for the extra variables
    df1 = len(extra_names)
    df2 = n - X_aug.shape[1]
    F_aug = ((R2_aug - R2_6) / df1) / ((1 - R2_aug) / df2)
    p_aug = 1 - sp_stats.f.cdf(F_aug, df1, df2)
    print(f"  F-test: F = {F_aug:.3f}, p = {p_aug:.4f}")

# Summary of the irreducible floor
print(f"\n  SUMMARY OF RESIDUAL BUDGET:")
print(f"  Total residual variance: {rms_6**2:.6f} dex²")
print(f"  Estimated noise: {rms_noise_est**2:.6f} dex² ({noise_fraction*100:.0f}%)")
print(f"  Physical residual: {rms_6**2 - rms_noise_est**2:.6f} dex² ({(1-noise_fraction)*100:.0f}%)")
print(f"  Physical RMS: {np.sqrt(max(0, rms_6**2 - rms_noise_est**2)):.4f} dex")

print("\n✓ Test 8 passed: irreducible floor estimated")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #509 SUMMARY")
print("=" * 70)
print(f"\nResidual RMS: {rms_6:.4f} dex")
print(f"Noise contribution: ~{noise_fraction*100:.0f}%")
print(f"NN residual correlation: {r_nn:+.4f}")
print(f"Shapiro-Wilk p: {w_p:.4f}")
if significant:
    print(f"Significant residual correlations:")
    for name, r, p in significant:
        print(f"  {name}: r = {r:+.4f}")
else:
    print(f"No significant residual correlations found")
print(f"\nAll 8 tests passed ✓")
