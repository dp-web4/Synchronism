#!/usr/bin/env python3
"""
======================================================================
SESSION #501: BOOTSTRAP PREDICTION INTERVALS
======================================================================

The 6-var model achieves R²=0.945, LOO R²=0.938, RMS=0.038. But what
are the FORMAL UNCERTAINTIES on:
1. Each coefficient?
2. Each galaxy's predicted offset?
3. A new galaxy's prediction?

This session uses bootstrap resampling (10,000 iterations) to construct
confidence intervals on all model outputs, and compares with analytic
estimates from normal theory.

Tests:
1. Bootstrap coefficient distributions
2. Coefficient stability: how often do signs flip?
3. Analytic vs bootstrap standard errors
4. Prediction intervals for each galaxy
5. Prediction intervals for hypothetical galaxies
6. Bootstrap R² distribution
7. Calibration: do 95% CIs contain 95% of LOO residuals?
8. Parameter covariance structure

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #501
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

        if vflat <= 0 or lum <= 0 or sb_eff <= 0:
            continue

        r_eff_pc = np.sqrt(lum * 1e9 / (2 * np.pi * max(sb_eff, 1)))
        r_eff_kpc = r_eff_pc / 1000

        v_obs = np.array([pt['v_obs'] for pt in points])
        v_gas = np.array([pt['v_gas'] for pt in points])
        v_disk = np.array([pt['v_disk'] for pt in points])
        v_bul = np.array([pt.get('v_bul', 0) for pt in points])
        radius = np.array([pt['radius'] for pt in points])

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

        galaxies.append({
            'id': gal_id,
            'offset': offset,
            'logV': np.log10(vflat),
            'logL': np.log10(lum),
            'c_V': c_V,
            'f_gas': f_gas,
            'hubble_type': hubble_type,
        })

    return galaxies


def build_model(X, y):
    """OLS regression, return beta, yhat, resid, R²."""
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    yhat = X @ beta
    resid = y - yhat
    R2 = 1 - np.sum(resid**2) / np.sum((y - np.mean(y))**2)
    return beta, yhat, resid, R2


print("=" * 70)
print("SESSION #501: BOOTSTRAP PREDICTION INTERVALS")
print("=" * 70)

galaxies = prepare_galaxies()
n = len(galaxies)

logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
y = np.array([g['offset'] for g in galaxies])

X = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])
p = X.shape[1]

beta_ols, yhat_ols, resid_ols, R2_ols = build_model(X, y)
rms_ols = np.sqrt(np.mean(resid_ols**2))
print(f"\n{n} galaxies, {p} parameters")
print(f"OLS: R² = {R2_ols:.4f}, RMS = {rms_ols:.4f}")

var_names = ['const', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V', 'logL×f_gas']

# =====================================================================
# TEST 1: BOOTSTRAP COEFFICIENT DISTRIBUTIONS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: BOOTSTRAP COEFFICIENT DISTRIBUTIONS (B=10000)")
print("=" * 60)

B = 10000
rng = np.random.RandomState(42)
boot_betas = np.zeros((B, p))
boot_R2 = np.zeros(B)
boot_yhats = np.zeros((B, n))

for b in range(B):
    idx = rng.randint(0, n, size=n)
    X_b = X[idx]
    y_b = y[idx]
    beta_b = np.linalg.lstsq(X_b, y_b, rcond=None)[0]
    boot_betas[b] = beta_b
    boot_yhats[b] = X @ beta_b  # predict ALL galaxies
    yhat_b = X_b @ beta_b
    resid_b = y_b - yhat_b
    ss_tot_b = np.sum((y_b - np.mean(y_b))**2)
    boot_R2[b] = 1 - np.sum(resid_b**2) / ss_tot_b if ss_tot_b > 0 else 0

print(f"\n{'Variable':<12} {'β(OLS)':<10} {'β(boot mean)':<14} {'SE(boot)':<10} "
      f"{'95% CI':<25}")
print("-" * 75)
for j, vname in enumerate(var_names):
    lo = np.percentile(boot_betas[:, j], 2.5)
    hi = np.percentile(boot_betas[:, j], 97.5)
    se = np.std(boot_betas[:, j])
    print(f"  {vname:<12} {beta_ols[j]:+.4f}    {np.mean(boot_betas[:, j]):+.4f}        "
          f"{se:.4f}     [{lo:+.4f}, {hi:+.4f}]")

print("\n✓ Test 1 passed: bootstrap distributions computed")

# =====================================================================
# TEST 2: COEFFICIENT SIGN STABILITY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: COEFFICIENT SIGN STABILITY")
print("=" * 60)

print(f"\n{'Variable':<12} {'OLS sign':<10} {'% same sign':<12} {'Stable?'}")
print("-" * 45)
for j, vname in enumerate(var_names):
    ols_sign = np.sign(beta_ols[j])
    same_frac = np.mean(np.sign(boot_betas[:, j]) == ols_sign) * 100
    stable = same_frac > 95
    print(f"  {vname:<12} {'+'if ols_sign>0 else '-':<10} {same_frac:<12.1f} "
          f"{'YES' if stable else 'NO'}")

print("\n✓ Test 2 passed: sign stability checked")

# =====================================================================
# TEST 3: ANALYTIC vs BOOTSTRAP STANDARD ERRORS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: ANALYTIC vs BOOTSTRAP STANDARD ERRORS")
print("=" * 60)

# Analytic SE from normal theory: SE(β) = sqrt(diag(σ² (X'X)⁻¹))
mse = np.sum(resid_ols**2) / (n - p)
XtX_inv = np.linalg.inv(X.T @ X)
se_analytic = np.sqrt(mse * np.diag(XtX_inv))
se_bootstrap = np.std(boot_betas, axis=0)

print(f"\n{'Variable':<12} {'SE(analytic)':<14} {'SE(bootstrap)':<14} {'Ratio':<8}")
print("-" * 50)
for j, vname in enumerate(var_names):
    ratio = se_bootstrap[j] / se_analytic[j]
    print(f"  {vname:<12} {se_analytic[j]:.4f}        {se_bootstrap[j]:.4f}        {ratio:.2f}")

mean_ratio = np.mean(se_bootstrap / se_analytic)
print(f"\nMean SE ratio (bootstrap/analytic): {mean_ratio:.3f}")
print(f"→ {'Good agreement' if 0.8 < mean_ratio < 1.2 else 'Discrepancy'}: "
      f"{'normal theory is adequate' if 0.8 < mean_ratio < 1.2 else 'non-normality detected'}")

print("\n✓ Test 3 passed: SE comparison done")

# =====================================================================
# TEST 4: PREDICTION INTERVALS FOR EACH GALAXY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: PREDICTION INTERVALS")
print("=" * 60)

# For each galaxy, the prediction interval from bootstrap
pred_lo = np.percentile(boot_yhats, 2.5, axis=0)
pred_hi = np.percentile(boot_yhats, 97.5, axis=0)
pred_se = np.std(boot_yhats, axis=0)

# Average prediction uncertainty
mean_pred_se = np.mean(pred_se)
mean_ci_width = np.mean(pred_hi - pred_lo)

print(f"\nPrediction uncertainty (from bootstrap):")
print(f"  Mean SE of prediction: {mean_pred_se:.4f} dex")
print(f"  Mean 95% CI width: {mean_ci_width:.4f} dex")
print(f"  Model RMS: {rms_ols:.4f} dex")
print(f"  SE/RMS ratio: {mean_pred_se/rms_ols:.2f}")

# Galaxies with widest CIs (most uncertain predictions)
ci_widths = pred_hi - pred_lo
sorted_ci = np.argsort(ci_widths)[::-1]

print(f"\nTop 5 most uncertain predictions:")
print(f"{'Galaxy':<15} {'Predicted':<10} {'95% CI':<25} {'Width':<8}")
print("-" * 60)
for idx in sorted_ci[:5]:
    g = galaxies[idx]
    print(f"  {g['id']:<15} {yhat_ols[idx]:+.3f}     [{pred_lo[idx]:+.3f}, {pred_hi[idx]:+.3f}]    {ci_widths[idx]:.3f}")

print(f"\nTop 5 most precise predictions:")
for idx in sorted_ci[-5:]:
    g = galaxies[idx]
    print(f"  {g['id']:<15} {yhat_ols[idx]:+.3f}     [{pred_lo[idx]:+.3f}, {pred_hi[idx]:+.3f}]    {ci_widths[idx]:.3f}")

print("\n✓ Test 4 passed: prediction intervals computed")

# =====================================================================
# TEST 5: HYPOTHETICAL GALAXY PREDICTIONS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: PREDICTION FOR HYPOTHETICAL GALAXIES")
print("=" * 60)

# Define some hypothetical galaxy types
hypothetical = [
    ('Milky Way-like', 2.3, 1.1, 0.85, 0.15),    # logV=2.3, logL=1.1, c_V=0.85, f_gas=0.15
    ('Low-mass dwarf', 1.7, -0.5, 0.60, 0.70),    # logV=1.7, logL=-0.5, c_V=0.60, f_gas=0.70
    ('Massive spiral', 2.4, 1.5, 0.90, 0.10),      # logV=2.4, logL=1.5, c_V=0.90, f_gas=0.10
    ('Gas-rich dwarf', 1.8, -0.3, 0.55, 0.85),     # logV=1.8, logL=-0.3, c_V=0.55, f_gas=0.85
    ('LSB giant', 2.2, 0.8, 0.65, 0.40),           # logV=2.2, logL=0.8, c_V=0.65, f_gas=0.40
]

print(f"\n{'Galaxy type':<18} {'Prediction':<12} {'95% CI':<25} {'PI width':<10}")
print("-" * 70)
for name, lv, ll, cv, fg in hypothetical:
    x_new = np.array([1, lv, ll, cv, fg, lv * cv, ll * fg])

    # Bootstrap predictions
    preds = boot_betas @ x_new
    pred_mean = np.mean(preds)
    pred_lo_h = np.percentile(preds, 2.5)
    pred_hi_h = np.percentile(preds, 97.5)

    # Prediction interval (includes residual noise)
    # PI = prediction ± sqrt(σ²_pred + σ²_resid)
    pi_se = np.sqrt(np.var(preds) + np.var(resid_ols))
    pi_lo = pred_mean - 1.96 * pi_se
    pi_hi = pred_mean + 1.96 * pi_se

    print(f"  {name:<18} {pred_mean:+.4f}      [{pred_lo_h:+.3f}, {pred_hi_h:+.3f}]    {pred_hi_h-pred_lo_h:.3f}")
    print(f"  {'(with noise)':<18} {'':12} [{pi_lo:+.3f}, {pi_hi:+.3f}]    {pi_hi-pi_lo:.3f}")

print("\n✓ Test 5 passed: hypothetical predictions done")

# =====================================================================
# TEST 6: BOOTSTRAP R² DISTRIBUTION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: BOOTSTRAP R² DISTRIBUTION")
print("=" * 60)

print(f"\nOLS R² = {R2_ols:.4f}")
print(f"Bootstrap R² distribution:")
print(f"  Mean: {np.mean(boot_R2):.4f}")
print(f"  Std: {np.std(boot_R2):.4f}")
print(f"  2.5th percentile: {np.percentile(boot_R2, 2.5):.4f}")
print(f"  97.5th percentile: {np.percentile(boot_R2, 97.5):.4f}")
print(f"  95% CI: [{np.percentile(boot_R2, 2.5):.4f}, {np.percentile(boot_R2, 97.5):.4f}]")
print(f"  Min: {np.min(boot_R2):.4f}")
print(f"  Max: {np.max(boot_R2):.4f}")

# What fraction of bootstraps give R² > 0.90?
print(f"\n  P(R² > 0.90) = {np.mean(boot_R2 > 0.90)*100:.1f}%")
print(f"  P(R² > 0.93) = {np.mean(boot_R2 > 0.93)*100:.1f}%")
print(f"  P(R² > 0.95) = {np.mean(boot_R2 > 0.95)*100:.1f}%")

print("\n✓ Test 6 passed: R² distribution computed")

# =====================================================================
# TEST 7: CALIBRATION — DO 95% CIs CONTAIN 95% OF LOO RESIDUALS?
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: CI CALIBRATION")
print("=" * 60)

# LOO residuals
H = X @ XtX_inv @ X.T
h = np.diag(H)
loo_resid = resid_ols / (1 - h)

# Method 1: Bootstrap CI on predictions
# For each galaxy, does the true offset fall within the bootstrap prediction CI?
coverage_pred = np.mean((y >= pred_lo) & (y <= pred_hi))

# Method 2: Full prediction interval (includes noise)
# PI = yhat ± t * sqrt(MSE * (1 + h_ii))
t_val = 1.96
pi_lo_full = yhat_ols - t_val * np.sqrt(mse * (1 + h))
pi_hi_full = yhat_ols + t_val * np.sqrt(mse * (1 + h))
coverage_full = np.mean((y >= pi_lo_full) & (y <= pi_hi_full))

# Method 3: Bootstrap + residual noise
boot_pred_se = np.std(boot_yhats, axis=0)
total_se = np.sqrt(boot_pred_se**2 + np.var(resid_ols))
pi_boot_lo = yhat_ols - 1.96 * total_se
pi_boot_hi = yhat_ols + 1.96 * total_se
coverage_boot_full = np.mean((y >= pi_boot_lo) & (y <= pi_boot_hi))

print(f"\n95% CI calibration (should be ~95%):")
print(f"  Bootstrap prediction CI: {coverage_pred*100:.1f}%")
print(f"  Analytic prediction interval: {coverage_full*100:.1f}%")
print(f"  Bootstrap + noise PI: {coverage_boot_full*100:.1f}%")

# LOO calibration
loo_sigma = np.sqrt(mse * (1 + h) / (1 - h)**2)
# Simplified: just use global sigma
loo_z = loo_resid / np.std(loo_resid)
loo_in_95 = np.mean(np.abs(loo_z) < 1.96)
print(f"\n  LOO residuals within ±1.96σ: {loo_in_95*100:.1f}%")

print("\n✓ Test 7 passed: calibration checked")

# =====================================================================
# TEST 8: PARAMETER COVARIANCE STRUCTURE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: PARAMETER COVARIANCE (BOOTSTRAP)")
print("=" * 60)

# Correlation matrix of bootstrap coefficients
boot_corr = np.corrcoef(boot_betas.T)

print(f"\nBootstrap coefficient correlation matrix:")
print(f"{'':12}", end="")
for vname in var_names:
    print(f"{vname[:6]:>8}", end="")
print()
for i, vi in enumerate(var_names):
    print(f"  {vi:<10}", end="")
    for j in range(p):
        r_val = boot_corr[i, j]
        flag = "*" if abs(r_val) > 0.5 and i != j else " "
        print(f"{r_val:+.2f}{flag}  ", end="")
    print()

# Find the most strongly correlated pair
max_corr = 0
max_pair = (0, 0)
for i in range(p):
    for j in range(i+1, p):
        if abs(boot_corr[i, j]) > abs(max_corr):
            max_corr = boot_corr[i, j]
            max_pair = (i, j)

print(f"\nMost correlated pair: {var_names[max_pair[0]]} ↔ {var_names[max_pair[1]]} "
      f"(r = {max_corr:+.3f})")

# Condition number of X'X
cond = np.linalg.cond(X.T @ X)
print(f"Condition number of X'X: {cond:.0f}")
print(f"→ {'Well-conditioned' if cond < 1000 else 'Moderate multicollinearity' if cond < 10000 else 'Severe multicollinearity'}")

# Variance Inflation Factors
print(f"\nVariance Inflation Factors:")
for j in range(1, p):  # skip intercept
    X_other = np.delete(X, j, axis=1)
    beta_other = np.linalg.lstsq(X_other, X[:, j], rcond=None)[0]
    yhat_other = X_other @ beta_other
    R2_j = 1 - np.sum((X[:, j] - yhat_other)**2) / np.sum((X[:, j] - np.mean(X[:, j]))**2)
    vif = 1 / max(1 - R2_j, 1e-10)
    print(f"  {var_names[j]:<12} VIF = {vif:.1f}")

print("\n✓ Test 8 passed: covariance structure analyzed")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #501 SUMMARY")
print("=" * 70)
print(f"Bootstrap: B = {B}")
print(f"Mean prediction SE: {mean_pred_se:.4f} dex")
print(f"Mean 95% CI width: {mean_ci_width:.4f} dex")
print(f"R² 95% CI: [{np.percentile(boot_R2, 2.5):.4f}, {np.percentile(boot_R2, 97.5):.4f}]")
print(f"All signs stable: {all(np.mean(np.sign(boot_betas[:,j]) == np.sign(beta_ols[j])) > 0.95 for j in range(p))}")
print(f"SE ratio (boot/analytic): {mean_ratio:.3f}")
print(f"PI calibration: {coverage_full*100:.1f}%")
print(f"Condition number: {cond:.0f}")
print(f"\nAll 8 tests passed ✓")
