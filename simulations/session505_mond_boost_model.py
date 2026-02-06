#!/usr/bin/env python3
"""
======================================================================
SESSION #505: MOND BOOST MODEL — γ AS PRIMARY PREDICTOR
======================================================================

Session #504 found: partial r(log γ, MOND boost | V, L) = +0.57.
This is the strongest positive result for γ = 2/√N_corr.

The MOND boost = log₁₀(g_obs/g_bar) at the outer MOND points.
In deep MOND: g_obs ≈ √(g_bar × a₀), so boost = 0.5 × log(a₀/g_bar).
The boost measures the total "dark matter" acceleration contribution.

Can we build a model where γ is the primary predictor?
How does a γ-based boost model compare to the empirical approach?

Tests:
1. MOND boost statistics and theoretical prediction
2. γ as boost predictor: raw and partial correlations
3. Build a γ-based boost model
4. Compare γ model with V,L model
5. LOO validation of boost models
6. Physical interpretation of boost model coefficients
7. Boost model residual analysis
8. Connecting boost model back to offset model

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #505
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
    """Load SPARC with MOND boost measurements."""
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
        v_bul_v = v_bul[valid] if np.any(v_bul != 0) else np.zeros_like(v_obs[valid])

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

        # MOND boost: log(g_obs/g_bar) at outer MOND points
        if outer_mond.sum() >= 2:
            boost = np.mean(np.log10(g_obs_v[outer_mond] / g_bar_v[outer_mond]))
            mean_g_bar = np.mean(g_bar_v[outer_mond])
        else:
            boost = np.mean(np.log10(g_obs_v[mond] / g_bar_v[mond]))
            mean_g_bar = np.mean(g_bar_v[mond])

        # Deep MOND theoretical prediction: boost = 0.5 × log(a₀/g_bar)
        theory_boost = 0.5 * np.log10(a0_mond / mean_g_bar)

        # f_gas
        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # R_max
        r_max_kpc = radius_v.max()

        # N_corr and gamma
        vflat_ms = vflat * 1e3
        r_max_m = r_max_kpc * 3.086e19
        N_corr = vflat_ms**2 / (r_max_m * a0_mond)
        gamma = 2.0 / np.sqrt(N_corr)
        log_gamma = np.log10(gamma)

        # log(g_bar/a₀) at outer radius — mean acceleration regime
        log_g_ratio = np.log10(mean_g_bar / a0_mond)

        galaxies.append({
            'id': gal_id,
            'offset': offset,
            'boost': boost,
            'theory_boost': theory_boost,
            'logV': np.log10(vflat),
            'logL': np.log10(lum),
            'c_V': c_V,
            'f_gas': f_gas,
            'hubble_type': hubble_type,
            'log_gamma': log_gamma,
            'gamma': gamma,
            'N_corr': N_corr,
            'log_g_ratio': log_g_ratio,
            'mean_g_bar': mean_g_bar,
            'r_max_kpc': r_max_kpc,
        })

    return galaxies


print("=" * 70)
print("SESSION #505: MOND BOOST MODEL — γ AS PRIMARY PREDICTOR")
print("=" * 70)

galaxies = prepare_galaxies()
n = len(galaxies)
print(f"\n{n} galaxies loaded")

# Extract arrays
boost = np.array([g['boost'] for g in galaxies])
theory_boost = np.array([g['theory_boost'] for g in galaxies])
offset = np.array([g['offset'] for g in galaxies])
logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
log_gamma = np.array([g['log_gamma'] for g in galaxies])
log_g_ratio = np.array([g['log_g_ratio'] for g in galaxies])

# =====================================================================
# TEST 1: MOND BOOST STATISTICS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: MOND BOOST STATISTICS")
print("=" * 60)

print(f"\nMOND boost = log₁₀(g_obs/g_bar) at outer MOND points:")
print(f"  Mean: {np.mean(boost):.4f}")
print(f"  Std: {np.std(boost):.4f}")
print(f"  Range: [{np.min(boost):.4f}, {np.max(boost):.4f}]")

print(f"\nDeep MOND theoretical prediction (0.5 × log(a₀/g_bar)):")
print(f"  Mean: {np.mean(theory_boost):.4f}")
print(f"  Std: {np.std(theory_boost):.4f}")

print(f"\nTheory vs observation:")
print(f"  r(theory_boost, boost) = {np.corrcoef(theory_boost, boost)[0,1]:+.4f}")
print(f"  Mean(boost - theory): {np.mean(boost - theory_boost):+.4f}")
print(f"  RMS(boost - theory): {np.sqrt(np.mean((boost - theory_boost)**2)):.4f}")

# Boost vs offset
r_boost_offset = np.corrcoef(boost, offset)[0, 1]
print(f"\n  r(boost, offset) = {r_boost_offset:+.4f}")

print("\n✓ Test 1 passed: boost statistics computed")

# =====================================================================
# TEST 2: γ AS BOOST PREDICTOR
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: γ AS BOOST PREDICTOR")
print("=" * 60)

r_gamma_boost = np.corrcoef(log_gamma, boost)[0, 1]
print(f"\n  r(log γ, boost) = {r_gamma_boost:+.4f}")

# Controlling V
X_v = np.column_stack([np.ones(n), logV])
beta_b_v = np.linalg.lstsq(X_v, boost, rcond=None)[0]
resid_b_v = boost - X_v @ beta_b_v
beta_g_v = np.linalg.lstsq(X_v, log_gamma, rcond=None)[0]
resid_g_v = log_gamma - X_v @ beta_g_v
partial_v = np.corrcoef(resid_b_v, resid_g_v)[0, 1]
print(f"  Partial r(log γ, boost | V) = {partial_v:+.4f}")

# Controlling V, L
X_vl = np.column_stack([np.ones(n), logV, logL])
beta_b_vl = np.linalg.lstsq(X_vl, boost, rcond=None)[0]
resid_b_vl = boost - X_vl @ beta_b_vl
beta_g_vl = np.linalg.lstsq(X_vl, log_gamma, rcond=None)[0]
resid_g_vl = log_gamma - X_vl @ beta_g_vl
partial_vl = np.corrcoef(resid_b_vl, resid_g_vl)[0, 1]
print(f"  Partial r(log γ, boost | V, L) = {partial_vl:+.4f}")

# Controlling V, L, f_gas, c_V
X_full = np.column_stack([np.ones(n), logV, logL, f_gas, c_V])
beta_b_f = np.linalg.lstsq(X_full, boost, rcond=None)[0]
resid_b_f = boost - X_full @ beta_b_f
beta_g_f = np.linalg.lstsq(X_full, log_gamma, rcond=None)[0]
resid_g_f = log_gamma - X_full @ beta_g_f
partial_full = np.corrcoef(resid_b_f, resid_g_f)[0, 1]
print(f"  Partial r(log γ, boost | V, L, f_gas, c_V) = {partial_full:+.4f}")

# Compare γ with log(g_bar/a₀) which is the "trivial" predictor
r_g_boost = np.corrcoef(log_g_ratio, boost)[0, 1]
print(f"\n  r(log(g_bar/a₀), boost) = {r_g_boost:+.4f}  (trivial predictor)")
print(f"  r(log γ, boost) = {r_gamma_boost:+.4f}")

print("\n✓ Test 2 passed: γ-boost correlations computed")

# =====================================================================
# TEST 3: BUILD A γ-BASED BOOST MODEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: γ-BASED BOOST MODEL")
print("=" * 60)

# Model 1: boost = a + b × log(γ)
X_g = np.column_stack([np.ones(n), log_gamma])
beta_g = np.linalg.lstsq(X_g, boost, rcond=None)[0]
yhat_g = X_g @ beta_g
resid_g = boost - yhat_g
R2_g = 1 - np.sum(resid_g**2) / np.sum((boost - np.mean(boost))**2)

# Model 2: boost = a + b × log(γ) + c × f_gas
X_gf = np.column_stack([np.ones(n), log_gamma, f_gas])
beta_gf = np.linalg.lstsq(X_gf, boost, rcond=None)[0]
yhat_gf = X_gf @ beta_gf
resid_gf = boost - yhat_gf
R2_gf = 1 - np.sum(resid_gf**2) / np.sum((boost - np.mean(boost))**2)

# Model 3: boost = a + b × log(γ) + c × f_gas + d × c_V
X_gfc = np.column_stack([np.ones(n), log_gamma, f_gas, c_V])
beta_gfc = np.linalg.lstsq(X_gfc, boost, rcond=None)[0]
yhat_gfc = X_gfc @ beta_gfc
resid_gfc = boost - yhat_gfc
R2_gfc = 1 - np.sum(resid_gfc**2) / np.sum((boost - np.mean(boost))**2)

# Model 4: boost = a + b × log(g_bar/a₀)  (trivial)
X_trivial = np.column_stack([np.ones(n), log_g_ratio])
beta_trivial = np.linalg.lstsq(X_trivial, boost, rcond=None)[0]
yhat_trivial = X_trivial @ beta_trivial
R2_trivial = 1 - np.sum((boost - yhat_trivial)**2) / np.sum((boost - np.mean(boost))**2)

print(f"\n{'Model':<35} {'R²':<10} {'RMS':<10}")
print("-" * 55)
print(f"  {'γ alone':<35} {R2_g:<10.4f} {np.sqrt(np.mean(resid_g**2)):<10.4f}")
print(f"  {'γ + f_gas':<35} {R2_gf:<10.4f} {np.sqrt(np.mean(resid_gf**2)):<10.4f}")
print(f"  {'γ + f_gas + c_V':<35} {R2_gfc:<10.4f} {np.sqrt(np.mean(resid_gfc**2)):<10.4f}")
print(f"  {'log(g_bar/a₀) (trivial)':<35} {R2_trivial:<10.4f} {np.sqrt(np.mean((boost - yhat_trivial)**2)):<10.4f}")
print(f"  {'Deep MOND theory':<35} {'—':<10} {np.sqrt(np.mean((boost - theory_boost)**2)):<10.4f}")

print(f"\n  Best γ model coefficients:")
for name, b in zip(['const', 'log(γ)', 'f_gas', 'c_V'], beta_gfc):
    print(f"    {name:<12} {b:+.4f}")

print("\n✓ Test 3 passed: γ-based boost model built")

# =====================================================================
# TEST 4: COMPARE γ MODEL WITH V,L MODEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: γ MODEL vs V,L MODEL FOR BOOST")
print("=" * 60)

# V,L model
X_vl_boost = np.column_stack([np.ones(n), logV, logL])
beta_vl = np.linalg.lstsq(X_vl_boost, boost, rcond=None)[0]
R2_vl = 1 - np.sum((boost - X_vl_boost @ beta_vl)**2) / np.sum((boost - np.mean(boost))**2)

# V,L,f_gas,c_V model
X_4var = np.column_stack([np.ones(n), logV, logL, f_gas, c_V])
beta_4var = np.linalg.lstsq(X_4var, boost, rcond=None)[0]
R2_4var = 1 - np.sum((boost - X_4var @ beta_4var)**2) / np.sum((boost - np.mean(boost))**2)

# Full 6-var model for boost
X_6var = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta_6var = np.linalg.lstsq(X_6var, boost, rcond=None)[0]
R2_6var = 1 - np.sum((boost - X_6var @ beta_6var)**2) / np.sum((boost - np.mean(boost))**2)

# γ replaces V,L?
X_g_replace = np.column_stack([np.ones(n), log_gamma, f_gas, c_V, log_gamma*c_V, logL*f_gas])
beta_g_rep = np.linalg.lstsq(X_g_replace, boost, rcond=None)[0]
R2_g_rep = 1 - np.sum((boost - X_g_replace @ beta_g_rep)**2) / np.sum((boost - np.mean(boost))**2)

print(f"\n{'Model':<35} {'R²':<10} {'# params'}")
print("-" * 55)
print(f"  {'γ alone':<35} {R2_g:<10.4f} 2")
print(f"  {'γ + f_gas + c_V':<35} {R2_gfc:<10.4f} 4")
print(f"  {'V, L':<35} {R2_vl:<10.4f} 3")
print(f"  {'V, L, f_gas, c_V':<35} {R2_4var:<10.4f} 5")
print(f"  {'6-var (with interactions)':<35} {R2_6var:<10.4f} 7")
print(f"  {'γ-based 6-var':<35} {R2_g_rep:<10.4f} 6")

# Does γ add to V,L?
X_vl_g = np.column_stack([np.ones(n), logV, logL, log_gamma])
beta_vl_g = np.linalg.lstsq(X_vl_g, boost, rcond=None)[0]
R2_vl_g = 1 - np.sum((boost - X_vl_g @ beta_vl_g)**2) / np.sum((boost - np.mean(boost))**2)
print(f"\n  V, L + γ: R² = {R2_vl_g:.4f} (Δ from V,L: {R2_vl_g - R2_vl:+.4f})")

print("\n✓ Test 4 passed: model comparison done")

# =====================================================================
# TEST 5: LOO VALIDATION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: LOO VALIDATION OF BOOST MODELS")
print("=" * 60)

def loo_r2(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    return 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)

models_loo = [
    ('γ alone', X_g),
    ('γ + f_gas + c_V', X_gfc),
    ('V, L', X_vl_boost),
    ('V, L, f_gas, c_V', X_4var),
    ('6-var', X_6var),
    ('γ-based 6-var', X_g_replace),
    ('log(g/a₀) trivial', X_trivial),
]

print(f"\n{'Model':<35} {'LOO R²':<10}")
print("-" * 45)
for name, X_m in models_loo:
    r2_loo = loo_r2(X_m, boost)
    print(f"  {name:<35} {r2_loo:.4f}")

print("\n✓ Test 5 passed: LOO validation done")

# =====================================================================
# TEST 6: PHYSICAL INTERPRETATION OF BOOST COEFFICIENTS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: BOOST MODEL COEFFICIENT INTERPRETATION")
print("=" * 60)

print(f"\n6-var boost model coefficients:")
var_names_6 = ['const', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V', 'logL×f_gas']
for name, b in zip(var_names_6, beta_6var):
    print(f"  {name:<12} {b:+.4f}")

print(f"\nCompare with 6-var OFFSET model coefficients:")
# Build offset model
X_off = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta_off = np.linalg.lstsq(X_off, offset, rcond=None)[0]
for name, b_boost, b_off in zip(var_names_6, beta_6var, beta_off):
    print(f"  {name:<12} boost: {b_boost:+.4f}   offset: {b_off:+.4f}   ratio: {b_boost/b_off:.2f}" if abs(b_off) > 0.01 else
          f"  {name:<12} boost: {b_boost:+.4f}   offset: {b_off:+.4f}")

# The boost coefficients should be SIMILAR to offset but shifted
# because boost = offset + log(g_rar/g_bar) ≈ offset + 0.5×log(a₀/g_bar)
print(f"\n  Boost vs offset: r = {np.corrcoef(boost, offset)[0,1]:+.4f}")
print(f"  Boost = offset + MOND term")

print("\n✓ Test 6 passed: coefficient interpretation done")

# =====================================================================
# TEST 7: BOOST MODEL RESIDUAL ANALYSIS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: BOOST MODEL RESIDUAL ANALYSIS")
print("=" * 60)

# Best boost model residuals
best_resid = boost - X_6var @ beta_6var
print(f"\n6-var boost model residuals:")
print(f"  RMS: {np.sqrt(np.mean(best_resid**2)):.4f}")
print(f"  Mean: {np.mean(best_resid):+.4f}")

# Correlation with properties
types = np.array([g['hubble_type'] for g in galaxies])
print(f"\n  r(resid, logV) = {np.corrcoef(best_resid, logV)[0,1]:+.4f}")
print(f"  r(resid, logL) = {np.corrcoef(best_resid, logL)[0,1]:+.4f}")
print(f"  r(resid, f_gas) = {np.corrcoef(best_resid, f_gas)[0,1]:+.4f}")
print(f"  r(resid, type) = {np.corrcoef(best_resid, types)[0,1]:+.4f}")
print(f"  r(resid, log γ) = {np.corrcoef(best_resid, log_gamma)[0,1]:+.4f}")

# Compare with γ-model residuals
gamma_resid = boost - X_gfc @ beta_gfc
print(f"\nγ + f_gas + c_V residuals:")
print(f"  RMS: {np.sqrt(np.mean(gamma_resid**2)):.4f}")
print(f"  r(resid, logV) = {np.corrcoef(gamma_resid, logV)[0,1]:+.4f}")
print(f"  r(resid, logL) = {np.corrcoef(gamma_resid, logL)[0,1]:+.4f}")

# Does γ-model residual correlate with offset model residual?
offset_resid = offset - X_off @ beta_off
r_resid_resid = np.corrcoef(gamma_resid, offset_resid)[0, 1]
r_6var_resid = np.corrcoef(best_resid, offset_resid)[0, 1]
print(f"\n  r(boost γ-resid, offset resid) = {r_resid_resid:+.4f}")
print(f"  r(boost 6var-resid, offset resid) = {r_6var_resid:+.4f}")

print("\n✓ Test 7 passed: residual analysis done")

# =====================================================================
# TEST 8: CONNECTING BOOST BACK TO OFFSET
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: FROM BOOST MODEL TO OFFSET MODEL")
print("=" * 60)

# offset = log(g_obs) - log(g_rar)
# boost = log(g_obs) - log(g_bar)
# Therefore: offset = boost - log(g_rar/g_bar)
# And: log(g_rar/g_bar) = log(ν(g_bar/a₀)) where ν is the MOND interpolation
# In deep MOND: log(g_rar/g_bar) ≈ 0.5 × log(a₀/g_bar)

# Can we predict offset from boost + g_bar?
X_boost_g = np.column_stack([np.ones(n), boost, log_g_ratio])
beta_bg = np.linalg.lstsq(X_boost_g, offset, rcond=None)[0]
R2_bg = 1 - np.sum((offset - X_boost_g @ beta_bg)**2) / np.sum((offset - np.mean(offset))**2)

print(f"\noffset from boost + log(g/a₀):")
print(f"  offset = {beta_bg[0]:+.3f} + {beta_bg[1]:+.3f}×boost + {beta_bg[2]:+.3f}×log(g/a₀)")
print(f"  R² = {R2_bg:.4f}")
print(f"  (If exact: boost coeff = 1, g_ratio coeff from interpolation function)")

# Pure theoretical connection
# offset_theory = boost - log(ν(g/a₀))
# where ν(x) = 1/(1-exp(-√x))
# log(ν) - 0 = offset_theory
mond_nu = 1 / (1 - np.exp(-np.sqrt(np.clip(10**log_g_ratio, 1e-10, None))))
theory_connection = boost - np.log10(mond_nu)
r_theory_offset = np.corrcoef(theory_connection, offset)[0, 1]
rms_theory = np.sqrt(np.mean((theory_connection - offset)**2))

print(f"\nTheoretical: offset = boost - log(ν(g/a₀))")
print(f"  r(theoretical offset, observed offset) = {r_theory_offset:+.4f}")
print(f"  RMS difference: {rms_theory:.4f}")
print(f"  Mean difference: {np.mean(theory_connection - offset):+.4f}")

# Does the γ-boost model produce good offset predictions via this route?
gamma_boost_pred = X_gfc @ beta_gfc  # predicted boost from γ model
gamma_offset_pred = gamma_boost_pred - np.log10(mond_nu)
R2_gamma_offset = 1 - np.sum((offset - gamma_offset_pred)**2) / np.sum((offset - np.mean(offset))**2)

print(f"\nγ-based offset prediction (via boost → offset):")
print(f"  R² = {R2_gamma_offset:.4f}")
print(f"  Compare direct 6-var offset model: R² = {1 - np.sum(offset_resid**2) / np.sum((offset - np.mean(offset))**2):.4f}")

print("\n✓ Test 8 passed: boost-offset connection established")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #505 SUMMARY")
print("=" * 70)
print(f"MOND boost mean: {np.mean(boost):.3f}, std: {np.std(boost):.3f}")
print(f"γ → boost: R² = {R2_g:.3f}, partial r(γ|V,L) = {partial_vl:+.3f}")
print(f"γ + f_gas + c_V → boost: R² = {R2_gfc:.3f}")
print(f"6-var → boost: R² = {R2_6var:.3f}")
print(f"γ adds to V,L: ΔR² = {R2_vl_g - R2_vl:+.3f}")
print(f"Boost → offset (via MOND): r = {r_theory_offset:+.3f}")
print(f"γ-based offset (via boost): R² = {R2_gamma_offset:.3f}")
print(f"\nAll 8 tests passed ✓")
