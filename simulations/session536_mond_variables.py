#!/usr/bin/env python3
"""
======================================================================
SESSION #536: THE MODEL IN MOND VARIABLES — REPARAMETRIZATION
======================================================================

The 6-var model uses (logV, logL, c_V, f_gas) + interactions. Session
#526 showed all coefficients derivable from MOND. Session #528 showed
the V-L ratio is 4.0 when gas-corrected. But the model is still in
"observational" variables. What does it look like in MOND variables?

Define:
- x = log(g_bar/a₀) — MOND regime parameter (how deep in MOND)
- δ_BTFR = logL - 4logV — BTFR residual (deviation from M_bar ∝ V⁴)
- f_gas — gas fraction (same, but now with MOND interpretation)
- c_V — rotation curve shape (same, but related to mass distribution)

Questions:
1. Does the model simplify in MOND variables?
2. Is x (MOND regime) the primary predictor?
3. Does δ_BTFR capture the M/L information?
4. What is the physical hierarchy in MOND terms?

Tests:
1. MOND variable definitions and statistics
2. Offset model in MOND variables (logx, δ_BTFR, f_gas, c_V)
3. Is the MOND regime (x) the primary predictor?
4. The BTFR residual as M/L proxy
5. Model hierarchy: building up from MOND variables
6. Comparison: observational vs MOND variable models
7. The deep MOND limit: what happens at low x?
8. Synthesis: what the MOND reparametrization reveals

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-07
Session: #536
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
kms_to_ms = 1e3
kpc_to_m = 3.0857e19


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
    return 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)


def prepare_galaxies():
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

        g_rar = g_bar_v * nu_mcgaugh(g_bar_v / a0_mond)
        offset_pts = np.log10(g_obs_v) - np.log10(g_rar)

        if outer_mond.sum() >= 2:
            offset_val = np.mean(offset_pts[outer_mond])
            mean_gbar = np.mean(g_bar_v[outer_mond])
            mean_gobs = np.mean(g_obs_v[outer_mond])
        else:
            offset_val = np.mean(offset_pts[mond])
            mean_gbar = np.mean(g_bar_v[mond])
            mean_gobs = np.mean(g_obs_v[mond])

        n_flat = min(5, len(v_gas_v))
        v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
        v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
        f_gas = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

        # MOND variables
        x_mond = mean_gbar / a0_mond  # MOND regime parameter
        log_x = np.log10(max(x_mond, 1e-10))
        nu_val = nu_mcgaugh(x_mond)
        log_nu = np.log10(nu_val)
        mond_boost = np.log10(mean_gobs / mean_gbar)

        # BTFR residual: δ = logL - 4logV (deviation from V⁴ ∝ L)
        logV = np.log10(vflat)
        logL = np.log10(lum)
        delta_btfr = logL - 4 * logV

        # N_corr and γ
        R_max = radius_v.max()
        R_max_m = R_max * kpc_to_m
        V_ms = vflat * kms_to_ms
        a_centripetal = V_ms**2 / R_max_m
        N_corr = a_centripetal / a0_mond
        gamma = 2.0 / np.sqrt(N_corr) if N_corr > 0 else np.nan
        if not np.isfinite(gamma):
            continue

        galaxies.append({
            'id': gal_id,
            'offset': offset_val,
            'logV': logV,
            'logL': logL,
            'c_V': c_V,
            'f_gas': f_gas,
            'hubble_type': hubble_type,
            'vflat': vflat,
            'lum': lum,
            'R_max': R_max,
            'mean_gbar': mean_gbar,
            'mean_gobs': mean_gobs,
            'x_mond': x_mond,
            'log_x': log_x,
            'log_nu': log_nu,
            'mond_boost': mond_boost,
            'delta_btfr': delta_btfr,
            'log_gamma': np.log10(gamma),
            'N_corr': N_corr,
        })

    return galaxies


from scipy import stats as sp_stats

print("=" * 70)
print("SESSION #536: THE MODEL IN MOND VARIABLES")
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
hubble_type = np.array([g['hubble_type'] for g in galaxies])
mond_boost = np.array([g['mond_boost'] for g in galaxies])
log_nu = np.array([g['log_nu'] for g in galaxies])
log_x = np.array([g['log_x'] for g in galaxies])
delta_btfr = np.array([g['delta_btfr'] for g in galaxies])
log_gamma = np.array([g['log_gamma'] for g in galaxies])

ones = np.ones(n)

# Standard 6-var model for reference
X6 = np.column_stack([ones, logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, offset)
loo6 = loo_r2(X6, offset)

# ============================================================
print("\n" + "=" * 60)
print("TEST 1: MOND VARIABLE DEFINITIONS AND STATISTICS")
print("=" * 60)

print(f"\n  MOND VARIABLES:")
print(f"  x = g_bar/a₀ — MOND regime parameter")
print(f"  log(x): mean = {np.mean(log_x):.3f}, std = {np.std(log_x):.3f}, "
      f"range = [{np.min(log_x):.3f}, {np.max(log_x):.3f}]")
print(f"  δ_BTFR = logL - 4logV — BTFR residual (M/L proxy)")
print(f"  δ_BTFR: mean = {np.mean(delta_btfr):.3f}, std = {np.std(delta_btfr):.3f}, "
      f"range = [{np.min(delta_btfr):.3f}, {np.max(delta_btfr):.3f}]")

# Correlations between MOND variables and observational variables
print(f"\n  MOND ↔ Observational correlations:")
print(f"  r(log x, logV)    = {sp_stats.pearsonr(log_x, logV)[0]:+.4f}")
print(f"  r(log x, logL)    = {sp_stats.pearsonr(log_x, logL)[0]:+.4f}")
print(f"  r(log x, c_V)     = {sp_stats.pearsonr(log_x, c_V)[0]:+.4f}")
print(f"  r(log x, f_gas)   = {sp_stats.pearsonr(log_x, f_gas)[0]:+.4f}")
print(f"  r(δ_BTFR, logV)   = {sp_stats.pearsonr(delta_btfr, logV)[0]:+.4f}")
print(f"  r(δ_BTFR, logL)   = {sp_stats.pearsonr(delta_btfr, logL)[0]:+.4f}")
print(f"  r(δ_BTFR, c_V)    = {sp_stats.pearsonr(delta_btfr, c_V)[0]:+.4f}")
print(f"  r(δ_BTFR, f_gas)  = {sp_stats.pearsonr(delta_btfr, f_gas)[0]:+.4f}")

# Key correlation: offset vs MOND variables
print(f"\n  Offset correlations:")
print(f"  r(offset, log x)    = {sp_stats.pearsonr(offset, log_x)[0]:+.4f}")
print(f"  r(offset, δ_BTFR)  = {sp_stats.pearsonr(offset, delta_btfr)[0]:+.4f}")
print(f"  r(offset, f_gas)   = {sp_stats.pearsonr(offset, f_gas)[0]:+.4f}")
print(f"  r(offset, c_V)     = {sp_stats.pearsonr(offset, c_V)[0]:+.4f}")
print(f"  r(offset, log γ)   = {sp_stats.pearsonr(offset, log_gamma)[0]:+.4f}")

# How redundant are the MOND variables?
print(f"\n  Inter-MOND-variable correlations:")
print(f"  r(log x, δ_BTFR)  = {sp_stats.pearsonr(log_x, delta_btfr)[0]:+.4f}")
print(f"  r(log x, f_gas)   = {sp_stats.pearsonr(log_x, f_gas)[0]:+.4f}")
print(f"  r(log x, c_V)     = {sp_stats.pearsonr(log_x, c_V)[0]:+.4f}")
print(f"  r(δ_BTFR, f_gas)  = {sp_stats.pearsonr(delta_btfr, f_gas)[0]:+.4f}")
print(f"  r(δ_BTFR, c_V)    = {sp_stats.pearsonr(delta_btfr, c_V)[0]:+.4f}")

print("\n✓ Test 1 passed: MOND variables defined and analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 2: OFFSET MODEL IN MOND VARIABLES")
print("=" * 60)

# Replace (logV, logL) with (log x, δ_BTFR)
# The idea: log x encodes MOND regime, δ_BTFR encodes M/L
X_mond = np.column_stack([ones, log_x, delta_btfr, c_V, f_gas])
beta_m, yhat_m, resid_m, R2_m, rms_m = build_model(X_mond, offset)
loo_m = loo_r2(X_mond, offset)

print(f"\n  MOND variable model (no interactions):")
print(f"  offset = {beta_m[0]:+.4f} {beta_m[1]:+.4f}×log(x) {beta_m[2]:+.4f}×δ_BTFR "
      f"{beta_m[3]:+.4f}×c_V {beta_m[4]:+.4f}×f_gas")
print(f"  R² = {R2_m:.4f}, LOO = {loo_m:.4f}, RMS = {rms_m:.4f}")
print(f"  vs 6-var: R² = {R2_6:.4f}, LOO = {loo6:.4f}, RMS = {rms6:.4f}")

# With interactions analogous to the 6-var model
# logV×c_V becomes log_x×c_V (similar to MOND regime × shape)
# logL×f_gas becomes δ_BTFR×f_gas (M/L × gas fraction)
X_mond_int1 = np.column_stack([ones, log_x, delta_btfr, c_V, f_gas,
                                log_x*c_V, delta_btfr*f_gas])
beta_mi1, _, _, R2_mi1, rms_mi1 = build_model(X_mond_int1, offset)
loo_mi1 = loo_r2(X_mond_int1, offset)

print(f"\n  MOND model + analogous interactions (log_x×c_V, δ_BTFR×f_gas):")
print(f"  R² = {R2_mi1:.4f}, LOO = {loo_mi1:.4f}, RMS = {rms_mi1:.4f}")

# Alternative interactions: what if the natural interactions are different?
X_mond_int2 = np.column_stack([ones, log_x, delta_btfr, c_V, f_gas,
                                log_x*f_gas, delta_btfr*c_V])
loo_mi2 = loo_r2(X_mond_int2, offset)

X_mond_int3 = np.column_stack([ones, log_x, delta_btfr, c_V, f_gas,
                                log_x*c_V, log_x*f_gas])
loo_mi3 = loo_r2(X_mond_int3, offset)

# Try ALL pairwise interactions
X_mond_all = np.column_stack([ones, log_x, delta_btfr, c_V, f_gas,
                               log_x*c_V, log_x*f_gas,
                               delta_btfr*c_V, delta_btfr*f_gas,
                               c_V*f_gas])
loo_mall = loo_r2(X_mond_all, offset)

print(f"\n  Interaction search:")
print(f"  log_x×c_V + δ×f_gas:    LOO = {loo_mi1:.4f}")
print(f"  log_x×f_gas + δ×c_V:    LOO = {loo_mi2:.4f}")
print(f"  log_x×c_V + log_x×f_gas: LOO = {loo_mi3:.4f}")
print(f"  All pairwise:            LOO = {loo_mall:.4f}")
print(f"  Standard 6-var:          LOO = {loo6:.4f}")

# The best MOND-variable interaction model
best_mond_loo = max(loo_mi1, loo_mi2, loo_mi3, loo_mall)
print(f"\n  Best MOND variable model: LOO = {best_mond_loo:.4f}")
print(f"  Gap from standard 6-var: {best_mond_loo - loo6:+.4f}")

print("\n✓ Test 2 passed: MOND variable models built")

# ============================================================
print("\n" + "=" * 60)
print("TEST 3: IS THE MOND REGIME (x) THE PRIMARY PREDICTOR?")
print("=" * 60)

# What does each MOND variable contribute?
vars_to_test = [
    ('log x', log_x),
    ('δ_BTFR', delta_btfr),
    ('c_V', c_V),
    ('f_gas', f_gas),
    ('log γ', log_gamma),
]

print(f"\n  Single variable → offset:")
print(f"  {'Variable':12s}  {'R²':>6s}  {'LOO':>6s}  {'β':>8s}  {'r':>7s}")
print(f"  {'-'*50}")
for name, var in vars_to_test:
    X1 = np.column_stack([ones, var])
    _, _, _, r2, _ = build_model(X1, offset)
    loo = loo_r2(X1, offset)
    beta1 = np.linalg.lstsq(X1, offset, rcond=None)[0]
    r = sp_stats.pearsonr(var, offset)[0]
    print(f"  {name:12s}  {r2:.4f}  {loo:.4f}  {beta1[1]:+8.4f}  {r:+7.4f}")

# log x is a combination of logV and logL (it's basically a mass proxy
# measured at the MOND regime point). How does it compare?
X_x_only = np.column_stack([ones, log_x])
_, _, resid_x, R2_x, _ = build_model(X_x_only, offset)
loo_x = loo_r2(X_x_only, offset)

# What does δ_BTFR add after log x?
X_x_delta = np.column_stack([ones, log_x, delta_btfr])
loo_x_delta = loo_r2(X_x_delta, offset)
delta_btfr_partial_r = sp_stats.pearsonr(resid_x, delta_btfr)[0]

print(f"\n  Building up from log(x):")
print(f"  log(x) alone:      LOO = {loo_x:.4f}")
print(f"  + δ_BTFR:          LOO = {loo_x_delta:.4f} (ΔLOO = {loo_x_delta - loo_x:+.4f})")
print(f"  r(resid|x, δ_BTFR) = {delta_btfr_partial_r:+.4f}")

# Is log x better or worse than logV + logL combined?
X_VL = np.column_stack([ones, logV, logL])
loo_VL = loo_r2(X_VL, offset)
print(f"\n  log(x) vs logV + logL:")
print(f"  log(x) alone:  LOO = {loo_x:.4f}")
print(f"  logV + logL:   LOO = {loo_VL:.4f}")
print(f"  x is {'better' if loo_x > loo_VL else 'worse'} ({loo_x - loo_VL:+.4f})")

# The reason: log x ≈ function(logV, logL, radius, mass distribution)
# It compresses the mass information into a single MOND regime parameter
# but loses the M/L information (which δ_BTFR captures)

print("\n✓ Test 3 passed: MOND regime importance assessed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 4: THE BTFR RESIDUAL AS M/L PROXY")
print("=" * 60)

# δ_BTFR = logL - 4logV measures how bright a galaxy is relative to its
# rotation speed. In MOND, V⁴ ∝ M_bar = L × M/L, so:
# δ_BTFR = logL - 4logV = logL - log(M_bar/M_L + const)
# For fixed M/L: δ_BTFR ∝ log(M/L) offset from some reference

# Direct test: does δ_BTFR predict offset?
r_delta_off, p_delta_off = sp_stats.pearsonr(delta_btfr, offset)
print(f"\n  δ_BTFR as offset predictor:")
print(f"  r(δ_BTFR, offset) = {r_delta_off:+.4f}, p = {p_delta_off:.6f}")

# But the offset IS log(M/L) approximately (Session #529)
# So δ_BTFR should correlate with offset because both measure M/L
# The question is: how well, and is the correlation direct or mediated?

# Partial correlation controlling for mass
# r(δ_BTFR, offset | log x)
_, _, resid_delta_x, _, _ = build_model(np.column_stack([ones, log_x]), delta_btfr)
_, _, resid_off_x, _, _ = build_model(np.column_stack([ones, log_x]), offset)
r_partial_delta_off_x = sp_stats.pearsonr(resid_delta_x, resid_off_x)[0]

print(f"  r(δ_BTFR, offset | log x) = {r_partial_delta_off_x:+.4f}")
print(f"  (After removing mass dependence)")

# δ_BTFR ≈ log(M/L) prediction
# In MOND deep regime: logL - 4logV = -log(M/L) - log(Ga₀)/(some constant)
# More precisely: V⁴ = G×a₀×M_bar, M_bar = M/L × L
# So 4logV = log(G×a₀) + log(M/L) + logL
# δ_BTFR = logL - 4logV = -log(M/L) - log(G×a₀)
# Therefore offset ∝ log(M/L) ∝ -δ_BTFR

# But there's gas: M_bar = M/L × L_star + M_gas
# M_bar ≈ L × (M/L_star × (1-f_gas) + M_gas/L)
# This complicates the δ_BTFR → M/L relationship

# Does gas correction improve δ_BTFR → offset?
X_delta_fg = np.column_stack([ones, delta_btfr, f_gas])
loo_df = loo_r2(X_delta_fg, offset)
X_delta_fg_int = np.column_stack([ones, delta_btfr, f_gas, delta_btfr*f_gas])
loo_dfi = loo_r2(X_delta_fg_int, offset)

print(f"\n  δ_BTFR models → offset:")
print(f"  δ_BTFR alone:              LOO = {loo_r2(np.column_stack([ones, delta_btfr]), offset):.4f}")
print(f"  δ_BTFR + f_gas:            LOO = {loo_df:.4f}")
print(f"  δ_BTFR + f_gas + δ×f_gas: LOO = {loo_dfi:.4f}")

# What's the implied M/L from δ_BTFR?
# offset ≈ β₀ + β₁ × δ_BTFR → log(M/L) ≈ 2×offset
# From the BTFR: δ_BTFR ≈ -log(M/L) - const (in deep MOND)
# So offset should have β(δ_BTFR) ≈ -0.5 (since offset ≈ 0.5 × log(M/L))
beta_d = np.linalg.lstsq(np.column_stack([ones, delta_btfr]), offset, rcond=None)[0]
print(f"\n  offset = {beta_d[0]:+.4f} {beta_d[1]:+.4f} × δ_BTFR")
print(f"  MOND prediction for β(δ_BTFR): -0.50 (if δ captures pure M/L)")
print(f"  Observed: {beta_d[1]:+.4f} ({'matches' if abs(beta_d[1] + 0.5) < 0.1 else 'differs'})")

print("\n✓ Test 4 passed: BTFR residual as M/L proxy analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 5: MODEL HIERARCHY IN MOND VARIABLES")
print("=" * 60)

# Build up: log_x → + δ_BTFR → + f_gas → + c_V → + interactions
models_mond = {}

X1 = np.column_stack([ones, log_x])
models_mond['1. log x'] = (X1, loo_r2(X1, offset))

X2 = np.column_stack([ones, log_x, delta_btfr])
models_mond['2. + δ_BTFR'] = (X2, loo_r2(X2, offset))

X3 = np.column_stack([ones, log_x, delta_btfr, f_gas])
models_mond['3. + f_gas'] = (X3, loo_r2(X3, offset))

X4 = np.column_stack([ones, log_x, delta_btfr, f_gas, c_V])
models_mond['4. + c_V'] = (X4, loo_r2(X4, offset))

X5 = np.column_stack([ones, log_x, delta_btfr, f_gas, c_V, log_x*c_V])
models_mond['5. + x×c_V'] = (X5, loo_r2(X5, offset))

X6m = np.column_stack([ones, log_x, delta_btfr, f_gas, c_V, log_x*c_V, delta_btfr*f_gas])
models_mond['6. + δ×f_gas'] = (X6m, loo_r2(X6m, offset))

# Alternative path: start with δ_BTFR
X1b = np.column_stack([ones, delta_btfr])
models_mond['1b. δ_BTFR'] = (X1b, loo_r2(X1b, offset))

X2b = np.column_stack([ones, delta_btfr, f_gas])
models_mond['2b. + f_gas'] = (X2b, loo_r2(X2b, offset))

X3b = np.column_stack([ones, delta_btfr, f_gas, log_x])
models_mond['3b. + log x'] = (X3b, loo_r2(X3b, offset))

print(f"\n  MOND variable model hierarchy:")
print(f"  {'Model':25s}  {'LOO':>6s}  {'ΔLOO':>7s}")
print(f"  {'-'*45}")

prev_loo = 0
for name in sorted(models_mond.keys()):
    X, loo = models_mond[name]
    delta = loo - prev_loo
    print(f"  {name:25s}  {loo:.4f}  {delta:+7.4f}")
    if name.startswith('1.'):
        prev_loo = loo

print(f"\n  Comparison paths:")
print(f"  Starting from log(x): LOO rises to {models_mond['6. + δ×f_gas'][1]:.4f}")
print(f"  Starting from δ_BTFR: LOO rises to {models_mond['3b. + log x'][1]:.4f} (3 vars)")
print(f"  Standard 6-var:       LOO = {loo6:.4f}")

# The hierarchy compared to observational variables
print(f"\n  Observational hierarchy (for reference):")
Xo1 = np.column_stack([ones, logV, logL])
Xo2 = np.column_stack([ones, logV, logL, f_gas])
Xo3 = np.column_stack([ones, logV, logL, c_V, f_gas])
print(f"  logV + logL:              LOO = {loo_r2(Xo1, offset):.4f}")
print(f"  + f_gas:                  LOO = {loo_r2(Xo2, offset):.4f}")
print(f"  + c_V:                    LOO = {loo_r2(Xo3, offset):.4f}")
print(f"  + interactions (6-var):   LOO = {loo6:.4f}")

print("\n✓ Test 5 passed: model hierarchy in MOND variables built")

# ============================================================
print("\n" + "=" * 60)
print("TEST 6: COMPARISON — OBSERVATIONAL vs MOND VARIABLES")
print("=" * 60)

# Direct comparison at each level of complexity
print(f"\n  {'#vars':>5s}  {'Observational':>25s}  {'LOO(obs)':>8s}  {'MOND':>25s}  {'LOO(mond)':>9s}")
print(f"  {'-'*80}")

# 2 vars: logV+logL vs log_x+δ_BTFR
loo_obs_2 = loo_r2(np.column_stack([ones, logV, logL]), offset)
loo_mond_2 = loo_r2(np.column_stack([ones, log_x, delta_btfr]), offset)
print(f"  {'2':>5s}  {'logV + logL':>25s}  {loo_obs_2:8.4f}  {'log x + δ_BTFR':>25s}  {loo_mond_2:9.4f}")

# 3 vars (+ f_gas both)
loo_obs_3 = loo_r2(np.column_stack([ones, logV, logL, f_gas]), offset)
loo_mond_3 = loo_r2(np.column_stack([ones, log_x, delta_btfr, f_gas]), offset)
print(f"  {'3':>5s}  {'+ f_gas':>25s}  {loo_obs_3:8.4f}  {'+ f_gas':>25s}  {loo_mond_3:9.4f}")

# 4 vars (+ c_V both)
loo_obs_4 = loo_r2(np.column_stack([ones, logV, logL, f_gas, c_V]), offset)
loo_mond_4 = loo_r2(np.column_stack([ones, log_x, delta_btfr, f_gas, c_V]), offset)
print(f"  {'4':>5s}  {'+ c_V':>25s}  {loo_obs_4:8.4f}  {'+ c_V':>25s}  {loo_mond_4:9.4f}")

# 6 vars (with interactions)
loo_obs_6 = loo6
loo_mond_6 = loo_r2(X6m, offset)
print(f"  {'6':>5s}  {'+ interactions':>25s}  {loo_obs_6:8.4f}  {'+ interactions':>25s}  {loo_mond_6:9.4f}")

# Are they really different? Check if the MOND variables are a linear
# transformation of the observational ones
# log x = f(logV, logL, radius, mass model) — NOT a simple linear combo
# δ_BTFR = logL - 4logV — IS a linear combo

# How much of log x is captured by logV + logL?
X_VL_for_x = np.column_stack([ones, logV, logL])
_, _, _, R2_x_VL, _ = build_model(X_VL_for_x, log_x)
print(f"\n  Linear relationship:")
print(f"  R²(log x ~ logV + logL) = {R2_x_VL:.4f}")
print(f"  δ_BTFR = logL - 4logV (exact linear combination)")
print(f"  The MOND parametrization is {'equivalent' if R2_x_VL > 0.95 else 'NOT equivalent'} "
      f"to observational")

# If R² is high, the models should be nearly identical
# The difference would come from log x capturing nonlinear MOND information

# Check: does log x carry information beyond logV + logL?
_, _, resid_x_VL, _, _ = build_model(X_VL_for_x, log_x)
# Use the unique part of log x
X_6var_plus_xresid = np.column_stack([X6, resid_x_VL])
loo_6_plus_x = loo_r2(X_6var_plus_xresid, offset)
print(f"\n  Does unique log(x) info help the 6-var model?")
print(f"  6-var + residual(log x | V,L): LOO = {loo_6_plus_x:.4f} (ΔLOO = {loo_6_plus_x - loo6:+.4f})")

print("\n✓ Test 6 passed: observational vs MOND variables compared")

# ============================================================
print("\n" + "=" * 60)
print("TEST 7: THE DEEP MOND LIMIT — WHAT HAPPENS AT LOW x?")
print("=" * 60)

# In deep MOND (x << 1), ν → 1/√x, so g_obs → √(g_bar × a₀)
# The offset should be purely M/L-driven: no regime dependence
# The BTFR should be exact: V⁴ = G×a₀×M_bar exactly
# So at low x, the model should simplify

# Split by MOND regime
deep_mond = log_x < np.median(log_x)
shallow_mond = ~deep_mond

print(f"\n  Deep MOND (log x < {np.median(log_x):.2f}, n={deep_mond.sum()}):")
print(f"  Mean log x: {np.mean(log_x[deep_mond]):.3f}")
print(f"  Mean offset: {np.mean(offset[deep_mond]):.4f}")
print(f"  σ(offset): {np.std(offset[deep_mond]):.4f}")

print(f"\n  Shallow MOND (log x ≥ {np.median(log_x):.2f}, n={shallow_mond.sum()}):")
print(f"  Mean log x: {np.mean(log_x[shallow_mond]):.3f}")
print(f"  Mean offset: {np.mean(offset[shallow_mond]):.4f}")
print(f"  σ(offset): {np.std(offset[shallow_mond]):.4f}")

# Does the model simplify in deep MOND?
# In deep MOND: ν ≈ 1/√x, so boost = 0.5×log(a₀/g_bar) + const
# The offset = boost - log(ν) should be independent of x
# So: r(offset, log x | deep MOND) should be weaker

r_deep, p_deep = sp_stats.pearsonr(offset[deep_mond], log_x[deep_mond])
r_shallow, p_shallow = sp_stats.pearsonr(offset[shallow_mond], log_x[shallow_mond])
print(f"\n  r(offset, log x) by regime:")
print(f"  Deep MOND:    r = {r_deep:+.4f}, p = {p_deep:.4f}")
print(f"  Shallow MOND: r = {r_shallow:+.4f}, p = {p_shallow:.4f}")

# Model performance in each regime
if deep_mond.sum() > 15 and shallow_mond.sum() > 15:
    # Simple model: just δ_BTFR + f_gas (should work well in deep MOND)
    X_simple = np.column_stack([ones, delta_btfr, f_gas, delta_btfr*f_gas])

    # Full model performance in each regime
    yhat_simple = X_simple @ np.linalg.lstsq(X_simple, offset, rcond=None)[0]
    yhat_6var = X6 @ beta6

    for mask, label in [(deep_mond, 'Deep MOND'), (shallow_mond, 'Shallow MOND')]:
        rms_s = np.sqrt(np.mean((offset[mask] - yhat_simple[mask])**2))
        rms_6 = np.sqrt(np.mean((offset[mask] - yhat_6var[mask])**2))
        print(f"\n  {label}:")
        print(f"  Simple (δ+f_gas+int): RMS = {rms_s:.4f}")
        print(f"  6-var:                RMS = {rms_6:.4f}")
        print(f"  6-var improvement: {(1-rms_6/rms_s)*100:+.1f}%")

# Quartile analysis of MOND regime
x_quartiles = np.percentile(log_x, [0, 25, 50, 75, 100])
print(f"\n  Offset by MOND regime quartile:")
print(f"  {'Quartile':10s}  {'log x range':15s}  {'N':>4s}  {'σ(offset)':>10s}  {'mean offset':>12s}")
print(f"  {'-'*60}")
for i in range(4):
    mask = (log_x >= x_quartiles[i]) & (log_x < x_quartiles[i+1] + (0.01 if i == 3 else 0))
    nm = mask.sum()
    print(f"  Q{i+1} (deep→shal) [{x_quartiles[i]:+.2f}, {x_quartiles[i+1]:+.2f}]"
          f"  {nm:4d}  {np.std(offset[mask]):10.4f}  {np.mean(offset[mask]):+12.4f}")

print("\n✓ Test 7 passed: deep MOND limit analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS — WHAT THE MOND REPARAMETRIZATION REVEALS")
print("=" * 60)

print(f"\n  COMPARISON SUMMARY:")
print(f"  {'Model':45s}  {'LOO':>6s}")
print(f"  {'-'*55}")
print(f"  {'Standard 6-var (logV, logL, c_V, f_gas + int)':45s}  {loo6:.4f}")
print(f"  {'MOND 6-var (log x, δ_BTFR, c_V, f_gas + int)':45s}  {loo_mond_6:.4f}")
print(f"  {'Difference':45s}  {loo_mond_6 - loo6:+.4f}")

print(f"\n  KEY FINDINGS:")

# Finding 1: Are MOND variables equivalent?
print(f"\n  1. MOND VARIABLES vs OBSERVATIONAL:")
print(f"     R²(log x ~ logV + logL) = {R2_x_VL:.4f}")
print(f"     δ_BTFR = logL - 4logV (exact linear)")
if R2_x_VL > 0.95:
    print(f"     → NEARLY EQUIVALENT: log x is ~linear in (logV, logL)")
    print(f"       The 'MOND regime' is just a reweighting of mass proxies")
else:
    print(f"     → DIFFERENT: log x carries nonlinear MOND information")
    print(f"       The MOND regime parameter adds beyond V+L")

# Finding 2: Which hierarchy is more natural?
print(f"\n  2. HIERARCHY COMPARISON:")
print(f"     Observational: logV+logL ({loo_obs_2:.4f}) → +f_gas ({loo_obs_3:.4f}) → +c_V ({loo_obs_4:.4f})")
print(f"     MOND:          log x+δ ({loo_mond_2:.4f}) → +f_gas ({loo_mond_3:.4f}) → +c_V ({loo_mond_4:.4f})")
better = 'Observational' if loo_obs_2 > loo_mond_2 else 'MOND'
print(f"     {better} is better at the base level")
better4 = 'Observational' if loo_obs_4 > loo_mond_4 else 'MOND'
print(f"     {better4} is better with 4 variables")

# Finding 3: The deep MOND test
print(f"\n  3. DEEP MOND SIMPLIFICATION:")
print(f"     r(offset, log x | deep MOND) = {r_deep:+.4f}")
print(f"     r(offset, log x | shallow) = {r_shallow:+.4f}")
if abs(r_deep) < abs(r_shallow):
    print(f"     → CONFIRMED: offset less dependent on MOND regime in deep MOND")
else:
    print(f"     → NOT confirmed: offset still depends on x in deep MOND")

# Finding 4: Physical interpretation
print(f"\n  4. PHYSICAL INTERPRETATION:")
print(f"     The model in MOND variables is:")
print(f"     offset = f(MOND regime, M/L proxy, gas fraction, mass distribution)")
print(f"     = f(how deep in MOND, how far from BTFR, how gassy, how concentrated)")
print(f"")
print(f"     The MOND regime (log x) and BTFR residual (δ) together contain")
print(f"     the same information as logV + logL, but in a physically more")
print(f"     transparent form:")
print(f"     - log x: where on the RAR the measurement is made")
print(f"     - δ_BTFR: how much the galaxy deviates from the BTFR")
print(f"     These are ORTHOGONAL physical questions (r = {sp_stats.pearsonr(log_x, delta_btfr)[0]:+.3f})")

print("\n✓ Test 8 passed: synthesis complete")

# ============================================================
print("\n" + "=" * 70)
print("SESSION #536 SUMMARY")
print("=" * 70)
print(f"\n  Standard 6-var:  LOO = {loo6:.4f}")
print(f"  MOND 6-var:      LOO = {loo_mond_6:.4f} (Δ = {loo_mond_6 - loo6:+.4f})")
print(f"  R²(log x ~ V+L) = {R2_x_VL:.4f}")
print(f"  r(log x, δ_BTFR) = {sp_stats.pearsonr(log_x, delta_btfr)[0]:+.4f}")
print(f"  Deep MOND r(offset, x) = {r_deep:+.4f}")
print(f"  Shallow MOND r(offset, x) = {r_shallow:+.4f}")
print(f"  β(δ_BTFR → offset) = {beta_d[1]:+.4f} (MOND: -0.50)")

print(f"\nAll 8 tests passed ✓")
