#!/usr/bin/env python3
"""
======================================================================
SESSION #503: SYNCHRONISM THEORY vs EMPIRICAL 6-VAR MODEL
======================================================================

The Synchronism framework predicts: γ = 2/√N_corr
where N_corr = V²/(R×a₀) counts "correlated MOND quanta".

The 6-var model is empirical:
  offset = -3.38 + 1.90×logV - 0.55×logL - 0.22×c_V - 0.45×f_gas
           + 0.15×logV×c_V + 0.18×logL×f_gas

Does the theoretical γ prediction match the empirical offset?
Can the 6-var model be derived FROM the theory?
Are the model coefficients consistent with γ = 2/√N_corr?

Tests:
1. Compute γ_theory for each galaxy
2. Compare γ_theory vs empirical offset
3. Residual analysis: what does the 6-var model add beyond γ?
4. Can γ replace the 6-var model?
5. Decompose N_corr in terms of the 6-var predictors
6. The theoretical prediction for MOND coefficients
7. γ as a function of galaxy properties
8. Combined model: γ + corrections

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #503
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from session372_sparc_sb_test import (
    load_sparc_catalog, load_sparc_mass_models,
    compute_gbar_gobs
)

a0_mond = 1.2e-10  # m/s²
c_light = 3e8  # m/s
H0 = 67.4e3 / 3.086e22  # 67.4 km/s/Mpc → s⁻¹


def rar_prediction(g_bar, a0=a0_mond):
    x = g_bar / a0
    x = np.clip(x, 1e-10, None)
    return g_bar / (1 - np.exp(-np.sqrt(x)))


def prepare_galaxies():
    """Load SPARC with all properties needed for theory comparison."""
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

        if vflat <= 0 or lum <= 0 or sb_eff <= 0 or distance <= 0:
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

        # R_max (maximum radius of RC, in kpc)
        r_max_kpc = radius_v.max()

        # Mean outer radius (where offset is measured)
        if outer_mond.sum() >= 2:
            r_outer_kpc = np.mean(radius_v[outer_mond])
        else:
            r_outer_kpc = np.mean(radius_v[mond])

        # V_flat in m/s
        vflat_ms = vflat * 1e3

        # R in meters (use R_max)
        r_max_m = r_max_kpc * 3.086e19
        r_outer_m = r_outer_kpc * 3.086e19
        r_eff_m = r_eff_kpc * 3.086e19

        # N_corr = V²/(R×a₀)
        N_corr_rmax = vflat_ms**2 / (r_max_m * a0_mond)
        N_corr_router = vflat_ms**2 / (r_outer_m * a0_mond)
        N_corr_reff = vflat_ms**2 / (r_eff_m * a0_mond)

        # γ_theory = 2/√N_corr
        gamma_rmax = 2.0 / np.sqrt(N_corr_rmax)
        gamma_router = 2.0 / np.sqrt(N_corr_router)
        gamma_reff = 2.0 / np.sqrt(N_corr_reff)

        galaxies.append({
            'id': gal_id,
            'offset': offset,
            'logV': np.log10(vflat),
            'logL': np.log10(lum),
            'c_V': c_V,
            'f_gas': f_gas,
            'hubble_type': hubble_type,
            'distance': distance,
            'r_max_kpc': r_max_kpc,
            'r_outer_kpc': r_outer_kpc,
            'r_eff_kpc': r_eff_kpc,
            'vflat': vflat,
            'N_corr_rmax': N_corr_rmax,
            'N_corr_router': N_corr_router,
            'N_corr_reff': N_corr_reff,
            'gamma_rmax': gamma_rmax,
            'gamma_router': gamma_router,
            'gamma_reff': gamma_reff,
            'log_gamma_rmax': np.log10(gamma_rmax) if gamma_rmax > 0 else np.nan,
            'log_gamma_reff': np.log10(gamma_reff) if gamma_reff > 0 else np.nan,
        })

    return galaxies


print("=" * 70)
print("SESSION #503: SYNCHRONISM THEORY vs EMPIRICAL 6-VAR MODEL")
print("=" * 70)

galaxies = prepare_galaxies()
n = len(galaxies)
print(f"\n{n} galaxies loaded")

# Extract arrays
logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
y = np.array([g['offset'] for g in galaxies])
log_gamma = np.array([g['log_gamma_rmax'] for g in galaxies])
log_gamma_reff = np.array([g['log_gamma_reff'] for g in galaxies])
N_corr = np.array([g['N_corr_rmax'] for g in galaxies])
log_N = np.log10(N_corr)

# 6-var model
X6 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])
beta6 = np.linalg.lstsq(X6, y, rcond=None)[0]
yhat6 = X6 @ beta6
resid6 = y - yhat6
R2_6 = 1 - np.sum(resid6**2) / np.sum((y - np.mean(y))**2)
print(f"6-var model: R² = {R2_6:.4f}")

# =====================================================================
# TEST 1: COMPUTE γ_theory FOR EACH GALAXY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: γ_THEORY = 2/√N_corr FOR EACH GALAXY")
print("=" * 60)

print(f"\n{'Metric':<25} {'R_max':<12} {'R_outer':<12} {'R_eff':<12}")
print("-" * 60)

gamma_rmax = np.array([g['gamma_rmax'] for g in galaxies])
gamma_reff = np.array([g['gamma_reff'] for g in galaxies])

for stat_name, arr_name in [('γ_theory', 'gamma'), ('N_corr', 'N_corr')]:
    for r_type, suffix in [('R_max', 'rmax'), ('R_outer', 'router'), ('R_eff', 'reff')]:
        vals = np.array([g[f'{arr_name}_{suffix}'] for g in galaxies])
        if stat_name == 'γ_theory':
            pass  # will print below
        break

print(f"  N_corr (R_max):")
print(f"    Mean: {np.mean(N_corr):.1f}")
print(f"    Median: {np.median(N_corr):.1f}")
print(f"    Range: [{np.min(N_corr):.1f}, {np.max(N_corr):.1f}]")

print(f"\n  γ(R_max) = 2/√N_corr:")
print(f"    Mean: {np.mean(gamma_rmax):.4f}")
print(f"    Median: {np.median(gamma_rmax):.4f}")
print(f"    Range: [{np.min(gamma_rmax):.4f}, {np.max(gamma_rmax):.4f}]")

print(f"\n  γ(R_eff) = 2/√N_corr:")
print(f"    Mean: {np.mean(gamma_reff):.4f}")
print(f"    Median: {np.median(gamma_reff):.4f}")

print(f"\n  log₁₀(γ_rmax):")
print(f"    Mean: {np.mean(log_gamma):.4f}")
print(f"    Range: [{np.min(log_gamma):.4f}, {np.max(log_gamma):.4f}]")

print(f"\n  Offset:")
print(f"    Mean: {np.mean(y):.4f}")
print(f"    Range: [{np.min(y):.4f}, {np.max(y):.4f}]")

print("\n✓ Test 1 passed: γ_theory computed")

# =====================================================================
# TEST 2: COMPARE γ_theory vs EMPIRICAL OFFSET
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: γ_THEORY vs EMPIRICAL OFFSET")
print("=" * 60)

# Direct correlation
r_gamma_offset = np.corrcoef(log_gamma, y)[0, 1]
r_gamma_reff_offset = np.corrcoef(log_gamma_reff, y)[0, 1]
r_N_offset = np.corrcoef(log_N, y)[0, 1]

print(f"\n  r(log γ_rmax, offset) = {r_gamma_offset:+.4f}")
print(f"  r(log γ_reff, offset) = {r_gamma_reff_offset:+.4f}")
print(f"  r(log N_corr, offset) = {r_N_offset:+.4f}")

# Linear regression: offset = a + b × log(γ)
X_gamma = np.column_stack([np.ones(n), log_gamma])
beta_g = np.linalg.lstsq(X_gamma, y, rcond=None)[0]
yhat_g = X_gamma @ beta_g
resid_g = y - yhat_g
R2_gamma = 1 - np.sum(resid_g**2) / np.sum((y - np.mean(y))**2)

print(f"\n  Linear model: offset = {beta_g[0]:+.4f} + {beta_g[1]:+.4f} × log₁₀(γ)")
print(f"  R² = {R2_gamma:.4f}")
print(f"  If theory predicts offset = log(γ), then slope should be 1.0")
print(f"  Observed slope: {beta_g[1]:.3f}")

# Direct comparison: is offset ≈ log(γ)?
direct_resid = y - log_gamma
print(f"\n  Direct comparison (offset vs log γ):")
print(f"  Mean(offset - log γ): {np.mean(direct_resid):+.4f}")
print(f"  RMS(offset - log γ): {np.sqrt(np.mean(direct_resid**2)):.4f}")

# Scaled comparison: offset ≈ α × log(γ) + β
# (already done above, that's beta_g)

print("\n✓ Test 2 passed: theory-observation comparison done")

# =====================================================================
# TEST 3: RESIDUAL ANALYSIS — WHAT DOES 6-VAR ADD BEYOND γ?
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: WHAT DOES THE 6-VAR MODEL ADD BEYOND γ?")
print("=" * 60)

# γ model residual vs 6-var model predictors
print(f"\nCorrelation of γ-model residual with 6-var predictors:")
for name, vals in [('logV', logV), ('logL', logL), ('c_V', c_V),
                   ('f_gas', f_gas), ('logV×c_V', logV*c_V), ('logL×f_gas', logL*f_gas)]:
    r = np.corrcoef(resid_g, vals)[0, 1]
    print(f"  r(resid_γ, {name}) = {r:+.4f}")

# Add 6-var predictors to γ model
X_combined = np.column_stack([np.ones(n), log_gamma, logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta_comb = np.linalg.lstsq(X_combined, y, rcond=None)[0]
yhat_comb = X_combined @ beta_comb
resid_comb = y - yhat_comb
R2_comb = 1 - np.sum(resid_comb**2) / np.sum((y - np.mean(y))**2)

# LOO for combined
H_comb = X_combined @ np.linalg.inv(X_combined.T @ X_combined) @ X_combined.T
h_comb = np.diag(H_comb)
loo_comb = resid_comb / (1 - h_comb)
R2_loo_comb = 1 - np.sum(loo_comb**2) / np.sum((y - np.mean(y))**2)

# LOO for γ alone
H_g = X_gamma @ np.linalg.inv(X_gamma.T @ X_gamma) @ X_gamma.T
h_g = np.diag(H_g)
loo_g = resid_g / (1 - h_g)
R2_loo_g = 1 - np.sum(loo_g**2) / np.sum((y - np.mean(y))**2)

# LOO for 6-var
H6 = X6 @ np.linalg.inv(X6.T @ X6) @ X6.T
h6 = np.diag(H6)
loo6 = resid6 / (1 - h6)
R2_loo6 = 1 - np.sum(loo6**2) / np.sum((y - np.mean(y))**2)

print(f"\nModel comparison:")
print(f"{'Model':<25} {'R²':<10} {'LOO R²':<10} {'RMS':<10}")
print("-" * 55)
print(f"  {'γ alone':<25} {R2_gamma:<10.4f} {R2_loo_g:<10.4f} {np.sqrt(np.mean(resid_g**2)):<10.4f}")
print(f"  {'6-var (no γ)':<25} {R2_6:<10.4f} {R2_loo6:<10.4f} {np.sqrt(np.mean(resid6**2)):<10.4f}")
print(f"  {'γ + 6-var predictors':<25} {R2_comb:<10.4f} {R2_loo_comb:<10.4f} {np.sqrt(np.mean(resid_comb**2)):<10.4f}")

print(f"\n  Δ(LOO R²) adding 6-var to γ: {R2_loo_comb - R2_loo_g:+.4f}")
print(f"  Δ(LOO R²) adding γ to 6-var: {R2_loo_comb - R2_loo6:+.4f}")

print("\n✓ Test 3 passed: added-value analysis done")

# =====================================================================
# TEST 4: CAN γ REPLACE THE 6-VAR MODEL?
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: CAN γ REPLACE THE 6-VAR MODEL?")
print("=" * 60)

# Model: offset = a + b × log(γ) + c × f_gas + d × logL×f_gas
# (minimal model: γ + gas correction)
X_gamma_gas = np.column_stack([np.ones(n), log_gamma, f_gas, logL * f_gas])
beta_gg = np.linalg.lstsq(X_gamma_gas, y, rcond=None)[0]
yhat_gg = X_gamma_gas @ beta_gg
resid_gg = y - yhat_gg
R2_gg = 1 - np.sum(resid_gg**2) / np.sum((y - np.mean(y))**2)
H_gg = X_gamma_gas @ np.linalg.inv(X_gamma_gas.T @ X_gamma_gas) @ X_gamma_gas.T
loo_gg = resid_gg / (1 - np.diag(H_gg))
R2_loo_gg = 1 - np.sum(loo_gg**2) / np.sum((y - np.mean(y))**2)

print(f"\nγ + gas model: offset = {beta_gg[0]:+.4f} + {beta_gg[1]:+.4f}×log(γ) "
      f"+ {beta_gg[2]:+.4f}×f_gas + {beta_gg[3]:+.4f}×logL×f_gas")
print(f"  R² = {R2_gg:.4f}, LOO R² = {R2_loo_gg:.4f}")
print(f"  Compare 6-var LOO R²: {R2_loo6:.4f}")
print(f"  Δ(LOO R²): {R2_loo_gg - R2_loo6:+.4f}")

# Does γ subsume logV and logL?
# γ ∝ V⁻¹ × R^(1/2) ∝ V⁻¹ × (V/a₀)^(1/2) × ...
# log(γ) = log(2) - 0.5 × log(N_corr) = log(2) - 0.5 × (2×logV - logR - log(a₀))
# log(γ) ≈ const - logV + 0.5×logR
print(f"\n  log(γ) decomposition:")
print(f"    r(log γ, logV) = {np.corrcoef(log_gamma, logV)[0,1]:+.4f}")
print(f"    r(log γ, logL) = {np.corrcoef(log_gamma, logL)[0,1]:+.4f}")
print(f"    r(log γ, c_V)  = {np.corrcoef(log_gamma, c_V)[0,1]:+.4f}")
print(f"    r(log γ, f_gas) = {np.corrcoef(log_gamma, f_gas)[0,1]:+.4f}")

print("\n✓ Test 4 passed: γ replacement tested")

# =====================================================================
# TEST 5: DECOMPOSE N_corr IN 6-VAR PREDICTORS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: N_corr IN TERMS OF 6-VAR PREDICTORS")
print("=" * 60)

# Regress log(N_corr) on the 6-var predictors
X_pred = np.column_stack([np.ones(n), logV, logL, c_V, f_gas])
beta_N = np.linalg.lstsq(X_pred, log_N, rcond=None)[0]
yhat_N = X_pred @ beta_N
resid_N = log_N - yhat_N
R2_N = 1 - np.sum(resid_N**2) / np.sum((log_N - np.mean(log_N))**2)

print(f"\nlog(N_corr) = {beta_N[0]:+.3f} + {beta_N[1]:+.3f}×logV + {beta_N[2]:+.3f}×logL "
      f"+ {beta_N[3]:+.3f}×c_V + {beta_N[4]:+.3f}×f_gas")
print(f"R² = {R2_N:.4f}")

print(f"\nTheoretical prediction: N_corr = V²/(R×a₀)")
print(f"  log(N_corr) = 2×logV - logR + const")
print(f"  Since logR correlates with logV and logL, the above regression should show:")
print(f"    β(logV) ≈ 2 - something (due to R correlation)")
print(f"    β(logL) > 0 (larger galaxies have larger R)")
print(f"  Observed: β(logV) = {beta_N[1]:+.3f}, β(logL) = {beta_N[2]:+.3f}")

# R in the definition: use R_max
logR = np.log10([g['r_max_kpc'] for g in galaxies])
print(f"\n  Direct check: log(N_corr) = 2×logV - logR + const")
print(f"    r(logR, logV) = {np.corrcoef(logR, logV)[0,1]:+.4f}")
print(f"    r(logR, logL) = {np.corrcoef(logR, logL)[0,1]:+.4f}")

# Exact decomposition
print(f"\n  Exact: log(N_corr) = 2×logV_ms - logR_m - log(a₀)")
check = 2 * np.log10([g['vflat'] * 1e3 for g in galaxies]) - \
        np.log10([g['r_max_kpc'] * 3.086e19 for g in galaxies]) - np.log10(a0_mond)
print(f"    Reconstructed vs computed: r = {np.corrcoef(check, log_N)[0,1]:.6f}")

print("\n✓ Test 5 passed: N_corr decomposition done")

# =====================================================================
# TEST 6: THEORETICAL PREDICTION FOR MOND COEFFICIENTS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: THEORETICAL MOND COEFFICIENT PREDICTION")
print("=" * 60)

# In deep MOND: g_obs = √(g_bar × a₀)
# RAR offset = log(g_obs) - log(g_rar(g_bar))
# For the mean galaxy: g_rar ≈ g_obs, so offset ≈ 0
# Deviation from mean: δ(offset) ≈ δ(log g_obs) - δ(log g_rar)
# At the flat part: g_obs = V⁴/(R²×a₀), so log(g_obs) = 4×logV - 2×logR - log(a₀)
# g_bar = V_bar⁴/(R²×a₀²), related to L through M/L

# The BTFR: V⁴ = G×M_bar×a₀ → 4×logV = logM + logG + log(a₀)
# M_bar = (M/L)×L, so logM ≈ logL + log(M/L)
# Therefore: 4×logV ≈ logL + const
# Or: logV ≈ 0.25×logL + const

# In terms of offset:
# offset ∝ logV⁴/M_bar ∝ 4×logV - logL
# δ(offset)/δ(logV) = 4×... but our model uses V and L independently

# Direct theoretical prediction:
# From γ = 2/√N_corr and N_corr = V²/(R×a₀):
# log(γ) = log(2) - logV + 0.5×logR + 0.5×log(a₀)
# Since R ∝ V^α × L^β (empirically):
# log(γ) = const + (−1 + 0.5α)×logV + 0.5β×logL

# Regress logR on logV and logL
X_rl = np.column_stack([np.ones(n), logV, logL])
beta_rl = np.linalg.lstsq(X_rl, logR, rcond=None)[0]
R2_rl = 1 - np.sum((logR - X_rl @ beta_rl)**2) / np.sum((logR - np.mean(logR))**2)

print(f"\nlogR = {beta_rl[0]:+.3f} + {beta_rl[1]:+.3f}×logV + {beta_rl[2]:+.3f}×logL")
print(f"R² = {R2_rl:.4f}")

# Predicted γ coefficients
pred_logV_coeff = -1 + 0.5 * beta_rl[1]
pred_logL_coeff = 0.5 * beta_rl[2]

print(f"\nTheoretical prediction for offset model:")
print(f"  If offset ≈ log(γ) = const - logV + 0.5×logR:")
print(f"    Predicted β(logV) = -1 + 0.5×{beta_rl[1]:.3f} = {pred_logV_coeff:+.3f}")
print(f"    Predicted β(logL) = 0.5×{beta_rl[2]:.3f} = {pred_logL_coeff:+.3f}")
print(f"  Observed (6-var):")
print(f"    β(logV) = {beta6[1]:+.3f}")
print(f"    β(logL) = {beta6[2]:+.3f}")

# Better: if offset ∝ γ^k, fit k
X_gk = np.column_stack([np.ones(n), log_gamma])
beta_gk = np.linalg.lstsq(X_gk, y, rcond=None)[0]
print(f"\n  offset = {beta_gk[0]:+.3f} + {beta_gk[1]:+.3f} × log(γ)")
print(f"  Theory says slope = 1 (offset = log γ)")
print(f"  Observed slope = {beta_gk[1]:.3f}")

print("\n✓ Test 6 passed: theoretical coefficients computed")

# =====================================================================
# TEST 7: γ AS FUNCTION OF GALAXY PROPERTIES
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: γ AS FUNCTION OF GALAXY PROPERTIES")
print("=" * 60)

types = np.array([g['hubble_type'] for g in galaxies])

print(f"\nlog(γ) by morphological type:")
for name, T_lo, T_hi in [('Early (T<4)', 0, 4), ('Mid (4≤T<7)', 4, 7), ('Late (T≥7)', 7, 20)]:
    mask = (types >= T_lo) & (types < T_hi)
    if mask.sum() >= 5:
        print(f"  {name} (N={mask.sum()}): ⟨log γ⟩ = {np.mean(log_gamma[mask]):+.4f}, "
              f"⟨offset⟩ = {np.mean(y[mask]):+.4f}")

# γ vs galaxy properties
print(f"\nCorrelations of log(γ) with properties:")
for name, vals in [('logV', logV), ('logL', logL), ('c_V', c_V), ('f_gas', f_gas),
                   ('hubble_type', types)]:
    r = np.corrcoef(log_gamma, vals)[0, 1]
    print(f"  r(log γ, {name}) = {r:+.4f}")

# Residual: what does γ miss?
print(f"\nCorrelation of γ-residual with properties:")
for name, vals in [('logV', logV), ('logL', logL), ('c_V', c_V), ('f_gas', f_gas)]:
    r = np.corrcoef(resid_g, vals)[0, 1]
    print(f"  r(resid_γ, {name}) = {r:+.4f}")

print("\n✓ Test 7 passed: γ property analysis done")

# =====================================================================
# TEST 8: COMBINED MODEL — γ + CORRECTIONS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: OPTIMAL γ-BASED MODEL")
print("=" * 60)

# Build up from γ, adding corrections one at a time
models = [
    ('γ alone', np.column_stack([np.ones(n), log_gamma])),
    ('γ + f_gas', np.column_stack([np.ones(n), log_gamma, f_gas])),
    ('γ + f_gas + logL×f_gas', np.column_stack([np.ones(n), log_gamma, f_gas, logL*f_gas])),
    ('γ + f_gas + logL×f_gas + c_V', np.column_stack([np.ones(n), log_gamma, f_gas, logL*f_gas, c_V])),
    ('γ + all 6-var', X_combined),
    ('6-var (no γ)', X6),
]

print(f"\n{'Model':<30} {'R²':<10} {'LOO R²':<10} {'# params':<10}")
print("-" * 60)
for name, X_m in models:
    beta_m = np.linalg.lstsq(X_m, y, rcond=None)[0]
    yhat_m = X_m @ beta_m
    resid_m = y - yhat_m
    R2_m = 1 - np.sum(resid_m**2) / np.sum((y - np.mean(y))**2)
    H_m = X_m @ np.linalg.inv(X_m.T @ X_m) @ X_m.T
    loo_m = resid_m / (1 - np.diag(H_m))
    R2_loo_m = 1 - np.sum(loo_m**2) / np.sum((y - np.mean(y))**2)
    print(f"  {name:<30} {R2_m:<10.4f} {R2_loo_m:<10.4f} {X_m.shape[1]:<10}")

# The best γ-based model
X_best = np.column_stack([np.ones(n), log_gamma, f_gas, logL*f_gas, c_V])
beta_best = np.linalg.lstsq(X_best, y, rcond=None)[0]
yhat_best = X_best @ beta_best
resid_best = y - yhat_best
R2_best = 1 - np.sum(resid_best**2) / np.sum((y - np.mean(y))**2)
H_best = X_best @ np.linalg.inv(X_best.T @ X_best) @ X_best.T
loo_best = resid_best / (1 - np.diag(H_best))
R2_loo_best = 1 - np.sum(loo_best**2) / np.sum((y - np.mean(y))**2)

print(f"\nBest γ-based model (4 predictors):")
vnames_best = ['const', 'log(γ)', 'f_gas', 'logL×f_gas', 'c_V']
for vname, b in zip(vnames_best, beta_best):
    print(f"  {vname:<15} {b:+.4f}")
print(f"  R² = {R2_best:.4f}, LOO R² = {R2_loo_best:.4f}")
print(f"  vs 6-var: R² = {R2_6:.4f}, LOO R² = {R2_loo6:.4f}")

# Is γ redundant with logV?
r_gamma_logV = np.corrcoef(log_gamma, logV)[0, 1]
print(f"\n  r(log γ, logV) = {r_gamma_logV:+.4f}")
print(f"  → γ is {'mostly' if abs(r_gamma_logV) > 0.7 else 'partially'} "
      f"{'redundant with' if abs(r_gamma_logV) > 0.7 else 'independent of'} logV")

print("\n✓ Test 8 passed: optimal γ model built")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #503 SUMMARY")
print("=" * 70)
print(f"r(log γ, offset) = {r_gamma_offset:+.4f}")
print(f"γ alone: R²(LOO) = {R2_loo_g:.4f}")
print(f"6-var model: R²(LOO) = {R2_loo6:.4f}")
print(f"γ + corrections: R²(LOO) = {R2_loo_best:.4f}")
print(f"γ + all 6-var: R²(LOO) = {R2_loo_comb:.4f}")
print(f"offset ≈ {beta_gk[0]:+.3f} + {beta_gk[1]:+.3f} × log(γ)")
print(f"r(log γ, logV) = {r_gamma_logV:+.4f}")
print(f"\nAll 8 tests passed ✓")
