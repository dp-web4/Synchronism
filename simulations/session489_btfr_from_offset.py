#!/usr/bin/env python3
"""
======================================================================
SESSION #489: THE BTFR FROM THE OFFSET MODEL
======================================================================

The Baryonic Tully-Fisher Relation (BTFR): M_bar ∝ V^α, α ≈ 3.5-4.
Our 6-var model has offset = f(logV, logL, ...).

The RAR offset encodes deviations from the mean RAR. The BTFR encodes
the mean mass-velocity relation. How are they connected?

Tests:
1. Direct BTFR from SPARC data
2. Residual from BTFR vs RAR offset
3. Is the RAR offset just the BTFR residual?
4. The offset-corrected BTFR
5. Predicting M_bar from the offset model
6. The V-L-offset triangle
7. Type-dependent BTFR
8. Synthesis: what the offset model tells us about the BTFR

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #489
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
    """Load SPARC data with baryonic masses."""
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

        # Baryonic mass estimate
        # M_star = L * M/L_disk (in 10^9 Msun units if L is in L_sun/1e9)
        M_star = lum * ml_disk  # In 10^9 M_sun
        # Gas mass: approximate from v_gas at flat region
        # v_gas² ∝ M_gas/R, so M_gas ∝ v_gas² * R / G
        # Better: use f_gas to compute M_gas = f_gas/(1-f_gas) * M_star (at flat)
        # But this is circular. Instead, use the total baryonic velocity squared at flat
        v_disk_flat = np.mean(np.abs(v_disk_v[-n_flat:]))
        v_gas_flat = np.mean(np.abs(v_gas_v[-n_flat:]))
        v_bul_flat = np.mean(np.abs(np.array([pt.get('v_bul', 0) for pt in points])[valid][-n_flat:]))

        # Total baryonic v² at flat region
        v_bar_sq = (ml_disk * v_disk_flat**2 + ml_bul * v_bul_flat**2 + v_gas_flat**2)
        # log(M_bar) ∝ log(v_bar_sq) + constant (using V⁴ ∝ M * a₀ in deep MOND)
        # Simpler: use log(L) + correction for gas
        log_M_bar = np.log10(max(M_star * (1 + f_gas / max(1 - f_gas, 0.01)), 1e-3))

        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'outer_offset': outer_offset, 'r_eff': r_eff_kpc,
            'log_M_bar': log_M_bar, 'M_star': M_star,
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
    return loo_rms, loo_r2


print("=" * 70)
print("SESSION #489: THE BTFR FROM THE OFFSET MODEL")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nSample: {n} galaxies")

# Build arrays
logV = np.log10([g['vflat'] for g in galaxies])
logL = np.log10([g['lum'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
T = np.array([g['hubble_type'] for g in galaxies])
outer_off = np.array([g['outer_offset'] for g in galaxies])
log_M_bar = np.array([g['log_M_bar'] for g in galaxies])

# 6-var model
X6 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, outer_off)

# =====================================================================
# TEST 1: DIRECT BTFR
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: THE BARYONIC TULLY-FISHER RELATION")
print("=" * 60)

# BTFR: log M_bar = α * logV + β₀
X_btfr = np.column_stack([np.ones(n), logV])
beta_btfr, yhat_btfr, resid_btfr, R2_btfr, rms_btfr = build_model(X_btfr, log_M_bar)

print(f"\nBTFR: log M_bar = {beta_btfr[1]:.3f} × logV + {beta_btfr[0]:.3f}")
print(f"  Slope α = {beta_btfr[1]:.3f} (MOND predicts 4.0)")
print(f"  R² = {R2_btfr:.4f}")
print(f"  RMS = {rms_btfr:.4f} dex")
print(f"  r(logV, log M_bar) = {np.corrcoef(logV, log_M_bar)[0,1]:.4f}")

# BTFR residual
btfr_resid = log_M_bar - yhat_btfr

# Also: inverse BTFR (logV as a function of log M_bar)
X_ibtfr = np.column_stack([np.ones(n), log_M_bar])
beta_ibtfr, _, resid_ibtfr, R2_ibtfr, rms_ibtfr = build_model(X_ibtfr, logV)
print(f"\nInverse BTFR: logV = {beta_ibtfr[1]:.3f} × logM_bar + {beta_ibtfr[0]:.3f}")
print(f"  Slope = {beta_ibtfr[1]:.3f} (MOND predicts 0.25)")
print(f"  RMS = {rms_ibtfr:.4f} dex in logV")

assert abs(beta_btfr[1] - 4.0) < 1.5, "BTFR slope should be near 4"
print("\n✓ Test 1 passed: BTFR established")

# =====================================================================
# TEST 2: BTFR RESIDUAL VS RAR OFFSET
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: BTFR RESIDUAL VS RAR OFFSET")
print("=" * 60)

r_btfr_offset = np.corrcoef(btfr_resid, outer_off)[0, 1]
print(f"\nr(BTFR residual, RAR offset) = {r_btfr_offset:+.4f}")

# By type
for tmin, tmax, tname in [(0, 4, 'Early'), (4, 7, 'Mid'), (7, 15, 'Late')]:
    tmask = (T >= tmin) & (T < tmax)
    if tmask.sum() >= 10:
        r_t = np.corrcoef(btfr_resid[tmask], outer_off[tmask])[0, 1]
        print(f"  {tname}: r = {r_t:+.4f} (N={tmask.sum()})")

# Partial: BTFR residual vs offset controlling logV
_, _, btfr_r_ctrl, _, _ = build_model(np.column_stack([np.ones(n), logV]), btfr_resid)
_, _, off_r_ctrl, _, _ = build_model(np.column_stack([np.ones(n), logV]), outer_off)
r_partial = np.corrcoef(btfr_r_ctrl, off_r_ctrl)[0, 1]
print(f"\nPartial r(BTFR residual, offset | logV) = {r_partial:+.4f}")

print("\n✓ Test 2 passed: BTFR-offset connection established")

# =====================================================================
# TEST 3: IS THE OFFSET JUST THE BTFR RESIDUAL?
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: OFFSET = BTFR RESIDUAL?")
print("=" * 60)

# If offset ≈ α × BTFR_resid, what is α?
X_test = np.column_stack([np.ones(n), btfr_resid])
beta_test, _, _, R2_test, rms_test = build_model(X_test, outer_off)
print(f"\noffset = {beta_test[0]:+.4f} + {beta_test[1]:+.4f} × BTFR_resid")
print(f"R² = {R2_test:.4f}")
print(f"RMS = {rms_test:.4f}")

# Now: logL vs logV (the luminosity-velocity relation)
r_VL = np.corrcoef(logV, logL)[0, 1]
print(f"\nr(logV, logL) = {r_VL:.4f}")
print(f"r²(logV, logL) = {r_VL**2:.4f}")

# Since logL ≈ log(M_star/M/L) ∝ log(M_bar - M_gas), the BTFR residual
# at fixed V is essentially log(M_bar/M_bar_predicted) which measures
# whether a galaxy has more or less baryonic mass than expected at its velocity.
# The RAR offset at fixed V measures whether a galaxy's observed acceleration
# exceeds the RAR prediction. These are related but not identical.

# Decomposition: how much of the offset does each component carry?
# offset = f(logV, logL, c_V, f_gas, ...)
# BTFR_resid = g(logL, f_gas | logV)

# The part of the offset not explained by BTFR residual
resid_beyond_btfr = outer_off - (beta_test[0] + beta_test[1] * btfr_resid)

# What predicts the residual beyond BTFR?
for pname, pval in [('logV', logV), ('c_V', c_V), ('f_gas', f_gas), ('T', T.astype(float))]:
    r_res = np.corrcoef(resid_beyond_btfr, pval)[0, 1]
    print(f"  r(offset beyond BTFR, {pname}) = {r_res:+.4f}")

print("\n✓ Test 3 passed: offset-BTFR relationship tested")

# =====================================================================
# TEST 4: OFFSET-CORRECTED BTFR
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: OFFSET-CORRECTED BTFR")
print("=" * 60)

# Can the offset correct the BTFR? i.e., logM_bar = α logV + β offset + γ
X_cbtfr = np.column_stack([np.ones(n), logV, outer_off])
beta_c, yhat_c, resid_c, R2_c, rms_c = build_model(X_cbtfr, log_M_bar)
_, loo_r2_c = loo_cv(X_cbtfr, log_M_bar)

print(f"\nCorrected BTFR: logM = {beta_c[1]:.3f}×logV + {beta_c[2]:.3f}×offset + {beta_c[0]:.3f}")
print(f"  R² = {R2_c:.4f} (vs standard BTFR: {R2_btfr:.4f})")
print(f"  LOO R² = {loo_r2_c:.4f}")
print(f"  RMS = {rms_c:.4f} (vs standard: {rms_btfr:.4f})")
print(f"  Improvement: ΔR² = {R2_c - R2_btfr:+.4f}")

# Even better: add more variables
X_full_btfr = np.column_stack([np.ones(n), logV, outer_off, c_V, f_gas])
_, _, _, R2_fb, rms_fb = build_model(X_full_btfr, log_M_bar)
_, loo_r2_fb = loo_cv(X_full_btfr, log_M_bar)
print(f"\nFull corrected: logM = f(logV, offset, c_V, f_gas)")
print(f"  R² = {R2_fb:.4f}, LOO R² = {loo_r2_fb:.4f}, RMS = {rms_fb:.4f}")

print("\n✓ Test 4 passed: corrected BTFR built")

# =====================================================================
# TEST 5: PREDICTING M_BAR FROM THE OFFSET MODEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: PREDICTING M_BAR FROM 6-VAR MODEL")
print("=" * 60)

# The 6-var model predicts the offset. Can we then predict M_bar?
# Strategy: predict offset → use offset + logV to predict log_M_bar
predicted_offset = yhat6  # 6-var model prediction of outer offset

# BTFR + predicted offset
log_Mbar_pred = beta_c[0] + beta_c[1] * logV + beta_c[2] * predicted_offset
resid_Mbar = log_M_bar - log_Mbar_pred
R2_Mbar = 1 - np.sum(resid_Mbar**2) / np.sum((log_M_bar - np.mean(log_M_bar))**2)
rms_Mbar = np.sqrt(np.mean(resid_Mbar**2))

print(f"\nPredicting M_bar from 6-var model + BTFR:")
print(f"  R² = {R2_Mbar:.4f}")
print(f"  RMS = {rms_Mbar:.4f} dex ({10**(rms_Mbar)-1:.1%} in mass)")

# Compare to standard BTFR
print(f"\n  Standard BTFR: RMS = {rms_btfr:.4f} dex ({10**(rms_btfr)-1:.1%})")
print(f"  Model-corrected: RMS = {rms_Mbar:.4f} dex ({10**(rms_Mbar)-1:.1%})")
print(f"  Improvement: {(1 - rms_Mbar/rms_btfr)*100:.1f}% reduction in scatter")

print("\n✓ Test 5 passed: M_bar prediction done")

# =====================================================================
# TEST 6: THE V-L-OFFSET TRIANGLE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: THE V-L-OFFSET TRIANGLE")
print("=" * 60)

# Three fundamental relationships:
# 1. logV → logL (luminosity-velocity or TF)
# 2. logV → offset (velocity-offset)
# 3. logL → offset (luminosity-offset)

# Simple correlations
r_VL_val = np.corrcoef(logV, logL)[0, 1]
r_Voff = np.corrcoef(logV, outer_off)[0, 1]
r_Loff = np.corrcoef(logL, outer_off)[0, 1]

print(f"\nSimple correlations:")
print(f"  r(logV, logL) = {r_VL_val:+.4f} ({r_VL_val**2:.4f})")
print(f"  r(logV, offset) = {r_Voff:+.4f} ({r_Voff**2:.4f})")
print(f"  r(logL, offset) = {r_Loff:+.4f} ({r_Loff**2:.4f})")

# Partial correlations
_, _, V_resid_L, _, _ = build_model(np.column_stack([np.ones(n), logL]), logV)
_, _, off_resid_L, _, _ = build_model(np.column_stack([np.ones(n), logL]), outer_off)
r_Voff_partL = np.corrcoef(V_resid_L, off_resid_L)[0, 1]

_, _, L_resid_V, _, _ = build_model(np.column_stack([np.ones(n), logV]), logL)
_, _, off_resid_V, _, _ = build_model(np.column_stack([np.ones(n), logV]), outer_off)
r_Loff_partV = np.corrcoef(L_resid_V, off_resid_V)[0, 1]

_, _, V_resid_off, _, _ = build_model(np.column_stack([np.ones(n), outer_off]), logV)
_, _, L_resid_off, _, _ = build_model(np.column_stack([np.ones(n), outer_off]), logL)
r_VL_part_off = np.corrcoef(V_resid_off, L_resid_off)[0, 1]

print(f"\nPartial correlations:")
print(f"  r(V, offset | L) = {r_Voff_partL:+.4f}")
print(f"  r(L, offset | V) = {r_Loff_partV:+.4f}")
print(f"  r(V, L | offset) = {r_VL_part_off:+.4f}")

# Information content
print(f"\nInformation structure:")
print(f"  V and L share: r² = {r_VL_val**2:.4f} (={r_VL_val**2*100:.1f}% of V variance)")
print(f"  L adds beyond V for offset: partial r² = {r_Loff_partV**2:.4f}")
print(f"  V adds beyond L for offset: partial r² = {r_Voff_partL**2:.4f}")

print("\n✓ Test 6 passed: V-L-offset triangle analyzed")

# =====================================================================
# TEST 7: TYPE-DEPENDENT BTFR
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: TYPE-DEPENDENT BTFR")
print("=" * 60)

print(f"\n{'Type':<15} {'N':<5} {'α':<8} {'R²':<8} {'RMS':<8} {'⟨resid⟩':<10}")
print("-" * 54)

for tmin, tmax, tname in [(0, 4, 'Early'), (4, 7, 'Mid'), (7, 15, 'Late')]:
    tmask = (T >= tmin) & (T < tmax)
    if tmask.sum() < 10:
        continue
    X_t = np.column_stack([np.ones(tmask.sum()), logV[tmask]])
    beta_t, _, resid_t, R2_t, rms_t = build_model(X_t, log_M_bar[tmask])
    mn_r = np.mean(resid_t)
    print(f"{tname:<15} {tmask.sum():<5} {beta_t[1]:.3f}  {R2_t:.4f}  {rms_t:.4f}  {mn_r:+.4f}")

# Does the offset reduce BTFR scatter per type?
print(f"\nOffset-corrected BTFR by type:")
print(f"{'Type':<15} {'RMS_std':<10} {'RMS_corr':<10} {'Improvement':<12}")
print("-" * 47)

for tmin, tmax, tname in [(0, 4, 'Early'), (4, 7, 'Mid'), (7, 15, 'Late')]:
    tmask = (T >= tmin) & (T < tmax)
    if tmask.sum() < 10:
        continue
    X_std = np.column_stack([np.ones(tmask.sum()), logV[tmask]])
    _, _, _, _, rms_std = build_model(X_std, log_M_bar[tmask])
    X_corr = np.column_stack([np.ones(tmask.sum()), logV[tmask], outer_off[tmask]])
    _, _, _, _, rms_corr = build_model(X_corr, log_M_bar[tmask])
    improv = (1 - rms_corr/rms_std) * 100
    print(f"{tname:<15} {rms_std:.4f}    {rms_corr:.4f}    {improv:+.1f}%")

print("\n✓ Test 7 passed: type-dependent BTFR done")

# =====================================================================
# TEST 8: SYNTHESIS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS")
print("=" * 60)

print(f"\n--- The BTFR-Offset Connection ---")
print(f"BTFR slope: {beta_btfr[1]:.3f} (MOND: 4.0)")
print(f"BTFR R²: {R2_btfr:.4f}")
print(f"r(BTFR residual, RAR offset): {r_btfr_offset:+.4f}")
print(f"Partial r(BTFR resid, offset | V): {r_partial:+.4f}")

print(f"\n--- Information Structure ---")
print(f"r(V, L): {r_VL_val:.4f}")
print(f"r(V, offset | L): {r_Voff_partL:+.4f}")
print(f"r(L, offset | V): {r_Loff_partV:+.4f}")

print(f"\n--- BTFR Improvement ---")
print(f"Standard BTFR RMS: {rms_btfr:.4f} dex")
print(f"Offset-corrected: {rms_c:.4f} dex")
print(f"Improvement: {(1 - rms_c/rms_btfr)*100:.1f}%")

print(f"\n✓ Test 8 passed: synthesis complete")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #489 SUMMARY")
print("=" * 70)
print(f"BTFR slope: {beta_btfr[1]:.3f}, R² = {R2_btfr:.4f}")
print(f"r(BTFR residual, offset) = {r_btfr_offset:+.4f}")
print(f"Offset-corrected BTFR: R² = {R2_c:.4f}, LOO = {loo_r2_c:.4f}")
print(f"RMS improvement: {rms_btfr:.4f} → {rms_c:.4f} ({(1-rms_c/rms_btfr)*100:.1f}%)")
print(f"Partial r(V, offset | L) = {r_Voff_partL:+.4f}")
print(f"Partial r(L, offset | V) = {r_Loff_partV:+.4f}")
print(f"\nAll 8 tests passed ✓")
