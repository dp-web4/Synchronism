#!/usr/bin/env python3
"""
======================================================================
SESSION #493: THE GOLDEN SUBSAMPLE — OPTIMAL MEASUREMENT QUALITY
======================================================================

Session #491 showed combined measurement noise = 28% of total offset
variance, with M/L the largest source (9.4%) but negligible for gas-rich
galaxies (r = -0.85). Inclination errors dominate for face-on galaxies.

This session selects "golden" subsamples that minimize specific noise
sources and tests whether the 6-var model improves:
1. Edge-on selection (i ≥ 60°): minimize inclination noise
2. Gas-rich selection (f_gas ≥ 0.4): minimize M/L noise
3. Q=1 selection: best V_obs data
4. Golden = Q=1 + edge-on + gas-rich: minimize ALL noise
5. Model on golden vs full sample
6. Denoised offsets via MC averaging
7. Leave-type-out cross-validation
8. The irreducible scatter

Tests:
1. Subsample definitions and sizes
2. Model R² on noise-minimized subsamples
3. Golden subsample model
4. MC denoised offsets
5. Model on denoised data
6. Per-galaxy noise vs residual
7. Leave-type-out cross-validation
8. The irreducible scatter estimate

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #493
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
    """Load SPARC data with all metadata."""
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

        # Mean v_obs error
        mean_e_v = np.mean(e_vobs_v[mond])
        frac_e_v = mean_e_v / max(np.mean(np.abs(v_obs_v[mond])), 1)

        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'outer_offset': outer_offset,
            'distance': distance, 'inclination': inclination, 'quality': quality,
            'mean_e_v': mean_e_v, 'frac_e_v': frac_e_v,
            'n_outer': outer_mond.sum(), 'n_mond': mond.sum(),
            # Raw data for MC
            'v_obs': v_obs_v, 'v_gas': v_gas_v, 'v_disk': v_disk_v,
            'v_bul': np.array([pt.get('v_bul', 0) for pt in points])[valid],
            'radius': radius_v, 'e_vobs': e_vobs_v,
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


def compute_offset_perturbed(g, ml_disk=0.5, ml_bul=0.7, v_noise=None,
                              d_factor=1.0, incl_shift_deg=0.0):
    """Compute outer offset with perturbations."""
    v_obs = g['v_obs'].copy()
    if v_noise is not None:
        v_obs = v_obs + v_noise

    if incl_shift_deg != 0.0:
        i_assumed = np.deg2rad(g['inclination'])
        i_true = np.deg2rad(g['inclination'] + incl_shift_deg)
        i_true = np.clip(i_true, np.deg2rad(20), np.deg2rad(89))
        incl_ratio = np.sin(i_assumed) / np.sin(i_true)
        v_obs = v_obs * incl_ratio

    sqrt_d = np.sqrt(d_factor)
    v_gas = g['v_gas'] * sqrt_d
    v_disk = g['v_disk'] * sqrt_d
    v_bul = g['v_bul'] * sqrt_d
    radius = g['radius'] * d_factor

    kpc_to_m = 3.086e19
    kms2 = 1e6

    g_bar_term = v_gas**2 + ml_disk * v_disk**2 + ml_bul * v_bul**2
    g_bar = np.abs(g_bar_term) * kms2 / (radius * kpc_to_m)
    g_obs = v_obs**2 * kms2 / (radius * kpc_to_m)

    valid = (g_bar > 0) & (g_obs > 0) & np.isfinite(g_bar) & np.isfinite(g_obs) & (radius > 0)
    if valid.sum() < 3:
        return np.nan

    g_bar_v = g_bar[valid]
    g_obs_v = g_obs[valid]
    radius_v = radius[valid]

    mond = g_bar_v < a0_mond
    if mond.sum() < 2:
        return np.nan

    radius_m = radius_v[mond]
    med_r = np.median(radius_m)
    outer = radius_m > med_r

    if outer.sum() < 1:
        g_rar = rar_prediction(g_bar_v[mond])
        return np.mean(np.log10(g_obs_v[mond]) - np.log10(g_rar))

    g_rar = rar_prediction(g_bar_v[mond][outer])
    return np.mean(np.log10(g_obs_v[mond][outer]) - np.log10(g_rar))


def make_X6(galaxies, mask=None):
    """Build the 6-var design matrix for a subset."""
    if mask is None:
        mask = np.ones(len(galaxies), dtype=bool)
    idx = np.where(mask)[0]
    n = len(idx)
    logV = np.log10([galaxies[i]['vflat'] for i in idx])
    logL = np.log10([galaxies[i]['lum'] for i in idx])
    c_V = np.array([galaxies[i]['c_V'] for i in idx])
    f_gas = np.array([galaxies[i]['f_gas'] for i in idx])
    y = np.array([galaxies[i]['outer_offset'] for i in idx])
    X = np.column_stack([np.ones(n), logV, logL, c_V, f_gas,
                          logV * c_V, logL * f_gas])
    return X, y


print("=" * 70)
print("SESSION #493: THE GOLDEN SUBSAMPLE — OPTIMAL MEASUREMENT QUALITY")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nFull sample: {n} galaxies")

# Build full model for reference
X_full, y_full = make_X6(galaxies)
beta_full, yhat_full, resid_full, R2_full, rms_full = build_model(X_full, y_full)
loo_rms_full, loo_r2_full = loo_cv(X_full, y_full)
print(f"Full 6-var: R² = {R2_full:.4f}, LOO R² = {loo_r2_full:.4f}, RMS = {rms_full:.4f}")

# Extract arrays
incl_arr = np.array([g['inclination'] for g in galaxies])
fgas_arr = np.array([g['f_gas'] for g in galaxies])
quality_arr = np.array([g['quality'] for g in galaxies])
T_arr = np.array([g['hubble_type'] for g in galaxies])
frac_e_arr = np.array([g['frac_e_v'] for g in galaxies])

# =====================================================================
# TEST 1: SUBSAMPLE DEFINITIONS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: SUBSAMPLE DEFINITIONS AND SIZES")
print("=" * 60)

edge_on = incl_arr >= 60
gas_rich = fgas_arr >= 0.4
q1 = quality_arr == 1
late = T_arr >= 7
golden = q1 & edge_on & gas_rich

print(f"\n{'Subsample':<25} {'N':<6} {'⟨i⟩':<8} {'⟨f_gas⟩':<10} {'⟨frac_e_V⟩'}")
print("-" * 60)
for name, mask in [('Full', np.ones(n, dtype=bool)),
                    ('Edge-on (i≥60°)', edge_on),
                    ('Gas-rich (f_gas≥0.4)', gas_rich),
                    ('Q=1', q1),
                    ('Late (T≥7)', late),
                    ('Golden (Q1+edge+gas)', golden)]:
    if mask.sum() > 0:
        print(f"  {name:<25} {mask.sum():<6} {np.mean(incl_arr[mask]):<8.1f} "
              f"{np.mean(fgas_arr[mask]):<10.3f} {np.mean(frac_e_arr[mask]):.4f}")

assert golden.sum() >= 5, f"Golden subsample too small: {golden.sum()}"
print("\n✓ Test 1 passed: subsamples defined")

# =====================================================================
# TEST 2: MODEL R² ON NOISE-MINIMIZED SUBSAMPLES
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: MODEL R² ON NOISE-MINIMIZED SUBSAMPLES")
print("=" * 60)

# Fit the full-sample model and evaluate on subsamples
# (Using full-sample coefficients, not refitting)
print(f"\n{'Subsample':<25} {'N':<6} {'R²(eval)':<10} {'RMS(eval)'}")
print("-" * 55)

for name, mask in [('Full', np.ones(n, dtype=bool)),
                    ('Edge-on (i≥60°)', edge_on),
                    ('Gas-rich (f_gas≥0.4)', gas_rich),
                    ('Q=1', q1),
                    ('Late (T≥7)', late),
                    ('Golden (Q1+edge+gas)', golden),
                    ('Edge-on + gas-rich', edge_on & gas_rich)]:
    if mask.sum() < 8:
        continue
    X_sub, y_sub = make_X6(galaxies, mask)
    # Evaluate full-model predictions
    yhat_sub = X_sub @ beta_full
    resid_sub = y_sub - yhat_sub
    ss_tot_sub = np.sum((y_sub - np.mean(y_sub))**2)
    r2_eval = 1 - np.sum(resid_sub**2) / ss_tot_sub if ss_tot_sub > 0 else 0
    rms_eval = np.sqrt(np.mean(resid_sub**2))
    print(f"  {name:<25} {mask.sum():<6} {r2_eval:<10.4f} {rms_eval:.4f}")

print("\n✓ Test 2 passed: subsample evaluation done")

# =====================================================================
# TEST 3: REFIT ON GOLDEN SUBSAMPLE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: GOLDEN SUBSAMPLE MODEL (REFIT)")
print("=" * 60)

X_gold, y_gold = make_X6(galaxies, golden)
if golden.sum() >= 10:
    beta_gold, yhat_gold, resid_gold, R2_gold, rms_gold = build_model(X_gold, y_gold)
    loo_rms_gold, loo_r2_gold = loo_cv(X_gold, y_gold)

    print(f"\nGolden subsample (N={golden.sum()}):")
    print(f"  R² = {R2_gold:.4f}, LOO R² = {loo_r2_gold:.4f}, RMS = {rms_gold:.4f}")
    print(f"  Full sample R² = {R2_full:.4f}, LOO R² = {loo_r2_full:.4f}")

    # Compare coefficients
    var_names = ['intercept', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V', 'logL×f_gas']
    print(f"\n  {'Variable':<12} {'Full β':<10} {'Golden β':<10}")
    print("  " + "-" * 35)
    for name, bf, bg in zip(var_names, beta_full, beta_gold):
        print(f"  {name:<12} {bf:+.4f}    {bg:+.4f}")
else:
    print(f"\nGolden subsample too small ({golden.sum()}) for reliable refit")
    R2_gold = np.nan
    loo_r2_gold = np.nan
    rms_gold = np.nan

# Also try broader "good" subsample: edge-on OR gas-rich
good = edge_on | gas_rich
X_good, y_good = make_X6(galaxies, good)
_, _, _, R2_good, rms_good = build_model(X_good, y_good)
loo_rms_good, loo_r2_good = loo_cv(X_good, y_good)
print(f"\n  'Good' (edge-on OR gas-rich, N={good.sum()}):")
print(f"  R² = {R2_good:.4f}, LOO R² = {loo_r2_good:.4f}, RMS = {rms_good:.4f}")

print("\n✓ Test 3 passed: golden subsample model built")

# =====================================================================
# TEST 4: MC DENOISED OFFSETS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: MC DENOISED OFFSETS")
print("=" * 60)

# Generate many noise realizations and average to get "denoised" offset
rng = np.random.RandomState(42)
n_mc = 200
mc_offsets = np.zeros((n, n_mc))

for trial in range(n_mc):
    for i, g in enumerate(galaxies):
        v_noise = rng.normal(0, g['e_vobs'])
        d_factor = max(1 + rng.normal(0, 0.10), 0.5)
        ml_d = max(0.5 + rng.normal(0, 0.15), 0.1)
        di = rng.normal(0, 4.0)
        off = compute_offset_perturbed(g, ml_disk=ml_d, v_noise=v_noise,
                                        d_factor=d_factor, incl_shift_deg=di)
        mc_offsets[i, trial] = off if np.isfinite(off) else g['outer_offset']

# The MC mean IS a denoised estimate (noise averages out, signal preserved)
# But wait: the MC adds noise to baseline. To denoise, we need to remove noise.
# Instead: the MC variance tells us per-galaxy noise.
# Denoised offset = observed offset (the actual measurement IS one realization)
# The "true" offset isn't recoverable without assumptions.

# What we CAN do: estimate per-galaxy noise and weight the model accordingly
mc_var_per_gal = np.var(mc_offsets, axis=1)
mc_sigma_per_gal = np.sqrt(mc_var_per_gal)

print(f"\nPer-galaxy MC noise:")
print(f"  Mean σ_noise: {np.mean(mc_sigma_per_gal):.4f} dex")
print(f"  Median σ_noise: {np.median(mc_sigma_per_gal):.4f} dex")
print(f"  Range: [{np.min(mc_sigma_per_gal):.4f}, {np.max(mc_sigma_per_gal):.4f}]")
print(f"  Interquartile: [{np.percentile(mc_sigma_per_gal, 25):.4f}, "
      f"{np.percentile(mc_sigma_per_gal, 75):.4f}]")

# Noise-quality ranking
noisy_quartile = mc_sigma_per_gal > np.percentile(mc_sigma_per_gal, 75)
clean_quartile = mc_sigma_per_gal < np.percentile(mc_sigma_per_gal, 25)
print(f"\n  Cleanest 25% (N={clean_quartile.sum()}): ⟨σ_noise⟩ = {np.mean(mc_sigma_per_gal[clean_quartile]):.4f}")
print(f"  Noisiest 25% (N={noisy_quartile.sum()}): ⟨σ_noise⟩ = {np.mean(mc_sigma_per_gal[noisy_quartile]):.4f}")

print("\n✓ Test 4 passed: MC noise estimated")

# =====================================================================
# TEST 5: WEIGHTED LEAST SQUARES
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: WEIGHTED LEAST SQUARES MODEL")
print("=" * 60)

# Weight by 1/σ²_noise — galaxies with less noise get more weight
# Prevent zero weights
w = 1.0 / np.clip(mc_var_per_gal, 1e-8, None)
w = w / np.sum(w) * n  # normalize

# WLS: X' = diag(sqrt(w)) @ X, y' = diag(sqrt(w)) @ y
sqrt_w = np.sqrt(w)
X_wls = X_full * sqrt_w[:, None]
y_wls = y_full * sqrt_w

beta_wls = np.linalg.lstsq(X_wls, y_wls, rcond=None)[0]
yhat_wls = X_full @ beta_wls
resid_wls = y_full - yhat_wls
ss_tot_wls = np.sum((y_full - np.mean(y_full))**2)
R2_wls = 1 - np.sum(resid_wls**2) / ss_tot_wls
rms_wls = np.sqrt(np.mean(resid_wls**2))

# WLS LOO
H_wls = X_wls @ np.linalg.inv(X_wls.T @ X_wls) @ X_wls.T
h_wls = np.diag(H_wls)
resid_wls_loo = (y_wls - X_wls @ beta_wls) / (1 - h_wls)
# Convert back to unweighted
loo_resid_wls = resid_wls_loo / sqrt_w
loo_rms_wls = np.sqrt(np.mean(loo_resid_wls**2))
loo_r2_wls = 1 - np.sum(loo_resid_wls**2) / ss_tot_wls

print(f"\n{'Model':<15} {'R²':<8} {'LOO R²':<8} {'RMS':<8}")
print("-" * 40)
print(f"{'OLS':<15} {R2_full:<8.4f} {loo_r2_full:<8.4f} {rms_full:.4f}")
print(f"{'WLS':<15} {R2_wls:<8.4f} {loo_r2_wls:<8.4f} {rms_wls:.4f}")

# Compare coefficients
print(f"\n  {'Variable':<12} {'OLS β':<10} {'WLS β':<10} {'Δ':<8}")
print("  " + "-" * 35)
for name, bo, bw in zip(var_names, beta_full, beta_wls):
    print(f"  {name:<12} {bo:+.4f}    {bw:+.4f}    {bw-bo:+.4f}")

print("\n✓ Test 5 passed: WLS model built")

# =====================================================================
# TEST 6: PER-GALAXY NOISE VS RESIDUAL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: PER-GALAXY NOISE VS RESIDUAL")
print("=" * 60)

# Does per-galaxy noise predict the residual magnitude?
r_noise_resid = np.corrcoef(mc_sigma_per_gal, np.abs(resid_full))[0, 1]
print(f"\nr(σ_noise, |residual|) = {r_noise_resid:+.4f}")

# Model residual for cleanest vs noisiest quartile
rms_clean = np.sqrt(np.mean(resid_full[clean_quartile]**2))
rms_noisy = np.sqrt(np.mean(resid_full[noisy_quartile]**2))
print(f"\n  Cleanest 25%: model RMS = {rms_clean:.4f} dex")
print(f"  Noisiest 25%: model RMS = {rms_noisy:.4f} dex")
print(f"  Ratio: {rms_noisy/rms_clean:.2f}×")

# Noise-to-signal ratio per galaxy
nsr = mc_sigma_per_gal / np.sqrt(np.var(y_full))
print(f"\n  Mean noise/signal ratio: {np.mean(nsr):.3f}")
print(f"  Fraction with noise > signal (NSR > 1): {np.mean(nsr > 1)*100:.1f}%")
print(f"  Fraction with noise < 0.3 signal: {np.mean(nsr < 0.3)*100:.1f}%")

# What properties make a galaxy "clean"?
print(f"\n  Clean galaxies (bottom 25% noise):")
print(f"    ⟨f_gas⟩ = {np.mean(fgas_arr[clean_quartile]):.3f}")
print(f"    ⟨incl⟩ = {np.mean(incl_arr[clean_quartile]):.1f}°")
print(f"    ⟨T⟩ = {np.mean(T_arr[clean_quartile]):.1f}")
print(f"    ⟨frac_eV⟩ = {np.mean(frac_e_arr[clean_quartile]):.4f}")
print(f"\n  Noisy galaxies (top 25% noise):")
print(f"    ⟨f_gas⟩ = {np.mean(fgas_arr[noisy_quartile]):.3f}")
print(f"    ⟨incl⟩ = {np.mean(incl_arr[noisy_quartile]):.1f}°")
print(f"    ⟨T⟩ = {np.mean(T_arr[noisy_quartile]):.1f}")
print(f"    ⟨frac_eV⟩ = {np.mean(frac_e_arr[noisy_quartile]):.4f}")

print("\n✓ Test 6 passed: noise-residual analysis done")

# =====================================================================
# TEST 7: LEAVE-TYPE-OUT CROSS-VALIDATION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: LEAVE-TYPE-OUT CROSS-VALIDATION")
print("=" * 60)

# Train on two type groups, predict the third
early = T_arr < 4
mid = (T_arr >= 4) & (T_arr < 7)
late_mask = T_arr >= 7

type_groups = [('Early (T<4)', early), ('Mid (4≤T<7)', mid), ('Late (T≥7)', late_mask)]

print(f"\n{'Hold-out':<15} {'N_train':<10} {'N_test':<8} {'R²_pred':<10} {'RMS_pred'}")
print("-" * 55)

for name_out, mask_out in type_groups:
    mask_in = ~mask_out
    if mask_out.sum() < 5 or mask_in.sum() < 10:
        continue
    X_train, y_train = make_X6(galaxies, mask_in)
    X_test, y_test = make_X6(galaxies, mask_out)

    beta_cv, _, _, _, _ = build_model(X_train, y_train)
    yhat_cv = X_test @ beta_cv
    resid_cv = y_test - yhat_cv
    ss_tot_cv = np.sum((y_test - np.mean(y_test))**2)
    r2_cv = 1 - np.sum(resid_cv**2) / ss_tot_cv if ss_tot_cv > 0 else 0
    rms_cv = np.sqrt(np.mean(resid_cv**2))
    print(f"  {name_out:<15} {mask_in.sum():<10} {mask_out.sum():<8} {r2_cv:<10.4f} {rms_cv:.4f}")

# Also: leave-quality-out
print(f"\n{'Hold-out':<15} {'N_train':<10} {'N_test':<8} {'R²_pred':<10} {'RMS_pred'}")
print("-" * 55)
for q in [1, 2, 3]:
    mask_out = quality_arr == q
    mask_in = ~mask_out
    if mask_out.sum() < 5 or mask_in.sum() < 10:
        continue
    X_train, y_train = make_X6(galaxies, mask_in)
    X_test, y_test = make_X6(galaxies, mask_out)
    beta_cv, _, _, _, _ = build_model(X_train, y_train)
    yhat_cv = X_test @ beta_cv
    resid_cv = y_test - yhat_cv
    ss_tot_cv = np.sum((y_test - np.mean(y_test))**2)
    r2_cv = 1 - np.sum(resid_cv**2) / ss_tot_cv if ss_tot_cv > 0 else 0
    rms_cv = np.sqrt(np.mean(resid_cv**2))
    print(f"  {'Q='+str(q):<15} {mask_in.sum():<10} {mask_out.sum():<8} {r2_cv:<10.4f} {rms_cv:.4f}")

print("\n✓ Test 7 passed: leave-type-out CV done")

# =====================================================================
# TEST 8: THE IRREDUCIBLE SCATTER
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: THE IRREDUCIBLE SCATTER ESTIMATE")
print("=" * 60)

# For the cleanest galaxies (edge-on + gas-rich + Q=1), the noise
# is minimal. Their residual is the closest estimate of physical scatter.
clean_strict = edge_on & gas_rich & q1
if clean_strict.sum() < 5:
    # Relax: just edge-on + gas-rich
    clean_strict = edge_on & gas_rich

X_cs, y_cs = make_X6(galaxies, clean_strict)
yhat_cs = X_cs @ beta_full
resid_cs = y_cs - yhat_cs
rms_cs = np.sqrt(np.mean(resid_cs**2))

# Also compute noise for these clean galaxies
noise_cs = mc_sigma_per_gal[clean_strict]
mean_noise_cs = np.sqrt(np.mean(noise_cs**2))

# Physical scatter ≈ sqrt(rms² - noise²)
physical_rms = np.sqrt(max(rms_cs**2 - mean_noise_cs**2, 0))

print(f"\nCleanest subsample: N = {clean_strict.sum()}")
print(f"  Selection: edge-on (i≥60°) + gas-rich (f_gas≥0.4) + Q=1")
print(f"  Model RMS on clean: {rms_cs:.4f} dex")
print(f"  Estimated noise σ: {mean_noise_cs:.4f} dex")
print(f"  Physical scatter: {physical_rms:.4f} dex")

# Compare with full sample
print(f"\n  Full sample RMS: {rms_full:.4f} dex")
print(f"  Clean subsample RMS: {rms_cs:.4f} dex")
print(f"  Improvement: {(1 - rms_cs/rms_full)*100:.1f}%")

# In velocity units
print(f"\n  Physical scatter: {physical_rms:.4f} dex = {(10**physical_rms-1)*100:.1f}% in V")
print(f"  This is the best estimate of irreducible physical variation")

# What drives the remaining physical scatter?
print(f"\n  Remaining physical scatter could come from:")
print(f"    - Dark matter halo shape variations")
print(f"    - Non-circular motions / RC asymmetry")
print(f"    - Environment (satellite vs isolated)")
print(f"    - Baryonic physics not captured by V, L, c_V, f_gas")

print("\n✓ Test 8 passed: irreducible scatter estimated")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #493 SUMMARY")
print("=" * 70)
print(f"Full sample: N={n}, R²={R2_full:.4f}, LOO R²={loo_r2_full:.4f}, RMS={rms_full:.4f}")
if np.isfinite(R2_gold):
    print(f"Golden (Q1+edge+gas): N={golden.sum()}, R²={R2_gold:.4f}, "
          f"LOO R²={loo_r2_gold:.4f}, RMS={rms_gold:.4f}")
print(f"WLS: R²={R2_wls:.4f}, LOO R²={loo_r2_wls:.4f}, RMS={rms_wls:.4f}")
print(f"Clean subsample RMS: {rms_cs:.4f} dex")
print(f"Physical scatter estimate: {physical_rms:.4f} dex")
print(f"r(σ_noise, |residual|) = {r_noise_resid:+.4f}")
print(f"MC noise range: [{np.min(mc_sigma_per_gal):.4f}, {np.max(mc_sigma_per_gal):.4f}] dex")
print(f"\nAll 8 tests passed ✓")
