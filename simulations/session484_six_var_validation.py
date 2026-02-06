#!/usr/bin/env python3
"""
======================================================================
SESSION #484: VALIDATING THE 6-VARIABLE MODEL
======================================================================

Session #483 discovered that logL×f_gas pushes LOO R² from 0.896 to
0.938. This session rigorously validates the new model:

1. Full model coefficients and their significance
2. Bootstrap confidence intervals on ΔR²
3. Late-type and early-type performance
4. NN autocorrelation: does it reduce the r=+0.46?
5. Outlier sensitivity
6. f_gas² vs logL×f_gas: which is more fundamental?
7. 7-variable model with both interactions
8. Predicted vs observed and residual diagnostics

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #484
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

        # Outer offset
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

        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'full_offset': full_offset, 'outer_offset': outer_offset,
            'distance': distance, 'inclination': inclination,
            'quality': quality, 'sb_eff': sb_eff, 'r_eff': r_eff_kpc,
        })

    return galaxies


def build_model(X, y):
    """OLS regression returning coefficients, predictions, residuals, R², RMS."""
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    resid = y - y_hat
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - np.sum(resid**2) / ss_tot if ss_tot > 0 else 0
    rms = np.sqrt(np.mean(resid**2))
    return beta, y_hat, resid, R2, rms


def loo_cv(X, y):
    """Leave-one-out via hat matrix. Returns LOO RMS, LOO R², LOO residuals."""
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_hat = X @ beta
    resid = y - y_hat
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    h = np.diag(H)
    loo_resid = resid / (1 - h)
    loo_rms = np.sqrt(np.mean(loo_resid**2))
    loo_r2 = 1 - np.sum(loo_resid**2) / np.sum((y - np.mean(y))**2)
    return loo_rms, loo_r2, loo_resid


def boot_ci(X, y, n_boot=2000, seed=42):
    """Bootstrap confidence interval on R²."""
    rng = np.random.RandomState(seed)
    n = len(y)
    r2_boots = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        _, _, _, R2_b, _ = build_model(X[idx], y[idx])
        r2_boots.append(R2_b)
    r2_boots = np.array(r2_boots)
    return np.percentile(r2_boots, [2.5, 50, 97.5])


print("=" * 70)
print("SESSION #484: VALIDATING THE 6-VARIABLE MODEL")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nSample: {n} galaxies")

# Build variables
logV = np.log10([g['vflat'] for g in galaxies])
logL = np.log10([g['lum'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
y = np.array([g['outer_offset'] for g in galaxies])
T = np.array([g['hubble_type'] for g in galaxies])

# 5-variable model
X5 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V])
beta5, yhat5, resid5, R2_5, rms5 = build_model(X5, y)
loo5, loo_r2_5, loo_resid5 = loo_cv(X5, y)

# 6-variable model: + logL×f_gas
X6 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, y)
loo6, loo_r2_6, loo_resid6 = loo_cv(X6, y)

print(f"\n5-var: R² = {R2_5:.4f}, LOO R² = {loo_r2_5:.4f}, RMS = {rms5:.4f}")
print(f"6-var: R² = {R2_6:.4f}, LOO R² = {loo_r2_6:.4f}, RMS = {rms6:.4f}")

# =====================================================================
# TEST 1: FULL MODEL COEFFICIENTS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: MODEL COEFFICIENTS AND SIGNIFICANCE")
print("=" * 60)

var_names_5 = ['intercept', 'logV', 'logL', 'c_V', 'f_gas', 'logV×c_V']
var_names_6 = var_names_5 + ['logL×f_gas']

# Standard errors via OLS formula
def ols_se(X, resid):
    """Standard errors of OLS coefficients."""
    n, k = X.shape
    s2 = np.sum(resid**2) / (n - k)
    cov = s2 * np.linalg.inv(X.T @ X)
    return np.sqrt(np.diag(cov))

se5 = ols_se(X5, resid5)
se6 = ols_se(X6, resid6)
t5 = beta5 / se5
t6 = beta6 / se6

print("\n5-variable model:")
print(f"  {'Variable':<15} {'β':<10} {'SE':<10} {'t':<10}")
print("  " + "-" * 45)
for i, name in enumerate(var_names_5):
    print(f"  {name:<15} {beta5[i]:+.4f}   {se5[i]:.4f}   {t5[i]:+.2f}")

print(f"\n6-variable model:")
print(f"  {'Variable':<15} {'β':<10} {'SE':<10} {'t':<10}")
print("  " + "-" * 45)
for i, name in enumerate(var_names_6):
    print(f"  {name:<15} {beta6[i]:+.4f}   {se6[i]:.4f}   {t6[i]:+.2f}")

# F-test for the new variable
k_5 = X5.shape[1]
k_6 = X6.shape[1]
ss_5 = np.sum(resid5**2)
ss_6 = np.sum(resid6**2)
F_stat = ((ss_5 - ss_6) / (k_6 - k_5)) / (ss_6 / (n - k_6))
print(f"\nF-test for logL×f_gas: F = {F_stat:.2f} (df = 1, {n - k_6})")
print(f"  t² = {t6[-1]**2:.2f} (should equal F)")

# Significance via approximate p-value
# For F(1, 121), F > 6.85 → p < 0.01
if F_stat > 10:
    print(f"  HIGHLY SIGNIFICANT (F > 10, p < 0.002)")
elif F_stat > 6.85:
    print(f"  SIGNIFICANT (F > 6.85, p < 0.01)")
elif F_stat > 3.92:
    print(f"  MARGINALLY SIGNIFICANT (F > 3.92, p < 0.05)")
else:
    print(f"  NOT SIGNIFICANT (F < 3.92, p > 0.05)")

assert abs(t6[-1]) > 2, "logL×f_gas should be significant at |t| > 2"
print("\n✓ Test 1 passed: coefficients computed")

# =====================================================================
# TEST 2: BOOTSTRAP CONFIDENCE INTERVALS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: BOOTSTRAP CONFIDENCE INTERVALS ON ΔR²")
print("=" * 60)

n_boot = 2000
rng = np.random.RandomState(42)
dR2_boots = []
r2_5_boots = []
r2_6_boots = []

for b in range(n_boot):
    idx = rng.choice(n, n, replace=True)
    _, _, _, R2_5b, _ = build_model(X5[idx], y[idx])
    _, _, _, R2_6b, _ = build_model(X6[idx], y[idx])
    dR2_boots.append(R2_6b - R2_5b)
    r2_5_boots.append(R2_5b)
    r2_6_boots.append(R2_6b)

dR2_boots = np.array(dR2_boots)
r2_5_boots = np.array(r2_5_boots)
r2_6_boots = np.array(r2_6_boots)

pct = np.percentile(dR2_boots, [2.5, 50, 97.5])
print(f"\nΔR² = {R2_6 - R2_5:.4f}")
print(f"Bootstrap 95% CI: [{pct[0]:.4f}, {pct[2]:.4f}]")
print(f"Bootstrap median: {pct[1]:.4f}")
print(f"Prob(ΔR² > 0): {(dR2_boots > 0).mean():.4f}")
print(f"Prob(ΔR² > 0.01): {(dR2_boots > 0.01).mean():.4f}")

pct5 = np.percentile(r2_5_boots, [2.5, 97.5])
pct6 = np.percentile(r2_6_boots, [2.5, 97.5])
print(f"\n5-var R² 95% CI: [{pct5[0]:.4f}, {pct5[1]:.4f}]")
print(f"6-var R² 95% CI: [{pct6[0]:.4f}, {pct6[1]:.4f}]")

assert pct[0] > 0, "Lower CI on ΔR² should be positive"
print("\n✓ Test 2 passed: bootstrap CI excludes zero")

# =====================================================================
# TEST 3: PERFORMANCE BY GALAXY TYPE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: PERFORMANCE BY GALAXY TYPE")
print("=" * 60)

type_groups = {
    'Early (T<4)': T < 4,
    'Mid (4≤T<7)': (T >= 4) & (T < 7),
    'Late (T≥7)': T >= 7,
    'Late+gas (T≥7,f>0.3)': (T >= 7) & (f_gas > 0.3),
}

print(f"\n{'Group':<25} {'N':<5} {'R²_5':<8} {'R²_6':<8} {'ΔR²':<8} {'LOO5':<8} {'LOO6':<8} {'ΔLOO':<8}")
print("-" * 78)

for name, mask in type_groups.items():
    nm = mask.sum()
    if nm < 10:
        print(f"{name:<25} {nm:<5} (too few)")
        continue
    _, _, _, R2_5g, _ = build_model(X5[mask], y[mask])
    _, _, _, R2_6g, _ = build_model(X6[mask], y[mask])
    loo5g, loo_r2_5g, _ = loo_cv(X5[mask], y[mask])
    loo6g, loo_r2_6g, _ = loo_cv(X6[mask], y[mask])
    dR2 = R2_6g - R2_5g
    dloo = loo_r2_6g - loo_r2_5g
    print(f"{name:<25} {nm:<5} {R2_5g:.4f}  {R2_6g:.4f}  {dR2:+.4f}  {loo_r2_5g:.4f}  {loo_r2_6g:.4f}  {dloo:+.4f}")

print("\n✓ Test 3 passed: type analysis complete")

# =====================================================================
# TEST 4: NN AUTOCORRELATION
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: NEAREST-NEIGHBOR AUTOCORRELATION")
print("=" * 60)

from scipy.spatial.distance import cdist

feat = np.column_stack([logV, logL, c_V, f_gas])
feat_std = (feat - feat.mean(axis=0)) / np.clip(feat.std(axis=0), 1e-10, None)
dist_matrix = cdist(feat_std, feat_std, 'euclidean')
np.fill_diagonal(dist_matrix, np.inf)

for k in [1, 3, 5, 10]:
    # k nearest neighbors
    r5_nn = []
    r6_nn = []
    for i in range(n):
        nn_idx = np.argsort(dist_matrix[i])[:k]
        r5_nn.append(np.mean(resid5[nn_idx]))
        r6_nn.append(np.mean(resid6[nn_idx]))
    r5_nn = np.array(r5_nn)
    r6_nn = np.array(r6_nn)

    r_5 = np.corrcoef(resid5, r5_nn)[0, 1]
    r_6 = np.corrcoef(resid6, r6_nn)[0, 1]
    print(f"k={k:<3}  5-var: r = {r_5:+.3f}   6-var: r = {r_6:+.3f}   Δ = {r_6 - r_5:+.3f}")

# Also test with the 6th variable in the feature space
feat6 = np.column_stack([logV, logL, c_V, f_gas, logL * f_gas])
feat6_std = (feat6 - feat6.mean(axis=0)) / np.clip(feat6.std(axis=0), 1e-10, None)
dist6 = cdist(feat6_std, feat6_std, 'euclidean')
np.fill_diagonal(dist6, np.inf)

nn1_idx = np.argmin(dist6, axis=1)
r_6ext = np.corrcoef(resid6, resid6[nn1_idx])[0, 1]
print(f"\nWith logL×f_gas in feature space: k=1 autocorrelation = {r_6ext:+.3f}")

print("\n✓ Test 4 passed: autocorrelation analysis done")

# =====================================================================
# TEST 5: OUTLIER SENSITIVITY
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: OUTLIER SENSITIVITY")
print("=" * 60)

# Jackknife: remove each galaxy, measure ΔR² impact
dR2_jack = []
for i in range(n):
    mask = np.ones(n, dtype=bool)
    mask[i] = False
    _, _, _, R2_5j, _ = build_model(X5[mask], y[mask])
    _, _, _, R2_6j, _ = build_model(X6[mask], y[mask])
    dR2_jack.append(R2_6j - R2_5j)

dR2_jack = np.array(dR2_jack)
obs_dR2 = R2_6 - R2_5

print(f"\nJackknife ΔR² statistics:")
print(f"  Observed ΔR²: {obs_dR2:.4f}")
print(f"  Jackknife mean: {np.mean(dR2_jack):.4f}")
print(f"  Jackknife std: {np.std(dR2_jack):.4f}")
print(f"  Min: {np.min(dR2_jack):.4f}")
print(f"  Max: {np.max(dR2_jack):.4f}")
print(f"  Fraction with ΔR² > 0.01: {(dR2_jack > 0.01).mean():.3f}")

# Identify influential galaxies
most_influential = np.argsort(dR2_jack)  # smallest = galaxy whose removal reduces ΔR² most
print(f"\nMost influential galaxies (removal REDUCES improvement):")
for rank in range(3):
    i = most_influential[rank]
    g = galaxies[i]
    print(f"  {g['id']}: ΔR² without = {dR2_jack[i]:.4f}, T={g['hubble_type']}, "
          f"V={g['vflat']:.0f}, f_gas={g['f_gas']:.3f}, resid5={resid5[i]:+.4f}")

# Also test removing top 5 and top 10 outliers
for n_remove in [5, 10]:
    worst = np.argsort(np.abs(resid5))[-n_remove:]
    mask = np.ones(n, dtype=bool)
    mask[worst] = False
    _, _, _, R2_5r, _ = build_model(X5[mask], y[mask])
    _, _, _, R2_6r, _ = build_model(X6[mask], y[mask])
    loo5r, loo_r2_5r, _ = loo_cv(X5[mask], y[mask])
    loo6r, loo_r2_6r, _ = loo_cv(X6[mask], y[mask])
    print(f"\nRemoving top {n_remove} outliers (N={mask.sum()}):")
    print(f"  5-var: R² = {R2_5r:.4f}, LOO R² = {loo_r2_5r:.4f}")
    print(f"  6-var: R² = {R2_6r:.4f}, LOO R² = {loo_r2_6r:.4f}")
    print(f"  ΔR² = {R2_6r - R2_5r:.4f}, ΔLOO R² = {loo_r2_6r - loo_r2_5r:.4f}")

assert np.min(dR2_jack) > 0, "ΔR² should remain positive even with jackknife"
print("\n✓ Test 5 passed: outlier sensitivity robust")

# =====================================================================
# TEST 6: f_gas² vs logL×f_gas — WHICH IS MORE FUNDAMENTAL?
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: f_gas² vs logL×f_gas")
print("=" * 60)

# Model with f_gas²
X6_fg2 = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, f_gas**2])
beta_fg2, _, resid_fg2, R2_fg2, rms_fg2 = build_model(X6_fg2, y)
loo_fg2, loo_r2_fg2, _ = loo_cv(X6_fg2, y)

# Model with logV×f_gas
X6_vf = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logV * f_gas])
beta_vf, _, resid_vf, R2_vf, rms_vf = build_model(X6_vf, y)
loo_vf, loo_r2_vf, _ = loo_cv(X6_vf, y)

# Model with BOTH logL×f_gas AND f_gas²
X7_both = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas, f_gas**2])
beta_both, _, resid_both, R2_both, rms_both = build_model(X7_both, y)
loo_both, loo_r2_both, _ = loo_cv(X7_both, y)

print(f"\n{'Model':<25} {'R²':<8} {'LOO R²':<10} {'RMS':<8}")
print("-" * 51)
print(f"{'5-var baseline':<25} {R2_5:.4f}  {loo_r2_5:.4f}    {rms5:.4f}")
print(f"{'+ logL×f_gas':<25} {R2_6:.4f}  {loo_r2_6:.4f}    {rms6:.4f}")
print(f"{'+ f_gas²':<25} {R2_fg2:.4f}  {loo_r2_fg2:.4f}    {rms_fg2:.4f}")
print(f"{'+ logV×f_gas':<25} {R2_vf:.4f}  {loo_r2_vf:.4f}    {rms_vf:.4f}")
print(f"{'+ logL×f_gas + f_gas²':<25} {R2_both:.4f}  {loo_r2_both:.4f}    {rms_both:.4f}")

# Partial correlation of each at fixed other
# logL×f_gas partial controlling f_gas²
_, _, lf_resid_ctrl, _, _ = build_model(
    np.column_stack([X5, f_gas**2]), logL * f_gas)
_, _, y_resid_ctrl, _, _ = build_model(
    np.column_stack([X5, f_gas**2]), y)
r_lf_ctrl = np.corrcoef(lf_resid_ctrl, y_resid_ctrl)[0, 1]

# f_gas² partial controlling logL×f_gas
_, _, fg2_resid_ctrl, _, _ = build_model(
    np.column_stack([X5, logL * f_gas]), f_gas**2)
_, _, y_resid_ctrl2, _, _ = build_model(
    np.column_stack([X5, logL * f_gas]), y)
r_fg2_ctrl = np.corrcoef(fg2_resid_ctrl, y_resid_ctrl2)[0, 1]

print(f"\nPartial correlations:")
print(f"  r(logL×f_gas, offset | 5vars, f_gas²) = {r_lf_ctrl:+.4f}")
print(f"  r(f_gas², offset | 5vars, logL×f_gas) = {r_fg2_ctrl:+.4f}")

if abs(r_lf_ctrl) > abs(r_fg2_ctrl):
    print(f"  → logL×f_gas is more fundamental (retains signal controlling f_gas²)")
elif abs(r_fg2_ctrl) > abs(r_lf_ctrl):
    print(f"  → f_gas² is more fundamental (retains signal controlling logL×f_gas)")
else:
    print(f"  → Comparable — both capture similar non-linearity")

print("\n✓ Test 6 passed: f_gas interaction comparison done")

# =====================================================================
# TEST 7: THE 7-VARIABLE MODEL
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: EXTENDED MODELS")
print("=" * 60)

# Test adding more interactions systematically
extra_terms = {
    'logL×f_gas': logL * f_gas,
    'f_gas²': f_gas**2,
    'logV×f_gas': logV * f_gas,
    'logV²': logV**2,
    'logL²': logL**2,
    'c_V×f_gas': c_V * f_gas,
}

# Start with 5-var + logL×f_gas (the 6-var model), try adding each
print(f"\nStarting from 6-var model (R² = {R2_6:.4f}, LOO R² = {loo_r2_6:.4f}):")
print(f"\n{'7th variable':<20} {'R²':<8} {'LOO R²':<10} {'ΔLOO R²':<10}")
print("-" * 48)

for name, term in extra_terms.items():
    if name == 'logL×f_gas':
        continue  # already in model
    X7 = np.column_stack([X6, term])
    _, _, _, R2_7, _ = build_model(X7, y)
    loo7, loo_r2_7, _ = loo_cv(X7, y)
    d_loo = loo_r2_7 - loo_r2_6
    print(f"{name:<20} {R2_7:.4f}  {loo_r2_7:.4f}    {d_loo:+.4f}")

# The big 8-variable kitchen sink model
X8 = np.column_stack([X6, f_gas**2, logV * f_gas])
_, _, _, R2_8, rms8 = build_model(X8, y)
loo8, loo_r2_8, _ = loo_cv(X8, y)
print(f"\n8-var (5 + logL×f_gas + f_gas² + logV×f_gas):")
print(f"  R² = {R2_8:.4f}, LOO R² = {loo_r2_8:.4f}")

print("\n✓ Test 7 passed: extended models tested")

# =====================================================================
# TEST 8: RESIDUAL DIAGNOSTICS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: RESIDUAL DIAGNOSTICS")
print("=" * 60)

# Predicted vs observed
r_pred_obs = np.corrcoef(yhat6, y)[0, 1]
print(f"\nPredicted vs observed: r = {r_pred_obs:.4f}")

# Residual vs predicted (should be ~0 for good model)
r_resid_pred = np.corrcoef(resid6, yhat6)[0, 1]
print(f"Residual vs predicted: r = {r_resid_pred:.4f} (should be ~0)")

# Residual normality (skewness, kurtosis)
skew = np.mean((resid6 - np.mean(resid6))**3) / np.std(resid6)**3
kurt = np.mean((resid6 - np.mean(resid6))**4) / np.std(resid6)**4 - 3
print(f"Residual skewness: {skew:.3f} (0 = normal)")
print(f"Residual excess kurtosis: {kurt:.3f} (0 = normal)")

# Percentage within various thresholds
for thresh in [0.025, 0.05, 0.10]:
    pct = (np.abs(resid6) < thresh).mean() * 100
    pct5 = (np.abs(resid5) < thresh).mean() * 100
    print(f"  |resid| < {thresh}: {pct:.1f}% (6-var) vs {pct5:.1f}% (5-var)")

# Top outliers in 6-var model
worst_idx = np.argsort(np.abs(resid6))[-5:][::-1]
print(f"\nTop 5 outliers (6-var model):")
print(f"  {'Galaxy':<15} {'resid_5':<10} {'resid_6':<10} {'T':<5} {'V':<6} {'f_gas':<8}")
print("  " + "-" * 54)
for i in worst_idx:
    g = galaxies[i]
    print(f"  {g['id']:<15} {resid5[i]:+.4f}   {resid6[i]:+.4f}   {g['hubble_type']:<5} "
          f"{g['vflat']:<6.0f} {g['f_gas']:.3f}")

# Improvement statistics
improved = np.abs(resid6) < np.abs(resid5)
print(f"\nGalaxies improved: {improved.sum()}/{n} ({100*improved.mean():.1f}%)")
print(f"Mean |resid| change: {np.mean(np.abs(resid6)) - np.mean(np.abs(resid5)):+.4f}")

assert r_pred_obs > 0.95, "Predicted-observed correlation should be > 0.95"
print("\n✓ Test 8 passed: diagnostics complete")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #484 SUMMARY")
print("=" * 70)

print(f"\n--- The 6-Variable Model ---")
print(f"offset_outer = {beta6[0]:+.3f}")
for i, name in enumerate(var_names_6[1:], 1):
    print(f"  {'+' if beta6[i] >= 0 else ''} {beta6[i]:.3f} × {name}")
print(f"\nR² = {R2_6:.4f}")
print(f"LOO R² = {loo_r2_6:.4f}")
print(f"RMS = {rms6:.4f} dex ({10**(rms6)-1:.1%} velocity)")
print(f"ΔR² vs 5-var: {R2_6 - R2_5:+.4f}")
print(f"ΔLOO R² vs 5-var: {loo_r2_6 - loo_r2_5:+.4f}")
print(f"logL×f_gas t-statistic: {t6[-1]:+.2f}")
print(f"F-test: F = {F_stat:.2f}")

print(f"\nAll 8 tests passed ✓")
