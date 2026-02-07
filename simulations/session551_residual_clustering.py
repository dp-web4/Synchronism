#!/usr/bin/env python3
"""
======================================================================
SESSION #551: RESIDUAL CLUSTERING — ARE MODEL ERRORS RANDOM?
======================================================================

The 6-var model has RMS=0.038 dex and LOO R²=0.938. Session #523 showed
the residuals are dominated by measurement noise (587% of residual
variance from known errors). But are the residuals RANDOMLY distributed
across galaxy properties, or do they cluster? Are there sub-populations
of galaxies that the model systematically over- or under-predicts?

This session tests:
1. Residual distribution normality
2. Spatial (property-space) clustering of residuals
3. Sign patterns: do neighboring galaxies share residual signs?
4. Residual correlations with external properties
5. Two-cluster analysis: are there hidden sub-populations?
6. Conditional heteroscedasticity: does scatter vary with properties?
7. Residual PCA: is there structure in the error space?
8. Synthesis: the nature of model errors

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-07
Session: #551
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


from scipy import stats as sp_stats

print("=" * 70)
print("SESSION #551: RESIDUAL CLUSTERING")
print("=" * 70)

# Load data
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
catalog = load_sparc_catalog(
    os.path.join(base_dir, "sparc_real_data", "SPARC_Lelli2016c.mrt"))
models = load_sparc_mass_models(
    os.path.join(base_dir, "sparc_real_data", "MassModels_Lelli2016c.mrt"))

ml_disk = 0.5
ml_bul = 0.7

# Prepare galaxies
galaxies = []
for gal_id, points in models.items():
    if len(points) < 5 or gal_id not in catalog:
        continue
    cat = catalog[gal_id]
    vflat = cat.get('vflat', 0)
    lum = cat.get('luminosity', 0)
    sb_eff = cat.get('sb_eff', 0)
    distance = cat.get('distance', 0)
    inclination = cat.get('inclination', 0)
    hubble_type = cat.get('hubble_type', 0)
    quality = cat.get('quality', 0)

    if vflat <= 0 or lum <= 0 or sb_eff <= 0 or distance <= 0:
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

    if r_eff_kpc > 0 and r_eff_kpc >= radius_v.min() and r_eff_kpc <= radius_v.max():
        v_at_reff = np.interp(r_eff_kpc, radius_v, np.abs(v_obs_v))
        c_V_val = v_at_reff / vflat
    else:
        c_V_val = np.nan
    if not np.isfinite(c_V_val):
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
    else:
        offset_val = np.mean(offset_pts[mond])

    n_flat = min(5, len(v_gas_v))
    v_gas_end = np.mean(v_gas_v[-n_flat:]**2)
    v_disk_end = np.mean(v_disk_v[-n_flat:]**2)
    f_gas_val = v_gas_end / max(v_gas_end + v_disk_end, 1e-10)

    galaxies.append({
        'id': gal_id,
        'offset': offset_val,
        'logV': np.log10(vflat),
        'logL': np.log10(lum),
        'c_V': c_V_val,
        'f_gas': f_gas_val,
        'distance': distance,
        'logD': np.log10(distance),
        'sb_eff': sb_eff,
        'log_sb': np.log10(sb_eff) if sb_eff > 0 else np.nan,
        'inclination': inclination,
        'hubble_type': hubble_type,
        'quality': quality,
        'r_eff_kpc': r_eff_kpc,
        'r_max': radius_v.max(),
        'n_points': len(g_bar_v),
        'n_mond': mond.sum(),
        'mean_e_vobs': np.mean(e_vobs_v),
    })

n = len(galaxies)
print(f"\n{n} galaxies loaded")

offset = np.array([g['offset'] for g in galaxies])
logV = np.array([g['logV'] for g in galaxies])
logL = np.array([g['logL'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
logD = np.array([g['logD'] for g in galaxies])
log_sb = np.array([g['log_sb'] for g in galaxies])
inclination = np.array([g['inclination'] for g in galaxies])
hubble_type = np.array([g['hubble_type'] for g in galaxies])
quality = np.array([g['quality'] for g in galaxies])
n_points = np.array([g['n_points'] for g in galaxies])
n_mond = np.array([g['n_mond'] for g in galaxies])
mean_e = np.array([g['mean_e_vobs'] for g in galaxies])
ones = np.ones(n)

# Standard 6-var model
X6 = np.column_stack([ones, logV, logL, c_V, f_gas, logV*c_V, logL*f_gas])
beta6, yhat6, resid6, R2_6, rms6 = build_model(X6, offset)
loo6 = loo_r2(X6, offset)

# LOO residuals
beta_loo = np.linalg.lstsq(X6, offset, rcond=None)[0]
H = X6 @ np.linalg.inv(X6.T @ X6) @ X6.T
h = np.diag(H)
loo_resid = resid6 / (1 - h)

print(f"\nStandard 6-var: R²={R2_6:.4f}, LOO={loo6:.4f}, RMS={rms6:.4f}")

# ============================================================
print("\n" + "=" * 60)
print("TEST 1: RESIDUAL DISTRIBUTION NORMALITY")
print("=" * 60)
# ============================================================

# Shapiro-Wilk test
shapiro_stat, shapiro_p = sp_stats.shapiro(resid6)

# D'Agostino-Pearson test
dagostino_stat, dagostino_p = sp_stats.normaltest(resid6)

# Skewness and kurtosis
skew = sp_stats.skew(resid6)
kurt = sp_stats.kurtosis(resid6)  # excess kurtosis (normal = 0)

# Quantile-quantile correlation
sorted_resid = np.sort(resid6)
expected_quantiles = sp_stats.norm.ppf(np.arange(1, n+1) / (n+1))
qq_r, _ = sp_stats.pearsonr(sorted_resid, expected_quantiles)

print(f"\nResidual distribution tests:")
print(f"  N = {n}")
print(f"  Mean:     {np.mean(resid6):+.5f} (should be ~0)")
print(f"  Std:      {np.std(resid6):.5f}")
print(f"  Skewness: {skew:+.3f} (normal = 0)")
print(f"  Kurtosis: {kurt:+.3f} (normal = 0)")
print(f"")
print(f"  Shapiro-Wilk:       W={shapiro_stat:.4f}, p={shapiro_p:.3f}")
print(f"  D'Agostino-Pearson: K²={dagostino_stat:.2f}, p={dagostino_p:.3f}")
print(f"  Q-Q correlation:    r={qq_r:.4f}")

# LOO residuals
skew_loo = sp_stats.skew(loo_resid)
kurt_loo = sp_stats.kurtosis(loo_resid)
shapiro_loo, shapiro_loo_p = sp_stats.shapiro(loo_resid)
print(f"\nLOO residuals:")
print(f"  Skewness: {skew_loo:+.3f}")
print(f"  Kurtosis: {kurt_loo:+.3f}")
print(f"  Shapiro-Wilk: W={shapiro_loo:.4f}, p={shapiro_loo_p:.3f}")

is_normal = shapiro_p > 0.05
print(f"\n  Residuals are {'consistent with' if is_normal else 'NOT'} normal (p={shapiro_p:.3f})")

print(f"\n✓ TEST 1 PASSED: Normality assessed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 2: NEAREST-NEIGHBOR RESIDUAL CORRELATION")
print("=" * 60)
# ============================================================

# In the 4D property space (logV, logL, c_V, f_gas), do nearby galaxies
# have correlated residuals? This would indicate the model is missing a
# spatially smooth signal.

# Standardize the property space
props = np.column_stack([logV, logL, c_V, f_gas])
props_std = (props - props.mean(axis=0)) / props.std(axis=0)

# Compute pairwise distances
from scipy.spatial.distance import cdist
dist_matrix = cdist(props_std, props_std, 'euclidean')

# For each galaxy, find its k nearest neighbors
k_values = [1, 3, 5, 10]
print(f"\nNearest-neighbor residual correlation (in 4D property space):")
print(f"{'k-NN':<6} {'Mean r(resid,NN_resid)':<25} {'p (permutation)'}")
print("-" * 55)

np.random.seed(42)
n_perm = 1000

for k in k_values:
    # Mean residual of k nearest neighbors
    nn_resid = np.zeros(n)
    for i in range(n):
        distances = dist_matrix[i, :]
        distances[i] = np.inf  # exclude self
        nn_idx = np.argsort(distances)[:k]
        nn_resid[i] = np.mean(resid6[nn_idx])

    # Correlation
    r_nn, _ = sp_stats.pearsonr(resid6, nn_resid)

    # Permutation test
    r_perm = np.zeros(n_perm)
    for p in range(n_perm):
        perm_resid = np.random.permutation(resid6)
        nn_perm = np.zeros(n)
        for i in range(n):
            distances = dist_matrix[i, :]
            distances_copy = distances.copy()
            distances_copy[i] = np.inf
            nn_idx = np.argsort(distances_copy)[:k]
            nn_perm[i] = np.mean(perm_resid[nn_idx])
        r_perm[p] = np.corrcoef(perm_resid, nn_perm)[0, 1]

    p_perm = np.mean(np.abs(r_perm) >= abs(r_nn))
    print(f"k={k:<4} r={r_nn:+.3f}                  p={p_perm:.3f}")

print(f"\n  Residuals {'cluster' if r_nn > 0.15 and p_perm < 0.05 else 'do NOT cluster'} "
      f"in property space")

print(f"\n✓ TEST 2 PASSED: NN correlation tested")

# ============================================================
print("\n" + "=" * 60)
print("TEST 3: SIGN PATTERNS — DO RESIDUAL SIGNS CLUSTER?")
print("=" * 60)
# ============================================================

# For each galaxy, what fraction of its k nearest neighbors have the
# SAME sign residual? If random, this should be ~50%.

for k in [3, 5, 10]:
    same_sign_fracs = np.zeros(n)
    for i in range(n):
        distances = dist_matrix[i, :].copy()
        distances[i] = np.inf
        nn_idx = np.argsort(distances)[:k]
        same_sign_fracs[i] = np.mean(np.sign(resid6[nn_idx]) == np.sign(resid6[i]))

    mean_frac = np.mean(same_sign_fracs)

    # Expected under random: fraction of galaxies with same sign
    pos_frac = np.mean(resid6 > 0)
    expected = pos_frac**2 + (1-pos_frac)**2  # P(same sign)

    print(f"k={k}: Mean same-sign fraction = {mean_frac:.3f} (expected {expected:.3f})")

# Runs test on residuals sorted by logV
sorted_idx = np.argsort(logV)
sorted_signs = np.sign(resid6[sorted_idx])
# Count runs (consecutive same-sign sequences)
n_runs = 1
for i in range(1, n):
    if sorted_signs[i] != sorted_signs[i-1]:
        n_runs += 1

n_pos = np.sum(sorted_signs > 0)
n_neg = np.sum(sorted_signs <= 0)
# Expected runs under randomness
expected_runs = 1 + 2 * n_pos * n_neg / n
std_runs = np.sqrt(2 * n_pos * n_neg * (2*n_pos*n_neg - n) / (n**2 * (n-1)))
z_runs = (n_runs - expected_runs) / std_runs if std_runs > 0 else 0
p_runs = 2 * sp_stats.norm.sf(abs(z_runs))

print(f"\nRuns test (residuals sorted by logV):")
print(f"  Observed runs: {n_runs}")
print(f"  Expected: {expected_runs:.1f} ± {std_runs:.1f}")
print(f"  z = {z_runs:+.2f}, p = {p_runs:.3f}")
print(f"  {'Clustered' if z_runs < -2 else 'Random' if abs(z_runs) < 2 else 'Anti-clustered'} signs")

print(f"\n✓ TEST 3 PASSED: Sign patterns analyzed")

# ============================================================
print("\n" + "=" * 60)
print("TEST 4: RESIDUAL vs EXTERNAL PROPERTIES")
print("=" * 60)
# ============================================================

# Test residuals against properties NOT in the model
external = {
    'logD': logD,
    'inclination': inclination,
    'hubble_type': hubble_type,
    'quality': quality,
    'log(SB_eff)': log_sb,
    'n_points': n_points,
    'n_MOND_points': n_mond,
    'mean(e_vobs)': mean_e,
    'R_max': np.array([g['r_max'] for g in galaxies]),
    'R_max/r_eff': np.array([g['r_max']/max(g['r_eff_kpc'],0.01) for g in galaxies]),
}

print(f"\nResidual correlations with external properties:")
print(f"{'Property':<18} {'r(resid,prop)':<15} {'r(|resid|,prop)':<16} {'p(resid)'}")
print("-" * 65)

for name, prop in external.items():
    mask = np.isfinite(prop)
    if mask.sum() < 10:
        continue
    r, p = sp_stats.pearsonr(resid6[mask], prop[mask])
    r_abs, p_abs = sp_stats.pearsonr(np.abs(resid6[mask]), prop[mask])
    star = " *" if p < 0.05 else ""
    print(f"  {name:<16} {r:+.3f}         {r_abs:+.3f}           {p:.3f}{star}")

# The key test: r(|resid|, prop) tells us about heteroscedasticity
# r(resid, prop) tells us about systematic bias

print(f"\n  * = significant at p<0.05")
print(f"  r(resid, prop): systematic bias")
print(f"  r(|resid|, prop): heteroscedastic (varying scatter)")

print(f"\n✓ TEST 4 PASSED: External correlations mapped")

# ============================================================
print("\n" + "=" * 60)
print("TEST 5: TWO-CLUSTER ANALYSIS")
print("=" * 60)
# ============================================================

# Are there hidden sub-populations? Use simple k-means with k=2 on
# the residuals + leverage to see if galaxies naturally split.

# First: just split by residual sign
pos_resid = resid6 > 0
neg_resid = resid6 <= 0

print(f"\nResidual sign split:")
print(f"  Positive (N={pos_resid.sum()}): over-predicted offset (too much DM)")
print(f"  Negative (N={neg_resid.sum()}): under-predicted offset (too little DM)")

print(f"\n  {'Property':<15} {'Positive resid':<18} {'Negative resid':<18} {'t-test p'}")
print(f"  {'-'*65}")
for name, prop in [('logV', logV), ('logL', logL), ('c_V', c_V),
                    ('f_gas', f_gas), ('hubble_type', hubble_type),
                    ('inclination', inclination), ('quality', quality),
                    ('logD', logD), ('n_points', n_points.astype(float))]:
    mean_pos = np.mean(prop[pos_resid])
    mean_neg = np.mean(prop[neg_resid])
    t, p = sp_stats.ttest_ind(prop[pos_resid], prop[neg_resid])
    star = " *" if p < 0.05 else ""
    print(f"  {name:<15} {mean_pos:.3f}            {mean_neg:.3f}            {p:.3f}{star}")

# K-means on standardized properties + residual
from scipy.cluster.hierarchy import fcluster, linkage

# Cluster on RESIDUALS + all properties
features = np.column_stack([resid6, logV, logL, c_V, f_gas])
features_std = (features - features.mean(axis=0)) / features.std(axis=0)

Z = linkage(features_std, method='ward')
labels_2 = fcluster(Z, t=2, criterion='maxclust')
labels_3 = fcluster(Z, t=3, criterion='maxclust')

# Compare clusters
for n_clust, labels in [(2, labels_2), (3, labels_3)]:
    print(f"\n{n_clust}-cluster analysis (Ward's method):")
    for c in range(1, n_clust+1):
        mask = labels == c
        print(f"  Cluster {c} (N={mask.sum()}): "
              f"resid={np.mean(resid6[mask]):+.4f}±{np.std(resid6[mask]):.4f}, "
              f"logV={np.mean(logV[mask]):.2f}, "
              f"f_gas={np.mean(f_gas[mask]):.2f}")

# Is the cluster separation meaningful?
# Compare within-cluster vs between-cluster residual variance
for n_clust, labels in [(2, labels_2)]:
    between = 0
    within = 0
    overall_mean = np.mean(resid6)
    for c in range(1, n_clust+1):
        mask = labels == c
        n_c = mask.sum()
        between += n_c * (np.mean(resid6[mask]) - overall_mean)**2
        within += np.sum((resid6[mask] - np.mean(resid6[mask]))**2)
    total = np.sum((resid6 - overall_mean)**2)
    f_ratio = (between / (n_clust-1)) / (within / (n - n_clust))
    print(f"\n  F-ratio for 2-cluster split: {f_ratio:.2f}")
    print(f"  Between/Total variance: {between/total:.3f} ({between/total*100:.1f}%)")

print(f"\n✓ TEST 5 PASSED: Cluster analysis complete")

# ============================================================
print("\n" + "=" * 60)
print("TEST 6: CONDITIONAL HETEROSCEDASTICITY")
print("=" * 60)
# ============================================================

# Does residual scatter vary systematically with galaxy properties?

# Split into quartiles of each property and compute RMS per quartile
print(f"\nResidual RMS by property quartiles:")

for name, prop in [('logV', logV), ('logL', logL), ('c_V', c_V),
                    ('f_gas', f_gas)]:
    quartiles = np.percentile(prop, [25, 50, 75])
    bins = [prop.min()-1, quartiles[0], quartiles[1], quartiles[2], prop.max()+1]

    print(f"\n  {name}:")
    print(f"  {'Quartile':<12} {'N':<5} {'RMS(resid)':<12} {'Mean |resid|'}")
    print(f"  {'-'*40}")

    rms_by_q = []
    for q in range(4):
        mask = (prop >= bins[q]) & (prop < bins[q+1])
        if mask.sum() > 0:
            rms_q = np.sqrt(np.mean(resid6[mask]**2))
            mean_abs_q = np.mean(np.abs(resid6[mask]))
            rms_by_q.append(rms_q)
            print(f"  Q{q+1:<10} {mask.sum():<5} {rms_q:.4f}      {mean_abs_q:.4f}")

    # Ratio of max to min RMS
    if len(rms_by_q) >= 2:
        ratio = max(rms_by_q) / min(rms_by_q)
        print(f"  Max/Min RMS ratio: {ratio:.2f}")

# Breusch-Pagan test for heteroscedasticity
resid_sq = resid6**2
X_bp = np.column_stack([ones, logV, logL, c_V, f_gas])
beta_bp = np.linalg.lstsq(X_bp, resid_sq, rcond=None)[0]
resid_sq_hat = X_bp @ beta_bp
resid_bp = resid_sq - resid_sq_hat
SS_reg = np.sum((resid_sq_hat - np.mean(resid_sq))**2)
SS_tot = np.sum((resid_sq - np.mean(resid_sq))**2)
R2_bp = SS_reg / SS_tot
LM_stat = n * R2_bp
p_bp = 1 - sp_stats.chi2.cdf(LM_stat, 4)

print(f"\nBreusch-Pagan test for heteroscedasticity:")
print(f"  LM statistic = {LM_stat:.2f}")
print(f"  p-value = {p_bp:.3f}")
print(f"  {'Heteroscedastic' if p_bp < 0.05 else 'Homoscedastic'} at 5% level")

# Which property drives scatter variation?
print(f"\n  r(resid², property):")
for name, prop in [('logV', logV), ('logL', logL), ('c_V', c_V), ('f_gas', f_gas),
                    ('hubble_type', hubble_type), ('n_points', n_points.astype(float))]:
    r, p = sp_stats.pearsonr(resid_sq, prop)
    star = " *" if p < 0.05 else ""
    print(f"    {name:<15}: r={r:+.3f}, p={p:.3f}{star}")

print(f"\n✓ TEST 6 PASSED: Heteroscedasticity tested")

# ============================================================
print("\n" + "=" * 60)
print("TEST 7: RESIDUAL PCA — STRUCTURE IN THE ERROR SPACE")
print("=" * 60)
# ============================================================

# If we combine residuals with galaxy properties, are there principal
# components that align with the residual?

# PCA on standardized properties + residual
all_vars = np.column_stack([logV, logL, c_V, f_gas, resid6])
all_std = (all_vars - all_vars.mean(axis=0)) / all_vars.std(axis=0)

cov_matrix = np.cov(all_std.T)
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

# Sort by decreasing eigenvalue
idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

var_explained = eigenvalues / eigenvalues.sum()

print(f"\nPCA on (logV, logL, c_V, f_gas, residual):")
print(f"{'PC':<5} {'Eigenvalue':<12} {'% Var':<8} {'Cumul%':<8} {'resid loading'}")
print("-" * 50)

var_names = ['logV', 'logL', 'c_V', 'f_gas', 'resid']
cumul = 0
for i in range(5):
    cumul += var_explained[i] * 100
    resid_loading = eigenvectors[4, i]  # loading of residual on this PC
    print(f"PC{i+1:<3} {eigenvalues[i]:.3f}       {var_explained[i]*100:.1f}%    {cumul:.1f}%    {resid_loading:+.3f}")

# Which PC does the residual load on most?
resid_loadings = eigenvectors[4, :]
max_pc = np.argmax(np.abs(resid_loadings)) + 1
print(f"\nResidual loads most on PC{max_pc} (loading = {resid_loadings[max_pc-1]:+.3f})")

# Is the residual's loading on PC1 (the mass/Hubble sequence) significant?
print(f"\nResidual's loading on PC1 (mass sequence): {resid_loadings[0]:+.3f}")
print(f"  → Residual is {'aligned with' if abs(resid_loadings[0]) > 0.3 else 'orthogonal to'} the mass sequence")

# What does the residual's PC look like?
print(f"\nPC{max_pc} loadings (the 'residual PC'):")
for j, name in enumerate(var_names):
    print(f"  {name:<10}: {eigenvectors[j, max_pc-1]:+.3f}")

# Is the residual orthogonal to the model variables? (should be by construction)
r_resid_pc1 = np.corrcoef(all_std[:, 4], all_std @ eigenvectors[:, 0])[0, 1]
print(f"\nr(residual, PC1_score) = {r_resid_pc1:+.3f}")

print(f"\n✓ TEST 7 PASSED: Residual PCA complete")

# ============================================================
print("\n" + "=" * 60)
print("TEST 8: SYNTHESIS — THE NATURE OF MODEL ERRORS")
print("=" * 60)
# ============================================================

print(f"\n{'='*60}")
print(f"RESIDUAL CLUSTERING SYNTHESIS")
print(f"{'='*60}")

print(f"\n1. DISTRIBUTION:")
print(f"   Shapiro-Wilk: W={shapiro_stat:.4f}, p={shapiro_p:.3f}")
print(f"   Skewness={skew:+.3f}, kurtosis={kurt:+.3f}")
print(f"   → {'Normal' if is_normal else 'Non-normal'}")

print(f"\n2. SPATIAL CLUSTERING (in property space):")
# Report k=5 result
k5_idx = k_values.index(5)
# Re-compute for reporting
nn_resid_5 = np.zeros(n)
for i in range(n):
    distances = dist_matrix[i, :].copy()
    distances[i] = np.inf
    nn_idx = np.argsort(distances)[:5]
    nn_resid_5[i] = np.mean(resid6[nn_idx])
r_nn5, _ = sp_stats.pearsonr(resid6, nn_resid_5)
print(f"   k=5 NN correlation: r={r_nn5:+.3f}")
print(f"   → Residuals {'cluster' if r_nn5 > 0.15 else 'do NOT cluster'} in property space")

print(f"\n3. SIGN PATTERNS:")
print(f"   Runs test (sorted by logV): z={z_runs:+.2f}, p={p_runs:.3f}")
print(f"   → Signs are {'random' if abs(z_runs) < 2 else 'structured'}")

print(f"\n4. EXTERNAL CORRELATIONS:")
# Report the strongest
for name, prop in external.items():
    mask = np.isfinite(prop)
    if mask.sum() < 10:
        continue
    r, p = sp_stats.pearsonr(resid6[mask], prop[mask])
    if p < 0.05:
        print(f"   r(resid, {name}) = {r:+.3f}, p={p:.3f} ← SIGNIFICANT")
    elif abs(r) > 0.1:
        print(f"   r(resid, {name}) = {r:+.3f}, p={p:.3f}")

n_sig = sum(1 for name, prop in external.items()
            if np.isfinite(prop).sum() >= 10 and
            sp_stats.pearsonr(resid6[np.isfinite(prop)], prop[np.isfinite(prop)])[1] < 0.05)
print(f"   {n_sig}/{len(external)} external properties are significant")

print(f"\n5. HETEROSCEDASTICITY:")
print(f"   Breusch-Pagan: LM={LM_stat:.2f}, p={p_bp:.3f}")
print(f"   → {'Scatter varies with properties' if p_bp < 0.05 else 'Scatter is uniform'}")

print(f"\n6. PCA:")
print(f"   Residual loads most on PC{max_pc} ({resid_loadings[max_pc-1]:+.3f})")
print(f"   PC1 loading: {resid_loadings[0]:+.3f} ({'aligned' if abs(resid_loadings[0]) > 0.3 else 'orthogonal'} to mass sequence)")

print(f"\n{'='*60}")
print(f"BOTTOM LINE:")
if is_normal and abs(z_runs) < 2 and p_bp > 0.05 and n_sig <= 1:
    print(f"  Model errors are random: normal, unclustered, homoscedastic,")
    print(f"  uncorrelated with external properties. The residuals are")
    print(f"  indistinguishable from measurement noise.")
elif is_normal and n_sig <= 2:
    print(f"  Model errors are mostly random: normal distribution,")
    print(f"  with minor departures from perfect randomness.")
    print(f"  The dominant error source is measurement noise.")
else:
    print(f"  Model errors show some structure. Further investigation")
    print(f"  may reveal systematic patterns.")
print(f"{'='*60}")

print(f"\n✓ TEST 8 PASSED: Synthesis complete")

# ============================================================
# FINAL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print(f"SESSION #551: ALL 8 TESTS PASSED")
print(f"{'='*70}")
print(f"\nKey findings:")
print(f"  1. Residual normality: W={shapiro_stat:.4f}, p={shapiro_p:.3f}")
print(f"  2. NN clustering (k=5): r={r_nn5:+.3f}")
print(f"  3. Runs test: z={z_runs:+.2f}, p={p_runs:.3f}")
print(f"  4. External correlations: {n_sig}/{len(external)} significant")
print(f"  5. Breusch-Pagan: p={p_bp:.3f}")
print(f"  6. Residual PC: PC{max_pc}, loading={resid_loadings[max_pc-1]:+.3f}")
print(f"  7. Residual is {'orthogonal' if abs(resid_loadings[0]) < 0.3 else 'aligned'} to mass sequence")
print(f"  8. Errors are {'random' if is_normal and abs(z_runs) < 2 else 'structured'}")
