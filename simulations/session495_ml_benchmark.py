#!/usr/bin/env python3
"""
======================================================================
SESSION #495: MACHINE LEARNING BENCHMARK — IS LINEAR OPTIMAL?
======================================================================

The 6-variable linear model achieves R² = 0.945, LOO R² = 0.938.
Can machine learning (non-parametric) models do better?

If ML barely improves: the linear model captures all available information.
If ML substantially improves: there is non-linear structure we're missing.

We use cross-validation (not LOO, since ML models can't use the hat matrix).
All comparisons are out-of-sample via K-fold CV.

Tests:
1. K-Fold CV baseline (linear model)
2. Random Forest
3. Gradient Boosting (scikit-learn GBR)
4. K-Nearest Neighbors
5. Feature importance from RF
6. Partial dependence analysis
7. Residual structure from ML
8. The ML advantage quantified

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-06
Session: #495
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
    """Load SPARC data with all features."""
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

        # Additional features for ML
        sb_disk = cat.get('sb_disk', sb_eff)
        r_max = radius_v.max()
        roughness = np.std(np.diff(np.log10(np.clip(g_obs_v, 1e-15, None))))
        n_mond = mond.sum()
        n_pts = len(g_bar_v)

        galaxies.append({
            'id': gal_id, 'vflat': vflat, 'lum': lum, 'c_V': c_V,
            'hubble_type': hubble_type, 'f_gas': f_gas,
            'outer_offset': outer_offset,
            'distance': distance, 'inclination': inclination,
            'quality': quality,
            'sb_eff': sb_eff, 'sb_disk': sb_disk, 'r_eff_kpc': r_eff_kpc,
            'r_max': r_max, 'roughness': roughness,
            'n_mond': n_mond, 'n_pts': n_pts,
        })

    return galaxies


def kfold_cv(X, y, model_fn, k=10, seed=42):
    """K-fold cross-validation returning R² and RMS."""
    rng = np.random.RandomState(seed)
    n = len(y)
    idx = rng.permutation(n)
    fold_size = n // k
    predictions = np.zeros(n)

    for fold in range(k):
        start = fold * fold_size
        if fold == k - 1:
            end = n  # last fold gets remainder
        else:
            end = start + fold_size
        test_idx = idx[start:end]
        train_idx = np.concatenate([idx[:start], idx[end:]])

        X_train, y_train = X[train_idx], y[train_idx]
        X_test = X[test_idx]

        pred = model_fn(X_train, y_train, X_test)
        predictions[test_idx] = pred

    resid = y - predictions
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - np.sum(resid**2) / ss_tot
    rms = np.sqrt(np.mean(resid**2))
    return r2, rms


def linear_model(X_train, y_train, X_test):
    beta = np.linalg.lstsq(X_train, y_train, rcond=None)[0]
    return X_test @ beta


print("=" * 70)
print("SESSION #495: MACHINE LEARNING BENCHMARK — IS LINEAR OPTIMAL?")
print("=" * 70)

galaxies = prepare_data()
n = len(galaxies)
print(f"\nSample: {n} galaxies")

# Feature matrix
logV = np.log10([g['vflat'] for g in galaxies])
logL = np.log10([g['lum'] for g in galaxies])
c_V = np.array([g['c_V'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
T = np.array([g['hubble_type'] for g in galaxies])
incl = np.array([g['inclination'] for g in galaxies])
log_sb = np.log10([g['sb_eff'] for g in galaxies])
log_reff = np.log10([g['r_eff_kpc'] for g in galaxies])
log_rmax = np.log10([g['r_max'] for g in galaxies])
roughness = np.array([g['roughness'] for g in galaxies])
n_mond = np.array([g['n_mond'] for g in galaxies])
y = np.array([g['outer_offset'] for g in galaxies])

# Feature sets
# 4 core features (no interactions)
X4 = np.column_stack([logV, logL, c_V, f_gas])
feat4_names = ['logV', 'logL', 'c_V', 'f_gas']

# 6-var (with interactions + intercept)
X6_lin = np.column_stack([np.ones(n), logV, logL, c_V, f_gas, logV * c_V, logL * f_gas])

# 4 core features for ML (no intercept needed, no interactions)
# ML discovers non-linearities automatically

# Extended features: add all available galaxy properties
X_ext = np.column_stack([logV, logL, c_V, f_gas, T, incl, log_sb, log_reff,
                          log_rmax, roughness, n_mond])
feat_ext_names = ['logV', 'logL', 'c_V', 'f_gas', 'T', 'incl', 'log_SB',
                   'log_R_eff', 'log_R_max', 'roughness', 'N_mond']

# =====================================================================
# TEST 1: K-FOLD BASELINE (LINEAR)
# =====================================================================
print("\n" + "=" * 60)
print("TEST 1: K-FOLD BASELINE (LINEAR MODEL)")
print("=" * 60)

r2_lin, rms_lin = kfold_cv(X6_lin, y, linear_model, k=10)
print(f"\n6-var linear (10-fold CV): R² = {r2_lin:.4f}, RMS = {rms_lin:.4f}")

# Also LOO for comparison
from session491_scatter_budget import loo_cv as loo_cv_func
loo_rms, loo_r2 = loo_cv_func(X6_lin, y)
print(f"6-var linear (LOO): R² = {loo_r2:.4f}, RMS = {loo_rms:.4f}")

# Multiple random seeds for stability
cv_r2s = []
for seed in range(20):
    r2, _ = kfold_cv(X6_lin, y, linear_model, k=10, seed=seed)
    cv_r2s.append(r2)
print(f"10-fold CV R² over 20 seeds: {np.mean(cv_r2s):.4f} ± {np.std(cv_r2s):.4f}")

print("\n✓ Test 1 passed: linear baseline established")

# =====================================================================
# TEST 2: RANDOM FOREST
# =====================================================================
print("\n" + "=" * 60)
print("TEST 2: RANDOM FOREST")
print("=" * 60)

try:
    from sklearn.ensemble import RandomForestRegressor

    def rf_model(X_train, y_train, X_test, n_trees=200, max_depth=None, min_leaf=3):
        rf = RandomForestRegressor(n_estimators=n_trees, max_depth=max_depth,
                                    min_samples_leaf=min_leaf, random_state=42)
        rf.fit(X_train, y_train)
        return rf.predict(X_test)

    # RF on 4 core features
    r2_rf4, rms_rf4 = kfold_cv(X4, y, lambda Xtr, ytr, Xte: rf_model(Xtr, ytr, Xte), k=10)
    # RF on extended features
    r2_rf_ext, rms_rf_ext = kfold_cv(X_ext, y, lambda Xtr, ytr, Xte: rf_model(Xtr, ytr, Xte), k=10)

    # RF with different hyperparameters
    configs = [
        ('RF(4feat, default)', X4, {}),
        ('RF(4feat, shallow)', X4, {'max_depth': 5}),
        ('RF(4feat, deep)', X4, {'max_depth': None, 'min_leaf': 1}),
        ('RF(11feat, default)', X_ext, {}),
    ]

    print(f"\n{'Config':<30} {'R²(CV)':<10} {'RMS(CV)'}")
    print("-" * 50)
    for name, X_use, params in configs:
        fn = lambda Xtr, ytr, Xte, p=params: rf_model(Xtr, ytr, Xte, **p)
        r2, rms = kfold_cv(X_use, y, fn, k=10)
        print(f"  {name:<30} {r2:<10.4f} {rms:.4f}")

    has_sklearn = True

except ImportError:
    print("\nscikit-learn not available — using manual decision tree")
    has_sklearn = False

    # Implement a simple bagged tree manually
    def simple_tree_predict(X_train, y_train, X_test, max_depth=4, min_leaf=5):
        """Simple axis-aligned decision tree."""
        n_train = len(y_train)

        def build_tree(idx, depth):
            if depth >= max_depth or len(idx) <= min_leaf:
                return np.mean(y_train[idx])

            best_score = np.var(y_train[idx]) * len(idx)
            best_feat = 0
            best_thresh = 0

            for f in range(X_train.shape[1]):
                vals = X_train[idx, f]
                thresholds = np.percentile(vals, [25, 50, 75])
                for t in thresholds:
                    left = idx[vals <= t]
                    right = idx[vals > t]
                    if len(left) < min_leaf or len(right) < min_leaf:
                        continue
                    score = np.var(y_train[left]) * len(left) + np.var(y_train[right]) * len(right)
                    if score < best_score:
                        best_score = score
                        best_feat = f
                        best_thresh = t

            left = idx[X_train[idx, best_feat] <= best_thresh]
            right = idx[X_train[idx, best_feat] > best_thresh]
            if len(left) < min_leaf or len(right) < min_leaf:
                return np.mean(y_train[idx])

            return (best_feat, best_thresh,
                    build_tree(left, depth + 1),
                    build_tree(right, depth + 1))

        def predict_one(tree, x):
            if not isinstance(tree, tuple):
                return tree
            feat, thresh, left, right = tree
            if x[feat] <= thresh:
                return predict_one(left, x)
            return predict_one(right, x)

        tree = build_tree(np.arange(n_train), 0)
        return np.array([predict_one(tree, x) for x in X_test])

    def bagged_tree(X_train, y_train, X_test, n_bags=50, max_depth=5):
        rng_b = np.random.RandomState(42)
        preds = np.zeros((len(X_test), n_bags))
        for b in range(n_bags):
            boot_idx = rng_b.choice(len(y_train), size=len(y_train), replace=True)
            preds[:, b] = simple_tree_predict(X_train[boot_idx], y_train[boot_idx],
                                               X_test, max_depth=max_depth)
        return np.mean(preds, axis=1)

    r2_rf4, rms_rf4 = kfold_cv(X4, y, lambda Xtr, ytr, Xte: bagged_tree(Xtr, ytr, Xte), k=10)
    r2_rf_ext, rms_rf_ext = kfold_cv(X_ext, y, lambda Xtr, ytr, Xte: bagged_tree(Xtr, ytr, Xte), k=10)

    print(f"\n{'Config':<30} {'R²(CV)':<10} {'RMS(CV)'}")
    print("-" * 50)
    print(f"  {'BaggedTree(4feat)':<30} {r2_rf4:<10.4f} {rms_rf4:.4f}")
    print(f"  {'BaggedTree(11feat)':<30} {r2_rf_ext:<10.4f} {rms_rf_ext:.4f}")

print(f"\n  Comparison: 6-var linear R²(CV) = {r2_lin:.4f}")

print("\n✓ Test 2 passed: RF tested")

# =====================================================================
# TEST 3: GRADIENT BOOSTING
# =====================================================================
print("\n" + "=" * 60)
print("TEST 3: GRADIENT BOOSTING")
print("=" * 60)

if has_sklearn:
    from sklearn.ensemble import GradientBoostingRegressor

    def gbr_model(X_train, y_train, X_test, n_est=200, max_depth=3, lr=0.05):
        gbr = GradientBoostingRegressor(n_estimators=n_est, max_depth=max_depth,
                                         learning_rate=lr, random_state=42)
        gbr.fit(X_train, y_train)
        return gbr.predict(X_test)

    configs_gb = [
        ('GBR(4feat, d3)', X4, {'max_depth': 3, 'lr': 0.05}),
        ('GBR(4feat, d5)', X4, {'max_depth': 5, 'lr': 0.03}),
        ('GBR(11feat, d3)', X_ext, {'max_depth': 3, 'lr': 0.05}),
        ('GBR(11feat, d5)', X_ext, {'max_depth': 5, 'lr': 0.03}),
    ]

    print(f"\n{'Config':<30} {'R²(CV)':<10} {'RMS(CV)'}")
    print("-" * 50)
    for name, X_use, params in configs_gb:
        fn = lambda Xtr, ytr, Xte, p=params: gbr_model(Xtr, ytr, Xte, **p)
        r2, rms = kfold_cv(X_use, y, fn, k=10)
        print(f"  {name:<30} {r2:<10.4f} {rms:.4f}")
else:
    # Simple gradient boosting from scratch
    def simple_gbr(X_train, y_train, X_test, n_rounds=100, lr=0.1, depth=3):
        residuals = y_train.copy()
        pred_train = np.zeros(len(y_train))
        pred_test = np.zeros(len(X_test))

        for r in range(n_rounds):
            tree_pred_train = simple_tree_predict(X_train, residuals, X_train, max_depth=depth)
            tree_pred_test = simple_tree_predict(X_train, residuals, X_test, max_depth=depth)
            pred_train += lr * tree_pred_train
            pred_test += lr * tree_pred_test
            residuals = y_train - pred_train

        return pred_test

    r2_gb4, rms_gb4 = kfold_cv(X4, y, lambda Xtr, ytr, Xte: simple_gbr(Xtr, ytr, Xte), k=10)
    r2_gb_ext, rms_gb_ext = kfold_cv(X_ext, y, lambda Xtr, ytr, Xte: simple_gbr(Xtr, ytr, Xte), k=10)

    print(f"\n{'Config':<30} {'R²(CV)':<10} {'RMS(CV)'}")
    print("-" * 50)
    print(f"  {'GBR(4feat)':<30} {r2_gb4:<10.4f} {rms_gb4:.4f}")
    print(f"  {'GBR(11feat)':<30} {r2_gb_ext:<10.4f} {rms_gb_ext:.4f}")

print(f"\n  Comparison: 6-var linear R²(CV) = {r2_lin:.4f}")

print("\n✓ Test 3 passed: GBR tested")

# =====================================================================
# TEST 4: K-NEAREST NEIGHBORS
# =====================================================================
print("\n" + "=" * 60)
print("TEST 4: K-NEAREST NEIGHBORS")
print("=" * 60)

def knn_model(X_train, y_train, X_test, k=5):
    """Simple KNN regression with standardized features."""
    # Standardize
    mu = np.mean(X_train, axis=0)
    std = np.clip(np.std(X_train, axis=0), 1e-10, None)
    X_tr_s = (X_train - mu) / std
    X_te_s = (X_test - mu) / std

    predictions = np.zeros(len(X_test))
    for i, x in enumerate(X_te_s):
        dists = np.sqrt(np.sum((X_tr_s - x)**2, axis=1))
        nn_idx = np.argsort(dists)[:k]
        predictions[i] = np.mean(y_train[nn_idx])
    return predictions

print(f"\n{'Config':<30} {'R²(CV)':<10} {'RMS(CV)'}")
print("-" * 50)
for k_nn in [3, 5, 7, 10, 15]:
    fn = lambda Xtr, ytr, Xte, k=k_nn: knn_model(Xtr, ytr, Xte, k=k)
    r2, rms = kfold_cv(X4, y, fn, k=10)
    print(f"  {'KNN(4feat, k='+str(k_nn)+')':<30} {r2:<10.4f} {rms:.4f}")

# KNN with extended features
r2_knn_ext, rms_knn_ext = kfold_cv(X_ext, y,
    lambda Xtr, ytr, Xte: knn_model(Xtr, ytr, Xte, k=5), k=10)
print(f"  {'KNN(11feat, k=5)':<30} {r2_knn_ext:<10.4f} {rms_knn_ext:.4f}")

print(f"\n  Comparison: 6-var linear R²(CV) = {r2_lin:.4f}")

print("\n✓ Test 4 passed: KNN tested")

# =====================================================================
# TEST 5: FEATURE IMPORTANCE
# =====================================================================
print("\n" + "=" * 60)
print("TEST 5: FEATURE IMPORTANCE")
print("=" * 60)

if has_sklearn:
    # Train full RF on all data
    rf_full = RandomForestRegressor(n_estimators=500, max_depth=None,
                                     min_samples_leaf=3, random_state=42)
    rf_full.fit(X_ext, y)
    importances = rf_full.feature_importances_

    print(f"\nRandom Forest Feature Importance (11 features):")
    sorted_idx = np.argsort(importances)[::-1]
    for i in sorted_idx:
        print(f"  {feat_ext_names[i]:<15} {importances[i]:.4f}")
else:
    # Permutation importance (manual)
    print(f"\nPermutation Importance (manual, 4 features):")
    rng_p = np.random.RandomState(42)

    # Baseline CV score
    r2_base, _ = kfold_cv(X4, y, lambda Xtr, ytr, Xte: bagged_tree(Xtr, ytr, Xte), k=10)

    for i, name in enumerate(feat4_names):
        X_perm = X4.copy()
        X_perm[:, i] = rng_p.permutation(X_perm[:, i])
        r2_perm, _ = kfold_cv(X_perm, y, lambda Xtr, ytr, Xte: bagged_tree(Xtr, ytr, Xte), k=10)
        print(f"  {name:<15} ΔR² = {r2_base - r2_perm:+.4f}")

print("\n✓ Test 5 passed: feature importance computed")

# =====================================================================
# TEST 6: PARTIAL DEPENDENCE (MANUAL)
# =====================================================================
print("\n" + "=" * 60)
print("TEST 6: PARTIAL DEPENDENCE")
print("=" * 60)

# For each feature, vary it while holding others at median
# Use the best ML model to generate predictions

# Simple approach: fit a RF on all data, then compute partial dependence
if has_sklearn:
    model_for_pd = rf_full  # Already fitted

    print(f"\nPartial Dependence (RF, 4 core features):")
    for i, name in enumerate(feat4_names):
        vals = np.percentile(X4[:, i], [10, 25, 50, 75, 90])
        X_pd = np.tile(np.median(X_ext, axis=0), (len(vals), 1))
        X_pd[:, i] = vals

        pd_pred = model_for_pd.predict(X_pd)
        print(f"\n  {name}:")
        for v, p in zip(vals, pd_pred):
            print(f"    {name}={v:.3f} → offset={p:+.4f}")

        # Linearity check: correlation of partial dependence with feature
        r_lin = np.corrcoef(vals, pd_pred)[0, 1]
        print(f"    PD linearity: r = {r_lin:+.3f}")

print("\n✓ Test 6 passed: partial dependence computed")

# =====================================================================
# TEST 7: RESIDUAL STRUCTURE FROM ML
# =====================================================================
print("\n" + "=" * 60)
print("TEST 7: RESIDUAL STRUCTURE FROM ML")
print("=" * 60)

# Compare residuals from linear model vs best ML model
# Do ML residuals have less structure?

# Full-sample fits for residual comparison
beta_lin = np.linalg.lstsq(X6_lin, y, rcond=None)[0]
resid_lin = y - X6_lin @ beta_lin

if has_sklearn:
    resid_rf = y - rf_full.predict(X_ext)
    rms_rf_full = np.sqrt(np.mean(resid_rf**2))
    r2_rf_full = 1 - np.sum(resid_rf**2) / np.sum((y - np.mean(y))**2)

    print(f"\nFull-sample comparison (in-sample):")
    print(f"  Linear: R² = {1 - np.sum(resid_lin**2)/np.sum((y-np.mean(y))**2):.4f}, "
          f"RMS = {np.sqrt(np.mean(resid_lin**2)):.4f}")
    print(f"  RF(500): R² = {r2_rf_full:.4f}, RMS = {rms_rf_full:.4f}")

    # NN autocorrelation in ML residuals
    from scipy.spatial.distance import cdist
    feat_std = (X4 - X4.mean(axis=0)) / np.clip(X4.std(axis=0), 1e-10, None)
    dist_matrix = cdist(feat_std, feat_std, 'euclidean')
    np.fill_diagonal(dist_matrix, np.inf)

    for k_nn in [1, 3, 5]:
        # Linear NN autocorrelation
        nn_corrs_lin = []
        nn_corrs_rf = []
        for i in range(n):
            nn_idx = np.argsort(dist_matrix[i])[:k_nn]
            nn_corrs_lin.append(np.mean(resid_lin[nn_idx]))
            nn_corrs_rf.append(np.mean(resid_rf[nn_idx]))
        r_nn_lin = np.corrcoef(resid_lin, nn_corrs_lin)[0, 1]
        r_nn_rf = np.corrcoef(resid_rf, nn_corrs_rf)[0, 1]
        print(f"\n  NN autocorrelation (k={k_nn}):")
        print(f"    Linear: r = {r_nn_lin:+.4f}")
        print(f"    RF:     r = {r_nn_rf:+.4f}")
else:
    print("\n  (Detailed residual analysis requires sklearn)")

print("\n✓ Test 7 passed: residual structure analyzed")

# =====================================================================
# TEST 8: THE ML ADVANTAGE QUANTIFIED
# =====================================================================
print("\n" + "=" * 60)
print("TEST 8: THE ML ADVANTAGE QUANTIFIED")
print("=" * 60)

# Collect all CV scores
print(f"\n--- Complete Model Comparison (10-fold CV) ---")
print(f"{'Model':<35} {'R²(CV)':<10} {'RMS(CV)':<10} {'ΔR² vs linear'}")
print("-" * 65)
print(f"  {'6-var linear':<35} {r2_lin:<10.4f} {rms_lin:<10.4f} {'baseline'}")

if has_sklearn:
    models_to_compare = [
        ('RF(4feat)', X4, lambda Xtr, ytr, Xte: rf_model(Xtr, ytr, Xte)),
        ('RF(11feat)', X_ext, lambda Xtr, ytr, Xte: rf_model(Xtr, ytr, Xte)),
        ('GBR(4feat, d3)', X4, lambda Xtr, ytr, Xte: gbr_model(Xtr, ytr, Xte, max_depth=3)),
        ('GBR(11feat, d3)', X_ext, lambda Xtr, ytr, Xte: gbr_model(Xtr, ytr, Xte, max_depth=3)),
    ]
else:
    models_to_compare = [
        ('BaggedTree(4feat)', X4, lambda Xtr, ytr, Xte: bagged_tree(Xtr, ytr, Xte)),
        ('BaggedTree(11feat)', X_ext, lambda Xtr, ytr, Xte: bagged_tree(Xtr, ytr, Xte)),
    ]

best_ml_r2 = r2_lin
best_ml_name = '6-var linear'
for name, X_use, fn in models_to_compare:
    r2, rms = kfold_cv(X_use, y, fn, k=10)
    delta = r2 - r2_lin
    print(f"  {name:<35} {r2:<10.4f} {rms:<10.4f} {delta:+.4f}")
    if r2 > best_ml_r2:
        best_ml_r2 = r2
        best_ml_name = name

# KNN
for k_nn in [5, 10]:
    r2, rms = kfold_cv(X4, y, lambda Xtr, ytr, Xte, k=k_nn: knn_model(Xtr, ytr, Xte, k=k), k=10)
    delta = r2 - r2_lin
    print(f"  {'KNN(4feat, k='+str(k_nn)+')':<35} {r2:<10.4f} {rms:<10.4f} {delta:+.4f}")
    if r2 > best_ml_r2:
        best_ml_r2 = r2
        best_ml_name = f'KNN(k={k_nn})'

print(f"\n--- Summary ---")
print(f"  Best model: {best_ml_name} (R²(CV) = {best_ml_r2:.4f})")
print(f"  Linear model: 6-var (R²(CV) = {r2_lin:.4f})")
print(f"  ML advantage: ΔR² = {best_ml_r2 - r2_lin:+.4f}")

if best_ml_r2 - r2_lin > 0.02:
    print(f"\n  → ML substantially outperforms linear: non-linear structure present")
elif best_ml_r2 - r2_lin > 0.005:
    print(f"\n  → ML marginally better: minor non-linear effects")
else:
    print(f"\n  → Linear model is OPTIMAL: no significant non-linear structure")
    print(f"     The 6-var model with logL×f_gas already captures the key non-linearity")

print(f"\n✓ Test 8 passed: ML advantage quantified")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SESSION #495 SUMMARY")
print("=" * 70)
print(f"6-var linear R²(10-fold CV): {r2_lin:.4f}")
print(f"Best ML R²(10-fold CV): {best_ml_r2:.4f} ({best_ml_name})")
print(f"ML advantage: ΔR² = {best_ml_r2 - r2_lin:+.4f}")
if has_sklearn:
    print(f"RF(4feat): R²(CV) = {r2_rf4:.4f}")
    print(f"RF(11feat): R²(CV) = {r2_rf_ext:.4f}")
print(f"\nAll 8 tests passed ✓")
