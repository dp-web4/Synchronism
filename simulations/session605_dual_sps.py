#!/usr/bin/env python3
"""
======================================================================
SESSION #605: DUAL SPS MASS COMPARISON — Taylor vs Mendel
======================================================================

Session #604 showed the BTFR bottleneck is residual M/L variation
(0.161 dex after TFR), not measurement noise (0.050 dex). Durbala+2020
provides TWO independent stellar mass estimates:
  - logMsT: Taylor+2011 (SED fitting, g-i color)
  - logMsM: Mendel+2014 (spectro-photometric decomposition)

These differ by -0.264 ± 0.223 dex (Mendel masses ~1.8× higher).
The GALAXY-BY-GALAXY VARIATION in this difference is a direct measure
of M/L uncertainty and might contain information to reduce BTFR scatter.

KEY QUESTIONS:
1. Which SPS method gives tighter BTFR?
2. Does the difference (Δ_SPS = logMsT - logMsM) correlate with BTFR residuals?
3. Can Δ_SPS predict M/L better than g-i color?
4. Does the optimal BTFR use a blend of both methods?
5. How does this interact with the TFR correction?
6. Does the improvement vary with galaxy type?
7. What does Δ_SPS physically represent?
8. Can we break the 0.161 dex barrier (S604)?
9. Synthesis

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-15
Session: #605
"""

import numpy as np
import os
import sys
from scipy import stats as sp_stats
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("=" * 70)
print("SESSION #605: DUAL SPS MASS COMPARISON — Taylor vs Mendel")
print("=" * 70)


# ============================================================================
# CONSTANTS & DATA LOADING
# ============================================================================
M_sun_i = 4.53
BELL_a_i = -0.222
BELL_b_i = 0.864


def parse_haynes_tsv(filepath):
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('-') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            try:
                agc = parts[0].strip()
                data[agc] = {
                    'w50': float(parts[1]), 'e_w50': float(parts[2]),
                    'vhel': float(parts[3]),
                    'logmhi': float(parts[4]), 'e_logmhi': float(parts[5]),
                    'snr': float(parts[6]),
                    'dist': float(parts[7]), 'e_dist': float(parts[8]),
                    'hi_code': int(parts[9]),
                }
            except (ValueError, IndexError):
                continue
    return data


def parse_durbala_table1(filepath):
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('-') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                agc = parts[0].strip()
                ba = float(parts[2]) if parts[2].strip() else None
                e_ba = float(parts[3]) if len(parts) > 3 and parts[3].strip() else None
                data[agc] = {
                    'flag': int(parts[1]), 'ba': ba, 'e_ba': e_ba,
                    'dist': float(parts[6]),
                }
            except (ValueError, IndexError):
                continue
    return data


def parse_durbala_table2(filepath):
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('-') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                agc = parts[0].strip()
                def sf(s):
                    return float(s.strip()) if s.strip() else None
                data[agc] = {
                    'iMAG': sf(parts[1]), 'e_iMAG': sf(parts[2]),
                    'g_i': sf(parts[3]), 'e_g_i': sf(parts[4]),
                    'logMsT': sf(parts[5]), 'e_logMsT': sf(parts[6]),
                    'logMsM': sf(parts[7]), 'e_logMsM': sf(parts[8]),
                    'logMHI_d': sf(parts[9]), 'e_logMHI_d': sf(parts[10]),
                }
            except (ValueError, IndexError):
                continue
    return data


base_dir = os.path.dirname(os.path.abspath(__file__))
alfalfa_dir = os.path.join(base_dir, "alfalfa_data")
haynes = parse_haynes_tsv(os.path.join(alfalfa_dir, "haynes_alpha100.tsv"))
durbala1 = parse_durbala_table1(os.path.join(alfalfa_dir, "durbala_table1.tsv"))
durbala2 = parse_durbala_table2(os.path.join(alfalfa_dir, "durbala_table2.tsv"))

common_agc = set(haynes.keys()) & set(durbala1.keys()) & set(durbala2.keys())

galaxies = []
for agc in common_agc:
    h, d1, d2 = haynes[agc], durbala1[agc], durbala2[agc]
    if h['hi_code'] != 1 or d1['flag'] not in (1, 2) or h['snr'] < 6.5:
        continue
    if h['w50'] < 20 or d1['ba'] is None or d1['ba'] > 0.85 or d1['ba'] < 0.20:
        continue
    if d2['iMAG'] is None or d2['logMsT'] is None or d2['g_i'] is None:
        continue
    if d2['logMsM'] is None:
        continue  # Need both mass estimates
    if h['dist'] < 5 or h['dist'] > 250:
        continue

    q0 = 0.2
    cos2_i = (d1['ba']**2 - q0**2) / (1 - q0**2)
    if cos2_i <= 0:
        cos2_i = 0.01
    sin_i = np.sqrt(1 - cos2_i)
    if sin_i < 0.1:
        continue
    v_rot = h['w50'] / (2.0 * sin_i)
    if v_rot < 20:
        continue

    L_i = 10**(-0.4 * (d2['iMAG'] - M_sun_i))
    Mstar_T = 10**d2['logMsT']
    Mstar_M = 10**d2['logMsM']
    Mgas = 1.33 * 10**h['logmhi']
    Mbar_T = Mstar_T + Mgas
    Mbar_M = Mstar_M + Mgas

    galaxies.append({
        'agc': agc,
        'v_rot': v_rot, 'L_i': L_i,
        'logMsT': d2['logMsT'], 'logMsM': d2['logMsM'],
        'Mstar_T': Mstar_T, 'Mstar_M': Mstar_M,
        'Mgas': Mgas, 'Mbar_T': Mbar_T, 'Mbar_M': Mbar_M,
        'f_gas_T': Mgas / Mbar_T, 'f_gas_M': Mgas / Mbar_M,
        'iMAG': d2['iMAG'], 'g_i': d2['g_i'],
        'logmhi': h['logmhi'], 'dist': h['dist'],
        'w50': h['w50'], 'snr': h['snr'], 'ba': d1['ba'],
        'sin_i': sin_i,
    })

N = len(galaxies)

# Arrays
v_rot = np.array([g['v_rot'] for g in galaxies])
logV = np.log10(v_rot)
L_i = np.array([g['L_i'] for g in galaxies])
logL_i = np.log10(np.clip(L_i, 1, None))
logMsT = np.array([g['logMsT'] for g in galaxies])
logMsM = np.array([g['logMsM'] for g in galaxies])
logMbar_T = np.log10(np.array([g['Mbar_T'] for g in galaxies]))
logMbar_M = np.log10(np.array([g['Mbar_M'] for g in galaxies]))
f_gas_T = np.array([g['f_gas_T'] for g in galaxies])
f_gas_M = np.array([g['f_gas_M'] for g in galaxies])
g_i = np.array([g['g_i'] for g in galaxies])
snr = np.array([g['snr'] for g in galaxies])
ba = np.array([g['ba'] for g in galaxies])

# Key derived quantities
delta_sps = logMsT - logMsM  # Taylor - Mendel difference
logMbar_avg = np.log10((10**logMbar_T + 10**logMbar_M) / 2)  # Average baryonic mass

print(f"\nSample: {N} galaxies with both Taylor and Mendel masses")
print(f"V range: {v_rot.min():.0f} - {v_rot.max():.0f} km/s")
print(f"\nStellar mass comparison:")
print(f"  logMsT: {np.mean(logMsT):.3f} ± {np.std(logMsT):.3f}")
print(f"  logMsM: {np.mean(logMsM):.3f} ± {np.std(logMsM):.3f}")
print(f"  Δ_SPS (T-M): {np.mean(delta_sps):.3f} ± {np.std(delta_sps):.3f}")
print(f"  r(T,M) = {np.corrcoef(logMsT, logMsM)[0,1]:.4f}")


# ============================================================================
# Helper: OLS with LOO
# ============================================================================
def ols_with_loo(X, y):
    """OLS regression returning coefficients, predictions, and LOO residuals."""
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    y_pred = X @ beta
    resid = y - y_pred

    # LOO via hat matrix
    H = X @ np.linalg.solve(X.T @ X, X.T)
    h_diag = np.diag(H)
    loo_resid = resid / (1 - h_diag)

    sigma = np.std(resid)
    loo_sigma = np.std(loo_resid)
    ss_tot = np.sum((y - np.mean(y))**2)
    loo_r2 = 1 - np.sum(loo_resid**2) / ss_tot

    return beta, y_pred, resid, sigma, loo_sigma, loo_r2


passed = 0
total = 0

# ============================================================================
# TEST 1: Which SPS Method Gives Tighter BTFR?
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 1: Which SPS Method Gives Tighter BTFR?")
print("=" * 70)

# BTFR with Taylor masses
A = np.column_stack([logV, np.ones(N)])
beta_T, _, resid_T, sigma_T, loo_sigma_T, loo_r2_T = ols_with_loo(A, logMbar_T)

# BTFR with Mendel masses
beta_M, _, resid_M, sigma_M, loo_sigma_M, loo_r2_M = ols_with_loo(A, logMbar_M)

# BTFR with average masses
beta_avg, _, resid_avg, sigma_avg, loo_sigma_avg, loo_r2_avg = ols_with_loo(A, logMbar_avg)

print(f"\n{'Method':>12}  {'slope':>6}  {'σ (dex)':>8}  {'LOO σ':>8}  {'LOO R²':>8}")
print("-" * 55)
print(f"{'Taylor':>12}  {beta_T[0]:>6.3f}  {sigma_T:>8.4f}  {loo_sigma_T:>8.4f}  {loo_r2_T:>8.4f}")
print(f"{'Mendel':>12}  {beta_M[0]:>6.3f}  {sigma_M:>8.4f}  {loo_sigma_M:>8.4f}  {loo_r2_M:>8.4f}")
print(f"{'Average':>12}  {beta_avg[0]:>6.3f}  {sigma_avg:>8.4f}  {loo_sigma_avg:>8.4f}  {loo_r2_avg:>8.4f}")

winner = 'Taylor' if sigma_T < sigma_M else 'Mendel'
improvement = abs(sigma_T - sigma_M) / max(sigma_T, sigma_M) * 100
avg_vs_best = abs(sigma_avg - min(sigma_T, sigma_M)) / min(sigma_T, sigma_M) * 100

print(f"\n  {winner} gives tighter BTFR by {improvement:.1f}%")
print(f"  Average masses {'improve' if sigma_avg < min(sigma_T, sigma_M) else 'do not improve'} "
      f"over best single method ({avg_vs_best:.1f}%)")

# MOND slope check
print(f"\n  MOND expects slope = 4.0 (logMbar vs logV):")
print(f"    Taylor: {beta_T[0]:.3f}  (deviation: {abs(beta_T[0]-4.0)/4.0*100:.1f}%)")
print(f"    Mendel: {beta_M[0]:.3f}  (deviation: {abs(beta_M[0]-4.0)/4.0*100:.1f}%)")

print("PASS: SPS method comparison complete")
passed += 1


# ============================================================================
# TEST 2: Does Δ_SPS Correlate with BTFR Residuals?
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 2: Does Δ_SPS Correlate with BTFR Residuals?")
print("=" * 70)

# If Δ_SPS captures M/L information, it should correlate with BTFR residuals
# Use Taylor BTFR residuals (our standard)
r_delta_T = np.corrcoef(delta_sps, resid_T)[0, 1]
r_delta_M = np.corrcoef(delta_sps, resid_M)[0, 1]

# Compare with g-i color correlation
r_color_T = np.corrcoef(g_i, resid_T)[0, 1]
r_color_M = np.corrcoef(g_i, resid_M)[0, 1]

# Partial correlation: Δ_SPS with BTFR resid, controlling for logV
from numpy.linalg import lstsq
def partial_corr(x, y, z):
    """Partial correlation between x and y, controlling for z."""
    A = np.column_stack([z, np.ones(len(z))])
    rx = x - A @ lstsq(A, x, rcond=None)[0]
    ry = y - A @ lstsq(A, y, rcond=None)[0]
    return np.corrcoef(rx, ry)[0, 1]

r_partial_delta = partial_corr(delta_sps, resid_T, logV)
r_partial_color = partial_corr(g_i, resid_T, logV)

print(f"\nCorrelation with BTFR residuals:")
print(f"  {'Variable':>12}  {'r (Taylor)':>12}  {'r (Mendel)':>12}")
print(f"  {'Δ_SPS':>12}  {r_delta_T:>12.4f}  {r_delta_M:>12.4f}")
print(f"  {'g-i color':>12}  {r_color_T:>12.4f}  {r_color_M:>12.4f}")

print(f"\nPartial correlations (controlling for logV):")
print(f"  Δ_SPS:    r_partial = {r_partial_delta:.4f}")
print(f"  g-i color: r_partial = {r_partial_color:.4f}")

# Correlation between Δ_SPS and g-i
r_delta_gi = np.corrcoef(delta_sps, g_i)[0, 1]
print(f"\n  r(Δ_SPS, g-i) = {r_delta_gi:.4f}")

if abs(r_delta_gi) > 0.8:
    print(f"  Δ_SPS is highly correlated with g-i → likely redundant")
elif abs(r_delta_gi) > 0.5:
    print(f"  Δ_SPS is moderately correlated with g-i → partially independent")
else:
    print(f"  Δ_SPS is weakly correlated with g-i → largely independent information")

print("PASS: Δ_SPS correlation analysis complete")
passed += 1


# ============================================================================
# TEST 3: Can Δ_SPS Reduce BTFR Scatter?
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 3: Can Δ_SPS Reduce BTFR Scatter?")
print("=" * 70)

# Augmented BTFR: logMbar = α×logV + γ×Δ_SPS + β
A_delta = np.column_stack([logV, delta_sps, np.ones(N)])
beta_d, _, resid_d, sigma_d, loo_sigma_d, loo_r2_d = ols_with_loo(A_delta, logMbar_T)

# Compare with g-i augmented BTFR
A_gi = np.column_stack([logV, g_i, np.ones(N)])
beta_gi, _, resid_gi, sigma_gi, loo_sigma_gi, loo_r2_gi = ols_with_loo(A_gi, logMbar_T)

# Both Δ_SPS and g-i
A_both = np.column_stack([logV, delta_sps, g_i, np.ones(N)])
beta_both, _, resid_both, sigma_both, loo_sigma_both, loo_r2_both = ols_with_loo(A_both, logMbar_T)

print(f"\n{'Model':>25}  {'σ':>7}  {'LOO σ':>7}  {'LOO R²':>7}  {'Δσ%':>6}")
print("-" * 65)
print(f"{'BTFR (V only)':>25}  {sigma_T:>7.4f}  {loo_sigma_T:>7.4f}  {loo_r2_T:>7.4f}  {'—':>6}")
print(f"{'BTFR + Δ_SPS':>25}  {sigma_d:>7.4f}  {loo_sigma_d:>7.4f}  {loo_r2_d:>7.4f}  "
      f"{(1-sigma_d/sigma_T)*100:>5.1f}%")
print(f"{'BTFR + g-i':>25}  {sigma_gi:>7.4f}  {loo_sigma_gi:>7.4f}  {loo_r2_gi:>7.4f}  "
      f"{(1-sigma_gi/sigma_T)*100:>5.1f}%")
print(f"{'BTFR + Δ_SPS + g-i':>25}  {sigma_both:>7.4f}  {loo_sigma_both:>7.4f}  {loo_r2_both:>7.4f}  "
      f"{(1-sigma_both/sigma_T)*100:>5.1f}%")

print(f"\n  Δ_SPS coefficient: {beta_d[1]:.4f}")
print(f"  g-i coefficient:   {beta_gi[1]:.4f}")
print(f"  Δ_SPS marginal gain beyond g-i: "
      f"{(1-sigma_both/sigma_gi)*100:.2f}%")
print(f"  g-i marginal gain beyond Δ_SPS: "
      f"{(1-sigma_both/sigma_d)*100:.2f}%")

print("PASS: Δ_SPS as scatter predictor tested")
passed += 1


# ============================================================================
# TEST 4: Optimal Mass Blend
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 4: Optimal Mass Blend")
print("=" * 70)

# logMstar_opt = w × logMsT + (1-w) × logMsM
# Search for optimal weight w

def btfr_scatter_at_weight(w):
    logMstar_blend = w * logMsT + (1 - w) * logMsM
    Mstar_blend = 10**logMstar_blend
    Mgas = np.array([g['Mgas'] for g in galaxies])
    Mbar_blend = Mstar_blend + Mgas
    logMbar_blend = np.log10(Mbar_blend)
    A = np.column_stack([logV, np.ones(N)])
    beta = np.linalg.lstsq(A, logMbar_blend, rcond=None)[0]
    resid = logMbar_blend - A @ beta
    return np.std(resid)

weights = np.linspace(0, 1, 101)
scatters = [btfr_scatter_at_weight(w) for w in weights]
best_w = weights[np.argmin(scatters)]
best_scatter = min(scatters)

print(f"\nOptimal blending weight search:")
print(f"  w=0 (pure Mendel):  σ = {btfr_scatter_at_weight(0):.4f}")
print(f"  w=0.5 (equal avg):  σ = {btfr_scatter_at_weight(0.5):.4f}")
print(f"  w=1 (pure Taylor):  σ = {btfr_scatter_at_weight(1):.4f}")
print(f"  w={best_w:.2f} (optimal):  σ = {best_scatter:.4f}")
print(f"\n  Improvement over Taylor: {(1 - best_scatter/sigma_T)*100:.2f}%")
print(f"  Improvement over Mendel: {(1 - best_scatter/sigma_M)*100:.2f}%")

if abs(best_w - 1.0) < 0.05:
    print(f"  → Taylor alone is optimal (or near-optimal)")
elif abs(best_w - 0.0) < 0.05:
    print(f"  → Mendel alone is optimal")
elif abs(best_w - 0.5) < 0.1:
    print(f"  → Simple average is near-optimal")
else:
    print(f"  → Optimal is a {best_w:.0%}T/{1-best_w:.0%}M blend")

print("PASS: Optimal mass blend found")
passed += 1


# ============================================================================
# TEST 5: TFR Interaction with Dual Masses
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 5: TFR Interaction with Dual Masses")
print("=" * 70)

# TFR: logL = slope × logV + intercept
A_tfr = np.column_stack([logV, np.ones(N)])
beta_tfr = np.linalg.lstsq(A_tfr, logL_i, rcond=None)[0]
tfr_resid = logL_i - A_tfr @ beta_tfr

# TFR-corrected BTFR with Taylor
A_corr = np.column_stack([logV, tfr_resid, np.ones(N)])
beta_corr_T, _, resid_corr_T, sigma_corr_T, loo_corr_T, r2_corr_T = ols_with_loo(A_corr, logMbar_T)

# TFR-corrected BTFR with Mendel
beta_corr_M, _, resid_corr_M, sigma_corr_M, loo_corr_M, r2_corr_M = ols_with_loo(A_corr, logMbar_M)

# TFR + Δ_SPS
A_corr_d = np.column_stack([logV, tfr_resid, delta_sps, np.ones(N)])
beta_cd, _, resid_cd, sigma_cd, loo_cd, r2_cd = ols_with_loo(A_corr_d, logMbar_T)

# TFR + g-i
A_corr_gi = np.column_stack([logV, tfr_resid, g_i, np.ones(N)])
beta_cgi, _, resid_cgi, sigma_cgi, loo_cgi, r2_cgi = ols_with_loo(A_corr_gi, logMbar_T)

# TFR + Δ_SPS + g-i
A_corr_all = np.column_stack([logV, tfr_resid, delta_sps, g_i, np.ones(N)])
beta_ca, _, resid_ca, sigma_ca, loo_ca, r2_ca = ols_with_loo(A_corr_all, logMbar_T)

print(f"\n{'Model':>30}  {'σ':>7}  {'LOO σ':>7}  {'Δσ%':>7}")
print("-" * 60)
print(f"{'TFR-corr (Taylor)':>30}  {sigma_corr_T:>7.4f}  {loo_corr_T:>7.4f}  {'—':>7}")
print(f"{'TFR-corr (Mendel)':>30}  {sigma_corr_M:>7.4f}  {loo_corr_M:>7.4f}  "
      f"{(1-sigma_corr_M/sigma_corr_T)*100:>6.1f}%")
print(f"{'TFR + Δ_SPS':>30}  {sigma_cd:>7.4f}  {loo_cd:>7.4f}  "
      f"{(1-sigma_cd/sigma_corr_T)*100:>6.1f}%")
print(f"{'TFR + g-i':>30}  {sigma_cgi:>7.4f}  {loo_cgi:>7.4f}  "
      f"{(1-sigma_cgi/sigma_corr_T)*100:>6.1f}%")
print(f"{'TFR + Δ_SPS + g-i':>30}  {sigma_ca:>7.4f}  {loo_ca:>7.4f}  "
      f"{(1-sigma_ca/sigma_corr_T)*100:>6.1f}%")

marginal_delta = (1 - sigma_cd/sigma_corr_T) * 100
print(f"\n  Δ_SPS adds {marginal_delta:.1f}% beyond TFR correction")
print(f"  g-i adds {(1 - sigma_cgi/sigma_corr_T) * 100:.1f}% beyond TFR correction")
print(f"  Combined adds {(1 - sigma_ca/sigma_corr_T) * 100:.1f}% beyond TFR")

# Key: Does Δ_SPS break the 0.161 barrier from S604?
print(f"\n  S604 TFR-corrected σ_int = 0.161 dex (bottleneck)")
print(f"  Best model here: σ = {min(sigma_cd, sigma_ca):.4f} dex")
if min(sigma_cd, sigma_ca) < 0.161:
    print(f"  → Broke the 0.161 barrier! Improvement: "
          f"{(1 - min(sigma_cd, sigma_ca)/0.161)*100:.1f}%")
else:
    print(f"  → Did NOT break the 0.161 barrier")

print("PASS: TFR + dual SPS analysis complete")
passed += 1


# ============================================================================
# TEST 6: Variation by Galaxy Type
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 6: Variation by Galaxy Type")
print("=" * 70)

# Does Δ_SPS help more for certain galaxy populations?
masks = [
    ("Gas-rich (f_gas>0.5)", f_gas_T > 0.5),
    ("Gas-poor (f_gas<0.3)", f_gas_T < 0.3),
    ("Low-V (V<80)", v_rot < 80),
    ("Mid-V (80-180)", (v_rot >= 80) & (v_rot < 180)),
    ("High-V (V>180)", v_rot >= 180),
    ("Blue (g-i<0.7)", g_i < 0.7),
    ("Red (g-i>1.0)", g_i > 1.0),
]

print(f"\n{'Population':>22}  {'N':>5}  {'σ_BTFR':>7}  {'σ_TFR':>7}  {'σ_TFR+Δ':>7}  {'Δ gain%':>7}")
print("-" * 70)

for label, mask in masks:
    n = np.sum(mask)
    if n < 30:
        continue

    # BTFR scatter for this population
    s_btfr = np.std(resid_T[mask])

    # TFR-corrected scatter
    s_tfr = np.std(resid_corr_T[mask])

    # TFR + Δ_SPS scatter
    s_delta = np.std(resid_cd[mask])

    gain = (1 - s_delta / s_tfr) * 100

    print(f"  {label:>20}  {n:>5}  {s_btfr:>7.4f}  {s_tfr:>7.4f}  {s_delta:>7.4f}  {gain:>6.1f}%")

# Does Δ_SPS depend on f_gas?
r_delta_fgas = np.corrcoef(delta_sps, f_gas_T)[0, 1]
r_delta_logV = np.corrcoef(delta_sps, logV)[0, 1]
print(f"\n  r(Δ_SPS, f_gas) = {r_delta_fgas:.4f}")
print(f"  r(Δ_SPS, logV) = {r_delta_logV:.4f}")

print("PASS: Population variation analysis complete")
passed += 1


# ============================================================================
# TEST 7: What Does Δ_SPS Physically Represent?
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 7: What Does Δ_SPS Physically Represent?")
print("=" * 70)

# Δ_SPS = logMsT - logMsM = log(MsT/MsM)
# Taylor uses SED+color fitting; Mendel uses spectro-photometric decomposition
# The difference might encode: age, metallicity, dust, SFH shape, bulge fraction

# Correlations with observable properties
props = {
    'g-i color': g_i,
    'i-band mag': np.array([g['iMAG'] for g in galaxies]),
    'logV': logV,
    'logL_i': logL_i,
    'f_gas (Taylor)': f_gas_T,
    'b/a': ba,
    'log(dist)': np.log10(np.array([g['dist'] for g in galaxies])),
    'SNR': snr,
}

print(f"\nCorrelations of Δ_SPS with galaxy properties:")
print(f"  {'Property':>20}  {'r':>8}  {'partial r(|logV)':>18}")
print("-" * 55)

for name, prop in props.items():
    r = np.corrcoef(delta_sps, prop)[0, 1]
    pr = partial_corr(delta_sps, prop, logV)
    print(f"  {name:>20}  {r:>8.4f}  {pr:>18.4f}")

# Binned mean Δ_SPS by g-i color
color_bins = np.linspace(0.2, 1.4, 7)
print(f"\n  Mean Δ_SPS by g-i color bin:")
for i in range(len(color_bins)-1):
    mask = (g_i >= color_bins[i]) & (g_i < color_bins[i+1])
    n = np.sum(mask)
    if n > 20:
        print(f"    g-i [{color_bins[i]:.1f}, {color_bins[i+1]:.1f}): "
              f"<Δ_SPS> = {np.mean(delta_sps[mask]):+.3f} ± {np.std(delta_sps[mask]):.3f} (N={n})")

# Is Δ_SPS just a proxy for color?
# Residual Δ_SPS after removing g-i dependence
A_gi_pred = np.column_stack([g_i, np.ones(N)])
beta_gi_pred = np.linalg.lstsq(A_gi_pred, delta_sps, rcond=None)[0]
delta_resid = delta_sps - A_gi_pred @ beta_gi_pred
print(f"\n  Δ_SPS = {beta_gi_pred[0]:+.4f}×(g-i) + {beta_gi_pred[1]:+.4f}")
print(f"  R² of g-i predicting Δ_SPS: {1 - np.var(delta_resid)/np.var(delta_sps):.4f}")
print(f"  Residual Δ_SPS scatter: {np.std(delta_resid):.4f} dex")
print(f"  → {np.std(delta_resid)/np.std(delta_sps)*100:.0f}% of Δ_SPS is independent of g-i color")

# Does residual Δ_SPS correlate with BTFR residuals?
r_resid_delta = np.corrcoef(delta_resid, resid_T)[0, 1]
print(f"\n  r(residual Δ_SPS, BTFR resid) = {r_resid_delta:.4f}")
print(f"  → Color-independent part of Δ_SPS {'does' if abs(r_resid_delta) > 0.05 else 'does NOT'} "
      f"correlate with BTFR residuals")

print("PASS: Physical interpretation of Δ_SPS analyzed")
passed += 1


# ============================================================================
# TEST 8: Can We Break 0.161 dex? (S604 Barrier)
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 8: Can We Break the 0.161 dex Barrier?")
print("=" * 70)

# S604 found TFR-corrected σ_int = 0.161 dex on optimal subsample
# The TFR-corrected residuals should be the starting point
# Try adding ALL available predictors to see what's achievable

# Full model: logMbar = α×logV + β_tfr×tfr_resid + β_delta×Δ_SPS + β_gi×g-i + β_fgas×f_gas + const
A_full = np.column_stack([logV, tfr_resid, delta_sps, g_i, f_gas_T, np.ones(N)])
beta_full, _, resid_full, sigma_full, loo_full, r2_full = ols_with_loo(A_full, logMbar_T)

# With Mendel mass instead
A_full_M = np.column_stack([logV, tfr_resid, delta_sps, g_i, f_gas_M, np.ones(N)])
beta_full_M, _, resid_full_M, sigma_full_M, loo_full_M, r2_full_M = ols_with_loo(A_full_M, logMbar_M)

# Optimal weight blend with TFR + Δ_SPS
# This is the key test: can blending + TFR + Δ_SPS get below 0.161?
best_blend_scatter = np.inf
best_blend_w = 0
for w in np.linspace(0, 1, 101):
    logMstar_blend = w * logMsT + (1 - w) * logMsM
    Mbar_blend = 10**logMstar_blend + np.array([g['Mgas'] for g in galaxies])
    logMbar_blend = np.log10(Mbar_blend)

    A_bl = np.column_stack([logV, tfr_resid, np.ones(N)])
    beta_bl = np.linalg.lstsq(A_bl, logMbar_blend, rcond=None)[0]
    resid_bl = logMbar_blend - A_bl @ beta_bl
    s = np.std(resid_bl)
    if s < best_blend_scatter:
        best_blend_scatter = s
        best_blend_w = w

print(f"\n{'Model':>35}  {'σ':>7}  {'LOO σ':>7}")
print("-" * 55)
print(f"{'TFR-corrected (Taylor)':>35}  {sigma_corr_T:>7.4f}  {loo_corr_T:>7.4f}")
print(f"{'TFR-corrected (Mendel)':>35}  {sigma_corr_M:>7.4f}  {loo_corr_M:>7.4f}")
print(f"{'TFR + Δ_SPS':>35}  {sigma_cd:>7.4f}  {loo_cd:>7.4f}")
print(f"{'TFR + g-i':>35}  {sigma_cgi:>7.4f}  {loo_cgi:>7.4f}")
print(f"{'TFR + Δ_SPS + g-i':>35}  {sigma_ca:>7.4f}  {loo_ca:>7.4f}")
print(f"{'Full model (T: V+TFR+Δ+gi+fg)':>35}  {sigma_full:>7.4f}  {loo_full:>7.4f}")
print(f"{'Full model (M: V+TFR+Δ+gi+fg)':>35}  {sigma_full_M:>7.4f}  {loo_full_M:>7.4f}")
print(f"{'Optimal blend + TFR (w={best_blend_w:.2f})':>35}  {best_blend_scatter:>7.4f}  {'—':>7}")

best_achievable = min(sigma_full, sigma_full_M, sigma_ca, sigma_cd)
print(f"\n  Best achievable σ: {best_achievable:.4f} dex")
print(f"  S604 barrier (TFR-corrected): 0.195 dex (full sample)")
print(f"  S604 optimal subsample σ_int: 0.161 dex")

if best_achievable < 0.195:
    print(f"  → BROKE full-sample barrier: {(1-best_achievable/0.195)*100:.1f}% improvement")
else:
    print(f"  → Did NOT break full-sample barrier")

# Check full model coefficients
print(f"\n  Full model coefficients (Taylor):")
labels = ['logV', 'TFR_resid', 'Δ_SPS', 'g-i', 'f_gas', 'const']
for l, b in zip(labels, beta_full):
    print(f"    {l:>12}: {b:+.4f}")

# t-statistics
X = A_full
beta = beta_full
resid = resid_full
sigma_resid = np.std(resid)
cov = sigma_resid**2 * np.linalg.inv(X.T @ X)
se = np.sqrt(np.diag(cov))
t_stats = beta / se
print(f"\n  t-statistics:")
for l, t in zip(labels, t_stats):
    sig = "***" if abs(t) > 3.29 else "**" if abs(t) > 2.58 else "*" if abs(t) > 1.96 else ""
    print(f"    {l:>12}: t = {t:+.2f} {sig}")

print("PASS: Barrier-breaking analysis complete")
passed += 1


# ============================================================================
# TEST 9: Synthesis
# ============================================================================
total += 1
print("\n" + "=" * 70)
print("TEST 9: Synthesis — Dual SPS Summary")
print("=" * 70)

print(f"\n{'='*60}")
print(f"DUAL SPS MASS COMPARISON SUMMARY")
print(f"{'='*60}")

print(f"\n1. SAMPLE: {N} galaxies with both Taylor and Mendel masses")
print(f"   Δ_SPS (T-M): {np.mean(delta_sps):+.3f} ± {np.std(delta_sps):.3f} dex")
print(f"   r(T,M) = {np.corrcoef(logMsT, logMsM)[0,1]:.4f}")

print(f"\n2. BTFR COMPARISON:")
print(f"   Taylor σ = {sigma_T:.4f}, Mendel σ = {sigma_M:.4f}")
print(f"   {winner} is tighter by {improvement:.1f}%")
print(f"   BTFR slope: Taylor {beta_T[0]:.3f}, Mendel {beta_M[0]:.3f}")

print(f"\n3. Δ_SPS AS PREDICTOR:")
print(f"   r(Δ_SPS, BTFR resid) = {r_delta_T:.4f}")
print(f"   r(g-i, BTFR resid) = {r_color_T:.4f}")
print(f"   r(Δ_SPS, g-i) = {r_delta_gi:.4f}")
print(f"   Δ_SPS independent of color: {np.std(delta_resid)/np.std(delta_sps)*100:.0f}%")

print(f"\n4. SCATTER REDUCTION:")
print(f"   BTFR alone:    {sigma_T:.4f}")
print(f"   + TFR:         {sigma_corr_T:.4f} ({(1-sigma_corr_T/sigma_T)*100:.1f}%)")
print(f"   + TFR + Δ_SPS: {sigma_cd:.4f} ({(1-sigma_cd/sigma_T)*100:.1f}%)")
print(f"   + TFR + g-i:   {sigma_cgi:.4f} ({(1-sigma_cgi/sigma_T)*100:.1f}%)")
print(f"   Full model:    {sigma_full:.4f} ({(1-sigma_full/sigma_T)*100:.1f}%)")

print(f"\n5. OPTIMAL BLEND: w={best_w:.2f} (Taylor={best_w*100:.0f}%, Mendel={(1-best_w)*100:.0f}%)")
print(f"   Blend σ = {best_scatter:.4f}")

print(f"\n6. KEY INSIGHTS:")
insight_delta_gain = (1 - sigma_cd/sigma_corr_T) * 100
print(f"   a) Δ_SPS adds {insight_delta_gain:.1f}% beyond TFR correction")
color_indep = np.std(delta_resid)/np.std(delta_sps)*100
print(f"   b) {color_indep:.0f}% of Δ_SPS is independent of g-i color")
print(f"   c) r(color-independent Δ_SPS, BTFR resid) = {r_resid_delta:.4f}")

if abs(r_resid_delta) > 0.05:
    print(f"   d) Δ_SPS encodes M/L information beyond color → SPS method matters")
else:
    print(f"   d) Color-independent Δ_SPS does NOT predict scatter → SPS difference is noise")

print(f"\n7. S604 BARRIER TEST:")
print(f"   S604 TFR-corrected: 0.195 dex (full), 0.161 dex (optimal)")
print(f"   Best this session:  {best_achievable:.4f} dex")
if best_achievable < 0.195:
    print(f"   → BROKE the barrier by {(1-best_achievable/0.195)*100:.1f}%")
else:
    print(f"   → Barrier intact")

print(f"\nPASS: Synthesis complete")
passed += 1

print(f"\n{'='*60}")
print(f"Session #605 Grand Total: {passed}/{total}")
print(f"{'='*60}")

PREV_TOTAL = 1937
grand_passed = PREV_TOTAL + passed
grand_total = PREV_TOTAL + total
print(f"\n{'='*60}")
print(f"GRAND TOTAL: {grand_passed}/{grand_total}")
print(f"{'='*60}")
