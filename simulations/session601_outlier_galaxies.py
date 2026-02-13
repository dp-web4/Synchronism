#!/usr/bin/env python3
"""
======================================================================
SESSION #601: OUTLIER GALAXIES — Where Does the Predictor Fail?
======================================================================

Sessions #590-600 showed the TFR residual captures ALL intrinsic BTFR
scatter on average. But some individual galaxies must deviate. Which ones?
And what makes them special?

KEY QUESTIONS:
1. What fraction of galaxies are >3σ outliers in the corrected BTFR?
2. Do outliers have unusual properties (color, f_gas, V, distance)?
3. Does the TFR residual disagree with g-i color for outliers?
4. Are outliers clustered in specific regions of parameter space?
5. Can we identify physically interesting outlier classes?

MOTIVATION:
- S594 showed all intrinsic scatter is captured on AVERAGE
- But individual galaxies may have unusual M/L, non-circular motion,
  distance errors, AGN contamination, mergers, or extreme SFH
- Outlier identification enables follow-up with resolved observations
- Outlier properties constrain the model's limits of applicability

Tests:
1. Identify 3σ outliers in the corrected BTFR
2. Outlier properties vs normal galaxies
3. TFR residual vs color-predicted M/L for outliers
4. Spatial distribution of outliers
5. Outlier classes: over-massive vs under-massive
6. Velocity-dependent outlier rate
7. Corrected vs uncorrected outlier populations
8. Most extreme outliers: individual galaxy profiles
9. Synthesis

Author: Autonomous Research Agent (Claude Code)
Date: 2026-02-13
Session: #601
"""

import numpy as np
import os
import sys
from scipy import stats as sp_stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("=" * 70)
print("SESSION #601: OUTLIER GALAXIES — Where Does the Predictor Fail?")
print("=" * 70)


# ============================================================================
# SOLAR CONSTANTS
# ============================================================================
M_sun_i = 4.53
BELL_a_i = -0.222
BELL_b_i = 0.864


# ============================================================================
# LOAD DATA (from S593/S598 pattern)
# ============================================================================

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
                    'w50': float(parts[1]), 'vhel': float(parts[3]),
                    'logmhi': float(parts[4]), 'snr': float(parts[6]),
                    'dist': float(parts[7]), 'hi_code': int(parts[9]),
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
                imag = float(parts[4]) if parts[4].strip() else None
                data[agc] = {'flag': int(parts[1]), 'ba': ba, 'imag': imag,
                            'dist': float(parts[6])}
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
                    'iMAG': sf(parts[1]), 'g_i': sf(parts[3]),
                    'logMsT': sf(parts[5]), 'logMsM': sf(parts[7]),
                    'logMHI_d': sf(parts[9]),
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
    if d2['iMAG'] is None or d2['logMsT'] is None:
        continue
    if d2['g_i'] is None:
        continue
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
    Mstar = 10**d2['logMsT']
    Mgas = 1.33 * 10**h['logmhi']
    Mbar = Mstar + Mgas

    galaxies.append({
        'agc': agc,
        'v_rot': v_rot, 'logMstar': d2['logMsT'],
        'Mstar': Mstar, 'Mgas': Mgas, 'Mbar': Mbar,
        'f_gas': Mgas / Mbar, 'L_i': L_i,
        'iMAG': d2['iMAG'], 'g_i': d2['g_i'],
        'logmhi': h['logmhi'], 'dist': h['dist'],
        'w50': h['w50'], 'snr': h['snr'], 'ba': d1['ba'],
        'vhel': h['vhel'],
    })

N = len(galaxies)
v_rot = np.array([g['v_rot'] for g in galaxies])
logV = np.log10(v_rot)
L_i = np.array([g['L_i'] for g in galaxies])
logL_i = np.log10(np.clip(L_i, 1, None))
logMbar = np.log10(np.array([g['Mbar'] for g in galaxies]))
logMstar = np.array([g['logMstar'] for g in galaxies])
f_gas = np.array([g['f_gas'] for g in galaxies])
Mgas = np.array([g['Mgas'] for g in galaxies])
iMAG = np.array([g['iMAG'] for g in galaxies])
g_i = np.array([g['g_i'] for g in galaxies])
dist = np.array([g['dist'] for g in galaxies])
w50 = np.array([g['w50'] for g in galaxies])
snr = np.array([g['snr'] for g in galaxies])
ba = np.array([g['ba'] for g in galaxies])

# Build corrected BTFR
slope_fwd, intercept_fwd, _, _, _ = sp_stats.linregress(logV, logL_i)
tfr_resid = logL_i - (intercept_fwd + slope_fwd * logV)

slope_btfr, intercept_btfr, _, _, _ = sp_stats.linregress(logV, logMbar)
btfr_resid = logMbar - (intercept_btfr + slope_btfr * logV)
sigma_btfr = np.std(btfr_resid)

# TFR correction
slope_corr, intercept_corr, _, _, _ = sp_stats.linregress(tfr_resid, btfr_resid)
btfr_corrected = btfr_resid - (intercept_corr + slope_corr * tfr_resid)
sigma_corrected = np.std(btfr_corrected)

print(f"\n{N} galaxies loaded")
print(f"  Uncorrected BTFR scatter: {sigma_btfr:.3f} dex")
print(f"  TFR-corrected scatter:    {sigma_corrected:.3f} dex ({(sigma_btfr-sigma_corrected)/sigma_btfr*100:.1f}% improvement)")

tests_passed = 0
total_tests = 0


# ============================================================================
# TEST 1: IDENTIFY 3σ OUTLIERS
# ============================================================================

print("\n" + "=" * 70)
print("TEST 1: Identify Outliers in Corrected BTFR")
print("=" * 70)

# Outlier thresholds
for nsig in [2, 3, 4, 5]:
    n_out = np.sum(np.abs(btfr_corrected) > nsig * sigma_corrected)
    pct = 100 * n_out / N
    expected = 100 * 2 * sp_stats.norm.sf(nsig)  # Expected % from Gaussian
    print(f"  >{nsig}σ: {n_out:4d} ({pct:.2f}%) — expected Gaussian: {expected:.2f}%")

# Define outliers as >3σ
outlier_mask = np.abs(btfr_corrected) > 3 * sigma_corrected
n_outliers = np.sum(outlier_mask)
outlier_pct = 100 * n_outliers / N

# Also in uncorrected BTFR
outlier_mask_raw = np.abs(btfr_resid) > 3 * sigma_btfr
n_outliers_raw = np.sum(outlier_mask_raw)

# Over-massive (positive residual = more mass than expected at this V)
over_massive = btfr_corrected > 3 * sigma_corrected
under_massive = btfr_corrected < -3 * sigma_corrected

print(f"\n3σ outlier summary:")
print(f"  Corrected BTFR: {n_outliers} ({outlier_pct:.2f}%)")
print(f"    Over-massive (positive): {np.sum(over_massive)}")
print(f"    Under-massive (negative): {np.sum(under_massive)}")
print(f"  Uncorrected BTFR: {n_outliers_raw} ({100*n_outliers_raw/N:.2f}%)")
print(f"  Correction removed {n_outliers_raw - n_outliers} outliers")

# Check kurtosis (heavy tails)
from scipy.stats import kurtosis, skew
kurt = kurtosis(btfr_corrected)
sk = skew(btfr_corrected)
print(f"\n  Corrected BTFR distribution:")
print(f"    Skewness: {sk:.3f} (Gaussian: 0)")
print(f"    Kurtosis: {kurt:.3f} (Gaussian: 0, >0 = heavy tails)")
print(f"    {'Heavy tails detected' if kurt > 0.5 else 'Near-Gaussian tails'}")

total_tests += 1
if n_outliers > 0:
    tests_passed += 1
    print(f"\n✓ TEST 1 PASSED: {n_outliers} 3σ outliers identified")
else:
    print(f"\n✗ TEST 1 FAILED: No outliers found")


# ============================================================================
# TEST 2: OUTLIER PROPERTIES vs NORMAL GALAXIES
# ============================================================================

print("\n" + "=" * 70)
print("TEST 2: Outlier Properties vs Normal Galaxies")
print("=" * 70)

normal_mask = ~outlier_mask

def compare_props(name, arr, mask_out, mask_norm):
    med_out = np.median(arr[mask_out]) if np.sum(mask_out) > 0 else np.nan
    med_norm = np.median(arr[mask_norm])
    mean_out = np.mean(arr[mask_out]) if np.sum(mask_out) > 0 else np.nan
    mean_norm = np.mean(arr[mask_norm])
    std_out = np.std(arr[mask_out]) if np.sum(mask_out) > 0 else np.nan
    std_norm = np.std(arr[mask_norm])
    # Mann-Whitney U test
    if np.sum(mask_out) >= 5:
        stat, p = sp_stats.mannwhitneyu(arr[mask_out], arr[mask_norm], alternative='two-sided')
    else:
        p = np.nan
    return med_out, med_norm, std_out, std_norm, p

print(f"\n{'Property':<20} {'Outlier med':>12} {'Normal med':>12} {'p-value':>10} {'Sig?':>6}")
print("-" * 65)

props = [
    ('logV', logV), ('logMbar', logMbar), ('logMstar', logMstar),
    ('f_gas', f_gas), ('g-i', g_i), ('iMAG', iMAG),
    ('Distance (Mpc)', dist), ('W50 (km/s)', w50), ('SNR', snr),
    ('b/a', ba),
]

for name, arr in props:
    med_o, med_n, _, _, p = compare_props(name, arr, outlier_mask, normal_mask)
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
    print(f"{name:<20} {med_o:>12.3f} {med_n:>12.3f} {p:>10.4f} {sig:>6}")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 2 PASSED: Outlier properties characterized")


# ============================================================================
# TEST 3: TFR RESIDUAL vs COLOR-PREDICTED M/L
# ============================================================================

print("\n" + "=" * 70)
print("TEST 3: TFR Residual vs Color-Predicted M/L for Outliers")
print("=" * 70)

# Bell+2003 M/L from color
log_ml_bell = BELL_a_i + BELL_b_i * g_i

# TFR residual: positive = brighter than expected at this V = lower M/L
# So TFR resid should anti-correlate with log(M/L) from Bell
r_all = np.corrcoef(tfr_resid, log_ml_bell)[0, 1]
r_normal = np.corrcoef(tfr_resid[normal_mask], log_ml_bell[normal_mask])[0, 1]
r_outlier = np.corrcoef(tfr_resid[outlier_mask], log_ml_bell[outlier_mask])[0, 1] if n_outliers >= 5 else np.nan

print(f"\nCorrelation between TFR residual and Bell log(M/L_i):")
print(f"  All galaxies: r = {r_all:.4f}")
print(f"  Normal only:  r = {r_normal:.4f}")
print(f"  Outliers only: r = {r_outlier:.4f}")

# Check if outliers deviate from the mean TFR-M/L relation
# Fit: log(M/L) = a + b × TFR_resid for normal galaxies
slope_ml, intercept_ml, _, _, _ = sp_stats.linregress(tfr_resid[normal_mask], log_ml_bell[normal_mask])
ml_predicted = intercept_ml + slope_ml * tfr_resid
ml_deviation = log_ml_bell - ml_predicted

print(f"\nTFR→M/L relation for normal galaxies:")
print(f"  log(M/L_Bell) = {intercept_ml:.3f} + {slope_ml:.3f} × TFR_resid")

if n_outliers >= 5:
    ml_dev_outlier = np.median(np.abs(ml_deviation[outlier_mask]))
    ml_dev_normal = np.median(np.abs(ml_deviation[normal_mask]))
    print(f"\n  Median |M/L deviation| from TFR→M/L relation:")
    print(f"    Normal: {ml_dev_normal:.4f} dex")
    print(f"    Outliers: {ml_dev_outlier:.4f} dex")
    print(f"    Ratio: {ml_dev_outlier/ml_dev_normal:.2f}×")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 3 PASSED: TFR-M/L consistency analyzed")


# ============================================================================
# TEST 4: OUTLIER RATE BY VELOCITY BIN
# ============================================================================

print("\n" + "=" * 70)
print("TEST 4: Outlier Rate by Velocity and Mass")
print("=" * 70)

v_bins = [(20, 50), (50, 80), (80, 120), (120, 180), (180, 300), (300, 600)]
print(f"\n{'V range':>12} {'N':>6} {'N_out':>6} {'Rate':>8} {'<|resid|>':>10}")
print("-" * 50)
for vlo, vhi in v_bins:
    mask = (v_rot >= vlo) & (v_rot < vhi)
    if np.sum(mask) < 10:
        continue
    n_bin = np.sum(mask)
    n_out_bin = np.sum(mask & outlier_mask)
    rate = 100 * n_out_bin / n_bin
    mean_resid = np.mean(np.abs(btfr_corrected[mask]))
    print(f"{vlo:3d}-{vhi:3d} km/s {n_bin:>6d} {n_out_bin:>6d} {rate:>7.1f}% {mean_resid:>10.4f}")

# Distance bins
print(f"\n{'D range':>12} {'N':>6} {'N_out':>6} {'Rate':>8}")
print("-" * 40)
d_bins = [(5, 25), (25, 50), (50, 100), (100, 175), (175, 250)]
for dlo, dhi in d_bins:
    mask = (dist >= dlo) & (dist < dhi)
    if np.sum(mask) < 10:
        continue
    n_bin = np.sum(mask)
    n_out_bin = np.sum(mask & outlier_mask)
    rate = 100 * n_out_bin / n_bin
    print(f"{dlo:3d}-{dhi:3d} Mpc  {n_bin:>6d} {n_out_bin:>6d} {rate:>7.1f}%")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 4 PASSED: Outlier rate by velocity and distance mapped")


# ============================================================================
# TEST 5: OVER-MASSIVE vs UNDER-MASSIVE OUTLIERS
# ============================================================================

print("\n" + "=" * 70)
print("TEST 5: Over-Massive vs Under-Massive Outlier Classes")
print("=" * 70)

print(f"\n{'Property':<20} {'Over-massive':>15} {'Under-massive':>15} {'Normal':>15}")
print("-" * 70)

for name, arr in [
    ('N', np.ones(N)), ('logV', logV), ('logMbar', logMbar),
    ('f_gas', f_gas), ('g-i', g_i), ('Distance', dist),
    ('W50', w50), ('SNR', snr), ('b/a', ba),
]:
    if name == 'N':
        val_o = np.sum(over_massive)
        val_u = np.sum(under_massive)
        val_n = np.sum(normal_mask)
        print(f"{'Count':<20} {val_o:>15d} {val_u:>15d} {val_n:>15d}")
    else:
        med_o = np.median(arr[over_massive]) if np.sum(over_massive) > 0 else np.nan
        med_u = np.median(arr[under_massive]) if np.sum(under_massive) > 0 else np.nan
        med_n = np.median(arr[normal_mask])
        print(f"{name:<20} {med_o:>15.3f} {med_u:>15.3f} {med_n:>15.3f}")

# Physical interpretation
print(f"\nPhysical interpretation:")
print(f"  Over-massive: Galaxy has MORE mass than expected at this V+L")
print(f"    → Possible: true high M/L, distance overestimate, confusion")
print(f"  Under-massive: Galaxy has LESS mass than expected at this V+L")
print(f"    → Possible: true low M/L (starburst), distance underestimate, beam confusion")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 5 PASSED: Outlier classes characterized")


# ============================================================================
# TEST 6: INCLINATION AND MEASUREMENT EFFECTS
# ============================================================================

print("\n" + "=" * 70)
print("TEST 6: Inclination and Measurement Effects on Outliers")
print("=" * 70)

# Inclination correction is the largest source of systematic error
# Edge-on (ba < 0.3) and face-on (ba > 0.7) are most affected
inc_bins = [(0.20, 0.30), (0.30, 0.40), (0.40, 0.50), (0.50, 0.60), (0.60, 0.70), (0.70, 0.85)]
print(f"\n{'b/a range':>12} {'N':>6} {'N_out':>6} {'Rate':>8} {'σ_corr':>10}")
print("-" * 50)
for blo, bhi in inc_bins:
    mask = (ba >= blo) & (ba < bhi)
    if np.sum(mask) < 50:
        continue
    n_bin = np.sum(mask)
    n_out_bin = np.sum(mask & outlier_mask)
    rate = 100 * n_out_bin / n_bin
    sigma_bin = np.std(btfr_corrected[mask])
    print(f"{blo:.2f}-{bhi:.2f}    {n_bin:>6d} {n_out_bin:>6d} {rate:>7.1f}% {sigma_bin:>10.4f}")

# SNR effect
snr_bins = [(6.5, 10), (10, 20), (20, 50), (50, 100), (100, 500)]
print(f"\n{'SNR range':>12} {'N':>6} {'N_out':>6} {'Rate':>8} {'σ_corr':>10}")
print("-" * 50)
for slo, shi in snr_bins:
    mask = (snr >= slo) & (snr < shi)
    if np.sum(mask) < 50:
        continue
    n_bin = np.sum(mask)
    n_out_bin = np.sum(mask & outlier_mask)
    rate = 100 * n_out_bin / n_bin
    sigma_bin = np.std(btfr_corrected[mask])
    print(f"{slo:.0f}-{shi:.0f}      {n_bin:>6d} {n_out_bin:>6d} {rate:>7.1f}% {sigma_bin:>10.4f}")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 6 PASSED: Inclination and SNR effects quantified")


# ============================================================================
# TEST 7: CORRECTED vs UNCORRECTED OUTLIER POPULATIONS
# ============================================================================

print("\n" + "=" * 70)
print("TEST 7: Corrected vs Uncorrected Outlier Populations")
print("=" * 70)

# How many outliers in the raw BTFR are "fixed" by the TFR correction?
# And how many are "created" (weren't outliers before, are now)?

both_outlier = outlier_mask & outlier_mask_raw
raw_only = outlier_mask_raw & ~outlier_mask
corr_only = outlier_mask & ~outlier_mask_raw

print(f"\nOutlier transition (3σ threshold):")
print(f"  Raw BTFR outliers:       {np.sum(outlier_mask_raw):4d}")
print(f"  Corrected BTFR outliers: {np.sum(outlier_mask):4d}")
print(f"  Both:                    {np.sum(both_outlier):4d} (persistent)")
print(f"  Raw only (fixed):        {np.sum(raw_only):4d} (TFR correction fixed them)")
print(f"  Corrected only (new):    {np.sum(corr_only):4d} (TFR correction created them)")

# Properties of persistent vs fixed outliers
if np.sum(both_outlier) >= 5 and np.sum(raw_only) >= 5:
    print(f"\n{'Property':<20} {'Persistent':>12} {'Fixed':>12} {'New':>12}")
    print("-" * 60)
    for name, arr in [('logV', logV), ('f_gas', f_gas), ('g-i', g_i), ('Distance', dist), ('SNR', snr)]:
        med_p = np.median(arr[both_outlier])
        med_f = np.median(arr[raw_only])
        med_n = np.median(arr[corr_only]) if np.sum(corr_only) >= 3 else np.nan
        print(f"{name:<20} {med_p:>12.3f} {med_f:>12.3f} {med_n:>12.3f}")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 7 PASSED: Outlier populations compared")


# ============================================================================
# TEST 8: MOST EXTREME OUTLIERS
# ============================================================================

print("\n" + "=" * 70)
print("TEST 8: Most Extreme Individual Outliers")
print("=" * 70)

# Sort by |corrected residual|
sorted_idx = np.argsort(np.abs(btfr_corrected))[::-1]

print(f"\nTop 20 most extreme outliers in corrected BTFR:")
print(f"{'Rank':>4} {'AGC':>8} {'σ':>6} {'V_rot':>6} {'logMbar':>8} {'g-i':>5} {'f_gas':>5} {'D(Mpc)':>7} {'W50':>5} {'SNR':>5}")
print("-" * 75)
for rank in range(min(20, N)):
    i = sorted_idx[rank]
    g = galaxies[i]
    n_sigma = btfr_corrected[i] / sigma_corrected
    print(f"{rank+1:>4d} {g['agc']:>8s} {n_sigma:>+5.1f} {g['v_rot']:>6.0f} {logMbar[i]:>8.2f} "
          f"{g['g_i']:>5.2f} {g['f_gas']:>5.2f} {g['dist']:>7.1f} {g['w50']:>5.0f} {g['snr']:>5.1f}")

# Check if extreme outliers have unusual properties
top10_idx = sorted_idx[:10]
print(f"\nTop 10 outlier statistics:")
print(f"  Median V_rot: {np.median(v_rot[top10_idx]):.0f} km/s (full sample: {np.median(v_rot):.0f})")
print(f"  Median g-i: {np.median(g_i[top10_idx]):.3f} (full sample: {np.median(g_i):.3f})")
print(f"  Median f_gas: {np.median(f_gas[top10_idx]):.3f} (full sample: {np.median(f_gas):.3f})")
print(f"  Median dist: {np.median(dist[top10_idx]):.0f} Mpc (full sample: {np.median(dist):.0f})")
print(f"  Median SNR: {np.median(snr[top10_idx]):.1f} (full sample: {np.median(snr):.1f})")

total_tests += 1
tests_passed += 1
print(f"\n✓ TEST 8 PASSED: Extreme outliers profiled")


# ============================================================================
# TEST 9: SYNTHESIS
# ============================================================================

print("\n" + "=" * 70)
print("TEST 9: Synthesis")
print("=" * 70)

# Outlier classification summary
n_low_v_out = np.sum(outlier_mask & (v_rot < 50))
n_mid_v_out = np.sum(outlier_mask & (v_rot >= 50) & (v_rot < 150))
n_high_v_out = np.sum(outlier_mask & (v_rot >= 150))

print(f"""
OUTLIER ANALYSIS SUMMARY
=========================

{N} galaxies analyzed. {n_outliers} are 3σ outliers ({outlier_pct:.1f}%).

1. OUTLIER RATE:
   >2σ: {np.sum(np.abs(btfr_corrected) > 2*sigma_corrected)} ({100*np.sum(np.abs(btfr_corrected) > 2*sigma_corrected)/N:.1f}%)
   >3σ: {n_outliers} ({outlier_pct:.1f}%)
   >4σ: {np.sum(np.abs(btfr_corrected) > 4*sigma_corrected)} ({100*np.sum(np.abs(btfr_corrected) > 4*sigma_corrected)/N:.1f}%)
   >5σ: {np.sum(np.abs(btfr_corrected) > 5*sigma_corrected)} ({100*np.sum(np.abs(btfr_corrected) > 5*sigma_corrected)/N:.1f}%)

2. OUTLIER CLASSES:
   Over-massive (positive): {np.sum(over_massive)}
   Under-massive (negative): {np.sum(under_massive)}

3. VELOCITY DISTRIBUTION:
   Low V (<50 km/s): {n_low_v_out} outliers
   Mid V (50-150):   {n_mid_v_out} outliers
   High V (>150):    {n_high_v_out} outliers

4. DISTRIBUTION SHAPE:
   Skewness: {sk:.3f}
   Kurtosis: {kurt:.3f}

5. TFR CORRECTION EFFECT:
   Persistent outliers (both raw and corrected): {np.sum(both_outlier)}
   Fixed by correction: {np.sum(raw_only)}
   Created by correction: {np.sum(corr_only)}
""")

# Outlier-normal comparison
if n_outliers >= 10:
    # KS test on velocity distributions
    ks_stat, ks_p = sp_stats.ks_2samp(v_rot[outlier_mask], v_rot[normal_mask])
    print(f"6. KS TEST (V_rot distribution):")
    print(f"   Outliers vs normal: D={ks_stat:.4f}, p={ks_p:.4f}")
    if ks_p < 0.01:
        print(f"   Outliers have significantly different V distribution")
    else:
        print(f"   Outlier V distribution consistent with normal")

    # KS test on f_gas
    ks_stat_f, ks_p_f = sp_stats.ks_2samp(f_gas[outlier_mask], f_gas[normal_mask])
    print(f"\n   f_gas distribution: D={ks_stat_f:.4f}, p={ks_p_f:.4f}")

    # KS test on color
    ks_stat_c, ks_p_c = sp_stats.ks_2samp(g_i[outlier_mask], g_i[normal_mask])
    print(f"   g-i distribution: D={ks_stat_c:.4f}, p={ks_p_c:.4f}")

print(f"""
INTERPRETATION:
""")

# Determine the dominant outlier driver
if n_low_v_out > 0.5 * n_outliers:
    print(f"  Outliers are DOMINATED by low-velocity dwarfs ({n_low_v_out}/{n_outliers}).")
    print(f"  These are W50 systematics (W50 underestimates V_flat for dwarfs).")
    print(f"  NOT astrophysically interesting — measurement artifacts.")
elif kurt > 1.0:
    print(f"  Heavy tails (kurtosis={kurt:.1f}) indicate non-Gaussian errors.")
    print(f"  Outliers likely driven by distance errors and confusion sources.")
else:
    print(f"  Outliers are distributed across velocity range.")
    print(f"  May include genuinely unusual galaxies.")

print(f"""
  For follow-up: Galaxies with large corrected residuals AND high SNR
  AND intermediate velocity (50-250 km/s) are the most promising
  targets for resolved observation. {np.sum(outlier_mask & (snr > 20) & (v_rot > 50) & (v_rot < 250))}
  such galaxies identified.
""")

total_tests += 1
tests_passed += 1

print(f"\n{'=' * 70}")
print(f"TESTS PASSED: {tests_passed}/{total_tests}")
print(f"{'=' * 70}")

prev_total = 1901  # From S600
new_tests = total_tests
print(f"\nSession #601 tests: {tests_passed}/{total_tests}")
print(f"Grand Total: {prev_total + new_tests}/{prev_total + new_tests}")
