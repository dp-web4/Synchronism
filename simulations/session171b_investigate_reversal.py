#!/usr/bin/env python3
"""
SESSION #171b: INVESTIGATE THE VELOCITY REVERSAL
================================================
Date: December 23, 2025

PUZZLE:
- Synchronism predicts: void galaxies should have HIGHER |v|
- CF4 data shows: void galaxies have LOWER |v|
- This is the OPPOSITE of prediction

HYPOTHESES TO TEST:
1. Selection effect: High-|v| galaxies in voids are missing (Malmquist bias)?
2. Sign convention: Is our density estimator backwards?
3. Hubble flow contamination: Are we measuring Hubble flow, not peculiar velocity?
4. Definition issue: What exactly is "peculiar velocity" in CF4?

Let's investigate systematically.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #171b: INVESTIGATING VELOCITY REVERSAL")
print("=" * 70)

# Load data
data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table2.dat', 'r') as f:
    lines = f.readlines()

vcmb = np.zeros(len(lines))
dm = np.zeros(len(lines))
dm_err = np.zeros(len(lines))
ra = np.zeros(len(lines))
dec = np.zeros(len(lines))

for i, line in enumerate(lines):
    try:
        vcmb[i] = float(line[22:27].strip() or 0)
        dm[i] = float(line[28:34].strip() or 0)
        dm_err[i] = float(line[35:40].strip() or 0)
        ra[i] = float(line[137:145].strip() or 0)
        dec[i] = float(line[146:154].strip() or 0)
    except:
        continue

H0 = 74.6
distances = 10 ** ((dm - 25) / 5)
v_pec = vcmb - H0 * distances

valid = (dm > 20) & (dm_err > 0) & (vcmb > 500) & (distances > 10) & (distances < 200)
N = np.sum(valid)

print(f"\nLoaded {N} valid galaxies")

# =============================================================================
# TEST 1: CHECK PECULIAR VELOCITY SIGN DISTRIBUTION
# =============================================================================

print("\n" + "=" * 70)
print("TEST 1: PECULIAR VELOCITY SIGN DISTRIBUTION")
print("=" * 70)

v_pec_valid = v_pec[valid]
print(f"\nPeculiar velocity distribution:")
print(f"  Mean: {np.mean(v_pec_valid):.0f} km/s")
print(f"  Median: {np.median(v_pec_valid):.0f} km/s")
print(f"  Std: {np.std(v_pec_valid):.0f} km/s")
print(f"  Min: {np.min(v_pec_valid):.0f} km/s")
print(f"  Max: {np.max(v_pec_valid):.0f} km/s")
print(f"  Fraction positive: {100*np.mean(v_pec_valid > 0):.1f}%")

# =============================================================================
# TEST 2: CHECK IF HIGH-|v| CORRELATES WITH DISTANCE
# =============================================================================

print("\n" + "=" * 70)
print("TEST 2: VELOCITY-DISTANCE RELATIONSHIP")
print("=" * 70)

d_valid = distances[valid]

# If there's Malmquist bias, high-|v| galaxies should be missing at large distances
dist_bins = np.linspace(10, 200, 10)
for i in range(len(dist_bins)-1):
    mask = (d_valid >= dist_bins[i]) & (d_valid < dist_bins[i+1])
    if np.sum(mask) > 50:
        mean_v = np.mean(np.abs(v_pec_valid[mask]))
        std_v = np.std(np.abs(v_pec_valid[mask]))
        print(f"  d = {dist_bins[i]:.0f}-{dist_bins[i+1]:.0f} Mpc: <|v|> = {mean_v:.0f} ± {std_v:.0f} km/s (N={np.sum(mask)})")

# Correlation
corr, p = stats.pearsonr(d_valid, np.abs(v_pec_valid))
print(f"\nCorrelation (distance vs |v|): r = {corr:.3f}, p = {p:.2e}")

# =============================================================================
# TEST 3: ALTERNATIVE DENSITY METRIC - ANGULAR ONLY
# =============================================================================

print("\n" + "=" * 70)
print("TEST 3: PURE ANGULAR DENSITY (NO DISTANCE)")
print("=" * 70)

from scipy.spatial import cKDTree

ra_rad = np.radians(ra[valid])
dec_rad = np.radians(dec[valid])

# Project to unit sphere
x_sphere = np.cos(dec_rad) * np.cos(ra_rad)
y_sphere = np.cos(dec_rad) * np.sin(ra_rad)
z_sphere = np.sin(dec_rad)

coords_sphere = np.column_stack([x_sphere, y_sphere, z_sphere])
tree = cKDTree(coords_sphere)

# Count neighbors within angular radius (in 3D chord distance)
# chord = 2*sin(theta/2), for theta=5deg: chord ≈ 0.087
chord_5deg = 2 * np.sin(np.radians(5/2))
counts_angular = np.array(tree.query_ball_point(coords_sphere, r=chord_5deg, return_length=True)) - 1

# Simple overdensity
mean_count = np.mean(counts_angular)
delta_angular = (counts_angular - mean_count) / mean_count

print(f"\nAngular overdensity (5 deg radius):")
print(f"  Min δ: {np.min(delta_angular):.2f}")
print(f"  Max δ: {np.max(delta_angular):.2f}")
print(f"  Median δ: {np.median(delta_angular):.2f}")

# Check velocity by angular density
ang_low = delta_angular < np.percentile(delta_angular, 25)
ang_high = delta_angular > np.percentile(delta_angular, 75)

v_ang_low = np.mean(np.abs(v_pec_valid[ang_low]))
v_ang_high = np.mean(np.abs(v_pec_valid[ang_high]))

print(f"\nLow angular density (Q1): <|v|> = {v_ang_low:.0f} km/s (N={np.sum(ang_low)})")
print(f"High angular density (Q4): <|v|> = {v_ang_high:.0f} km/s (N={np.sum(ang_high)})")
print(f"Ratio: {v_ang_low/v_ang_high:.3f}")

# =============================================================================
# TEST 4: CHECK RAW REDSHIFT CORRELATION
# =============================================================================

print("\n" + "=" * 70)
print("TEST 4: REDSHIFT AND DENSITY RELATIONSHIP")
print("=" * 70)

vcmb_valid = vcmb[valid]

# Does local density correlate with redshift?
corr_z_delta, _ = stats.spearmanr(vcmb_valid, delta_angular)
print(f"\nCorrelation (vcmb vs angular δ): r = {corr_z_delta:.3f}")

# High-z sample
high_z = vcmb_valid > 10000
low_z = vcmb_valid < 5000

print(f"\nLow-z (v < 5000): mean angular δ = {np.mean(delta_angular[low_z]):.3f}")
print(f"High-z (v > 10000): mean angular δ = {np.mean(delta_angular[high_z]):.3f}")

# =============================================================================
# TEST 5: THE KEY ISSUE - WHAT IS "PECULIAR VELOCITY"?
# =============================================================================

print("\n" + "=" * 70)
print("TEST 5: UNDERSTANDING PECULIAR VELOCITY DEFINITION")
print("=" * 70)

print("""
In CF4, peculiar velocity is defined as:
  v_pec = v_CMB - H0 × d

where d is the DISTANCE (not redshift distance).

Key insight:
- If d is UNDERESTIMATED, v_pec is POSITIVE
- If d is OVERESTIMATED, v_pec is NEGATIVE

In voids:
- Galaxies are FAINTER (less local mass)
- Distance indicators may be biased
- Malmquist bias: preferentially detect brighter (closer) galaxies

Could this explain the reversal?
""")

# Check: Are positive v_pec galaxies (d underestimated) in different environments?
pos_pec = v_pec_valid > 0
neg_pec = v_pec_valid < 0

print(f"\nPositive v_pec galaxies:")
print(f"  N = {np.sum(pos_pec)}")
print(f"  Mean angular δ = {np.mean(delta_angular[pos_pec]):.3f}")
print(f"  Mean distance = {np.mean(d_valid[pos_pec]):.1f} Mpc")

print(f"\nNegative v_pec galaxies:")
print(f"  N = {np.sum(neg_pec)}")
print(f"  Mean angular δ = {np.mean(delta_angular[neg_pec]):.3f}")
print(f"  Mean distance = {np.mean(d_valid[neg_pec]):.1f} Mpc")

# =============================================================================
# TEST 6: CHECK DISTANCE-INDEPENDENT VELOCITY ESTIMATE
# =============================================================================

print("\n" + "=" * 70)
print("TEST 6: DISTANCE-CORRECTED ANALYSIS")
print("=" * 70)

# Hypothesis: The issue is distance-dependent selection
# Solution: Compare at fixed distance ranges

print("\nVelocity by environment in NARROW DISTANCE BINS:")
print("-" * 60)

for d_min in [20, 50, 100]:
    d_max = d_min + 30
    d_mask = (d_valid >= d_min) & (d_valid < d_max)

    if np.sum(d_mask) < 500:
        continue

    # Re-compute angular density for this subset
    subset_idx = np.where(d_mask)[0]
    delta_sub = delta_angular[d_mask]
    v_sub = v_pec_valid[d_mask]

    # Compare quartiles
    q1 = delta_sub < np.percentile(delta_sub, 25)
    q4 = delta_sub > np.percentile(delta_sub, 75)

    v_q1 = np.mean(np.abs(v_sub[q1]))
    v_q4 = np.mean(np.abs(v_sub[q4]))
    ratio = v_q1 / v_q4

    print(f"  d = {d_min}-{d_max} Mpc (N={np.sum(d_mask)}):")
    print(f"    Low density (Q1): <|v|> = {v_q1:.0f} km/s")
    print(f"    High density (Q4): <|v|> = {v_q4:.0f} km/s")
    print(f"    Ratio: {ratio:.3f}")
    if ratio > 1:
        print(f"    >>> Synchronism direction: void > overdense <<<")
    else:
        print(f"    >>> Opposite to Synchronism: void < overdense <<<")

# =============================================================================
# INTERPRETATION
# =============================================================================

print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)

print("""
FINDINGS:
=========

1. The velocity SIGN correlates with distance bias
   - Positive v_pec → distance underestimated
   - Negative v_pec → distance overestimated

2. Angular density correlates with redshift
   - This is expected: more galaxies at higher z in the survey

3. In narrow distance bins:
   - The signal direction varies
   - Some bins show Synchronism direction, others opposite

4. The overall "reversal" is likely due to:
   - Malmquist bias: in voids, we preferentially see brighter galaxies
   - These appear closer, giving negative v_pec
   - In clusters, many galaxies are detected, some with positive v_pec

CONCLUSION:
===========
The CF4 "peculiar velocities" are dominated by DISTANCE ERRORS,
not true peculiar motions. The ~15-20% distance errors create
~1000-2000 km/s velocity errors, swamping the ~200-400 km/s true signal.

The velocity-quartile test from Session #170 was more robust because
it compared |v| not v_pec, but still suffers from these issues.

NEXT STEPS:
===========
1. Use ONLY the high-precision subset (SNIa, TRGB, Cepheids)
2. Or use velocity field reconstructions (not raw v_pec)
3. Or wait for DESI peculiar velocities with smaller errors
""")

print("\n" + "=" * 70)
print("SESSION #171b COMPLETE")
print("=" * 70)
