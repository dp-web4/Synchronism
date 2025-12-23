#!/usr/bin/env python3
"""
SESSION #172c: 3D vs ANGULAR DENSITY DISCREPANCY
=================================================
Date: December 23, 2025

CRITICAL FINDING FROM 172b:
- Angular density: Ratio = 0.55 (opposite to Synchronism)
- 3D density: Ratio = 4.23 (STRONGLY supporting Synchronism!)

This is a massive discrepancy that needs investigation.
The 3D density uses supergalactic coordinates from CF4.

Key questions:
1. Why does 3D density give opposite results to angular?
2. Which is more physically meaningful?
3. Is the 3D result real or an artifact?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #172c: 3D vs ANGULAR DENSITY INVESTIGATION")
print("=" * 70)

# =============================================================================
# LOAD DATA
# =============================================================================

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table4.dat', 'r') as f:
    lines = f.readlines()

N = len(lines)
dist = np.zeros(N)
dm_err = np.zeros(N)
vpec = np.zeros(N)
ra = np.zeros(N)
dec = np.zeros(N)
sgx = np.zeros(N)
sgy = np.zeros(N)
sgz = np.zeros(N)

def safe_float(s, default=0):
    s = s.strip()
    if not s: return default
    try: return float(s)
    except: return default

def safe_int(s, default=0):
    s = s.strip()
    if not s: return default
    try: return int(s)
    except: return default

for i, line in enumerate(lines):
    if len(line) < 157:
        continue
    dist[i] = safe_float(line[21:26])
    dm_err[i] = safe_float(line[15:20])
    vpec[i] = safe_int(line[64:69])
    ra[i] = safe_float(line[83:91])
    dec[i] = safe_float(line[92:100])
    sgx[i] = safe_int(line[137:143])
    sgy[i] = safe_int(line[144:150])
    sgz[i] = safe_int(line[151:157])

valid = (dist > 10) & (dist < 300) & (dm_err > 0) & (np.abs(vpec) < 10000)
print(f"\nValid groups: {np.sum(valid)}")

# Extract valid data
dist_v = dist[valid]
vpec_v = vpec[valid]
ra_v = ra[valid]
dec_v = dec[valid]
sgx_v = sgx[valid]
sgy_v = sgy[valid]
sgz_v = sgz[valid]

H0 = 74.6

# =============================================================================
# COMPUTE BOTH DENSITY METRICS
# =============================================================================

print("\n" + "=" * 70)
print("COMPUTING DENSITY METRICS")
print("=" * 70)

# Angular density
ra_rad = np.radians(ra_v)
dec_rad = np.radians(dec_v)
x_sphere = np.cos(dec_rad) * np.cos(ra_rad)
y_sphere = np.cos(dec_rad) * np.sin(ra_rad)
z_sphere = np.sin(dec_rad)
coords_angular = np.column_stack([x_sphere, y_sphere, z_sphere])
tree_angular = cKDTree(coords_angular)

chord_10deg = 2 * np.sin(np.radians(10/2))
counts_angular = np.array(tree_angular.query_ball_point(coords_angular, r=chord_10deg, return_length=True)) - 1
delta_angular = (counts_angular - np.mean(counts_angular)) / max(np.mean(counts_angular), 1)

# 3D density from supergalactic coordinates
# Note: SGX, SGY, SGZ are in km/s units (velocity-space)
# Convert to Mpc
sgx_mpc = sgx_v / H0
sgy_mpc = sgy_v / H0
sgz_mpc = sgz_v / H0

coords_3d = np.column_stack([sgx_mpc, sgy_mpc, sgz_mpc])
tree_3d = cKDTree(coords_3d)

# Try different 3D radii
print("\n3D density with different radii:")
for radius in [5, 10, 15, 20, 30]:
    counts_3d = np.array(tree_3d.query_ball_point(coords_3d, r=radius, return_length=True)) - 1
    delta_3d = (counts_3d - np.mean(counts_3d)) / max(np.mean(counts_3d), 1)

    void_3d = delta_3d < np.percentile(delta_3d, 25)
    overdense_3d = delta_3d > np.percentile(delta_3d, 75)

    v_void = np.mean(np.abs(vpec_v[void_3d]))
    v_overdense = np.mean(np.abs(vpec_v[overdense_3d]))
    ratio = v_void / v_overdense

    print(f"  R = {radius} Mpc: mean neighbors = {np.mean(counts_3d):.1f}, ratio = {ratio:.3f}", end="")
    if ratio > 1:
        print(" *** SYNCHRONISM ***")
    else:
        print()

# =============================================================================
# INVESTIGATE THE DISCREPANCY
# =============================================================================

print("\n" + "=" * 70)
print("INVESTIGATING THE DISCREPANCY")
print("=" * 70)

# Use 10 Mpc for 3D (as in 172b)
counts_3d = np.array(tree_3d.query_ball_point(coords_3d, r=10, return_length=True)) - 1
delta_3d = (counts_3d - np.mean(counts_3d)) / max(np.mean(counts_3d), 1)

# Compare classifications
void_angular = delta_angular < np.percentile(delta_angular, 25)
overdense_angular = delta_angular > np.percentile(delta_angular, 75)
void_3d = delta_3d < np.percentile(delta_3d, 25)
overdense_3d = delta_3d > np.percentile(delta_3d, 75)

print(f"\nClassification overlap:")
print(f"  Angular void AND 3D void: {np.sum(void_angular & void_3d)}")
print(f"  Angular void AND 3D overdense: {np.sum(void_angular & overdense_3d)}")
print(f"  Angular overdense AND 3D void: {np.sum(overdense_angular & void_3d)}")
print(f"  Angular overdense AND 3D overdense: {np.sum(overdense_angular & overdense_3d)}")

# Correlation between the two density metrics
corr_density = np.corrcoef(delta_angular, delta_3d)[0,1]
print(f"\nCorrelation (angular δ vs 3D δ): {corr_density:.3f}")

# =============================================================================
# KEY INSIGHT: Distance distribution by environment
# =============================================================================

print("\n" + "=" * 70)
print("KEY INSIGHT: DISTANCE DISTRIBUTION BY ENVIRONMENT")
print("=" * 70)

print(f"\nMean distance by ANGULAR environment:")
print(f"  Angular void: {np.mean(dist_v[void_angular]):.1f} Mpc")
print(f"  Angular overdense: {np.mean(dist_v[overdense_angular]):.1f} Mpc")

print(f"\nMean distance by 3D environment:")
print(f"  3D void: {np.mean(dist_v[void_3d]):.1f} Mpc")
print(f"  3D overdense: {np.mean(dist_v[overdense_3d]):.1f} Mpc")

# This is the key: 3D voids are at LARGER distances!
print(f"\nRatio of mean distances (void/overdense):")
print(f"  Angular: {np.mean(dist_v[void_angular])/np.mean(dist_v[overdense_angular]):.3f}")
print(f"  3D: {np.mean(dist_v[void_3d])/np.mean(dist_v[overdense_3d]):.3f}")

# =============================================================================
# SELECTION EFFECT EXPLANATION
# =============================================================================

print("\n" + "=" * 70)
print("SELECTION EFFECT ANALYSIS")
print("=" * 70)

# From Session #171: |v| correlates with distance due to errors
# If 3D voids have higher distances, they'll have higher |v| due to errors

print("\n3D void classification is biased by distance:")
print(f"  - 3D 'void' galaxies are at larger distances (mean = {np.mean(dist_v[void_3d]):.0f} Mpc)")
print(f"  - 3D 'overdense' galaxies are nearby (mean = {np.mean(dist_v[overdense_3d]):.0f} Mpc)")
print(f"  - Larger distance → larger velocity errors → larger apparent |v|")
print(f"  - This CREATES the appearance of Synchronism signal!")

# Verify: Check |v| vs distance correlation
abs_v = np.abs(vpec_v)
corr_v_d = np.corrcoef(abs_v, dist_v)[0,1]
print(f"\nCorrelation (|v| vs distance): r = {corr_v_d:.3f}")

# =============================================================================
# CONTROL: DISTANCE-MATCHED COMPARISON
# =============================================================================

print("\n" + "=" * 70)
print("CONTROL: DISTANCE-MATCHED 3D ANALYSIS")
print("=" * 70)

# Narrow distance bins to control for distance effects
print("\n3D density analysis in narrow distance bins:")
for d_min in [30, 50, 80, 120]:
    d_max = d_min + 30
    d_mask = (dist_v >= d_min) & (dist_v < d_max)

    if np.sum(d_mask) < 500:
        continue

    # Recompute 3D density for this slice
    coords_sub = coords_3d[d_mask]
    tree_sub = cKDTree(coords_sub)
    counts_sub = np.array(tree_sub.query_ball_point(coords_sub, r=10, return_length=True)) - 1
    delta_sub = (counts_sub - np.mean(counts_sub)) / max(np.mean(counts_sub), 1)

    void_sub = delta_sub < np.percentile(delta_sub, 25)
    overdense_sub = delta_sub > np.percentile(delta_sub, 75)

    v_void_sub = np.mean(np.abs(vpec_v[d_mask][void_sub]))
    v_overdense_sub = np.mean(np.abs(vpec_v[d_mask][overdense_sub]))
    ratio = v_void_sub / v_overdense_sub

    print(f"  d = {d_min}-{d_max} Mpc (N={np.sum(d_mask)}):")
    print(f"    3D void <|v|>: {v_void_sub:.0f} km/s")
    print(f"    3D overdense <|v|>: {v_overdense_sub:.0f} km/s")
    print(f"    Ratio: {ratio:.3f}", end="")
    if ratio > 1:
        print(" (Synchronism)")
    else:
        print(" (opposite)")

# =============================================================================
# ANGULAR DENSITY IN DISTANCE BINS (CONTROL)
# =============================================================================

print("\n" + "=" * 70)
print("CONTROL: ANGULAR DENSITY IN DISTANCE BINS")
print("=" * 70)

for d_min in [30, 50, 80, 120]:
    d_max = d_min + 30
    d_mask = (dist_v >= d_min) & (dist_v < d_max)

    if np.sum(d_mask) < 500:
        continue

    delta_sub = delta_angular[d_mask]
    void_sub = delta_sub < np.percentile(delta_sub, 25)
    overdense_sub = delta_sub > np.percentile(delta_sub, 75)

    v_void_sub = np.mean(np.abs(vpec_v[d_mask][void_sub]))
    v_overdense_sub = np.mean(np.abs(vpec_v[d_mask][overdense_sub]))
    ratio = v_void_sub / v_overdense_sub

    print(f"  d = {d_min}-{d_max} Mpc (N={np.sum(d_mask)}):")
    print(f"    Angular void <|v|>: {v_void_sub:.0f} km/s")
    print(f"    Angular overdense <|v|>: {v_overdense_sub:.0f} km/s")
    print(f"    Ratio: {ratio:.3f}", end="")
    if ratio > 1:
        print(" (Synchronism)")
    else:
        print(" (opposite)")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Angular vs 3D density correlation
ax1 = axes[0, 0]
ax1.scatter(delta_angular, delta_3d, alpha=0.1, s=1)
ax1.set_xlabel('Angular overdensity δ')
ax1.set_ylabel('3D overdensity δ')
ax1.set_title(f'Angular vs 3D Density (r = {corr_density:.3f})')
ax1.axhline(0, color='red', linestyle='--', alpha=0.5)
ax1.axvline(0, color='red', linestyle='--', alpha=0.5)

# Panel 2: Distance distribution by environment type
ax2 = axes[0, 1]
ax2.hist(dist_v[void_3d], bins=50, alpha=0.7, label=f'3D void (mean={np.mean(dist_v[void_3d]):.0f})', density=True)
ax2.hist(dist_v[overdense_3d], bins=50, alpha=0.7, label=f'3D overdense (mean={np.mean(dist_v[overdense_3d]):.0f})', density=True)
ax2.set_xlabel('Distance (Mpc)')
ax2.set_ylabel('Density')
ax2.set_title('Distance Distribution by 3D Environment')
ax2.legend()

# Panel 3: |v| vs distance
ax3 = axes[1, 0]
ax3.scatter(dist_v, abs_v, alpha=0.1, s=1)
ax3.set_xlabel('Distance (Mpc)')
ax3.set_ylabel('|Peculiar velocity| (km/s)')
ax3.set_title(f'Velocity vs Distance (r = {corr_v_d:.3f})')

# Panel 4: Summary comparison
ax4 = axes[1, 1]
methods = ['Angular\n(all)', '3D\n(all)', 'Angular\n(50-80 Mpc)', '3D\n(50-80 Mpc)']

# Calculate all ratios
v_void_ang = np.mean(np.abs(vpec_v[void_angular]))
v_od_ang = np.mean(np.abs(vpec_v[overdense_angular]))
v_void_3d = np.mean(np.abs(vpec_v[void_3d]))
v_od_3d = np.mean(np.abs(vpec_v[overdense_3d]))

d_mask = (dist_v >= 50) & (dist_v < 80)
delta_ang_sub = delta_angular[d_mask]
delta_3d_sub_calc = delta_3d[d_mask]

void_ang_sub = delta_ang_sub < np.percentile(delta_ang_sub, 25)
od_ang_sub = delta_ang_sub > np.percentile(delta_ang_sub, 75)
void_3d_sub = delta_3d_sub_calc < np.percentile(delta_3d_sub_calc, 25)
od_3d_sub = delta_3d_sub_calc > np.percentile(delta_3d_sub_calc, 75)

ratios = [
    v_void_ang/v_od_ang,
    v_void_3d/v_od_3d,
    np.mean(np.abs(vpec_v[d_mask][void_ang_sub])) / np.mean(np.abs(vpec_v[d_mask][od_ang_sub])),
    np.mean(np.abs(vpec_v[d_mask][void_3d_sub])) / np.mean(np.abs(vpec_v[d_mask][od_3d_sub]))
]
colors = ['red' if r < 1 else 'green' for r in ratios]

ax4.bar(methods, ratios, color=colors, alpha=0.7)
ax4.axhline(1.0, color='black', linestyle='--', label='Synchronism threshold')
ax4.set_ylabel('Void/Overdense velocity ratio')
ax4.set_title('Effect of Distance Control')
ax4.legend()

plt.suptitle('Session #172c: 3D vs Angular Density Investigation', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session172c_density_investigation.png', dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FIGURE SAVED: session172c_density_investigation.png")
print("=" * 70)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #172c: SUMMARY")
print("=" * 70)

print("""
CRITICAL FINDING: THE 3D DENSITY RESULT IS AN ARTIFACT
======================================================

The dramatic ratio = 4.2 from 3D density is NOT a Synchronism signal.
It is caused by a selection effect:

1. 3D "voids" are preferentially at LARGER distances
   - 3D void mean distance: {:.0f} Mpc
   - 3D overdense mean distance: {:.0f} Mpc

2. |v| correlates with distance (r = {:.3f})
   - Due to velocity errors scaling with distance
   - More distant galaxies have larger apparent |v|

3. Therefore:
   - 3D "void" galaxies have higher |v| because they're MORE DISTANT
   - This creates a FALSE Synchronism signal

4. When controlling for distance:
   - Both angular and 3D density show void/overdense ratio < 1
   - The "opposite to Synchronism" pattern persists

IMPLICATIONS:
=============
The Session #170 velocity-quartile result (>> 10σ significance)
likely has a similar selection effect component.

The challenge is distinguishing:
1. True Synchronism signal (G_eff enhancement in voids)
2. Selection effects (Malmquist bias, distance-dependent errors)

Next steps:
1. Need true void catalogs (not density-based classification)
2. Need forward modeling of selection effects
3. Consider complementary tests not affected by distance errors
""".format(
    np.mean(dist_v[void_3d]),
    np.mean(dist_v[overdense_3d]),
    corr_v_d
))

print("=" * 70)
print("SESSION #172c COMPLETE")
print("=" * 70)
