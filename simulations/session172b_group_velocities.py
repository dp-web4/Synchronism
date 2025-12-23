#!/usr/bin/env python3
"""
SESSION #172b: CF4 GROUP PECULIAR VELOCITIES ANALYSIS
======================================================
Date: December 23, 2025

Following Session #171 recommendation:
"Use velocity field reconstructions (not raw v_pec)"

CF4 table4.dat contains GROUP peculiar velocities:
- Vpec: Peculiar velocity using ramp in Equation 11
- These are averaged over group members, reducing noise
- May be more robust than individual galaxy v_pec

Also contains:
- Vpds: Davis & Scrimgeour (2014) formula
- Vpwf: Watkins & Feldman (2015) formula
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #172b: CF4 GROUP PECULIAR VELOCITIES")
print("=" * 70)

# =============================================================================
# LOAD GROUP DATA (table4.dat)
# =============================================================================

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table4.dat', 'r') as f:
    lines = f.readlines()

print(f"\nLoaded {len(lines)} groups from CF4 table4.dat")

# Parse according to ReadMe:
# DMzp: 9-14, e_DMzp: 16-20
# Dist: 22-26
# V3k: 40-44 (CMB velocity)
# Vpds: 52-57 (Davis & Scrimgeour peculiar velocity)
# Vpwf: 59-63 (Watkins & Feldman peculiar velocity)
# Vpec: 65-69 (ramp-corrected peculiar velocity)
# RAdeg: 84-91, DEdeg: 93-100

N = len(lines)
pgc = np.zeros(N, dtype=int)
dm = np.zeros(N)
dm_err = np.zeros(N)
dist = np.zeros(N)
v3k = np.zeros(N)
vpds = np.zeros(N)
vpwf = np.zeros(N)
vpec = np.zeros(N)
ra = np.zeros(N)
dec = np.zeros(N)
sgx = np.zeros(N)
sgy = np.zeros(N)
sgz = np.zeros(N)

def safe_float(s, default=0):
    s = s.strip()
    if not s or s == '?' or s == '-':
        return default
    try:
        return float(s)
    except:
        return default

def safe_int(s, default=0):
    s = s.strip()
    if not s or s == '?' or s == '-':
        return default
    try:
        return int(s)
    except:
        return default

for i, line in enumerate(lines):
    if len(line) < 157:
        continue

    pgc[i] = safe_int(line[0:7])
    dm[i] = safe_float(line[8:14])
    dm_err[i] = safe_float(line[15:20])
    dist[i] = safe_float(line[21:26])
    v3k[i] = safe_int(line[39:44])
    vpds[i] = safe_int(line[51:57])
    vpwf[i] = safe_int(line[58:63])
    vpec[i] = safe_int(line[64:69])
    ra[i] = safe_float(line[83:91])
    dec[i] = safe_float(line[92:100])
    sgx[i] = safe_int(line[137:143])
    sgy[i] = safe_int(line[144:150])
    sgz[i] = safe_int(line[151:157])

# Quality cuts
valid = (dist > 10) & (dist < 300) & (dm_err > 0) & (np.abs(vpec) < 10000)
N_valid = np.sum(valid)

print(f"Groups after quality cuts: {N_valid}")

# =============================================================================
# COMPARE DIFFERENT PECULIAR VELOCITY DEFINITIONS
# =============================================================================

print("\n" + "=" * 70)
print("PECULIAR VELOCITY DEFINITIONS COMPARISON")
print("=" * 70)

print(f"\nVpec (ramp, Eq. 11):")
print(f"  Mean: {np.mean(vpec[valid]):.0f} km/s")
print(f"  Std: {np.std(vpec[valid]):.0f} km/s")
print(f"  Range: [{np.min(vpec[valid]):.0f}, {np.max(vpec[valid]):.0f}] km/s")

print(f"\nVpds (Davis & Scrimgeour 2014):")
print(f"  Mean: {np.mean(vpds[valid]):.0f} km/s")
print(f"  Std: {np.std(vpds[valid]):.0f} km/s")

print(f"\nVpwf (Watkins & Feldman 2015):")
print(f"  Mean: {np.mean(vpwf[valid]):.0f} km/s")
print(f"  Std: {np.std(vpwf[valid]):.0f} km/s")

# Correlation between methods
corr_pds_pec = np.corrcoef(vpds[valid], vpec[valid])[0,1]
corr_pwf_pec = np.corrcoef(vpwf[valid], vpec[valid])[0,1]
print(f"\nCorrelation Vpds vs Vpec: {corr_pds_pec:.3f}")
print(f"Correlation Vpwf vs Vpec: {corr_pwf_pec:.3f}")

# =============================================================================
# ENVIRONMENT FROM ANGULAR DENSITY
# =============================================================================

print("\n" + "=" * 70)
print("ENVIRONMENT ESTIMATION (ANGULAR DENSITY)")
print("=" * 70)

ra_v = ra[valid]
dec_v = dec[valid]
vpec_v = vpec[valid]
dist_v = dist[valid]
dm_err_v = dm_err[valid]

# Angular density
ra_rad = np.radians(ra_v)
dec_rad = np.radians(dec_v)
x_sphere = np.cos(dec_rad) * np.cos(ra_rad)
y_sphere = np.cos(dec_rad) * np.sin(ra_rad)
z_sphere = np.sin(dec_rad)

coords_sphere = np.column_stack([x_sphere, y_sphere, z_sphere])
tree = cKDTree(coords_sphere)

# 10-degree radius
chord_10deg = 2 * np.sin(np.radians(10/2))
counts_angular = np.array(tree.query_ball_point(coords_sphere, r=chord_10deg, return_length=True)) - 1

mean_count = np.mean(counts_angular)
delta_angular = (counts_angular - mean_count) / max(mean_count, 1)

print(f"\nAngular overdensity (10° radius):")
print(f"  Mean neighbors: {mean_count:.1f}")
print(f"  δ range: [{np.min(delta_angular):.2f}, {np.max(delta_angular):.2f}]")

# =============================================================================
# TEST 1: VELOCITY-QUARTILE TEST (GROUP VELOCITIES)
# =============================================================================

print("\n" + "=" * 70)
print("TEST 1: VELOCITY-QUARTILE TEST (GROUPS)")
print("=" * 70)

abs_v = np.abs(vpec_v)
q1_thresh = np.percentile(abs_v, 25)
q4_thresh = np.percentile(abs_v, 75)

q1_mask = abs_v < q1_thresh
q4_mask = abs_v > q4_thresh

delta_q1 = delta_angular[q1_mask]
delta_q4 = delta_angular[q4_mask]

print(f"\n|Vpec| thresholds:")
print(f"  Q1 (low): < {q1_thresh:.0f} km/s (N={np.sum(q1_mask)})")
print(f"  Q4 (high): > {q4_thresh:.0f} km/s (N={np.sum(q4_mask)})")

print(f"\nEnvironment density by velocity quartile:")
print(f"  Q1 (low |v|): mean δ = {np.mean(delta_q1):.3f}")
print(f"  Q4 (high |v|): mean δ = {np.mean(delta_q4):.3f}")

u_stat, p_mw = stats.mannwhitneyu(delta_q1, delta_q4, alternative='greater')
print(f"\nMann-Whitney U test:")
print(f"  p-value: {p_mw:.2e}")

if np.mean(delta_q1) > np.mean(delta_q4):
    print(f"\n>>> High-|v| groups are in LOWER density environments")
    print(f">>> CONSISTENT with Synchronism prediction!")
else:
    print(f"\n>>> High-|v| groups are in HIGHER density environments")
    print(f">>> OPPOSITE to Synchronism prediction")

# =============================================================================
# TEST 2: VELOCITY BY ENVIRONMENT
# =============================================================================

print("\n" + "=" * 70)
print("TEST 2: VELOCITY BY ENVIRONMENT (GROUPS)")
print("=" * 70)

void_mask = delta_angular < np.percentile(delta_angular, 25)
overdense_mask = delta_angular > np.percentile(delta_angular, 75)

v_void = np.abs(vpec_v[void_mask])
v_overdense = np.abs(vpec_v[overdense_mask])

print(f"\n|Vpec| by environment:")
print(f"  Void (Q1, N={np.sum(void_mask)}): <|v|> = {np.mean(v_void):.0f} ± {np.std(v_void)/np.sqrt(len(v_void)):.0f} km/s")
print(f"  Overdense (Q4, N={np.sum(overdense_mask)}): <|v|> = {np.mean(v_overdense):.0f} ± {np.std(v_overdense)/np.sqrt(len(v_overdense)):.0f} km/s")
print(f"  Ratio (void/overdense): {np.mean(v_void)/np.mean(v_overdense):.3f}")

print(f"\n  Synchronism predicts: ratio > 1")
if np.mean(v_void) > np.mean(v_overdense):
    print(f"  >>> OBSERVED: Void > Overdense (SYNCHRONISM DIRECTION!)")
else:
    print(f"  >>> OBSERVED: Void < Overdense (opposite)")

# =============================================================================
# TEST 3: DISTANCE-BINNED ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("TEST 3: DISTANCE-BINNED ANALYSIS (GROUPS)")
print("=" * 70)

print("\nVelocity by environment in distance bins:")
print("-" * 60)

results_by_dist = []

for d_min in [10, 30, 50, 70, 100, 150, 200]:
    d_max = d_min + 50
    d_mask = (dist_v >= d_min) & (dist_v < d_max)

    if np.sum(d_mask) < 100:
        continue

    delta_sub = delta_angular[d_mask]
    v_sub = vpec_v[d_mask]

    void_sub = delta_sub < np.percentile(delta_sub, 25)
    overdense_sub = delta_sub > np.percentile(delta_sub, 75)

    if np.sum(void_sub) < 10 or np.sum(overdense_sub) < 10:
        continue

    v_void_sub = np.mean(np.abs(v_sub[void_sub]))
    v_overdense_sub = np.mean(np.abs(v_sub[overdense_sub]))
    ratio = v_void_sub / v_overdense_sub

    results_by_dist.append((d_min, d_max, np.sum(d_mask), ratio))

    print(f"  d = {d_min}-{d_max} Mpc (N={np.sum(d_mask)}):")
    print(f"    Void (N={np.sum(void_sub)}): <|v|> = {v_void_sub:.0f} km/s")
    print(f"    Overdense (N={np.sum(overdense_sub)}): <|v|> = {v_overdense_sub:.0f} km/s")
    print(f"    Ratio: {ratio:.3f}", end="")
    if ratio > 1:
        print(" *** SYNCHRONISM ***")
    else:
        print()

# =============================================================================
# TEST 4: 3D DENSITY FROM SUPERGALACTIC COORDINATES
# =============================================================================

print("\n" + "=" * 70)
print("TEST 4: 3D DENSITY FROM SUPERGALACTIC COORDINATES")
print("=" * 70)

# Use SGX, SGY, SGZ (in Mpc units from velocities) for true 3D density
sgx_v = sgx[valid]
sgy_v = sgy[valid]
sgz_v = sgz[valid]

# Convert to proper distance using H0
H0 = 74.6
sgx_mpc = sgx_v / H0
sgy_mpc = sgy_v / H0
sgz_mpc = sgz_v / H0

coords_3d = np.column_stack([sgx_mpc, sgy_mpc, sgz_mpc])
tree_3d = cKDTree(coords_3d)

# Count neighbors within 10 Mpc
counts_3d = np.array(tree_3d.query_ball_point(coords_3d, r=10, return_length=True)) - 1
mean_3d = np.mean(counts_3d)
delta_3d = (counts_3d - mean_3d) / max(mean_3d, 1)

print(f"\n3D overdensity (10 Mpc radius):")
print(f"  Mean neighbors: {mean_3d:.1f}")
print(f"  δ range: [{np.min(delta_3d):.2f}, {np.max(delta_3d):.2f}]")

# Velocity by 3D environment
void_3d = delta_3d < np.percentile(delta_3d, 25)
overdense_3d = delta_3d > np.percentile(delta_3d, 75)

v_void_3d = np.mean(np.abs(vpec_v[void_3d]))
v_overdense_3d = np.mean(np.abs(vpec_v[overdense_3d]))

print(f"\n|Vpec| by 3D environment:")
print(f"  Void (Q1): <|v|> = {v_void_3d:.0f} km/s")
print(f"  Overdense (Q4): <|v|> = {v_overdense_3d:.0f} km/s")
print(f"  Ratio: {v_void_3d/v_overdense_3d:.3f}")

# =============================================================================
# TEST 5: DISTANCE ERROR ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("TEST 5: VELOCITY ERROR ESTIMATION")
print("=" * 70)

# σ_v ≈ H0 * d * σ_dm * ln(10)/5
v_err = H0 * dist_v * dm_err_v * 0.461

print(f"\nVelocity error statistics:")
print(f"  Mean σ_v: {np.mean(v_err):.0f} km/s")
print(f"  Median σ_v: {np.median(v_err):.0f} km/s")
print(f"  Mean |Vpec|: {np.mean(np.abs(vpec_v)):.0f} km/s")
print(f"  Signal-to-noise: {np.mean(np.abs(vpec_v))/np.mean(v_err):.2f}")

# Error-weighted analysis
weights = 1 / v_err**2
weights = weights / np.sum(weights)

w_void = weights[void_mask] / np.sum(weights[void_mask])
w_overdense = weights[overdense_mask] / np.sum(weights[overdense_mask])

weighted_v_void = np.sum(np.abs(vpec_v[void_mask]) * w_void)
weighted_v_overdense = np.sum(np.abs(vpec_v[overdense_mask]) * w_overdense)

print(f"\nError-weighted analysis:")
print(f"  Void: {weighted_v_void:.0f} km/s")
print(f"  Overdense: {weighted_v_overdense:.0f} km/s")
print(f"  Weighted ratio: {weighted_v_void/weighted_v_overdense:.3f}")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Velocity distributions
ax1 = axes[0, 0]
ax1.hist(vpec_v, bins=100, alpha=0.7, label='Vpec (ramp)')
ax1.axvline(0, color='red', linestyle='--')
ax1.set_xlabel('Peculiar velocity (km/s)')
ax1.set_ylabel('Count')
ax1.set_title(f'Group Peculiar Velocity Distribution (N={N_valid})')
ax1.legend()
ax1.set_xlim(-5000, 5000)

# Panel 2: Distance-binned ratios
ax2 = axes[0, 1]
if results_by_dist:
    d_centers = [(r[0]+r[1])/2 for r in results_by_dist]
    ratios = [r[3] for r in results_by_dist]
    ax2.plot(d_centers, ratios, 'o-', markersize=10)
    ax2.axhline(1.0, color='red', linestyle='--', label='Synchronism threshold')
    ax2.set_xlabel('Distance (Mpc)')
    ax2.set_ylabel('Void/Overdense velocity ratio')
    ax2.set_title('Velocity Ratio by Distance')
    ax2.legend()
    ax2.set_ylim(0.5, 1.5)

# Panel 3: Velocity-environment scatter
ax3 = axes[1, 0]
scatter = ax3.scatter(delta_angular, np.abs(vpec_v), c=dist_v, alpha=0.3, s=5, cmap='viridis')
plt.colorbar(scatter, ax=ax3, label='Distance (Mpc)')
ax3.set_xlabel('Angular overdensity δ')
ax3.set_ylabel('|Peculiar velocity| (km/s)')
ax3.set_title('Velocity vs Environment')

# Panel 4: Velocity quartile comparison
ax4 = axes[1, 1]
quartiles = ['Q1 (low |v|)', 'Q4 (high |v|)']
densities = [np.mean(delta_q1), np.mean(delta_q4)]
errors = [np.std(delta_q1)/np.sqrt(len(delta_q1)), np.std(delta_q4)/np.sqrt(len(delta_q4))]
colors = ['red', 'blue']
ax4.bar(quartiles, densities, yerr=errors, capsize=5, color=colors, alpha=0.7)
ax4.axhline(0, color='black', linestyle='-', linewidth=0.5)
ax4.set_ylabel('Mean angular overdensity δ')
ax4.set_title(f'Environment by Velocity Quartile')

plt.suptitle('Session #172b: CF4 Group Peculiar Velocities', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session172b_group_velocities.png', dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FIGURE SAVED: session172b_group_velocities.png")
print("=" * 70)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #172b: SUMMARY")
print("=" * 70)

print(f"""
CF4 GROUP PECULIAR VELOCITY ANALYSIS
=====================================

Sample: {N_valid} galaxy groups with quality cuts

1. VELOCITY-QUARTILE TEST:
   - High-|v| groups: mean δ = {np.mean(delta_q4):.3f}
   - Low-|v| groups: mean δ = {np.mean(delta_q1):.3f}
   - Direction: {"SYNCHRONISM" if np.mean(delta_q1) > np.mean(delta_q4) else "Opposite"}

2. VELOCITY BY ENVIRONMENT (Angular):
   - Void <|v|>: {np.mean(v_void):.0f} km/s
   - Overdense <|v|>: {np.mean(v_overdense):.0f} km/s
   - Ratio: {np.mean(v_void)/np.mean(v_overdense):.3f}

3. VELOCITY BY ENVIRONMENT (3D):
   - Void <|v|>: {v_void_3d:.0f} km/s
   - Overdense <|v|>: {v_overdense_3d:.0f} km/s
   - Ratio: {v_void_3d/v_overdense_3d:.3f}

4. KEY FINDING:
   The group peculiar velocities show {"SIMILAR" if np.mean(v_void)/np.mean(v_overdense) < 1 else "REVERSED"}
   behavior to individual galaxies: void groups have
   {"LOWER" if np.mean(v_void)/np.mean(v_overdense) < 1 else "HIGHER"} peculiar velocities.

INTERPRETATION:
===============
The apparent "reversal" from Synchronism predictions persists in the group
velocity data. This suggests either:
1. The angular density metric doesn't capture true voids
2. Selection effects persist in the group averaging
3. The Synchronism effect is smaller than expected at these scales
4. Alternative explanations for the observed velocity-density correlations
""")

print("=" * 70)
print("SESSION #172b COMPLETE")
print("=" * 70)
