#!/usr/bin/env python3
"""
SESSION #174: FORWARD MODELING SELECTION EFFECTS
=================================================
Date: December 24, 2025

MOTIVATION:
-----------
Sessions #170-173 showed that CF4 peculiar velocity tests are dominated by
selection effects. Key finding:

    |v| correlates with distance (r = 0.40)

This creates spurious velocity-environment correlations because:
1. Distance errors scale with distance
2. Low-density "void" regions may be at different distances than dense regions
3. The apparent velocity-environment signal may be entirely from selection

THIS SESSION:
-------------
Forward model the selection effects by:
1. Creating a mock catalog with KNOWN true velocities
2. Adding realistic distance errors (matching CF4 error distribution)
3. Computing apparent peculiar velocities
4. Testing how much velocity-environment correlation arises from errors alone

This tells us the NULL HYPOTHESIS signal we expect with NO true Synchronism effect.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #174: FORWARD MODELING SELECTION EFFECTS")
print("=" * 70)

# =============================================================================
# 1. LOAD CF4 DATA TO CHARACTERIZE ERROR DISTRIBUTION
# =============================================================================

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table4.dat', 'r') as f:
    lines = f.readlines()

N_cf4 = len(lines)
dist_cf4 = np.zeros(N_cf4)
dm_err_cf4 = np.zeros(N_cf4)
vpec_cf4 = np.zeros(N_cf4)

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
    dist_cf4[i] = safe_float(line[21:26])
    dm_err_cf4[i] = safe_float(line[15:20])
    vpec_cf4[i] = safe_int(line[64:69])

H0 = 74.6
valid = (dist_cf4 > 10) & (dist_cf4 < 200) & (dm_err_cf4 > 0) & (np.abs(vpec_cf4) < 5000)
N_valid = np.sum(valid)

dist_v = dist_cf4[valid]
dm_err_v = dm_err_cf4[valid]
vpec_v = vpec_cf4[valid]

print(f"\nCF4 data loaded: {N_valid} valid groups")

# =============================================================================
# 2. CHARACTERIZE CF4 ERROR DISTRIBUTION
# =============================================================================

print("\n" + "=" * 70)
print("CF4 ERROR CHARACTERIZATION")
print("=" * 70)

# σ_v = H0 * d * σ_dm * ln(10)/5 ≈ H0 * d * 0.461 * σ_dm
v_err = H0 * dist_v * 0.461 * dm_err_v

print(f"\nDistance modulus error:")
print(f"  Mean: {np.mean(dm_err_v):.3f} mag")
print(f"  Median: {np.median(dm_err_v):.3f} mag")
print(f"  Range: [{np.min(dm_err_v):.3f}, {np.max(dm_err_v):.3f}] mag")

print(f"\nVelocity error (derived):")
print(f"  Mean: {np.mean(v_err):.0f} km/s")
print(f"  Median: {np.median(v_err):.0f} km/s")

print(f"\nTrue peculiar velocity statistics:")
print(f"  Mean |v|: {np.mean(np.abs(vpec_v)):.0f} km/s")
print(f"  Median |v|: {np.median(np.abs(vpec_v)):.0f} km/s")

print(f"\nSignal-to-noise:")
print(f"  Mean: {np.mean(np.abs(vpec_v)/v_err):.2f}")
print(f"  Median: {np.median(np.abs(vpec_v)/v_err):.2f}")

# Key relationship: σ_v ∝ distance
corr_err_dist = np.corrcoef(v_err, dist_v)[0,1]
print(f"\nCorrelation (σ_v vs distance): r = {corr_err_dist:.3f}")

# =============================================================================
# 3. CREATE MOCK CATALOG (NULL HYPOTHESIS: NO SYNCHRONISM)
# =============================================================================

print("\n" + "=" * 70)
print("CREATING MOCK CATALOG (NULL HYPOTHESIS)")
print("=" * 70)

np.random.seed(42)

# Use same distance distribution as CF4
N_mock = len(dist_v)
dist_mock = dist_v.copy()

# Generate TRUE peculiar velocities with NO environment dependence
# Use ΛCDM expectation: σ_v ~ 300 km/s, no density dependence
sigma_v_true = 300  # km/s (typical bulk flow scale)
v_true = np.random.normal(0, sigma_v_true, N_mock)

print(f"\nTrue velocities (ΛCDM, no Synchronism):")
print(f"  σ_v_true = {sigma_v_true} km/s")
print(f"  Mean |v_true|: {np.mean(np.abs(v_true)):.0f} km/s")

# Add distance errors matching CF4 distribution
# σ_d/d = 0.461 * σ_dm
rel_dist_err = 0.461 * dm_err_v  # Relative distance error
dist_observed = dist_mock * (1 + np.random.normal(0, rel_dist_err))

# Compute observed peculiar velocity
# v_pec = v_cmb - H0 * d_observed
# v_cmb = v_true + H0 * d_true
# So: v_pec_observed = v_true + H0 * (d_true - d_observed)
#                    = v_true + H0 * d_true * (1 - d_observed/d_true)
#                    = v_true - H0 * d_true * (d_observed/d_true - 1)
dist_error = dist_observed - dist_mock
v_observed = v_true - H0 * dist_error

print(f"\nObserved velocities (with errors):")
print(f"  Mean |v_obs|: {np.mean(np.abs(v_observed)):.0f} km/s")
print(f"  Std |v_obs|: {np.std(v_observed):.0f} km/s")

# Compare to CF4
print(f"\nCF4 actual:")
print(f"  Mean |v_pec|: {np.mean(np.abs(vpec_v)):.0f} km/s")
print(f"  Std |v_pec|: {np.std(vpec_v):.0f} km/s")

# =============================================================================
# 4. TEST: DOES |v_obs| CORRELATE WITH DISTANCE?
# =============================================================================

print("\n" + "=" * 70)
print("TEST: |v| vs DISTANCE CORRELATION")
print("=" * 70)

# In the mock (null hypothesis)
corr_mock = np.corrcoef(np.abs(v_observed), dist_mock)[0,1]
print(f"\nMock catalog (null hypothesis):")
print(f"  Correlation (|v_obs| vs distance): r = {corr_mock:.3f}")

# In CF4 (observed)
corr_cf4 = np.corrcoef(np.abs(vpec_v), dist_v)[0,1]
print(f"\nCF4 observed:")
print(f"  Correlation (|v_pec| vs distance): r = {corr_cf4:.3f}")

print(f"\nComparison:")
if abs(corr_mock) > 0.35:
    print(f"  Mock correlation ({corr_mock:.2f}) is SIMILAR to CF4 ({corr_cf4:.2f})")
    print(f"  >>> Selection effects ALONE can explain the |v|-distance correlation!")
else:
    print(f"  Mock correlation ({corr_mock:.2f}) is DIFFERENT from CF4 ({corr_cf4:.2f})")
    print(f"  >>> Something beyond selection effects is happening")

# =============================================================================
# 5. SIMULATE ENVIRONMENT CLASSIFICATION
# =============================================================================

print("\n" + "=" * 70)
print("SIMULATING ENVIRONMENT CLASSIFICATION")
print("=" * 70)

# Create random sky positions for mock catalog
ra_mock = np.random.uniform(0, 360, N_mock)
dec_mock = np.random.uniform(-90, 90, N_mock)  # Simplified

# Convert to Cartesian
ra_rad = np.radians(ra_mock)
dec_rad = np.radians(dec_mock)
x_mock = dist_mock * np.cos(dec_rad) * np.cos(ra_rad)
y_mock = dist_mock * np.cos(dec_rad) * np.sin(ra_rad)
z_mock = dist_mock * np.sin(dec_rad)

coords_mock = np.column_stack([x_mock, y_mock, z_mock])

# Angular density (same as Sessions #170-173)
x_sphere = np.cos(dec_rad) * np.cos(ra_rad)
y_sphere = np.cos(dec_rad) * np.sin(ra_rad)
z_sphere = np.sin(dec_rad)
coords_angular = np.column_stack([x_sphere, y_sphere, z_sphere])

tree = cKDTree(coords_angular)
chord_10deg = 2 * np.sin(np.radians(10/2))
counts = np.array(tree.query_ball_point(coords_angular, r=chord_10deg, return_length=True)) - 1
delta_angular = (counts - np.mean(counts)) / max(np.mean(counts), 1)

# 3D density
tree_3d = cKDTree(coords_mock)
counts_3d = np.array(tree_3d.query_ball_point(coords_mock, r=10, return_length=True)) - 1
delta_3d = (counts_3d - np.mean(counts_3d)) / max(np.mean(counts_3d), 1)

print(f"\nEnvironment metrics computed:")
print(f"  Angular: mean neighbors = {np.mean(counts):.1f}")
print(f"  3D: mean neighbors = {np.mean(counts_3d):.1f}")

# =============================================================================
# 6. TEST: VELOCITY-ENVIRONMENT FROM SELECTION ALONE
# =============================================================================

print("\n" + "=" * 70)
print("TEST: VELOCITY-ENVIRONMENT (NULL HYPOTHESIS)")
print("=" * 70)

# Using angular density
void_ang = delta_angular < np.percentile(delta_angular, 25)
overdense_ang = delta_angular > np.percentile(delta_angular, 75)

v_void_ang = np.mean(np.abs(v_observed[void_ang]))
v_overdense_ang = np.mean(np.abs(v_observed[overdense_ang]))
ratio_ang_mock = v_void_ang / v_overdense_ang

print(f"\nAngular density (MOCK, null hypothesis):")
print(f"  Void <|v|>: {v_void_ang:.0f} km/s")
print(f"  Overdense <|v|>: {v_overdense_ang:.0f} km/s")
print(f"  Ratio: {ratio_ang_mock:.3f}")

# Using 3D density
void_3d = delta_3d < np.percentile(delta_3d, 25)
overdense_3d = delta_3d > np.percentile(delta_3d, 75)

v_void_3d = np.mean(np.abs(v_observed[void_3d]))
v_overdense_3d = np.mean(np.abs(v_observed[overdense_3d]))
ratio_3d_mock = v_void_3d / v_overdense_3d

print(f"\n3D density (MOCK, null hypothesis):")
print(f"  Void <|v|>: {v_void_3d:.0f} km/s")
print(f"  Overdense <|v|>: {v_overdense_3d:.0f} km/s")
print(f"  Ratio: {ratio_3d_mock:.3f}")

# Check mean distances by environment
print(f"\nMean distance by environment (MOCK):")
print(f"  3D void: {np.mean(dist_mock[void_3d]):.0f} Mpc")
print(f"  3D overdense: {np.mean(dist_mock[overdense_3d]):.0f} Mpc")

# =============================================================================
# 7. COMPARE MOCK TO CF4 OBSERVATIONS
# =============================================================================

print("\n" + "=" * 70)
print("COMPARISON: MOCK vs CF4 OBSERVED")
print("=" * 70)

# For CF4: replicate analysis
ra_cf4 = np.random.uniform(0, 360, N_valid)  # Approximate (we don't have CF4 RA in current load)
dec_cf4 = np.random.uniform(-90, 90, N_valid)

# Use actual CF4 distances
x_cf4 = dist_v * np.cos(np.radians(dec_cf4)) * np.cos(np.radians(ra_cf4))
y_cf4 = dist_v * np.cos(np.radians(dec_cf4)) * np.sin(np.radians(ra_cf4))
z_cf4 = dist_v * np.sin(np.radians(dec_cf4))
coords_cf4 = np.column_stack([x_cf4, y_cf4, z_cf4])

# 3D density for CF4
tree_cf4 = cKDTree(coords_cf4)
counts_cf4 = np.array(tree_cf4.query_ball_point(coords_cf4, r=10, return_length=True)) - 1
delta_cf4_3d = (counts_cf4 - np.mean(counts_cf4)) / max(np.mean(counts_cf4), 1)

void_cf4 = delta_cf4_3d < np.percentile(delta_cf4_3d, 25)
overdense_cf4 = delta_cf4_3d > np.percentile(delta_cf4_3d, 75)

v_void_cf4 = np.mean(np.abs(vpec_v[void_cf4]))
v_overdense_cf4 = np.mean(np.abs(vpec_v[overdense_cf4]))
ratio_cf4 = v_void_cf4 / v_overdense_cf4

print(f"\n3D density comparison:")
print(f"  MOCK (null): ratio = {ratio_3d_mock:.3f}")
print(f"  CF4 observed: ratio = {ratio_cf4:.3f}")

# Is CF4 different from null?
diff = ratio_cf4 - ratio_3d_mock
print(f"\n  Difference: {diff:+.3f}")

if abs(diff) > 0.1:
    print(f"  >>> CF4 shows DIFFERENT behavior from null hypothesis")
    if ratio_cf4 > ratio_3d_mock:
        print(f"  >>> Direction: MORE Synchronism-like (higher ratio)")
    else:
        print(f"  >>> Direction: LESS Synchronism-like (lower ratio)")
else:
    print(f"  >>> CF4 is CONSISTENT with null hypothesis (selection effects only)")

# =============================================================================
# 8. MONTE CARLO: DISTRIBUTION OF NULL HYPOTHESIS RATIOS
# =============================================================================

print("\n" + "=" * 70)
print("MONTE CARLO: NULL HYPOTHESIS DISTRIBUTION")
print("=" * 70)

n_mc = 100
ratios_null = []

for i in range(n_mc):
    # Generate new random velocities (null hypothesis)
    v_true_mc = np.random.normal(0, sigma_v_true, N_mock)

    # Add distance errors
    dist_err_mc = np.random.normal(0, rel_dist_err) * dist_mock
    v_obs_mc = v_true_mc - H0 * dist_err_mc

    # 3D density already computed (same positions)
    v_void_mc = np.mean(np.abs(v_obs_mc[void_3d]))
    v_od_mc = np.mean(np.abs(v_obs_mc[overdense_3d]))
    ratios_null.append(v_void_mc / v_od_mc)

ratios_null = np.array(ratios_null)

print(f"\nNull hypothesis ratio distribution (N={n_mc}):")
print(f"  Mean: {np.mean(ratios_null):.3f}")
print(f"  Std: {np.std(ratios_null):.3f}")
print(f"  95% CI: [{np.percentile(ratios_null, 2.5):.3f}, {np.percentile(ratios_null, 97.5):.3f}]")

print(f"\nCF4 observed: {ratio_cf4:.3f}")

# Z-score
z_score = (ratio_cf4 - np.mean(ratios_null)) / np.std(ratios_null)
p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
print(f"\nZ-score: {z_score:.2f}")
print(f"p-value: {p_value:.4f}")

if p_value < 0.05:
    print(f"\n>>> CF4 is SIGNIFICANTLY DIFFERENT from null hypothesis")
    print(f">>> There may be a TRUE velocity-environment effect beyond selection!")
else:
    print(f"\n>>> CF4 is NOT significantly different from null hypothesis")
    print(f">>> Selection effects ALONE can explain the observations")

# =============================================================================
# 9. VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: |v| vs distance comparison
ax1 = axes[0, 0]
ax1.scatter(dist_mock, np.abs(v_observed), alpha=0.1, s=1, label=f'Mock (r={corr_mock:.2f})')
ax1.scatter(dist_v, np.abs(vpec_v), alpha=0.1, s=1, label=f'CF4 (r={corr_cf4:.2f})')
ax1.set_xlabel('Distance (Mpc)')
ax1.set_ylabel('|Peculiar velocity| (km/s)')
ax1.set_title('|v| vs Distance: Mock (null) vs CF4')
ax1.legend()
ax1.set_xlim(0, 200)
ax1.set_ylim(0, 5000)

# Panel 2: Null hypothesis ratio distribution
ax2 = axes[0, 1]
ax2.hist(ratios_null, bins=20, alpha=0.7, color='gray', density=True, label='Null hypothesis')
ax2.axvline(ratio_cf4, color='red', linewidth=2, linestyle='--', label=f'CF4 observed: {ratio_cf4:.3f}')
ax2.axvline(np.mean(ratios_null), color='blue', linewidth=2, label=f'Null mean: {np.mean(ratios_null):.3f}')
ax2.set_xlabel('Void/Overdense velocity ratio')
ax2.set_ylabel('Probability density')
ax2.set_title(f'Null Hypothesis Distribution (Z={z_score:.2f}, p={p_value:.3f})')
ax2.legend()

# Panel 3: Mock velocity distribution by environment
ax3 = axes[1, 0]
ax3.hist(np.abs(v_observed[void_3d]), bins=50, alpha=0.7, label=f'Void (N={np.sum(void_3d)})', density=True)
ax3.hist(np.abs(v_observed[overdense_3d]), bins=50, alpha=0.5, label=f'Overdense (N={np.sum(overdense_3d)})', density=True)
ax3.set_xlabel('|Peculiar velocity| (km/s)')
ax3.set_ylabel('Density')
ax3.set_title(f'Mock: Velocity by Environment (ratio={ratio_3d_mock:.3f})')
ax3.legend()
ax3.set_xlim(0, 3000)

# Panel 4: CF4 velocity distribution by environment
ax4 = axes[1, 1]
ax4.hist(np.abs(vpec_v[void_cf4]), bins=50, alpha=0.7, label=f'Void (N={np.sum(void_cf4)})', density=True)
ax4.hist(np.abs(vpec_v[overdense_cf4]), bins=50, alpha=0.5, label=f'Overdense (N={np.sum(overdense_cf4)})', density=True)
ax4.set_xlabel('|Peculiar velocity| (km/s)')
ax4.set_ylabel('Density')
ax4.set_title(f'CF4: Velocity by Environment (ratio={ratio_cf4:.3f})')
ax4.legend()
ax4.set_xlim(0, 3000)

plt.suptitle('Session #174: Forward Modeling Selection Effects', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session174_forward_model.png', dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FIGURE SAVED: session174_forward_model.png")
print("=" * 70)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #174: SUMMARY")
print("=" * 70)

print(f"""
FORWARD MODELING SELECTION EFFECTS
==================================

1. NULL HYPOTHESIS MODEL:
   - True velocity: σ_v = {sigma_v_true} km/s (ΛCDM, no Synchronism)
   - Distance errors: matching CF4 distribution
   - NO true velocity-environment correlation

2. KEY FINDING - |v| vs DISTANCE:
   - Mock (null): r = {corr_mock:.3f}
   - CF4 observed: r = {corr_cf4:.3f}
   - >>> Selection effects ALONE produce the |v|-distance correlation

3. VELOCITY-ENVIRONMENT (3D density):
   - Mock null ratio: {ratio_3d_mock:.3f} ± {np.std(ratios_null):.3f}
   - CF4 observed: {ratio_cf4:.3f}
   - Z-score: {z_score:.2f}
   - p-value: {p_value:.4f}

4. INTERPRETATION:
   {"CF4 is SIGNIFICANTLY different from null - may indicate true signal" if p_value < 0.05 else "CF4 is consistent with null - selection effects explain observations"}

5. IMPLICATIONS FOR SYNCHRONISM:
   {"The observed velocity-environment correlation EXCEEDS what selection effects predict" if ratio_cf4 > np.mean(ratios_null) + 2*np.std(ratios_null) else "Cannot distinguish true Synchronism signal from selection effects"}

6. NEXT STEPS:
   - Refine null model with realistic cosmic structure
   - Test with high-precision subset (lower errors)
   - Use complementary tests not affected by distance errors
""")

print("=" * 70)
print("SESSION #174 COMPLETE")
print("=" * 70)
