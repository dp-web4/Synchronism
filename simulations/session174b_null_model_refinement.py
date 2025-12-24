#!/usr/bin/env python3
"""
SESSION #174b: INTERPRETING THE NULL MODEL RESULT
==================================================
Date: December 24, 2025

SURPRISING FINDING FROM 174:
----------------------------
Null hypothesis (pure selection effects) predicts: ratio = 2.82
CF4 observed: ratio = 2.15

CF4 shows LESS velocity-environment correlation than selection effects alone!

This means:
1. True velocities may be ANTI-correlated with environment
2. OR our error model is too extreme
3. OR there's a physical suppression effect

This session investigates the null model and refines it.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #174b: NULL MODEL REFINEMENT")
print("=" * 70)

# =============================================================================
# 1. LOAD CF4 DATA
# =============================================================================

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table4.dat', 'r') as f:
    lines = f.readlines()

N_cf4 = len(lines)
dist_cf4 = np.zeros(N_cf4)
dm_err_cf4 = np.zeros(N_cf4)
vpec_cf4 = np.zeros(N_cf4)
ra_cf4 = np.zeros(N_cf4)
dec_cf4 = np.zeros(N_cf4)

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
    ra_cf4[i] = safe_float(line[83:91])
    dec_cf4[i] = safe_float(line[92:100])

H0 = 74.6
valid = (dist_cf4 > 10) & (dist_cf4 < 200) & (dm_err_cf4 > 0) & (np.abs(vpec_cf4) < 5000)
N_valid = np.sum(valid)

dist_v = dist_cf4[valid]
dm_err_v = dm_err_cf4[valid]
vpec_v = vpec_cf4[valid]
ra_v = ra_cf4[valid]
dec_v = dec_cf4[valid]

print(f"\nLoaded {N_valid} valid CF4 groups")

# =============================================================================
# 2. USE ACTUAL CF4 POSITIONS (not random)
# =============================================================================

print("\n" + "=" * 70)
print("USING ACTUAL CF4 POSITIONS")
print("=" * 70)

# Convert to Cartesian using actual RA/Dec
ra_rad = np.radians(ra_v)
dec_rad = np.radians(dec_v)
x_cf4 = dist_v * np.cos(dec_rad) * np.cos(ra_rad)
y_cf4 = dist_v * np.cos(dec_rad) * np.sin(ra_rad)
z_cf4 = dist_v * np.sin(dec_rad)
coords_cf4 = np.column_stack([x_cf4, y_cf4, z_cf4])

# 3D density with actual positions
tree_3d = cKDTree(coords_cf4)
counts_3d = np.array(tree_3d.query_ball_point(coords_cf4, r=10, return_length=True)) - 1
delta_3d = (counts_3d - np.mean(counts_3d)) / max(np.mean(counts_3d), 1)

void_3d = delta_3d < np.percentile(delta_3d, 25)
overdense_3d = delta_3d > np.percentile(delta_3d, 75)

print(f"\n3D density with actual positions:")
print(f"  Void (N={np.sum(void_3d)}): mean d = {np.mean(dist_v[void_3d]):.0f} Mpc")
print(f"  Overdense (N={np.sum(overdense_3d)}): mean d = {np.mean(dist_v[overdense_3d]):.0f} Mpc")

# CF4 actual ratio
v_void_cf4 = np.mean(np.abs(vpec_v[void_3d]))
v_od_cf4 = np.mean(np.abs(vpec_v[overdense_3d]))
ratio_cf4 = v_void_cf4 / v_od_cf4

print(f"\nCF4 observed:")
print(f"  Void <|v|>: {v_void_cf4:.0f} km/s")
print(f"  Overdense <|v|>: {v_od_cf4:.0f} km/s")
print(f"  Ratio: {ratio_cf4:.3f}")

# =============================================================================
# 3. ANALYZE THE DISTANCE BIAS IN 3D DENSITY
# =============================================================================

print("\n" + "=" * 70)
print("DISTANCE BIAS IN 3D DENSITY CLASSIFICATION")
print("=" * 70)

# The issue: 3D density is computed in observed (erroneous) distance space
# "Voids" in 3D tend to be at large distances (fewer neighbors due to sparse sampling)
# "Overdense" regions are nearby (more neighbors)

dist_ratio = np.mean(dist_v[void_3d]) / np.mean(dist_v[overdense_3d])
print(f"\nDistance ratio (void/overdense): {dist_ratio:.2f}")

# This creates spurious |v|-environment correlation because:
# |v_obs| ∝ distance (due to errors)
# 3D void ↔ high distance ↔ high |v_obs|

# The NULL model should use SAME distance-environment relationship as CF4
# Not random positions!

# =============================================================================
# 4. REFINED NULL MODEL: SAME POSITIONS, RANDOM VELOCITIES
# =============================================================================

print("\n" + "=" * 70)
print("REFINED NULL MODEL")
print("=" * 70)

np.random.seed(42)

# Use ACTUAL CF4 positions (same 3D density classification)
# But replace velocities with random + errors

sigma_v_true = 300  # km/s

# Velocity error model
v_err = H0 * dist_v * 0.461 * dm_err_v

n_mc = 100
ratios_null_refined = []

for i in range(n_mc):
    # True velocity: random with no environment dependence
    v_true = np.random.normal(0, sigma_v_true, N_valid)

    # Add distance-dependent errors
    rel_dist_err = 0.461 * dm_err_v
    dist_error = np.random.normal(0, rel_dist_err) * dist_v

    # Observed velocity
    v_obs = v_true - H0 * dist_error

    # Compute ratio using SAME environment classification
    v_void = np.mean(np.abs(v_obs[void_3d]))
    v_od = np.mean(np.abs(v_obs[overdense_3d]))
    ratios_null_refined.append(v_void / v_od)

ratios_null_refined = np.array(ratios_null_refined)

print(f"\nRefined null hypothesis (same positions as CF4):")
print(f"  Mean ratio: {np.mean(ratios_null_refined):.3f}")
print(f"  Std: {np.std(ratios_null_refined):.3f}")
print(f"  95% CI: [{np.percentile(ratios_null_refined, 2.5):.3f}, {np.percentile(ratios_null_refined, 97.5):.3f}]")

print(f"\nCF4 observed: {ratio_cf4:.3f}")

z_score = (ratio_cf4 - np.mean(ratios_null_refined)) / np.std(ratios_null_refined)
p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

print(f"\nZ-score: {z_score:.2f}")
print(f"p-value: {p_value:.4f}")

# =============================================================================
# 5. INTERPRETATION
# =============================================================================

print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)

if ratio_cf4 < np.mean(ratios_null_refined):
    print(f"\nCF4 shows LESS velocity-environment correlation than null prediction")
    print(f"\nPossible explanations:")
    print(f"  1. True peculiar velocities are ANTI-correlated with 3D density")
    print(f"  2. Dense regions have higher true |v| than voids (opposite to Synchronism)")
    print(f"  3. Malmquist bias or other selection effects we haven't modeled")
    print(f"  4. CF4 v_pec includes corrections that reduce distance-error effects")
else:
    print(f"\nCF4 shows MORE velocity-environment correlation than null prediction")
    print(f"  >>> This would suggest a TRUE Synchronism-like signal!")

# =============================================================================
# 6. TEST DIFFERENT TRUE VELOCITY SCENARIOS
# =============================================================================

print("\n" + "=" * 70)
print("TESTING TRUE VELOCITY SCENARIOS")
print("=" * 70)

# Scenario 1: Synchronism - enhanced velocities in voids
# Scenario 2: Anti-Synchronism - suppressed velocities in voids
# Scenario 3: Null - no environment dependence

scenarios = [
    ("Null (σ=300, no env)", 300, 1.0),
    ("Synchronism (+30% in voids)", 300, 1.3),
    ("Anti-Synchronism (-30% in voids)", 300, 0.7),
    ("Higher true σ (500 km/s)", 500, 1.0),
]

print("\nScenario testing:")
print("-" * 60)

for name, sigma, void_enhancement in scenarios:
    ratios_scenario = []

    for _ in range(50):
        # Generate true velocities
        v_true = np.random.normal(0, sigma, N_valid)

        # Apply environment dependence
        if void_enhancement != 1.0:
            enhancement = np.ones(N_valid)
            enhancement[void_3d] = void_enhancement
            v_true = v_true * enhancement

        # Add errors
        rel_dist_err = 0.461 * dm_err_v
        dist_error = np.random.normal(0, rel_dist_err) * dist_v
        v_obs = v_true - H0 * dist_error

        # Compute ratio
        v_void = np.mean(np.abs(v_obs[void_3d]))
        v_od = np.mean(np.abs(v_obs[overdense_3d]))
        ratios_scenario.append(v_void / v_od)

    mean_ratio = np.mean(ratios_scenario)
    std_ratio = np.std(ratios_scenario)

    match = abs(mean_ratio - ratio_cf4) < 2 * std_ratio
    print(f"  {name}: ratio = {mean_ratio:.3f} ± {std_ratio:.3f}", end="")
    if match:
        print(" *** MATCHES CF4 ***")
    else:
        print()

# =============================================================================
# 7. WHAT ENHANCEMENT MATCHES CF4?
# =============================================================================

print("\n" + "=" * 70)
print("FINDING BEST-FIT ENHANCEMENT")
print("=" * 70)

# Binary search for the void enhancement that matches CF4
enhancements = np.linspace(0.5, 1.5, 21)
mean_ratios = []

for enh in enhancements:
    ratios_test = []
    for _ in range(30):
        v_true = np.random.normal(0, 300, N_valid)
        enhancement = np.ones(N_valid)
        enhancement[void_3d] = enh
        v_true = v_true * enhancement

        rel_dist_err = 0.461 * dm_err_v
        dist_error = np.random.normal(0, rel_dist_err) * dist_v
        v_obs = v_true - H0 * dist_error

        v_void = np.mean(np.abs(v_obs[void_3d]))
        v_od = np.mean(np.abs(v_obs[overdense_3d]))
        ratios_test.append(v_void / v_od)

    mean_ratios.append(np.mean(ratios_test))

mean_ratios = np.array(mean_ratios)

# Find enhancement that best matches CF4
best_idx = np.argmin(np.abs(mean_ratios - ratio_cf4))
best_enhancement = enhancements[best_idx]

print(f"\nBest-fit void enhancement: {best_enhancement:.2f}")
print(f"  Null prediction (enh=1.0): ratio = {mean_ratios[10]:.3f}")
print(f"  Best-fit ratio: {mean_ratios[best_idx]:.3f}")
print(f"  CF4 observed: {ratio_cf4:.3f}")

if best_enhancement < 1.0:
    reduction = (1.0 - best_enhancement) * 100
    print(f"\n>>> CF4 is consistent with {reduction:.0f}% LOWER velocities in voids")
    print(f">>> This is OPPOSITE to Synchronism prediction!")
elif best_enhancement > 1.0:
    increase = (best_enhancement - 1.0) * 100
    print(f"\n>>> CF4 is consistent with {increase:.0f}% HIGHER velocities in voids")
    print(f">>> This SUPPORTS Synchronism prediction!")
else:
    print(f"\n>>> CF4 is consistent with null hypothesis (no enhancement)")

# =============================================================================
# 8. VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Enhancement vs predicted ratio
ax1 = axes[0, 0]
ax1.plot(enhancements, mean_ratios, 'b-o', markersize=8)
ax1.axhline(ratio_cf4, color='red', linestyle='--', linewidth=2, label=f'CF4 observed: {ratio_cf4:.3f}')
ax1.axvline(1.0, color='gray', linestyle=':', label='Null (no enhancement)')
ax1.axvline(best_enhancement, color='green', linestyle='--', label=f'Best fit: {best_enhancement:.2f}')
ax1.set_xlabel('Void velocity enhancement factor')
ax1.set_ylabel('Expected void/overdense ratio')
ax1.set_title('Finding True Void Enhancement')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel 2: Null distribution comparison
ax2 = axes[0, 1]
ax2.hist(ratios_null_refined, bins=20, alpha=0.7, color='gray', density=True, label='Null hypothesis')
ax2.axvline(ratio_cf4, color='red', linewidth=2, linestyle='--', label=f'CF4: {ratio_cf4:.3f}')
ax2.axvline(np.mean(ratios_null_refined), color='blue', linewidth=2, label=f'Null mean: {np.mean(ratios_null_refined):.3f}')
ax2.set_xlabel('Void/Overdense ratio')
ax2.set_ylabel('Density')
ax2.set_title(f'Null vs Observed (Z={z_score:.2f})')
ax2.legend()

# Panel 3: Distance distribution by environment
ax3 = axes[1, 0]
ax3.hist(dist_v[void_3d], bins=30, alpha=0.7, label=f'3D void (mean={np.mean(dist_v[void_3d]):.0f} Mpc)', density=True)
ax3.hist(dist_v[overdense_3d], bins=30, alpha=0.5, label=f'3D overdense (mean={np.mean(dist_v[overdense_3d]):.0f} Mpc)', density=True)
ax3.set_xlabel('Distance (Mpc)')
ax3.set_ylabel('Density')
ax3.set_title('Distance Distribution by 3D Environment')
ax3.legend()

# Panel 4: Velocity by environment
ax4 = axes[1, 1]
ax4.hist(np.abs(vpec_v[void_3d]), bins=50, alpha=0.7, label=f'Void: <|v|>={v_void_cf4:.0f}', density=True)
ax4.hist(np.abs(vpec_v[overdense_3d]), bins=50, alpha=0.5, label=f'Overdense: <|v|>={v_od_cf4:.0f}', density=True)
ax4.set_xlabel('|Peculiar velocity| (km/s)')
ax4.set_ylabel('Density')
ax4.set_title(f'CF4 Velocity Distribution (ratio={ratio_cf4:.3f})')
ax4.legend()
ax4.set_xlim(0, 3000)

plt.suptitle('Session #174b: Null Model Refinement', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session174b_null_refinement.png', dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FIGURE SAVED: session174b_null_refinement.png")
print("=" * 70)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #174b: SUMMARY")
print("=" * 70)

print(f"""
NULL MODEL REFINEMENT
=====================

1. REFINED NULL HYPOTHESIS:
   - Uses ACTUAL CF4 positions (same 3D density classification)
   - Random velocities with σ = 300 km/s
   - Same distance error distribution as CF4

2. KEY RESULT:
   - Null prediction: ratio = {np.mean(ratios_null_refined):.3f} ± {np.std(ratios_null_refined):.3f}
   - CF4 observed: ratio = {ratio_cf4:.3f}
   - Z-score: {z_score:.2f}
   - Direction: CF4 shows LESS correlation than null

3. BEST-FIT TRUE ENHANCEMENT:
   - Best match: void enhancement = {best_enhancement:.2f}
   - This means: {100*abs(1-best_enhancement):.0f}% {"LOWER" if best_enhancement < 1 else "HIGHER"}
     true velocities in voids

4. INTERPRETATION:
   The 3D density test shows that AFTER accounting for selection effects,
   CF4 galaxies in 3D "void" regions have {"LOWER" if best_enhancement < 1 else "HIGHER"}
   true peculiar velocities than those in overdense regions.

   This is {"OPPOSITE to" if best_enhancement < 1 else "CONSISTENT with"} Synchronism predictions.

5. CAVEAT:
   The 3D "void" classification is based on OBSERVED distances, not true
   structure. The true cosmic voids may have different velocity properties.

6. NEXT STEPS:
   - Test with angular density (less distance-dependent)
   - Use true void catalogs
   - Consider that CF4 v_pec may already be bias-corrected
""")

print("=" * 70)
print("SESSION #174b COMPLETE")
print("=" * 70)
