#!/usr/bin/env python3
"""
SESSION #171: REFINED ENVIRONMENT ANALYSIS
===========================================
Date: December 23, 2025
Focus: Improve environment classification and quantify Synchronism signature

CONTEXT FROM SESSION #170:
- Velocity-quartile test: >> 10σ (high-|v| galaxies in lower density)
- But angular overdensity metric had issues
- Need better environment proxies

THIS SESSION:
1. Use redshift-space density field
2. Identify void-like and cluster-like environments
3. Quantify velocity enhancement properly
4. Compare to Synchronism predictions
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #171: REFINED ENVIRONMENT ANALYSIS")
print("=" * 70)
print("Improving density estimation and Synchronism signature quantification")
print("=" * 70)

# =============================================================================
# LOAD CF4 DATA
# =============================================================================

print("\n" + "=" * 70)
print("PART 1: LOADING CF4 DATA")
print("=" * 70)

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table2.dat', 'r') as f:
    lines = f.readlines()

N_total = len(lines)
print(f"Loading {N_total} galaxies...")

vcmb = np.zeros(N_total)
dm = np.zeros(N_total)
dm_err = np.zeros(N_total)
ra = np.zeros(N_total)
dec = np.zeros(N_total)

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
frac_dist_err = 0.4605 * dm_err
v_pec = vcmb - H0 * distances
v_err = H0 * distances * frac_dist_err

# Quality cuts
valid = (dm > 20) & (dm_err > 0) & (vcmb > 500) & (distances > 10) & (distances < 200)
valid &= (v_err < 3000)  # Reasonable errors

N_valid = np.sum(valid)
print(f"Valid galaxies: {N_valid}")

# Extract valid data
d_valid = distances[valid]
ra_valid = ra[valid]
dec_valid = dec[valid]
v_pec_valid = v_pec[valid]
v_err_valid = v_err[valid]
vcmb_valid = vcmb[valid]

# =============================================================================
# PART 2: REDSHIFT-SPACE DENSITY FIELD
# =============================================================================

print("\n" + "=" * 70)
print("PART 2: REDSHIFT-SPACE DENSITY ESTIMATION")
print("=" * 70)

print("""
Strategy: Use redshift (vcmb) as distance proxy for density estimation.
This avoids distance error issues while still providing 3D structure.

Steps:
1. Convert (RA, Dec, v_cmb) to pseudo-3D coordinates
2. Compute local density via fixed comoving volume
3. Compare to mean density at that redshift
""")

# Convert to pseudo-Cartesian using vcmb as distance proxy
# (In km/s units, so effectively r = vcmb/100 gives ~Mpc)
r_redshift = vcmb_valid / 100  # Pseudo-distance in 100 km/s units

ra_rad = np.radians(ra_valid)
dec_rad = np.radians(dec_valid)

x_z = r_redshift * np.cos(dec_rad) * np.cos(ra_rad)
y_z = r_redshift * np.cos(dec_rad) * np.sin(ra_rad)
z_z = r_redshift * np.sin(dec_rad)

coords_z = np.column_stack([x_z, y_z, z_z])

# Build k-d tree in redshift space
tree_z = cKDTree(coords_z)

# Compute density with fixed physical radius
# Use radius ~ 10 Mpc equivalent (in our units: ~10)
radius_fixed = 10.0  # In 100 km/s units ~ 10 Mpc
print(f"Computing density with radius = {radius_fixed} (~ 10 Mpc)...")

counts = tree_z.query_ball_point(coords_z, r=radius_fixed, return_length=True)
counts = np.array(counts) - 1  # Subtract self

# Density normalization: compare to expected at each redshift shell
# Define redshift shells
v_shells = np.linspace(500, 20000, 20)
v_centers = 0.5 * (v_shells[:-1] + v_shells[1:])
mean_counts = np.zeros(len(v_centers))

for i in range(len(v_centers)):
    mask = (vcmb_valid >= v_shells[i]) & (vcmb_valid < v_shells[i+1])
    if np.sum(mask) > 100:
        mean_counts[i] = np.mean(counts[mask])

# Interpolate expected count for each galaxy
expected_counts = np.interp(vcmb_valid, v_centers, mean_counts)
expected_counts = np.maximum(expected_counts, 1)

# Overdensity
delta_z = (counts - expected_counts) / expected_counts

print(f"\nRedshift-space overdensity:")
print(f"  Min δ: {np.min(delta_z):.2f}")
print(f"  Max δ: {np.max(delta_z):.2f}")
print(f"  Median δ: {np.median(delta_z):.2f}")
print(f"  Mean δ: {np.mean(delta_z):.2f}")

# Environment classification
n_void = np.sum(delta_z < -0.5)
n_underdense = np.sum((delta_z >= -0.5) & (delta_z < -0.2))
n_mean = np.sum((delta_z >= -0.2) & (delta_z < 0.5))
n_overdense = np.sum((delta_z >= 0.5) & (delta_z < 2))
n_cluster = np.sum(delta_z >= 2)

print(f"\nEnvironment breakdown:")
print(f"  Void (δ < -0.5): {n_void} ({100*n_void/N_valid:.1f}%)")
print(f"  Underdense: {n_underdense} ({100*n_underdense/N_valid:.1f}%)")
print(f"  Mean: {n_mean} ({100*n_mean/N_valid:.1f}%)")
print(f"  Overdense: {n_overdense} ({100*n_overdense/N_valid:.1f}%)")
print(f"  Cluster: {n_cluster} ({100*n_cluster/N_valid:.1f}%)")

# =============================================================================
# PART 3: SYNCHRONISM PREDICTIONS
# =============================================================================

print("\n" + "=" * 70)
print("PART 3: SYNCHRONISM PREDICTIONS")
print("=" * 70)

Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2

def C_coherence(rho_ratio):
    """Coherence function"""
    x = np.maximum(rho_ratio, 0.01)
    x_phi = x ** (1/phi)
    return Omega_m + (1 - Omega_m) * x_phi / (1 + x_phi)

def velocity_enhancement(delta):
    """v_sync / v_LCDM"""
    rho_ratio = np.maximum(1 + delta, 0.01)
    return 0.97 * np.sqrt(1.0 / C_coherence(rho_ratio))

print("\nSynchronism predictions:")
for delta in [-0.8, -0.5, -0.3, 0.0, 0.5, 1.0, 3.0]:
    v_enh = velocity_enhancement(delta)
    print(f"  δ = {delta:+5.1f}: v_enhancement = {v_enh:.3f} ({100*(v_enh-1):+.1f}%)")

# =============================================================================
# PART 4: VELOCITY ANALYSIS BY ENVIRONMENT
# =============================================================================

print("\n" + "=" * 70)
print("PART 4: ERROR-WEIGHTED VELOCITY ANALYSIS")
print("=" * 70)

def weighted_mean(v, e):
    w = 1.0 / (e**2 + 1e-10)
    return np.sum(v * w) / np.sum(w)

def weighted_err(e):
    w = 1.0 / (e**2 + 1e-10)
    return 1.0 / np.sqrt(np.sum(w))

# Environment bins
env_bins = [
    ('Deep Void', -1.0, -0.5),
    ('Void', -0.5, -0.3),
    ('Underdense', -0.3, -0.1),
    ('Mean', -0.1, 0.3),
    ('Overdense', 0.3, 1.0),
    ('Cluster', 1.0, 10.0)
]

print("\nEnvironment-binned velocity analysis:")
print("-" * 80)
print(f"{'Environment':12s} {'N':>6s} {'<δ>':>8s} {'<|v|>':>12s} {'Error':>8s} {'Pred v_enh':>10s}")
print("-" * 80)

results = []
for env_name, delta_min, delta_max in env_bins:
    mask = (delta_z >= delta_min) & (delta_z < delta_max)
    n = np.sum(mask)
    if n < 100:
        continue

    mean_delta = np.mean(delta_z[mask])
    mean_v = weighted_mean(np.abs(v_pec_valid[mask]), v_err_valid[mask])
    mean_err = weighted_err(v_err_valid[mask])
    pred_enh = velocity_enhancement(mean_delta)

    print(f"{env_name:12s} {n:6d} {mean_delta:+8.3f} {mean_v:12.1f} {mean_err:8.1f} {pred_enh:10.3f}")

    results.append({
        'env': env_name,
        'n': n,
        'delta': mean_delta,
        'v': mean_v,
        'err': mean_err,
        'pred': pred_enh
    })

# =============================================================================
# PART 5: SYNCHRONISM SIGNATURE - VOID vs OVERDENSE RATIO
# =============================================================================

print("\n" + "=" * 70)
print("PART 5: SYNCHRONISM SIGNATURE TEST")
print("=" * 70)

void_mask = delta_z < -0.3
overdense_mask = delta_z > 0.5

n_void_test = np.sum(void_mask)
n_od_test = np.sum(overdense_mask)

print(f"\nVoid sample (δ < -0.3): N = {n_void_test}")
print(f"Overdense sample (δ > 0.5): N = {n_od_test}")

if n_void_test > 200 and n_od_test > 200:
    # Weighted velocities
    v_void = weighted_mean(np.abs(v_pec_valid[void_mask]), v_err_valid[void_mask])
    v_void_err = weighted_err(v_err_valid[void_mask])

    v_od = weighted_mean(np.abs(v_pec_valid[overdense_mask]), v_err_valid[overdense_mask])
    v_od_err = weighted_err(v_err_valid[overdense_mask])

    # Ratio
    ratio = v_void / v_od
    ratio_err = ratio * np.sqrt((v_void_err/v_void)**2 + (v_od_err/v_od)**2)

    # Predictions
    mean_delta_void = np.mean(delta_z[void_mask])
    mean_delta_od = np.mean(delta_z[overdense_mask])

    pred_void = velocity_enhancement(mean_delta_void)
    pred_od = velocity_enhancement(mean_delta_od)
    pred_ratio = pred_void / pred_od

    # Significance
    signif_lcdm = (ratio - 1.0) / ratio_err
    signif_sync = (ratio - pred_ratio) / ratio_err

    print(f"\nRESULTS:")
    print("-" * 50)
    print(f"Void <|v|>: {v_void:.1f} ± {v_void_err:.1f} km/s (mean δ = {mean_delta_void:.3f})")
    print(f"Overdense <|v|>: {v_od:.1f} ± {v_od_err:.1f} km/s (mean δ = {mean_delta_od:.3f})")
    print()
    print(f"Observed ratio: {ratio:.4f} ± {ratio_err:.4f}")
    print(f"Expected (Synchronism): {pred_ratio:.4f}")
    print(f"Expected (ΛCDM): 1.0000")
    print()
    print(f"Deviation from ΛCDM: {signif_lcdm:.2f}σ")
    print(f"Deviation from Synchronism: {signif_sync:.2f}σ")

    if signif_lcdm > 3:
        print("\n>>> SIGNIFICANT DEVIATION FROM ΛCDM <<<")
        if abs(signif_sync) < 2:
            print(">>> CONSISTENT WITH SYNCHRONISM <<<")
        else:
            print(">>> MAGNITUDE DIFFERS FROM PREDICTION <<<")

# =============================================================================
# PART 6: VELOCITY-ENVIRONMENT CORRELATION
# =============================================================================

print("\n" + "=" * 70)
print("PART 6: VELOCITY-ENVIRONMENT CORRELATION")
print("=" * 70)

# Bin by environment
delta_bins = np.linspace(-0.8, 3, 20)
delta_centers = 0.5 * (delta_bins[:-1] + delta_bins[1:])

v_binned = np.zeros(len(delta_centers))
v_binned_err = np.zeros(len(delta_centers))
n_binned = np.zeros(len(delta_centers))

for i in range(len(delta_centers)):
    mask = (delta_z >= delta_bins[i]) & (delta_z < delta_bins[i+1])
    n = np.sum(mask)
    n_binned[i] = n
    if n >= 50:
        v_binned[i] = weighted_mean(np.abs(v_pec_valid[mask]), v_err_valid[mask])
        v_binned_err[i] = weighted_err(v_err_valid[mask])

# Normalize to mean
valid_bins = n_binned >= 50
mean_v_norm = np.mean(v_binned[valid_bins])
v_norm = v_binned / mean_v_norm
v_norm_err = v_binned_err / mean_v_norm

# Synchronism prediction curve
delta_theory = np.linspace(-0.8, 3, 100)
v_theory = velocity_enhancement(delta_theory)
v_theory_norm = v_theory / velocity_enhancement(0)

# Fit quality
if np.sum(valid_bins) >= 5:
    # Chi-squared for Synchronism vs data
    theory_at_bins = velocity_enhancement(delta_centers[valid_bins])
    theory_at_bins_norm = theory_at_bins / velocity_enhancement(0)

    residuals = v_norm[valid_bins] - theory_at_bins_norm
    chi2_sync = np.sum((residuals / v_norm_err[valid_bins])**2)

    # Chi-squared for ΛCDM (flat)
    chi2_lcdm = np.sum(((v_norm[valid_bins] - 1.0) / v_norm_err[valid_bins])**2)

    dof = np.sum(valid_bins) - 1
    delta_chi2 = chi2_lcdm - chi2_sync

    print(f"\nModel comparison ({np.sum(valid_bins)} bins):")
    print(f"  χ² (Synchronism): {chi2_sync:.1f}")
    print(f"  χ² (ΛCDM flat): {chi2_lcdm:.1f}")
    print(f"  Δχ²: {delta_chi2:.1f}")
    print(f"  Preference for Synchronism: {np.sqrt(max(delta_chi2, 0)):.1f}σ")

# =============================================================================
# PART 7: BOOTSTRAP CONFIDENCE INTERVALS
# =============================================================================

print("\n" + "=" * 70)
print("PART 7: BOOTSTRAP UNCERTAINTY ESTIMATION")
print("=" * 70)

n_bootstrap = 1000
boot_ratios = np.zeros(n_bootstrap)

void_indices = np.where(void_mask)[0]
od_indices = np.where(overdense_mask)[0]

print(f"Running {n_bootstrap} bootstrap iterations...")

for b in range(n_bootstrap):
    # Resample with replacement
    boot_void = np.random.choice(void_indices, len(void_indices), replace=True)
    boot_od = np.random.choice(od_indices, len(od_indices), replace=True)

    v_void_b = weighted_mean(np.abs(v_pec_valid[boot_void]), v_err_valid[boot_void])
    v_od_b = weighted_mean(np.abs(v_pec_valid[boot_od]), v_err_valid[boot_od])

    boot_ratios[b] = v_void_b / v_od_b

# Confidence intervals
boot_mean = np.mean(boot_ratios)
boot_std = np.std(boot_ratios)
boot_16 = np.percentile(boot_ratios, 16)
boot_84 = np.percentile(boot_ratios, 84)

print(f"\nBootstrap results (N={n_bootstrap}):")
print(f"  Mean ratio: {boot_mean:.4f}")
print(f"  Std: {boot_std:.4f}")
print(f"  68% CI: [{boot_16:.4f}, {boot_84:.4f}]")
print(f"  Deviation from 1.0: {(boot_mean - 1)/boot_std:.2f}σ")

# =============================================================================
# PART 8: VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("PART 8: GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Environment distribution
ax1 = axes[0, 0]
ax1.hist(delta_z, bins=60, color='forestgreen', edgecolor='black', alpha=0.7)
ax1.set_xlabel('Redshift-space Overdensity δ', fontsize=12)
ax1.set_ylabel('Number of Galaxies', fontsize=12)
ax1.set_title('Environment Distribution (CF4)', fontsize=14)
ax1.axvline(0, color='gray', linestyle='--', alpha=0.5, label='Mean density')
ax1.axvline(-0.3, color='blue', linestyle=':', label='Void threshold')
ax1.axvline(0.5, color='red', linestyle=':', label='Overdense threshold')
ax1.legend()
ax1.set_xlim(-1, 4)

# Panel 2: Velocity vs environment
ax2 = axes[0, 1]
ax2.errorbar(delta_centers[valid_bins], v_norm[valid_bins],
             yerr=v_norm_err[valid_bins], fmt='o', color='navy',
             capsize=3, label='CF4 data (weighted)', markersize=6)
ax2.plot(delta_theory, v_theory_norm, 'r-', linewidth=2, label='Synchronism')
ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='ΛCDM')
ax2.set_xlabel('Overdensity δ', fontsize=12)
ax2.set_ylabel('Normalized |v| (weighted)', fontsize=12)
ax2.set_title('Velocity Enhancement vs Environment', fontsize=14)
ax2.legend()
ax2.set_xlim(-1, 3)
ax2.set_ylim(0.7, 1.5)

# Panel 3: Bootstrap distribution
ax3 = axes[1, 0]
ax3.hist(boot_ratios, bins=40, color='steelblue', edgecolor='black', alpha=0.7)
ax3.axvline(boot_mean, color='red', linestyle='-', linewidth=2, label=f'Mean = {boot_mean:.3f}')
ax3.axvline(1.0, color='gray', linestyle='--', linewidth=2, label='ΛCDM = 1.0')
ax3.axvline(pred_ratio, color='green', linestyle=':', linewidth=2, label=f'Synchronism = {pred_ratio:.3f}')
ax3.set_xlabel('Void/Overdense Velocity Ratio', fontsize=12)
ax3.set_ylabel('Bootstrap Count', fontsize=12)
ax3.set_title('Bootstrap Distribution of Ratio', fontsize=14)
ax3.legend()

# Panel 4: Summary
ax4 = axes[1, 1]
ax4.text(0.5, 0.95, 'SESSION #171 RESULTS', fontsize=16, fontweight='bold',
         ha='center', transform=ax4.transAxes)
ax4.text(0.1, 0.80, f'Void galaxies (δ < -0.3): N = {n_void_test}', fontsize=11, transform=ax4.transAxes)
ax4.text(0.1, 0.72, f'Overdense galaxies (δ > 0.5): N = {n_od_test}', fontsize=11, transform=ax4.transAxes)
ax4.text(0.1, 0.60, f'Observed ratio: {ratio:.4f} ± {ratio_err:.4f}', fontsize=12, transform=ax4.transAxes)
ax4.text(0.1, 0.52, f'Synchronism prediction: {pred_ratio:.4f}', fontsize=12, transform=ax4.transAxes)
ax4.text(0.1, 0.44, f'ΛCDM prediction: 1.0000', fontsize=12, transform=ax4.transAxes)
ax4.text(0.1, 0.30, f'Deviation from ΛCDM: {signif_lcdm:.1f}σ',
         fontsize=14, fontweight='bold', color='red' if signif_lcdm > 3 else 'black',
         transform=ax4.transAxes)
ax4.text(0.1, 0.20, f'Bootstrap mean: {boot_mean:.4f} ± {boot_std:.4f}',
         fontsize=11, transform=ax4.transAxes)

status = "POTENTIAL DETECTION" if signif_lcdm > 3 else "INCONCLUSIVE"
ax4.text(0.1, 0.05, f'Status: {status}', fontsize=14, fontweight='bold',
         color='green' if signif_lcdm > 3 else 'gray', transform=ax4.transAxes)
ax4.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session171_refined_analysis.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session171_refined_analysis.png")

# =============================================================================
# PART 9: SESSION SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #171 SUMMARY")
print("=" * 70)

summary = f"""
SESSION #171: REFINED ENVIRONMENT ANALYSIS
==========================================

METHODOLOGY IMPROVEMENTS:
-------------------------
• Redshift-space density (avoids distance errors)
• Shell-normalized overdensity (fair comparison)
• Error-weighted statistics throughout
• Bootstrap uncertainty quantification

DATA:
-----
• CF4 galaxies: {N_valid} (after quality cuts)
• Void sample (δ < -0.3): {n_void_test}
• Overdense sample (δ > 0.5): {n_od_test}

KEY RESULTS:
------------
• Void <|v|>: {v_void:.1f} ± {v_void_err:.1f} km/s
• Overdense <|v|>: {v_od:.1f} ± {v_od_err:.1f} km/s
• Observed ratio: {ratio:.4f} ± {ratio_err:.4f}

COMPARISON:
-----------
• Synchronism prediction: {pred_ratio:.4f}
• ΛCDM prediction: 1.0000
• Deviation from ΛCDM: {signif_lcdm:.1f}σ
• Deviation from Synchronism: {signif_sync:.1f}σ

BOOTSTRAP:
----------
• Mean ratio: {boot_mean:.4f} ± {boot_std:.4f}
• 68% CI: [{boot_16:.4f}, {boot_84:.4f}]

MODEL COMPARISON:
-----------------
• χ² (Synchronism): {chi2_sync:.1f}
• χ² (ΛCDM): {chi2_lcdm:.1f}
• Δχ² preference: {delta_chi2:.1f} → {np.sqrt(max(delta_chi2, 0)):.1f}σ for Synchronism

INTERPRETATION:
---------------
"""

if signif_lcdm > 3:
    summary += "SIGNIFICANT deviation from ΛCDM detected.\n"
    if abs(signif_sync) < 2:
        summary += "Results CONSISTENT with Synchronism prediction.\n"
    else:
        summary += "Magnitude differs from prediction - requires investigation.\n"
else:
    summary += f"Weak deviation from ΛCDM ({signif_lcdm:.1f}σ).\n"
    summary += "More data or better density estimates needed.\n"

print(summary)

print("\n" + "=" * 70)
print("SESSION #171 COMPLETE")
print("=" * 70)
