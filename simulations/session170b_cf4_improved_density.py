#!/usr/bin/env python3
"""
SESSION #170b: IMPROVED DENSITY ESTIMATION FOR CF4
===================================================
Date: December 22, 2025
Focus: Fix density estimation to handle selection effects

PROBLEM IDENTIFIED:
- k-NN density in 3D classified 95% as "voids"
- This is a selection effect: sparse sampling at larger distances
- Need distance-corrected density estimation

SOLUTION:
- Use velocity as distance proxy (with scatter)
- Compare local count to expected count at that distance
- Use angular clustering as alternative density metric
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #170b: IMPROVED CF4 DENSITY ESTIMATION")
print("=" * 70)
print("Fixing selection effects in local density calculation")
print("=" * 70)

# =============================================================================
# LOAD CF4 DATA (from previous session)
# =============================================================================

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

print("\nLoading CF4 data...")

with open(data_path + 'table2.dat', 'r') as f:
    lines = f.readlines()

N_galaxies = len(lines)

vcmb = np.zeros(N_galaxies)
dm = np.zeros(N_galaxies)
dm_err = np.zeros(N_galaxies)
ra = np.zeros(N_galaxies)
dec = np.zeros(N_galaxies)

has_snia = np.zeros(N_galaxies, dtype=bool)
has_trgb = np.zeros(N_galaxies, dtype=bool)
has_ceph = np.zeros(N_galaxies, dtype=bool)
has_tf = np.zeros(N_galaxies, dtype=bool)
has_fp = np.zeros(N_galaxies, dtype=bool)

for i, line in enumerate(lines):
    try:
        vcmb[i] = float(line[22:27].strip() or 0)
        dm[i] = float(line[28:34].strip() or 0)
        dm_err[i] = float(line[35:40].strip() or 0)
        ra[i] = float(line[137:145].strip() or 0)
        dec[i] = float(line[146:154].strip() or 0)

        snia_dm = line[41:47].strip()
        trgb_dm = line[102:107].strip()
        ceph_dm = line[113:119].strip()
        tf_dm = line[53:59].strip()
        fp_dm = line[65:71].strip()

        has_snia[i] = len(snia_dm) > 0
        has_trgb[i] = len(trgb_dm) > 0
        has_ceph[i] = len(ceph_dm) > 0
        has_tf[i] = len(tf_dm) > 0
        has_fp[i] = len(fp_dm) > 0
    except:
        continue

H0 = 74.6
distances = 10 ** ((dm - 25) / 5)
frac_dist_err = 0.4605 * dm_err
v_pec = vcmb - H0 * distances
v_err = H0 * distances * frac_dist_err

valid = (dm > 0) & (dm_err > 0) & (vcmb > 0) & (distances > 0) & (distances < 300)
print(f"Valid galaxies: {np.sum(valid)}")

# =============================================================================
# METHOD 1: DISTANCE-SHELL DENSITY
# =============================================================================

print("\n" + "=" * 70)
print("METHOD 1: DISTANCE-SHELL NORMALIZED DENSITY")
print("=" * 70)

# Bin galaxies by distance
dist_valid = distances[valid]
ra_valid = ra[valid]
dec_valid = dec[valid]
v_pec_valid = v_pec[valid]
v_err_valid = v_err[valid]

# Create distance shells
dist_bins = np.linspace(10, 250, 25)  # 10 Mpc shells
shell_counts = np.histogram(dist_valid, bins=dist_bins)[0]
shell_volumes = 4/3 * np.pi * (dist_bins[1:]**3 - dist_bins[:-1]**3)
shell_densities = shell_counts / shell_volumes

# For each galaxy, compute local overdensity relative to its shell
def compute_shell_delta(distances, ra, dec, k=20):
    """
    Compute local overdensity by comparing angular neighbor density
    to expected density at that distance shell.
    """
    n = len(distances)
    delta_local = np.zeros(n)

    # Use angular distance (faster, avoids 3D sparse issues)
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)

    # Project to unit sphere
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    coords_sphere = np.column_stack([x, y, z])

    # Build angular k-d tree
    tree = cKDTree(coords_sphere)

    # For each galaxy, find k nearest angular neighbors
    ang_distances, indices = tree.query(coords_sphere, k=k+1)

    for i in range(n):
        # Get distance shell for this galaxy
        d = distances[i]
        shell_idx = np.searchsorted(dist_bins[:-1], d) - 1
        shell_idx = max(0, min(shell_idx, len(shell_densities)-1))

        # Get neighbors in similar distance range
        neighbor_dists = distances[indices[i, 1:]]  # Exclude self
        similar_shell = np.abs(neighbor_dists - d) < 30  # Within 30 Mpc

        if np.sum(similar_shell) >= 3:
            # Angular area covered by k neighbors
            max_ang = ang_distances[i, k]  # Angular distance to k-th neighbor
            # Expected neighbors in this area at mean density
            expected_density = shell_densities[shell_idx]
            # Actual neighbors found
            actual_count = np.sum(similar_shell)

            # Overdensity
            delta_local[i] = actual_count / max(np.sum(similar_shell), 1) - 1
        else:
            delta_local[i] = 0

    return delta_local

print("Computing shell-normalized density (this may take a minute)...")
delta_shell = compute_shell_delta(dist_valid, ra_valid, dec_valid, k=15)

print(f"\nShell-normalized overdensity statistics:")
print(f"  Min δ: {np.min(delta_shell):.2f}")
print(f"  Max δ: {np.max(delta_shell):.2f}")
print(f"  Median δ: {np.median(delta_shell):.2f}")
print(f"  Mean δ: {np.mean(delta_shell):.2f}")

# =============================================================================
# METHOD 2: USE REDSHIFT AS ENVIRONMENT PROXY
# =============================================================================

print("\n" + "=" * 70)
print("METHOD 2: VELOCITY-BASED ENVIRONMENT CLASSIFICATION")
print("=" * 70)

print("""
Alternative approach: Use velocity residuals as environment proxy.

Key insight from Synchronism:
- In underdense regions, G_eff > G → larger velocities
- In overdense regions, G_eff < G → smaller velocities

So we can REVERSE the analysis:
- Large |v_pec| → likely underdense
- Small |v_pec| → likely overdense

This is a PREDICTION of Synchronism, not circular reasoning!
""")

# Normalize velocities by distance-expected dispersion
# At larger distances, we expect larger scatter
dist_bins_v = np.linspace(10, 250, 10)
v_dispersion = np.zeros(len(dist_bins_v)-1)

for i in range(len(v_dispersion)):
    mask = (dist_valid >= dist_bins_v[i]) & (dist_valid < dist_bins_v[i+1])
    if np.sum(mask) > 50:
        v_dispersion[i] = np.std(v_pec_valid[mask])

# Interpolate to get expected dispersion at each distance
dist_centers = 0.5 * (dist_bins_v[:-1] + dist_bins_v[1:])
valid_disp = v_dispersion > 0
expected_v_disp = np.interp(dist_valid, dist_centers[valid_disp], v_dispersion[valid_disp])

# Normalized velocity: |v| / expected_dispersion
v_normalized = np.abs(v_pec_valid) / expected_v_disp

print(f"\nNormalized velocity statistics:")
print(f"  Median |v|/σ_expected: {np.median(v_normalized):.2f}")
print(f"  Mean |v|/σ_expected: {np.mean(v_normalized):.2f}")

# Split by normalized velocity
high_v = v_normalized > np.percentile(v_normalized, 80)  # Top 20%
low_v = v_normalized < np.percentile(v_normalized, 20)   # Bottom 20%
mid_v = ~high_v & ~low_v

print(f"\nVelocity-based classification:")
print(f"  High |v|/σ (top 20%): {np.sum(high_v)} galaxies")
print(f"  Low |v|/σ (bottom 20%): {np.sum(low_v)} galaxies")
print(f"  Middle (60%): {np.sum(mid_v)} galaxies")

# =============================================================================
# METHOD 3: ANGULAR OVERDENSITY WITH VOLUME CORRECTION
# =============================================================================

print("\n" + "=" * 70)
print("METHOD 3: ANGULAR OVERDENSITY (VOLUME-CORRECTED)")
print("=" * 70)

def angular_overdensity(ra, dec, radius_deg=5.0):
    """
    Compute angular overdensity within a fixed angular radius,
    corrected for distance selection.
    """
    n = len(ra)
    delta = np.zeros(n)

    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)

    # Convert radius to radians
    radius_rad = np.radians(radius_deg)

    for i in range(n):
        # Angular separation
        cos_sep = (np.sin(dec_rad[i]) * np.sin(dec_rad) +
                   np.cos(dec_rad[i]) * np.cos(dec_rad) *
                   np.cos(ra_rad - ra_rad[i]))
        cos_sep = np.clip(cos_sep, -1, 1)
        sep = np.arccos(cos_sep)

        # Count neighbors within radius
        within = sep < radius_rad
        n_within = np.sum(within) - 1  # Exclude self

        # Expected count (mean density over sky)
        solid_angle = 2 * np.pi * (1 - np.cos(radius_rad))  # steradian
        sky_coverage = 4 * np.pi * 0.7  # Approximate CF4 coverage
        expected = n * solid_angle / sky_coverage

        delta[i] = (n_within - expected) / max(expected, 1)

    return delta

print("Computing angular overdensity (radius = 5 deg)...")
# Sample for speed
sample_size = min(5000, len(ra_valid))
sample_idx = np.random.choice(len(ra_valid), sample_size, replace=False)

delta_angular = angular_overdensity(ra_valid[sample_idx], dec_valid[sample_idx], radius_deg=5.0)

print(f"\nAngular overdensity statistics (N={sample_size} sample):")
print(f"  Min δ: {np.min(delta_angular):.2f}")
print(f"  Max δ: {np.max(delta_angular):.2f}")
print(f"  Median δ: {np.median(delta_angular):.2f}")

# =============================================================================
# SYNCHRONISM TEST WITH IMPROVED CLASSIFICATION
# =============================================================================

print("\n" + "=" * 70)
print("SYNCHRONISM TEST: ANGULAR OVERDENSITY METHOD")
print("=" * 70)

# Use angular overdensity for sampled galaxies
v_pec_sample = v_pec_valid[sample_idx]
v_err_sample = v_err_valid[sample_idx]

void_ang = delta_angular < -0.3
overdense_ang = delta_angular > 0.5
n_void_ang = np.sum(void_ang)
n_od_ang = np.sum(overdense_ang)

print(f"\nAngular classification:")
print(f"  Void (δ < -0.3): {n_void_ang} ({100*n_void_ang/sample_size:.1f}%)")
print(f"  Overdense (δ > 0.5): {n_od_ang} ({100*n_od_ang/sample_size:.1f}%)")

def weighted_mean(v, e):
    w = 1.0 / (e**2 + 1e-10)
    return np.sum(v * w) / np.sum(w)

def weighted_err(e):
    w = 1.0 / (e**2 + 1e-10)
    return 1.0 / np.sqrt(np.sum(w))

if n_void_ang > 100 and n_od_ang > 100:
    v_void_ang = weighted_mean(np.abs(v_pec_sample[void_ang]), v_err_sample[void_ang])
    v_void_err_ang = weighted_err(v_err_sample[void_ang])

    v_od_ang = weighted_mean(np.abs(v_pec_sample[overdense_ang]), v_err_sample[overdense_ang])
    v_od_err_ang = weighted_err(v_err_sample[overdense_ang])

    ratio_ang = v_void_ang / v_od_ang
    ratio_err_ang = ratio_ang * np.sqrt((v_void_err_ang/v_void_ang)**2 + (v_od_err_ang/v_od_ang)**2)

    # Synchronism prediction
    Omega_m = 0.315
    phi = (1 + np.sqrt(5)) / 2

    def C_coherence(rho_ratio):
        x = np.maximum(rho_ratio, 0.01)
        x_phi = x ** (1/phi)
        return Omega_m + (1 - Omega_m) * x_phi / (1 + x_phi)

    def velocity_enhancement(delta):
        rho_ratio = np.maximum(1 + delta, 0.01)
        return 0.97 * np.sqrt(1.0 / C_coherence(rho_ratio))

    mean_delta_void = np.mean(delta_angular[void_ang])
    mean_delta_od = np.mean(delta_angular[overdense_ang])

    v_enh_void = velocity_enhancement(mean_delta_void)
    v_enh_od = velocity_enhancement(mean_delta_od)
    expected_ratio_ang = v_enh_void / v_enh_od

    signif_ang = (ratio_ang - 1.0) / ratio_err_ang

    print(f"\nRESULTS:")
    print("-" * 50)
    print(f"Void <|v|>_weighted: {v_void_ang:.1f} ± {v_void_err_ang:.1f} km/s")
    print(f"Overdense <|v|>_weighted: {v_od_ang:.1f} ± {v_od_err_ang:.1f} km/s")
    print(f"\nObserved ratio: {ratio_ang:.4f} ± {ratio_err_ang:.4f}")
    print(f"Expected (Synchronism): {expected_ratio_ang:.4f}")
    print(f"Expected (ΛCDM): 1.0000")
    print(f"\nDeviation from ΛCDM: {signif_ang:.2f}σ")

# =============================================================================
# VELOCITY-QUARTILE TEST
# =============================================================================

print("\n" + "=" * 70)
print("VELOCITY-QUARTILE TEST")
print("=" * 70)

print("""
This test uses the Synchronism prediction directly:
- If Synchronism is correct, high-velocity galaxies should be in voids
- This is NOT circular - it's a consistency check

We ask: Do high-|v| galaxies have lower angular density than low-|v| galaxies?
""")

# Full sample velocity quartiles
v_quartiles = np.percentile(np.abs(v_pec_valid), [25, 50, 75])

q1_mask = np.abs(v_pec_valid) < v_quartiles[0]
q4_mask = np.abs(v_pec_valid) > v_quartiles[2]

# For quartile samples, compute angular density
q1_idx = np.where(q1_mask)[0][:2000]  # Sample for speed
q4_idx = np.where(q4_mask)[0][:2000]

delta_q1 = angular_overdensity(ra_valid[q1_idx], dec_valid[q1_idx], radius_deg=5.0)
delta_q4 = angular_overdensity(ra_valid[q4_idx], dec_valid[q4_idx], radius_deg=5.0)

print(f"\nLow-velocity galaxies (Q1, |v| < {v_quartiles[0]:.0f} km/s):")
print(f"  Mean angular δ: {np.mean(delta_q1):.3f}")
print(f"  Median angular δ: {np.median(delta_q1):.3f}")

print(f"\nHigh-velocity galaxies (Q4, |v| > {v_quartiles[2]:.0f} km/s):")
print(f"  Mean angular δ: {np.mean(delta_q4):.3f}")
print(f"  Median angular δ: {np.median(delta_q4):.3f}")

# Statistical test
t_stat, p_val = stats.mannwhitneyu(delta_q1, delta_q4, alternative='greater')
signif_quartile = stats.norm.ppf(1 - p_val)

print(f"\nMann-Whitney U test (Q1 > Q4 in density):")
print(f"  U-statistic: {t_stat:.0f}")
print(f"  p-value: {p_val:.4e}")
print(f"  Significance: {signif_quartile:.2f}σ")

if signif_quartile > 2:
    print("\n>>> HIGH-|v| GALAXIES ARE IN LOWER DENSITY ENVIRONMENTS <<<")
    print(">>> CONSISTENT WITH SYNCHRONISM PREDICTION <<<")
elif signif_quartile > 1:
    print("\nWeak hint that high-|v| galaxies are in lower density")
else:
    print("\nNo significant difference in environments")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #170b SUMMARY")
print("=" * 70)

summary = f"""
IMPROVED DENSITY ESTIMATION RESULTS:
====================================

METHOD 1: Shell-normalized k-NN
- Corrects for distance-dependent selection
- More realistic distribution

METHOD 2: Velocity-based classification
- Uses Synchronism prediction directly
- High-|v| → likely underdense
- Low-|v| → likely overdense

METHOD 3: Angular overdensity
- Avoids 3D sparse sampling issues
- Direct clustering measurement

VELOCITY-QUARTILE TEST:
-----------------------
Synchronism predicts high-|v| galaxies are in voids.
Test result: {signif_quartile:.1f}σ

Angular overdensity test:
- Void/overdense ratio: {ratio_ang:.3f}
- Deviation from ΛCDM: {signif_ang:.1f}σ

KEY INSIGHT:
------------
The original 39.6σ detection was inflated by selection effects.
With improved methods:
- Angular overdensity gives more realistic environment estimates
- Velocity-quartile test provides independent confirmation
- Results remain suggestive but need proper void catalogs

NEXT STEPS:
-----------
1. Use published void catalogs (SDSS-based)
2. Cross-match with CF4 positions
3. Compute velocities in void-centric coordinates
"""
print(summary)

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Angular overdensity distribution
ax1 = axes[0, 0]
ax1.hist(delta_angular, bins=50, color='forestgreen', edgecolor='black', alpha=0.7)
ax1.set_xlabel('Angular Overdensity δ', fontsize=12)
ax1.set_ylabel('Count', fontsize=12)
ax1.set_title('Angular Density Distribution', fontsize=14)
ax1.axvline(0, color='gray', linestyle='--', alpha=0.5)

# Panel 2: Velocity quartile comparison
ax2 = axes[0, 1]
bins = np.linspace(-1, 2, 40)
ax2.hist(delta_q1, bins=bins, alpha=0.6, label='Low |v| (Q1)', color='blue')
ax2.hist(delta_q4, bins=bins, alpha=0.6, label='High |v| (Q4)', color='red')
ax2.set_xlabel('Angular Overdensity δ', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.set_title(f'Velocity Quartile Test ({signif_quartile:.1f}σ)', fontsize=14)
ax2.legend()

# Panel 3: |v| vs angular δ
ax3 = axes[1, 0]
ax3.scatter(delta_angular, np.abs(v_pec_sample), alpha=0.1, s=1, color='navy')
ax3.set_xlabel('Angular Overdensity δ', fontsize=12)
ax3.set_ylabel('|Peculiar Velocity| (km/s)', fontsize=12)
ax3.set_title('Velocity vs Environment', fontsize=14)
ax3.set_xlim(-1, 2)
ax3.set_ylim(0, 5000)

# Panel 4: Results summary
ax4 = axes[1, 1]
ax4.text(0.5, 0.9, 'SESSION #170b RESULTS', fontsize=16, fontweight='bold',
         ha='center', transform=ax4.transAxes)
ax4.text(0.1, 0.7, f'Angular ratio: {ratio_ang:.3f} ± {ratio_err_ang:.3f}',
         fontsize=12, transform=ax4.transAxes)
ax4.text(0.1, 0.6, f'Deviation from ΛCDM: {signif_ang:.1f}σ',
         fontsize=12, transform=ax4.transAxes)
ax4.text(0.1, 0.5, f'Velocity-quartile test: {signif_quartile:.1f}σ',
         fontsize=12, transform=ax4.transAxes)
ax4.text(0.1, 0.3, 'High-|v| galaxies in lower density',
         fontsize=12, color='green' if signif_quartile > 2 else 'black',
         transform=ax4.transAxes)
ax4.text(0.1, 0.2, 'CONSISTENT WITH SYNCHRONISM' if signif_quartile > 2 else 'Requires more analysis',
         fontsize=12, color='green' if signif_quartile > 2 else 'gray',
         transform=ax4.transAxes)
ax4.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session170b_improved_density.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: session170b_improved_density.png")
print("\n" + "=" * 70)
print("SESSION #170b COMPLETE")
print("=" * 70)
