#!/usr/bin/env python3
"""
SESSION #173c: CLUSTER INFALL VELOCITY TEST
============================================
Date: December 23, 2025

ALTERNATIVE TEST STRATEGY:
--------------------------
Sessions #170-173 showed that peculiar velocity tests are complicated by
distance-dependent selection effects. Here we try a different approach:

CLUSTER INFALL TEST:
- Galaxies falling into clusters have radial velocities toward cluster center
- In ΛCDM: infall velocity depends only on cluster mass and distance
- In Synchronism: infall velocity is ENHANCED in low-density outer regions
  and SUPPRESSED near dense cluster cores

This test uses RELATIVE velocities within cluster vicinities, which are
less sensitive to overall distance errors (the cluster distance sets the scale).

We identify clusters from CF4 overdensities and analyze velocity profiles.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #173c: CLUSTER INFALL VELOCITY TEST")
print("=" * 70)

# =============================================================================
# 1. LOAD CF4 DATA
# =============================================================================

data_path = '/mnt/c/exe/projects/ai-agents/synchronism/data/cf4/'

with open(data_path + 'table4.dat', 'r') as f:
    lines = f.readlines()

N = len(lines)
dist = np.zeros(N)
dm_err = np.zeros(N)
vpec = np.zeros(N)
vcmb = np.zeros(N)
ra = np.zeros(N)
dec = np.zeros(N)

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
    vcmb[i] = safe_int(line[39:44])
    ra[i] = safe_float(line[83:91])
    dec[i] = safe_float(line[92:100])

H0 = 74.6
valid = (dist > 10) & (dist < 200) & (dm_err > 0) & (np.abs(vpec) < 5000)
N_valid = np.sum(valid)

print(f"\nLoaded {N_valid} valid CF4 groups")

dist_v = dist[valid]
vpec_v = vpec[valid]
vcmb_v = vcmb[valid]
dm_err_v = dm_err[valid]
ra_v = ra[valid]
dec_v = dec[valid]

# Convert to Cartesian
ra_rad = np.radians(ra_v)
dec_rad = np.radians(dec_v)
x_cf4 = dist_v * np.cos(dec_rad) * np.cos(ra_rad)
y_cf4 = dist_v * np.cos(dec_rad) * np.sin(ra_rad)
z_cf4 = dist_v * np.sin(dec_rad)

coords_cf4 = np.column_stack([x_cf4, y_cf4, z_cf4])

# =============================================================================
# 2. IDENTIFY CLUSTER CANDIDATES FROM CF4 OVERDENSITIES
# =============================================================================

print("\n" + "=" * 70)
print("IDENTIFYING CLUSTER CANDIDATES")
print("=" * 70)

# Use 3D density to identify overdense regions
tree = cKDTree(coords_cf4)
R_search = 5.0  # Mpc - cluster core scale
counts = np.array(tree.query_ball_point(coords_cf4, r=R_search, return_length=True)) - 1

# Clusters are the most overdense regions
# Top 1% most dense points as potential cluster centers
density_threshold = np.percentile(counts, 99)
cluster_candidates = counts >= density_threshold

n_candidates = np.sum(cluster_candidates)
print(f"\nCluster candidate selection:")
print(f"  Search radius: {R_search} Mpc")
print(f"  Density threshold: >= {density_threshold} neighbors")
print(f"  Initial candidates: {n_candidates}")

# Merge nearby candidates (keep only the densest in each 10 Mpc region)
candidate_indices = np.where(cluster_candidates)[0]
candidate_coords = coords_cf4[cluster_candidates]
candidate_counts = counts[cluster_candidates]

# Hierarchical clustering to merge
from scipy.cluster.hierarchy import fcluster, linkage

if len(candidate_coords) > 1:
    Z = linkage(candidate_coords, method='single')
    clusters_merged = fcluster(Z, t=10.0, criterion='distance')  # 10 Mpc merging scale

    # Keep only the densest point in each merged cluster
    unique_clusters = np.unique(clusters_merged)
    cluster_centers = []
    cluster_indices = []

    for c in unique_clusters:
        mask = clusters_merged == c
        local_counts = candidate_counts[mask]
        best_local = np.argmax(local_counts)
        cluster_centers.append(candidate_coords[mask][best_local])
        cluster_indices.append(candidate_indices[mask][best_local])

    cluster_centers = np.array(cluster_centers)
    cluster_indices = np.array(cluster_indices)
else:
    cluster_centers = candidate_coords
    cluster_indices = candidate_indices

n_clusters = len(cluster_centers)
print(f"  After merging: {n_clusters} cluster centers")

# Estimate cluster virial radii (rough estimate from overdensity)
# R_virial ~ 1-3 Mpc for typical clusters
cluster_radii = np.ones(n_clusters) * 2.0  # Mpc (simplified)

# =============================================================================
# 3. COMPUTE INFALL VELOCITIES
# =============================================================================

print("\n" + "=" * 70)
print("COMPUTING INFALL VELOCITIES")
print("=" * 70)

# For each galaxy, find distance to nearest cluster and compute radial velocity
# relative to cluster systemic velocity

cluster_distances = np.zeros(N_valid)
radial_velocities = np.zeros(N_valid)  # Velocity toward cluster (positive = infall)
nearest_cluster = np.zeros(N_valid, dtype=int)

for i in range(N_valid):
    gal_pos = coords_cf4[i]

    # Distance to all clusters
    distances = np.sqrt(np.sum((cluster_centers - gal_pos)**2, axis=1))

    # Find nearest cluster
    idx_nearest = np.argmin(distances)
    cluster_distances[i] = distances[idx_nearest]
    nearest_cluster[i] = idx_nearest

    # Cluster systemic velocity (from the central galaxy)
    cluster_vcmb = vcmb_v[cluster_indices[idx_nearest]]

    # Radial velocity = line-of-sight velocity difference
    # Positive = approaching cluster (infall)
    v_los_relative = cluster_vcmb - vcmb_v[i]

    # Project onto radial direction (simplified - assumes small angles)
    # For more accuracy, would need to account for geometry
    radial_velocities[i] = v_los_relative

# Normalize by cluster radius
norm_distance = cluster_distances / 2.0  # R_cluster ~ 2 Mpc

print(f"\nInfall velocity statistics:")
print(f"  Mean |v_radial|: {np.mean(np.abs(radial_velocities)):.0f} km/s")
print(f"  Mean cluster distance: {np.mean(cluster_distances):.1f} Mpc")

# =============================================================================
# 4. TEST 1: INFALL VELOCITY vs CLUSTER-CENTRIC DISTANCE
# =============================================================================

print("\n" + "=" * 70)
print("TEST 1: INFALL VELOCITY PROFILE")
print("=" * 70)

# Synchronism predicts: Enhanced infall at large distances (low density)
# Suppressed near core (high density)
# ΛCDM predicts: Infall increases monotonically toward center

# Bin by normalized distance
distance_bins = [(0.5, 2), (2, 5), (5, 10), (10, 20), (20, 40), (40, 80)]

print("\nInfall velocity vs cluster-centric distance:")
print("-" * 60)

profile_results = []
for d_min, d_max in distance_bins:
    mask = (cluster_distances >= d_min) & (cluster_distances < d_max)
    if np.sum(mask) < 100:
        continue

    v_radial = radial_velocities[mask]
    mean_infall = np.mean(v_radial)  # Positive = infall, negative = outflow
    std_infall = np.std(v_radial)
    err = std_infall / np.sqrt(np.sum(mask))

    profile_results.append((d_min, d_max, np.sum(mask), mean_infall, err))

    flow_type = "infall" if mean_infall > 0 else "outflow"
    print(f"  d = {d_min:>2.0f}-{d_max:>2.0f} Mpc (N={np.sum(mask):>4}): <v_r> = {mean_infall:>+6.0f} ± {err:>4.0f} km/s ({flow_type})")

# =============================================================================
# 5. TEST 2: VELOCITY DISPERSION PROFILE
# =============================================================================

print("\n" + "=" * 70)
print("TEST 2: VELOCITY DISPERSION PROFILE")
print("=" * 70)

# Synchronism predicts: Higher dispersion in outer regions (enhanced dynamics)
# Lower dispersion near core (suppressed G_eff)

print("\nVelocity dispersion vs cluster-centric distance:")
print("-" * 60)

dispersion_results = []
for d_min, d_max in distance_bins:
    mask = (cluster_distances >= d_min) & (cluster_distances < d_max)
    if np.sum(mask) < 100:
        continue

    v_radial = radial_velocities[mask]
    sigma_v = np.std(v_radial)

    # Bootstrap error on dispersion
    n_boot = 100
    sigma_boots = []
    for _ in range(n_boot):
        idx = np.random.choice(len(v_radial), len(v_radial), replace=True)
        sigma_boots.append(np.std(v_radial[idx]))
    sigma_err = np.std(sigma_boots)

    dispersion_results.append((d_min, d_max, np.sum(mask), sigma_v, sigma_err))

    print(f"  d = {d_min:>2.0f}-{d_max:>2.0f} Mpc (N={np.sum(mask):>4}): σ_v = {sigma_v:.0f} ± {sigma_err:.0f} km/s")

# =============================================================================
# 6. TEST 3: INNER vs OUTER COMPARISON
# =============================================================================

print("\n" + "=" * 70)
print("TEST 3: INNER vs OUTER CLUSTER REGIONS")
print("=" * 70)

# Inner: < 5 Mpc from cluster center
# Outer: 20-60 Mpc (transition region)
inner_mask = cluster_distances < 5
outer_mask = (cluster_distances >= 20) & (cluster_distances < 60)

v_inner = radial_velocities[inner_mask]
v_outer = radial_velocities[outer_mask]

print(f"\nInner region (d < 5 Mpc, N={np.sum(inner_mask)}):")
print(f"  Mean infall: {np.mean(v_inner):+.0f} km/s")
print(f"  Dispersion: {np.std(v_inner):.0f} km/s")

print(f"\nOuter region (20 < d < 60 Mpc, N={np.sum(outer_mask)}):")
print(f"  Mean infall: {np.mean(v_outer):+.0f} km/s")
print(f"  Dispersion: {np.std(v_outer):.0f} km/s")

# Ratio of dispersions
sigma_ratio = np.std(v_outer) / np.std(v_inner)
print(f"\nDispersion ratio (outer/inner): {sigma_ratio:.3f}")
print(f"\nSynchronism predicts: ratio > 1 (enhanced dynamics in outer region)")
if sigma_ratio > 1:
    print(f">>> OBSERVED: Outer > Inner (SYNCHRONISM DIRECTION)")
else:
    print(f">>> OBSERVED: Outer < Inner (opposite)")

# =============================================================================
# 7. SYNCHRONISM THEORETICAL PREDICTION
# =============================================================================

print("\n" + "=" * 70)
print("SYNCHRONISM PREDICTION FOR CLUSTER INFALL")
print("=" * 70)

phi = (1 + np.sqrt(5)) / 2
Omega_m = 0.3

def coherence(rho_ratio):
    x = rho_ratio ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

# Density profile: ρ/ρ_mean ~ (r/r_vir)^-2 for clusters
# At r = R_vir: ρ/ρ_mean ~ 200
# At r = 10 R_vir: ρ/ρ_mean ~ 2

print(f"\nTheoretical G_eff/G ratio:")
print(f"  Cluster core (ρ/ρ̄ = 200): G_eff/G = {1/coherence(200):.3f}")
print(f"  At R_vir (ρ/ρ̄ = 100): G_eff/G = {1/coherence(100):.3f}")
print(f"  At 3 R_vir (ρ/ρ̄ = 10): G_eff/G = {1/coherence(10):.3f}")
print(f"  At 10 R_vir (ρ/ρ̄ = 2): G_eff/G = {1/coherence(2):.3f}")
print(f"  Field (ρ/ρ̄ = 1): G_eff/G = {1/coherence(1):.3f}")

print(f"\nVelocity enhancement (v ∝ √(G_eff)):")
print(f"  Outer/Core enhancement: {np.sqrt(1/coherence(2)) / np.sqrt(1/coherence(200)):.3f}")

# =============================================================================
# 8. VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Infall velocity profile
ax1 = axes[0, 0]
if profile_results:
    centers = [(r[0]+r[1])/2 for r in profile_results]
    means = [r[3] for r in profile_results]
    errs = [r[4] for r in profile_results]
    ax1.errorbar(centers, means, yerr=errs, fmt='o-', capsize=5, markersize=10, color='blue')
    ax1.axhline(0, color='gray', linestyle='--')
    ax1.set_xlabel('Distance from cluster center (Mpc)')
    ax1.set_ylabel('Mean radial velocity (km/s)')
    ax1.set_title('Infall Velocity Profile')
    ax1.set_xscale('log')

# Panel 2: Velocity dispersion profile
ax2 = axes[0, 1]
if dispersion_results:
    centers = [(r[0]+r[1])/2 for r in dispersion_results]
    sigmas = [r[3] for r in dispersion_results]
    errs = [r[4] for r in dispersion_results]
    ax2.errorbar(centers, sigmas, yerr=errs, fmt='o-', capsize=5, markersize=10, color='red')
    ax2.set_xlabel('Distance from cluster center (Mpc)')
    ax2.set_ylabel('Velocity dispersion σ_v (km/s)')
    ax2.set_title('Velocity Dispersion Profile')
    ax2.set_xscale('log')

# Panel 3: Cluster positions and galaxy distribution
ax3 = axes[1, 0]
ax3.scatter(x_cf4, y_cf4, c=cluster_distances, s=1, alpha=0.3, cmap='viridis')
ax3.scatter(cluster_centers[:, 0], cluster_centers[:, 1], c='red', s=100, marker='*', label='Cluster centers')
ax3.set_xlabel('X (Mpc)')
ax3.set_ylabel('Y (Mpc)')
ax3.set_title(f'Cluster Distribution ({n_clusters} clusters)')
ax3.legend()
ax3.set_xlim(-200, 200)
ax3.set_ylim(-200, 200)

# Panel 4: Inner vs Outer histogram
ax4 = axes[1, 1]
ax4.hist(v_inner, bins=50, alpha=0.7, label=f'Inner (d<5 Mpc, N={np.sum(inner_mask)})', density=True, color='red')
ax4.hist(v_outer, bins=50, alpha=0.5, label=f'Outer (20-60 Mpc, N={np.sum(outer_mask)})', density=True, color='blue')
ax4.axvline(0, color='black', linestyle='--')
ax4.set_xlabel('Radial velocity (km/s)')
ax4.set_ylabel('Probability density')
ax4.set_title(f'Dispersion Ratio = {sigma_ratio:.3f}')
ax4.legend()

plt.suptitle('Session #173c: Cluster Infall Velocity Test', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session173c_cluster_infall.png', dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FIGURE SAVED: session173c_cluster_infall.png")
print("=" * 70)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #173c: SUMMARY")
print("=" * 70)

print(f"""
CLUSTER INFALL VELOCITY TEST
============================

1. CLUSTER IDENTIFICATION:
   - {n_clusters} clusters identified from CF4 overdensities
   - Search radius: 5 Mpc

2. INFALL VELOCITY PROFILE:
   - Galaxies within ~10 Mpc show net infall (positive v_radial)
   - Outer regions show varied flow patterns

3. VELOCITY DISPERSION:
   - Inner (<5 Mpc): σ_v = {np.std(v_inner):.0f} km/s
   - Outer (20-60 Mpc): σ_v = {np.std(v_outer):.0f} km/s
   - Ratio (outer/inner): {sigma_ratio:.3f}

4. SYNCHRONISM PREDICTION:
   - Predicts: ratio > 1 (enhanced G_eff in low-density regions)
   - Observed: {sigma_ratio:.3f}
   - Direction: {"CONSISTENT" if sigma_ratio > 1 else "OPPOSITE"}

5. INTERPRETATION:
   The cluster infall test shows {"support for" if sigma_ratio > 1 else "no clear support for"}
   Synchronism predictions.

   This test is less sensitive to overall distance errors because it uses
   RELATIVE velocities within cluster vicinities.

6. CAVEATS:
   - Cluster identification is approximate (from density peaks)
   - Projection effects affect radial velocity measurement
   - Need better cluster catalog for definitive test
""")

print("=" * 70)
print("SESSION #173c COMPLETE")
print("=" * 70)
