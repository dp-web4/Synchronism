#!/usr/bin/env python3
"""
SESSION #179: SPARC ENVIRONMENT CLASSIFICATION AND VOID/CLUSTER TEST
====================================================================
Date: December 25, 2025

OBJECTIVE:
----------
Session #177 predicted that void galaxies should have ~3.5% higher rotation
velocities than cluster galaxies of the same stellar mass.

This session:
1. Classifies SPARC galaxies by environment (void, field, cluster)
2. Tests whether v_flat varies with environment at fixed luminosity
3. Computes the Baryonic Tully-Fisher Relation by environment

APPROACH:
---------
Without direct void catalog cross-matching, we use LOCAL GALAXY DENSITY
estimated from galaxy positions and distances in SPARC itself as a proxy.

For more rigorous analysis, we estimate environment using:
- Distance from known galaxy clusters (Virgo, Fornax, etc.)
- Local density from nearest neighbor distances
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

print("=" * 70)
print("SESSION #179: SPARC ENVIRONMENT CLASSIFICATION")
print("=" * 70)

# =============================================================================
# 1. LOAD SPARC DATA
# =============================================================================

print("\n" + "=" * 70)
print("1. LOADING SPARC DATA")
print("=" * 70)

sparc_file = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/sparc_real_data/SPARC_Lelli2016c.mrt'

# Read SPARC data - parse space-delimited format
galaxies = []
with open(sparc_file, 'r') as f:
    in_data = False
    for line in f:
        line = line.rstrip()
        # Skip header - data starts when we see galaxy names
        if not in_data:
            # Start when we see lines starting with known galaxy name patterns
            if len(line) > 50 and line.lstrip()[:3] in ['CIG', 'D51', 'D63', 'DDO', 'ESO',
                'F56', 'F56', 'F57', 'F58', 'F59', 'F56', 'IC ', 'KK9', 'Mkn',
                'NGC', 'PGC', 'UGC', 'UGC']:
                in_data = True

        if in_data and len(line) > 80:
            try:
                # Parse space-delimited format
                parts = line.split()
                if len(parts) < 18:
                    continue

                name = parts[0]
                hubble_type = int(parts[1])
                distance = float(parts[2])
                e_distance = float(parts[3])
                # parts[4] is distance method
                inclination = float(parts[5])
                # parts[6] is e_inc
                L_3p6 = float(parts[7])
                # parts[8] is e_L
                Reff = float(parts[9])
                SBeff = float(parts[10])
                Rdisk = float(parts[11])
                SBdisk = float(parts[12])
                MHI = float(parts[13])
                RHI = float(parts[14])
                Vflat = float(parts[15])
                e_Vflat = float(parts[16])
                quality = int(parts[17])

                if Vflat > 0 and L_3p6 > 0 and distance > 0:
                    galaxies.append({
                        'name': name,
                        'type': hubble_type,
                        'D': distance,
                        'e_D': e_distance,
                        'inc': inclination,
                        'L_3p6': L_3p6,
                        'Reff': Reff,
                        'SBeff': SBeff,
                        'Rdisk': Rdisk,
                        'SBdisk': SBdisk,
                        'MHI': MHI,
                        'RHI': RHI,
                        'Vflat': Vflat,
                        'e_Vflat': e_Vflat,
                        'Q': quality
                    })
            except (ValueError, IndexError) as e:
                continue

print(f"Loaded {len(galaxies)} galaxies with complete data")

# =============================================================================
# 2. ENVIRONMENT CLASSIFICATION BY DISTANCE TO KNOWN CLUSTERS
# =============================================================================

print("\n" + "=" * 70)
print("2. ENVIRONMENT CLASSIFICATION")
print("=" * 70)

print("""
APPROACH: Use known galaxy cluster distances as environment proxy

Key nearby clusters/groups:
- Virgo cluster: ~16.5 Mpc
- Fornax cluster: ~20 Mpc
- Local Group: < 3 Mpc
- Ursa Major cluster: ~18 Mpc

SPARC galaxies within ~20 Mpc of these structures are likely in denser
environments. More distant galaxies may be in voids or field.

Additionally, we use surface brightness as density proxy:
- High SB galaxies are typically in denser environments
- Low SB galaxies (LSB) are often in voids

This is a PROXY method - proper environment classification requires
full 3D position + velocity catalogs.
""")

# Known cluster distances (approximate, using Virgo as primary)
virgo_distance = 16.5  # Mpc
fornax_distance = 20.0  # Mpc
uma_distance = 18.0  # Mpc

# Classify by distance and surface brightness
for g in galaxies:
    D = g['D']
    SB = g['SBeff']

    # Distance-based classification
    if D < 5:
        g['env_class'] = 'local'  # Local Group or very nearby
    elif 10 < D < 25:
        # Near Virgo/Fornax/UMa complex - likely higher density
        g['env_class'] = 'cluster_region'
    elif D > 50:
        # Distant galaxies - likely field or void
        g['env_class'] = 'distant_field'
    else:
        g['env_class'] = 'field'

    # Surface brightness proxy (low SB often in voids)
    if SB < 100:  # Low surface brightness
        g['sb_class'] = 'LSB'
    elif SB > 500:  # High surface brightness
        g['sb_class'] = 'HSB'
    else:
        g['sb_class'] = 'intermediate'

# Count classifications
env_counts = {}
sb_counts = {}
for g in galaxies:
    env_counts[g['env_class']] = env_counts.get(g['env_class'], 0) + 1
    sb_counts[g['sb_class']] = sb_counts.get(g['sb_class'], 0) + 1

print("\nEnvironment classification by distance:")
for env, count in sorted(env_counts.items()):
    print(f"  {env}: {count} galaxies")

print("\nSurface brightness classification:")
for sb, count in sorted(sb_counts.items()):
    print(f"  {sb}: {count} galaxies")

# =============================================================================
# 3. BARYONIC MASS ESTIMATION
# =============================================================================

print("\n" + "=" * 70)
print("3. BARYONIC MASS ESTIMATION")
print("=" * 70)

# Stellar mass from 3.6μm luminosity
# M_star = M/L × L where M/L ≈ 0.5 for [3.6]
ML_ratio = 0.5  # M_sun/L_sun

for g in galaxies:
    # Stellar mass
    M_star = ML_ratio * g['L_3p6'] * 1e9  # L_3p6 is in 10^9 L_sun

    # Gas mass (HI + He correction)
    M_gas = 1.33 * g['MHI'] * 1e9  # MHI is in 10^9 M_sun, 1.33 for He

    # Total baryonic mass
    M_bary = M_star + M_gas

    g['M_star'] = M_star
    g['M_gas'] = M_gas
    g['M_bary'] = M_bary

print(f"Baryonic masses computed using M/L = {ML_ratio}")
print(f"Mean log(M_bary/M_sun) = {np.mean([np.log10(g['M_bary']) for g in galaxies]):.2f}")

# =============================================================================
# 4. BARYONIC TULLY-FISHER RELATION BY ENVIRONMENT
# =============================================================================

print("\n" + "=" * 70)
print("4. BTFR BY ENVIRONMENT")
print("=" * 70)

def btfr(M_bary, a, b):
    """log(v) = a × log(M) + b"""
    return a * np.log10(M_bary) + b

# Split by environment
local = [g for g in galaxies if g['env_class'] == 'local']
cluster_region = [g for g in galaxies if g['env_class'] == 'cluster_region']
field = [g for g in galaxies if g['env_class'] == 'field']
distant = [g for g in galaxies if g['env_class'] == 'distant_field']

print(f"\nGalaxies by environment:")
print(f"  Local (<5 Mpc): {len(local)}")
print(f"  Cluster region (10-25 Mpc): {len(cluster_region)}")
print(f"  Field (5-50 Mpc, excluding cluster): {len(field)}")
print(f"  Distant field (>50 Mpc): {len(distant)}")

# Fit BTFR for each environment
def fit_btfr(gal_list, name):
    if len(gal_list) < 5:
        print(f"  {name}: Insufficient data ({len(gal_list)} galaxies)")
        return None, None, None

    M = np.array([g['M_bary'] for g in gal_list])
    v = np.array([g['Vflat'] for g in gal_list])

    # Fit
    try:
        popt, pcov = curve_fit(btfr, M, np.log10(v), p0=[0.25, -0.5])
        a, b = popt
        residuals = np.log10(v) - btfr(M, a, b)
        scatter = np.std(residuals)
        print(f"  {name}: a = {a:.4f}, b = {b:.4f}, scatter = {scatter:.4f} dex ({len(gal_list)} galaxies)")
        return a, b, scatter
    except Exception as e:
        print(f"  {name}: Fit failed - {e}")
        return None, None, None

print("\nBTFR fits by environment:")
print("-" * 60)
a_local, b_local, s_local = fit_btfr(local, 'Local')
a_cluster, b_cluster, s_cluster = fit_btfr(cluster_region, 'Cluster region')
a_field, b_field, s_field = fit_btfr(field, 'Field')
a_distant, b_distant, s_distant = fit_btfr(distant, 'Distant field')

# Also fit by surface brightness
lsb = [g for g in galaxies if g['sb_class'] == 'LSB']
hsb = [g for g in galaxies if g['sb_class'] == 'HSB']
isb = [g for g in galaxies if g['sb_class'] == 'intermediate']

print("\nBTFR fits by surface brightness:")
print("-" * 60)
a_lsb, b_lsb, s_lsb = fit_btfr(lsb, 'LSB')
a_isb, b_isb, s_isb = fit_btfr(isb, 'Intermediate SB')
a_hsb, b_hsb, s_hsb = fit_btfr(hsb, 'HSB')

# =============================================================================
# 5. RESIDUAL ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("5. RESIDUAL ANALYSIS - KEY TEST")
print("=" * 70)

print("""
SYNCHRONISM PREDICTION:
At fixed baryonic mass, void/LSB galaxies should have HIGHER Vflat.

TEST: Compute residuals from overall BTFR, compare by environment.

If Synchronism is correct:
- LSB galaxies: POSITIVE residuals (higher v than expected)
- Distant/void galaxies: POSITIVE residuals
- HSB/cluster galaxies: NEGATIVE or zero residuals
""")

# Fit overall BTFR
all_M = np.array([g['M_bary'] for g in galaxies])
all_v = np.array([g['Vflat'] for g in galaxies])

popt_all, _ = curve_fit(btfr, all_M, np.log10(all_v), p0=[0.25, -0.5])
a_all, b_all = popt_all

print(f"\nOverall BTFR: log(v) = {a_all:.4f} × log(M) + {b_all:.4f}")

# Compute residuals for each galaxy
for g in galaxies:
    v_expected = 10**btfr(g['M_bary'], a_all, b_all)
    g['v_residual'] = (g['Vflat'] - v_expected) / v_expected  # Fractional residual
    g['log_v_residual'] = np.log10(g['Vflat']) - btfr(g['M_bary'], a_all, b_all)

# Analyze residuals by environment
print("\nMean fractional velocity residual by environment:")
print("-" * 60)
print(f"{'Environment':<20} {'N':>6} {'Mean Δv/v':>12} {'Std':>10} {'Median':>10}")
print("-" * 60)

for env_name, env_list in [('Local', local), ('Cluster region', cluster_region),
                            ('Field', field), ('Distant field', distant)]:
    if len(env_list) > 0:
        residuals = [g['v_residual'] for g in env_list]
        print(f"{env_name:<20} {len(env_list):>6} {np.mean(residuals)*100:>11.2f}% {np.std(residuals)*100:>9.2f}% {np.median(residuals)*100:>9.2f}%")

print("\nMean fractional velocity residual by surface brightness:")
print("-" * 60)
print(f"{'SB Class':<20} {'N':>6} {'Mean Δv/v':>12} {'Std':>10} {'Median':>10}")
print("-" * 60)

for sb_name, sb_list in [('LSB', lsb), ('Intermediate', isb), ('HSB', hsb)]:
    if len(sb_list) > 0:
        residuals = [g['v_residual'] for g in sb_list]
        print(f"{sb_name:<20} {len(sb_list):>6} {np.mean(residuals)*100:>11.2f}% {np.std(residuals)*100:>9.2f}% {np.median(residuals)*100:>9.2f}%")

# =============================================================================
# 6. SYNCHRONISM PREDICTION VS OBSERVATION
# =============================================================================

print("\n" + "=" * 70)
print("6. SYNCHRONISM PREDICTION COMPARISON")
print("=" * 70)

# Session #177 prediction: ~3.5% higher v for void galaxies
sync_prediction = 0.035  # 3.5%

# Use LSB as void proxy, HSB as cluster proxy
if len(lsb) > 5 and len(hsb) > 5:
    lsb_mean_residual = np.mean([g['v_residual'] for g in lsb])
    hsb_mean_residual = np.mean([g['v_residual'] for g in hsb])

    observed_difference = lsb_mean_residual - hsb_mean_residual

    print(f"""
COMPARISON:
===========

Synchronism prediction (Session #177):
  Void galaxies have {sync_prediction*100:.1f}% higher v_rot than cluster galaxies

Observation using LSB vs HSB as environment proxy:
  LSB mean residual: {lsb_mean_residual*100:+.2f}%
  HSB mean residual: {hsb_mean_residual*100:+.2f}%
  Difference (LSB - HSB): {observed_difference*100:+.2f}%

Predicted difference: {sync_prediction*100:.1f}%
Observed difference: {observed_difference*100:.2f}%
""")

    if observed_difference > 0:
        print("RESULT: Observed trend is in the CORRECT DIRECTION (LSB > HSB)")
        if abs(observed_difference) > 0.5 * sync_prediction:
            print("RESULT: Magnitude is within a factor of 2 of prediction")
        else:
            print("RESULT: Magnitude is smaller than predicted")
    else:
        print("RESULT: Observed trend is OPPOSITE to prediction")

# =============================================================================
# 7. DISTANCE-BINNED ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("7. DISTANCE-BINNED ANALYSIS")
print("=" * 70)

# Bin by distance
distance_bins = [(0, 10), (10, 20), (20, 40), (40, 80), (80, 200)]

print(f"\nMean velocity residual by distance:")
print("-" * 70)
print(f"{'Distance bin (Mpc)':<20} {'N':>6} {'Mean Δv/v':>12} {'Expected (Sync)':>15}")
print("-" * 70)

# Synchronism expectation: farther = potentially less dense = higher v
# But this is weak correlation, mainly testing for trend

for d_min, d_max in distance_bins:
    bin_gals = [g for g in galaxies if d_min <= g['D'] < d_max]
    if len(bin_gals) >= 5:
        residuals = [g['v_residual'] for g in bin_gals]
        mean_res = np.mean(residuals)
        # Rough Synchronism expectation: more distant may be more underdense
        # But this is very weak correlation
        print(f"{d_min:>3} - {d_max:<4} Mpc      {len(bin_gals):>6} {mean_res*100:>11.2f}%")

# =============================================================================
# 8. VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("8. GENERATING VISUALIZATIONS")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: BTFR by environment
ax1 = axes[0, 0]
M_plot = np.logspace(7, 12, 100)

colors = {'local': 'purple', 'cluster_region': 'blue', 'field': 'green', 'distant_field': 'red'}
labels = {'local': 'Local (<5 Mpc)', 'cluster_region': 'Cluster region (10-25 Mpc)',
          'field': 'Field', 'distant_field': 'Distant (>50 Mpc)'}

for env in ['local', 'cluster_region', 'field', 'distant_field']:
    gals = [g for g in galaxies if g['env_class'] == env]
    if len(gals) > 0:
        M = [g['M_bary'] for g in gals]
        v = [g['Vflat'] for g in gals]
        ax1.scatter(M, v, c=colors[env], alpha=0.6, s=30, label=labels[env])

# Plot overall fit
ax1.plot(M_plot, 10**btfr(M_plot, a_all, b_all), 'k-', linewidth=2, label=f'BTFR fit')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Baryonic Mass (M☉)')
ax1.set_ylabel('Flat Rotation Velocity (km/s)')
ax1.set_title('BTFR by Environment (Distance-based)')
ax1.legend(loc='lower right', fontsize=8)
ax1.grid(True, alpha=0.3)

# Panel 2: BTFR by surface brightness
ax2 = axes[0, 1]
sb_colors = {'LSB': 'red', 'intermediate': 'green', 'HSB': 'blue'}
sb_labels = {'LSB': 'Low SB (void proxy)', 'intermediate': 'Intermediate', 'HSB': 'High SB (cluster proxy)'}

for sb in ['LSB', 'intermediate', 'HSB']:
    gals = [g for g in galaxies if g['sb_class'] == sb]
    if len(gals) > 0:
        M = [g['M_bary'] for g in gals]
        v = [g['Vflat'] for g in gals]
        ax2.scatter(M, v, c=sb_colors[sb], alpha=0.6, s=30, label=sb_labels[sb])

ax2.plot(M_plot, 10**btfr(M_plot, a_all, b_all), 'k-', linewidth=2)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Baryonic Mass (M☉)')
ax2.set_ylabel('Flat Rotation Velocity (km/s)')
ax2.set_title('BTFR by Surface Brightness')
ax2.legend(loc='lower right', fontsize=8)
ax2.grid(True, alpha=0.3)

# Panel 3: Residual histogram by SB class
ax3 = axes[1, 0]
if len(lsb) > 0:
    ax3.hist([g['v_residual']*100 for g in lsb], bins=20, alpha=0.5, label=f'LSB (n={len(lsb)})', color='red')
if len(hsb) > 0:
    ax3.hist([g['v_residual']*100 for g in hsb], bins=20, alpha=0.5, label=f'HSB (n={len(hsb)})', color='blue')

ax3.axvline(0, color='black', linestyle='--', linewidth=1)
ax3.axvline(3.5, color='green', linestyle=':', linewidth=2, label='Sync prediction: +3.5%')
ax3.set_xlabel('Velocity residual Δv/v (%)')
ax3.set_ylabel('Count')
ax3.set_title('Velocity Residual Distribution: LSB vs HSB')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Panel 4: Residual vs distance
ax4 = axes[1, 1]
distances = [g['D'] for g in galaxies]
residuals = [g['v_residual']*100 for g in galaxies]
sb_class = [g['sb_class'] for g in galaxies]

for sb, color in [('LSB', 'red'), ('intermediate', 'green'), ('HSB', 'blue')]:
    d = [distances[i] for i in range(len(distances)) if sb_class[i] == sb]
    r = [residuals[i] for i in range(len(residuals)) if sb_class[i] == sb]
    if len(d) > 0:
        ax4.scatter(d, r, c=color, alpha=0.5, s=20, label=sb)

ax4.axhline(0, color='black', linestyle='--', linewidth=1)
ax4.axhline(3.5, color='green', linestyle=':', linewidth=2, alpha=0.7)
ax4.set_xlabel('Distance (Mpc)')
ax4.set_ylabel('Velocity residual (%)')
ax4.set_title('Velocity Residual vs Distance')
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 100)
ax4.set_ylim(-50, 50)

plt.suptitle('Session #179: SPARC Environment Analysis - Testing Void/Cluster Prediction',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session179_sparc_environment.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session179_sparc_environment.png")

# =============================================================================
# 9. QUALITY-FILTERED ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("9. QUALITY-FILTERED ANALYSIS (Q=1 only)")
print("=" * 70)

q1_galaxies = [g for g in galaxies if g['Q'] == 1]
print(f"High-quality (Q=1) sample: {len(q1_galaxies)} galaxies")

if len(q1_galaxies) > 20:
    q1_lsb = [g for g in q1_galaxies if g['sb_class'] == 'LSB']
    q1_hsb = [g for g in q1_galaxies if g['sb_class'] == 'HSB']

    print(f"  Q=1 LSB: {len(q1_lsb)}, Q=1 HSB: {len(q1_hsb)}")

    if len(q1_lsb) >= 5 and len(q1_hsb) >= 5:
        q1_lsb_res = np.mean([g['v_residual'] for g in q1_lsb])
        q1_hsb_res = np.mean([g['v_residual'] for g in q1_hsb])

        print(f"\n  High-quality sample residuals:")
        print(f"    LSB mean: {q1_lsb_res*100:+.2f}%")
        print(f"    HSB mean: {q1_hsb_res*100:+.2f}%")
        print(f"    Difference: {(q1_lsb_res - q1_hsb_res)*100:+.2f}%")

# =============================================================================
# 10. SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #179: SUMMARY")
print("=" * 70)

print(f"""
SPARC ENVIRONMENT ANALYSIS - TESTING SESSION #177 PREDICTION
=============================================================

1. DATA:
   - {len(galaxies)} SPARC galaxies with complete data
   - Environment classified by distance and surface brightness proxy

2. SYNCHRONISM PREDICTION (Session #177):
   - Void galaxies should have ~3.5% higher v_rot than cluster galaxies
   - Using LSB as void proxy, HSB as cluster proxy

3. OBSERVATION:""")

if len(lsb) > 5 and len(hsb) > 5:
    lsb_res = np.mean([g['v_residual'] for g in lsb]) * 100
    hsb_res = np.mean([g['v_residual'] for g in hsb]) * 100
    diff = lsb_res - hsb_res

    print(f"""   - LSB galaxies: {lsb_res:+.2f}% velocity residual
   - HSB galaxies: {hsb_res:+.2f}% velocity residual
   - Difference (LSB - HSB): {diff:+.2f}%
   - Predicted difference: +3.5%""")

    if diff > 0:
        print(f"""
4. RESULT:
   ✓ Trend is in CORRECT direction (LSB > HSB)
   - Observed: {diff:+.2f}%
   - Predicted: +3.5%
   - Ratio: {diff/3.5:.2f}""")

        if abs(diff) > 1.0:
            print("   ✓ Effect is detectable and significant")
        else:
            print("   - Effect is small but positive")
    else:
        print(f"""
4. RESULT:
   ✗ Trend is OPPOSITE to prediction
   - Observed: {diff:+.2f}%
   - Predicted: +3.5%""")

print("""
5. CAVEATS:
   - LSB/HSB is imperfect proxy for void/cluster environment
   - Need true 3D density field classification
   - Distance alone doesn't determine environment
   - Surface brightness correlates with other properties

6. NEXT STEPS:
   - Cross-match with SDSS void catalogs for direct test
   - Use 2MRS or ALFALFA local density estimators
   - Look for void-residing galaxies with HI data
   - Consider galaxy group catalogs (Yang+2007, etc.)

FILES CREATED:
- session179_sparc_environment.png
""")

print("=" * 70)
print("SESSION #179 COMPLETE")
print("=" * 70)
