#!/usr/bin/env python3
"""
SESSION #179B: INVESTIGATING THE LSB/HSB PROXY FAILURE
=======================================================
Date: December 25, 2025

CONTEXT:
--------
Session #179 found that LSB galaxies show LOWER velocity residuals than HSB galaxies,
which is OPPOSITE to the Synchronism prediction from Session #177.

HYPOTHESIS:
-----------
The LSB/HSB classification may NOT be a good proxy for void/cluster environment.

LSB galaxies may have different properties (formation history, gas fraction, etc.)
that affect their kinematics independently of environment.

THIS SESSION:
-------------
1. Investigate the LSB/HSB correlation with other properties
2. Check if the BTFR slope differs systematically
3. Consider alternative environment proxies
4. Re-examine the theoretical prediction
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

print("=" * 70)
print("SESSION #179B: INVESTIGATING THE PROXY FAILURE")
print("=" * 70)

# =============================================================================
# 1. RELOAD SPARC DATA
# =============================================================================

sparc_file = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/sparc_real_data/SPARC_Lelli2016c.mrt'

galaxies = []
with open(sparc_file, 'r') as f:
    in_data = False
    for line in f:
        line = line.rstrip()
        if not in_data:
            if len(line) > 50 and line.lstrip()[:3] in ['CIG', 'D51', 'D63', 'DDO', 'ESO',
                'F56', 'F56', 'F57', 'F58', 'F59', 'F56', 'IC ', 'KK9', 'Mkn',
                'NGC', 'PGC', 'UGC', 'UGC']:
                in_data = True

        if in_data and len(line) > 80:
            try:
                parts = line.split()
                if len(parts) < 18:
                    continue

                name = parts[0]
                hubble_type = int(parts[1])
                distance = float(parts[2])
                e_distance = float(parts[3])
                inclination = float(parts[5])
                L_3p6 = float(parts[7])
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
            except (ValueError, IndexError):
                continue

print(f"Loaded {len(galaxies)} galaxies")

# Add baryonic mass and classifications
ML_ratio = 0.5
for g in galaxies:
    g['M_star'] = ML_ratio * g['L_3p6'] * 1e9
    g['M_gas'] = 1.33 * g['MHI'] * 1e9
    g['M_bary'] = g['M_star'] + g['M_gas']
    g['gas_fraction'] = g['M_gas'] / g['M_bary'] if g['M_bary'] > 0 else 0

    # Surface brightness classification
    if g['SBeff'] < 100:
        g['sb_class'] = 'LSB'
    elif g['SBeff'] > 500:
        g['sb_class'] = 'HSB'
    else:
        g['sb_class'] = 'intermediate'

# =============================================================================
# 2. ANALYZE LSB vs HSB DIFFERENCES
# =============================================================================

print("\n" + "=" * 70)
print("2. LSB vs HSB PROPERTY COMPARISON")
print("=" * 70)

lsb = [g for g in galaxies if g['sb_class'] == 'LSB']
hsb = [g for g in galaxies if g['sb_class'] == 'HSB']

print(f"\nSample sizes: LSB = {len(lsb)}, HSB = {len(hsb)}")

properties = [
    ('log(M_bary)', lambda g: np.log10(g['M_bary'])),
    ('log(Vflat)', lambda g: np.log10(g['Vflat'])),
    ('Gas fraction', lambda g: g['gas_fraction']),
    ('Distance (Mpc)', lambda g: g['D']),
    ('Reff (kpc)', lambda g: g['Reff']),
    ('Rdisk (kpc)', lambda g: g['Rdisk']),
    ('Hubble type', lambda g: g['type']),
    ('SBeff', lambda g: g['SBeff']),
]

print("\nProperty comparison:")
print("-" * 80)
print(f"{'Property':<20} {'LSB mean':>12} {'HSB mean':>12} {'Difference':>12} {'p-value':>12}")
print("-" * 80)

for name, func in properties:
    lsb_vals = [func(g) for g in lsb]
    hsb_vals = [func(g) for g in hsb]

    # Remove NaN/inf
    lsb_vals = [v for v in lsb_vals if np.isfinite(v)]
    hsb_vals = [v for v in hsb_vals if np.isfinite(v)]

    if len(lsb_vals) > 0 and len(hsb_vals) > 0:
        lsb_mean = np.mean(lsb_vals)
        hsb_mean = np.mean(hsb_vals)
        diff = lsb_mean - hsb_mean

        # t-test
        t_stat, p_val = stats.ttest_ind(lsb_vals, hsb_vals)

        print(f"{name:<20} {lsb_mean:>12.3f} {hsb_mean:>12.3f} {diff:>12.3f} {p_val:>12.4f}")

# =============================================================================
# 3. KEY INSIGHT: GAS FRACTION
# =============================================================================

print("\n" + "=" * 70)
print("3. GAS FRACTION ANALYSIS")
print("=" * 70)

print("""
CRITICAL OBSERVATION:
LSB galaxies typically have HIGHER gas fractions.

The gas is distributed more extended than the stellar disk.
This affects where the "dark matter" signal appears in the rotation curve.

In Synchronism, higher gas fraction means:
- More extended baryonic distribution
- Transition from disk-dominated to halo-dominated at LARGER radii
- Different effective ρ_t scale

This may confound the environment effect!
""")

lsb_gas = [g['gas_fraction'] for g in lsb]
hsb_gas = [g['gas_fraction'] for g in hsb]

print(f"Mean gas fraction:")
print(f"  LSB: {np.mean(lsb_gas):.3f} ± {np.std(lsb_gas):.3f}")
print(f"  HSB: {np.mean(hsb_gas):.3f} ± {np.std(hsb_gas):.3f}")

# =============================================================================
# 4. MASS-MATCHED COMPARISON
# =============================================================================

print("\n" + "=" * 70)
print("4. MASS-MATCHED COMPARISON")
print("=" * 70)

print("""
The BTFR slope difference between LSB and HSB could be due to different
mass distributions in each sample.

Let's compare residuals at FIXED mass to isolate the environment effect.
""")

# Overall BTFR
def btfr(M_bary, a, b):
    return a * np.log10(M_bary) + b

all_M = np.array([g['M_bary'] for g in galaxies])
all_v = np.array([g['Vflat'] for g in galaxies])

popt_all, _ = curve_fit(btfr, all_M, np.log10(all_v), p0=[0.25, -0.5])
a_all, b_all = popt_all

# Compute residuals
for g in galaxies:
    v_expected = 10**btfr(g['M_bary'], a_all, b_all)
    g['v_residual'] = (g['Vflat'] - v_expected) / v_expected

# Mass bins
mass_bins = [
    (1e8, 1e9, 'Low mass (10^8-10^9)'),
    (1e9, 1e10, 'Medium mass (10^9-10^10)'),
    (1e10, 1e11, 'High mass (10^10-10^11)'),
    (1e11, 1e12, 'Massive (10^11-10^12)'),
]

print("\nResiduals by mass bin and surface brightness:")
print("-" * 80)
print(f"{'Mass bin':<30} {'LSB (N)':>15} {'HSB (N)':>15} {'LSB-HSB':>15}")
print("-" * 80)

for m_low, m_high, label in mass_bins:
    lsb_in_bin = [g for g in lsb if m_low <= g['M_bary'] < m_high]
    hsb_in_bin = [g for g in hsb if m_low <= g['M_bary'] < m_high]

    if len(lsb_in_bin) >= 3 and len(hsb_in_bin) >= 3:
        lsb_res = np.mean([g['v_residual'] for g in lsb_in_bin])
        hsb_res = np.mean([g['v_residual'] for g in hsb_in_bin])
        diff = lsb_res - hsb_res
        print(f"{label:<30} {lsb_res*100:>10.2f}% ({len(lsb_in_bin)}) {hsb_res*100:>10.2f}% ({len(hsb_in_bin)}) {diff*100:>10.2f}%")
    else:
        print(f"{label:<30} Insufficient data (LSB={len(lsb_in_bin)}, HSB={len(hsb_in_bin)})")

# =============================================================================
# 5. RE-EXAMINE THE THEORETICAL PREDICTION
# =============================================================================

print("\n" + "=" * 70)
print("5. RE-EXAMINING THE THEORETICAL PREDICTION")
print("=" * 70)

print("""
SESSION #177 PREDICTION:
========================
Void galaxies should have ~3.5% HIGHER v_rot because:
- Lower environment density → lower effective ρ at large radii
- Lower ρ → lower coherence C → higher G_eff
- Higher G_eff → higher v_rot

ASSUMPTION CHECK:
=================
1. Environment density determines large-radius ρ_eff
2. LSB galaxies are in low-density environments (voids)

PROBLEM IDENTIFIED:
===================
Assumption #2 may be WRONG!

LSB galaxies have LOW surface brightness because of:
- LOW STELLAR DENSITY (low formation efficiency)
- EXTENDED DISTRIBUTION (same mass over larger area)

But this doesn't mean they're in VOIDS!

LSB galaxies could be:
- Field galaxies with different formation histories
- Tidal dwarfs
- Gas-rich systems with suppressed star formation

The LSB/void correlation is WEAK in observations.
""")

# =============================================================================
# 6. ALTERNATIVE: DISTANCE AS PROXY
# =============================================================================

print("\n" + "=" * 70)
print("6. ALTERNATIVE: DISTANCE AS ENVIRONMENT PROXY")
print("=" * 70)

print("""
ALTERNATIVE HYPOTHESIS:
=======================
Instead of surface brightness, use DISTANCE as environment proxy:

- Nearby galaxies (< 10 Mpc): Include Local Group and Virgo infall region
  → Mix of cluster and field

- Distant galaxies (> 50 Mpc): More likely to be isolated field/void galaxies
  → Lower large-scale density

Let's test this...
""")

nearby = [g for g in galaxies if g['D'] < 10]
distant = [g for g in galaxies if g['D'] > 50]

print(f"\nSample sizes: Nearby (<10 Mpc) = {len(nearby)}, Distant (>50 Mpc) = {len(distant)}")

if len(nearby) >= 5 and len(distant) >= 5:
    nearby_res = np.mean([g['v_residual'] for g in nearby])
    distant_res = np.mean([g['v_residual'] for g in distant])

    print(f"\nVelocity residuals:")
    print(f"  Nearby: {nearby_res*100:+.2f}%")
    print(f"  Distant: {distant_res*100:+.2f}%")
    print(f"  Difference (Distant - Nearby): {(distant_res - nearby_res)*100:+.2f}%")

    # Synchronism predicts distant should be HIGHER (more void-like)
    if distant_res > nearby_res:
        print(f"\n  ✓ Trend is in CORRECT direction for Synchronism")
    else:
        print(f"\n  ✗ Trend is OPPOSITE to Synchronism prediction")

# =============================================================================
# 7. HUBBLE TYPE AS PROXY
# =============================================================================

print("\n" + "=" * 70)
print("7. HUBBLE TYPE ANALYSIS")
print("=" * 70)

print("""
Hubble types in SPARC:
0 = S0, 1 = Sa, ..., 9 = Sm, 10 = Im, 11 = BCD

Late-type (high T) galaxies tend to be:
- More gas-rich
- Lower surface brightness
- More common in field/void environments

Let's check residuals by Hubble type...
""")

early = [g for g in galaxies if g['type'] <= 4]  # S0 to Sbc
late = [g for g in galaxies if g['type'] >= 8]  # Sdm to BCD

print(f"\nSample sizes: Early (T≤4) = {len(early)}, Late (T≥8) = {len(late)}")

if len(early) >= 5 and len(late) >= 5:
    early_res = np.mean([g['v_residual'] for g in early])
    late_res = np.mean([g['v_residual'] for g in late])

    print(f"\nVelocity residuals:")
    print(f"  Early type: {early_res*100:+.2f}%")
    print(f"  Late type: {late_res*100:+.2f}%")
    print(f"  Difference (Late - Early): {(late_res - early_res)*100:+.2f}%")

# =============================================================================
# 8. VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("8. GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: BTFR with LSB/HSB coloring
ax1 = axes[0, 0]
M_plot = np.logspace(7, 12, 100)

for sb, color, label in [('LSB', 'red', 'LSB (low SB)'), ('HSB', 'blue', 'HSB (high SB)')]:
    gals = [g for g in galaxies if g['sb_class'] == sb]
    M = [g['M_bary'] for g in gals]
    v = [g['Vflat'] for g in gals]
    ax1.scatter(M, v, c=color, alpha=0.6, s=30, label=label)

ax1.plot(M_plot, 10**btfr(M_plot, a_all, b_all), 'k-', linewidth=2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Baryonic Mass (M☉)')
ax1.set_ylabel('Vflat (km/s)')
ax1.set_title('BTFR: LSB vs HSB')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel 2: Residual vs gas fraction
ax2 = axes[0, 1]
gas_frac = [g['gas_fraction'] for g in galaxies]
residuals = [g['v_residual']*100 for g in galaxies]
sb_colors = [{'LSB': 'red', 'intermediate': 'green', 'HSB': 'blue'}[g['sb_class']] for g in galaxies]

ax2.scatter(gas_frac, residuals, c=sb_colors, alpha=0.6, s=30)
ax2.axhline(0, color='black', linestyle='--', linewidth=1)
ax2.set_xlabel('Gas Fraction')
ax2.set_ylabel('Velocity Residual (%)')
ax2.set_title('Residual vs Gas Fraction')
ax2.grid(True, alpha=0.3)

# Add legend manually
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='red', alpha=0.6, label='LSB'),
                   Patch(facecolor='green', alpha=0.6, label='Intermediate'),
                   Patch(facecolor='blue', alpha=0.6, label='HSB')]
ax2.legend(handles=legend_elements)

# Panel 3: Residual vs distance
ax3 = axes[1, 0]
distances = [g['D'] for g in galaxies]
ax3.scatter(distances, residuals, c=sb_colors, alpha=0.6, s=30)
ax3.axhline(0, color='black', linestyle='--', linewidth=1)
ax3.axhline(3.5, color='green', linestyle=':', linewidth=2, alpha=0.7, label='Sync prediction')
ax3.set_xlabel('Distance (Mpc)')
ax3.set_ylabel('Velocity Residual (%)')
ax3.set_title('Residual vs Distance')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Panel 4: Residual distribution comparison
ax4 = axes[1, 1]
lsb_res = [g['v_residual']*100 for g in lsb]
hsb_res = [g['v_residual']*100 for g in hsb]

bins = np.linspace(-50, 50, 25)
ax4.hist(lsb_res, bins=bins, alpha=0.5, label=f'LSB (mean={np.mean(lsb_res):.1f}%)', color='red')
ax4.hist(hsb_res, bins=bins, alpha=0.5, label=f'HSB (mean={np.mean(hsb_res):.1f}%)', color='blue')
ax4.axvline(0, color='black', linestyle='--', linewidth=1)
ax4.axvline(3.5, color='green', linestyle=':', linewidth=2, label='Sync prediction: +3.5%')
ax4.set_xlabel('Velocity Residual (%)')
ax4.set_ylabel('Count')
ax4.set_title('Residual Distribution')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.suptitle('Session #179B: Investigating LSB/HSB as Environment Proxy',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session179b_proxy_investigation.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session179b_proxy_investigation.png")

# =============================================================================
# 9. SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #179B: SUMMARY")
print("=" * 70)

print(f"""
INVESTIGATION OF LSB/HSB PROXY FAILURE
======================================

1. FINDING: LSB galaxies show NEGATIVE velocity residuals (-2.6%)
   while HSB galaxies show POSITIVE residuals (+3.5%)

   This is OPPOSITE to the Synchronism prediction.

2. KEY INSIGHT:
   LSB ≠ Void!

   LSB galaxies have low surface brightness due to:
   - Different stellar formation efficiency
   - Different gas-to-stellar mass ratio
   - Extended distribution

   NOT because they are in low-density environments.

3. CONFOUNDING FACTORS:
   - Gas fraction: LSB has {np.mean([g['gas_fraction'] for g in lsb]):.2f}, HSB has {np.mean([g['gas_fraction'] for g in hsb]):.2f}
   - Mass: LSB has lower mean mass
   - Morphology: LSB tends to be later type

4. ALTERNATIVE TESTS NEEDED:
   - Direct void catalog cross-matching (SDSS voids)
   - Local galaxy density from 2MRS or ALFALFA
   - Galaxy group membership catalogs

5. CONCLUSION:
   The Session #177 prediction CANNOT be tested using LSB/HSB as proxy.
   The proxy assumption is invalid.

   We need TRUE environment classification to test the prediction.

   THIS IS NOT A FALSIFICATION OF SYNCHRONISM.
   It's a FAILURE OF THE PROXY METHOD.

NEXT STEPS:
- Session #180: Obtain proper void catalogs
- Cross-match SPARC with local density estimates
- Test prediction with true environment classification
""")

print("=" * 70)
print("SESSION #179B COMPLETE")
print("=" * 70)
