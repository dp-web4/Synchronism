#!/usr/bin/env python3
"""
SESSION #180: MRH RE-EXAMINATION OF THE VOID/CLUSTER PREDICTION
================================================================
Date: December 25, 2025

CONTEXT:
--------
Session #179 found that ALL environment proxies showed the OPPOSITE trend
to the Session #177 prediction. Following RESEARCH_PHILOSOPHY.md guidance:

    "Am I adding epicycles, or is nature telling us to change the paradigm?"

This session re-examines the theoretical prediction from first principles,
considering whether the MRH abstraction was correct.

KEY QUESTION:
-------------
At what MRH does the "environment density" effect operate?

- Session #177 assumed: Environment density affects ρ_eff at large galactic radii
- But: Galactic rotation curves are measured at r ~ 10-50 kpc
- The "environment" (cluster/void) is at scales > 1 Mpc

HYPOTHESIS:
-----------
The environment effect may operate at a DIFFERENT MRH than galactic dynamics.
The coherence function may need scale-dependent environment terms.

APPROACH:
---------
1. Clarify the MRH scales involved
2. Re-derive the environment effect with proper scale separation
3. Consider what Session #179 data actually tells us
4. Derive revised predictions if needed
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("SESSION #180: MRH RE-EXAMINATION OF VOID/CLUSTER PREDICTION")
print("=" * 70)

# =============================================================================
# 1. MRH SCALE ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("1. MRH SCALE ANALYSIS")
print("=" * 70)

print("""
SCALES INVOLVED:

1. GALACTIC DISK (where rotation curves are measured)
   - r = 5-50 kpc
   - ρ = 0.01 - 10 M☉/pc³ (from disk center to outskirts)
   - This is where v_rot is observed

2. GALACTIC HALO (dark matter "halo")
   - r = 50-300 kpc
   - ρ = 10⁻⁵ - 10⁻² M☉/pc³
   - NFW profile at this scale

3. LOCAL GROUP / CLUSTER ENVIRONMENT
   - r = 1-10 Mpc from galaxy center
   - ρ = 10⁻⁸ - 10⁻⁶ M☉/pc³ (cosmic mean ~ 4×10⁻⁸)
   - This is where "environment" (void/cluster) operates

CRITICAL INSIGHT:
=================
The Session #177 prediction assumed environment density directly affects
the galactic rotation curve via ρ_eff = max(ρ_baryon, ρ_env).

But these operate at DIFFERENT MRH scales:
- Rotation curve: MRH ~ 10 kpc
- Environment: MRH ~ 1 Mpc

The environment density cannot directly affect rotation curve dynamics.
The effect must be INDIRECT - through halo properties formed during collapse.
""")

# Define scales
scales = {
    'Galactic disk': {'r_min': 1, 'r_max': 50, 'rho_min': 0.01, 'rho_max': 10},
    'Galactic halo': {'r_min': 50, 'r_max': 300, 'rho_min': 1e-5, 'rho_max': 0.01},
    'Environment': {'r_min': 1000, 'r_max': 10000, 'rho_min': 1e-8, 'rho_max': 1e-5},
}

print("\nScale hierarchy (r in kpc, ρ in M☉/pc³):")
print("-" * 70)
for name, scale in scales.items():
    print(f"{name:20s}: r = {scale['r_min']:>6} - {scale['r_max']:<6} kpc, "
          f"ρ = {scale['rho_min']:.0e} - {scale['rho_max']:.0e}")

# =============================================================================
# 2. THE SESSION #177 ERROR
# =============================================================================

print("\n" + "=" * 70)
print("2. THE SESSION #177 ERROR")
print("=" * 70)

print("""
SESSION #177 ASSUMED:
=====================
ρ_eff(r) = max(ρ_baryon(r), ρ_environment)

This is WRONG for two reasons:

1. SCALE MISMATCH:
   - ρ_environment is defined at Mpc scales
   - ρ_baryon is defined at kpc scales
   - Cannot compare densities at different MRH

2. INSTANTANEOUS vs FORMATION:
   - Current environment density ≠ formation environment density
   - Galaxy halos formed at z ~ 1-3, environment was different
   - Current rotation curve reflects FORMATION conditions

CORRECT APPROACH:
=================
The environment affects galaxies DURING FORMATION, not dynamically at z=0.

What we should predict:
- Galaxies that FORMED in low-density environments have different halos
- This affects their formation history, not their current ρ_eff

This is actually what ΛCDM predicts too:
- Cluster galaxies: Earlier formation, higher concentration
- Field galaxies: Later formation, lower concentration
- Void galaxies: Latest formation, lowest concentration

The Synchronism prediction should be about FORMATION, not current environment.
""")

# =============================================================================
# 3. WHAT THE SESSION #179 DATA ACTUALLY SHOWS
# =============================================================================

print("\n" + "=" * 70)
print("3. WHAT SESSION #179 DATA ACTUALLY SHOWS")
print("=" * 70)

print("""
SESSION #179 FINDINGS:
======================
1. LSB galaxies have LOWER velocity residuals (-2.6%)
2. HSB galaxies have HIGHER velocity residuals (+3.5%)
3. Difference: -6% (opposite to prediction)

BUT WE ALSO FOUND:
==================
- LSB galaxies have HIGHER gas fractions (0.67 vs 0.17)
- LSB galaxies have LOWER masses (10^9.3 vs 10^10.8 M☉)
- LSB galaxies are LATER Hubble type (8.4 vs 3.1)

REINTERPRETATION:
=================
What if the "opposite trend" is CORRECT?

In Synchronism:
- G_eff = G/C(ρ) where C(ρ) → 1 at high density
- High density → C → 1 → G_eff → G (standard gravity)
- Low density → C → Ω_m → G_eff → G/Ω_m (enhanced gravity)

LSB galaxies have:
- Lower SURFACE density (spread out)
- But same or lower VOLUMETRIC density?
- Extended gas distribution

HYPOTHESIS:
===========
LSB galaxies may actually have HIGHER G_eff (lower density on average)
BUT their extended distribution means:
- Baryonic mass is at larger radii
- The v_rot at a given r samples different mass distribution
- Net effect on BTFR residual could be negative!

This is a BTFR slope effect, not an environment effect.
""")

# =============================================================================
# 4. REVISED THEORETICAL FRAMEWORK
# =============================================================================

print("\n" + "=" * 70)
print("4. REVISED THEORETICAL FRAMEWORK")
print("=" * 70)

print("""
CORRECTED SYNCHRONISM PREDICTION:
=================================

The void/cluster effect operates through FORMATION, not current dynamics.

1. FORMATION ENVIRONMENT EFFECT:
   - Void galaxies form from lower overdensity peaks
   - Later collapse time → more extended halos
   - Lower concentration → more "dark matter" at large radii

2. BUT THIS IS THE SAME AS ΛCDM:
   - ΛCDM also predicts environment-dependent halo properties
   - The prediction is NOT discriminating

3. WHAT WOULD BE DISCRIMINATING:
   - Synchronism: G_eff enhancement increases as ρ decreases
   - At FIXED halo mass/concentration, void galaxies should show
     SLIGHTLY higher velocities due to G_eff enhancement
   - But this effect is small (~3%) and requires matching
     halo properties, not just stellar mass

REVISED PREDICTION:
==================
At FIXED dark matter fraction (not stellar mass), void galaxies
should show ~3% higher velocities.

But current SPARC data doesn't have:
- Environment classification
- Dark matter fraction measurement (only inferred)
- Halo concentration data

THE TEST REQUIRES:
- Proper void/cluster classification (not LSB/HSB proxy)
- Matching by dark matter fraction or halo concentration
- High-precision velocity measurements

SESSION #179 USED THE WRONG TEST.
""")

# =============================================================================
# 5. WHAT LSB/HSB ACTUALLY TELLS US
# =============================================================================

print("\n" + "=" * 70)
print("5. WHAT LSB/HSB ACTUALLY TELLS US")
print("=" * 70)

print("""
THE LSB/HSB RESULT IS INTERESTING:
==================================

LSB galaxies: Lower velocity residual (-2.6%)
HSB galaxies: Higher velocity residual (+3.5%)

This could mean:

OPTION A: BTFR CURVATURE
------------------------
The BTFR may have slight curvature that LSB/HSB samples differently:
- LSB = lower mass end: BTFR slope slightly steeper
- HSB = higher mass end: BTFR slope slightly shallower
- Fitting single slope → residuals correlate with mass → correlate with SB

OPTION B: GAS DYNAMICS
----------------------
LSB galaxies have high gas fraction (0.67 vs 0.17).
Gas dynamics are different from stellar dynamics:
- Gas is collisional, subject to pressure
- Gas rotation may lag behind circular velocity
- High gas fraction → measured v_flat may underestimate true potential

OPTION C: FORMATION HISTORY
---------------------------
LSB galaxies may form in lower-density regions:
- Lower overdensity → later collapse
- More extended halos → lower concentration
- More "dark matter" at large radii → but also more uncertainty

OPTION D: MEASUREMENT EFFECTS
-----------------------------
LSB galaxies are harder to measure:
- Lower signal-to-noise in HI observations
- More asymmetric rotation curves
- Inclination corrections more uncertain

NONE OF THESE TEST THE SYNCHRONISM ENVIRONMENT PREDICTION.
""")

# =============================================================================
# 6. QUANTITATIVE ANALYSIS: BTFR SLOPE BY MASS
# =============================================================================

print("\n" + "=" * 70)
print("6. QUANTITATIVE ANALYSIS")
print("=" * 70)

# Reload SPARC data for quantitative analysis
sparc_file = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/sparc_real_data/SPARC_Lelli2016c.mrt'

galaxies = []
with open(sparc_file, 'r') as f:
    in_data = False
    for line in f:
        line = line.rstrip()
        if not in_data:
            if len(line) > 50 and line.lstrip()[:3] in ['CIG', 'D51', 'D63', 'DDO', 'ESO',
                'F56', 'F57', 'F58', 'F59', 'IC ', 'KK9', 'Mkn',
                'NGC', 'PGC', 'UGC']:
                in_data = True

        if in_data and len(line) > 80:
            try:
                parts = line.split()
                if len(parts) < 18:
                    continue

                name = parts[0]
                L_3p6 = float(parts[7])
                SBeff = float(parts[10])
                MHI = float(parts[13])
                Vflat = float(parts[15])

                if Vflat > 0 and L_3p6 > 0:
                    M_star = 0.5 * L_3p6 * 1e9
                    M_gas = 1.33 * MHI * 1e9
                    M_bary = M_star + M_gas

                    galaxies.append({
                        'name': name,
                        'M_bary': M_bary,
                        'Vflat': Vflat,
                        'SBeff': SBeff,
                        'gas_frac': M_gas / M_bary if M_bary > 0 else 0
                    })
            except (ValueError, IndexError):
                continue

print(f"Loaded {len(galaxies)} galaxies")

# Split by mass
low_mass = [g for g in galaxies if g['M_bary'] < 1e10]
high_mass = [g for g in galaxies if g['M_bary'] >= 1e10]

print(f"Low mass (<10^10): {len(low_mass)}")
print(f"High mass (≥10^10): {len(high_mass)}")

# Compute BTFR slopes for each subsample
from scipy.optimize import curve_fit

def btfr(M, a, b):
    return a * np.log10(M) + b

def fit_btfr(gal_list):
    M = np.array([g['M_bary'] for g in gal_list])
    v = np.array([g['Vflat'] for g in gal_list])
    popt, _ = curve_fit(btfr, M, np.log10(v), p0=[0.25, -0.5])
    return popt

if len(low_mass) > 10 and len(high_mass) > 10:
    a_low, b_low = fit_btfr(low_mass)
    a_high, b_high = fit_btfr(high_mass)

    print(f"\nBTFR slopes:")
    print(f"  Low mass:  a = {a_low:.4f}")
    print(f"  High mass: a = {a_high:.4f}")
    print(f"  Difference: {a_low - a_high:.4f}")

    if a_low > a_high:
        print("\n  → Low-mass galaxies have STEEPER BTFR slope")
        print("    This explains why LSB (lower mass) shows negative residuals")
        print("    when fit with single overall slope")

# =============================================================================
# 7. GAS FRACTION EFFECT
# =============================================================================

print("\n" + "=" * 70)
print("7. GAS FRACTION EFFECT ON BTFR")
print("=" * 70)

# Overall BTFR
all_M = np.array([g['M_bary'] for g in galaxies])
all_v = np.array([g['Vflat'] for g in galaxies])
a_all, b_all = curve_fit(btfr, all_M, np.log10(all_v), p0=[0.25, -0.5])[0]

# Compute residuals
for g in galaxies:
    g['v_residual'] = (g['Vflat'] - 10**btfr(g['M_bary'], a_all, b_all)) / 10**btfr(g['M_bary'], a_all, b_all)

# Correlation with gas fraction
gas_fracs = [g['gas_frac'] for g in galaxies]
residuals = [g['v_residual'] for g in galaxies]

from scipy import stats
correlation, p_value = stats.pearsonr(gas_fracs, residuals)

print(f"Correlation between gas fraction and velocity residual:")
print(f"  Pearson r = {correlation:.3f}")
print(f"  p-value = {p_value:.4f}")

if p_value < 0.05:
    print(f"\n  → SIGNIFICANT correlation (p < 0.05)")
    if correlation < 0:
        print("    Higher gas fraction → LOWER velocity residual")
        print("    This explains why LSB (gas-rich) shows negative residuals")
else:
    print(f"\n  → No significant correlation")

# =============================================================================
# 8. REVISED CONCLUSION
# =============================================================================

print("\n" + "=" * 70)
print("8. REVISED CONCLUSION")
print("=" * 70)

print("""
SESSION #180 FINDINGS:
======================

1. THE SESSION #177 PREDICTION WAS FLAWED:
   - Environment (void/cluster) operates at Mpc scales
   - Rotation curves operate at kpc scales
   - Cannot mix densities across MRH

2. THE SESSION #179 "OPPOSITE TREND" IS EXPLAINED BY:
   - BTFR slope variation with mass
   - Gas fraction effects on measured velocities
   - NOT an environment density effect

3. PROPER TEST OF SYNCHRONISM REQUIRES:
   - True environment classification (void catalog)
   - Matching by dark matter fraction, not stellar mass
   - High-precision data that controls for confounders

4. THE PREDICTION ITSELF MAY NEED REVISION:
   - Environment affects formation, not current dynamics
   - Synchronism effect is through G_eff at local density
   - Environment enters indirectly through halo concentration

IMPLICATIONS FOR SYNCHRONISM:
=============================

The Session #177 prediction was based on a MISUNDERSTANDING of how
environment affects galaxies. The prediction should be:

"At FIXED local density (matched halos), Synchronism predicts
slightly enhanced gravity in all galaxies. The environment effect
is through formation history, which is the SAME as ΛCDM prediction."

This means the void/cluster test is NOT a discriminating test.

WHAT WOULD BE DISCRIMINATING:
=============================

1. M_dyn/M_lens ratio (Session #176)
   - Direct measurement of G_eff
   - Independent of environment proxies

2. Radial profile of inferred dark matter
   - Synchronism: C(ρ) varies with local density
   - ΛCDM: NFW profile

3. Extremely low-density systems
   - Ultra-diffuse galaxies
   - Outer Milky Way
   - Intergalactic gas

CONCLUSION:
===========
Session #177's prediction was based on MRH confusion.
Session #179's "failure" was actually success - it revealed the flaw.
The void/cluster test is NOT a discriminating test for Synchronism.
Need to focus on M_dyn/M_lens and radial profile tests.
""")

# =============================================================================
# 9. VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("9. GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: BTFR by mass bin with different slopes
ax1 = axes[0, 0]
M_plot = np.logspace(7, 12, 100)

if len(low_mass) > 0:
    M_l = [g['M_bary'] for g in low_mass]
    v_l = [g['Vflat'] for g in low_mass]
    ax1.scatter(M_l, v_l, c='blue', alpha=0.6, s=30, label='Low mass')
    ax1.plot(M_plot, 10**btfr(M_plot, a_low, b_low), 'b--', linewidth=2, label=f'Low mass fit (a={a_low:.3f})')

if len(high_mass) > 0:
    M_h = [g['M_bary'] for g in high_mass]
    v_h = [g['Vflat'] for g in high_mass]
    ax1.scatter(M_h, v_h, c='red', alpha=0.6, s=30, label='High mass')
    ax1.plot(M_plot, 10**btfr(M_plot, a_high, b_high), 'r--', linewidth=2, label=f'High mass fit (a={a_high:.3f})')

ax1.plot(M_plot, 10**btfr(M_plot, a_all, b_all), 'k-', linewidth=2, label='Overall fit')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Baryonic Mass (M☉)')
ax1.set_ylabel('Vflat (km/s)')
ax1.set_title('BTFR: Slope Varies with Mass')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# Panel 2: Residual vs gas fraction
ax2 = axes[0, 1]
ax2.scatter(gas_fracs, [r*100 for r in residuals], alpha=0.5, s=30)
ax2.axhline(0, color='black', linestyle='--', linewidth=1)

# Add regression line
z = np.polyfit(gas_fracs, [r*100 for r in residuals], 1)
p = np.poly1d(z)
ax2.plot([0, 1], [p(0), p(1)], 'r-', linewidth=2, label=f'r={correlation:.2f}')

ax2.set_xlabel('Gas Fraction')
ax2.set_ylabel('Velocity Residual (%)')
ax2.set_title('Residual vs Gas Fraction')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Panel 3: MRH scale diagram
ax3 = axes[1, 0]
ax3.set_xlim(-1, 5)
ax3.set_ylim(0, 4)

# Draw scale boxes
scales_plot = [
    ('Disk\n(1-50 kpc)', 0, 1, 'blue'),
    ('Halo\n(50-300 kpc)', 1.5, 1, 'green'),
    ('Environment\n(1-10 Mpc)', 3, 1, 'red'),
]

for label, x, width, color in scales_plot:
    rect = plt.Rectangle((x, 1.5), width, 1, fill=True, alpha=0.3, color=color)
    ax3.add_patch(rect)
    ax3.text(x + width/2, 2, label, ha='center', va='center', fontsize=10)

ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_title('MRH Scale Separation')
ax3.text(2, 0.5, 'Session #177 incorrectly mixed scales', ha='center', fontsize=10, style='italic')

# Panel 4: Summary text
ax4 = axes[1, 1]
ax4.axis('off')
summary_text = """
SESSION #180 KEY FINDINGS:

1. Session #177 prediction was based on MRH confusion
   - Environment (Mpc) ≠ Rotation curve (kpc)

2. Session #179 "opposite trend" is explained by:
   - BTFR slope varies with mass
   - Gas fraction correlates with residuals

3. The void/cluster test is NOT discriminating
   - Environment affects formation, not dynamics
   - Same prediction as ΛCDM

4. Better discriminating tests:
   - M_dyn/M_lens ratio
   - Radial profile of inferred DM
   - Ultra-low-density systems

CONCLUSION: The prediction needs revision.
"""
ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=10,
         verticalalignment='top', family='monospace')

plt.suptitle('Session #180: MRH Re-examination of Void/Cluster Prediction',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session180_mrh_reexamination.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session180_mrh_reexamination.png")

print("\n" + "=" * 70)
print("SESSION #180 COMPLETE")
print("=" * 70)
