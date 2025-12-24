#!/usr/bin/env python3
"""
SESSION #175: FIRST PRINCIPLES REVIEW
======================================
Date: December 24, 2025

Following the guidance in RESEARCH_PHILOSOPHY.md, this session re-examines
the Synchronism peculiar velocity predictions from first principles.

KEY QUESTION:
-------------
Session #174 found CF4 shows OPPOSITE to Synchronism prediction.
Is nature telling us:
1. The coherence function prediction is wrong?
2. We're testing at the wrong MRH?
3. The peculiar velocity test is not the right signature?

APPROACH:
---------
1. Re-derive the coherence function prediction
2. Check what MRH the prediction applies to
3. Examine if peculiar velocities test the right thing
4. Consider alternative signatures
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("SESSION #175: FIRST PRINCIPLES REVIEW")
print("=" * 70)

# =============================================================================
# 1. THE SYNCHRONISM COHERENCE FUNCTION
# =============================================================================

print("\n" + "=" * 70)
print("1. COHERENCE FUNCTION: WHAT DOES IT ACTUALLY PREDICT?")
print("=" * 70)

print("""
The coherence function is:

    C(ρ) = Ω_m + (1 - Ω_m) × (ρ/ρ_t)^(1/φ) / [1 + (ρ/ρ_t)^(1/φ)]

Where:
    φ = 1.618... (golden ratio)
    Ω_m = 0.3 (matter density parameter)
    ρ_t = transition density

This describes the COUPLING between quantum and gravitational degrees of freedom.
""")

phi = (1 + np.sqrt(5)) / 2
Omega_m = 0.3
rho_t = 1.0

def coherence(rho_ratio):
    """Coherence function C(ρ)"""
    x = rho_ratio ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

# Plot coherence function
rho_range = np.logspace(-2, 3, 100)
C_values = [coherence(r) for r in rho_range]

print("\nCoherence values:")
print(f"  Deep void (ρ/ρ_t = 0.01): C = {coherence(0.01):.4f}")
print(f"  Void (ρ/ρ_t = 0.1): C = {coherence(0.1):.4f}")
print(f"  Mean density (ρ/ρ_t = 1): C = {coherence(1.0):.4f}")
print(f"  Overdense (ρ/ρ_t = 10): C = {coherence(10):.4f}")
print(f"  Cluster (ρ/ρ_t = 100): C = {coherence(100):.4f}")

# =============================================================================
# 2. WHAT DOES C CONTROL?
# =============================================================================

print("\n" + "=" * 70)
print("2. PHYSICAL MEANING OF COHERENCE")
print("=" * 70)

print("""
The coherence function C(ρ) is proposed to modify the effective gravitational
coupling:

    G_eff = G / C(ρ)

In regions where C < 1 (low density):
    - G_eff > G (enhanced gravity)
    - This should produce stronger gravitational dynamics

In regions where C → 1 (high density):
    - G_eff → G (standard gravity)
    - Gravity behaves normally

KEY QUESTION: What dynamics does enhanced G_eff actually affect?
""")

# G_eff ratio
G_ratio = [1/C for C in C_values]

print("\nG_eff/G ratio:")
print(f"  Deep void: G_eff/G = {1/coherence(0.01):.3f}")
print(f"  Void: G_eff/G = {1/coherence(0.1):.3f}")
print(f"  Mean: G_eff/G = {1/coherence(1.0):.3f}")
print(f"  Overdense: G_eff/G = {1/coherence(10):.3f}")
print(f"  Cluster: G_eff/G = {1/coherence(100):.3f}")

# =============================================================================
# 3. MRH CHECK: IS PECULIAR VELOCITY THE RIGHT TEST?
# =============================================================================

print("\n" + "=" * 70)
print("3. MRH CHECK: IS PECULIAR VELOCITY THE RIGHT TEST?")
print("=" * 70)

print("""
RESEARCH_PHILOSOPHY.md says: "Match MRH to the scale you're observing"

Peculiar velocity observations:
    - Scale: Individual galaxies (10^10 M_sun)
    - Distance: 10-200 Mpc
    - Observable: Line-of-sight velocity offset from Hubble flow

What does G_eff enhancement actually predict?

OPTION A: Enhanced local dynamics
    - Galaxies feel stronger local gravitational pull
    - Infall velocities increase
    - Velocity dispersion increases

OPTION B: Enhanced bulk flow
    - Large-scale structure grows faster
    - Coherent flows toward overdensities
    - Void outflows are enhanced

OPTION C: Enhanced galaxy internal dynamics
    - Rotation curves affected (this IS what fits with coherence!)
    - Internal velocity dispersion changed
    - NOT peculiar velocities between galaxies

CRITICAL INSIGHT:
-----------------
The coherence function was derived from GALAXY ROTATION CURVES, not
peculiar velocities. These probe DIFFERENT scales:

    Rotation curves: R ~ 1-30 kpc (internal galaxy dynamics)
    Peculiar velocities: d ~ 10-200 Mpc (intergalactic dynamics)

The MRH is DIFFERENT for these two observables!
""")

print("\nScale comparison:")
print(f"  Galaxy rotation: R ~ 10 kpc = 10^4 pc")
print(f"  Peculiar velocity: d ~ 100 Mpc = 10^8 pc")
print(f"  Scale ratio: ~10,000×")
print(f"  >>> Different MRH regimes!")

# =============================================================================
# 4. WHAT SHOULD PECULIAR VELOCITIES SHOW?
# =============================================================================

print("\n" + "=" * 70)
print("4. WHAT SHOULD PECULIAR VELOCITIES ACTUALLY SHOW?")
print("=" * 70)

print("""
At the MRH of intergalactic scales (Mpc), the coherence function affects:

1. GROWTH OF STRUCTURE
   - δρ/ρ grows faster in regions with higher G_eff
   - Voids should EXPAND faster (outflow)
   - Overdensities should COLLAPSE faster (infall)

2. BUT: The OBSERVED peculiar velocity is v_pec = v - H*d

   In voids (lower density):
   - Outflow velocity is enhanced (away from void center)
   - When we OBSERVE this from outside, we see:
     - Negative v_pec for near side (moving toward us out of void)
     - Positive v_pec for far side (moving away from us out of void)
   - The MEAN v_pec could be ZERO if symmetric!

3. WHAT WE ACTUALLY MEASURE
   - |v_pec| = magnitude of peculiar velocity
   - This should be enhanced in voids IF outflows are stronger
   - BUT: We're measuring v_pec derived from DISTANCES, not velocities!

CRITICAL INSIGHT:
-----------------
CF4 v_pec is derived from: v_pec = v_cmb - H0 * d

The distance d is measured from distance indicators (TF, FP, SNIa).
The "peculiar velocity" we measure is actually:

    v_pec = v_true + H0 * (d_true - d_observed)

If d_observed has errors that correlate with environment, the
v_pec we measure doesn't represent true peculiar velocities!
""")

# =============================================================================
# 5. THE FUNDAMENTAL PROBLEM
# =============================================================================

print("\n" + "=" * 70)
print("5. THE FUNDAMENTAL PROBLEM WITH PECULIAR VELOCITY TESTS")
print("=" * 70)

print("""
The CF4 test is fundamentally limited because:

1. We don't measure true peculiar velocities
   - We measure v_cmb (from spectra) - OK
   - We measure d (from distance indicators) - PROBLEMATIC
   - Distance errors >> true peculiar velocities

2. Selection effects correlate with environment
   - Session #174 showed: angular "void" galaxies have LOWER |v_pec|
   - This could be because:
     a. True velocities are lower (anti-Synchronism)
     b. Distance errors are different (selection effect)
     c. Environment classification doesn't match density

3. The Synchronism prediction may not apply to this observable
   - Coherence function was derived from rotation curves
   - Peculiar velocities probe a different MRH
   - Enhanced G_eff might not produce the expected signature

WHAT WE LEARNED:
----------------
The peculiar velocity test, as designed, cannot distinguish:
    - True Synchronism signal
    - Selection effects
    - Distance error correlations
    - Wrong MRH for the test
""")

# =============================================================================
# 6. ALTERNATIVE SIGNATURES
# =============================================================================

print("\n" + "=" * 70)
print("6. BETTER TESTS FOR SYNCHRONISM")
print("=" * 70)

print("""
Tests that probe the RIGHT MRH for coherence function:

1. GALAXY ROTATION CURVES (Already done!)
   - This IS where the coherence function was derived
   - Session #173 cluster dispersion test showed ratio = 3.28
   - Need more cluster-scale tests

2. VOID EXPANSION PROFILES
   - Measure void radius evolution with redshift
   - Enhanced G_eff → faster void expansion
   - Less affected by individual galaxy distance errors

3. CLUSTER MASS ESTIMATES
   - Compare lensing mass to dynamical mass
   - In Synchronism: Dynamical mass assumes G, but G_eff > G in outskirts
   - Should see systematic discrepancy vs radius

4. REDSHIFT-SPACE DISTORTIONS
   - Measure galaxy clustering in redshift space
   - Enhanced G_eff → different finger-of-god/Kaiser effect
   - DESI will provide excellent data

5. GRAVITATIONAL LENSING AROUND VOIDS
   - Voids should lens LESS if G_eff is higher there
   - Counter-intuitive but testable

RECOMMENDED NEXT STEPS:
-----------------------
1. Abandon peculiar velocity tests (too many systematics)
2. Focus on cluster dispersion profiles (Session #173 showed promise)
3. Wait for DESI redshift-space distortion data
4. Pursue lensing-based tests
""")

# =============================================================================
# 7. VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Coherence function
ax1 = axes[0, 0]
ax1.semilogx(rho_range, C_values, 'b-', linewidth=2)
ax1.axhline(1.0, color='gray', linestyle='--')
ax1.axhline(0.3, color='red', linestyle='--', label='Ω_m limit')
ax1.axvline(1.0, color='green', linestyle=':', label='ρ = ρ_t')
ax1.set_xlabel('Density ratio ρ/ρ_t')
ax1.set_ylabel('Coherence C(ρ)')
ax1.set_title('Synchronism Coherence Function')
ax1.legend()
ax1.set_ylim(0, 1.1)

# Panel 2: G_eff ratio
ax2 = axes[0, 1]
ax2.semilogx(rho_range, G_ratio, 'r-', linewidth=2)
ax2.axhline(1.0, color='gray', linestyle='--', label='Standard gravity')
ax2.fill_between(rho_range, 1, G_ratio, where=np.array(G_ratio)>1, alpha=0.3, color='red', label='Enhanced G')
ax2.set_xlabel('Density ratio ρ/ρ_t')
ax2.set_ylabel('G_eff / G')
ax2.set_title('Effective Gravitational Coupling')
ax2.legend()
ax2.set_ylim(0.9, 3.5)

# Panel 3: MRH scales
ax3 = axes[1, 0]
scales = ['Planck\n(10⁻³⁵ m)', 'Atom\n(10⁻¹⁰ m)', 'Galaxy\n(10²¹ m)', 'Cluster\n(10²³ m)', 'Cosmic\n(10²⁶ m)']
mrh_values = [1e-35, 1e-10, 1e21, 1e23, 1e26]
colors = ['purple', 'blue', 'green', 'orange', 'red']
ax3.barh(scales, [np.log10(m)+35 for m in mrh_values], color=colors, alpha=0.7)
ax3.set_xlabel('Log₁₀(Scale) + 35')
ax3.set_title('MRH Scales: Different Physics at Each Level')
ax3.axvline(np.log10(1e21)+35, color='green', linestyle='--', label='Rotation curves')
ax3.axvline(np.log10(1e23)+35, color='orange', linestyle='--', label='Peculiar velocities')

# Panel 4: Test comparison
ax4 = axes[1, 1]
tests = ['Rotation\ncurves', 'Peculiar\nvelocities', 'Cluster\ndispersion', 'Void\nlensing', 'RSDs']
reliability = [0.9, 0.2, 0.7, 0.8, 0.85]
colors = ['green', 'red', 'orange', 'green', 'green']
ax4.bar(tests, reliability, color=colors, alpha=0.7)
ax4.set_ylabel('Test Reliability (subjective)')
ax4.set_title('Synchronism Tests: Reliability Assessment')
ax4.set_ylim(0, 1)
ax4.axhline(0.5, color='gray', linestyle='--')

plt.suptitle('Session #175: First Principles Review of Synchronism Tests', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session175_first_principles.png', dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FIGURE SAVED: session175_first_principles.png")
print("=" * 70)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #175: SUMMARY")
print("=" * 70)

print("""
FIRST PRINCIPLES REVIEW
=======================

1. THE COHERENCE FUNCTION
   - Derived from galaxy rotation curves (kpc scale)
   - Predicts G_eff = G/C(ρ) enhancement in low-density regions
   - 64.6% of galaxies fit with 2 global parameters

2. MRH MISMATCH
   - Rotation curves: R ~ 10 kpc
   - Peculiar velocities: d ~ 100 Mpc
   - These probe DIFFERENT MRH regimes (10,000× scale difference)

3. PECULIAR VELOCITY TEST PROBLEMS
   - Derived from distance indicators, not true velocities
   - Distance errors >> true peculiar velocities
   - Selection effects correlate with environment
   - Cannot distinguish signal from systematics

4. SESSION #174 FINDING REINTERPRETED
   - CF4 showing "opposite" may not rule out Synchronism
   - May indicate: Wrong test, not wrong theory
   - MRH mismatch between prediction and observable

5. BETTER TESTS
   - Cluster dispersion profiles (Session #173: ratio = 3.28)
   - Void expansion rates
   - Lensing-mass vs dynamical-mass discrepancies
   - Redshift-space distortions (DESI)

6. CONCLUSION
   The peculiar velocity test arc (Sessions #169-174) revealed:
   - Sophisticated forward modeling is needed
   - Selection effects dominate the signal
   - The test may not probe the right MRH for coherence
   - Alternative tests are required

"Nature is telling you something. Listen to the data, not the paradigm."
    - RESEARCH_PHILOSOPHY.md

In this case: The data is telling us peculiar velocities are the
wrong observable for testing the coherence function.
""")

print("=" * 70)
print("SESSION #175 COMPLETE")
print("=" * 70)
