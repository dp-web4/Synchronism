#!/usr/bin/env python3
"""
Session #200: Abell 520 Analysis - The "Train Wreck" Cluster
============================================================

Abell 520 is a merging cluster that presents an interesting puzzle.
Initial studies (Mahdavi+ 2007) found a "dark core" - a lensing peak
that coincided with the hot gas rather than the galaxies.

This was initially seen as contradicting ΛCDM (dark matter should
follow galaxies, not gas, as in the Bullet Cluster).

Later work found more complex structure with multiple mass peaks.

How does Synchronism interpret Abell 520?

Date: December 30, 2025
Machine: CBP
"""

import numpy as np

print("="*70)
print("SESSION #200: ABELL 520 - THE 'TRAIN WRECK' CLUSTER")
print("="*70)

print("""
ABELL 520 OBSERVATIONS
======================

Basic Properties:
- Redshift: z = 0.199
- Total mass: ~1 × 10^15 M_sun
- Status: Major merger in progress
- Nickname: "The Train Wreck"

Key Observations (Mahdavi+ 2007, 2008):
---------------------------------------
1. Multiple substructures visible in X-ray and lensing
2. Central "dark core" - lensing peak with:
   - High inferred mass (~10^14 M_sun)
   - Low galaxy density
   - Coincident with hot X-ray gas
3. Two galaxy concentrations OFFSET from main lensing peaks

The Puzzle:
-----------
In the Bullet Cluster:
- DM (lensing) follows galaxies (collisionless)
- Gas (X-ray) is offset (collisional)
- This SUPPORTS CDM being collisionless

In Abell 520 (initially):
- Lensing peak at GAS location
- Galaxies offset from main mass
- This seems to CONTRADICT CDM!

If DM is collisionless, why isn't it with the galaxies?

RESOLVED? (Jee+ 2012, 2014):
----------------------------
Deeper HST lensing found MORE substructure:
- Multiple mass peaks, not just one "dark core"
- Some peaks DO correlate with galaxies
- The "dark core" may be real but interpretation changed

Current understanding:
- Complex merger geometry (possibly 3+ subclusters)
- Line-of-sight projection effects
- "Dark core" may be overlap region of multiple halos
""")

print("\n" + "="*70)
print("SYNCHRONISM INTERPRETATION OF ABELL 520")
print("="*70)

print("""
SYNCHRONISM FRAMEWORK:
======================

Recall from Session #196-197:
- Total mass = Baryons + Indifferent Patterns
- Indifferent patterns are collisionless (like CDM)
- They follow the galaxies, not the gas
- G_eff enhancement adds apparent mass

In a complex merger like Abell 520:

1. INDIFFERENT PATTERNS
   - Follow the galaxy subclusters (collisionless)
   - Each subcluster brings its indifferent mass
   - Complex geometry creates complex lensing

2. G_EFF ENHANCEMENT
   - a/a₀ varies across the cluster
   - Different regions have different G_eff
   - This adds position-dependent "extra mass"

3. PROJECTION EFFECTS
   - Multiple subclusters along line of sight
   - Lensing integrates through entire structure
   - Can create "phantom" mass concentrations

THE "DARK CORE" IN SYNCHRONISM:
-------------------------------

Possibility 1: Real indifferent mass concentration
- Some indifferent patterns DID stay with gas
- Self-interaction cross-section not exactly zero
- But this contradicts Bullet Cluster...

Possibility 2: Projection artifact
- Multiple mass peaks overlap in projection
- Creates apparent high-mass, low-galaxy region
- More likely given complex geometry

Possibility 3: G_eff geometry
- Region of low acceleration (a << a₀)
- High G_eff inflates apparent mass
- But lensing should see true mass...

Most likely: Combination of projection + complex structure
""")

# =============================================================================
# QUANTITATIVE ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("QUANTITATIVE ANALYSIS")
print("="*70)

# Physical constants
G = 6.674e-11
c = 2.998e8
H0 = 70 * 1e3 / 3.086e22
Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2
a0 = c * H0 * Omega_m**phi

def coherence(a):
    if a <= 0:
        return Omega_m
    x = (a / a0) ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def G_eff_over_G(a):
    return 1.0 / coherence(a)

# Abell 520 parameters
M_total = 1e15 * 1.989e30  # kg
rho_crit = 3 * H0**2 / (8 * np.pi * G)
r_200 = (3 * M_total / (4 * np.pi * 200 * rho_crit))**(1/3)

print(f"Abell 520 estimated parameters:")
print(f"M_total ~ 10^15 M_sun")
print(f"r_200 ~ {r_200/3.086e22:.2f} Mpc")

# Acceleration at different radii
radii_kpc = [100, 300, 500, 1000, 2000]
print(f"\nAcceleration profile:")
print(f"{'r (kpc)':<12} {'a/a₀':<12} {'G_eff/G':<12}")
print("-"*40)

for r_kpc in radii_kpc:
    r = r_kpc * 3.086e19  # kpc to m
    # Simplified: assume NFW-ish enclosed mass
    r_s = r_200 / 4  # typical concentration
    x = r / r_s
    f_x = np.log(1 + x) - x / (1 + x)
    f_c = np.log(1 + 4) - 4 / 5
    M_enc = M_total * f_x / f_c

    a = G * M_enc / r**2
    print(f"{r_kpc:<12} {a/a0:<12.3f} {G_eff_over_G(a):<12.3f}")

print("""
OBSERVATIONS:
- Throughout Abell 520, a ~ a₀ (transition regime)
- G_eff/G ~ 1.4-1.7
- This is similar to other clusters

NO SPECIAL G_eff EFFECT in Abell 520 specifically.
The "dark core" must have other explanations.
""")

# =============================================================================
# COMPARISON: BULLET vs ABELL 520
# =============================================================================

print("\n" + "="*70)
print("COMPARISON: BULLET CLUSTER vs ABELL 520")
print("="*70)

print("""
                    BULLET CLUSTER          ABELL 520
                    ==============          =========
Merger geometry     Simple, binary          Complex, multiple subclusters
Collision stage     Post-collision          Mid-collision?
Viewing angle       Near plane-of-sky       Uncertain, possibly L.O.S.

Lensing peaks       2 main peaks            Multiple peaks, complex
Peak locations      On galaxies             Mixed (some on gas?)
Gas location        Between peaks           Between some peaks

ΛCDM interpretation:
Bullet: Clear success - DM collisionless, follows galaxies
A520:   Initially puzzling, later "explained" by complex geometry

SYNCHRONISM interpretation:
Bullet: Success - Indifferent patterns collisionless (Session #197)
A520:   Consistent - Complex geometry, projection effects

KEY INSIGHT:
Both frameworks invoke similar explanations for Abell 520:
- Complex merger geometry
- Projection effects
- Multiple overlapping structures

Neither framework requires modification to explain A520.
""")

# =============================================================================
# PREDICTIONS FOR ABELL 520
# =============================================================================

print("\n" + "="*70)
print("SYNCHRONISM PREDICTIONS FOR ABELL 520")
print("="*70)

print("""
SPECIFIC PREDICTIONS:

1. MASS RATIOS
   M_dyn/M_lens should follow same G_eff pattern as other clusters
   - No special physics in A520
   - Same ~1.1-1.3 ratio expected (with anisotropy)

2. RADIAL TREND
   If radial M_dyn/M_lens data available:
   - Should show increasing trend with radius
   - Same as other clusters (Session #199 prediction)

3. "DARK CORE" INTERPRETATION
   If real (not projection artifact):
   - Some indifferent mass CAN be stripped
   - But this requires self-interaction
   - Constrains interaction cross-section

4. SUBSTRUCTURE
   Each subcluster should show:
   - Indifferent mass following galaxies
   - G_eff enhancement in low-a regions
   - Mass ratio consistent with f_indifferent ~ 4

FALSIFICATION SCENARIO:
If Abell 520 shows QUALITATIVELY different M_dyn/M_lens
than other clusters, this would require explanation.

Current data: No such anomaly reported.
""")

# =============================================================================
# CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("CONCLUSIONS: ABELL 520")
print("="*70)

print("""
SUMMARY:

1. INITIAL PUZZLE LARGELY RESOLVED
   - "Dark core" was overinterpreted
   - Complex geometry + projection effects
   - Multiple mass peaks found with deeper data

2. NO CHALLENGE TO SYNCHRONISM
   - Same framework as Bullet Cluster applies
   - Indifferent patterns are collisionless
   - G_eff enhancement is present but not special

3. USEFUL AS CONSISTENCY CHECK
   - If A520 showed anomalous M_dyn/M_lens, would be concern
   - No such anomaly reported
   - Consistent with framework

4. REMAINING QUESTION
   - Is there ANY real "dark core" without galaxies?
   - If so, implies some indifferent mass can separate
   - Would constrain self-interaction cross-section
   - Current data: Unclear if truly anomalous

NEXT:
- Look for M_dyn/M_lens radial data in A520 specifically
- Compare to other merger systems
- Consider self-interaction constraints
""")
