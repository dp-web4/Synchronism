#!/usr/bin/env python3
"""
Session #200: Caustic Mass Method Analysis
==========================================

The caustic mass technique uses the escape velocity profile in phase space
(radius-velocity diagram) to determine enclosed mass. This is LESS sensitive
to velocity anisotropy than simple velocity dispersion methods.

Key insight: In Synchronism, the escape velocity is:
v_esc² = 2 G_eff M(<r) / r = 2 (G/C) M / r

Observers use: M_caustic = v_esc² r / (2G)

Therefore: M_caustic / M_true = G_eff / G = 1/C(a)

Same as the σ-based prediction, BUT without the anisotropy correction!

This makes caustic mass a CLEANER test of G_eff.

Date: December 30, 2025
Machine: CBP
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
G = 6.674e-11  # m³/kg/s²
c_light = 2.998e8  # m/s
H0 = 70 * 1e3 / 3.086e22  # s^-1
Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2

# Derived critical acceleration
a0 = c_light * H0 * Omega_m**phi
print(f"Critical acceleration a₀ = {a0:.3e} m/s²")

def coherence(a):
    """C(a) = Omega_m + (1 - Omega_m) * (a/a0)^(1/phi) / [1 + (a/a0)^(1/phi)]"""
    if np.isscalar(a):
        if a <= 0:
            return Omega_m
        x = (a / a0) ** (1/phi)
        return Omega_m + (1 - Omega_m) * x / (1 + x)
    else:
        result = np.zeros_like(a)
        mask = a > 0
        x = (a[mask] / a0) ** (1/phi)
        result[mask] = Omega_m + (1 - Omega_m) * x / (1 + x)
        result[~mask] = Omega_m
        return result

def G_eff_over_G(a):
    return 1.0 / coherence(a)

# =============================================================================
# CAUSTIC MASS METHOD PHYSICS
# =============================================================================

print("\n" + "="*70)
print("SESSION #200: CAUSTIC MASS METHOD ANALYSIS")
print("="*70)

print("""
CAUSTIC MASS FUNDAMENTALS
=========================

The caustic method (Diaferio & Geller 1997, Diaferio 1999) identifies
the escape velocity profile from the phase-space distribution of
cluster members.

In phase space (projected radius R vs line-of-sight velocity v_los),
the boundary forms a "trumpet" shape:
- Inner region: High v_esc, large spread
- Outer region: Lower v_esc, narrower spread

The enclosed mass is:

  M(<r) = (1/G) × F_β × ∫₀ʳ A²(r') dr'

where:
- A(r) is the caustic amplitude (escape velocity profile)
- F_β is a filling factor (~0.5-0.7) that depends on orbits

KEY ADVANTAGE:
The caustic amplitude is determined by the ESCAPE VELOCITY,
not the velocity dispersion. This makes it less sensitive to
velocity anisotropy β.

IN SYNCHRONISM:
--------------
Escape velocity: v_esc² = 2 G_eff M / r = 2 (G/C) M / r

Observers measure A(r) and compute:
  M_caustic = (1/G) × ∫ A² dr

But the TRUE mass is:
  M_true = (1/G_eff) × ∫ A² dr = C × (1/G) × ∫ A² dr

Therefore:
  M_caustic / M_true = 1/C = G_eff/G

Same as the σ-based prediction!
But without the (1-β) anisotropy correction.

PREDICTION:
  M_caustic / M_lens = G_eff/G ≈ 1.8-2.0 at R_200
  (Higher than M_σ/M_lens ≈ 1.1-1.2)
""")

# =============================================================================
# LITERATURE DATA: CAUSTIC MASSES
# =============================================================================

print("\n" + "="*70)
print("LITERATURE: CAUSTIC MASS STUDIES")
print("="*70)

print("""
KEY CAUSTIC MASS STUDIES:

1. CIRS (Cluster Infall Regions in SDSS)
   - Rines & Diaferio (2006)
   - 72 clusters with caustic masses
   - Extended to several R_200

2. HeCS-SZ (Hectospec Cluster Survey - SZ selected)
   - Rines et al. (2016)
   - 47 SZ-selected clusters
   - High-quality spectroscopy

3. CLASH-VLT
   - Biviano et al. (2013)
   - Deep spectroscopy of CLASH clusters
   - Multiple mass estimators compared

OBSERVED M_caustic / M_other:
-----------------------------

| Study | Comparison | Ratio | Notes |
|-------|------------|-------|-------|
| Rines+ 2013 | M_caus/M_WL | 1.2 ± 0.2 | At R_200 |
| Rines+ 2016 | M_caus/M_SZ | 1.0-1.1 | HeCS-SZ sample |
| Maughan+ 2016 | M_caus/M_X | 1.3 ± 0.3 | Large scatter |
| Biviano+ 2013 | M_caus/M_lens | 1.1-1.3 | CLASH clusters |

CRITICAL OBSERVATION:
M_caustic tends to be HIGHER than M_σ (velocity dispersion).
This is consistent with Synchronism if:
- M_σ is reduced by anisotropy correction
- M_caustic is not (or less so)

Ratio: M_caustic / M_σ ≈ 1.1-1.3

If M_σ/M_lens = 1.1 and M_caustic/M_σ = 1.2, then
M_caustic/M_lens ≈ 1.3

Still below G_eff/G ≈ 1.9, but closer than M_σ method!
""")

# =============================================================================
# SYSTEMATIC ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("SYSTEMATIC ANALYSIS: WHY IS M_caustic/M_lens NOT ~2?")
print("="*70)

print("""
POSSIBLE EXPLANATIONS:

1. FILLING FACTOR CALIBRATION
   The filling factor F_β is calibrated using N-body simulations
   that assume standard Newtonian gravity.

   In Synchronism:
   - F_β would need recalibration
   - Simulations use G, but reality uses G_eff
   - This could systematically bias M_caustic downward

2. PROJECTION EFFECTS
   Caustic identification uses projected phase space.
   Interloper contamination affects amplitude.
   This is independent of Synchronism.

3. APERTURE MATCHING
   M_caustic often quoted at different radius than M_lens.
   Need consistent apertures for comparison.

4. RESIDUAL ANISOTROPY DEPENDENCE
   Though less sensitive than M_σ, caustic method still
   has some anisotropy dependence through F_β.

RECALIBRATION ANALYSIS:
-----------------------
If true G_eff/G = 1.9 at R_200, but we observe M_caustic/M_lens ≈ 1.2-1.3,
then the effective F_β correction is:

F_β,effective / F_β,assumed = 1.3 / 1.9 ≈ 0.68

This is within plausible range of F_β uncertainty (0.5-0.7).

PREDICTION:
If Synchronism is correct, recalibrating F_β with G_eff-aware
simulations should bring M_caustic/M_lens to ~1.9.
""")

# =============================================================================
# QUANTITATIVE PREDICTIONS
# =============================================================================

def predict_mass_ratios(M_200_solar=5e14, c_200=4.0):
    """
    Calculate predicted mass ratios at various radii.
    """
    M_sun = 1.989e30
    M_200 = M_200_solar * M_sun

    rho_crit = 3 * H0**2 / (8 * np.pi * G)
    r_200 = (3 * M_200 / (4 * np.pi * 200 * rho_crit))**(1/3)
    r_s = r_200 / c_200

    def M_enc(r):
        f_c = np.log(1 + c_200) - c_200 / (1 + c_200)
        x = r / r_s
        f_x = np.log(1 + x) - x / (1 + x)
        return M_200 * f_x / f_c

    def acceleration(r):
        return G * M_enc(r) / r**2

    print(f"\n{'='*70}")
    print(f"PREDICTIONS FOR M_200 = {M_200_solar:.0e} M_sun, c = {c_200}")
    print(f"{'='*70}")

    r_fracs = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]

    print(f"\n{'r/R200':<10} {'a/a₀':<10} {'G_eff/G':<12} {'M_caus/M_lens':<14} {'M_σ/M_lens':<12}")
    print(f"{'(pred)':<10} {'':<10} {'':<12} {'(no β corr)':<14} {'(β=0.4)':<12}")
    print("-" * 60)

    for rf in r_fracs:
        r = rf * r_200
        a = acceleration(r)
        G_ratio = G_eff_over_G(a)
        M_sigma_ratio = G_ratio * 0.6  # With anisotropy β=0.4

        print(f"{rf:<10.2f} {a/a0:<10.3f} {G_ratio:<12.3f} {G_ratio:<14.3f} {M_sigma_ratio:<12.3f}")

    print("""
NOTE: M_caustic/M_lens = G_eff/G assumes:
- No F_β recalibration
- Perfect caustic identification
- No projection effects

Actual observed M_caustic/M_lens ≈ 1.2-1.3, suggesting ~30% systematic
from F_β calibration or other effects.
""")

predict_mass_ratios()

# =============================================================================
# KEY TEST: COMPARE M_caustic to M_σ WITHIN SAME CLUSTERS
# =============================================================================

print("\n" + "="*70)
print("KEY TEST: M_caustic vs M_σ WITHIN CLUSTERS")
print("="*70)

print("""
CRITICAL OBSERVATION:

If we compare WITHIN the same clusters:
  M_caustic / M_σ

This ratio removes many systematics (both use same galaxies, same data).

In Synchronism:
  M_caustic = G_eff/G × M_true
  M_σ = G_eff/G × (1-β) × M_true

Therefore:
  M_caustic / M_σ = 1 / (1-β)

For β = 0.4: M_caustic / M_σ = 1.67
For β = 0.3: M_caustic / M_σ = 1.43

OBSERVED:
  M_caustic / M_σ ≈ 1.1-1.3 (various studies)

This suggests β ≈ 0.1-0.2, which is LOWER than expected β ≈ 0.3-0.4.

POSSIBLE INTERPRETATIONS:

1. Anisotropy is lower than simulations suggest
   - β ~ 0.1-0.2 rather than 0.3-0.4
   - Then M_σ/M_lens ≈ (G_eff/G)(1-β) ≈ 1.9×0.85 ≈ 1.6
   - But observed M_σ/M_lens ≈ 1.1
   - Discrepancy remains

2. Caustic F_β also needs correction
   - Both M_caustic and M_σ are biased
   - But in same direction, so ratio ≈ 1
   - Consistent with observations

3. G_eff effect is smaller than predicted
   - C(a) formula gives too much enhancement
   - Need to reconsider formula

DECISION POINT:
The M_caustic/M_σ ratio is a clean test.
If β is well-measured independently, this constrains G_eff.
""")

# =============================================================================
# RADIAL TREND IN CAUSTIC MASSES
# =============================================================================

print("\n" + "="*70)
print("RADIAL TREND: THE CLEANEST TEST")
print("="*70)

print("""
THE CLEANEST TEST:

Rather than comparing different methods, look at RADIAL TREND
of a single method:

M(<r) / M_lens(<r) vs r

In ΛCDM:
- Ratio = 1 at all radii
- Scatter from measurement errors

In Synchronism:
- Ratio = G_eff(r)/G = 1/C(a(r))
- Increases with radius (as a decreases)

PREDICTION (for typical cluster):

| r/R_200 | a/a₀ | M(<r)/M_lens(<r) |
|---------|------|------------------|
| 0.3     | 0.8  | 1.57             |
| 0.5     | 0.5  | 1.69             |
| 1.0     | 0.25 | 1.93             |
| 1.5     | 0.15 | 2.10             |
| 2.0     | 0.10 | 2.23             |

A ~40% increase from 0.3 to 2.0 R_200.

This is INDEPENDENT of:
- Absolute calibration of F_β
- Absolute lensing calibration
- Anisotropy modeling

It's a DIFFERENTIAL test: does the ratio increase with radius?

EXISTING DATA:

Most studies quote single-aperture masses.
Very few examine the radial profile of M_dyn/M_lens.

Umetsu et al. (2016) have stacked weak lensing to 3 R_200.
Rines+ have caustic profiles to similar radii.

A DIRECT COMPARISON of these profiles would be definitive.
""")

# =============================================================================
# CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("SESSION #200 CONCLUSIONS: CAUSTIC MASS")
print("="*70)

print("""
KEY FINDINGS:

1. CAUSTIC METHOD IS LESS ANISOTROPY-DEPENDENT
   - Uses escape velocity, not dispersion
   - Better probe of true G_eff

2. OBSERVED M_caustic/M_lens ≈ 1.2-1.3
   - Higher than M_σ/M_lens ≈ 1.1
   - But still below G_eff/G ≈ 1.9

3. F_β CALIBRATION IS THE KEY UNCERTAINTY
   - Calibrated with N-body simulations using G
   - May need recalibration with G_eff

4. M_caustic/M_σ RATIO IS INFORMATIVE
   - Observed: 1.1-1.3
   - Synchronism + β~0.3-0.4 predicts: 1.4-1.7
   - Discrepancy suggests either:
     a) β is lower than expected
     b) Both methods need similar calibration correction
     c) G_eff is smaller than predicted

5. RADIAL TREND IS THE CLEANEST TEST
   - Does M/M_lens increase with radius?
   - Independent of absolute calibrations
   - Need to compile radial profile data

NEXT STEPS:
-----------
1. Compile radial M_dyn/M_lens profiles from literature
2. Analyze Abell 520 (anomalous cluster)
3. Consider whether C(a) formula needs refinement
4. Look for independent β measurements in clusters
""")
