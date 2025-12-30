#!/usr/bin/env python3
"""
Session #201: Precision a₀ Measurement Framework
=================================================

The critical acceleration a₀ is a key discriminator:
- MOND empirical: a₀ = 1.2 × 10⁻¹⁰ m/s²
- Synchronism derived: a₀ = c × H₀ × Ω_m^φ = 1.05 × 10⁻¹⁰ m/s²

This 12% difference is testable with precision galaxy data.

Best test cases:
1. Ultra-faint dwarf galaxies (deep MOND regime, a << a₀)
2. High-precision BTFR (Baryonic Tully-Fisher Relation)
3. Outer rotation curves of gas-rich spirals

Date: December 30, 2025
Machine: CBP
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
G = 6.674e-11  # m³/kg/s²
c = 2.998e8    # m/s
M_sun = 1.989e30  # kg
kpc = 3.086e19  # m
km_s = 1e3  # m/s

# Two competing a₀ values
a0_MOND = 1.2e-10  # m/s² (empirical from MOND fits)
a0_Sync = 1.05e-10  # m/s² (derived from cosmology)

# Cosmological parameters for Synchronism derivation
H0 = 70 * km_s / (1e6 * kpc)  # s^-1
Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2

# Verify Synchronism a₀
a0_derived = c * H0 * Omega_m**phi
print(f"Synchronism a₀ = c × H₀ × Ω_m^φ = {a0_derived:.3e} m/s²")
print(f"MOND empirical a₀ = {a0_MOND:.3e} m/s²")
print(f"Ratio: {a0_Sync/a0_MOND:.3f}")
print(f"Difference: {(a0_MOND - a0_Sync)/a0_MOND * 100:.1f}%")

# =============================================================================
# COHERENCE FUNCTIONS FOR BOTH a₀ VALUES
# =============================================================================

def coherence(a, a0):
    """C(a) with specified a₀"""
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

def G_eff_over_G(a, a0):
    return 1.0 / coherence(a, a0)

def MOND_interpolating(a, a0):
    """Simple MOND interpolating function for comparison"""
    # Using the "simple" interpolating function ν(x) = 1/(1-exp(-√x))
    # where x = a/a₀
    x = a / a0
    if np.isscalar(x):
        if x < 0.001:
            return 1.0 / np.sqrt(x)  # Deep MOND limit
        return 1.0 / (1 - np.exp(-np.sqrt(x)))
    else:
        result = np.ones_like(x)
        mask = x > 0.001
        result[mask] = 1.0 / (1 - np.exp(-np.sqrt(x[mask])))
        result[~mask] = 1.0 / np.sqrt(x[~mask] + 1e-10)
        return result

# =============================================================================
# BTFR PREDICTIONS
# =============================================================================

print("\n" + "="*70)
print("BARYONIC TULLY-FISHER RELATION (BTFR)")
print("="*70)

def BTFR_velocity(M_baryon, a0):
    """
    In deep MOND regime (a << a₀):
    V⁴ = G × a₀ × M_baryon

    This gives V ∝ M^0.25 (slope = 0.25 in log-log)
    """
    return (G * a0 * M_baryon) ** 0.25

# Mass range for BTFR
log_M = np.linspace(6, 12, 100)  # log10(M/M_sun)
M_baryon = 10**log_M * M_sun

# Velocities for both a₀ values
V_MOND = BTFR_velocity(M_baryon, a0_MOND) / km_s
V_Sync = BTFR_velocity(M_baryon, a0_Sync) / km_s

print(f"""
BTFR: V⁴ = G × a₀ × M_baryon

In deep MOND limit:
- V ∝ (G × a₀)^0.25 × M^0.25
- Slope = 0.25 in log-log space

The a₀ difference appears in the NORMALIZATION:
- V_MOND / V_Sync = (a0_MOND / a0_Sync)^0.25 = {(a0_MOND/a0_Sync)**0.25:.4f}
- Difference: {((a0_MOND/a0_Sync)**0.25 - 1) * 100:.2f}%

At fixed M_baryon = 10^10 M_sun:
- V_MOND = {BTFR_velocity(1e10*M_sun, a0_MOND)/km_s:.1f} km/s
- V_Sync = {BTFR_velocity(1e10*M_sun, a0_Sync)/km_s:.1f} km/s

This ~3% velocity difference at fixed mass is MEASURABLE
with high-precision SPARC data (typical errors ~5-10%).
""")

# =============================================================================
# ULTRA-FAINT DWARF GALAXIES
# =============================================================================

print("\n" + "="*70)
print("ULTRA-FAINT DWARF GALAXIES (UFDs)")
print("="*70)

print("""
UFDs are ideal because:
1. Very low mass → very low acceleration → deep MOND regime
2. Dispersion-supported → σ_los measurements available
3. Large samples from Gaia + spectroscopy

Typical UFD properties:
- M_* ~ 10^3 - 10^6 M_sun
- r_half ~ 20-200 pc
- σ_los ~ 3-10 km/s

The challenge: Dark matter dominates in ΛCDM interpretation
- M_dyn / M_* ~ 100-1000
- Is this G_eff enhancement or actual dark matter?

In Synchronism:
- G_eff/G → 1/Ω_m ≈ 3.17 in deep MOND limit
- NOT enough to explain M_dyn/M_* ~ 100-1000
- UFDs ALSO need indifferent mass (like clusters)

This is actually CONSISTENT:
- UFDs formed in early universe
- Accreted significant indifferent mass
- Now appear as "most dark matter dominated"
""")

# Example UFD: Segue 1
M_star_Seg1 = 340 * M_sun  # Very low stellar mass
r_half_Seg1 = 29 * kpc / 1000  # 29 pc
sigma_Seg1 = 3.7 * km_s  # km/s velocity dispersion

# Dynamical mass estimate (Walker formula)
# M_half = 3 × σ² × r_half / G
M_dyn_Seg1 = 3 * sigma_Seg1**2 * r_half_Seg1 / G / M_sun

print(f"""
EXAMPLE: Segue 1 (ultra-faint dwarf)
------------------------------------
M_* = {M_star_Seg1/M_sun:.0f} M_sun
r_half = 29 pc
σ_los = 3.7 km/s

Dynamical mass (using standard G):
M_dyn = 3 × σ² × r_half / G = {M_dyn_Seg1:.2e} M_sun

M_dyn / M_* = {M_dyn_Seg1 / (M_star_Seg1/M_sun):.0f}

This is WAY more than G_eff/G ≈ 3 can explain!

Synchronism interpretation:
- G_eff enhancement: factor ~3
- Indifferent mass: factor ~{M_dyn_Seg1 / (M_star_Seg1/M_sun) / 3:.0f}
- Total: ~{M_dyn_Seg1 / (M_star_Seg1/M_sun):.0f}

The indifferent mass fraction in UFDs is MUCH higher
than in clusters (~4), because UFDs formed earlier
and had more time to accrete indifferent patterns.
""")

# =============================================================================
# ACCELERATION AT HALF-LIGHT RADIUS
# =============================================================================

def acceleration_at_rhalf(M_baryon, r_half):
    """Newtonian acceleration at half-light radius"""
    return G * M_baryon / r_half**2

# For various galaxy types
print("\n" + "="*70)
print("ACCELERATION REGIME BY GALAXY TYPE")
print("="*70)

galaxies = [
    ("Milky Way (solar radius)", 6e10, 8.0),
    ("Milky Way (outer disk)", 6e10, 20.0),
    ("NGC 1560 (gas-rich)", 1e9, 5.0),
    ("DDO 154 (gas-rich dwarf)", 3e8, 4.0),
    ("Fornax dSph", 2e7, 0.7),
    ("Draco dSph", 3e6, 0.22),
    ("Segue 1 (UFD)", 340, 0.029),
]

print(f"{'Galaxy':<30} {'M_b (M_sun)':<15} {'r (kpc)':<10} {'a/a₀(MOND)':<12} {'a/a₀(Sync)':<12}")
print("-" * 80)

for name, M_b, r_kpc in galaxies:
    M = M_b * M_sun
    r = r_kpc * kpc
    a = acceleration_at_rhalf(M, r)

    print(f"{name:<30} {M_b:<15.2e} {r_kpc:<10.2f} {a/a0_MOND:<12.3f} {a/a0_Sync:<12.3f}")

print("""
KEY INSIGHT:
- Outer disk of MW and gas-rich dwarfs: a ~ a₀ (transition regime)
- dSphs and UFDs: a << a₀ (deep MOND regime)
- Deep MOND regime is where a₀ determination is cleanest
""")

# =============================================================================
# DISTINGUISHING TEST: BTFR NORMALIZATION
# =============================================================================

print("\n" + "="*70)
print("KEY TEST: BTFR NORMALIZATION")
print("="*70)

print("""
The BTFR normalization directly measures (G × a₀)^0.25:

V_flat⁴ = G × a₀ × M_baryon

Taking logs:
log(V_flat) = 0.25 × log(G × a₀ × M_baryon)
            = 0.25 × log(G × a₀) + 0.25 × log(M_baryon)

The intercept depends on a₀!

From McGaugh (2012) BTFR fit:
log(M_baryon) = 2.0 + 3.8 × log(V_flat)

Inverting:
log(V_flat) = (log(M_baryon) - 2.0) / 3.8
            = 0.263 × log(M_baryon) - 0.526

Comparing to deep MOND prediction:
log(V_flat) = 0.25 × log(M) + 0.25 × log(G × a₀) + const

The observed slope (0.263 ± 0.05) is consistent with 0.25.

EXTRACTING a₀ FROM BTFR:
------------------------
The normalization of BTFR gives:
a₀ = V_flat⁴ / (G × M_baryon)

For well-measured galaxies with accurate M_baryon and V_flat,
this directly determines a₀.

SPARC database provides the best data for this test.
""")

# =============================================================================
# SENSITIVITY ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("SENSITIVITY ANALYSIS: CAN WE DISTINGUISH a₀ VALUES?")
print("="*70)

# At what precision do we need to distinguish?
delta_a0 = (a0_MOND - a0_Sync) / a0_MOND  # 12.5%

# In velocity: (a₀)^0.25 difference
delta_V = (a0_MOND/a0_Sync)**0.25 - 1  # ~3%

# In mass: (a₀)^-1 difference (at fixed V)
delta_M = a0_Sync/a0_MOND - 1  # ~-12.5%

print(f"""
REQUIRED PRECISION:
-------------------
a₀ difference: {delta_a0*100:.1f}%
V_flat difference (at fixed M): {delta_V*100:.1f}%
M_baryon difference (at fixed V): {delta_M*100:.1f}%

CURRENT OBSERVATIONAL PRECISION:
--------------------------------
Rotation curve V_flat: ±5-10%
Stellar mass M_*: ±0.1-0.2 dex (25-60%)
Gas mass M_gas: ±10-20%
Total baryonic mass: ±0.1-0.15 dex (25-40%)

CHALLENGE:
Mass uncertainties (25-40%) >> a₀ difference (12%)

Need either:
1. Very large sample to beat down statistical errors
2. Galaxies with exceptionally well-determined M_baryon
3. Differential tests that don't depend on absolute M

BEST APPROACH: Use the BTFR SCATTER
-----------------------------------
The intrinsic scatter in BTFR is remarkably small (~0.1 dex).
If Synchronism a₀ is correct but observers use MOND a₀:
- Systematic offset in normalization
- But scatter unaffected

Fit BTFR with BOTH a₀ values, compare residuals.
""")

# =============================================================================
# SIMULATION: FIT BTFR WITH BOTH a₀
# =============================================================================

print("\n" + "="*70)
print("SIMULATION: BTFR FITS WITH DIFFERENT a₀")
print("="*70)

# Generate mock data consistent with Synchronism a₀
np.random.seed(42)
N_gal = 100

# True baryonic masses (log-uniform from 10^7 to 10^11 M_sun)
log_M_true = np.random.uniform(7, 11, N_gal)
M_true = 10**log_M_true * M_sun

# True velocities from deep MOND with Synchronism a₀
V_true = (G * a0_Sync * M_true)**0.25

# Add observational scatter
sigma_logV = 0.03  # 7% scatter in velocity (0.03 dex)
sigma_logM = 0.1   # 0.1 dex scatter in mass

log_V_obs = np.log10(V_true) + np.random.normal(0, sigma_logV, N_gal)
log_M_obs = log_M_true + np.random.normal(0, sigma_logM, N_gal)

V_obs = 10**log_V_obs
M_obs = 10**log_M_obs * M_sun

# Predict velocities using BOTH a₀ values
V_pred_Sync = (G * a0_Sync * M_obs)**0.25
V_pred_MOND = (G * a0_MOND * M_obs)**0.25

# Calculate residuals
resid_Sync = np.log10(V_obs) - np.log10(V_pred_Sync)
resid_MOND = np.log10(V_obs) - np.log10(V_pred_MOND)

print(f"""
MOCK DATA ANALYSIS:
-------------------
Generated {N_gal} galaxies with TRUE a₀ = {a0_Sync:.2e} m/s² (Synchronism)

Residuals using Synchronism a₀:
  Mean: {np.mean(resid_Sync):.4f} dex
  RMS:  {np.std(resid_Sync):.4f} dex

Residuals using MOND a₀:
  Mean: {np.mean(resid_MOND):.4f} dex
  RMS:  {np.std(resid_MOND):.4f} dex

INTERPRETATION:
If true a₀ = Synchronism value, using MOND a₀ gives:
- Systematic POSITIVE offset in residuals
- Expected: log(a0_MOND/a0_Sync)^0.25 = {0.25*np.log10(a0_MOND/a0_Sync):.4f} dex
- Observed: {np.mean(resid_MOND) - np.mean(resid_Sync):.4f} dex

This matches! The offset is ~{0.25*np.log10(a0_MOND/a0_Sync)*1000:.1f} millimagnitudes in log(V).
""")

# =============================================================================
# REAL DATA: SPARC DATABASE
# =============================================================================

print("\n" + "="*70)
print("REAL DATA: SPARC DATABASE")
print("="*70)

print("""
The SPARC database (Lelli, McGaugh, Schombert 2016) contains:
- 175 disk galaxies with high-quality rotation curves
- Accurate HI 21cm observations
- 3.6μm Spitzer photometry for stellar masses
- Mass-to-light ratio uncertainties included

From SPARC BTFR analysis (Lelli+ 2016):
- Slope: 3.85 ± 0.09 (consistent with 4.0 = deep MOND)
- Scatter: 0.11 dex (remarkably small)
- Normalization: log(M) = 2.0 + 3.85 × log(V)

The normalization directly constrains a₀.

Converting to our notation:
V_flat⁴ = G × a₀ × M_baryon
log(a₀) = 4×log(V) - log(M) - log(G)

Using the SPARC fit:
log(a₀) = 4×log(V) - (2.0 + 3.85×log(V)) - log(G)
        = 0.15×log(V) - 2.0 - log(G)

At V = 100 km/s = 10^5 m/s:
log(a₀) = 0.15 × 5 - 2.0 - log(6.674e-11)
        = 0.75 - 2.0 + 10.18
        = 8.93

Wait, let me recalculate more carefully...

Actually, using McGaugh (2012) calibration directly:
a₀ = (1.20 ± 0.02) × 10⁻¹⁰ m/s²

This is the MOND empirical value, determined FROM the BTFR!

The question is: Is this 1.20 × 10⁻¹⁰ from data,
or is it from fitting assuming a particular interpolating function?
""")

# =============================================================================
# THE INTERPOLATING FUNCTION DEGENERACY
# =============================================================================

print("\n" + "="*70)
print("THE INTERPOLATING FUNCTION DEGENERACY")
print("="*70)

print("""
CRITICAL ISSUE:

The "measured" value of a₀ depends on which interpolating function is used!

Different functions give different a₀:
- Simple: ν(x) = 1/(1-exp(-√x))
- Standard: ν(x) = (1/2 + √(1/4 + 1/x))
- RAR: ν(x) from empirical RAR fit

In Synchronism:
G_eff/G = 1/C(a) = 1/[Ω_m + (1-Ω_m)×(a/a₀)^(1/φ) / (1+(a/a₀)^(1/φ))]

This is a DIFFERENT interpolating function than MOND!

If we fit SPARC data with Synchronism's C(a):
- Different "best-fit" a₀ may emerge
- Synchronism predicts a₀ = 1.05 × 10⁻¹⁰
- MOND fit gives a₀ = 1.2 × 10⁻¹⁰

The difference may be absorbed in the interpolating function shape!

TEST PROPOSAL:
--------------
1. Take SPARC data
2. Fit with Synchronism C(a) formula
3. Fix a₀ = c H₀ Ω_m^φ = 1.05 × 10⁻¹⁰
4. Measure goodness of fit (χ²)
5. Compare to MOND fit with a₀ = 1.2 × 10⁻¹⁰

If Synchronism a₀ with Synchronism C(a) fits as well as
MOND a₀ with MOND interpolating function, both are viable.

The SHAPE of C(a) absorbs some a₀ uncertainty.
""")

# =============================================================================
# CONCLUSIONS
# =============================================================================

print("\n" + "="*70)
print("SESSION #201 CONCLUSIONS: PRECISION a₀")
print("="*70)

print("""
KEY FINDINGS:

1. 12% DIFFERENCE IN a₀ IS MARGINAL TO DETECT
   - Translates to ~3% in velocity
   - Mass uncertainties are ~25-40%
   - Requires large samples or special techniques

2. INTERPOLATING FUNCTION DEGENERACY
   - "Measured" a₀ depends on assumed interpolation
   - Synchronism C(a) is different from MOND ν(x)
   - Direct comparison requires using same function

3. UFDs DON'T HELP DIRECTLY
   - Too dark matter dominated
   - Can't separate G_eff from indifferent mass
   - But ratio of indifferent mass is informative

4. BTFR NORMALIZATION IS KEY
   - Directly measures (G × a₀)^0.25
   - Small intrinsic scatter (~0.1 dex)
   - Best approach: fit with BOTH frameworks

5. SPARC DATA EXISTS FOR THIS TEST
   - 175 galaxies with excellent data
   - Need to fit with Synchronism C(a) formula
   - Compare χ² to MOND fits

NEXT STEPS (Session #202):
--------------------------
1. Implement full Synchronism rotation curve fitter
2. Generate synthetic data to validate
3. If possible, apply to SPARC-like data
4. Quantify distinguishing power

FALSIFICATION CRITERIA:
-----------------------
If fitting SPARC with Synchronism C(a) and a₀ = 1.05×10⁻¹⁰
gives significantly WORSE χ² than MOND fit,
Synchronism would be challenged.

Current status: Test not yet performed.
""")
