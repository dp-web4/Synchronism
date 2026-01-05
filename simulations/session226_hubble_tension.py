#!/usr/bin/env python3
"""
Session #226: The Hubble Tension from Coherence Physics

Building on Session #224's void-dominated cosmology, this session derives
quantitative predictions for the Hubble tension.

Key Hypothesis:
- Local H₀ measurements probe higher-density environments (clusters, filaments)
- CMB H₀ measurements probe the void-dominated cosmic average
- Different environments experience different effective dark energy → different H₀

The Hubble Tension:
- Local (Cepheids): H₀ = 73.0 ± 1.0 km/s/Mpc (Riess et al. 2022)
- CMB (Planck): H₀ = 67.4 ± 0.5 km/s/Mpc (Planck 2018)
- Tension: ~9% difference, >5σ significance

Date: January 5, 2026
Machine: CBP
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq

# =============================================================================
# PART 1: CONSTANTS AND COHERENCE FUNCTION
# =============================================================================

print("=" * 70)
print("SESSION #226: THE HUBBLE TENSION FROM COHERENCE PHYSICS")
print("=" * 70)

# Physical constants
c = 2.998e8          # m/s
G = 6.674e-11        # m³/(kg·s²)

# Planck 2018 cosmological parameters (CMB-derived)
H_0_cmb = 67.4       # km/s/Mpc (CMB value)
Omega_m_cmb = 0.315  # Matter density (CMB value)
Omega_Lambda_cmb = 0.685  # Dark energy density (CMB value)

# Local measurements
H_0_local = 73.0     # km/s/Mpc (Cepheid value, Riess et al.)

# Golden ratio
phi = (1 + np.sqrt(5)) / 2

# Critical acceleration
a_0 = 1.2e-10  # m/s² (MOND scale)

def coherence_function(a, alpha=1/phi):
    """
    Coherence function C(a) from Sessions #217-224.

    C(a) = Ω_m + (1 - Ω_m) × (a/a₀)^α / [1 + (a/a₀)^α]
    """
    if a <= 0:
        return Omega_m_cmb

    x = (a / a_0) ** alpha
    return Omega_m_cmb + (1 - Omega_m_cmb) * x / (1 + x)

def dark_energy_ratio(a):
    """
    Dark energy to matter ratio from coherence.

    ρ_dark/ρ_m = 1/C(a) - 1
    """
    C = coherence_function(a)
    return 1/C - 1

def effective_omega_lambda(a):
    """
    Effective dark energy density parameter at acceleration a.

    Ω_Λ,eff = (ρ_dark/ρ_m) / (1 + ρ_dark/ρ_m)
    """
    ratio = dark_energy_ratio(a)
    return ratio / (1 + ratio)

def effective_omega_m(a):
    """
    Effective matter density parameter at acceleration a.

    Ω_m,eff = 1 - Ω_Λ,eff
    """
    return 1 - effective_omega_lambda(a)


# =============================================================================
# PART 2: ENVIRONMENT CHARACTERIZATION
# =============================================================================

print("\n" + "=" * 70)
print("PART 2: COSMIC ENVIRONMENT CHARACTERIZATION")
print("=" * 70)

# Environment properties
# Based on cosmic web classification (voids, filaments, walls, clusters)

# Typical accelerations in different environments
# Using gravitational acceleration at characteristic scale

def environment_acceleration(rho_over_rho_mean, scale_Mpc):
    """
    Characteristic gravitational acceleration in an environment.

    Parameters:
        rho_over_rho_mean: density relative to cosmic mean
        scale_Mpc: characteristic scale in Mpc

    Returns:
        acceleration in m/s²
    """
    # Mean cosmic density today
    H_0_SI = H_0_cmb * 1e3 / 3.086e22  # s⁻¹
    rho_mean = 3 * H_0_SI**2 / (8 * np.pi * G) * Omega_m_cmb

    # Environment density
    rho = rho_over_rho_mean * rho_mean

    # Scale in meters
    R = scale_Mpc * 3.086e22

    # Mass within scale
    M = rho * (4/3) * np.pi * R**3

    # Gravitational acceleration
    a = G * M / R**2

    return a

# Define environments
environments = {
    'Deep void': {
        'rho_over_mean': 0.1,    # 10% of mean density
        'scale_Mpc': 30,         # 30 Mpc typical void
        'volume_fraction': 0.50  # 50% of volume in deep voids
    },
    'Underdense': {
        'rho_over_mean': 0.5,
        'scale_Mpc': 20,
        'volume_fraction': 0.20
    },
    'Mean density': {
        'rho_over_mean': 1.0,
        'scale_Mpc': 10,
        'volume_fraction': 0.15
    },
    'Filament': {
        'rho_over_mean': 3.0,
        'scale_Mpc': 5,
        'volume_fraction': 0.10
    },
    'Cluster': {
        'rho_over_mean': 100.0,
        'scale_Mpc': 2,
        'volume_fraction': 0.05
    }
}

print("\nCosmic Environment Properties:")
print("-" * 85)
print(f"{'Environment':<15} {'ρ/⟨ρ⟩':>8} {'Scale':>8} {'a (m/s²)':>12} {'C(a)':>8} {'Ω_Λ,eff':>10}")
print("-" * 85)

env_data = {}
for name, props in environments.items():
    a = environment_acceleration(props['rho_over_mean'], props['scale_Mpc'])
    C = coherence_function(a)
    Omega_Lambda_eff = effective_omega_lambda(a)
    env_data[name] = {
        'a': a,
        'C': C,
        'Omega_Lambda_eff': Omega_Lambda_eff,
        'volume_fraction': props['volume_fraction']
    }
    print(f"{name:<15} {props['rho_over_mean']:>8.1f} {props['scale_Mpc']:>6} Mpc "
          f"{a:>12.3e} {C:>8.4f} {Omega_Lambda_eff:>10.4f}")


# =============================================================================
# PART 3: FRIEDMANN EQUATION WITH ENVIRONMENT-DEPENDENT Ω_Λ
# =============================================================================

print("\n" + "=" * 70)
print("PART 3: FRIEDMANN EQUATION AND H₀ DERIVATION")
print("=" * 70)

print("""
The Friedmann equation at z = 0:

    H² = H₀² × [Ω_m + Ω_Λ]

If Ω_Λ depends on environment, then H(environment) varies.

For a given H₀ (from CMB), the expansion rate in different environments:

    H_env = H₀_cmb × √[Ω_m + Ω_Λ,env]

But this isn't quite right - the CMB measures the GLOBAL H₀, which is a
volume-weighted average. Let's derive the relationship properly.
""")


# =============================================================================
# PART 4: VOLUME-WEIGHTED VS LOCAL H₀
# =============================================================================

print("\n" + "=" * 70)
print("PART 4: VOLUME-WEIGHTED VS LOCAL H₀")
print("=" * 70)

# Global (CMB) measurement: volume-weighted average
# H₀_cmb² = ⟨H²⟩ = Σ f_i × H_i² where f_i is volume fraction

# Local measurement: typically in higher-density regions (Cepheids in galaxies)
# H₀_local measures H in regions with ρ > ⟨ρ⟩

# The key insight: if Ω_Λ,eff varies with environment, then
# H² ∝ (Ω_m + Ω_Λ,eff) varies with environment

# Method 1: Direct calculation
# H_local / H_cmb = √[(Ω_m + Ω_Λ,local) / (Ω_m + Ω_Λ,global)]

# Volume-weighted average Ω_Λ
Omega_Lambda_global = sum(env_data[name]['Omega_Lambda_eff'] * env_data[name]['volume_fraction']
                          for name in environments)

print(f"Volume-weighted average Ω_Λ,eff: {Omega_Lambda_global:.4f}")
print(f"Standard ΛCDM Ω_Λ: {Omega_Lambda_cmb:.4f}")

# For local measurements, assume we're probing filament + mean density regions
# (Cepheids are in galaxies, which are in filaments and overdense regions)

local_weights = {
    'Deep void': 0.0,      # No Cepheids in voids
    'Underdense': 0.1,     # Few galaxies
    'Mean density': 0.3,   # Some galaxies
    'Filament': 0.5,       # Many Cepheid hosts
    'Cluster': 0.1         # Cluster ellipticals (fewer Cepheids)
}

# Normalize weights
total_weight = sum(local_weights.values())
local_weights = {k: v/total_weight for k, v in local_weights.items()}

Omega_Lambda_local = sum(env_data[name]['Omega_Lambda_eff'] * local_weights[name]
                         for name in environments)

print(f"\nLocal measurement weighted Ω_Λ,eff: {Omega_Lambda_local:.4f}")

# Calculate H₀ ratio
H_ratio = np.sqrt((Omega_m_cmb + Omega_Lambda_local) / (Omega_m_cmb + Omega_Lambda_global))

print(f"\nH₀_local / H₀_global = {H_ratio:.4f}")
print(f"Predicted H₀_local = {H_0_cmb * H_ratio:.2f} km/s/Mpc")
print(f"Observed H₀_local = {H_0_local:.2f} km/s/Mpc")


# =============================================================================
# PART 5: REFINED CALCULATION WITH LUMINOSITY DISTANCE
# =============================================================================

print("\n" + "=" * 70)
print("PART 5: LUMINOSITY DISTANCE AND H₀ MEASUREMENT")
print("=" * 70)

print("""
H₀ is measured from luminosity distance-redshift relation:

    d_L = c × (1+z) × ∫₀ᶻ dz' / H(z')

The key is that H(z) depends on the environment along the line of sight.

For local measurements (z < 0.1), the luminosity distance is approximately:

    d_L ≈ c × z / H₀ × (1 + z/2 × (1 - q₀))

where q₀ = Ω_m/2 - Ω_Λ is the deceleration parameter.

If the line of sight passes through different environments, the effective H₀
is environment-weighted.
""")

def hubble_parameter(z, Omega_m, Omega_Lambda):
    """Hubble parameter H(z) in km/s/Mpc."""
    return H_0_cmb * np.sqrt(Omega_m * (1+z)**3 + Omega_Lambda)

def luminosity_distance(z, Omega_m, Omega_Lambda):
    """
    Luminosity distance in Mpc.

    d_L = c × (1+z) × ∫₀ᶻ dz' / H(z')
    """
    c_kms = 299792.458  # km/s

    def integrand(zp):
        return 1 / hubble_parameter(zp, Omega_m, Omega_Lambda)

    integral, _ = quad(integrand, 0, z)
    return c_kms * (1 + z) * integral

# Compare luminosity distances in different environments at z = 0.05
z_test = 0.05

print(f"\nLuminosity distance at z = {z_test}:")
print("-" * 60)

for name, data in env_data.items():
    Omega_m = 1 - data['Omega_Lambda_eff']
    Omega_Lambda = data['Omega_Lambda_eff']
    d_L = luminosity_distance(z_test, Omega_m, Omega_Lambda)
    print(f"{name:<15}: Ω_m={Omega_m:.3f}, Ω_Λ={Omega_Lambda:.3f}, d_L = {d_L:.2f} Mpc")

# Global average
d_L_global = luminosity_distance(z_test, Omega_m_cmb, Omega_Lambda_cmb)
print(f"{'Global (CMB)':<15}: Ω_m={Omega_m_cmb:.3f}, Ω_Λ={Omega_Lambda_cmb:.3f}, d_L = {d_L_global:.2f} Mpc")


# =============================================================================
# PART 6: DERIVING H₀ FROM DISTANCE MEASUREMENTS
# =============================================================================

print("\n" + "=" * 70)
print("PART 6: INVERTING FOR H₀ FROM DISTANCE MEASUREMENTS")
print("=" * 70)

print("""
When we measure d_L and z, we infer H₀ assuming ΛCDM.

If the TRUE cosmology has environment-dependent Ω_Λ, the inferred H₀ differs.

Method:
1. Calculate TRUE d_L in local environment
2. Infer H₀ assuming standard ΛCDM
3. Compare with CMB H₀
""")

def infer_H0_from_distance(d_L_measured, z, Omega_m_assumed, Omega_Lambda_assumed):
    """
    Infer H₀ from measured d_L assuming ΛCDM cosmology.

    Given a measured d_L and z, find H₀ that reproduces d_L under ΛCDM.
    """
    c_kms = 299792.458

    def integrand(zp, H0):
        return 1 / (H0 * np.sqrt(Omega_m_assumed * (1+zp)**3 + Omega_Lambda_assumed))

    def d_L_model(H0):
        integral, _ = quad(lambda zp: integrand(zp, H0), 0, z)
        return c_kms * (1 + z) * integral

    # Find H₀ that gives the measured d_L
    # Use simple approximation for small z
    return c_kms * z * (1 + z/2) / d_L_measured


# Calculate for local environment
# Assume local measurements probe filament/mean density environments

# Weighted average Omega values for local measurements
Omega_Lambda_local_eff = sum(env_data[name]['Omega_Lambda_eff'] * local_weights[name]
                              for name in environments)
Omega_m_local_eff = 1 - Omega_Lambda_local_eff

print(f"\nLocal environment effective parameters:")
print(f"  Ω_m,local = {Omega_m_local_eff:.4f}")
print(f"  Ω_Λ,local = {Omega_Lambda_local_eff:.4f}")

# True luminosity distance in local environment
z_cepheid = 0.03  # Typical Cepheid host redshift
d_L_local_true = luminosity_distance(z_cepheid, Omega_m_local_eff, Omega_Lambda_local_eff)

print(f"\nAt z = {z_cepheid}:")
print(f"  True d_L (local env) = {d_L_local_true:.2f} Mpc")

# If we interpret this with global ΛCDM, what H₀ do we infer?
H0_inferred = infer_H0_from_distance(d_L_local_true, z_cepheid, Omega_m_cmb, Omega_Lambda_cmb)

print(f"  Inferred H₀ (assuming ΛCDM) = {H0_inferred:.2f} km/s/Mpc")
print(f"  True CMB H₀ = {H_0_cmb:.2f} km/s/Mpc")
print(f"  Observed local H₀ = {H_0_local:.2f} km/s/Mpc")


# =============================================================================
# PART 7: SYSTEMATIC ANALYSIS
# =============================================================================

print("\n" + "=" * 70)
print("PART 7: SYSTEMATIC ANALYSIS OF H₀ TENSION")
print("=" * 70)

# More careful analysis: the issue is that H₀ affects the distance scale

# The Hubble tension is:
# H₀_local / H₀_cmb = 73.0 / 67.4 = 1.083 (8.3% higher)

# Our coherence model predicts different expansion rates in different environments
# Let's calculate what fraction of the tension is explained

# Method: Use the fact that d_L ∝ 1/H₀ at low z
# If local d_L is smaller due to less dark energy, the inferred H₀ is larger

# The ratio of effective Ω_Λ
ratio_Omega_Lambda = Omega_Lambda_local_eff / Omega_Lambda_global

print(f"\nΩ_Λ ratio (local/global): {ratio_Omega_Lambda:.4f}")

# At z << 1, d_L ≈ cz/H₀, so H₀ ∝ 1/d_L
# But Ω_Λ affects the expansion history

# More careful: H(z)² = H₀² × [Ω_m(1+z)³ + Ω_Λ]
# At z = 0: H₀² = H₀_cmb² × [Ω_m + Ω_Λ]

# If local Ω_Λ < global Ω_Λ, then for the same physical H at z=0,
# the inferred H₀ would be LOWER, not higher.

# Wait - let me reconsider...

print("""
RECONSIDERATION:

The Hubble tension is H₀_local > H₀_cmb.

If local environments have LESS dark energy (lower Ω_Λ), they expand SLOWER.
This would mean H₀_local < H₀_cmb, opposite to observations!

Let me check the sign...
""")

# In voids (low ρ, low a), C(a) is SMALL → more dark energy (Ω_Λ,eff is HIGH)
# In clusters (high ρ, high a), C(a) is LARGE → less dark energy (Ω_Λ,eff is LOW)

# So LOCAL measurements (in overdense regions) see LESS dark energy → SLOWER expansion
# This is the WRONG DIRECTION for the Hubble tension!

print("\nCONTRADICTION IDENTIFIED:")
print("-" * 60)
print("Coherence physics predicts:")
print("  - Voids: high Ω_Λ,eff → fast expansion")
print("  - Clusters: low Ω_Λ,eff → slow expansion")
print("  - Local (overdense): H₀_local < H₀_cmb")
print("")
print("But observations show:")
print("  - H₀_local > H₀_cmb")
print("")
print("This is the OPPOSITE of what we need!")


# =============================================================================
# PART 8: RESOLVING THE CONTRADICTION
# =============================================================================

print("\n" + "=" * 70)
print("PART 8: RESOLVING THE CONTRADICTION")
print("=" * 70)

print("""
Wait - let me reconsider the physics more carefully.

The issue is HOW H₀ is measured locally:

1. CEPHEID METHOD:
   - Measure d_L to Cepheid host galaxies
   - Compare with redshift z
   - Infer H₀ = cz/d_L (approximately)

2. If local Ω_Λ is LOWER:
   - The universe expanded SLOWER locally
   - Objects are CLOSER than ΛCDM predicts
   - Measured d_L is SMALLER
   - Inferred H₀ = cz/d_L is LARGER!

So the sign might actually work out...

Let me redo the calculation more carefully.
""")

# The key is: d_L depends on the integral of 1/H(z)
# If Ω_Λ is lower, H(z) at z > 0 is lower (slower expansion in past)
# This means ∫dz/H(z) is larger → d_L is LARGER
# This gives SMALLER inferred H₀

# But wait - the expansion rate today H₀ is what we're measuring
# If Ω_Λ is lower and Ω_m is higher (to keep Ω_total = 1),
# then H₀² = H₀_cmb² × (Ω_m + Ω_Λ) is the SAME!

# The issue is the HISTORY, not the present rate.

# Let me think about this differently...

print("""
ALTERNATIVE INTERPRETATION:

The Hubble tension might NOT be explained by environment-dependent Ω_Λ
in the way I first thought.

Instead, consider:

1. The CMB measures H₀ × r_s (sound horizon × Hubble parameter)
2. If r_s is affected by coherence at recombination (Session #225!),
   the inferred H₀ could be systematically different

3. Local measurements are distance-ladder based:
   - Parallax → Cepheids → SNe Ia
   - Each rung could be affected by environment-dependent physics

Let me explore the SOUND HORIZON interpretation...
""")


# =============================================================================
# PART 9: SOUND HORIZON AND CMB H₀
# =============================================================================

print("\n" + "=" * 70)
print("PART 9: SOUND HORIZON AND CMB H₀ INFERENCE")
print("=" * 70)

print("""
The CMB doesn't measure H₀ directly. It measures:

    θ_* = r_s / D_A(z_*)

where:
- θ_* = acoustic scale angle
- r_s = sound horizon at recombination
- D_A(z_*) = angular diameter distance to recombination

Then H₀ is inferred from the full ΛCDM model.

If the sound horizon r_s is DIFFERENT due to coherence physics,
the inferred H₀ would be systematically shifted.
""")

# From Session #225: perturbation accelerations at recombination are in MOND regime
# This modifies the acoustic oscillation physics

# Sound horizon at recombination:
# r_s = ∫₀^(t_*) c_s dt / a(t)
# where c_s = c / √[3(1 + R)] is the sound speed
# R = 3ρ_b / 4ρ_γ

# If gravitational enhancement G_eff > G at the acoustic scale,
# the matter overdensities grow faster
# This could modify the relationship between θ_* and H₀

# From Session #225:
# Peak 1 (ℓ=220, scale~150 Mpc): G_eff/G = 1.017
# Peak 2 (ℓ=530, scale~60 Mpc): G_eff/G = 1.030

# The sound horizon is set by the balance of gravity and pressure
# If G_eff > G, the acoustic oscillation frequency changes

# ω_acoustic ∝ √(G_eff × ρ) → changes by √(G_eff/G)
# This affects the sound horizon

G_eff_acoustic = 1.02  # Typical enhancement at acoustic scale

print(f"Gravitational enhancement at acoustic scale: G_eff/G = {G_eff_acoustic:.3f}")

# The sound horizon scales as r_s ∝ c_s × t_rec
# where t_rec ∝ 1/√(G_eff × ρ)

# So r_s ∝ 1/√G_eff

r_s_modification = 1 / np.sqrt(G_eff_acoustic)
print(f"Sound horizon modification: r_s(coherence) / r_s(ΛCDM) = {r_s_modification:.4f}")

# If r_s is smaller, but θ_* is fixed by observation,
# then D_A must be smaller → H₀ must be larger!

# θ_* = r_s / D_A
# If r_s decreases by factor f, D_A must decrease by factor f
# But D_A ∝ 1/H₀ at z >> 1, so H₀ must increase by factor 1/f

H0_correction_factor = 1 / r_s_modification
print(f"\nH₀ correction factor: {H0_correction_factor:.4f}")

H0_cmb_corrected = H_0_cmb * H0_correction_factor
print(f"Corrected CMB H₀: {H0_cmb_corrected:.2f} km/s/Mpc")
print(f"This is {'higher' if H0_cmb_corrected > H_0_cmb else 'lower'} than raw CMB value of {H_0_cmb:.2f}")

# But wait - this goes the WRONG way!
# If coherence makes r_s smaller, the CMB H₀ should be inferred HIGHER
# This would INCREASE the tension with local measurements


# =============================================================================
# PART 10: THE REAL RESOLUTION - LATE TIME EFFECTS
# =============================================================================

print("\n" + "=" * 70)
print("PART 10: LATE-TIME DARK ENERGY VARIATION")
print("=" * 70)

print("""
Perhaps the resolution lies in LATE-TIME effects, not early universe:

1. Local H₀ measurements are at z ~ 0-0.1
2. At these redshifts, dark energy dominates
3. If Ω_Λ varies spatially, expansion rates vary spatially

KEY INSIGHT:
The distance ladder probes lines of sight that pass through specific environments.
If those environments systematically differ from cosmic average, H₀ will differ.
""")

# Consider the Cepheid host galaxies:
# - They're in spiral galaxies
# - Spirals are typically in filaments or sheets
# - NOT in deep voids

# The line of sight to z ~ 0.03 passes through:
# - 50% void (on average)
# - 30% filament/sheet
# - 15% cluster outskirts
# - 5% cluster cores

# But if Cepheid hosts are specifically in overdense regions,
# the sampling is biased

print("Line of sight analysis:")
print("-" * 60)

# For a random line of sight, the distance integral samples:
# d_L ∝ ∫ dz / H(z, environment)

# The EFFECTIVE H₀ along a line of sight depends on environment weighting

def effective_H0_line_of_sight(void_fraction, filament_fraction, cluster_fraction):
    """
    Calculate effective H₀ for a line of sight passing through
    different environment fractions.
    """
    # Environment-specific H₀
    # H² ∝ (Ω_m + Ω_Λ,eff)
    # At z ≈ 0, H = H₀

    H0_void = H_0_cmb * np.sqrt((Omega_m_cmb + env_data['Deep void']['Omega_Lambda_eff']) /
                                 (Omega_m_cmb + Omega_Lambda_cmb))
    H0_filament = H_0_cmb * np.sqrt((Omega_m_cmb + env_data['Filament']['Omega_Lambda_eff']) /
                                     (Omega_m_cmb + Omega_Lambda_cmb))
    H0_cluster = H_0_cmb * np.sqrt((Omega_m_cmb + env_data['Cluster']['Omega_Lambda_eff']) /
                                    (Omega_m_cmb + Omega_Lambda_cmb))

    # Effective H₀ is distance-weighted average of 1/H
    # d = ∫ dz/H(z) ≈ z × (1/H₀_eff)

    # 1/H₀_eff = Σ f_i / H₀_i
    H0_eff_inv = (void_fraction / H0_void +
                  filament_fraction / H0_filament +
                  cluster_fraction / H0_cluster)

    return 1 / H0_eff_inv

# Random line of sight
H0_random = effective_H0_line_of_sight(0.50, 0.30, 0.20)
print(f"Random line of sight H₀: {H0_random:.2f} km/s/Mpc")

# Cepheid-selected line of sight (biased toward overdense)
H0_cepheid = effective_H0_line_of_sight(0.20, 0.50, 0.30)
print(f"Cepheid-selected line of sight H₀: {H0_cepheid:.2f} km/s/Mpc")

# The difference
print(f"\nDifference: {(H0_cepheid/H0_random - 1)*100:.2f}%")


# =============================================================================
# PART 11: QUANTITATIVE PREDICTION
# =============================================================================

print("\n" + "=" * 70)
print("PART 11: QUANTITATIVE PREDICTION FOR HUBBLE TENSION")
print("=" * 70)

# Let me recalculate with proper weighting

# Environment properties with corrected H₀ values
for name, data in env_data.items():
    Omega_total = Omega_m_cmb + data['Omega_Lambda_eff']
    H0_env = H_0_cmb * np.sqrt(Omega_total / (Omega_m_cmb + Omega_Lambda_cmb))
    env_data[name]['H0'] = H0_env
    print(f"{name:<15}: Ω_Λ,eff = {data['Omega_Lambda_eff']:.4f}, H₀ = {H0_env:.2f} km/s/Mpc")

# Volume-weighted average H₀ (what CMB measures)
H0_global = sum(env_data[name]['H0'] * env_data[name]['volume_fraction']
                for name in environments)
print(f"\nVolume-weighted global H₀: {H0_global:.2f} km/s/Mpc")

# Local (Cepheid-weighted) H₀
H0_local_pred = sum(env_data[name]['H0'] * local_weights[name]
                     for name in environments)
print(f"Cepheid-weighted local H₀: {H0_local_pred:.2f} km/s/Mpc")

# Predicted tension
predicted_tension = (H0_local_pred - H0_global) / H0_global * 100
observed_tension = (H_0_local - H_0_cmb) / H_0_cmb * 100

print(f"\n" + "=" * 60)
print("HUBBLE TENSION COMPARISON")
print("=" * 60)
print(f"Predicted tension: {predicted_tension:.2f}%")
print(f"Observed tension:  {observed_tension:.2f}%")
print(f"Coherence explains: {abs(predicted_tension/observed_tension)*100:.1f}% of tension")


# =============================================================================
# PART 12: VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("PART 12: CREATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Environment-dependent Ω_Λ
ax1 = axes[0, 0]
env_names = list(environments.keys())
Omega_Lambda_values = [env_data[name]['Omega_Lambda_eff'] for name in env_names]
colors = ['blue', 'cyan', 'green', 'orange', 'red']

bars = ax1.bar(range(len(env_names)), Omega_Lambda_values, color=colors, edgecolor='black')
ax1.axhline(y=Omega_Lambda_cmb, color='black', linestyle='--', linewidth=2,
            label=f'ΛCDM: Ω_Λ = {Omega_Lambda_cmb}')
ax1.axhline(y=Omega_Lambda_global, color='purple', linestyle=':', linewidth=2,
            label=f'Global avg: {Omega_Lambda_global:.3f}')
ax1.set_xticks(range(len(env_names)))
ax1.set_xticklabels(env_names, rotation=45, ha='right')
ax1.set_ylabel('Effective Ω_Λ', fontsize=12)
ax1.set_title('Environment-Dependent Dark Energy', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel 2: Environment-dependent H₀
ax2 = axes[0, 1]
H0_values = [env_data[name]['H0'] for name in env_names]

bars = ax2.bar(range(len(env_names)), H0_values, color=colors, edgecolor='black')
ax2.axhline(y=H_0_cmb, color='black', linestyle='--', linewidth=2,
            label=f'CMB: H₀ = {H_0_cmb}')
ax2.axhline(y=H_0_local, color='red', linestyle='--', linewidth=2,
            label=f'Local: H₀ = {H_0_local}')
ax2.axhline(y=H0_global, color='purple', linestyle=':', linewidth=2,
            label=f'Global avg: {H0_global:.1f}')
ax2.set_xticks(range(len(env_names)))
ax2.set_xticklabels(env_names, rotation=45, ha='right')
ax2.set_ylabel('H₀ (km/s/Mpc)', fontsize=12)
ax2.set_title('Environment-Dependent Hubble Parameter', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Panel 3: Coherence function vs acceleration
ax3 = axes[1, 0]
a_range = np.logspace(-14, -8, 100)
C_values = [coherence_function(a) for a in a_range]

ax3.semilogx(a_range, C_values, 'b-', linewidth=2)
ax3.axhline(y=Omega_m_cmb, color='r', linestyle='--', label=f'C_min = Ω_m = {Omega_m_cmb}')
ax3.axhline(y=1, color='g', linestyle='--', label='C_max = 1')
ax3.axvline(x=a_0, color='purple', linestyle=':', label=f'a₀ = {a_0:.1e} m/s²')

# Mark environment accelerations
for name, data in env_data.items():
    ax3.scatter([data['a']], [data['C']], s=100, zorder=5, label=f'{name}')

ax3.set_xlabel('Acceleration (m/s²)', fontsize=12)
ax3.set_ylabel('Coherence C(a)', fontsize=12)
ax3.set_title('Coherence Function Across Environments', fontsize=14)
ax3.legend(loc='upper left', fontsize=8)
ax3.grid(True, alpha=0.3)

# Panel 4: Hubble tension summary
ax4 = axes[1, 1]
categories = ['CMB\n(Planck)', 'Predicted\nLocal', 'Observed\nLocal']
H0_comparison = [H_0_cmb, H0_local_pred, H_0_local]
errors = [0.5, 1.0, 1.0]  # Approximate errors
colors_bar = ['blue', 'green', 'red']

bars = ax4.bar(categories, H0_comparison, yerr=errors, color=colors_bar,
               edgecolor='black', capsize=5)

ax4.set_ylabel('H₀ (km/s/Mpc)', fontsize=12)
ax4.set_title('Hubble Tension: Prediction vs Observation', fontsize=14)
ax4.set_ylim(65, 76)
ax4.grid(True, alpha=0.3, axis='y')

# Add text annotations
for i, (h, err) in enumerate(zip(H0_comparison, errors)):
    ax4.text(i, h + err + 0.5, f'{h:.1f}', ha='center', fontsize=11, fontweight='bold')

# Add tension annotation
ax4.annotate('', xy=(2, H_0_local - 0.5), xytext=(0, H_0_cmb + 0.5),
             arrowprops=dict(arrowstyle='<->', color='black', lw=2))
ax4.text(1, (H_0_cmb + H_0_local)/2, f'Tension:\n{observed_tension:.1f}%',
         ha='center', va='center', fontsize=10, fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session226_hubble_tension.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session226_hubble_tension.png")


# =============================================================================
# PART 13: CONCLUSIONS
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #226: CONCLUSIONS")
print("=" * 70)

print(f"""
KEY FINDINGS:

1. ENVIRONMENT-DEPENDENT DARK ENERGY:
   - Deep voids: Ω_Λ,eff = {env_data['Deep void']['Omega_Lambda_eff']:.4f}
   - Mean density: Ω_Λ,eff = {env_data['Mean density']['Omega_Lambda_eff']:.4f}
   - Clusters: Ω_Λ,eff = {env_data['Cluster']['Omega_Lambda_eff']:.4f}

2. ENVIRONMENT-DEPENDENT H₀:
   - Deep voids: H₀ = {env_data['Deep void']['H0']:.2f} km/s/Mpc
   - Mean density: H₀ = {env_data['Mean density']['H0']:.2f} km/s/Mpc
   - Clusters: H₀ = {env_data['Cluster']['H0']:.2f} km/s/Mpc

3. HUBBLE TENSION PREDICTION:
   - CMB (volume-weighted): H₀ = {H0_global:.2f} km/s/Mpc
   - Local (Cepheid-weighted): H₀ = {H0_local_pred:.2f} km/s/Mpc
   - Predicted tension: {predicted_tension:.2f}%
   - Observed tension: {observed_tension:.2f}%

4. ASSESSMENT:
   - Coherence physics predicts OPPOSITE direction from observed tension!
   - Local (overdense) regions have LESS dark energy → SLOWER expansion
   - This would give H₀_local < H₀_cmb, not H₀_local > H₀_cmb

5. POSSIBLE RESOLUTIONS:
   a) The Hubble tension is NOT due to environment-dependent dark energy
   b) Local measurements are affected by other systematic effects
   c) The coherence model needs modification at low accelerations
   d) Early universe effects (sound horizon) are more important

6. NEXT STEPS:
   - Investigate sound horizon modification more carefully
   - Check if distance ladder calibration is affected
   - Consider time-dependent (not just environment-dependent) dark energy
   - Test against BAO measurements (which also show tension)

IMPORTANT NEGATIVE RESULT:
The naive coherence model does NOT explain the Hubble tension in terms of
environment-dependent expansion. This is valuable - it means the tension
is likely due to other physics (early universe, calibration, or new physics).
""")

print("\n" + "=" * 70)
print("SESSION #226 COMPLETE")
print("=" * 70)
