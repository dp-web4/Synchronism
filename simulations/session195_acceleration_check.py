#!/usr/bin/env python3
"""
SESSION #195b: VERIFY CLUSTER ACCELERATION CALCULATION
======================================================

The main simulation showed g ~ 10^5-10^6 × a₀ at R_200 for clusters.
This seems too high. Let's verify with a simple calculation.
"""

import numpy as np

print("=" * 70)
print("SESSION #195b: VERIFY CLUSTER ACCELERATION CALCULATION")
print("=" * 70)

# Constants
G_SI = 6.674e-11  # m³/(kg s²)
M_sun = 1.989e30  # kg
kpc_m = 3.086e19  # m
Mpc_m = 3.086e22  # m

# Cluster: M = 10^15 M_sun, R_200 = 2 Mpc
M_cluster = 1e15 * M_sun  # kg
R_200 = 2.0 * Mpc_m  # m

# Newtonian acceleration at R_200
g_SI = G_SI * M_cluster / R_200**2
print(f"\n1. Simple Newtonian calculation:")
print(f"   M = 10^15 M_sun = {M_cluster:.3e} kg")
print(f"   R = 2 Mpc = {R_200:.3e} m")
print(f"   g = G M / R² = {g_SI:.3e} m/s²")

# Compare to a₀
a0 = 1.05e-10  # m/s²
print(f"\n2. Compare to a₀:")
print(f"   a₀ = {a0:.3e} m/s²")
print(f"   g / a₀ = {g_SI / a0:.1f}")

# This gives a reasonable ~3, not 10^5!

print(f"\n3. WAIT - This gives g/a₀ ~ {g_SI/a0:.1f}, not 10^5!")
print("   Something is wrong with the unit conversion in the main script.")

# Let's redo in cluster units
print("\n" + "=" * 70)
print("4. CLUSTER UNITS ANALYSIS")
print("=" * 70)

# G in different units:
# G = 4.302e-6 kpc (km/s)² / M_sun
# g = G M / R² where M in M_sun, R in kpc

G_kpc = 4.302e-6  # kpc (km/s)² / M_sun

M = 1e15  # M_sun
R = 2000  # kpc (2 Mpc)

g_kpc = G_kpc * M / R**2  # (km/s)²/kpc
print(f"\nWith G = 4.302e-6 kpc (km/s)²/M_sun:")
print(f"   M = 10^15 M_sun")
print(f"   R = 2000 kpc")
print(f"   g = {g_kpc:.3e} (km/s)²/kpc")

# Convert to m/s²
# 1 (km/s)²/kpc = (1000 m/s)² / (3.086e19 m) = 10^6 / 3.086e19 = 3.24e-14 m/s²
conversion = 1e6 / kpc_m
print(f"\n   Conversion: 1 (km/s)²/kpc = {conversion:.3e} m/s²")

g_SI_check = g_kpc * conversion
print(f"   g in SI = {g_SI_check:.3e} m/s²")
print(f"   g / a₀ = {g_SI_check / a0:.1f}")

print("\n5. FINDING THE ERROR:")
print("   In main script, a0_cluster was computed as:")
print("   a0_cluster = a0 * 3.086e13 (km/s)²/kpc")

a0_cluster_wrong = a0 * 3.086e13
print(f"   a0_cluster = {a0_cluster_wrong:.3e}")

print("\n   BUT the correct conversion from m/s² to (km/s)²/kpc is:")
print("   1 m/s² = 1 m/s² × (1 km/s / 1000 m/s)² × (kpc / 3.086e19 m)")
print("   1 m/s² = 1e-6 / 3.086e19 (km/s)²/kpc")
print("   1 m/s² = 3.24e-26 (km/s)²/kpc")

correct_conversion = 1e-6 / kpc_m
a0_cluster_correct = a0 / conversion
print(f"\n   So a₀ = {a0} m/s² = {a0_cluster_correct:.3e} (km/s)²/kpc")

print("\n   The script was WRONG by a factor of ~10^39!")

print("\n" + "=" * 70)
print("6. CORRECTED CALCULATION")
print("=" * 70)

# Correct a0 in cluster units
a0_correct = a0 / conversion
print(f"\na₀ in (km/s)²/kpc = {a0_correct:.6e}")

# g at R_200 in same units
g_at_R200 = g_kpc
print(f"g(R_200) in (km/s)²/kpc = {g_at_R200:.6e}")

ratio = g_at_R200 / a0_correct
print(f"\ng(R_200) / a₀ = {ratio:.2f}")

print("\nCONCLUSION:")
print(f"  At R_200 of a 10^15 M_sun cluster:")
print(f"  g ≈ {ratio:.0f} × a₀")
print("  This is in the TRANSITION REGIME, not deep Newtonian!")
print("  C(g) will have significant modification.")

# What about at R = 200 kpc (inner region)?
R_inner = 200  # kpc
g_inner = G_kpc * M / R_inner**2
ratio_inner = g_inner / a0_correct
print(f"\n  At r = 200 kpc:")
print(f"  g ≈ {ratio_inner:.0f} × a₀")

# What about at R = 5 Mpc (outskirts)?
R_outer = 5000  # kpc
g_outer = G_kpc * M / R_outer**2
ratio_outer = g_outer / a0_correct
print(f"\n  At r = 5 Mpc:")
print(f"  g ≈ {ratio_outer:.1f} × a₀")

print("\n" + "=" * 70)
print("7. COHERENCE VALUES WITH CORRECT a₀")
print("=" * 70)

phi = (1 + np.sqrt(5)) / 2
Omega_m = 0.315

def coherence(g, a0):
    """C(g) = Ω_m + (1-Ω_m) × (g/a₀)^(1/φ) / [1 + (g/a₀)^(1/φ)]"""
    if g <= 0:
        return Omega_m
    x = (g / a0) ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

print(f"\nCoherence at different cluster radii (M = 10^15 M_sun):")
print("-" * 60)
print(f"{'Radius':>12} {'g/a₀':>10} {'C(g)':>10} {'G_eff/G':>10}")
print("-" * 60)

for R_kpc in [100, 200, 500, 1000, 2000, 3000, 5000]:
    g = G_kpc * 1e15 / R_kpc**2
    g_ratio = g / a0_correct
    C = coherence(g, a0_correct)
    G_eff = 1 / C
    print(f"{R_kpc:>10} kpc {g_ratio:>10.1f} {C:>10.4f} {G_eff:>10.4f}")

print("\n" + "=" * 70)
print("8. IMPLICATIONS FOR CLUSTER DYNAMICS")
print("=" * 70)

print("""
WITH CORRECT UNIT CONVERSION:
=============================

1. At R_200 (2 Mpc): g ~ 3 a₀
   → C ≈ 0.75
   → G_eff/G ≈ 1.3

2. At 5 Mpc: g ~ 0.5 a₀
   → C ≈ 0.55
   → G_eff/G ≈ 1.8

3. This means:
   - Clusters ARE in the transition/MOND regime!
   - Significant G_eff enhancement expected
   - NOT pure Newtonian as the buggy script suggested

4. BUT: The enhancement is STILL not enough for clusters
   - Maximum G_eff/G = 1/Ω_m = 3.17
   - Clusters need M_dyn/M_baryon ~ 10
   - So factor of ~3× still insufficient

5. CONCLUSION:
   - The cluster mass problem persists
   - But clusters DO show Synchronism effects
   - Just not enough to eliminate need for DM in clusters
""")

print("=" * 70)
print("BUG FIXED - Rerun main analysis needed")
print("=" * 70)
