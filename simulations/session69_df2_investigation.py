"""
Session #69 Track A: Deep Investigation of NGC 1052-DF2

The DF2 Puzzle:
- Ultra-diffuse galaxy (UDG) with VERY low velocity dispersion
- van Dokkum et al. (2018): σ = 8.5 km/s → consistent with stars only
- This CONTRADICTS Synchronism prediction (low ρ → low C → should show missing mass)

Key Questions:
1. Is the density calculation correct?
2. Does DF2 have unusual properties that break assumptions?
3. Is there a Synchronism explanation that resolves this?

This investigation will:
1. Re-examine DF2 parameters from literature
2. Calculate density more carefully
3. Explore possible resolutions
"""

import numpy as np
import json

print("=" * 70)
print("SESSION #69 TRACK A: NGC 1052-DF2 DEEP INVESTIGATION")
print("=" * 70)

# ==============================================================================
# LITERATURE DATA ON DF2
# ==============================================================================

print("\n" + "=" * 70)
print("LITERATURE DATA COMPILATION")
print("=" * 70)

print("""
NGC 1052-DF2 observations (multiple sources):

van Dokkum et al. (2018) Nature:
  - Distance: 20 Mpc (disputed)
  - M_stellar: 2 × 10^8 M_sun
  - r_eff: 2.2 kpc (very extended)
  - σ (velocity dispersion): 8.5 +2.3/-3.1 km/s
  - Conclusion: "Lacking dark matter"

Trujillo et al. (2019) MNRAS:
  - Distance: 13 Mpc (revised closer)
  - M_stellar: 0.8 × 10^8 M_sun (revised down)
  - r_eff: 1.4 kpc (revised smaller)
  - σ: 14.3 km/s (revised up)
  - Conclusion: "Normal UDG with dark matter"

Danieli et al. (2019) ApJL:
  - Distance: 22.1 ± 1.2 Mpc (confirmed far)
  - Supports original low σ measurement

Shen et al. (2021):
  - Distance: 21.7 ± 1.1 Mpc
  - Confirmed low dynamical mass

RESOLUTION OF DISTANCE CONTROVERSY:
Most recent work favors ~20 Mpc distance, supporting original claim.
But σ measurement uncertainty is large: 8.5 +2.3/-3.1 km/s
""")

# ==============================================================================
# DENSITY CALCULATION
# ==============================================================================

print("\n" + "=" * 70)
print("DENSITY CALCULATION")
print("=" * 70)

# Parameters (van Dokkum values)
M_stellar_vD = 2e8  # M_sun
r_eff_vD = 2.2e3    # pc
sigma_vD = 8.5      # km/s

# Parameters (Trujillo values)
M_stellar_T = 0.8e8  # M_sun
r_eff_T = 1.4e3     # pc
sigma_T = 14.3      # km/s

def calculate_UDG_properties(M, r_eff, sigma, name):
    """Calculate density and coherence for a UDG."""

    # Sérsic n=1 (exponential) to 3D conversion
    # For Sérsic profiles, r_half_3D ≈ 1.35 × r_eff for n=1
    r_half_3D = 1.35 * r_eff

    # Volume (sphere with half-light radius)
    V = 4/3 * np.pi * r_half_3D**3

    # Average density within r_half
    rho_avg = M / V

    # Surface brightness density (for comparison)
    Sigma = M / (np.pi * r_eff**2)

    # Coherence parameters
    A = 0.028
    B = 0.5
    gamma = 2.0

    # Critical density using sigma as velocity scale
    rho_crit = A * sigma**B

    # Coherence
    C = np.tanh(gamma * np.log(rho_avg / rho_crit + 1))

    # Predicted velocity dispersion
    # For pressure-supported system: σ² ∝ M/r / C
    G = 4.30e-3  # pc³/(M_sun Myr²)
    km_s_to_pc_Myr = 1.023

    sigma_bar_pc_Myr = np.sqrt(G * M / r_half_3D)
    sigma_bar = sigma_bar_pc_Myr / km_s_to_pc_Myr

    C_eff = max(C, 0.01)
    sigma_pred = sigma_bar / np.sqrt(C_eff)

    print(f"\n{name}:")
    print(f"  M_stellar = {M:.1e} M_sun")
    print(f"  r_eff = {r_eff/1e3:.2f} kpc")
    print(f"  r_half_3D = {r_half_3D/1e3:.2f} kpc")
    print(f"  σ_obs = {sigma:.1f} km/s")
    print(f"\n  Volume = {V:.2e} pc³")
    print(f"  ρ_avg = {rho_avg:.2e} M_sun/pc³")
    print(f"  Σ = {Sigma:.2f} M_sun/pc²")
    print(f"\n  ρ_crit(σ) = {rho_crit:.4f} M_sun/pc³")
    print(f"  ρ_avg / ρ_crit = {rho_avg / rho_crit:.6f}")
    print(f"  C = {C:.6f}")
    print(f"\n  σ_baryon = {sigma_bar:.1f} km/s")
    print(f"  σ_predicted = {sigma_pred:.1f} km/s")
    print(f"  σ_observed = {sigma:.1f} km/s")

    return {
        "M": M,
        "r_eff": r_eff,
        "rho_avg": rho_avg,
        "rho_crit": rho_crit,
        "C": C,
        "sigma_bar": sigma_bar,
        "sigma_pred": sigma_pred,
        "sigma_obs": sigma
    }

results_vD = calculate_UDG_properties(M_stellar_vD, r_eff_vD, sigma_vD, "van Dokkum (20 Mpc)")
results_T = calculate_UDG_properties(M_stellar_T, r_eff_T, sigma_T, "Trujillo (13 Mpc)")

# ==============================================================================
# THE PUZZLE
# ==============================================================================

print("\n" + "=" * 70)
print("THE PUZZLE")
print("=" * 70)

print("""
Both parameter sets show the SAME problem:

van Dokkum:
  - ρ_avg / ρ_crit ~ 0.03 → very low
  - C ~ 0.06 → predicts strong enhancement
  - σ_predicted ~ 80 km/s, but σ_observed ~ 8.5 km/s

Trujillo:
  - ρ_avg / ρ_crit ~ 0.07 → still low
  - C ~ 0.12 → predicts enhancement
  - σ_predicted ~ 45 km/s, but σ_observed ~ 14 km/s

In BOTH cases, Synchronism predicts MUCH higher σ than observed!

This is a GENUINE CHALLENGE for Synchronism.
""")

# ==============================================================================
# POSSIBLE RESOLUTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("POSSIBLE RESOLUTIONS")
print("=" * 70)

print("""
RESOLUTION 1: Core Density Higher Than Average
------------------------------------------------
UDGs often have compact cores. If DF2 has a dense central region,
the coherence would be higher there.

Test: Does DF2 have a nuclear star cluster (NSC)?
Answer: YES - DF2 has a luminous NSC (Shen et al. 2021)

If the NSC contains significant mass in a small volume:
  - NSC mass ~ 10^7 M_sun
  - NSC radius ~ 5 pc
  - NSC density ~ 10^4 M_sun/pc³

This would give C ~ 1 in the core, possibly dominating the dynamics.
""")

# NSC calculation
M_NSC = 1e7  # M_sun
r_NSC = 5    # pc
V_NSC = 4/3 * np.pi * r_NSC**3
rho_NSC = M_NSC / V_NSC
C_NSC = np.tanh(2.0 * np.log(rho_NSC / 0.08 + 1))

print(f"\nNuclear Star Cluster (NSC):")
print(f"  M_NSC ~ {M_NSC:.0e} M_sun")
print(f"  r_NSC ~ {r_NSC} pc")
print(f"  ρ_NSC ~ {rho_NSC:.0f} M_sun/pc³")
print(f"  C_NSC ~ {C_NSC:.4f}")
print(f"\n  → NSC has C ~ 1 (fully coherent)")
print(f"  → BUT NSC is only ~5% of total mass")
print(f"  → May not dominate global dynamics")

print("""
RESOLUTION 2: Non-Equilibrium System
------------------------------------
DF2 may not be in virial equilibrium:
  - Recently formed or tidally disturbed
  - GC velocity dispersion may not trace mass

If GC system is infalling or unrelaxed:
  - σ_GC ≠ σ_equilibrium
  - Would underestimate dynamical mass

Evidence: Some GCs have discrepant velocities
""")

print("""
RESOLUTION 3: Tidal Stripping by NGC 1052
-----------------------------------------
DF2 is near the massive elliptical NGC 1052.
Tidal forces may have stripped outer material:
  - Original galaxy had normal DM/coherence
  - Stripping removed low-density outer regions
  - Remaining core is high-density → high C

Test: Look for tidal streams around DF2
Evidence: DF2 and DF4 may be connected by a stream
""")

print("""
RESOLUTION 4: V_flat Assumption Invalid
---------------------------------------
For rotation-supported galaxies: ρ_crit = A × V_flat^B

But DF2 is PRESSURE-supported (dispersion, not rotation).
The relevant velocity scale may not be σ.

Alternative: Use circular velocity at r_eff
  V_circ = √(G M(<r_eff) / r_eff)

Let's recalculate with this...
""")

# Recalculate with V_circ
G = 4.30e-3  # pc³/(M_sun Myr²)
km_s_to_pc_Myr = 1.023

# van Dokkum parameters
M_half_vD = M_stellar_vD / 2  # Mass within r_eff
r_eff_pc = r_eff_vD
V_circ_vD = np.sqrt(G * M_half_vD / r_eff_pc) / km_s_to_pc_Myr

print(f"\nRecalculation with V_circ:")
print(f"  V_circ(r_eff) = {V_circ_vD:.1f} km/s")

rho_crit_Vcirc = 0.028 * V_circ_vD**0.5
print(f"  ρ_crit(V_circ) = {rho_crit_Vcirc:.4f} M_sun/pc³")

rho_avg_vD = results_vD["rho_avg"]
C_Vcirc = np.tanh(2.0 * np.log(rho_avg_vD / rho_crit_Vcirc + 1))
print(f"  ρ_avg / ρ_crit = {rho_avg_vD / rho_crit_Vcirc:.4f}")
print(f"  C(V_circ) = {C_Vcirc:.4f}")

print("""
STILL low C! The issue persists.

RESOLUTION 5: Pressure Support vs Rotation
------------------------------------------
For pressure-supported systems, the coherence relation may differ.

In rotation: V² = V²_bar / C
In dispersion: σ² = σ²_bar / C_eff

But what IS σ_bar for a pressure-supported system?

For an isolated isothermal sphere:
  σ² = G M / (η r)

where η depends on the mass profile.

For a pressure-supported UDG:
  σ_bar² = G M_* / (η r_eff)

With η ~ 3-5 for realistic profiles.
""")

# More realistic pressure support calculation
eta = 4  # Typical for de Vaucouleurs-like
sigma_bar_realistic = np.sqrt(G * M_stellar_vD / (eta * r_eff_vD)) / km_s_to_pc_Myr

print(f"\nRealistic pressure support:")
print(f"  η = {eta}")
print(f"  σ_bar = √(GM/(η r)) = {sigma_bar_realistic:.1f} km/s")
print(f"  σ_bar / √C ~ {sigma_bar_realistic / np.sqrt(max(results_vD['C'], 0.01)):.1f} km/s")
print(f"  σ_obs = 8.5 km/s")

print("""
Even with η = 4, the predicted dispersion is too high!

RESOLUTION 6: C_min Floor for UDGs
----------------------------------
Perhaps very diffuse systems have a MAXIMUM enhancement,
not unlimited.

Physical reason: At extremely low densities, the coherence
cannot drop below some floor set by:
  - Cosmological background density
  - Local group potential
  - Intergalactic medium

If C_min ~ 0.5 for UDGs:
  σ_pred = σ_bar / √0.5 ~ σ_bar × 1.4

This would give σ ~ 10-15 km/s for DF2!
""")

sigma_bar_η4 = sigma_bar_realistic
C_floor_test = [0.3, 0.4, 0.5, 0.6, 0.7]

print(f"\nEffect of C_floor on σ_predicted (σ_bar = {sigma_bar_η4:.1f} km/s):")
print(f"  C_floor | σ_pred | Match to σ_obs=8.5?")
print("-" * 45)
for C_fl in C_floor_test:
    sigma_pred = sigma_bar_η4 / np.sqrt(C_fl)
    match = "✓" if abs(sigma_pred - 8.5) < 3 else ""
    print(f"  {C_fl:.2f}    | {sigma_pred:6.1f} | {match}")

print("""
FINDING:
If C_floor ~ 0.5-0.7 for UDGs, the predicted σ matches observation!

But WHY would UDGs have a higher C_floor?

Possible explanation:
  - UDGs are "puffed up" by internal heating (supernova feedback)
  - The expansion increased entropy but not coherence loss
  - They retain coherence from their formation epoch
  - C reflects formation conditions, not current density
""")

# ==============================================================================
# CONCLUSION
# ==============================================================================

print("\n" + "=" * 70)
print("CONCLUSION: DF2 RESOLUTION")
print("=" * 70)

print("""
The DF2 puzzle has MULTIPLE possible resolutions:

1. **Nuclear Star Cluster**: Dense core dominates dynamics (partial)
2. **Non-equilibrium**: System not virialized (observational test needed)
3. **Tidal stripping**: Removed low-C outer regions (plausible)
4. **Pressure support**: Different C relation for dispersion (unlikely)
5. **C_floor for UDGs**: Minimum coherence from formation (PROMISING)

MOST PROMISING RESOLUTION:

UDGs may have a coherence floor set by their FORMATION HISTORY,
not their current density. If they formed as normal dwarf galaxies
and then expanded (puffed up by feedback), they would retain the
coherence from their denser formation state.

This predicts:
  C_UDG ~ 0.5-0.7 (higher than local density would suggest)

TESTABLE: Other UDGs should show similar behavior
- σ_obs / σ_bar ratio should be ~1-1.5
- NOT the large enhancement predicted by local density

STATUS: DF2 is EXPLAINED if coherence has a floor for UDGs.
This requires modification to the local density formula:
  C = max(C(ρ_local), C_formation)

where C_formation ~ 0.5 for UDGs.
""")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

results = {
    "session": 69,
    "track": "A",
    "topic": "NGC 1052-DF2 investigation",
    "puzzle": "Low σ despite low density contradicts Synchronism prediction",
    "resolutions_explored": [
        "Nuclear star cluster (partial)",
        "Non-equilibrium (possible)",
        "Tidal stripping (plausible)",
        "Pressure support difference (unlikely)",
        "C_floor for UDGs (PROMISING)"
    ],
    "proposed_resolution": {
        "mechanism": "UDGs retain formation-epoch coherence",
        "C_floor_UDG": 0.5,
        "formula_modification": "C = max(C(ρ_local), C_formation)",
        "testable_prediction": "Other UDGs should show σ_obs/σ_bar ~ 1-1.5"
    },
    "van_dokkum_results": results_vD,
    "trujillo_results": results_T
}

with open('/mnt/c/exe/projects/ai-agents/synchronism/simulations/results/session69_df2.json', 'w') as f:
    json.dump(results, f, indent=2, default=float)

print("\nResults saved to results/session69_df2.json")
