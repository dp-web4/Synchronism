#!/usr/bin/env python3
"""
Session 686: verify the load-bearing claims of
`Research/proposals/gamma_ncorr_sign_inversion_sharpness.md`
(maintainer + explorer adjudication, 2026-06-06):

  1. Galaxy fixed point: 2/sqrt(N) = 2*sqrt(N) at N=1 (trivial but worth
     showing the calibrated SPARC inputs sit at the swap-identity point).

  2. Effect of the proposed flip (gamma = 2*sqrt(N_corr)) on the C(rho)
     shape for BCS / BEC / collective systems: does the transition move
     from C ~ 0 everywhere (current) to sharp near a low-rho transition
     point (proposed)?

  3. Apply the flip to S676's ladder (electron -> molecule -> nanoparticle
     -> BCS -> BEC at fixed rho/rho_crit). Under the original formula
     gamma = 2/sqrt(N_corr), C was MONOTONE DECREASING in N_corr (S676
     "naming inversion"). Under the flip, C should be MONOTONE INCREASING
     in N_corr -- which is what fixing a sign error looks like.

These checks are all under Reading A (universal coherence scalar, presets
comparable on one C axis). Under Reading B (density-response function,
gamma as inverse effective T) neither S676 nor this proposal's claim is
even statable -- the foundational fork (A vs B) is the actual research
question.
"""
import numpy as np


def C(rho_over_rhocrit, gamma):
    return np.tanh(gamma * np.log(rho_over_rhocrit + 1.0))


# ----------------------------------------------------------------------------
# (1) Galaxy fixed point
# ----------------------------------------------------------------------------
print("=" * 74)
print("(1) Galaxy fixed point: N_corr = 1 -> 2/sqrt(1) = 2 = 2*sqrt(1)")
print("=" * 74)
for label in ["original 2/sqrt(N)", "flip 2*sqrt(N)"]:
    if "original" in label:
        gamma_galaxy = 2.0 / np.sqrt(1)
    else:
        gamma_galaxy = 2.0 * np.sqrt(1)
    print(f"   {label}: gamma(N=1) = {gamma_galaxy:.6f}")
print("   Both forms give gamma = 2 at the galaxy regime. SPARC fits,")
print("   rho_crit = A * V_flat^2, and any other N_corr=1 calibration are")
print("   IDENTICAL under either sign. The flip is free at galaxy scale.")

# ----------------------------------------------------------------------------
# (2) BCS / BEC shape under the flip
# ----------------------------------------------------------------------------
print()
print("=" * 74)
print("(2) BCS / BEC transition shape under the proposed flip")
print("=" * 74)
N_BCS = 1e7
gamma_orig_BCS = 2.0 / np.sqrt(N_BCS)
gamma_flip_BCS = 2.0 * np.sqrt(N_BCS)
print(f"   N_corr_BCS = {N_BCS:.0e}")
print(f"   gamma_original = 2/sqrt(N) = {gamma_orig_BCS:.3e}")
print(f"   gamma_flip     = 2*sqrt(N) = {gamma_flip_BCS:.3e}")
print()
print("   Sample the tanh shape at rho/rho_crit = {1e-6, 1e-5, 1e-4, 0.1, 1, 10}:")
print(f"   {'rho/rho_crit':>14} {'C original':>14} {'C flip':>14}")
for r in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0]:
    C_orig = C(r, gamma_orig_BCS)
    C_flip = C(r, gamma_flip_BCS)
    print(f"   {r:>14.0e} {C_orig:>14.6f} {C_flip:>14.6f}")
print()
print("   Reading:")
print("   - Original: C ~ 0 across the whole rho range (BCS pinned at C~0).")
print("   - Flip: C transitions sharply between rho/rho_crit ~ 1e-7 and 1e-4,")
print("     saturating to ~1 thereafter. The proposed flip produces the qualitative")
print("     'sharp transition at a low rho threshold' that BCS-with-a-real-Tc")
print("     demands, while preserving C ~ 1 above it. The half-transition point")
half_point = (np.exp(np.arctanh(0.5) / gamma_flip_BCS) - 1.0)
print(f"     under the flip is at rho/rho_crit ~ {half_point:.2e} (very low rho).")

# ----------------------------------------------------------------------------
# (3) S676 ladder under the flip
# ----------------------------------------------------------------------------
print()
print("=" * 74)
print("(3) S676 ladder under original vs flipped formula  (at rho = rho_crit)")
print("=" * 74)
print("   Quantum-coherence ladder: increasing N_corr should correspond to")
print("   INCREASING quantum coherence / ODLRO / collective ordering.")
print()
print(f"   {'System':<26} {'N_corr':>9} {'C orig':>10} {'C flip':>10}")
ladder = [
    ("Lone electron", 1),
    ("Diatomic molecule", 2),
    ("Small macromolecule", 100),
    ("Mesoscale nanoparticle", 1000),
    ("BCS superconductor", int(1e7)),
    ("BEC", int(1e8)),
]
for label, N in ladder:
    gamma_orig = 2.0 / np.sqrt(N)
    gamma_flip = 2.0 * np.sqrt(N)
    C_orig = C(1.0, gamma_orig)          # rho = rho_crit -> ratio = 1
    C_flip = C(1.0, gamma_flip)
    print(f"   {label:<26} {N:>9d} {C_orig:>10.4f} {C_flip:>10.4f}")
print()
print("   Original column: monotone DECREASING with N_corr. S676 'naming inversion'")
print("   -- the framework's namesake variable is anti-correlated with quantum")
print("   coherence / synchronization across the ladder.")
print()
print("   Flip column: monotone INCREASING with N_corr. The flip qualitatively")
print("   reverses the inversion -- which is what fixing a sign error does.")
print("   S676's 'C anti-correlated with quantum coherence' becomes 'C correlated")
print("   with quantum coherence' under the flip, at zero cost to galaxy fits.")
print()
print("   But both readings assume the cross-system C ladder is even well-defined")
print("   (Reading A). Under Reading B (C is system-specific density-response,")
print("   gamma is inverse effective T), neither S676 nor this proposal's claim")
print("   type-checks -- you cannot ladder C across systems at all. The flip and")
print("   the original formula are then equally defensible at zero empirical cost,")
print("   but the entire family of presets / visualizer / cross-system comparisons")
print("   that the framework's tooling produces becomes invalid.")

# ----------------------------------------------------------------------------
# (4) The foundational tension surfaced
# ----------------------------------------------------------------------------
print()
print("=" * 74)
print("(4) The foundational tension: framework's tools vs framework's formula")
print("=" * 74)
print("   - The visualizer, gamma-calculator, and preset table all commit to")
print("     Reading A (single C axis, presets comparable). This is the tooling.")
print("   - Reading A produces TWO simultaneous sign inversions (S676 magnitude;")
print("     this proposal sharpness). The proposed flip fixes both at zero cost")
print("     to galaxy fits, but has no first-principles derivation.")
print("   - Reading B (density-response with gamma as inverse effective T) saves")
print("     the formula's directionality but invalidates the cross-system tooling.")
print("   - There is no choice that keeps both the tooling and the original")
print("     formula direction. The framework is in a tools-vs-formula tension at")
print("     the C ontology layer, not at the parameter layer.")
print()
print("   The actual research question is NOT 'which sign' but 'commit C to A or")
print("   B' -- a frame question, identical to the explorer's bottom line. The")
print("   research-repo contribution: under either commitment, the framework")
print("   loses something already deployed (tooling vs formula direction).")
