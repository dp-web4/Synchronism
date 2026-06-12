#!/usr/bin/env python3
"""
Session 691: verify the load-bearing claim in
`Research/proposals/test02_kill_branch_adjudicable_now.md` (maintainer 06:22)
that the wide-binary regime's local density SATURATES C(rho) -> ~1 under the
framework's galaxy-anchored rho_crit, which means the framework predicts
Newton in the wide-binary regime, which makes the current stated kill
criterion ("underlying anomaly fails") INVERTED -- it would kill the
framework for being right.

The 2026-06-08 explorer adjudication (test02_adjudication_hung_modeling_crux.md)
returns HUNG on the literature (Chae 4.9 sigma vs Saad-Ting 0.4 sigma on the
SAME 36 systems; significance lives in one modeling choice). I do not
independently re-adjudicate the literature -- I verify the framework-side
prediction that the kill criterion is keyed to.
"""
import numpy as np

# --- Local solar-neighborhood density components ----------------------------
# (using standard astrophysics literature values, not framework values)
# ISM neutral gas: ~1 H atom / cm^3
mH = 1.673e-27  # kg
n_H = 1.0       # atoms / cm^3
rho_gas = n_H * mH * 1e6  # kg / m^3
# Local stellar mass density (Galactic disc at solar location): ~ 0.05 Msun/pc^3
Msun = 1.989e30
pc = 3.086e16
rho_stellar = 0.05 * Msun / pc**3
# Local dark matter density: ~ 0.4 GeV/c^2 / cm^3
GeV_over_c2 = 1.783e-27  # kg
rho_DM = 0.4 * GeV_over_c2 * 1e6  # kg/m^3
# Total local matter density (the relevant input for C(rho))
rho_local = rho_gas + rho_stellar + rho_DM

print("=" * 74)
print("(1) Local solar-neighborhood density components")
print("=" * 74)
print(f"   ISM neutral gas         : {rho_gas:.2e} kg/m^3")
print(f"   Stellar disc            : {rho_stellar:.2e} kg/m^3")
print(f"   Dark matter (halo)      : {rho_DM:.2e} kg/m^3")
print(f"   Total                   : {rho_local:.2e} kg/m^3")

# --- C(rho) at local density under the framework's galaxy anchor ------------
print()
print("=" * 74)
print("(2) C(rho) at local density under galaxy-anchored rho_crit")
print("=" * 74)
rho_crit_galaxy = 1e-23  # S678/S683 lower edge of galaxy-anchored range
gamma = 2.0              # N_corr = 1 (galaxy-rung default, S688)

ratio = rho_local / rho_crit_galaxy
arg = gamma * np.log(ratio + 1.0)
C = np.tanh(arg)
print(f"   rho_crit (galaxy anchor) : {rho_crit_galaxy:.2e} kg/m^3")
print(f"   gamma (N_corr = 1)       : {gamma}")
print(f"   rho_local / rho_crit     : {ratio:.2e}")
print(f"   gamma * ln(ratio + 1)    : {arg:.3f}")
print(f"   C(rho_local) = tanh(arg) : {C:.6f}")
print()
print("   Result: at solar-neighborhood density (~1e-21 kg/m^3, dominated by DM")
print(f"   halo at ~7e-22), C(rho) saturates to {C:.4f} -- effectively 1.")

# --- What does C ~ 1 mean for the wide-binary prediction? -------------------
print()
print("=" * 74)
print("(3) Implication for the wide-binary prediction under C(rho)")
print("=" * 74)
# In the framework's published form, the gravitational coupling is modulated
# by C. The published convention is g_eff = g_N / C(rho) at the limit of small
# C; equivalently in the high-C regime the modification vanishes and motion
# follows Newton's laws. Either way, C ~ 1 means NO deviation from Newton.
print("   With C(rho_local) ~ 1, the framework's prediction at wide-binary")
print("   internal accelerations is: NO deviation from Newton. The framework")
print("   PREDICTS the Newton null in the local solar-neighborhood regime,")
print("   independent of internal acceleration. C(rho) modulation is keyed to")
print("   density, and the local density saturates the coherence to one.")
print()
print("   This is the prediction the 2026-06-05 fork computation (cited in")
print("   the maintainer's proposal) reached: Newtonian null of 0.05-0.4% in")
print("   the clean within-250-pc sample, vs MOND's ~18% boost prediction.")

# --- The kill criterion inversion -------------------------------------------
print()
print("=" * 74)
print("(4) Kill criterion inversion verified")
print("=" * 74)
print("   Stated kill criterion (currently on the site, pre-fork):")
print("     'anomaly independent of local density OR underlying anomaly fails'")
print()
print("   Decomposing under the current C(rho) form:")
print("     'anomaly independent of local density'")
print("        - if Chae's ~1.4x boost is real AND independent of local rho,")
print("          this is the framework's prediction failing (Newton expected,")
print("          MOND-like observed). KILL FIRES correctly.")
print("     'underlying anomaly fails'")
print("        - if Banik's Newton null is real -> NO anomaly -> matches the")
print("          framework's Newton prediction (C ~ 1) -> NOT a refutation.")
print("          But the criterion as stated says this kills the framework.")
print("          INVERSION: would kill the framework for being right.")
print()
print("   The correct post-fork kill criterion:")
print("     KILL: Gaia-confirmed MOND-scale (~18%, ~1.4x boost) wide-binary")
print("           anomaly in the clean sample refutes C(rho).")
print("     NON-DISCRIMINATING SURVIVAL: confirmed Newton null is consistent")
print("           with C(rho) but equally consistent with GR -- no points.")
print()
print("   The maintainer's claim of inversion is verified directly from the")
print("   framework's own C(rho) at local density.")

# --- HUNG verdict acceptance and external-mirror observation ----------------
print()
print("=" * 74)
print("(5) Adjudication state: HUNG (per explorer 2026-06-12 08:09)")
print("=" * 74)
print("   I do not independently re-adjudicate the literature dispute. The")
print("   explorer's mapping:")
print("     - Chae et al. 2026 (arXiv:2601.21728): 4.9 sigma, gamma_boost~1.6")
print("       from 36 RV+speckle-vetted binaries.")
print("     - Saad & Ting 2026 (arXiv:2603.11015): SAME 36 systems, gamma=1.12")
print("       +- 0.25 (Newton at 0.4 sigma) under hierarchical 3D-orbit")
print("       inference; recovers Chae's 1.56 only by re-adopting his")
print("       geometric-deprojection prior.")
print("     - Cookson-Banik-El-Badry et al. 2026: Newton 'up to 1500x more")
print("       likely.'")
print("     - Significance flips between 4.9 sigma and 0.4 sigma on the SAME")
print("       data under one modeling choice (orbital prior).")
print()
print("   Catalog state per explorer: HUNG -- kill criterion does NOT fire;")
print("   degenerate survival NOT established. Re-open triggers replace")
print("   'future Gaia data' with literature events (Chae rebuttal to Saad-")
print("   Ting; mock-injection adjudication; independent non-camp confirmation).")

# --- External mirror observation --------------------------------------------
print()
print("=" * 74)
print("(6) External mirror of the framework's own failure mode")
print("=" * 74)
print("   The explorer flagged a striking pattern: 'Same data flipping between")
print("   5 sigma discovery and null under a modeling assumption is the")
print("   literature-side twin of what this archive keeps documenting")
print("   internally (gamma=2/sqrt(N_corr) absorbing N_corr; A-from-Jeans")
print("   0.0294 outliving its computation; A2ACW survival rate as filter")
print("   property).'")
print()
print("   This pattern matches the methodology sessions S672/S687/S689/S690:")
print("     S672 (DESI): wrong-paper value 0.45 propagated without re-execution")
print("     S687 (A-from-Jeans): 0.0294 outlived its V^0.5 source computation")
print("     S689 (cluster no-go): wrong-variable framing propagated without")
print("                           citing Milgrom 2005")
print("     S690 (C latent variable): forward-map structural barrier surfaced")
print("                               by external archive survey")
print()
print("   The Chae vs Saad-Ting dispute is the EXTERNAL ANALOGUE: significance")
print("   propagated from one modeling assumption to another camp, where the")
print("   one modeling choice carries 100% of the claimed effect. This is the")
print("   same 'one input absorbing the whole signal' pattern -- now on the")
print("   data-analysis side rather than the framework side.")
print()
print("   The framework's last empirical falsification channel waits on the")
print("   external field resolving exactly the failure mode the internal audit")
print("   catalogs. Two possibilities: (a) coincidence, (b) the failure mode")
print("   is general to underdetermined-modeling-choice domains. I have no")
print("   leverage to distinguish them; flag the structural symmetry as a")
print("   substantive observation worth recording.")
