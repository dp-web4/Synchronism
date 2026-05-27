#!/usr/bin/env python3
"""
Session 673: TEST-15 (GW170817 coherence-coupling) has no derived amplitude.
Adjudicates the proposal gw170817_test15_resolution_no_derived_amplitude.md by
re-grounding against the two primary sources (verified this session):

  Session 59 (Session59_GW_Coherence_Theory.md), VERIFIED quotes:
    line 206: c_g/c = 1 + alpha*(1 - <C>_LOS)        [Case-1 propagation claim]
    line 259: c_g/c = 1 + alpha*(1 - C) = 1 + alpha*f_DM
    line 149: "If alpha ~ 10^-15 (from GW170817 constraint)"   [READ OFF, not derived]
    line 150: "Dc_g/c ~ alpha * 0.9 ~ 10^-15"                  [<1-C> ~ 0.9, C_avg ~ 0.1]
    line 312: GW170817: |c_g - c|/c < 3e-15

  Session 642 (Session642_..._Field_Or_Parameterization.md), VERIFIED:
    "no Lagrangian, no action, no equation of motion for C(rho)... GW170817 does
     not directly constrain Synchronism... Case 3 - a parameterization."

These two are mutually exclusive (Session 59 IS the propagation claim Session 642
says does not exist). This script shows TEST-15 is non-discriminating either way,
and quantifies the "natural scale is dead" coupling hierarchy.
"""

# GW170817 bound (Session 59 line 312)
gw_bound = 3e-15                      # |c_g - c|/c < 3e-15
one_minus_C_path = 0.9               # Session 59's own LOS estimate (C_avg ~ 0.1)

# (A) The amplitude alpha is bounded, not derived
alpha_max = gw_bound / one_minus_C_path
print("=" * 72)
print("(A) alpha is read off GW170817, not derived")
print("=" * 72)
print(f"   c_g/c - 1 = alpha * (1 - <C>_LOS); Session 59 uses (1-<C>) ~ {one_minus_C_path}")
print(f"   GW170817: |c_g/c - 1| < {gw_bound:g}  ->  alpha < {alpha_max:.2g}")
print("   Session 59 line 149 literally sets 'alpha ~ 10^-15 (from GW170817 constraint)'.")
print("   The single free parameter is read off the very experiment TEST-15 performs.")
print("   => no independent prediction, no derived floor. NOT a Tier-1 quantity")
print("      (S670): there is no derived central object, only a data-fit upper bound.")

# (B) The natural scale is dead: same (1-C) drives DM and GW, but couplings must
#     differ by >=15 orders of magnitude
print()
print("=" * 72)
print("(B) The natural scale is dead -- required coupling hierarchy")
print("=" * 72)
# Dark-matter mechanism: (1-C) ~ O(1) coupling O(1) -> O(1) mass discrepancy
dm_coupling = 1.0                    # O(1) by construction (DM effect is order unity)
# If GW coupling alpha tracked the SAME coherence physics, alpha ~ O(1):
alpha_natural = 1.0
hierarchy = alpha_natural / alpha_max
print(f"   Session 59 line 259: c_g/c = 1 + alpha*(1-C) = 1 + alpha*f_DM")
print(f"   The DM mechanism uses (1-C)~O(1) with coupling ~O(1) to get O(1) mass excess.")
print(f"   If alpha tracked that same coherence physics: alpha ~ O(1) = {alpha_natural}.")
print(f"   GW170817 forces alpha < {alpha_max:.2g}.")
print(f"   Required decoupling: coherence couples to GW propagation >= {hierarchy:.0e}x")
print(f"   more weakly than to galactic dynamics -- an unstated, underived hierarchy.")
print("   alpha ~ O(1) is excluded by ~15 orders of magnitude.")

# (C) The dichotomy
print()
print("=" * 72)
print("(C) TEST-15 is non-discriminating either way")
print("=" * 72)
print("   Branch alpha ~ O(1) (natural scale, tracks DM coupling): REFUTED by GW170817")
print("                                                            at ~15 OOM.")
print("   Branch alpha read off data (<= 3e-15): GR-equivalent in the allowed range;")
print("        no derived floor -> unfalsifiable as a POSITIVE prediction.")
print("   Either branch: TEST-15 discriminates nothing. Novel discriminating-test")
print("   count unchanged at 0.")

print()
print("=" * 72)
print("(D) Internal contradiction (both primary sources verified this session)")
print("=" * 72)
print("   Session 59 : Case-1 propagation claim, TEST-15 a 'novel test'.")
print("   Session 642: Case-3, 'framework makes NO GW propagation claim'.")
print("   Mutually exclusive. Session 642 (later, more careful) is the honest position,")
print("   but it overlooked that Session 59 had already filed the Case-1 claim. Holding")
print("   Case-3 REQUIRES removing TEST-15 from the catalog -- Case-3 IS the statement")
print("   that TEST-15's prediction does not exist.")
print()
print("   GW sector joins the import-of-predictive-content pattern (S671): the only")
print("   GW parameter (alpha) is imported from GW170817, as the i was imported for QM")
print("   (S666), theta_D for chemistry (S669), and the MOND interpolating function for")
print("   galaxies (S661). The framework's GW 'prediction' is the GW170817 bound,")
print("   relabeled.")
