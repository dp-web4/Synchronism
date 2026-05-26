#!/usr/bin/env python3
"""
Session 670: A specificity test of the DEMOLITION method itself.

S662 measured the specificity of the A2ACW *vocabulary rule* (R1: "does a modern-
register restatement trigger prior art?") -> 0% on genuine discoveries. But the
method I actually used in S665/S666/S669 is different: REDUCTIO-TO-KNOWN-PHYSICS
("show the claim reduces to / relabels established physics"). Its specificity was
never measured. This script operationalizes and measures it.

Two tiers of reductio (criterion stated BEFORE application, applied uniformly):

  TIER 1 (exact): there exists a pre-existing named quantity/theorem Q such that the
    claim's central object equals Q by definition or exact derivation, with NO
    additional empirical content (no new confirmed prediction beyond Q's).
    Example fired correctly: gamma_phonon === 2T/theta_D (Debye temp), Delta_r=0 (S669).

  TIER 2 (weak): the claim's ingredients are all individually known, OR it resembles
    an existing framework. ("electrons+phonons+mean-field are all known", "C(rho)
    resembles MOND".)

Test set:
  - CONTROLS (must NOT be flagged -- these are genuine, confirmed discoveries):
    BCS, Noether's theorem, the Higgs mechanism, Dirac equation.
  - SYNCHRONISM verdicts I issued (should be flagged IF the method is sound):
    chemistry r=0.98 (S669), CFD irrotational (S665), dissipative-vs-unitary (S666),
    plus a resemblance-only verdict (C(rho) ~ Verlinde, S664) for contrast.

Specificity = fraction of genuine discoveries NOT flagged. High = good (few false
positives). The question: does each tier preserve the genuine discoveries?
"""

# case: (name, is_genuine_discovery,
#        tier1_exact_identity_to_known_Q, tier2_ingredients_known_or_resembles,
#        note)
cases = [
    # --- CONTROLS: genuine confirmed discoveries ---
    ("BCS superconductivity (1957)", True, False, True,
     "central object = energy gap Delta=2*hbar*wD*exp(-1/N0V); NOT identical to any "
     "pre-existing quantity; predicts NEW confirmed effects (isotope Tc~M^-1/2, "
     "2Delta/kTc=3.52, specific-heat jump). Ingredients (e-, phonon, mean-field) all "
     "known -> Tier2 FIRES (false positive). No exact identity -> Tier1 clean."),
    ("Noether's theorem (1918)", True, False, True,
     "symmetry<->conservation law; ingredients (Lagrangian mechanics, continuous "
     "groups) all known -> Tier2 fires. But the theorem is a NEW general result, not "
     "an identity relabeling one prior quantity -> Tier1 clean."),
    ("Higgs mechanism (1964)", True, False, True,
     "spontaneous symmetry breaking + gauge fields, both known -> Tier2 fires. "
     "Predicts a NEW particle (found 2012) and mass relations -> not a relabeling -> "
     "Tier1 clean."),
    ("Dirac equation (1928)", True, False, True,
     "combine SR + QM, both known -> Tier2 fires. Predicts antimatter + spin + g=2 "
     "(all NEW, confirmed) -> Tier1 clean."),
    # --- SYNCHRONISM verdicts I issued ---
    ("Sync chemistry r=0.98 (S669)", False, True, True,
     "central object gamma_phonon === 2T/theta_D BY DEFINITION; theta_D pre-exists "
     "(Debye 1912); Delta_r=0 to machine precision; NO new prediction -> Tier1 FIRES "
     "correctly (true positive)."),
    ("Sync CFD irrotational (S665)", False, True, True,
     "v=-g(I)grad(I) => curl v === 0 by vector identity (exact theorem); the claimed "
     "vortex phenomenology cannot exist -> Tier1 fires correctly."),
    ("Sync dissipative!=unitary (S666)", False, True, True,
     "real diffusion has real eigenvalues (decay), Schrodinger needs i (oscillation); "
     "the 'derivation' inserts i by hand -> structural identity, Tier1 fires."),
    ("Sync C(rho) ~ Verlinde (S664)", False, False, True,
     "RESEMBLANCE claim only (galaxy-regime reparametrization by reduction chain); no "
     "exact identity collapsing C(rho) to Verlinde -> Tier1 clean, Tier2 fires. This "
     "is the WEAK kind of verdict -- same tier that false-positives on BCS."),
]

def specificity(tier_index):
    # specificity = TN / (TN+FP) over the genuine-discovery controls
    controls = [c for c in cases if c[1]]            # is_genuine_discovery
    flagged = [c for c in controls if c[2+tier_index]]
    tn = len(controls) - len(flagged)
    return tn / len(controls), [c[0] for c in flagged]

def sensitivity(tier_index):
    # sensitivity = TP / (TP+FN) over the Synchronism verdicts (known reparametrizations)
    syncs = [c for c in cases if not c[1]]
    flagged = [c for c in syncs if c[2+tier_index]]
    return len(flagged) / len(syncs), [c[0] for c in syncs if not c[2+tier_index]]

print("=" * 74)
print("Specificity of the reductio method, by tier (on genuine-discovery controls)")
print("=" * 74)
for ti, label in [(0, "TIER 1 (exact identity to known quantity)"),
                  (1, "TIER 2 (ingredients known / resembles known framework)")]:
    spec, fps = specificity(ti)
    print(f"\n  {label}")
    print(f"    specificity (genuine discoveries NOT flagged) = {spec*100:.0f}%")
    if fps:
        print(f"    FALSE POSITIVES (real discoveries wrongly demoted): {', '.join(fps)}")
    else:
        print(f"    false positives: none")

print()
print("=" * 74)
print("Sensitivity of each tier on the Synchronism verdicts (known reparametrizations)")
print("=" * 74)
for ti, label in [(0, "TIER 1"), (1, "TIER 2")]:
    sens, missed = sensitivity(ti)
    print(f"  {label}: sensitivity = {sens*100:.0f}%"
          + (f"  (Tier1-clean: {', '.join(missed)})" if missed else ""))

print()
print("=" * 74)
print("VERDICT")
print("=" * 74)
print("  TIER 1 (exact identity/theorem, zero new content): specificity 100% -- fires")
print("  on NONE of BCS/Noether/Higgs/Dirac. These genuine discoveries synthesize")
print("  known ingredients into NEW confirmed predictions; no definitional identity")
print("  collapses them. Tier1 also fires on S665/S666/S669 (correctly).")
print()
print("  TIER 2 (ingredients-known / resemblance): specificity 0% -- would demote")
print("  EVERY genuine discovery (all of physics builds on known ingredients).")
print("  Same instrument flags 'C(rho) ~ Verlinde' (S664) -- so resemblance-based")
print("  reparametrization verdicts are NOT trustworthy on their own.")
print()
print("  => The program's robust results are the TIER-1 ones (S665 exact theorem,")
print("     S666 structural, S669 Delta_r=0 identity, + S661 executed refutation).")
print("     Resemblance-only reparametrization claims (some landscape positioning)")
print("     share a tier that false-positives on BCS and must be held weakly.")
print("     This operationalizes the novelty judgment S662 said was never automated:")
print("     'exact identity to a pre-existing quantity, with no new prediction.'")
