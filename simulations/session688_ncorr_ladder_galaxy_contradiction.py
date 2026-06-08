#!/usr/bin/env python3
"""
Session 688: verify the load-bearing 'both-directions contradiction' at the
galaxy rung in the ncorr_ladder_never_anchored proposal (2026-06-08).

Claim:
  - Framework asserts galaxy stars are independent => N_corr = 1 => gamma = 2
    (under original formula gamma = 2/sqrt(N_corr)).
  - SPARC RAR rejects gamma = 2 at Delta BIC = +184 (S661); the data-preferred
    free fit gives gamma ~ 0.49 (matches McGaugh's empirical interpolating fn).
  - Inverting the formula at gamma = 0.49 gives N_corr = (2/gamma)^2 ~ 17,
    contradicting the N_corr = 1 premise.
  - Therefore the galaxy rung is internally inconsistent in the framework's own
    law: asserted side gives rejected gamma; fitted side voids the N_corr
    assignment that licenses applying C(rho) at galaxies in the first place.

Also check: does the proposed S686 sign flip gamma = 2*sqrt(N_corr) repair this
contradiction? Inverting that at gamma = 0.49 gives N_corr = (gamma/2)^2 ~ 0.06,
which is <1 star -- nonphysical. So the contradiction is sign-independent: it
survives both the original formula and the S686-proposed flip.
"""
import numpy as np

print("=" * 74)
print("(1) Both-directions contradiction at the galaxy rung")
print("=" * 74)
print("   ASSERTED SIDE:")
print("     'Galaxy stars are independent' -> N_corr = 1")
N_corr_asserted = 1
gamma_orig_from_assertion = 2.0 / np.sqrt(N_corr_asserted)
print(f"     gamma = 2/sqrt(N_corr) -> gamma = 2/sqrt({N_corr_asserted}) = "
      f"{gamma_orig_from_assertion:.3f}")
print(f"     SPARC verdict (S661): gamma = 2 REJECTED at Delta BIC = +184")
print(f"     SPARC free-fit: gamma_best ~ 0.49 (= McGaugh interpolating function)")
print()
print("   FITTED SIDE:")
gamma_fit = 0.49
print(f"     gamma_best = {gamma_fit}")
N_corr_from_fit_orig = (2.0 / gamma_fit) ** 2
print(f"     gamma = 2/sqrt(N_corr) inverted: N_corr = (2/gamma)^2 = "
      f"(2/{gamma_fit})^2 = {N_corr_from_fit_orig:.2f}")
print(f"     -> ~17 stars correlated at the galaxy rung")
print(f"     But the framework asserts N_corr = 1 (stars independent) at this rung.")
print(f"     The fitted side voids the N_corr = 1 premise that licenses applying")
print(f"     C(rho) at galaxies in the first place.")
print()
print(f"   RESULT: asserted N_corr = 1 gives gamma = 2 (rejected by data);")
print(f"           fitted gamma = 0.49 forces N_corr = 17 (contradicts the")
print(f"           independence premise). Self-contained internal inconsistency.")

print()
print("=" * 74)
print("(2) Does the S686 sign flip  gamma = 2*sqrt(N_corr)  repair this?")
print("=" * 74)
print("   Under flip:")
gamma_flip_from_assertion = 2.0 * np.sqrt(N_corr_asserted)
print(f"     N_corr = 1 -> gamma_flip = 2*sqrt(1) = {gamma_flip_from_assertion:.3f}")
print(f"     (galaxy fixed point: same as original at N_corr = 1)")
N_corr_from_fit_flip = (gamma_fit / 2.0) ** 2
print(f"     gamma_fit = 0.49 inverted under flip: N_corr = (gamma/2)^2 = "
      f"({gamma_fit}/2)^2 = {N_corr_from_fit_flip:.4f}")
print(f"     -> N_corr ~ {N_corr_from_fit_flip:.3f} (less than ONE star correlated)")
print(f"     Nonphysical: N_corr must be a positive integer >= 1.")
print()
print("   RESULT: the sign flip does NOT repair the galaxy-rung contradiction.")
print("   Under the original formula the fitted gamma gives N_corr ~ 17 (too")
print("   many); under the flip it gives N_corr ~ 0.06 (less than one). Either")
print("   direction violates the N_corr = 1 premise. The contradiction is")
print("   sign-INDEPENDENT.")
print()
print("   This is a new constraint on the S686 fork: the 'A + flip' choice that")
print("   qualitatively repaired the BCS/BEC ladder still leaves the galaxy")
print("   rung internally inconsistent. The flip is not a universal fix.")

print()
print("=" * 74)
print("(3) What N_corr could the formula even support?")
print("=" * 74)
print("   gamma must satisfy: gamma_data = 2/sqrt(N_corr_asserted)")
print("   For galaxies to be a CONSISTENT rung (asserted N_corr = fitted N_corr):")
print(f"     N_corr_data = (2/gamma_data)^2 = (2/0.49)^2 = {N_corr_from_fit_orig:.2f}")
print("     The framework would have to assert 'galaxy stars correlate in groups")
print("     of ~17' instead of 'galaxy stars are independent'.")
print()
print("   This re-assertion has no independent physical motivation -- there is no")
print("   17-star natural correlation length in galactic dynamics. It would be")
print("   chosen specifically to satisfy the fit, which is the back-fit pattern")
print("   the proposal flags ('N_corr asserted, not measured, on every rung').")
print()
print("   The galaxy rung therefore offers no escape: independence assertion +")
print("   formula gives rejected gamma; data fit forces N_corr to a value with")
print("   no independent motivation; either choice nullifies the universal-law")
print("   claim for this rung.")

print()
print("=" * 74)
print("(4) Integration: 4 sub-sectors of structural fork pattern + this")
print("=" * 74)
print("   S661 galaxy RAR shape: gamma=2 refuted; gamma-free = MOND")
print("   S678/S683 cluster: wrong-variable; C(a) restoration = MOND")
print("   S684 EFE: bounded refuted; B_max -> inf = MOND")
print("   S686 N_corr ladder under flip: repairs sign inversions qualitatively")
print()
print("   S688 (this): the galaxy rung -- the ONLY rung where the framework")
print("   has any quantitative success at all -- is internally inconsistent")
print("   in its OWN ladder formula, sign-independent. The 'one equation 80")
print("   OOM' framing requires at least one rung where N_corr is counted")
print("   independently and the predicted gamma survives a like-for-like test;")
print("   the proposal's audit (4 confronted rungs, 0 surviving) leaves zero.")
print()
print("   Substance is consistent with prior sessions (S661 refutation at this")
print("   rung was already established). The proposal's NEW contribution is")
print("   framing the refutation as a both-directions internal contradiction")
print("   in the ladder's OWN law, not just a data-vs-prediction mismatch.")
