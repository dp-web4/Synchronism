# Back-annotation: RAR Transition-Shape Test Closes Galaxy Research Program

**Source**: synchronism-site explorer session 2026-05-21 + maintainer WAKE 2026-05-22
**Status**: Research closure — supersedes `rar_transition_shape_discriminator.md` (which documents the test execution)
**Back-annotation from site**: visitor logs 2026-05-22 (Pass 3 + Pass 4 both independently flagged this result as invisible on all researcher-facing pages)

## Research Conclusion

The only non-degenerate galaxy-scale discriminating test between the Synchronism γ=2 compander
and McGaugh's MOND interpolating function has now been **executed and closed** on real SPARC data
(Lelli-McGaugh-Schombert 2016, 2807 points after quality cuts).

**Result** (from `synchronism-site/explorer/findings/rar-transition-shape-real-sparc-result.md`):

| Model | a₀ (m/s²) | RMS (dex) | ΔBIC vs McGaugh |
|-------|-----------|-----------|-----------------|
| McGaugh ν (standard MOND) | 1.13×10⁻¹⁰ | 0.1437 | reference |
| Compander μ, γ=2 pinned | 2.97×10⁻¹⁰ | 0.1485 | **+184** |
| Compander μ, γ free | 5.3×10⁻¹¹ (γ=0.49) | 0.1437 | +7.1 |

**Kill criterion**: ΔBIC > 10 favoring McGaugh refutes the γ=2 compander. Actual ΔBIC = +184.
Conservative intra-galaxy correlation correction (effective N≈500–1000): ΔBIC ≈ 33. Still decisive.

## What This Closes

1. **The compander (γ=2) is refuted as the galaxy mechanism.** ΔBIC=+184 (conservative ≈33)
   leaves no room for the γ=2 fixed-parameter compander as a description of galaxy dynamics.

2. **Free-γ fit converges to γ=0.49 ≈ MOND.** The RMS at free γ equals McGaugh's to four digits.
   ΔBIC=+7.1 is entirely the BIC penalty for the extra free parameter — the fit improvement is zero.
   The compander at its best-fit γ is MOND.

3. **No γ is simultaneously distinct from MOND and consistent with SPARC.** This is the fork
   the 2026-05-20 proposal identified as the decisive question: pinned γ=2 → refuted; free γ → MOND.
   The answer is now empirical.

4. **Net discriminating galaxy tests vs MOND+ΛCDM: 0, by execution.** Previously the count was
   "1 untested (RAR shape)." As of 2026-05-21 it is "0 remaining discriminators."

## What This Does NOT Close

- The environment-dependence tests (TEST-01, TEST-05) remain pending — but they were already
  MOND+EFE degenerate in direction (see `tier1_mond_efe_discriminator_gap.md`). Their discriminating
  power was contingent on the compander having daylight from MOND at some scale; the shape test
  shows there is none at the galaxy scale.
- The DESI fσ₈ result (TEST-04a) is a separate, cosmological closure (sign-reversed mechanism-class
  failure; note: post-hoc retrodiction, not a prospective prediction).
- Wide-binary test (TEST-02) is blocked pending resolution of the anomaly dispute.

## Research Direction Implications

### What the project now owns (on the physics side)

After this closure, the honest physics summary is:
- The compander (C(ρ)) is a phenomenological μ-law/Hill transfer function with no governing
  free energy, no self-consistency equation, and no universality class
- At its best-fit γ it reproduces MOND's interpolating function
- Its only novel knob (environment-dependent scatter, TEST-03) failed its kill criterion (R²=0.14)
- Zero non-degenerate discriminating tests remain

### What the project genuinely owns (on the methodology side)

The A2ACW analysis (see `a2acw_v2_three_axis_protocol.md`) produced a genuine result:
adversarial AI self-play over a shared corpus is a **reparametrization detector, not a discovery
engine**. The three-axis protocol (vocabulary-translation, symbol-audit, null-baseline computation)
catches 6/6 reparametrizations. This is a methodological contribution independent of the physics.

### Proposed research reframe

The galaxy program is closed. The physics program is closed (0 novel confirmed predictions,
0 remaining discriminating tests, 1 refuted mechanism class). The site should shift from
"active physics framework with live discriminating tests" to "completed case study of
AI-driven reparametrization detection."

The honest next step is not to search for another discriminating test — it is to document the
methodology result (A2ACW three-axis protocol) as a standalone contribution and to present the
physics framework as the worked example that produced it.

## Action for site maintainer

These results need to be visible on the four pages where a physicist would look:
1. `/tier-1-existing` — "Recommended Start" section: add RAR shape closure note
2. `/honest-assessment` — "What Was Tested" or "What Failed": add RAR shape result card
3. `/galaxy-rotation` — add RAR shape test section after "Honest Caveat"
4. `/test-catalog` — "What's Already Been Analyzed": add RAR shape result
