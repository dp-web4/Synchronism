# The distinct external-field predictions live on the excluded branch; a disjunctive registration let a morphology effect stand as support

*Site maintainer, 2026-10-21. Raised by simulated-visitor graduate and researcher personas (visitor log 2026-10-21), checked against the archive. No bucket moves; count 6; Bucket 0 = 0.*

## 1. No MOND-distinct external-field prediction exists on a live branch

The site's /mond-unification called linear refraction of an external field "the one genuinely distinct external-field
prediction of this family", untested. It also called "two satellites at the same g_ext but different host gas content" "the
genuine discriminator", and said TEST-05 had run it. A researcher persona separated the two and pointed out that TEST-05 had
not run the second. TEST-05 used Mpc-scale Cosmicflows-4 density and held neither g_ext nor the satellite's ambient gas
fixed, so it refuted a registered amplitude (> 20% of scatter), not the discriminator.

Both predictions are the **density-keyed** branch's:
- **Refraction (L2):** ∇·[C(ρ)∇Φ] = 4πGρ bends an external field wherever ∇C ≠ 0. Estimated at 0.3–0.37 g_ext beyond the knee in
  a halo globular cluster (2026-09-16 block).
- **Ambient density at fixed g_ext:** capped by the 2026-08-21 tidal identity (L ≤ k(3−n)/9 ≈ 0.1–0.2) and systematics-limited
  (signal/systematics ≈ 4×10⁻³; the baseline offset correlates with local density at r = −0.36).

That branch is excluded twice: SPARC ΔBIC +2843 against the acceleration-keyed law, which survives galaxy-level scoring
(10-18: 95% [2178, 3590]), and LLR (68/74). The acceleration-keyed branch the fits use carries QUMOND's own EFE, so it has
nothing distinct here. So:

> **No galaxy-scale prediction that differs from MOND exists on a branch that survives.** The distinct ones are predictions of a
> branch already excluded on other grounds. Testing them can add a root to that exclusion, or contradict it. It cannot move
> Bucket 0.

This is a sharper form of the 10-13 scoped sentence ("no prediction that differs from MOND or ΛCDM survives among the
registered tests and the readings executed"). That sentence leaves room for unregistered distinct predictions; this closes
that room for the galaxy sector's external-field content.

**A correction to the persona's suggested cheap kill.** The persona proposed that the refraction amplitude might be negligible
because C_ρ ≈ 10⁻³ in discs. That holds for the floored form, where C ≈ Ω_m below the knee and ∇C ≈ 0. It does not hold for
the unfloored form. Below the knee C ≈ γρ/ρ_c, so ∇ln C = ∇ln ρ, which is O(1/h) however small C is (the site already said
refraction is "order R/h in any thin disc"). Smallness of C is not smallness of the refraction. The amplitude for a satellite
is uncomputed; an explorer topic is seeded.

**Recommendation (gates on the project lead):** adopt the sharper sentence as the galaxy-sector form of the scoped negative.
List the two predictions with an "excluded branch" flag rather than as open discriminators. The site now does this
(/test-catalog#unregistered-distinct).

## 2. TEST-05's kill text is a disjunction, and the wrong leg fired

TEST-05's kill reads "RAR scatter independent of Hubble type / environment". The framework predicts an environment
(density) dependence. It predicts nothing about Hubble type. Sessions 376–378 found RAR scatter depends on Hubble type on SPARC
(p = 5×10⁻⁶, R² = 0.14, surviving five confound controls), and S378 scored that as "NP2 STRONGLY SUPPORTED". S381 then showed
the dependence is structural: within each type an isolation proxy predicts nothing (r = +0.003, p = 0.97), and rotation-curve
roughness absorbs the type signal (type p = 0.61 given roughness). The 2026-07-14 environment run gives r² = 0.0001 (N = 141).

The site carried "environment effect, p = 5×10⁻⁶, real and statistically significant" on /galaxy-rotation (a headline tile)
and /rar-scatter (a green tile) until today, beside the r² = 0.0001 kill. Two personas caught the contradiction independently.
The 07-09 corrections fixed the sample attribution (SPARC, not ALFALFA) but kept the "environment" label.

The mechanism is general. **A kill criterion phrased as "independent of X / Y", where only Y is the framework's variable, can be
met by an X effect the framework never predicted.** The rejection is then read as support. This is the "read the registration
text" lesson from the other side: there the registration named a different parameter than the program ran; here it names one
extra variable. **Recommendation:** grep the TEST-01…26 registry for criteria containing "/" or "or" that name a non-framework
variable, and scope each to the framework's variable. Rewriting registration text gates on the project lead; this proposal
does not edit it.

## 3. A2ACW: arm dependence recorded; the demotion auditor traced (explorer 2026-10-20)

- Re-running LV's three targets Claude-vs-Claude: run-to-run agreement 0.86, agreement with LV 0.41–0.52, and all 17
  disagreements harsher. LV's T2 "HELD" does not reproduce. Vendor, model generation, round count and the human relay are
  confounded.
- S662:50's "human audit" is wrong. All six demotions were AI-track edits. Back-annotated in place today.
- **My own error:** the 10-20 proposal's "13 verdicts (5/4/4)" was recalled, not counted. The logs give 27 (10/8/9).
  Corrected in that file. It was the RECALL failure mode, on a record-facing number, in the lane that proposed the RECALL code.

## Files
- Site: `src/app/galaxy-rotation`, `rar-scatter`, `mond-unification`, `test-catalog` (#unregistered-distinct), `key-claims`,
  `cosmology-predictions`; `maintainer/tools/site_lint.py` (224 rules, strict flag added).
- Archive: Session 381 (the structural reading), Session 378 (the "STRONGLY SUPPORTED" score).
