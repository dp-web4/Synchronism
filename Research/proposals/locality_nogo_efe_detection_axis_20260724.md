# Proposal: Engage the EFE-Detection Literature as an Independent Evidence Axis for the Locality No-Go Writeup

**Date**: 2026-07-24
**Origin**: Site visitor track, Pass 4 (modified-gravity researcher persona), 2026-07-24 —
"Has the local-density no-go been checked against Chae's EFE-detection methodology? If the
EFE detection stands, it's independent evidence for non-local (acceleration-keyed)
phenomenology — which *strengthens* the no-go's triage value and would sharpen the
standalone writeup."
**Status**: Proposed (explorer execution candidate; feeds the preprint decision that gates on dp)

## The gap

The locality no-go (ρ_crit ∝ V⁻² required vs V⁺² asserted; ~1.7 dex cross-system offset;
cluster ρ_crit off by 10⁴–10⁶×) cleared its Milgrom prior-art gate on 2026-07-23 and is the
program's #1 transferable null. But the writeup-in-waiting argues *only from refutation*:
local-density coupling fails on data. There is a second, currently unengaged argument from
**positive detection**: Chae et al. (2020, ApJ 904, 51; 2021) claim a ~4σ detection of the
External Field Effect in SPARC — RAR deviations correlated with *external acceleration*, a
strictly non-local variable. The archive already touches this defensively (TEST-08's run
"does not contradict Chae+2020 — different estimator, different claim"; TEST-05's
adjudication quotes Chae's e_N median), but the no-go has never used it offensively:

> If environment enters galaxy dynamics at all, it enters through external *acceleration*
> (detected, disputed) and not through ambient *density* (registered run: r² = 0.0001).
> Both halves of that sentence point the same direction — the organizing variable is
> non-local.

## What to execute (explorer-scale, no new data needed)

1. **Adversarial read of the detection**: Chae et al. 2020/2021 vs the non-EFE readings
   (Freundlich et al. 2022; Paranjape & Sheth 2022 — ΛCDM-expected correlations mimicking
   the signal). The no-go writeup must not overweight a disputed 4σ; the honest form is
   "detected-and-disputed, but *whichever way it resolves*, the candidate variable is
   acceleration, not density."
2. **Sharpen the both-ways structure**: if EFE is real → non-local coupling confirmed →
   no-go's triage sharpened (density-keyed proposals excluded twice over). If EFE is an
   ΛCDM artifact → environment enters through neither variable → the framework's registered
   >20% environment claim was doubly wrong but the no-go itself is untouched. Verify there
   is no branch in which a *local-density* coupling is rescued.
3. **Site propagation**: an EFE paragraph on /for-researchers (the visitor called its
   absence a "missed connection"), with the TEST-05/TEST-08 numbers already on the site.

## Why this is research-direction, not site polish

The preprint strategy (stable fixed point, gate on dp) currently rests the citability case
on prior-art absence. A second independent axis — the live EFE-detection debate keying on
exactly the non-local variable class the no-go predicts must win — converts the writeup from
"we quantified a failure" to "we quantified a failure that the field's sharpest live
discriminator independently corroborates." That is the difference between a null result and
a null result with a forward-looking triage claim.

## Refutation criterion for the proposal itself

If the Chae methodology, on inspection, does not actually discriminate external-acceleration
coupling from ambient-density coupling (e.g., if e_N and ambient ρ are too collinear in
SPARC for the estimators to separate), then the "independent corroboration" framing dies and
the writeup should stay refutation-only. That check is cheap and should be step 0.

---

## STATUS UPDATE 2026-07-24 (explorer): step-0 EXECUTED — SEPARABLE, framing survives scoped

The pre-stated refutation criterion was run same-day: Chae+2020 erratum-corrected Table 2
(arXiv v2 source) × TEST-08 per-galaxy artifacts, N = 141/141. Primary Pearson
r(log e_env, dist-corr density) = +0.432, **r² = 0.187** — below the 0.25 separability
boundary (VIF 1.23; 19–27% shared variance across scales). The framing does NOT die, but it
survives only scoped: the acceleration-vs-density contrast is estimator-dependent (under the
whole-galaxy offset estimator neither variable signals; Chae's 5σ is mean-level,
low-acceleration-weighted, per-galaxy r(e, e_env) ≈ 0), and the detection remains disputed
(Paranjape & Sheth 2022 ΛCDM assembly-bias mimicry; Freundlich 2022; Sargent+ 2025
inconclusive). Both-ways verification complete: **no dispute branch rescues local-density
coupling.** Canonical scoped paragraph for the writeup:
`synchronism-site/explorer/findings/efe-step0-collinearity-e-env-vs-ambient-density.md`
§Implications. Execution record: `explorations/2026-07-24-efe-step0-collinearity-separable.md`.
Erratum trap: golden-galaxy e_env 0.094/0.102 are pre-erratum → corrected 0.040/0.050.
