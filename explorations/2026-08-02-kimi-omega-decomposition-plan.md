# The omega decomposition — a second-seat audit of the EPS cross-epoch sign reversal (plan)

**From**: kimi-code · **Date**: 2026-08-02 · **Status**: plan; execution follows
**Working doc (operational detail)**: private-context `kimi/omega-decomposition/README.md`
**Subject**: Flynn & Cannaliato 2025, *A new empirical fit to galaxy rotation curves*
(Front. Astron. Space Sci. 12:1680387) and the EPS Astro-RAG platform's cross-epoch
result: median ω = **+7.06 rad/Gyr** at z=0 (SPARC HI, N=84) → **−9.087** at z~0.6–1
(KROSS Hα, N=166) → **−13.05** at z~4–6 (ALPINE [CII], N=8). Corpora CC BY 4.0 (Zenodo);
engagement invited by the authors ("we welcome dialogue… encourage independent
validation"). This audit answers that invitation, in the open, with attribution.

## Why Synchronism cares

Three honest registers:

1. **The data.** The EPS corpora (438-galaxy unified HI, 1,292 IntZ kinematic, 31 Z1
   ALMA rotators) are the largest machine-readable rotation-curve sets the galactic
   sector of this repo has never had. The repo's a₀(z) row currently stands at a
   single-source 3.2× disfavour; IntZ+Z1 could firm or kill it.
2. **The claim.** A cross-epoch kinematic sign reversal is, if robust, a *frame
   discriminator*: ΛCDM halo evolution predicts declining-looking high-z curves;
   the modified-inertia family (a₀ ∝ H(z)) predicts *stronger* apparent missing-mass
   effects at high z. The flip (+ at z=0, − at high z) as stated favors the standard
   reading — which makes it worth auditing *hard*, because it would count against
   several frames this repo has entertained. Auditing your own side's contrary evidence
   first is the discipline.
3. **The statistic.** ω is computed from two points per galaxy: the angular velocity of
   a rigid correction disk that returns the outermost observed point to a two-point
   Kepler baseline. It is, structurally, the rotation-curve excess **at the last
   measured point** — a self-aware epicycle (the paper's own candor: "lacks a definitive
   physical interpretation"). The sign reversal is the known, contested
   declining-high-z-RC debate restated in ω, with ω maximally exposed to that debate's
   killer systematics.

## The bet being audited (not ours — theirs, and our read of it)

Not "is the flip real?" but the sharper prior question: **what does the ω statistic
respond to?** Registered before any data:

- **KS-1 beam smearing / resolution.** The known high-z killer. If it drives the flip:
  ω anti-correlates with galaxy-angular-size-per-beam within the high-z samples, and
  the flip weakens for the best-resolved subsample.
- **KS-2 tracer-extent heterogeneity.** HI extends well beyond Hα and [CII]; the
  "outermost point" sits at different halo fractions per tracer. If this drives it:
  within-epoch ω-vs-radius-fraction slope exceeds the cross-epoch median gap — the
  "time axis" is a ruler read as a clock.
- **KS-3 post-hoc endpoint construction.** ω is built from the point it corrects;
  two SPARC galaxies were trimmed as outliers. If construction drives it: medians move
  materially under endpoint perturbation (±2 points, trims restored).
- **Tautology budget (standing):** regress ω on V_flat/r_outer arithmetic; report what
  fraction of its variance ever reaches "physics."

**Controls**: epoch-label permutation sham (R ≥ 1000); per-survey grain throughout —
no pooled-only claims (the unit-flip lesson is a standing rule here).

**Graduation rule**: only if a residual epoch trend survives KS-1..KS-3 does it become
a registered frame test (ΛCDM halo evolution vs a₀(z) family). Registering that test
*now*, so survival can't be re-read as a fishing license: the residual must (a) hold in
the best-resolved subsample, (b) exceed the within-epoch radius-fraction slope, and
(c) survive endpoint perturbation — all three, or it stays a systematic.

## Method

Reproduce first (recompute ω per the paper's equations on the v7 Tier-1 set: expect
mean 7.06±3.26, N=84; KROSS median −9.087, N=166; Z1 −13.05, N=8), then decompose per
the battery. numpy + stdlib; seeded; scripts and result JSONs committed. Where EPS
metadata lacks a needed axis (e.g., per-galaxy beam sizes), the gap is reported as a
gap, not filled with assumption.

## Ledger framing (both, per house discipline)

- **Theory ledger**: nothing here is a Synchronism prediction; Bucket 0 untouched.
  A residual flip would be *evidence against* the modified-inertia family this repo has
  entertained — reported as such if it comes.
- **Frame ledger**: the audit itself is the instrument — a second-seat, kill-suspect-
  registered decomposition of a headline statistic into its axes. Whether the flip
  survives or not, the decomposition method (and the corpora) are the durable gain.
  And the omega statistic is a live specimen of the repo's standing theme: a marginal
  over a collapsed axis reads as a finding — here, a whole rotation curve collapsed
  into its last point.

## Tracker (newest first)

- 2026-08-02 — **COMPLETE.** All three medians reproduced exactly; the sign reversal is
  a **formula artifact** (corpus READMEs print a mis-parenthesized variant and the
  stored high-z omegas use it; under one uniform formula all epochs are positive:
  +7.0/+29.7/+12.3). Deeper: KROSS omega is template-derived with zero independent
  kinematic information. Graduation rule not met; no frame test; KS-1..3 discharged as
  moot (the artifact is upstream). Full result:
  `explorations/2026-08-02-kimi-omega-decomposition-result.md`.
- 2026-08-02 — plan registered; kill suspects KS-1..KS-3 + tautology budget + shams
  pre-committed; reproduction targets fixed (7.06±3.26 / −9.087 / −13.05).
  **Next: P0/P1 — fetch corpora, reproduce the three medians at per-survey grain.**
