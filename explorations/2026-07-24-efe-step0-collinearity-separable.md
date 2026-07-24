# EFE evidence-axis step-0 EXECUTED: e_env vs ambient density — SEPARABLE (r² = 0.19), framing survives scoped

**Date**: 2026-07-24 (site explorer track; register→execute gap: 1 day)
**Registered in**: `Research/proposals/locality_nogo_efe_detection_axis_20260724.md` (refutation criterion pre-stated 2026-07-24 by maintainer)
**Script**: `synchronism-site/explorer/scripts/efe_step0_collinearity.py`
**Data**: Chae et al. 2020 (ApJ 904, 51) Table 2 **erratum-corrected**, recovered from the arXiv:2009.11525v2 source tarball (archived `synchronism-site/explorer/data/chae2020_ms_r2.tex`) × TEST-08 registered per-galaxy artifacts (`simulations/test08_per_galaxy_results.json`). Join: **N = 141/141**.
**Full finding**: `synchronism-site/explorer/findings/efe-step0-collinearity-e-env-vs-ambient-density.md`

## Pre-declaration (written to file before computing)

Primary: Pearson r between log₁₀(e_env) and TEST-08's registered dist-corrected log(1+N_cyl).
Rule: r² ≥ 0.5 → collinear, corroboration framing DIES; r² < 0.25 → separable, survives; else gray.

## Result

**r = +0.432, r² = 0.187 → SEPARABLE.** VIF 1.23; identical in the low-acc subset (r² = 0.19, N = 106); 5 Mpc sphere secondary reaches r² = 0.27 → honest headline: **19–27% shared variance — separable, not orthogonal**. Spearman 0.54 (had Spearman been primary, verdict would be GRAY — sensitivity disclosed).

Instrument validated both directions: TEST-08's r² = 0.0001 null reproduces exactly on the join; Chae's headline reproduces exactly from the parse (N = 113, median e = 0.052, median e_env = 0.033, e>0 count 77≈78).

## What the writeup may and may not claim

**May**: the two variables are statistically separable in SPARC; the registered density arm is dead by execution (r² = 0.0001 vs >20% claimed); the live detection debate keys on external acceleration; **no branch of the dispute rescues local-density coupling** (real EFE → acceleration; ΛCDM mimicry [Paranjape & Sheth 2022, MNRAS 517, 130] → halo concentration via assembly bias; inconclusive [Sargent et al. 2025, arXiv:2511.03839] → density null stands).

**May NOT**: the unscoped "environment enters through acceleration, not density" as an equal-estimator variable contrast. Diagnostic 2×2: under the whole-galaxy mean-offset estimator **neither** variable signals (e_env r = −0.11 ns, density r = +0.01); Chae's 5σ is a mean-level result under low-acceleration weighting, and per-galaxy r(e, e_env) ≈ 0 (medians-only; individual e uncertainties span the e_env range). The contrast is estimator-dependent and ships scoped or not at all.

**Bonus sign structure**: r(e_env, RAR offset) = −0.11 (ns) carries the MOND-EFE sign — opposite Synchronism's registered void-high prediction. Both environment arms point anti-Synchronism. P&S 2022's sign discriminator (MOND vs ΛCDM mocks predict opposite signs) is testable on this pipeline with a low-acceleration-weighted offset — seeded as a follow-up (literature contribution, not a Synchronism test).

## Erratum trap (site-facing)

The widely quoted golden-galaxy e_env ≈ 0.094 (NGC5055) / 0.102 (NGC5033) are **pre-erratum**; corrected values are 0.040 / 0.050. Any site or archive text citing per-galaxy e_env must use the corrected Table 2. (The comparison literature has its own numbers-outliving-computation failure mode.)

*Run executed as registered; result recorded either way, per the standing kill-criterion rule.*
