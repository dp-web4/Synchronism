# Session 639: TEST-03 Disambiguation — Two Different Numbers, Two Different Tests

**Date**: 2026-04-30
**Type**: Site-Archive-Audit (9th instance, post-arc-closure)
**Trigger**: 2026-04-30 proposal `test03_kill_criterion_self_trigger.md`
**Grade**: A- (specific, resolves a self-trigger ambiguity, exposes a metric conflation)

---

## Setup

The Framework Stress Test arc was declared CLOSED at 22 sessions earlier today (2026-04-30 03:36 Publisher report). At 06:19 same day, a new visitor proposal arrived flagging that TEST-03's reported R²=0.14 may have triggered its own <20% kill criterion. The proposal asked: which of three interpretations applies (literal trigger, mis-specified criterion, no real derivation)?

S639 traces both numbers in the archive and finds the answer: **the site conflates two distinct metrics**. The kill criterion is being compared against the wrong quantity.

## Two Numbers, Two Different Tests

The site's TEST-03 mixes a derived prediction with a separate environmental ansatz:

| Metric | Source | What it measures | Status |
|--------|--------|------------------|--------|
| **51% improvement** | Session 593 (Grade A) | TFR residual reducing BTFR (Baryonic Tully-Fisher) scatter, σ: 0.402 → 0.195 dex on 14,437 galaxies | Derived; would *not* trigger <20% |
| **R² = 0.14** | Site `/galaxy-rotation` | Environmental density explaining RAR (Radial Acceleration Relation) scatter on 14,585 galaxies | Separate ansatz; *would* trigger <20% if measured against this |

### Verbatim from Session 593

> "The 58.3% drops to 51.4% when removing L_i from Mbar, confirming some algebraic circularity. But 51.4% is still massive — the i-band TFR residual genuinely predicts BTFR scatter on 14,437 galaxies, even when Mbar is constructed entirely from SPS masses."

### Verbatim from `/galaxy-rotation`

> "The environment-dependent scatter is real and statistically significant (p = 5×10⁻⁶ with N = 14,585), but it explains only 14% of the total RAR scatter (R² = 0.14)."

The 51% is BTFR; the 14% is RAR. The 51% predictor is TFR residual (M/L proxy); the 14% predictor is environmental density. They are independent, additive contributions — Session 594 line 58 shows TFR + environment gives 55.1% combined improvement (the 4% environmental adds onto the 51% TFR).

## Resolving the Three Interpretations

The proposal listed three possibilities:

- **A. Kill criterion is triggered** (TEST-03 is Failed)
  → Only if the kill criterion compares against the 14% RAR environmental result. But the site's 51% prediction is for BTFR, not RAR. **Not triggered against the actual derived prediction.**

- **B. Kill criterion was misspecified** (the criterion's denominator is wrong)
  → Closest to correct. The kill criterion is anchored to the wrong number. The 51% prediction is for BTFR scatter (TFR residual is the predictor); the 14% report is for RAR environmental ansatz. **The site needs to separate them and apply each kill criterion to its proper metric.**

- **C. Prediction was never derived** (51% is post-hoc)
  → False. The 51% is derived in Session 593 with explicit SPS-mass baseline (σ=0.402 dex) and post-TFR residual (σ=0.195 dex). Sessions 594 (decomposition) and 596 (synthesis) confirm. **The prediction has a clean derivation chain.**

**Verdict: Interpretation B with refinement.** The 51% prediction (BTFR scatter via TFR residual) survives any reasonable <20% threshold by 2.5×. The 14% result (RAR environmental ansatz) falls below threshold but is a different, weaker claim that should not be conflated with TEST-03's headline number.

## What the Site Should Do

This is operator/site-side action, but the recommended structure:

1. **Split TEST-03 into two clearly-labeled tests:**
   - TEST-03A: TFR residual reduces BTFR scatter — prediction 51%, kill <20%, **passing** (derived in Session 593, sample 14,437)
   - TEST-03B: Environmental density explains RAR scatter — prediction "significant signal," kill <20%, **below threshold** (R²=0.14, sample 14,585)
2. **Reclassify TEST-03B** as Speculative-Below-Threshold or move to honest-assessment per the proposal's recommendation.
3. **Compute ΔBIC vs MOND-only** for the environmental term (the site already flags this as missing). This is the proper test of whether the 14% adds anything over MOND's RAR fit.

## Why This Matters Beyond Bookkeeping

The conflation pattern is the same one S631–S638 exposed at the foundation level: **a site-level metric and an archive-level derivation share a label but measure different things**. In S631 the site claimed n≈2.2 was the BTFR exponent (full sample) when archive showed it was only the "baryonic component" subset. In S639 the site reports R²=0.14 next to a "<20%" criterion when the 14% measures RAR environmental scatter and the threshold is anchored to the BTFR TFR-residual claim.

These conflations are how a sub-threshold result coexists with an above-threshold derivation under the same TEST-ID. The audit channel is now identifying not just isolated errors but the *mechanism* by which the site/archive divergence accumulates: shared labels, distinct measurements, no enforced alignment.

## Audit-Channel Taxonomy Now 9 Modes (Post-Closure Extension)

| # | Type | Session |
|---|------|---------|
| 1 | Quantitative refutation + mislabeling | S631 |
| 2 | Dimensional inconsistency | S632 |
| 3 | Structural overclaim | S633 |
| 4 | Count discrepancy | S634 |
| 5 | Domain-level badge overclaim | S635 |
| 6 | Category error | S636 |
| 7 | Derivation succeeds but predicts undetectable signal | S637 |
| 8 | External-track derivation independently verified | S638 |
| 9 | **Metric conflation under a shared TEST-ID (NEW, post-arc)** | **S639** |

The arc was formally closed at 22 sessions earlier today. S639 arrives within hours of closure with a 9th distinct mode. This suggests the visitor channel continues to be productive even after the framework's predictive content was declared characterized — the audit pattern continues finding new error types as long as there are public claims to check against.

## Note on Arc Status

The Publisher's 2026-04-30 closure was based on "predictive content fully characterized." That claim still stands — S639 doesn't re-open foundational physics questions. But it does show the *correction process* is incomplete: the site has structural errors that haven't propagated to the public surface. Whether to reopen or treat S639 as a coda is the operator's call. From here it looks like an extension of the same methodology, applicable as long as proposals continue arriving.

## Files

- `Research/Session639_TEST03_Kill_Disambiguation.md` (this document)
- No simulation needed — this is a definitional clarification

## So What?

TEST-03 is not self-triggered. The 51% prediction (TFR residual on BTFR) is derived and passes its kill threshold by 2.5×. The 14% report (environmental density on RAR) is below threshold but is a different test that has been mislabeled under the same ID. The site should split the two and apply each kill criterion to its proper metric.

The pattern S639 exposes — *metric conflation under a shared TEST-ID* — is the mechanism by which site/archive divergence accumulates without anyone noticing: same label, different measurement, asymmetric correction propagation. Catching it took looking at both numbers and asking what each actually measures. The visitor channel did that; this session traced the answer in the archive.
