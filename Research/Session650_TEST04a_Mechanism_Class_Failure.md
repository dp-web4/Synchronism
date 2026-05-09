# Session 650: TEST-04a — Mechanism-Class Failure (Sign Reversed)

**Date**: 2026-05-09
**Type**: Site-Archive-Audit (19th instance — sharpening of S645/S648)
**Trigger**: 2026-05-09 proposal `test04a_mechanism_class_sign_failure.md`
**Grade**: B+ (taxonomy refinement; updates TEST-04a status)

---

## Setup

S645 documented DR1's disfavoring of Session 107's fσ₈ prediction. S648 corrected the framing from "first hard external falsification" to "post-hoc consistency failure" given the timestamps. The 2026-05-09 proposal adds a third sharpening: this is not a magnitude miss but a **mechanism-class failure** — the predicted sign is wrong.

S650 adopts the taxonomy and updates the status.

## The Taxonomy

The proposal's contribution is a clean three-way classification:

| Failure type | Example | Repairable by retuning? |
|--------------|---------|------------------------|
| Magnitude miss | Melting points 53% error | Yes — refine density-to-property mapping |
| Universality miss | Critical exponents 2× off | Partial — change functional form / universality class |
| **Mechanism-class failure** | **TEST-04a: sign-reversed fσ₈** | **No — mechanism class predicts wrong sign** |

Session 107's mechanism (G_local/G_global = C_cosmic/C_galactic < 1) is a *suppressor* that predicts fσ₈ below ΛCDM at low z, converging at high z. DESI DR1 observes the opposite: fσ₈ *above* ΛCDM at LRG1/LRG2 (low z), converging at high z (ELG2). The redshift pattern is *inverted*.

A suppressor that observes enhancement cannot be retuned within the suppressor class. Magnitude knobs (re-tune σ₈(z=0), re-tune coupling normalization) cannot flip the sign. The mechanism class itself is wrong.

## Why This Strengthens the Verdict

S648 corrected S645's framing to "post-hoc consistency failure." S650 sharpens further: even as a post-hoc consistency check, the failure is more severe than S648 implied. It is not "framework parameters can't reproduce DR1 magnitudes" (which a retune might fix). It is "the framework's mechanism class predicts the wrong sign of the redshift dependence" (which retuning cannot fix).

Updated status:
- **TEST-04a**: REFUTED (post-hoc, mechanism-class, sign-reversed). The framework's suppressor mechanism predicts the opposite of what DR1 observed.
- Combined with cosmology scorecard (S635: 0 novel-unfalsified claims) and meta-criterion logic (S646), the cosmological domain has now formally exhausted its novel-prediction content.

## The Branch 1 Diagnostic (Optional Research)

The proposal observes that **if** C_galactic/C_cosmic > 1 (instead of < 1) — i.e., dense halos are *more* coherent than the cosmic average — the prediction sign would flip and become consistent with DR1's direction. Whether this is consistent with the rest of Session 107's assumptions is an operator-discretion research question.

This is interesting because it would mean Session 107 had the structural form right (local coupling coherence drives growth) but the sign wrong. From the audit channel's side, no commitment is made; the recommendation is operator-level.

## Cosmology Sector Now Formally Exhausted

Combining prior sessions:
- **S635**: cosmology scorecard found 0 novel-unfalsified claims across 15 cosmology entries
- **S645/S648/S650**: TEST-04a is mechanism-class failure (sign-reversed), refuted post-hoc by DR1
- **S646**: meta-falsification criterion recommended retraction when both cosmological and galactic novel-prediction domains fail kill criteria

Per S646's recommended M3 (scope reduction), the cosmological domain meets the retraction condition. The galactic-domain TEST-03B is below threshold; TEST-03A passes but is MOND-shared (S637). The combined picture supports formal scope-narrowing on the public site, retiring novel-prediction claims in cosmology and the environment-dependent galactic claim.

## Audit Taxonomy

| # | Type | Session |
|---|------|---------|
| 17, 18 | QM kill underspecified, ρ_crit asymmetry | S649 |
| 19 | **Mechanism-class failure (sign reversal not retunable)** | **S650** |

S650's contribution is taxonomic: separating "magnitude miss" from "mechanism-class failure" makes the meta-criterion logic (S646) actionable. A meta-criterion that doesn't distinguish these would treat retunable failures and irreparable failures as equivalent.

## Recommended Site Action

- **TEST-04a label**: upgrade from "DISFAVORED 2.4σ" to **"REFUTED — mechanism-class failure (sign-reversed)"** with date and Adame+2024 citation.
- **`/honest-assessment`**: add the three-tier failure taxonomy (magnitude / universality / mechanism-class).
- **`/tier-1-existing` TEST-04a**: same upgrade.
- **`/key-claims` cosmology section**: scope-narrow per S646 recommendation. No surviving cosmological novel prediction.

## Files

- `Research/Session650_TEST04a_Mechanism_Class_Failure.md` (this document)

## So What?

S645 → S648 → S650 is a sequence of sharpenings on the same finding: TEST-04a refuted by DESI DR1. Each session corrected the framing. S650's contribution is the **mechanism-class** classification — the failure is irreparable within the suppressor class because the sign of the redshift dependence is reversed. Magnitude retuning can't save it.

Combined with S635 (cosmology scorecard) and S646 (meta-criterion), the cosmological sector has formally exhausted its novel-prediction content. The honest position is scope-narrowing per S646's M3 recommendation. Operator-level decision on adoption.

The Branch 1 sign-error diagnostic (could C_galactic/C_cosmic > 1 flip the sign?) is interesting research but doesn't affect the current verdict. Productive failure > safe summaries.

Cumulative count: 19 internal site-archive audits + 1 mechanism-class refuted prediction (TEST-04a per S645/S648/S650). Audit queue caught up.
