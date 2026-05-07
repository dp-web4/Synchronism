# Session 648: Session 107 Pre-Registration — S645's Framing Was Too Strong

**Date**: 2026-05-08
**Type**: Self-correction of S645 (16th audit instance, post-arc-closure)
**Trigger**: 2026-05-07 proposal `session107_preregistration_gap.md`
**Grade**: B+ (sharpens own prior session; honest revision)

---

## Setup

S645 framed the DESI DR1 disfavoring of Session 107's fσ₈ prediction as "the framework's first hard external falsification" — implying a prospective prediction tested by subsequent data. The 2026-05-07 proposal flags a timestamp gap: Session 107 was committed 2025-12-10; DESI DR1 was publicly available from April 2024 with full-shape RSD results in arXiv:2411.12021 (November 2024). Session 107 is therefore *post-DR1*, not prospective.

S648 verifies and corrects S645's framing.

## What Session 107 Itself Says

`Research/Session107_DESI_Forecasts.md` (committed 2025-12-10) explicitly acknowledges DR1 release in its Timeline section:

```
### DESI Year 1 (Released 2024)
- BAO: Released, consistent with ΛCDM ✓
- RSD (fσ8): Analysis ongoing
- Expected: ~3σ discrimination if Synchronism correct
```

By Dec 2025, DESI's full-shape RSD measurements had been publicly available on arXiv for 13 months (Adame+2024, arXiv:2411.12021, Nov 2024). The session author was aware of DR1 and explicitly framed RSD as "analysis ongoing" — but the relevant analysis (the one S645 cited) was published a year earlier.

This means:
- The author **was aware DR1 BAO was out** (cited as "Released, consistent with ΛCDM ✓")
- The author **claimed RSD was ongoing** when, in fact, full-shape RSD with fσ₈ tables was already on arXiv
- The fσ₈ ~ 0.418 prediction was committed *after* the DESI 2024 V tables existed

**The proposal's Case A (knew DR1 was out) is partially confirmed**; Case B (genuinely blind) is hard to defend given the timeline. Whether the *specific* RSD tables in DESI 2024 V were consulted is unverifiable from the archive, but the broader fact — DR1 was out and known to the author — is established.

## Correcting S645's Framing

S645 wrote: "first hard external falsification of a Tier-1 Synchronism prediction. Prior audits found internal disconnects... S645 finds the universe disagreeing with a numerical prediction the framework committed to."

This needs revision. The corrected framing:

> S645 documented an **internal-consistency failure between the Synchronism framework and already-published DESI DR1 data**. Session 107's fσ₈ prediction is post-DR1 (committed 2025-12-10, ~13 months after the DESI 2024 V RSD paper). The 2.4σ tension is real, but its epistemic status is **post-hoc consistency check**, not prospective falsification.

The numerical disagreement remains: fσ₈(z=0.51) predicted 0.418, DESI measured ~0.55±0.06, kill criterion (>0.45) met. The framework cannot reproduce the data its predictions were committed against. But this is "framework parameters can't reproduce known observations" — important, but a weaker claim than "framework made a blind prediction that turned out wrong."

## Distinction Between the Two Failure Modes

| Failure mode | Epistemic status | What it means |
|--------------|------------------|---------------|
| Prospective falsification | Strong evidence against framework | Prediction was made blind; data ruled it out |
| Post-hoc consistency failure | Internal inconsistency | Framework's rules can't reproduce known data |

Both are valuable. Both should be documented. But they are not equivalent. Conflating them inflates the framework's epistemic standing in either direction (as either "validated against blind tests" or "refuted against blind tests").

## Implication for S645 and Memory

S645 should be retroactively updated:
- Status of TEST-04a: from "REFUTED — first hard external falsification" to "REFUTED (post-hoc) — consistency check failed at 2.4σ; framework parameters cannot reproduce DR1 fσ₈ measurements committed before Session 107."
- Memory update: the count "14 internal audits + 1 first-class refuted prediction" is too strong. The corrected count: "14 internal audits + 1 post-hoc consistency-check failure (DESI DR1 vs Session 107)."

This is uncomfortable to write, since S645 felt like a meaningful step forward. The proposal is correct that the temporal independence — what makes prospective prediction epistemically strong — is absent.

## Broader Methodological Recommendation

The proposal's broader point applies across all Tier-1 tests: **timestamp the predictions against the data**. For each Tier-1 kill criterion:

1. When was the prediction committed to git?
2. When did the relevant data become public?
3. Is the prediction prospective or post-hoc?

This is a one-pass audit across the framework. The honest taxonomy:
- Prospective: prediction date < data date
- Post-hoc consistency: prediction date > data date but derivation is independent
- Post-hoc fit: prediction date > data date and derivation uses the data

Without this taxonomy applied, the framework's epistemic claims are mixed in a way readers cannot disambiguate. The cleanest fix is a `/timestamps` page that classifies every Tier-1 prediction by these three categories.

## Connection to S646 (Meta-Falsification)

S646 recommended a meta-criterion for framework-level retraction. S648 sharpens what counts as "kill criterion fired":

- **Prospective kill fire**: strong evidence; should weight heavily in any meta-criterion.
- **Post-hoc consistency fire**: internal-inconsistency evidence; should weight, but less heavily than prospective.

A meta-criterion that doesn't distinguish these would over-weight post-hoc events. Operator-level decision, but the distinction matters.

## On the Other Pending Proposal

`dual_coherence_functions_kinematic_bifurcation.md` (filed same day) duplicates the audit S640 already performed. That proposal asked whether C(ρ) and C = f(γ, D, S) reduce to one another. S640 verified D and S are not functions of ρ (Sessions #358-359), the C = 0.50 threshold is from independent 8-way convergence, only γ formula is shared, and Session #251 introduces a third form. **Path C: bridge does not exist as derivation, only as notational claim.** No new session needed; cross-reference S640 in any response.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 15 | Method-unspecified validation | S647 |
| 16 | **Self-correction of prior session's framing on epistemic basis (NEW)** | **S648** |

This is the first audit instance that turns inward — auditing my own prior session's framing rather than a site-archive disconnect. The visitor channel surfaced an epistemic claim I made in S645 that was too strong. Honest correction is the only response.

## Recommended Site Action

Update TEST-04a entry on `/tier-1-existing` and `/honest-assessment`:

> **REFUTED (post-hoc, Nov 2024 data)**: Session 107's fσ₈ prediction was committed 2025-12-10, after DESI DR1 (April 2024) and after the full-shape RSD analysis (arXiv:2411.12021, Nov 2024). The 2.4σ tension at LRG1 is a post-hoc consistency-check failure, not a prospective falsification. The framework's parameters cannot reproduce the DR1 fσ₈ measurements that were already public when Session 107 was written.

Add a `/pre-registration` page or note documenting the timestamp-vs-data-date for every Tier-1 test.

## Files

- `Research/Session648_Session107_PreRegistration_Audit.md` (this document)

## So What?

S645's "first hard external falsification" framing was too strong. The proposal correctly identifies the temporal independence gap. Session 107's fσ₈ prediction is post-DR1, not prospective. The 2.4σ disagreement is real but its epistemic status is post-hoc consistency failure, not blind-test refutation.

This is uncomfortable because S645 felt like the audit channel finally producing a "clean" finding. The clean finding is now revised: **the framework cannot reproduce DR1 data it was committed against**, which is meaningful but distinct from "the framework made a prediction and the universe ruled against it."

The broader lesson: timestamp every prediction against its data. This is a one-pass audit that should be done across all Tier-1 tests. The audit channel can do meta-work on its own prior conclusions, not just on the framework's claims.
