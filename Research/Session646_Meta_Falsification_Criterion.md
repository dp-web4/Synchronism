# Session 646: Framework Meta-Falsification Criterion — Methodology Recommendation

**Date**: 2026-05-07
**Type**: Methodology (governance-level, not audit-channel)
**Trigger**: 2026-05-06 proposal `framework_meta_falsification_criterion.md`
**Grade**: B+ (concurs with proposal; clarifies state with S639/S645 detail)

---

## Setup

The proposal flags a methodological gap: Synchronism has per-test kill criteria but no **pre-registered framework-level retraction criterion**. Without one, accumulated failures can be treated as isolated, never adding up to anything that would retract the framework. The proposal asks for one of three responses (Branch A: register a meta-criterion now; Branch B: wait for DESI DR2 then decide; Branch C: scope-narrow now).

S646 confirms the state with sharper detail from prior sessions, agrees with the proposal's diagnosis, and clarifies what the worker channel can and cannot do here.

## Confirming the Test State

The proposal's test-state table needs one correction and one update from prior sessions:

| Test | Proposal Status | Corrected Status | Source |
|------|----------------|-------------------|--------|
| TEST-03 | "Kill criterion TRIGGERED" | **Two tests under one ID**: TEST-03A (TFR residual on BTFR, 51%, **passing** by 2.5×) + TEST-03B (RAR environmental ansatz, R²=0.14, **below threshold**) | S639 disambiguation |
| TEST-04 | WITHDRAWN | WITHDRAWN | as proposal |
| TEST-04a | "DISFAVORED at 2.4σ" | **REFUTED** — DESI DR1 fired the kill criterion (fσ₈(z=0.51) ~0.55 vs threshold 0.45) | S645 (today) |
| TEST-02 | Premised on disputed anomaly | Same | as proposal |
| TEST-01/05 | Not run | Same | as proposal |
| TEST-07 | Not yet a prediction | Same | S632 |

**Updated framing**: The cosmological arc (TEST-04a) is now at the kill threshold the proposal anticipated for DR2 — DR1 already fires it at LRG1 (kill = fσ₈(z=0.5) > 0.45; observed = ~0.55 ± 0.06, met at 2.14σ). The proposal's Branch B (wait for DR2) is therefore *less needed* than the proposal anticipated; the kill criterion is met now, not pending DR2.

The galactic arc has TEST-03B failing (RAR environmental ansatz). TEST-03A passes — but it is the M/L-driven TFR-residual prediction, which is shared with MOND and not Synchronism-distinctive (per S637 cumulative finding that the framework reduces to MOND in the testable regime).

## The Proposal's Diagnosis Is Correct

A framework that treats every per-test kill firing as recoverable is a framework with no retraction condition — which is to say, it cannot be falsified at the framework level by accumulating evidence. This is the same self-sealing structure S621 identified ("Intent is pre-mathematical, every fix IS known physics"), now visible from a different angle: at the meta-level, there is no rule for when "the framework has lost."

Pre-registered meta-criteria are standard in clinical trials and registered replication projects for exactly this reason. Without one:
- Each test becomes individually disposable.
- The framework can persist through any number of failures.
- The "productive failure > safe summaries" stated value applies to individual tests but not to the program.

The proposal's three meta-criteria (M1 N-of-M, M2 primary-test, M3 scope-reduction) are all defensible. M3 (scope reduction triggered when both cosmological and galactic domains fail) is the most honest given the current state: TEST-04a refuted (cosmological domain), TEST-03B below threshold (galactic environment domain).

## What the Worker Channel Can and Cannot Do

This is a **governance-level recommendation**, not an audit-internal finding. The worker session can:

- Confirm the state (done above).
- Endorse the proposal's reasoning.
- Recommend adoption.

The worker session cannot:

- Unilaterally adopt a meta-criterion. That requires operator-level commitment.
- Update `/research-philosophy` or other site governance pages (those are operator-domain).
- Retract the framework's claims at scope.

The recommendation is for the Publisher / operator track to act on. From the audit-channel side, S646 affirms the diagnosis and the recommended action; the actual adoption is upstream.

## Recommended Action (Operator Decision)

The proposal's hybrid Branch B + Branch C is sound, but DR1 already fires the kill (S645). The cleanest action sequence:

1. **Register a meta-criterion now** (Branch A): explicit text in `/research-philosophy` stating that if Tier-1 tests in both the cosmological and galactic novel-prediction domains fire their kill criteria, the framework retracts novel physics claims in those domains. Include date stamps and audit trail.

2. **Implement scope narrowing** (Branch C, partial): given that TEST-04a is refuted (S645) and TEST-03B is below threshold (S639):
   - `/key-claims`: note that no cosmological novel-prediction is currently surviving.
   - `/honest-assessment`: TEST-04a → REFUTED with date.
   - Galactic-dynamics novel claim: scope to TEST-03B status (below threshold for environmental component).
   - Retain TEST-03A (TFR-residual, M/L-driven, MOND-compatible) as a useful tool, but not as a Synchronism-distinctive novel prediction.

3. **Surviving framework content**: A2ACW methodology, entity criterion (Γ < m), and the cross-track audit/perseveration-recognition pattern as methodology contributions. None of these are novel physics predictions; all are honest about what they are.

This converges with the picture from S637 (framework reduces to MOND in testable regime), S638 (Curie response in CM regime), and S645 (one external falsification). The meta-criterion makes the program-level conclusion explicit rather than implicit across 14 audits.

## Audit-Channel vs Methodology

The audit channel (S631-S644) finds and resolves site-archive disconnects. S645 reports the first external falsification. S646 sits at a different layer — it is the **methodology recommendation that the audit findings imply**: with TEST-04a refuted and TEST-03B below threshold, the framework's novel-physics scope should be formally narrowed.

This is not a 16th audit instance. It is the methodology-level synthesis that 14 audits + 1 falsification motivate.

## Files

- `Research/Session646_Meta_Falsification_Criterion.md` (this document)

## So What?

The framework lacks a pre-registered retraction criterion. This is a methodological gap that becomes load-bearing now that one Tier-1 prediction is refuted (TEST-04a, S645) and another component is below threshold (TEST-03B, S639). The proposal's Branch A (register meta-criterion) + Branch C (scope-narrow) is the honest action sequence.

The worker channel's role here is to confirm the state and recommend the action; the action itself is operator-level. S646 logs the recommendation. Whether and how it is implemented is upstream.

One proposal remains pending: chemistry validation N_corr method unspecified (2026-05-06).
