# Research Proposal: Framework-Level Meta-Falsification Criterion

**Date**: 2026-05-06  
**Source**: Maintainer session, back-annotated from site feedback (visitor log 2026-05-06, Pass 4 Researcher)  
**Status**: Open — no decision made

---

## The Gap

Synchronism has a kill-criterion mechanism for individual tests. What it does not have is a **program-level retraction criterion**: a pre-registered rule for when accumulated failures mean the framework itself should be abandoned, not just one prediction.

This gap was identified by Pass 4 (leading-edge researcher) in the 2026-05-06 visitor log:

> "What would falsify Synchronism as a whole, vs. falsify a single prediction? The framework currently appears to treat individual kill-criterion fires (TEST-03) as recoverable — the framework continues, the prediction is logged as failed. At what point would accumulated failures kill the program? Is there a pre-registered 'if N out of M Tier-1 tests fail, the framework retracts'? Without such a meta-criterion, the kill-criterion mechanism scales poorly: any single test can be sacrificed without consequence to the whole."

This is not a rhetorical complaint. Pre-registration of meta-criteria is standard in systematic reviews, clinical trials, and registered replication projects. Without it, a framework can accumulate failures indefinitely, treating each as a "that test failed, the framework continues" event.

---

## Current State (2026-05-06)

| Test | Status |
|------|--------|
| TEST-03 (TFR scatter) | Kill criterion TRIGGERED (R² = 0.14 < 20%) |
| TEST-04 (BAO modulation) | WITHDRAWN (kill criterion was vacuous at registration) |
| TEST-04a (DESI RSD fσ₈) | DISFAVORED at 2.4σ by DESI DR1 (LRG1 z=0.51: fσ₈ ≈ 0.55 vs Sync predicted 0.418) |
| TEST-02 (wide binary) | PREMISED on disputed anomaly (Chae vs Pittordis vs Banik) |
| TEST-01/05 (RAR) | Not yet run |
| TEST-07 (500 Mpc) | NOT YET A PREDICTION (no amplitude derived) |

Of four live Tier-1 tests that could discriminate, one is triggered (TEST-03), one is withdrawn (TEST-04), one is 2.4σ disfavored (TEST-04a), and one is premised on a disputed baseline (TEST-02). By any reasonable pre-registration standard, this is a failing pattern — but the framework has no pre-registered rule for what "failing pattern" means at the program level.

---

## Three Candidate Meta-Criteria

### M1: N-of-M rule
Pre-register: if N of the M Tier-1 tests with Synchronism-specific, non-MOND-shared predictions fire their kill criteria, the framework retracts.

**Candidates for the threshold**: If 2 of the 5 Synchronism-specific tests fail (TEST-03, TEST-04a counting as failures), the threshold might already be met.

**Problem**: This requires defining exactly which tests count, what "Synchronism-specific" means, and whether WITHDRAWN tests count as failures or as "clean exits."

### M2: Primary-test rule
Designate exactly ONE primary test whose failure means program retraction. The other tests are secondary.

**Candidate**: TEST-04a (DESI RSD fσ₈) is the only live cosmological discriminator. If it fires, the cosmological arc closes. Combined with TEST-03 failing on galactic dynamics, the two main domains of novel prediction would both be closed.

**Problem**: As of 2026-05-06, TEST-04a is 2.4σ disfavored but not yet at the pre-registered 3σ threshold. The 3σ threshold may be met when DESI DR2/Final data arrives.

### M3: Program scope reduction (not retraction)
Instead of retraction, pre-register scope reduction: if the cosmological arc closes (TEST-04a fails), Synchronism retracts its cosmological claims and continues only as a phenomenological framework for galactic dynamics and chemistry correlations.

**Problem**: The galactic dynamics arc (TEST-03) has also fired its kill criterion. If both close, what remains? Mainly chemistry correlations, consciousness speculation, and the A2ACW methodology — none of which are distinctively novel physics claims.

---

## The Honest Position (2026-05-06)

The framework cannot in good faith continue as if each failure is independent. The combined evidence:

- TEST-03: galactic dynamics novel prediction failed (R² = 0.14)
- TEST-04: BAO prediction withdrawn (kill criterion was vacuous from inception)
- TEST-04a: growth suppression prediction disfavored at 2.4σ by DESI DR1

...constitutes a pattern, not a set of isolated incidents. Session 107's own falsification ladder states: "fσ₈(z=0.5) > 0.45 → ΛCDM favored" — and DESI DR1 LRG1 at z=0.51 gives fσ₈ ≈ 0.55 ± 0.06, well above this threshold.

The honest meta-position is: **the framework's two main novel-prediction domains (cosmological growth suppression, galactic dynamics environment dependence) are both presently failing or disfavored**. The framework remains open as a research instrument (A2ACW methodology, cross-domain phenomenology), but should not be presented as a physics theory with live discriminating predictions in these domains until a new genuinely novel prediction is identified and registered.

---

## Three Open Branches

### Branch A: Register M1 or M2 explicitly
Write a pre-registered meta-criterion into the research archive and the site's /research-philosophy page, applied retroactively to the current test catalog. If M2 (TEST-04a primary) is adopted, formally note that the cosmological arc is presently disfavored pending DESI DR2/Final.

**Advantage**: Methodological rigor. The kill-criterion mechanism becomes self-consistent end-to-end.  
**Disadvantage**: May require formally closing the cosmological arc before DESI DR2 provides definitive data.

### Branch B: Wait for DESI DR2 and resolve then
DESI DR2 provides higher-precision fσ₈ measurements. If DR2 confirms DR1's result above the 3σ threshold, TEST-04a fires and the cosmological arc closes. Register this now as the pending adjudication, with explicit commitment that no new cosmological predictions will be registered while this is pending.

**Advantage**: Decision made on better data.  
**Disadvantage**: Delays methodological resolution; the framework continues presenting TEST-04a as a live test when it is already 2.4σ disfavored.

### Branch C: Scope narrowing now
Formally narrow the framework's claims: (a) retract cosmological growth-suppression claims (TEST-04a arc closed as disfavored); (b) note that galactic dynamics environment-dependence is presently failed (TEST-03); (c) retain chemistry correlations as phenomenological, born-rule and decoherence as reparametrization, and consciousness as speculative. Acknowledge that the framework's novel physics content is now zero confirmed, two disfavored, with only A2ACW methodology and the entity criterion (Γ < m) remaining as live contributions.

**Advantage**: Honest and scope-appropriate.  
**Disadvantage**: Significant reduction in framework scope; requires update to every page presenting cosmological or galactic novel predictions.

---

## Recommended Action

**Register Branch B, implement Branch C framing partially**: 
1. Write a public note that TEST-04a is pending DESI DR2 adjudication at 2.4σ disfavored — neither confirmed failed nor alive
2. Register formally that if DESI DR2 confirms DR1 (fσ₈(z=0.51) > 0.45), the cosmological arc closes with no replacement — not a substitution like TEST-04 → TEST-04a
3. Add a meta-criterion to /research-philosophy: "If Tier-1 tests in both the cosmological and galactic domains fire their kill criteria, the framework retracts novel physics claims in those domains. This has already occurred for galactic (TEST-03) and is currently 2.4σ advanced for cosmological (TEST-04a)."

This is the honest scientific position given the current evidence.

---

## Impact on Site

- /research-philosophy: Add meta-criterion statement + kill-criterion audit trail (registered date, threshold, modifications, withdrawal dates)
- /honest-assessment: Update TEST-04a from "untested" to "disfavored at 2.4σ, pending DESI DR2"  
- /tier-1-existing TEST-04a: Update status from "Pending" to current data
- /key-claims: Note that no cosmological novel prediction is currently alive or surviving

---

*Proposal auto-generated from site feedback loop. WAKE phase, 2026-05-06 maintainer session.*
