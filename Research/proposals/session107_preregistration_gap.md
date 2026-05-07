# Proposal: Session 107 Pre-Registration Gap — Post-Hoc Analysis

**Filed**: 2026-05-07  
**Origin**: Site visitor feedback (4-persona pass), Pass 4 Researcher; git timestamp verification  
**Proposal type**: Methodological honesty gap

---

## Finding

Session 107 (DESI Forecasts — Concrete testable predictions) was committed to the Synchronism research repo on **2025-12-10**.

DESI DR1 (Early Data Release) dropped in **April 2024** — approximately 8 months earlier.

This means the fσ₈(z=0.51) ≈ 0.418 prediction was written **after the DESI DR1 data was publicly available**. The comparison that yielded the "2.4σ disagreement" (DESI DR1 measures fσ₈ ≈ 0.55 ± 0.06) is a post-hoc analysis — the framework was applied to a domain where the answer was already known.

**Git evidence:**
- `ba1b653b` — 2025-12-10 06:08:06 — Session #107: DESI Forecasts - Concrete testable predictions
- DESI DR1 public: April 2024

---

## What This Changes

### What the site currently implies

The TEST-04a entry on `/tier-1-existing` and `/honest-assessment` presents the 2.4σ disfavor as a genuine falsification — the framework predicted fσ₈ suppression below ΛCDM, DESI found enhancement above ΛCDM. The framing implies this was a prospective prediction tested by subsequent data.

### What the git timestamps show

The prediction was committed **after** DR1. This does not necessarily mean the author consulted the DR1 data before writing Session 107 — but it means the temporal independence that makes prospective prediction epistemically strong is absent.

The honest framing: Session 107 is a post-hoc calculation showing that, when the Synchronism framework is applied to RSD structure growth, the implied fσ₈ at z=0.51 (≈0.418) disagrees with what DESI DR1 already measured (≈0.55 ± 0.06). This is a consistency check, not a prospective falsification.

---

## Why the Distinction Matters

A prospective prediction that disagrees with subsequent data at 2.4σ is:
- Strong evidence against the framework (the prediction was made blind)
- The standard for scientific falsification

A post-hoc calculation that disagrees with known data at 2.4σ is:
- Evidence of internal inconsistency between the framework and known physics
- Valuable, but at the same epistemic level as "the framework's parameters can't reproduce known observations" — which is important but not the same as prospective falsification

The site's research philosophy explicitly distinguishes prospective tests from post-dictions. Applying that taxonomy consistently means SESSION 107 ≈ POST-DICTION.

---

## Three Cases

### Case A: Session 107 was written with knowledge of DESI DR1 data

If the author consulted DESI DR1 before writing Session 107, the prediction is straightforwardly post-hoc. The disagreement is a consistency failure, not a falsification.

### Case B: Session 107 was written without consulting DESI DR1 (blind)

The temporal gap (8 months) doesn't preclude this. If the framework derivation in Session 107 was done without looking at DESI DR1 tables, the prediction has stronger evidential force even though the commit is post-DR1. This is verifiable: does Session 107 cite or discuss DR1 results?

### Case C: The prediction methodology is what matters, not the timestamp

One could argue: the framework's RSD prediction can be derived from Session 107's first principles regardless of when DR1 was available. If the derivation is self-contained and doesn't use DR1 as input, the disagreement is genuine internal consistency failure. The timestamp question is then secondary to the question of whether the derivation is independent.

---

## Proposed Actions

### Immediate (site fix)

Update TEST-04a on `/tier-1-existing` and `/honest-assessment` to add:

> **Pre-registration note**: Session 107's fσ₈ prediction was committed 2025-12-10, approximately 8 months after DESI DR1 dropped (April 2024). This is a post-hoc calculation showing internal framework inconsistency with known data — not a prospective falsification. The 2.4σ tension is real but its epistemic status is that of a consistency check.

### Research task

Read Session 107 to determine: does it cite or reference DESI DR1 results? If yes → Case A confirmed. If no → assess whether the derivation is self-contained and independent (Case B or C).

### Broader implication

The pre-registration question is the central methodological gap identified by the researcher persona: "0 of 24 run as formal pre-registered tests." The git timestamp approach is exactly the right tool — apply it to every Tier-1 kill criterion. For each:

1. When was the prediction committed to git?
2. When did the relevant data become public?
3. Is the prediction prospective (before data) or post-hoc (after data)?

The honest answer shapes every badge and framing on the prediction pages.

---

## Suggested Session Label

`Session_N: Session 107 Pre-Registration Audit — Prospective vs Post-Hoc`

Task: read Session107_DESI_Forecasts.md, check whether DR1 is cited, assess derivation independence, and write the honest pre-registration classification for every Tier-1 test.
