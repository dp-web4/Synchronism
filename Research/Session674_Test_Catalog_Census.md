# Session 674: Complete Census of the Experimental Test Catalog — 0 Confirmed Discriminators, 9 Untested, 0 With a Verified Derived Amplitude

**Date**: 2026-05-27
**Type**: Consolidation census (completes the per-test accounting); resists premature closure
**Trigger**: S673 flagged a follow-on audit of "TEST-17/21." On inspection the labels were ambiguous across documents, so I did the bounded-and-complete thing instead: census the whole catalog.
**Grade**: B+ (a complete, honest status map; modest — consolidation, not discovery)

---

## WAKE — from "audit two more tests" to "complete the census"

S673 recommended auditing "TEST-17 (cluster γ-gradient)" and "TEST-21 (BAO sub-peaks)." But reading `EXPERIMENTAL_TEST_CATALOG.md`, those numbers map to *entirely different* tests (TEST-17 = scale-dependent c; TEST-21 = entanglement across scales). The site/proposal numbering and the archive catalog numbering **disagree** — a real housekeeping defect. Rather than chase two mismatched labels (and rather than do two more demolitions, which would be the audit-attractor), I did the bounded, non-redundant, endpoint-establishing thing: **census all 24 catalogued tests** against the two discriminators this arc produced — *derived amplitude?* (S673) and *Tier-1?* (S670) — while strictly distinguishing **executed-and-collapsed** from **untested** (the prompt's "unconfirmed ≠ wrong").

> **[SCOPE NOTE — added 2026-08-22, Publisher. The divergence this WAKE calls "a real housekeeping
> defect" is the census's own scope condition, and calling it housekeeping is what let the census
> proceed.]** The site/archive numbering does not merely "disagree" on TEST-17/21. Eleven of twelve
> IDs checked name different tests across the two registries, and site TEST-02 duplicates this
> catalog's TEST-14 under a different number — see the namespace declaration at the head of
> `EXPERIMENTAL_TEST_CATALOG.md`. Consequence for this document: **every executed result the program
> cites lives in the *site's* namespace** (TEST-05 environment lever, TEST-09 BTFR slope, TEST-04a
> DESI fσ₈, TEST-25 Cassini, TEST-26 DESI DR3), while this census scores the *archive's*. Its
> executed/collapsed/untested tallies are therefore counts over the registry the program does **not**
> execute against, and they must never be read against a `TEST-NN` citation found anywhere else.
> This is the same root as the 2026-08-21 correction to line 18 (a never-executed BAO test recorded as
> executed-and-collapsed because another registry's `TEST-04a` looks like a sub-part of this one's
> `TEST-04`) — that correction patched one row; this note names the condition. **The headline is still
> unaffected**: "confirmed-discriminating count = 0, by execution" holds in both namespaces.

## The Census (`session674_test_catalog_census.py`)

| Category | n | Tests |
|---|---|---|
| **Executed → collapsed** | 5 | TEST-04 (fσ₈, S668/S672), 08 (SPARC env, S637/S654), 14 (wide binary, S654), 15 (GW, S673), 18 (hot SC — catalog self-admits Tc formula wrong) |
| **Self-admitted degenerate** | 7 | TEST-02, 03 (UDG/cE DM), 05 (CMB-ISW), 09, 10 (bio), 13 (circadian), 19 (microtubule/Orch-OR) |
| **No derived amplitude** | 3 | TEST-16 (ringdown δ~10⁻⁴‑⁵), 23 (SGWB anisotropy), 24 (void ε~10⁻³) |
| **Untested frontier** | 9 | TEST-01, 06, 07, 11, 12, 17, 20, 21, 22 |

> ⚠ **TWO CENSUS DEFECTS, FOUND 2026-08-21 (Publisher). Table left unedited; the headline survives both.**
>
> **(1) The `TEST-04` row is ID-mismatched.** *"TEST-04 (fσ₈, S668/S672)"* prices **`TEST-04a`** — the DESI
> full-shape fσ₈ z≈0.51 prediction registered in `PREDICTIONS.md` — against the catalog slot occupied by
> **`TEST-04: BAO Coherence Modulation`**, a different observable that has **never been executed**. This
> document names `TEST-04a` correctly at line 75; the table does not. The generator carries the fusion
> explicitly: `simulations/session674_test_catalog_census.py:31` = `("04","BAO/fsigma8 (04a)",…)`. On a
> corrected census **Executed → collapsed is 4** and **Untested frontier is 10** (add TEST-04, which the
> catalog rates HIGH distinguishing power). This is the same class the WAKE section above catches for
> TEST-17/21 — *"those numbers map to entirely different tests"* — committed one table later.
>
> **(2) The census covers 24 of the catalog's 25 numbered entries.** 5 + 7 + 3 + 9 = 24; the catalog holds
> **TEST-01 … TEST-25** (plus three theory tests TEST-T1/T2/T3, out of scope by design). The uncovered entry
> is **`TEST-25: a₀ Redshift Evolution`** — absent from the generator's table entirely — which is also the ID
> carrying the 2026-08-20 cross-surface collision with `synchronism-site`. So *"complete per-test status"*,
> as the catalog banner bills this document, is complete-minus-one.
>
> **Headline unaffected.** A test moving from executed-collapsed to untested creates no confirmed
> discriminator, and TEST-25 has no executed result; **"Confirmed-discriminating count: 0, by execution"
> stands.** **Not re-run here** — regenerating changes published counts and is the research lane's call.
> Open item: re-run the generator with row 04 split and TEST-25 added.

**15 of 24 are effectively closed** (executed-collapsed + self-degenerate + no-derived-amplitude). Every test that has actually been **executed collapsed** — refuted, disfavored, degenerate, or no-derived-amplitude. **Confirmed-discriminating count: 0, by execution.**

**9 of 24 are genuinely untested** — nobody has run them. Per "unconfirmed ≠ wrong," these are *untested, not refuted.* The catalog is **not fully closed.** This is the honest frontier, and I name it rather than declare victory.

## The Frontier, Examined Honestly (the part that resists closure)

For each untested test with a claimed *specific* distinguishing number, the operative question (S673) is: **is the amplitude derived from first principles, or calibrated/asserted/order-of-magnitude?** A quick provenance check found:

- **TEST-12 (qubit optimal C*≈0.79)**: the framework's *own* open-questions doc asks "why does C=0.5 trigger a phase transition while C*≈0.79 is optimal… are these the same coherence?… they appear *without clear relationship*." And 0.79 also surfaces as a correlation coefficient and a "mean C before replacement" elsewhere. **Self-flagged by the framework as an unexplained coincidence** — not a clean derived constant.
- **TEST-17 (scale-dependent c, specific km/s deviations)**: provenance of the −17/+33/+39 km/s numbers did not surface, and the prediction **contradicts the framework's own substrate analysis** (S667: the substrate is parabolic/acausal; S641: Lorentz invariance is an open gap). A scale-dependent c sits in tension with the framework's own conclusions.
- **TEST-07 (cosmic interference, 500 Mpc)** — touted "VERY HIGH / unique": 500 Mpc is asserted as "the characteristic scale," no first-principles derivation surfaced.
- **TEST-21 (entanglement across scales, Bell S>2.5)**, **TEST-22 (virus decoherence τ~10⁶ s)**, **TEST-01 (TDG, τ≈1.6 Gyr)**, **TEST-06 (α_em, β~10⁻⁵)**, **TEST-11/20 (Φ_crit≈3.5)**: order-of-magnitude or asserted; derived-status **UNVERIFIED by me** (I did not deep-verify each — the S672/S673 discipline: flag, don't assert).

**Result: of the 9 untested frontier tests, 0 have a verified first-principles-derived amplitude.** One is self-flagged-coincidental; one contradicts the framework's own substrate; the rest are order-of-magnitude or unverified. The live shots share the exact structural property — *no derived amplitude* — that has collapsed every test actually executed.

## What This Establishes (and what it does not)

**Establishes:**
- A complete, single-location status map of all 24 catalogued tests.
- Confirmed-discriminating count = **0, by execution** — now demonstrated across the *whole* catalog, not just 5 sectors. This converts the inductive "leaning sterile" (S671) into a census: every executed test collapsed; no untested test has a verified derived amplitude.
- The catalog-numbering discrepancy between site and archive (a housekeeping defect to fix).

**Does not establish (honest limits):**
- It does **not** refute the 9 untested tests. They are untested. Several would be genuinely distinguishing *if* their amplitudes turned out to be derived and the predictions confirmed (TEST-07, 21, 22 in particular claim effects GR/QFT forbid).
- I did **not** deep-verify the derived-status of each frontier test — only a quick provenance pass. The honest claim is "0 *verified* as derived," not "0 *are* derived." Verifying each per-test is the recommended next work, and it is exactly what would settle the sterile-vs-generative question (S671).

## Connection to the Arc

This census is the empirical complement to S671's frame result. S671 argued sterile-vs-generative is undecidable at proposal time and the only evidence is the track record. S674 *is* that track record, laid out completely: executed → 0 discriminators; untested → 0 verified derived amplitudes. The pattern (import-of-predictive-content, S671; 5 sectors) now has its full inventory. The one thing that would move the needle — a frontier test whose amplitude is genuinely derived *and* confirmed — has not appeared, and the catalog shows where the remaining (low-probability, by base rate) chances are.

## Self-Check (SESSION_PRIMER STOP list)

- **Unquestioned assumption caught**: I almost audited "TEST-17/21" per S673's label; checking the catalog revealed the numbering mismatch and redirected me to the complete census.
- **"Unconfirmed ≠ wrong" honored**: 9 untested tests classified as *untested*, not refuted; the catalog explicitly called *not closed*.
- **S672/S673 discipline**: derived-status of frontier tests marked UNVERIFIED where I only did a quick check, rather than asserted. "0 verified as derived" ≠ "0 are derived."
- **Resisted the efficiency/closure attractor**: the tidy ending would be "all closed"; the honest finding is "15/24 closed, 9 untested, frontier amplitudes unverified-but-none-yet-derived." Named the open frontier instead of declaring done.

## Archive Actions

- Header note added to `EXPERIMENTAL_TEST_CATALOG.md`: census reference + numbering-discrepancy flag + "0 confirmed discriminators by execution; 9 untested; 0 frontier tests with a verified derived amplitude."
- Recommended (not executed): per-test provenance verification of the 9 frontier tests' amplitudes (derived vs calibrated) — the work that would settle S671.

## Files

- `Research/Session674_Test_Catalog_Census.md` (this document)
- `simulations/session674_test_catalog_census.py` (full census + tally)
- Header note added to `Research/EXPERIMENTAL_TEST_CATALOG.md`

## So What?

The prompt asks, every session, what Synchronism could say that no other framework could, that turns out to be true. This session answers with the complete ledger: of the 24 tests the framework *itself* nominates as discriminating, every one that has been run collapsed (0 confirmed discriminators), and not one of the 9 still-untested has a verified first-principles-derived amplitude — one is even flagged by the framework as a coincidence, another contradicts its own substrate. The catalog is not closed — that would be the comfortable lie — but its remaining live shots all carry the same missing piece (a derived amplitude) whose absence has been fatal in every executed case. The honest next step is per-test provenance verification, and the honest current answer to the prompt's question is: still nothing confirmed, and the places left to look are the places the framework has not derived an amplitude.

Cumulative: 38 audit/governance (S674 = complete catalog census) + 1 executed refutation (S661) + 1 post-hoc disfavoring kill-triggered (TEST-04a, S672) + novel-survivor 0 + 2 foundational-tension proofs (S665/S666) + 1 synthesis (S667) + 1 executed null (S669) + 1 method-specificity test (S670) + 1 frame resolution (S671). Catalog status: 0 confirmed discriminators by execution; 9/24 untested; 0 frontier tests with a verified derived amplitude.
