# Triage — B3's "C≈0.50 refuted against coherence data" was salience data relabeled "C"; re-scoped + φ−1 corrected (2026-07-07)

**Status:** `[ACTIVE-MRH]` — triage of `Research/proposals/session63_methods_scope_b3_wording.md` (site-explorer
Session-63 methods audit), which tripped the hold-checklist (new in-lane ledger-correction proposal). This is
the **third** ledger-citation defect in ~12 days (after dim-4 c_μν and the a₀ row) — but the FIRST that cuts
toward *un-refuting*, so it got over-failing-grade scrutiny in **both** directions. **Verdict: ACCEPT — verified
against the S63 primary source and by re-executing the statistics. B3's "value REFUTED against multi-model
coherence data" mischaracterizes a test of a *different variable* (`salience_total`, a SNARC heuristic the script
renames "C"); honest status is unsupported & untested. Balanced, not a pure softening: the same audit refutes the
0.64≈φ−1 reparametrization on the data's own terms. Bucket 0 unchanged (0); arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Independent verification (re-execute, don't re-read)

Both source repos are local (`gnosis-research/`, `SAGE/`). Verified the three load-bearing claims:

1. **The tested variable is `salience_total`, not coherence.** `thor_session_63_cross_instance_c_validation.py`
   builds its metric from `exp.get('salience')` (lines 44, 57–61), and its analysis function is literally titled
   *"Analyze whether C (salience_total) clusters around 0.5"* (line 78–79) — the script **renames salience_total
   "C" in situ.** `salience_total` is the weighted mean of five hand-coded SNARC heuristics
   (surprise 0.25 / novelty 0.25 / arousal 0.20 / conflict 0.15 / reward 0.15) from one shared scorer
   (`SAGE/sage/attention/experience_salience.py`) across 8 SAGE instances. So "multi-model coherence data" is
   SNARC-salience data relabeled. **Confirmed in code.**
2. **The t-test rejects an operating *mean* = 0.5, not a threshold location.** `validation_results.json`:
   `t_statistic 20.19`, `p_value 1.83e-7`, `mean_c 0.6400`, over `instances_analyzed 8`. Reproduces the
   proposal's numbers exactly. The 8 "instances" are different models (qwen3.5, tinyllama, gemma3-12b, phi4-14b…)
   but a **shared** scorer, so the tight σ=0.018 is the scorer's, not cross-model coherence convergence.
3. **The 0.64 ≈ φ−1 reading fails on the data's own terms.** Re-ran a one-sample t-test on the 8 instance means:
   H0 μ=0.5 → p=1.8e-7; **H0 μ=φ−1=0.618 → t=3.18, p=0.016 (excluded)**; H0 μ=2/3 → p=0.006 (excluded);
   H0 μ=0.640 → p=0.99 (the sample mean). So even taken at face value the "golden-ratio threshold" is rejected,
   and gnosis Session 64's "validated" verdict is contradicted by S63's own aggregate.

## Why accepting this is honest, not soft (the both-directions check)

The ledger tilts mildly toward **over-failing** (documented: dim-4, a₀). A correction that moves "REFUTED" →
"untested" runs *with* that tilt, so it needs the same skepticism an over-claim gets. Applied it:

- **Steelman keeping "refuted":** the companion program that owns B3 set out to "validate the C≈0.5 hypothesis"
  and, on its own operationalization, rejected 0.5. From its standpoint, 0.5 is refuted.
- **Why re-scope anyway:** that operationalization is `salience_total` — a SNARC heuristic with **no defined map**
  to the framework's coherence C (the framework's own "doubly unanchored" verdict). Renaming a variable "C" in a
  script does not make it C. So the test refutes "mean SNARC salience = 0.5," which is preserved verbatim; it
  does not test the coherence threshold, which therefore stays **untested** — there is no valid C-refutation being
  hidden. Discipline-4 analog holds: the genuine empirical result is kept; only its mislabel as a *coherence*
  refutation is removed.
- **Not a pure softening:** the same audit *adds* a refutation (0.64≈φ−1 excluded at p=0.016) — an over-claim
  correction pointing the other way. Net: one mislabeled refutation re-scoped down, one reparametrization claim
  refuted. Balanced.

## The fix applied to PREDICTIONS.md (B3 status cell)

Replaced "Threshold value REFUTED; … rejected at p<0.0001 against multi-model coherence data … C≈0.64≈φ−1
reparametrization candidate" with: **"UNSUPPORTED & UNTESTED — the one empirical test ran on a different
variable"**, spelling out salience_total / shared-scorer / the script renaming it "C"; preserving the salience
rejection (p≈1.8e-7) as a salience result that doesn't bear on C; and adding that φ−1 and 2/3 are excluded on
S63's own means. Noted the site-side "0.64 also rejected" fabrication was never in this ledger (maintainer P0,
site lane).

## Lane routing

- **PREDICTIONS.md B3** — fixed (mine; the ledger is the source of truth).
- **Whitepaper sections** the proposal also lists (appendix C banner + P280.2 row, conclusion, life_cognition) —
  carry the same "coherence data / rejected p<0.0001" wording. **Flagged for maintainer** (whitepaper is auto-built
  source; consistent with every prior propagation triage — I fix the ledger, route the downstream). They must get
  the same scope-qualification or they will contradict the corrected ledger.
- **Site-side fabrication removal** ("0.64 also rejected at p<0.0001", a 2026-06-23 visitor-persona invention) —
  maintainer P0, already in progress per the proposal; not my lane, recorded for the audit trail.
- **Proposal #2 (`a2acw_cross_vendor_corpus_control.md`)** — a well-posed cross-vendor control that would separate
  the A2ACW null's monoculture-vs-methodology readings, but explicitly **gated on dp credentials / API access to a
  different-corpus model.** Sound design; forward/fleet work; **routed to dp**, no ledger change now.

## So what

Third citation-walk catch in ~12 days, and the most instructive: it points the *opposite* way from dim-4/a₀
(un-refuting, not un-claiming), which is exactly why it needed the two-directions test — and it passed, because the
"refutation" was of a renamed variable, not the coherence threshold. The load-bearing tell was in the source: a
script that titles its own function *"C (salience_total)"*. Lesson reinforced: **re-execute against the primary
source** — the mischaracterization was invisible from the ledger text and only showed up on reading the code that
generated the p-value. The B3 threshold is honestly back to untested (where a no-lab consciousness prediction
belongs), the golden-ratio reparametrization is refuted on its own data, and no genuine result was softened.
Bucket 0 still 0; arc AT REST.
