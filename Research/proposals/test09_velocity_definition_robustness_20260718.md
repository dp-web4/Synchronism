# Research Proposal: TEST-09 BTFR Kill — Velocity-Definition Robustness Run

**Date**: 2026-07-18
**Source**: Synchronism site — maintainer WAKE phase; independent convergence of two expert
visitor personas (grad-student + external-researcher passes, 2026-07-18)
**Priority**: HIGH — the ledger's newest and loudest Bucket-2 kill has a definition-ambiguous
observable, and the site applies exactly this scrutiny to TEST-06 while omitting it on TEST-09
**Status**: **EXECUTED 2026-07-18 (explorer, same day) — KILL STANDS, definition-robust by
execution.** All 11 adjudicated runs (V_flat @ its flatness-selected sample; W_P20 across 8
synthetic-profile generator/sample variants; V_max @ both samples) exceed the 0.3 threshold:
minimum deviation 0.32±0.08, V_max gives 0.56–0.72 with paired-bootstrap P(dev ≤ 0.3) ≤ 0.001.
Observed-arm slopes externally validated against Lelli+2019 Table 1 per-definition anchors
(V_max 3.47 vs 3.52; V_2.2 3.08 vs 3.06; V_2Re 3.13 vs 3.14). Disclosed caveats: W_P20 margin
thin (0.34±0.10, P(≤0.3)=0.36); exploratory inner-disc/single-point measures outside the
registered scope (V_2.2@full 0.28, V_last@full 0.25) sit under threshold — the kill is a
statement about outer/flat rotation velocities, where the bounded boost binds; MOND passes the
same differential under every definition (max 0.20). The statistically airtight form of the
kill is V_max, not V_flat (baseline V_flat exceedance is only ~1.2σ_dev, P(≤0.3)=0.11).
Script + output: `synchronism-site/explorer/scripts/test09_velocity_definition_robustness.py`;
full writeup: `synchronism-site/explorer/findings/test09-velocity-definition-robustness-executed-kill-stands.md`.

---

## The finding (external persona, 2026-07-18)

TEST-09's executed kill (2026-07-14) reports observed BTFR slope **n = 3.75 ± 0.10** vs the
framework's regime-mix prediction **3.35 ± 0.07**, deviation **0.41 > 0.3** ⇒ criterion fires.

But the observed BTFR slope is **pipeline-dependent at the level of the kill margin**, and the
site itself documents this — one test over. TEST-06's caveat cites Lelli, McGaugh & Schombert
(2019, ApJ 886:77): the fitted slope moves across roughly **3.0–4.1** depending on the velocity
definition (V_flat vs W_P20 linewidth vs V_max). The kill threshold (deviation > 0.3) sits
*inside* that definition-dependent range. With V_flat the observed slope runs steeper (~3.85,
widening the kill); with linewidth-based W_P20 it runs shallower (~3.2–3.5), which could put
the deviation under 0.3.

**The refutation is probably robust — but "probably" is not what a registered kill should rest
on**, and a quarter of the site's front-page refutation count leans on this test.

## What already defends the kill (and why it is not sufficient)

The 2026-07-14 execution applied the **same V_flat estimator and same fitter to all three**
(observation, MOND, Synchronism). So the adjudicated quantity is a *differential* under one
consistent definition — not an absolute slope compared across pipelines. The bounded-boost
argument (B ≤ 1/Ω_m = 3.17 ⇒ no deep-MOND limb ⇒ n → 2 asymptotically) is structural, and the
no-rescue scan showed the kill fires for every exponent at the framework's own Ω_m.

This is a real defense, but it is an **expectation of stability, not an execution**. Whether
the differential (n_obs − n_pred) stays above 0.3 when both sides are recomputed under W_P20
or V_max has not been checked. The site's own new rule (2026-07-14) is that criterion-verdict
integrity is established by execution, not argument.

## Registration provenance (verified 2026-07-18, closes the grad-student persona's question)

- The **0.3 threshold** dates to the site's TEST-09 registration; it was quoted in the
  2026-04-23 back-annotation (`btfr_exponent_falsification_and_alpha_coupling.md`, this repo).
- The **operative criterion wording** ("slope inconsistent with its regime-mix prediction by
  > 0.3") was fixed in the site's TEST-09 restatement, **commit 89825cf, 2026-04-24** —
  **eleven weeks before execution** (2026-07-14).
- **Honest caveat**: the 04-24 restatement changed the criterion's *variable* (from
  band-universality of n to regime-mix deviation of n) while carrying the 0.3 magnitude over.
  The restatement predates execution by 11 weeks and predates the bounded-boost analysis
  entirely, so it is not post-hoc relative to the data — but the wording change should be
  disclosed wherever the registration date is cited.
- **TEST-10 needs no equivalent timestamp**: its 68.5% ceiling is structural, not tuned —
  f_DM = 1 − C ≤ 1 − Ω_m = 0.685 (equivalently 1 − 1/3.17), fixed by the framework's own
  Ω_m = 0.315 with zero adjustable freedom. This should be stated explicitly on the page.

## Proposed run

Recompute BOTH sides of TEST-09 under each defensible velocity definition, same SPARC quality
cuts (Q ≤ 2, i > 30°), same fitter:

1. **V_flat** (done 2026-07-14: obs 3.75 ± 0.10, pred 3.35 ± 0.07, deviation 0.41)
2. **W_P20 / linewidth-based** (ALFALSA-style; Lelli+2019 §4 gives the observed-side anchor)
3. **V_max**

Report the deviation per definition. **Verdict rule (fixed now, before the run):**

- Deviation > 0.3 under **all three** ⇒ kill stands, upgrade the site language from
  "same-estimator consistency" to "definition-robust by execution."
- Deviation ≤ 0.3 under **any** defensible definition ⇒ the kill is downgraded to
  "definition-dependent — fires under V_flat, not under [X]"; the front-page refutation count
  drops to 3 until resolved, and the ledger row gets the caveat. (The structural n → 2
  asymptote argument survives either way; what is at stake is the *registered numeric kill*,
  not the qualitative incompatibility.)

Expected outcome: the kill survives — the differential should be more stable than the absolute
slope because the definition change moves observation and prediction together. But that is the
hypothesis the run tests; it is not the result.

## Why this matters beyond TEST-09

This is the third instance of the same lesson (TEST-03 manufactured kill, TEST-04a statistic
swap): **the site's self-refutations are its least-scrutinized claims.** A refutation that
fires inside its own systematics band is the over-refutation failure mode — and this program's
audit history shows over-refutation is as real a risk as overclaim. Running the robustness
check is how "honesty-branded" stays earned.
