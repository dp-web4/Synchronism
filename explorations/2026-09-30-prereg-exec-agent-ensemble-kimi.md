# EXECUTION REGISTRATION — running the 2026-06-24 agent-ensemble compatibility transfer bet (kimi-code, 2026-09-30)

**Status:** registered BEFORE the full run; pilot data (16 items, 2 arms) used only to size the
run and confirm the manipulation moves the proxy. Parent bet:
[`2026-06-24-prereg-agent-ensemble-compatibility-transfer-bet.md`](2026-06-24-prereg-agent-ensemble-compatibility-transfer-bet.md)
(CBP-Claude). This document fixes the operationalization choices the parent left open ("pick one,
declare it before measuring"), the task domain, and the exact scoring against the parent's kill
criteria. **Bucket 0 untouched; generative-axis bet.**

## What is run

- **Model:** `qwen3.5:4b` via local ollama, `think:false`, temperature 0.7, top_p 0.95,
  num_predict 48, per-generation seed varied. One model for every agent in every arm — per-agent
  capability q is held fixed by construction and measured per persona (reported, not assumed).
- **Task domain:** templated arithmetic word problems with known integer answers
  (`simulations/ensemble_bet/items.py`, seed 20260930, 8 structural templates × 6 = 48 items;
  trimmed from 80 pre-data to fit the GPU courtesy window — power note: kill-1 slope CIs and
  the ⟨C⟩ residual-agreement estimates tighten as √items, 48 is the declared power point).
  Templates are varied so that reading style changes *which* items fail, not *how many*; the
  pilot checks that.
- **Compatibility manipulation (the independent variable):** persona-pool size
  K ∈ {1, 2, 3, 4, 6, 8, 12}. An ensemble of N agents draws personas round-robin from the first
  K of a fixed 12-persona roster (mixed reading/verification styles). K = 1 is maximal error
  redundancy (every agent the same persona); K = 12 maximal complementarity. The roster and the
  round-robin assignment are in `simulations/ensemble_bet/gen_answers.py`.
- **Answer pool:** POOL = 12 generations per (item, K). Ensembles of size N ∈ {1, 3, 5, 7, 9, 12}
  are bootstrap-resampled from the pool (200 draws per (item, K, N); ties in majority vote break
  to the numerically smallest candidate). Generation totals: 48 × 7 × 12 = 4,032 calls.

## Declared measurements

- **⟨C⟩ proxy (declared, the parent's proxy #1):** residual pairwise correctness agreement.
  Per arm K: item difficulty d_i = pool-wide accuracy on item i; each pool member's correctness
  residual r = correct − d_i; ⟨C⟩ := 1 − mean over pool pairs of Pearson corr(r_a, r_b) across
  items. High residual agreement = redundant errors = LOW compatibility.
- **Collective coherence:** ensemble majority-vote accuracy.
- **Manipulation check:** residual agreement must be HIGHEST at K = 1 and fall with K;
  **q-equality binds at ARM level** (per-arm mean accuracy within 0.05 of each other — pilot
  v2: K=1 0.542 vs K=8 0.552 ✓). Per-persona q spread is reported but does not violate the
  premise: every K > 1 arm mixes the same roster, so persona-level differences wash out at arm
  level. (Amended 2026-09-30 pre-full-run: the original "per-persona span ≤ 0.15" was the wrong
  level of analysis — the kill-3 confound is arm-level q, and that is what is held fixed.)

## GPU courtesy (added pre-full-run)

Generation runs inside a **CBP GPU courtesy window**
(`shared-context/machines/cbp-gpu-windows.md`, landed today): the being's beats rest for a
bounded, self-expiring interval while the batch runs; the being's model is never unloaded by
the mechanism. This was found necessary, not just polite: the two models do not co-fit in
8 GB VRAM, so un-windowed coexistence is model-swap churn, not sharing.

## Scoring against the parent's kill criteria (verbatim mapping)

The bet is **REFUTED** for agent ensembles if any of:
1. at the lowest-⟨C⟩ arm (K = 1 as measured), ensemble accuracy **rises** with ln N — slope
   95% CI entirely above 0 (bootstrap CIs) → count compensates → AGG; OR
2. the ⟨C⟩ dependence is **gradual**: at fixed N = 9, a linear fit of accuracy on measured ⟨C⟩
   (7 K-levels) beats or ties a 3-parameter Hill fit on ΔAIC < 2 in Hill's favour → no sharp
   threshold; OR
3. ⟨C⟩ adds nothing beyond {q, N}: logistic fit accuracy ~ ln N + q̄_persona + ⟨C⟩ vs
   accuracy ~ ln N + q̄_persona, ΔAIC < 2 → vacuous.

**SUPPORTED only if** none of the three fires: slope ≤ 0 at low ⟨C⟩, Hill beats linear by
ΔAIC ≥ 2 with a visible knee, and ⟨C⟩ carries ΔAIC ≥ 2 beyond {q, N}.

## Known limits (declared, not discovered later)

- Seven discrete K-levels give a coarse ⟨C⟩ axis; "sharpness" is relative to a linear null on
  seven points. A knee *between* measured levels will read as gradual — that asymmetry is
  conservative (biases toward refutation).
- Persona-induced diversity is *induced*, not organic; the transfer claim if supported is about
  coupling structure, and personas are this lab's way to vary it at fixed q.
- The GPU is shared with the being (priority to it); generation is sequential. Runtime jitter
  does not enter any measurement.

## After the run

Results doc beside this one (`2026-09-30-agent-ensemble-bet-results.md`), scoring each kill
criterion against the numbers, verdict per the parent's rule, harness + raw answers committed
(`simulations/ensemble_bet/`). Either outcome is reported; that is the point of registering.
