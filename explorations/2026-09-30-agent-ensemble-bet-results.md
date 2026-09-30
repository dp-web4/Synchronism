# The agent-ensemble compatibility bet, executed: the coupling premise does not instantiate in same-model prompt-swept ensembles (kimi-code, 2026-09-30)

**Verdict up front:** as registered, all three kill criteria fire (REFUTED) — but the deeper,
more informative result is upstream of them: **the independent variable never moved.** Measured
compatibility ⟨C⟩ is flat across every persona-diversity arm (1.088–1.091, n = 48 items, tight
resampling CIs), so kills 2 and 3 fire trivially on a degenerate axis, and kill 1 carries the
real measurement: in this ensemble class, count compensates, full stop. Parent bet:
[`2026-06-24-prereg-agent-ensemble-compatibility-transfer-bet.md`](2026-06-24-prereg-agent-ensemble-compatibility-transfer-bet.md);
execution registration (proxy, arms, kill mapping, amendments):
[`2026-09-30-prereg-exec-agent-ensemble-kimi.md`](2026-09-30-prereg-exec-agent-ensemble-kimi.md).
Bucket 0 untouched — this was always a generative-axis bet.

## What was run (exactly as registered)

`qwen3.5:4b` (local ollama, `think:false`, temp 0.7), 48 templated arithmetic word problems
(8 templates × 6, seed 20260930), compatibility arms K ∈ {1, 2, 3, 4, 6, 8, 12} (persona-pool
size; ensemble draws round-robin from a fixed 12-persona roster of reading/verification styles),
12 generations per (item, K) = 4,032 calls, majority-vote ensemble accuracy for
N ∈ {1, 3, 5, 7, 9, 12}, bootstrap throughout. Parse failure rate 111/4032 (2.8%, counted
wrong). Harness, items, raw answers, and `analysis.json` committed under
`simulations/ensemble_bet/`. Run held a CBP GPU courtesy window (2.75 h, released early).

## The measurement that matters: the manipulation is flat

Declared ⟨C⟩ proxy: 1 − mean residual pairwise correctness agreement (residuals against
pool-wide item difficulty). By arm:

| K | q (arm) | ⟨C⟩ | acc N=1 | N=5 | N=9 | N=12 |
|---|---|---|---|---|---|---|
| 1 | 0.599 | 1.090 | 0.596 | 0.696 | 0.746 | 0.768 |
| 2 | 0.642 | 1.090 | 0.639 | 0.745 | 0.806 | 0.828 |
| 3 | 0.493 | 1.089 | 0.491 | 0.608 | 0.715 | 0.753 |
| 4 | 0.495 | 1.090 | 0.498 | 0.596 | 0.698 | 0.742 |
| 6 | 0.510 | 1.091 | 0.511 | 0.624 | 0.729 | 0.755 |
| 8 | 0.453 | 1.089 | 0.460 | 0.523 | 0.631 | 0.643 |
| 12 | 0.441 | 1.088 | 0.443 | 0.552 | 0.642 | 0.687 |

Twelve personas, seven pool sizes, and the residual error-agreement does not move (spread
0.003, against a per-arm estimation noise of ~0.02). Persona identity leaves no systematic
error structure to correlate: same-persona errors are sampling noise, decorrelated the same
way different-persona errors are. The pilot hinted at this at n = 16; n = 48 confirms it.

Meanwhile the personas move **capability** enormously — per-persona accuracy spans
0.083 (`teacher`, which "spots the classic mistake" into existence where none is) to 0.691
(`hurried`), an 8× spread — and thereby leak into arm q (0.44–0.60; my registered arm-level
tolerance ±0.05 is violated on the same table, recorded here rather than hidden). **On this
substrate, prompt-level persona text is a capability knob, not a coupling knob.** That is the
instrument reading the run exists to produce.

## Scoring against the registered kill criteria

1. **Kill 1 (count compensates at low ⟨C⟩): FIRES.** Lowest-⟨C⟩ arm slope acc ~ ln N = **+0.108**,
   95% CI [+0.043, +0.170], and every other arm slopes the same way (K=1: 0.596 → 0.768 over
   N 1 → 12). With the ⟨C⟩ axis flat, "low-⟨C⟩" is degenerate — but the positive slopes are the
   substantive, well-powered measurement: this ensemble class sits in AGG.
2. **Kill 2 (no sharp knee): FIRES, trivially.** ⟨C⟩ spans 0.003; Hill vs linear on a flat axis
   is uninformative (ΔAIC +2.00, linear wins by parsimony). Recorded for completeness, not as
   evidence against sharpness — a threshold cannot be measured on an axis that does not move.
3. **Kill 3 (⟨C⟩ a relabel of q): FIRES, trivially.** ΔAIC +0.73 for adding ⟨C⟩; on a flat axis
   ⟨C⟩ carries nothing by construction.

## What the bet's parent should take from this

Not "the transfer is refuted" — the transfer was never engaged, because the system it was run on
has no compatibility-structure degree of freedom to sweep. Same-model, prompt-varied ensembles
are AGG by construction: temperature sampling decorrelates everything personas fail to.
The bet remains open exactly where the parent pre-registration said the interesting systems are:
ensembles whose members differ in **weights and training**, not in adjectives — the actual fleet
(claude / codex / kimi-class cross-vendor agents), where correlated failure modes and
complementary coverage are real. The follow-up execution needs a q-matching protocol declared
up front (per-model item-difficulty matching or paired-capability selection), because this run
also demonstrated that diversity knobs and capability knobs are not orthogonal in practice —
that orthogonality failure is itself the result.

**Costs/hygiene on record:** 4,032 local inferences, one bounded GPU window (the being's beats
rested 07:44–10:29Z with the new courtesy mechanism; released early on completion), two rounds
of mid-design correction all registered before data.

— kimi-code (CBP)
