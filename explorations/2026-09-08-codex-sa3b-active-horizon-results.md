# SA-3B results: stopping can preserve the ignorance that caused it

Codex, 2026-09-08. Synthetic instrument experiment; no new physics claim.
[Protocol](2026-09-08-codex-sa3b-active-horizon-charter.md) committed and pushed
at `a09f253a` before implementation/evaluation. All registered policies,
worlds, prices, seeds, and horizons were retained without outcome-driven tuning.

## Bottom line

Two distinct limits appeared. Planning about future learning can improve
acquisition, but a longer planner cannot select a relationship absent from
its model. Supplying that relationship helps, yet does not prevent an early
run of misleading evidence from making stopping self-perpetuating.

**A rational purchase decision under the observer's current model is not a
certificate that the world contains no useful excluded information.** This
is the operational counterpart of SA-3A's identifiability control.

## What was actually tested

An observer predicts a binary target after choosing stop, history retrieval,
noisy sensing, or both. It receives target feedback afterward. It starts
without the world identity, but WITH known candidate likelihoods, noise rates,
costs, and a uniform prior. History is a stored bit, not a dynamical process.
The budget is three acquisition credits per round, not a shared episode bank.

Four persistent worlds: memory matters, sensor matters, neither matters, or
their XOR matters (parity). The restricted model omits parity; the expanded
model includes it from the start. Greedy and exact two-round Bayes planners
are compared with four fixed acquisition baselines. Each cell uses 24 paired
seeds × 64 rounds. This is 768 seed-policy runs / 49,152 decisions, not that
many independent environments or scientific findings.

## Full primary-metric comparison

Mean prediction error + acquisition price per round; lower is better.
All values rounded to four decimals. Oracle entries are analytic expected
utilities for a regime-aware purchaser, not fitted sample scores.

| Policy | Memory | Sensor | Noise | Parity |
|---|---:|---:|---:|---:|
| Restricted, myopic | 0.1797 | 0.3608 | 0.5098 | 0.5100 |
| Restricted, two-step | 0.1804 | 0.3139 | 0.5118 | 0.5162 |
| Expanded, myopic | 0.1693 | 0.3732 | 0.5070 | 0.4186 |
| Expanded, two-step | 0.1588 | 0.3492 | 0.5142 | 0.3742 |
| Always stop | 0.5033 | 0.5293 | 0.4961 | 0.4889 |
| Always retrieve | 0.1396 | 0.5426 | 0.5231 | 0.5543 |
| Always sense | 0.5676 | 0.2655 | 0.5852 | 0.5787 |
| Always buy both | 0.2307 | 0.3212 | 0.6115 | 0.3238 |
| Known-regime oracle | 0.1400 | 0.2600 | 0.5000 | 0.3000 |

No adaptive policy beats the best regime-specific fixed action in this run.
The adaptive task is to find that action without knowing the regime. The
expanded model does not dominate the restricted model: it improves parity
substantially but pays for additional uncertainty in the sensor world.

Paired differences below are mean ± one standard error across 24 seed
differences, NOT confidence intervals or registered significance decisions.

| Comparison (first minus second) | Memory | Sensor | Noise | Parity |
|---|---:|---:|---:|---:|
| Restricted two-step − myopic | +0.0007 ± 0.0007 | −0.0469 ± 0.0199 | +0.0020 ± 0.0021 | +0.0062 ± 0.0026 |
| Expanded two-step − myopic | −0.0106 ± 0.0101 | −0.0239 ± 0.0135 | +0.0072 ± 0.0029 | −0.0444 ± 0.0181 |
| Expanded − restricted two-step | −0.0217 ± 0.0168 | +0.0353 ± 0.0168 | +0.0024 ± 0.0020 | −0.1420 ± 0.0214 |

## The three useful distinctions

### 1. Immediate usefulness is not the value of learning

The registered exact probe starts with probabilities (memory 0.1, sensor 0.1,
noise 0.8). Retrieval is break-even for the current prediction, so myopic
selection stops. Two-round planning retrieves: expected cumulative cost
0.9856 versus 1.0000 for stopping then acting optimally next round.

This can be checked by hand. Retrieval plus current prediction costs 0.5.
Its bit agrees with the subsequent label with probability 0.54, increasing
memory's posterior probability to 1/6. Retrieving next round then costs
0.04 + 0.5 − 0.4/6 = 0.473333…. On disagreement, stopping costs 0.5.
Thus 0.5 + 0.54×0.473333… + 0.46×0.5 = 0.9856.

The benefit is learning what will be worth buying later, not information
needed for the current target alone. It is standard Bayesian decision theory,
not a newly discovered MRH law.

### 2. More planning does not repair a missing hypothesis

In parity, restricted two-step never purchases both and all 24 runs stop
throughout the final 16 rounds. Mean final probability assigned to its noise
model is 0.832. Expanded two-step purchases both in 17 of 24 runs throughout
the tail, reducing mean total utility from 0.5162 to 0.3742.

Parity was supplied in the expanded family before the run. This demonstrates
the consequence of representational support, not spontaneous discovery of
an interaction. Single-channel observations genuinely cannot distinguish
parity from independent noise in these worlds; their joint observation can.

### 3. Even the right family can stop prematurely

Expanded two-step still stops throughout the tail in 7/24 parity runs, 7/24
sensor runs, and 1/24 memory runs. All noise runs eventually stop, appropriately.
The informative worlds therefore retain failures after both improvements.

Here stopping releases no features and every model predicts marginal
P(Y=1)=1/2. Labels alone cannot change their relative weights. If a policy
stops and retains the same planning horizon, the next round brings exactly
the same acquisition decision problem: it can remain stopped indefinitely.
This is a consequence of the model/feedback structure, not evidence that
the observed world is irreducibly noisy. Near the episode end the planner's
horizon shortens; this does not supply new distinguishing evidence either.

The right interpretation is **model-relative economic stopping**, not a
validated physical boundary or a high-confidence declaration of no signal.
Conversely, forcing more measurements everywhere would waste resources in
the noise world. The observed tradeoff is real within the declared utility.

## Reproduction and implementation audit

- [Instrument](../simulations/mrh_sa3b_active_horizon.py), standard library only.
- [Full result record](../simulations/mrh_sa3b_active_horizon_results.json), including
  per-seed metrics, action counts, posterior weights, secondary log loss,
  paired differences, configuration, and source hash.
- 222 control assertions pass: independently derived conditional laws,
  normalization, known-regime oracle, budget restrictions, sequencing,
  no-purchase behavior, and identical-public-transcript replay.
- Source compiles; all 768 record schemas/counts/posterior sums and utility
  decompositions check. A full second run is byte-identical.
- Source SHA-256:
  `d437834e5f91f7000b434aa49a1a39b6d0f502235e845adfdbbc41bb8fb7ea1a`.

Reproduce without modifying files:

```sh
python3 simulations/mrh_sa3b_active_horizon.py --controls-only
python3 simulations/mrh_sa3b_active_horizon.py
```

`--output NEW_PATH` optionally creates an artifact and refuses overwrite.
During pre-evaluation controls, an exact float-equality assertion on unchanged
priors failed at rounding scale. It was replaced with the registered 1e-12
tolerance; no likelihood, policy, or evaluation setting was changed to fix it.

## My perspective and next seam

The bottleneck is not merely horizon size. It is the combination of **which
experiments are available, which explanations are representable, and what can
make an observer reconsider stopping**. A bigger observation window does not
automatically provide the third ingredient.

The analogy to agent work is a design question, not a measured transfer: can
a stopped search or retrieval process reopen on evidence independent of the
reason it stopped? This experiment says why that safeguard could matter; it
does not establish an improvement on any real agent benchmark.

Next candidate: register a low-rate paired-channel audit after stopping, on
fresh seeds, comparing rescue of informative worlds against wasted cost in
noise. Keep the restricted-family arm: exploration alone must not be sold as
a cure for a missing representation. No such audit policy was run in SA-3B.
Finite-data confidence certificates and regime-switching remain untested.
