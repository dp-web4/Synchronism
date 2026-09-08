# SA-3A result: the price depends on the loss, and can be unidentifiable

**Author:** Codex (Astra)  
**Date:** 2026-09-08  
**Status:** first registered rung complete; six exact controls pass  
**Registration:** [test cards](2026-09-08-codex-sa3a-loss-and-identifiability-charter.md), committed and pushed as `5e82a9f4` before implementation/execution  
**Instrument:** [mrh_sa3a_loss_identifiability.py](../simulations/mrh_sa3a_loss_identifiability.py)  
**Machine record:** [exact-controls.json](../simulations/mrh_sa3a_results/2026-09-08-exact-controls.json)

## Outcome

**61 acceptance checks passed across the six registered controls and instrument
validation.** The second execution was byte-identical. Compilation passed.
No acceptance criteria changed and no failed run preceded these results.

This establishes the distinctions in the constructed finite worlds. It is not
a new physical law, a fitted validity estimator, or a confirmation of the
Synchronism substrate. The expected values were known analytically and explicitly
registered as controls, not prospective discoveries.

## Measured population quantities

All information is in nats; square risks use the target's numerical units.
Numbers below are rounded for readability; the JSON retains numerical precision.

| Control | Excluded-information CMI | Square risk: restricted → extended | What the instrument must distinguish |
|---|---:|---:|---|
| T1: closed but noisy | 0 | 0.09 → 0.09 | No exclusion cost does not mean no prediction error. |
| T2: excluded variable changes magnitude, not mean | 0.693147 | 2.5 → 2.5 | Distribution information can improve log loss without improving mean prediction. |
| T3: excluded variable predicts noisy binary target | 0.368064 | 0.25 → 0.09 | A positive control with benefits under both losses. |
| T3: same target multiplied by 10 | 0.368064 | 25 → 9 | Information is invariant while squared-error cost changes by 100×. |
| T4: previous state added to present observation | 0.368064 | 0.25 → 0.09 | Retained history can recover predictive state. |
| T5: direct world / hidden-innovation world | 0 / 0.325083 | 0.09 → 0.09 / 0.09 → 0 | Identical visible distributions can conceal different exclusion costs. |
| T6: predictor offset by 10 | not a CMI case | predictor MSE = 100 | Centered residual variance is zero and the E-style index is 1 despite large error. |

For every joint-table case, independently evaluated optimal log-risk reduction
equals CMI within the registered absolute tolerance of 1e-12. In T2 the actual
MSE ratio is **1**, while the Gaussian formula misapplied to it would return
**4**. That is a direct counterexample to extending the Gaussian ratio identity
to arbitrary non-Gaussian targets, not a failure of the identity within its scope.

## The sharpest result: a channel cannot identify its own omitted information

T5's two worlds have exactly the same visible distribution:

| X | Y | Probability in either world |
|---:|---:|---:|
| 0 | 0 | 0.45 |
| 0 | 1 | 0.05 |
| 1 | 0 | 0.05 |
| 1 | 1 | 0.45 |

In the direct world, a hidden Z is irrelevant and independent of the innovation
that produces Y. In the other world, Z **is** that innovation. Seeing Z buys
nothing in one world and removes all uncertainty in the other. Seeing only X,Y
cannot distinguish them, even with unlimited observations of this channel.

The observer function accepts only P(X,Y), not Z or a world identifier. It
returns the same answer for both:

```text
status: not_identifiable
compatible exclusion-log-loss interval: [0, 0.3250829733914482] nats
```

This is the sharp population interval `[0,H(Y|X)]` over unrestricted finite
hidden extensions. The two registered worlds attain its endpoints. The output
is derived from that model-class assumption and the visible distribution; it
is not an algorithm discovering which hidden world is present. Restricting the
admissible hidden mechanisms or acquiring a new channel can change the answer.

**Consequence for the MRH program:** observable residuals can assess an observable
predictor, but cannot generally reveal how much an unmeasured variable would
improve it. An honest validity record needs an identifiability field, not just
a numerical residual or information score.

## The memory control

T4 uses `X_(t+1)=X_(t-1) XOR N_t` with 10% independent flip probability. The
four-state pair transition matrix was constructed explicitly and its uniform
stationary distribution checked. The present X_t alone carries no prediction
advantage for the next state; the previous X_(t-1) does. Retaining that history
reduces MSE from 0.25 to 0.09.

This is not evidence that every hidden state is recoverable from history. It
demonstrates why a present-only MRH is a different restriction from a
history-aware MRH. It builds on the earlier target/horizon-relative
[Markov Phase 3](2026-08-17-markov-phase3-causal-vs-relevant-horizon.md), rather
than replacing its definition.

## Scope, limitations, and handoff

- These are exact finite-support probability tables evaluated in floating-point
  arithmetic. There is no sampling uncertainty, fitting, or empirical p-value.
- The five malformed-table checks and deterministic-target observer check test
  the instrument. They are not extra scientific results or independent findings.
- No finite-data CMI estimator or uncertainty-calibration method was implemented.
  The observable entropy is assumed known; a plug-in empirical entropy would not
  automatically give a valid finite-sample upper confidence bound.
- No SPARC files, Kimi instruments, or prediction-bucket statuses were changed.
  The [overview](2026-09-08-codex-mrh-validity-overview.md) separately records
  review findings for Kimi's SA-2 lane; this suite does not recompute that lane.
- The refusal is a result under the specified observation channel and model
  class, not a claim that excluded physics is intrinsically unknowable.

Per the registered stopping rule, this rung stops here. The next candidate is
a **separately registered finite-data/instrument-noise calibration**: retain the
same counterexample worlds, test observable prediction-loss bounds for coverage,
and require that added samples never manufacture identification where the
population distributions are identical. That follow-up is proposed, not run.

## Reproduction and provenance

```bash
python3 simulations/mrh_sa3a_loss_identifiability.py
```

Stdlib only; stdout JSON; no file writes or external datasets. Implementation
SHA-256, also recorded in the machine output:

`586891ed25cae5db5d274d2fcb37584b80dd01e97afb0b5fcffd6579738548e2`

Self-check before committing: known ideal distributions are the untested
real-world assumption; exact enumeration is appropriate for these finite
controls but not a claim about real-data uncertainty; preserving the MRH
program while narrowing its claims respects the operator's frame; no physical
dimension, conservation law, or substrate axiom is modeled by these tables.
