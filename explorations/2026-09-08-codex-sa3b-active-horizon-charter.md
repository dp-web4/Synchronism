# SA-3B registration: buying evidence, learning its value, and stopping

Codex, 2026-09-08. Status: registered before implementation and evaluation.
Authorization: dp asked to pursue the active-horizon direction. This is a
bounded synthetic decision experiment, not a physics prediction or a SAGE run.
SA-3A's exact identifiability controls are complete; its proposed finite-data
confidence-interval calibration remains a separate, unexecuted candidate.

## Question and disclosed assumptions

Can an observer learn whether to retrieve a historical record, buy a fresh
measurement, buy both, or stop? Separate two obstacles: insufficient planning
about the value of learning, and absence of the useful relation from the
observer's hypothesis family. The observer does NOT discover arbitrary laws:
we supply the candidate likelihoods, noise rates, costs, and uniform prior.

The historical record here is a previously stored bit associated with the
current target, not a simulated dynamical trajectory. Retrieval and sensing
differ in provenance, price, and measurement noise; this is an acquisition
analogue of memory versus extra state, not a Mori–Zwanzig derivation.

Each world persists for 64 independent rounds. H and S are independent fair
bits each round. In the memory world Y = H xor N; in the sensor world
Y = S xor N; in the noise world Y is independent fair; in the parity world
Y = H xor S xor N. Task noise N has probability 0.1. Retrieval reveals H
exactly; sensing reveals S xor M with independent measurement noise 0.1.
No useful bit is initially visible. The world identity never reaches policy
code. Y is revealed only after acquisition and prediction, as free feedback.
This feedback assumption is essential, not a claim about deployment costs.

Actions are stop, history, sensor, both. Prices in prediction-error utility
units are 0, 0.04, 0.08, 0.12. Resource charges are 0, 1, 2, 3 credits;
the per-round cap is 3 credits. There is no bankable cross-round budget.
The utility is 0–1 prediction error plus acquisition price. Labels tied at
probability 0.5 predict zero. These prices are declared design choices, not
measured constants. Acquisitions within a round are bundled: no second query
conditional on the first query's value in this rung.

## Frozen policies and evaluation

Two hypothesis families: restricted {memory, sensor, noise}; expanded adds
parity. Both start uniform and update by Bayes' rule using ONLY acquired
values and subsequent target feedback. No oracle regime label, unqueried
values, or future targets are supplied. Expanded-family testing is registered
now, not introduced after seeing restricted-family failure.

For each family, evaluate (1) myopic Bayes action selection and (2) exact
two-round Bayes lookahead, receding at each round and truncated to one at the
last round. The latter integrates over acquired observations AND later label
feedback before valuing the next round. It is not a 64-round optimal policy.
Action ties within 1e-12 favor stop, history, sensor, both in that order.

Four static acquisition baselines always stop/retrieve/sense/buy both, with
the expanded family's same Bayesian predictor and feedback updates. Report
all eight policies on all four worlds, 24 seeds (0 through 23), 64 rounds.
Each seed/world provides an action-independent stream of latent bits and
noises shared across policies. No seed search or post-evaluation tuning.

Primary metric: mean error plus acquisition price per round, with per-seed
values retained. Report error, price, action frequencies, last-16-round action
frequencies and utility, and posterior weights. Paired differences compare
two-step versus myopic within each family and expanded versus restricted
two-step; report mean and standard error across 24 independent seeds, not
24×64 independent observations. No binary significance claim is registered.
Secondary prediction metric: log loss in nats, without optimizing for it.

Reference oracle knows the regime but sees only purchased measurements. Its
analytic per-round expected utilities are memory 0.14 (history), sensor 0.26
(sensor), noise 0.50 (stop), parity 0.30 (both). Realized finite samples can
beat an expectation; that is not an oracle violation.

## Controls and failure criteria

- All action/model event likelihoods normalize and agree with independent
  closed-form predictive probabilities for history, noisy sensor, and parity.
- Known-regime policy actions/values match the analytic oracle. A one-credit
  cap makes sensing/both illegal; sufficiently high prices force stopping.
- Public-history replay gives identical decisions and predictions; stopping
  yields no acquired values and labels alone leave the uniform prior intact.
- Predictions precede feedback; double feedback or feedback before prediction
  must be rejected. Invalid actions and over-budget requests must fail.
- Exact planning probe: restricted prior (0.1, 0.1, 0.8), where myopic
  retrieval is break-even. Report whether two-step planning purchases useful
  learning. Absence of a difference is a result, not permission to tune costs.
- Parity is a deliberate misspecification test for the restricted family.
  Report its failures without claiming that low loss on in-family worlds
  certifies closure. Expanded-family success would demonstrate supplied
  representational capacity, not spontaneous invention of parity.
- Run deterministic repeat checks, compile validation, and preserve source
  hash plus raw seed-level results. Failed implementation controls block
  interpretation; disappointing policy performance does not.

## Prior art, scope, and stop

This is standard finite Bayesian sequential decision-making applied to the
MRH question. Exploration versus immediate reward is established prior art;
see [Russo and Van Roy, Learning to Optimize via Information-Directed Sampling](https://arxiv.org/abs/1403.5556).
We do not implement IDS or claim its guarantees. SA-3A and the repository's
[Markov Phase 3](2026-08-17-markov-phase3-causal-vs-relevant-horizon.md) supply
the local predictive-horizon context.

Unquestioned assumption to expose: the action menu and hypothesis family
already contain the right experiment and explanation. Practice to avoid:
using evaluator truth as an observer's stopping certificate. Operator
pushback anticipated: hand-built worlds and prices make this an instrument
test, not evidence of fleet capability. No substrate axioms, conservation
laws, dimensional reduction, or new physics predictions are tested here.

Stop after this registered run and interpretation. Do not silently add
regime switches, confidence certificates, new families, longer planners,
or a SAGE integration. A new rung requires an explicit recorded design.
