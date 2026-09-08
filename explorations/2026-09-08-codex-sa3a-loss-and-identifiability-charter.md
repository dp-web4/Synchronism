# SA-3A claim: loss-aware predictive validity and identifiability

**Owner:** Codex (Astra)  
**Date:** 2026-09-08  
**Status at registration:** claimed; computation suite not yet implemented or run  
**Parent:** [MRH-validity SA-3](2026-09-08-kimi-mrh-validity-subarc-map.md)  
**Authorization:** dp requested an overview and choice of a complementary arc

**Execution update:** the [first-rung result](2026-09-08-codex-sa3a-exact-controls-results.md)
passes all six controls (61 checks). The registration below is preserved;
implementation and execution followed its publication in commit `5e82a9f4`.

## Boundary and prior work

Kimi has executed SA-2 rung 2 and an SA-1b first cut in
[results 02](2026-09-08-kimi-mrh-validity-results-02-sparc-closure.md), and is
continuing with SPARC rungs 3 and 4. This claim does not take those tasks or
change Kimi's instruments, charter, or kill criteria. SA-3's broader Fano/bound
work remains open; SA-3A is the bounded calibration and identifiability component.

The earlier [Markov Phase 3](2026-08-17-markov-phase3-causal-vs-relevant-horizon.md)
already defines predictive MRH by target, horizon and CMI tolerance. This work
adds explicit loss comparisons and an observation-only identifiability test.
Entropy/Bayes-risk identities, finite Markov chains, and observational
equivalence are established mathematics; no theorem or physics novelty is claimed.

## Question and hypothesis

Can one honest validity record distinguish:

1. no predictive benefit from an excluded variable;
2. benefit depending on the prediction loss;
3. benefit recovered by retaining history; and
4. benefit that is not identifiable from the available observations?

**Hypothesis:** a loss-aware record with an explicit observation channel can
represent these cases without labeling all residual reduction as closure.
An observation-only estimate cannot identify excluded-information value across
a model class containing observationally equivalent worlds with different values.

## Procedure fixed before implementation

Implement one stdlib-only, deterministic exact-enumeration script:
`simulations/mrh_sa3a_loss_identifiability.py`. Use natural-log units (nats).
Compute CMI from its conditional probability-ratio definition, and expected
optimal log loss directly from predictive probabilities. Compute square-loss
Bayes risks independently from conditional means. Validate table probability
mass, support, and nonnegative probabilities before scoring.

All results are population quantities on finite support, not estimates from
samples. Acceptance absolute tolerance is **1e-12**. No parameter search, model
training, external datasets, environment source access, or substrate simulation.
Known analytic expectations below are controls, not prospectively unknown discoveries.

### T1: closed but noisy channel (negative control)

X and Z are independent fair bits. Y=X XOR N, N~Bernoulli(0.1), independent.
Test CMI=0 and identical full/restricted risks, while log and square losses
remain positive. This must not be reported as zero prediction error.

### T2: non-Gaussian distribution information with no mean improvement

X is constant; Z is a fair bit. Conditional on Z=0, Y is equally likely -1,+1;
conditional on Z=1, Y is equally likely -2,+2.
Test CMI=ln(2), equal square risks **2.5**, and log-risk reduction ln(2).
This falsifies applying the Gaussian MSE-ratio identity outside its assumptions:
the MSE ratio is 1, not exp(2*CMI)=4.

### T3: excluded mean information and target scaling (positive control)

X is constant; Z is fair; Y=Z XOR N, N~Bernoulli(0.1).
Test square risk **0.25 -> 0.09**, and log-risk reduction
`ln(2)-h(0.1)`. Repeat with numeric target values multiplied by 10:
CMI unchanged, square risks and their difference multiplied by 100.

### T4: history recovers an omitted state

Use the stationary binary process `X_(t+1)=X_(t-1) XOR N_t`, with independent
Bernoulli(0.1) innovations. Verify the uniform distribution of the four pair
states is invariant under the explicit transition matrix before measuring.
Present-only observation is X_t; extended observation includes X_(t-1).
Expected square risk **0.25 -> 0.09**; expected log-risk gain
`ln(2)-h(0.1)`. This is a stochastic-process control, not a physical model.

### T5: observation-only identifiability fails (adversarial control)

Both worlds have fair X and Z~Bernoulli(0.1), independent of X.

- **Direct world:** Y=X XOR N, with N~Bernoulli(0.1) independent of X,Z.
- **Hidden-innovation world:** Y=X XOR Z.

Verify identical observable P(X,Y) in both worlds, but CMI **0** versus
**h(0.1)** respectively. The observer sees only X,Y, never Z or world labels.
An observation-only summary must therefore not issue a unique exclusion-loss
value: report **not identifiable**, with compatible witnesses and the tight
range **[0,h(0.1)]** for a class containing these endpoint worlds. That range
also follows from `0 <= I(Y;Z|X) <= H(Y|X)` for this finite target.
Do not feed hidden truth into the observer's summary and call it discovery.

### T6: centered variance misses bias (metric control)

For X uniform on {-1,+1}, target Y=X, use predictor X+10.
Test centered residual variance=0, an E-style index=1 against Var(Y), and
MSE=100. This tests the interpretation of a variance score, not the SPARC fit.

## Falsifier / failure handling

- Any incorrect analytic value, stationarity check, probability normalization,
  or log-loss/CMI equality fails the suite and must be reported before repair.
- If T5's observable distributions differ, it is not an identifiability
  counterexample; reject that construction rather than relax the criterion.
- Passing controls establishes the intended distinctions in these examples,
  not a general learned validity estimator, a physical bound, or a galaxy claim.
- No fitted thresholds and no empirical significance claims. Tests of formulas
  are not independent confirmations of Synchronism's physical hypotheses.

## Deliverables and stopping rule

Commit and push this claim before executing the new suite. Then deliver the
script, a machine-readable report, and a dated result note with the exact command
and implementation revision. Stop after these controls. A finite-data estimator,
sensor-selection policy, or application to SPARC requires a separately stated
next rung; do not silently broaden this test into one.

Pre-commit self-check: ideal known distributions are an explicit assumption;
no 3-D physics is approximated by these finite probability tables; no
conservation axiom is tested; a negative inference result does not dismiss the
MRH program. The useful result is a more precise boundary on what it can claim.
