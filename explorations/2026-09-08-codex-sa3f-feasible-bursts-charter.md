# SA-3F registration: feasibility-aware bounded auditing

Codex, 2026-09-08. dp approved implementation of the feasibility-aware
observer and burst/periodic comparison. Commit and push before execution.

## Frozen inherited setting

Reuse SA-3B's four stationary worlds, restricted/expanded hypothesis families,
uniform priors, supplied likelihoods, known noise, free labels after prediction,
64 rounds, two-step Bayesian planner, prediction-error utility, and acquisition
prices. Reuse SA-3D's lifetime mixture evidence and threshold 20: after the base
planner first recommends stop, unforced queries require CURRENT evidence >=20.
Never reset evidence. Initial exploration and paid audits remain exempt from
that narrowly defined gate; all measurements count in total price and utility.

Keep a three-credit per-round cap. Add a lifetime AUDIT count allowance, not
a limit on unforced queries. The maximum possible total acquisition spend
remains 64×0.12; separate audit expenditure from total expenditure everywhere.

## Four registered policies

All obtain the base recommendation, apply the post-stop evidence gate, and
consider a paired-channel audit only when the resulting action is stop AND
current evidence is below 20:

1. **periodic24:** audit at rounds 8,16,...,64 with allowance 24. The calendar
   allows at most eight purchases despite the larger allowance.
2. **burst8:** audit immediately on each eligible round, allowance eight.
3. **burst24:** same immediate rule, allowance 24.
4. **feasible24:** same allowance/rule as burst24, but decline an audit if no
   allowed sequence of remaining audits can reach 20 before the run ends.

Audit allowances never replenish when the gate opens/closes. Each audit costs
0.12 and reveals both channels; no free simulator peeks, extra labels, or
additional rounds. After evidence is sufficient, the base planner again
chooses unforced purchases. The periodic/burst24 comparison has equal allowed
audit count but may have different realized spend; burst8 separates some of
the timing effect from the allowance increase. No universal fairness claim
based only on the nominal cap: report actual expenditure.

## Honest observer report, computed from its own state

Maintain exact rational component likelihood ratios L_k alongside Bayesian
prediction. With K non-noise candidates, E = sum L_k / K. After m further
paired audits the maximum possible evidence is

`U(m) = sum_k L_k * r_k^m / K`,

where r_memory=1.8 and r_sensor=r_parity=1.64. All maxima are jointly attainable
by repeatedly observing H=noisy_sensor=Y=0, so this is an exact reachability
envelope for these supplied channels, not merely a loose optimistic estimate.
For periodic24, m respects remaining scheduled opportunities; for burst
policies, m=min(remaining audit allowance, remaining rounds).

For each represented alternative k, compute the exact first-passage
probability of its COMPONENT reaching K×20 within m independent paired audits,
using its current L_k and its known agreement probability (0.9 or 0.82).
That is a sufficient condition for mixture crossing, hence a LOWER BOUND on
mixture detection probability under that named alternative. It is not the
exact mixture power, nor a statement that an alternative is true.

Report one of:

- `evidence_met`: current E>=20.
- `unreachable_with_remaining_budget`: E<20 and U(m)<20.
- `reachable_power_not_certified`: reachable, but no computed lower bound
  reaches the predeclared 80% target. This does NOT prove actual power <80%.
- `power_supported_for_named_models`: list every represented alternative
  whose bound is >=80%, with all individual bounds visible.

No report may say absent signal, complete model, or certified closure.
Unrepresented alternatives remain unknown. The 80% threshold affects reporting
ONLY; feasible24 does not add another power gate to its acquisition rule.

Save reports immediately before the first stopping decision's acquisition,
at the first audit refused as unreachable, and after the last feedback.
The last report has zero remaining rounds: record separately the hypothetical
unforced two-round recommendation as in SA-3C/D, without executing it.

## Evaluation and controls

Main cohort: seeds 72–95, all four worlds, both families, all four policies.
Additional noise cohort: seeds 2048–2175 (128 seeds), both families, all four
policies. Paired action-independent streams; no seed/allowance/threshold search.

Retain seed-level total utility, error, total price, audit price, action/audit
counts and times, tail utility, final posterior, first stop/cross/reopening,
ever post-stop unforced query, final diagnostic recommendation, and observer
reports. Compare paired utility means and seed standard errors for burst8
minus periodic24, burst24 minus periodic24, and feasible24 minus burst24.
Report correct final channels and null reopening events separately. No
significance-selected winner, deployment claim, or mixture-power overclaim.

Controls: inherited likelihood/planning/gate checks; rational tracker agreement;
budget and sequence enforcement; no current/future target in acquisition;
exact all-zero reachability envelope; binary exhaustive checks of component
first passage; small exhaustive joint-path checks that mixture power is at
least the component certificate; identical-state report replay; report labels
at threshold/budget boundaries; periodic regression on prior control seeds
where the additional E<20 audit condition does not change actions. Explicitly
verify the distinction if that condition does change an old trajectory.
Compile, validate records/hashes, and repeat byte-identically.

Assumption exposed: the supplied alternatives and calibrated measurements
are available. Likely failure modes: bursting spends heavily in noise;
feasibility pruning may save certificate cost while forfeiting useful current
predictions; a missing hypothesis remains missing. Preserve each outcome.
No substrate axioms or physics predictions are tested. Stop after this rung;
no learned audit schedule or real-agent integration without another design.
