# SA-3F results: spend on reachable evidence, not an impossible reopening promise

Codex, 2026-09-08. [Registration](2026-09-08-codex-sa3f-feasible-bursts-charter.md)
was committed and pushed at `cf65fd7a` before execution. All declared policies,
caps, thresholds, seeds, and costs were retained. This is a synthetic observer
implementation, not a physical discovery or measured agent-benchmark transfer.

**The feasibility-aware observer is implemented.** It distinguishes an
unreachable evidence requirement from a reachable one whose detection power
is not certified, and names the supplied models for which its power lower
bound is sufficient. It never calls insufficient budget “no signal.”

The bounded burst recovers some cases that periodic auditing cannot. Pruning
unreachable audits preserves the full burst's reopening outcomes while
reducing its noise-world audit count by roughly half. Neither finding makes
the method uniformly cheaper or solves the missing-hypothesis problem.

## Policies, budgets, and scope

All policies share the previous two-step Bayesian predictor and the lifetime
E>=20 gate on unforced queries after first stopping. Measurements are paid,
labels arrive after prediction, and each round has the same three-credit cap.
There are 64 rounds and no additional observations outside the scored run.

- **periodic24:** 24-audit allowance, but only every eighth round is eligible;
  the calendar permits at most eight audits.
- **burst8:** eligible stopped rounds are audited immediately, up to eight.
- **burst24:** the same immediate rule, up to 24.
- **feasible24:** burst24, except refuse audits once the exact remaining-
  budget envelope cannot reach E=20.

An audit buys both channels and costs 0.12. The audit allowance is NOT a cap
on ordinary purchases; total acquisition price includes both. Caps never
replenish. The 80% power target is a reporting threshold, not another action
gate. Report/planner computation itself has no simulated fee in this model.

## Main cohort: timing alone was insufficient here

24 new paired seeds per world/family/policy. Mean error + total acquisition
price per round, lower better:

| Family | World | Periodic24 | Burst8 | Burst24 | Feasible24 |
|---|---|---:|---:|---:|---:|
| Restricted | Memory | 0.1618 | 0.1606 | 0.1453 | 0.1453 |
| Expanded | Memory | 0.1605 | 0.1587 | 0.1438 | 0.1438 |
| Restricted | Sensor | 0.3174 | 0.3186 | 0.2992 | 0.2981 |
| Expanded | Sensor | 0.3359 | 0.3392 | 0.3116 | 0.3115 |
| Restricted | Noise | 0.5429 | 0.5444 | 0.5816 | 0.5542 |
| Expanded | Noise | 0.5450 | 0.5399 | 0.5680 | 0.5447 |
| Restricted | Parity | 0.5098 | 0.5348 | 0.5608 | 0.5384 |
| Expanded | Parity | 0.4002 | 0.4102 | 0.3790 | 0.3816 |

With the expanded family, correct final autonomous choices change as follows:

| World | Periodic24 / Burst8 | Burst24 / Feasible24 |
|---|---:|---:|
| Memory | 23/24 | 24/24 |
| Sensor | 18/24 | 21/24 |
| Parity | 13/24 | 16/24 |

No periodic24 or burst8 main-cohort run reopens autonomously after stopping.
Both 24-audit burst variants reopen in five expanded sensor runs and five
expanded parity runs. Some subsequently lose sufficient current evidence and
stop again, so reopening once is not the same as a correct final choice.

For expanded parity, burst24 minus periodic24 utility is −0.021198 ± 0.016134;
for expanded sensor it is −0.024323 ± 0.015438. These are paired means ± one
seed-level standard error, not confidence intervals or declared significance.
The first-stop reports classify every stopped main-cohort periodic24/burst8
case as unreachable within its remaining audits. Faster use of the same
eight-audit allowance therefore cannot repair those cases; the larger
allowance, not merely faster timing, matters in this cohort.

The restricted family still has zero correct final parity choices under all
four policies. Feasible24 spends less than burst24 there, but remains worse
in utility than periodic24. It cannot invent the missing relationship.

## Larger noise cohort: pruning saves audits, not all overhead

128 additional paired noise seeds per family/policy:

| Family / policy | Mean audits | Mean total utility | Ever post-stop unforced query |
|---|---:|---:|---:|
| Restricted / periodic24 | 6.852 | 0.52893 | 0/128 |
| Restricted / burst8 | 7.898 | 0.52991 | 0/128 |
| Restricted / burst24 | 23.547 | 0.55925 | 0/128 |
| Restricted / feasible24 | 10.727 | 0.53485 | 0/128 |
| Expanded / periodic24 | 6.734 | 0.53065 | 0/128 |
| Expanded / burst8 | 7.891 | 0.53233 | 0/128 |
| Expanded / burst24 | 23.523 | 0.56125 | 1/128 |
| Expanded / feasible24 | 11.383 | 0.53910 | 1/128 |

Feasible24 minus burst24 utility is −0.024404 ± 0.003862 (restricted) and
−0.022153 ± 0.003667 (expanded). Audit-count savings are about 54% and 52%
respectively; those percentages do not describe savings in total utility.
Feasible24 still spends more and has worse mean noise utility than periodic24.

One expanded burst trajectory reopens falsely, and pruning preserves that
event too. The inherited specified-null 5% per-record bound is not a promise
of zero errors; initial exploration and paid audits remain outside its event.
These empirical counts do not prove the bound, whose likelihood-ratio
assumptions are unchanged. They also are not a joint fleet-wide error bound.

## Why pruning preserves reopening—and what it can still sacrifice

For current component likelihood ratios L_k and m remaining paired audits,
the exact maximum evidence is

`U(m) = (sum_k L_k * r_k^m) / K`,

with r_memory=1.8 and r_sensor=r_parity=1.64. Repeated H=sensor=Y=0 attains
every maximum simultaneously. While the gate is closed, non-audit stops
provide no distinguishing evidence. Therefore U(m)<20 means that no allowed
future audit outcomes can reopen it before the deadline.

Feasible24 and burst24 have identical histories up to the first refusal. At
that point burst24 cannot reach the threshold using its remaining audits
either. Thus pruning preserves gate-crossing/reopening outcomes and the final
gated recommendation, not just on average but path by path in this setting.
All 448 paired burst24/feasible24 trajectories across both cohorts agree on
first stopping, first crossing, first reopening, and final unforced action;
feasible24 never uses more audits.

This is an analytic consequence of the registered envelope, additionally
checked on the results. It does NOT imply dominance in prediction utility:
even an audit incapable of earning a reopening certificate can improve its
own target prediction. Expanded parity illustrates that: feasible24's mean
utility is 0.002526 ± 0.005464 worse than burst24 despite identical final
choices and fewer audits. Stochastic outcomes and current-prediction value
both remain visible in the scored comparison.

## The observer's report is conditional, not an assertion that a signal exists

The implemented statuses are `evidence_met`, `unreachable_with_remaining_budget`,
`reachable_power_not_certified`, and `power_supported_for_named_models`.
The report carries exact current/maximal evidence, remaining time and audits,
component ratios, per-model power lower bounds, and the explicit statement
that unrepresented alternatives are unknown.

For K alternatives, one component reaching K×20 is sufficient for the mixture
to reach 20. Its exact binary first-passage probability is therefore a valid
LOWER BOUND on mixture power under that named model. Failure of this sufficient
condition does not prove that mixture power is low. Reaching 80% does not
mean an 80% probability that the model is true.

An illustrative record is the first main-cohort seed, 72, in the restricted
observer's parity world. At round 14 its report includes:

```json
{
  "status": "power_supported_for_named_models",
  "supported_models": ["memory"],
  "power_lower_bounds": {"memory": 0.8531, "sensor": 0.5946},
  "unrepresented_alternatives": "unknown",
  "absence_of_signal_certified": false
}
```

Display probabilities are rounded here; the artifact retains exact fractions.
The report is correct about what it could test IF memory were the true
mechanism. It cannot promise to detect the actual parity relation, which its
family omits. This is why a feasibility report must name its alternatives
rather than summarize itself as “enough evidence available.”

## Deliverable and limits

The [implementation](../simulations/mrh_sa3f_feasible_bursts.py) provides the
observer, reusable report/envelope/component-power functions, regression
controls, and the frozen evaluation. [Full results](../simulations/mrh_sa3f_feasible_bursts_results.json)
retain per-seed costs, timing, posterior/evidence, and 3,527 observer reports.

There are 768 main-cohort and 1,024 noise-cohort runs, 114,688 decisions total.
876 new assertions plus 1,994 inherited checks pass. Independent binary/joint
path enumeration checks the power certificates and envelope; all source
hashes, budget invariants, report schemas, and paired pruning properties
check. A full second run is byte-identical. No controls failed or policies
were tuned after execution.

- Source SHA-256:
  `0d1f805d0f1f8af8d8a74b7ac8d74ae5f5f59c8a814b52b42c3b06815e9d1534`.
- Result SHA-256:
  `48fbe25870b69ca7c31c69d266111df7930cba00d1fb11770688572e16f9b1ff`.

```sh
python3 simulations/mrh_sa3f_feasible_bursts.py --controls-only
python3 simulations/mrh_sa3f_feasible_bursts.py
```

Output is stdout; `--output NEW_PATH` creates a new artifact and refuses to
overwrite. The finite model and its likelihoods are supplied, computation has
no simulated fee, and feedback is free. No learned schedule, new hypothesis
discovery, changing-world guarantee, or real-agent integration was attempted.

My assessment: this is a useful recovery-contract component, not a universal
query policy. It can refuse an impossible promise with a checkable reason,
retain the full burst's attainable reopenings, and expose what it does not
know. Choosing whether those reopenings are worth their total cost remains
a workload-dependent decision. This registered rung is complete.
