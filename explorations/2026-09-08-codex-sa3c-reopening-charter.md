# SA-3C registration: periodic audits after stopping

Codex, 2026-09-08. Prospective extension of completed SA-3B, authorized within
dp's request to pursue this arc. Register and push before implementation/run.

## Question

Can a small, externally scheduled evidence purchase rescue an observer from
SA-3B's self-perpetuating stopping state? Does the benefit depend on already
representing the useful relationship? This is an acquisition safeguard, not
a confidence certificate, new physics prediction, or agent-benchmark claim.

## Frozen design

Keep SA-3B's four worlds, two hypothesis families, priors, likelihoods, noise,
prediction loss, acquisition prices, per-round credit cap, and free label
feedback. Use its exact two-round planner, truncated at the last round.
Compare four policies: restricted and expanded family, each with and without
the following audit override.

At rounds 8, 16, 24, 32, 40, 48, 56, and 64, if the base planner recommends
stop, purchase BOTH history and the noisy sensor instead. Otherwise take
the base action. Predict and update through the same Bayesian model. The
override sees only the round number and recommended action, not truth,
unbought observations, prediction error, or world identity. Prices and
three-credit cap still apply. Maximum forced-audit expenditure is 8×0.12
per run; total policy expenditure can change further through learned actions.
The final-round audit is retained even though it has no subsequent learning
benefit within the run. It may improve that prediction or simply cost money.

Evaluate all four policies × four worlds × 24 NEW seeds (24 through 47),
64 rounds each, with paired, action-independent streams. No interval search,
new likelihoods, different seeds, or cost changes after evaluation.

## Metrics, controls, and interpretation

Primary: mean error + acquisition price. Retain per-seed error, cost, log
loss, action counts, audit counts, tail metrics, and final posterior. Compare
audit minus no-audit within each family using paired seed differences and
their standard errors. No significance threshold or universal-win claim.

Distinguish an actual reopening from a forced measurement: after round 64's
feedback, record the model's UNFORCED two-step recommendation for a hypothetical
fresh two-round continuation. This diagnostic is not an extra executed round
or additional scored data. For each seed, report whether the baseline would
stop and the audited policy would autonomously select the world's known
optimal channel (evaluation-only oracle label). Also report worsening cases.

Controls: reuse all SA-3B tests; no override off schedule or when a non-stop
action is preferred; exact scheduled audit counts for a known-noise observer;
known-noise expected extra cost 0.015/round and no predictive benefit; replay
consistency; unaudited implementation matches SA-3B baseline on an old seed.
The latter is a regression control only, excluded from new-seed estimates.
Compile and repeat full output byte-for-byte; record both source hashes.

A rescue is not proof that arbitrary missing variables are discoverable.
In parity the restricted family still cannot express the rule, even after
buying both inputs. The expanded observer was GIVEN that candidate rule.
If auditing fails to help, or harms more than it helps, publish that result
without tuning the schedule. If it helps informative worlds but wastes cost
in noise, that is the declared tradeoff, not a contradiction.

Unquestioned assumption: an external audit budget is available. Potential
misapplied practice: counting forced purchases as voluntary recovery. The
operator's fair objection: the schedule is arbitrary; this rung tests a
mechanism, not its optimum. No foundational substrate axioms are tested.

Stop after this rung and interpretation. Do not add representation expansion,
regime switching, learned schedules, or integrations without a new design.
