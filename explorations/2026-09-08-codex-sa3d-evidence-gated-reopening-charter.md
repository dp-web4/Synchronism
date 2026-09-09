# SA-3D registration: evidence-gated reopening

Codex, 2026-09-08. dp requested continuation after SA-3C. Commit and push
this design before implementation or evaluation. No new physics claim.

## Question and frozen scope

SA-3C's periodic audits rescued informative worlds but also provoked renewed
querying in pure noise. Can an anytime-valid evidence gate limit that latter
event, and what recovery, delay, and utility does it sacrifice?

Reuse SA-3B's worlds, supplied model families, priors, prices, known noise
rates, free post-prediction labels, 64-round duration, and three-credit
per-round cap. Reuse two-round Bayesian planning and SA-3C's every-eighth-
round paired-channel audit. No schedule, threshold, or world tuning.

## Evidence process and precise guarantee

The null P0 is the specific existing noise world: fair independent H, S,
and target, with the specified sensor noise. For each supplied non-noise
candidate k, initialize L_k(0)=1. After each purchased observation O and
subsequent target Y, multiply by

`P_k(O,Y | selected action) / P_0(O,Y | selected action)`.

E_t is the fixed uniform average of these cumulative likelihood ratios.
All acquired evidence counts from the start; never reset or select a
favorable subset. Unpurchased values and evaluator identity are unavailable.
Actions depend only on prior public evidence, so under P0 each likelihood
ratio is a nonnegative mean-one martingale, as is their fixed mixture.
Ville's inequality gives `P0(sup_t E_t >= 20) <= 0.05`.

This is an application of established sequential testing, not a new theorem;
see [Ramdas et al., Game-theoretic statistics and safe anytime-valid inference](https://arxiv.org/abs/2210.01948).

Once the base planner first recommends stop, mark that fact permanently.
Thereafter, in the gated policy, allow an unforced non-stop recommendation
only if the CURRENT E is at least 20. Otherwise veto it to stop. Apply the
periodic audit to the resulting recommendation: a vetoed stop can still buy
both on rounds 8,16,...,64. All bought values enter the same Bayesian update
and likelihood-ratio update. The evidence threshold is evaluated BEFORE
buying the current measurement, never using the target being predicted.

Thus any unforced non-stop action after first stopping requires a preceding
threshold crossing. Under P0 the probability of EVER such an action is at
most 5%, even with adaptive acquisition. This does NOT bound initial
exploration, scheduled audits, total spend, selection of a wrong informative
channel, or false alarms under an incorrectly specified null. Failing to
cross is NOT evidence of closure or a certificate of irreducible noise.

## Frozen comparisons and evaluation

Main cohort: fresh seeds 48–71; four worlds × two families × three policies:
no audit (SA-3B), periodic audit (SA-3C), evidence-gated periodic audit.
64 rounds per run. Use paired action-independent streams across policies.
Restricted-family parity remains a missing-representation control.

Null-event cohort: fresh seeds 1024–1535 (512 seeds), noise world only,
restricted/expanded families × periodic/gated audits. No no-audit arm in
this additional cohort. Record empirical event counts and proportions; these
do not prove or invalidate a mathematical 5% bound by themselves. A random
empirical fraction can exceed its population bound. The proof's assumptions
and instrument checks, not a favorable Monte Carlo count, justify validity.

Primary metrics: prediction error + acquisition price, event ANY unforced
query after first stopping, correct final autonomous recommendation, and
number/timing of audits and gate vetoes. Retain per-seed errors, prices,
log loss, action counts, tail metrics, posteriors, final and maximum E,
first stop/crossing/reopening rounds. The final recommendation is the same
unexecuted two-round-continuation diagnostic as SA-3C, with the evidence gate
applied and its underlying ungated recommendation also recorded.

Report paired utility differences and standard errors across seeds: gated
minus periodic, periodic minus no-audit, gated minus no-audit within family
in the main cohort; gated minus periodic in the null cohort. Report correct
final channel counts and false reopening counts separately. No threshold
search, significance-based winner selection, or pooled real-workload claim.

## Controls, risks, and stop

- Reuse SA-3B/SA-3C controls and regression-check unchanged baseline behavior.
- Verify likelihood-ratio normalization/conditional martingale identities
  with independent exact rational probabilities for every action and family.
- Check no observations means no evidence change; no double feedback; the
  current target cannot authorize its own acquisition.
- Test threshold boundary, post-stop-only gating, audit eligibility after
  veto, and no evidence resets across repeated stopping/reopening.
- Verify every allowed post-stop unforced query was preceded by E>=20;
  count scheduled audits separately and enforce the original credit cap.
- Compile, check schemas/source hashes, repeat full output byte-for-byte.

The assumption exposed is a fully specified stationary null and available
candidate alternatives. The likely cost is delayed or missed recovery from
requiring strong evidence. Audits themselves cannot be made free by a test.
Operator objection: a 5% gate is a chosen risk budget, not automatically the
right business or scientific utility. No substrate axiom is tested.

Stop after this rung's interpretation. Unknown/composite noise, nonstationary
worlds, representation discovery, optimized risk budgets, and deployment
require another design; do not quietly broaden this guarantee to them.
