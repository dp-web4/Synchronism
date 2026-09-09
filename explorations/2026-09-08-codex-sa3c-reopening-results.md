# SA-3C results: reopening works when the observer can use what it finds

Codex, 2026-09-08. Registered synthetic experiment, not new physics or a
demonstrated improvement on a real agent benchmark.

**Periodic audits rescued premature stopping, but did not repair a missing
representation. They also caused some unnecessary reopening in pure noise.**

[Registration](2026-09-08-codex-sa3c-reopening-charter.md) was committed and
pushed at `77c79f3e` before implementation/evaluation. The design was not
tuned after observing results. It extends [SA-3B](2026-09-08-codex-sa3b-active-horizon-results.md)
on fresh seeds 24–47, not a rerun of the original seed cohort.

## Intervention and what counts as recovery

Keep the two-round Bayesian planner, noise, likelihood families, costs,
feedback, and three-credit per-round cap. On every eighth round, override a
recommendation to stop with a paid purchase of both channels. Otherwise leave
the planner alone. The audit sees only the clock and proposed action.

There are 24 runs per policy/world, each 64 rounds. Forced measurements do
NOT count as autonomous recovery. After the final feedback, ask what the
unforced planner would purchase for a hypothetical two-round continuation.
That diagnostic receives no additional observations and is not a scored round.

## Main results

Utility is prediction error plus measurement price per round; lower is better.
Differences are audit minus baseline, reported as mean ± one paired-seed
standard error, not confidence intervals or significance declarations.

| Family | World | Baseline | Audited | Paired difference |
|---|---|---:|---:|---:|
| Restricted | Memory | 0.1787 | 0.1660 | −0.0127 ± 0.0093 |
| Restricted | Sensor | 0.3355 | 0.2849 | −0.0506 ± 0.0192 |
| Restricted | Noise | 0.5012 | 0.5218 | +0.0207 ± 0.0061 |
| Restricted | Parity | 0.4974 | 0.5147 | +0.0173 ± 0.0073 |
| Expanded | Memory | 0.1494 | 0.1494 | 0.0000 ± 0.0000 |
| Expanded | Sensor | 0.3266 | 0.2978 | −0.0288 ± 0.0124 |
| Expanded | Noise | 0.5016 | 0.5191 | +0.0175 ± 0.0059 |
| Expanded | Parity | 0.4232 | 0.3622 | −0.0610 ± 0.0184 |

The parity result makes the distinction cleanest. With parity in the model,
12/24 baseline runs autonomously choose both channels at the final diagnostic;
auditing raises that to 24/24. Every one of the 12 stopped baseline cases is
rescued by this diagnostic. Eleven have lower cumulative utility, one has
higher utility: final recovery is not identical to having recouped its cost
within 64 rounds. Mean forced audits in the expanded parity arm: 1.5 per run.

Without parity in the model, auditing produces **zero** correct autonomous
recoveries in that world. The final recommendations are stop in 22/24 runs
and sensor in 2/24, never both. Extra measurements increase mean utility
from 0.4974 to 0.5147. The information is purchased, but the predictor lacks
the relation needed to exploit it.

## The costs and remaining failures are part of the result

- Expanded sensor: correct final recommendations increase from 18/24 to
  23/24; one run still does not recover. Restricted sensor improves 17/24
  to 24/24. Neither implies guaranteed finite-time recovery.
- Restricted memory improves 22/24 to 23/24. Expanded memory is already
  correct in all 24 new-seed baseline runs, so no audits occur and its
  audited output is unchanged. This does not erase SA-3B's failures on a
  different seed cohort.
- In pure noise every baseline ends with the correct stop recommendation.
  Auditing makes 4/24 expanded observers recommend both and 2/24 restricted
  observers recommend sensing at the final diagnostic. Those purchases are
  unnecessary under evaluator truth. More observations can produce finite-
  sample false leads; a periodic audit is not a false-alarm guarantee.
- For a known-noise observer, the exact control purchases eight audits,
  adds 0.015 utility per round, and gains no predictive benefit. Unknown-
  regime runs can pay more through changed subsequent decisions and errors.
- Expanded parity's audited utility 0.3622 still exceeds the regime-aware
  oracle expectation 0.3000. The whole episode includes learning and mistakes;
  a correct final recommendation is not oracle-level cumulative performance.

The schedule was fixed, not optimized. The relative frequency of informative
and noise worlds in a real workload would determine whether this insurance
is worth its cost; these four hand-built worlds do not estimate that frequency.

## What this adds to the MRH question

The original stopping problem was endogenous: no measurement meant no
distinguishing feedback, which preserved the belief that justified stopping.
A clock-triggered audit can break that loop without inspecting hidden truth.
But its success requires a useful available experiment AND a model capable
of interpreting its outcome. It can also create false leads.

That gives three separable engineering questions rather than one instruction
to "expand the horizon":

1. **Acquisition:** can the observer access a distinguishing experiment?
2. **Representation:** can its candidate explanations express the distinction?
3. **Reopening:** can it gather disconfirming evidence after it decides to stop?

SA-3A showed why current observations alone do not identify all excluded
information. SA-3B showed that choosing information has both immediate and
future learning value. SA-3C shows a limited way to reopen an observation
policy, including its costs and failure modes. These are applications of
existing statistical decision ideas, not discoveries of physical laws.

My transfer proposal is an interface requirement, not an implementation claim:
a stop record should carry its assumed model family, unavailable/unbought
channels, cost-based reason, and explicit reopening condition. "No purchase
is currently worth its price under this model" should never silently become
"there is no useful information outside this boundary."

For real agent work the relevant test would preserve the same separation:
change retrieval/probing while holding representation fixed, then change
representation while holding evidence access fixed, and measure the extra
cost. No private repository was read or modified for this experiment, and no
benchmark transfer was attempted.

## Artifacts and verification

- [Instrument](../simulations/mrh_sa3c_reopening_audit.py), reusing the frozen
  [SA-3B implementation](../simulations/mrh_sa3b_active_horizon.py).
- [Full record](../simulations/mrh_sa3c_reopening_audit_results.json): 384
  seed-policy runs / 24,576 decisions, per-seed metrics and audit timing,
  posterior weights, autonomous diagnostics, paired differences, source hashes.
- 222 inherited SA-3B assertions plus 1,028 SA-3C assertions pass. Many are
  exhaustive schedule checks, not independent scientific tests.
- Compile, both source hashes, all seed/count/audit-budget/posterior/utility
  checks pass. Full repeat output is byte-identical.
- SA-3C source SHA-256:
  `ef044bf08ca09fd55d4b077d7d3ab043256da5cf235631aae0d5264bef459550`.
- SA-3B dependency SHA-256:
  `d437834e5f91f7000b434aa49a1a39b6d0f502235e845adfdbbc41bb8fb7ea1a`.

```sh
python3 simulations/mrh_sa3c_reopening_audit.py --controls-only
python3 simulations/mrh_sa3c_reopening_audit.py
```

Both commands leave result files unchanged; `--output NEW_PATH` optionally
creates a new JSON artifact and refuses overwrite.

## Stop point

This rung is complete. The next nontrivial question is calibrated reopening:
how to balance missed signal against audit expense and false leads without
privileged knowledge of the world's class. That needs a new registered loss,
workload assumption, and evaluation—not another tuned schedule on these seeds.
No confidence certificate, arbitrary hypothesis discovery, changing-regime
adaptation, substrate result, or deployed-agent improvement has been established.
