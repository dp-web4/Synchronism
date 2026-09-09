# SA-3G results: safety, power, and reachability fail differently

Codex, 2026-09-08. [Protocol](2026-09-08-codex-sa3g-calibration-charter.md)
committed and pushed at `d274f576` before implementation/evaluation.

The clean result is a separation, not a blanket statement that "uncalibrated
models are unsafe":

- A wrong **alternative** can badly overstate recovery power while leaving
  the specified-null false-alarm bound intact.
- A wrong **null** can break an incorrectly extended false-alarm claim.
- **Reachability** remains exactly checkable for the chosen update rule,
  even when its probability model is wrong.
- A predeclared **null uncertainty guard** restores a broader safety
  guarantee, but does not restore alternative-power calibration for free.

The quantitative warning: at a 24-audit budget, a nominal **85.28%** detection
probability becomes **60.27%** when true agreement falls from 82% to 75%.
Under a separate null perturbation to 55% agreement, the same nominal rule
crosses **6.44%** of the time; the guarded rule lowers that to **3.33%**, while
also reducing detection at 75% true agreement to **50.14%**.

These are exact synthetic probabilities, not sampling estimates or real-agent
benchmark results. No confidence interval is needed for the finite calculation;
the uncertainty is whether the supplied model class describes an application.

## What was actually tested

This rung isolates SA-3E's binary audit stream rather than rerunning SA-3F's
planner. A paid paired measurement yields a parity prediction; the later target
gives an agreement bit Z. True IID agreement is p; the evidence rule uses an
assumed alternative q and null ceiling u:

`agreement: multiply by q/u; disagreement: multiply by (1-q)/(1-u)`.

Evidence starts at 1 or 1/10, threshold is 20, each audit costs 0.12, and
auditing ends at crossing or the cap. q is one of {0.65, 0.82, 0.90}; p is one
of {0.50, 0.52, 0.55, 0.60, 0.65, 0.75, 0.82, 0.90}. Nominal u=0.50;
guarded u=0.55. All 96 exact curves, caps 0–64, are retained: **6,240 records**.
The primary cap is 24; 8 and 64 are registered reference/stress horizons.

The guarded null intentionally includes weak positive association. It is
**not independent fair noise**. Its guarantee covers conditional agreement
probabilities <=0.55 on each audit, including history-dependent ones. The
power curves concern IID alternatives. The guard is supplied, not learned.

There is no full-planner utility comparison here: reported spend is audit
expenditure, not prediction error plus all acquisition costs. A crossing is an
evidence event, not necessarily a correct downstream action or a reopening
the planner would choose. The existing SA-3F observer and JSON interface are
unchanged; no guarded evidence was substituted into their old records.

## 1. A wrong alternative is not automatically an unsafe test

For any q>u, the mean one-step factor at p=u is exactly 1. At p<=u it is
at most 1. Thus the product is a nonnegative supermartingale under the
conditional null, with crossing probability bounded by starting evidence / 20.
This is an application of established
[anytime-valid sequential inference](https://arxiv.org/abs/2210.01948), not a
new statistical theorem.

In particular, under fair noise the nominal gate retains its bound for all
three assumed alternatives. Its numerator need not describe an actually true
alternative to be a normalized betting rule. All tested fair-noise curves
respect their bound. What fails when q is wrong is the power claim obtained
by pretending that actual future agreement p equals q.

At start=1 and cap=24:

| Assumed q | Rule | Power if p=q | Actual power at p=0.75 | Crossing at p=0.55 |
|---|---|---:|---:|---:|
| 0.65 | Nominal | 12.95% | 46.17% | 1.92% |
| 0.65 | Guarded | 0.46% | 5.07% | 0.02% |
| 0.82 | Nominal | 85.28% | 60.27% | 6.44% |
| 0.82 | Guarded | 79.37% | 50.14% | 3.33% |
| 0.90 | Nominal | 98.29% | 56.60% | 6.90% |
| 0.90 | Guarded | 94.95% | 40.52% | 3.10% |

The q=0.90 guarded case is a particularly clear boundary: its broader-null
safety is valid, yet its nominal 94.95% power would be a gross overstatement
at p=0.75, where actual power is 40.52%. A safety guard is not a power guard.
Conversely, q=0.65 understates power when actual agreement is stronger; errors
need not all be optimistic.

## 2. A wrong null can break the extended numerical promise

For q=0.82 and start=1:

| Actual background agreement | Nominal crossing by 24 | Guarded crossing by 24 | Nominal crossing by 64 | Guarded crossing by 64 |
|---|---:|---:|---:|---:|
| 0.50 | 2.82% | 1.29% | 3.40% | 1.42% |
| 0.52 | 3.97% | 1.91% | 5.04% | 2.19% |
| 0.55 | 6.44% | 3.33% | 8.88% | 4.12% |

Two percentage points of background agreement drift suffice to exceed 5%
within the longer horizon: the exact p=0.52 nominal crossing probability is
**0.05035316… by 64** (first exceeds 5% at cap 56). The excess is small, not
catastrophic. At p=0.55 the violation is already visible at the primary cap.

These do **not** refute SA-3D's fair-noise theorem: p>0.50 is outside that
theorem's null. They refute carrying its numerical guarantee into a broader
environment without changing the rule or declaring the lost guarantee.
All guarded curves inside their declared null respect the bound; the
supermartingale argument extends beyond the tested IID grid to conditional
p<=0.55. Outside that class the guarded rule also loses its justification.

With starting evidence 1/10, the corresponding conditional bound is 0.5%,
not 5%. At p=0.55 and cap=24, nominal crossing is 0.9980%, guarded 0.1591%.
Again, the nominal continuation bound cannot be carried to the broader null.
These conditional starts require evidence inherited under the same rule;
switching nulls requires recomputing a valid record, not relabeling a ratio.

## 3. A supplied alternative interval prices the power uncertainty

Crossing is coordinatewise increasing in the agreement bits. For a declared
IID alternative interval p in [0.75,0.90], its exact worst-case power is
therefore the p=0.75 curve. No fit or arbitrary "confidence discount" is needed.
The interval itself is still a supplied assumption, not a measured guarantee.

For q=0.82 and start=1:

| Requirement | Nominal minimum cap | Guarded minimum cap |
|---|---:|---:|
| Any possible crossing | 7 | 8 |
| 80% power if p=0.82 exactly | 22 | 28 |
| 80% power throughout p in [0.75,0.90] | 40 | 64 |

Thus the old 24-audit allowance supports the nominal point-model 80% claim,
but neither the guarded point-model claim nor the interval-wide claim.
**None of the six registered q/rule combinations certifies 80% across that
interval within 24 audits.** That is an honest budget insufficiency for these
rules, not a proof that every conceivable test must fail at 24.

At p=0.82, guarding raises expected audits under cap 24 from 13.76 to 16.08
(spend 1.6514 to 1.9300), while lowering crossing from 85.28% to 79.37%.
At p=0.75, expected audits rise from 17.45 to 19.56 while crossing falls from
60.27% to 50.14%. Under this stop-on-crossing policy, safety costs both time
and recovery. This is not a universally acceptable trade; the workload must
price the benefit of recovery and the tolerated weak association.

Prior adverse evidence also matters. At start=1/10, p=0.75, q=0.82, neither
rule reaches 80% even by 64 (nominal 78.63%, guarded 57.98%). No reset is
authorized to make that result look better.

## 4. Reachability survives probability misspecification—within its scope

For a fixed rule, the best-case envelope remains `start*(q/u)^m`. It depends
on allowed updates and remaining audits, not on their true probabilities.
All registered true p values give both outcomes positive probability, so the
all-agreement path remains possible. A wrong p changes its probability, not
its existence. This preserves the rationale for refusing an impossible
crossing promise, even when the nominal power calculation is unreliable.

Changing the update rule is different: guarding changes the envelope. For
q=0.82 and start=1/10, the minimum possible crossing moves from 11 to 14
audits. The artifact's envelope is the maximum under hypothetical continued
updates; actual acquisition stops once crossing occurs.

This does not extend SA-3F's pruning theorem to arbitrary real agents:
additional channels, nonzero evidence from "stop," or a changed likelihood
family can invalidate that theorem's pathwise premises.

## Verification and research judgment

[Source](../simulations/mrh_sa3g_calibration.py),
[full artifact](../simulations/mrh_sa3g_calibration_results.json), and
[additional boundary/artifact tests](../simulations/test_mrh_sa3g_calibration.py).

47,046 instrument assertions pass: exact mass conservation, costs, conditional
null bounds, stopped moments, monotonicity, exhaustive small binary paths,
history-dependent null adversaries, and SA-3E/SA-3F regressions. Four additional
test groups pass, exercising nontrivial near-threshold paths for every rule
and replaying all 6,240 saved records and 780 power-report entries. A full
second execution is byte-identical. No failed controls or post-result policy
tuning; added near-threshold software checks did not change the experiment.

- Source SHA-256: `49b273ad2edbf1e0582ae02a9db12c097d0a686a023f6a905a1a3f1ab92cdf6a`.
- Artifact SHA-256: `55676e558e675b75d1da951bc7007d4685aa2d50b993423184bf29513826a110`.

```sh
python3 -B simulations/mrh_sa3g_calibration.py --controls-only
python3 -B simulations/mrh_sa3g_calibration.py --output NEW_PATH.json
python3 -B -m unittest discover -s simulations -p test_mrh_sa3g_calibration.py -v
```

The output option refuses overwrite. No data fetching, new dependency,
learned calibration, planner change, or physics-ledger change.

My perspective: the bottleneck is no longer just buying enough evidence.
It is deciding **which claims the evidence-generating process actually
licenses**. The observer needs separate records for null coverage,
alternative-power coverage, and rule-level feasibility. A mathematically
consistent JSON report cannot supply any missing coverage assumption.

The next substantive question is how to obtain those uncertainty bounds
from a finite, paid calibration record without reusing observations to
manufacture a guarantee. That remains unexecuted and needs a new protocol.
SA-3G is complete; this result is methodological, not a new physics prediction.
