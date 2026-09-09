# SA-3H results: calibration becomes the dominant acquisition budget

Codex, 2026-09-08. [Protocol](2026-09-08-codex-sa3h-paid-calibration-charter.md)
committed and pushed at `1242022f` before implementation/evaluation.

Finite calibration can justify the previously supplied bounds, but under
this conservative protocol **authorization becomes the bottleneck**. A good
test after authorization is not necessarily an effective end-to-end workflow.

With true background agreement 0.50 and signal agreement 0.82, buying **512
samples from each source** authorizes deployment only **37.85%** of the time.
At the 96-audit deployment cap, the fresh test's conditional detection is
**99.71%**, but end-to-end detection is only **37.74%**. Buying **2,048 per
source** raises end-to-end detection to **96.99%** at the price of **4,096
calibration measurements** before deployment begins.

This is not an optimal sample-complexity bound, nor a real-agent benchmark.
It is an exact accounting of one registered procedure that refuses to make
its calibration assumptions free.

## What the observer is allowed to know

Two separate labeled IID calibration sources provide background and signal
agreement bits. Their identities are supplied; their rates are unknown to the
observer. Each sample costs 0.12 in the synthetic accounting. Both n-sample
batches are paid in full, even when calibration later fails.

Only calibration counts enter the certificate decisions:

- A sufficiently low background count certifies the predeclared ceiling
  p0<=0.55, using a one-sided exact-binomial tail budget of 0.005.
- A sufficiently high signal count certifies the predeclared floor p1>=0.75,
  with its own 0.005 budget.
- Without the ceiling certificate, abstain and buy no deployment audits.
- With it, start one **fresh independent** test at evidence 1, assumed
  alternative q=0.82, factors q/0.55 and (1-q)/0.45, threshold **25**.
- The signal certificate affects the power report only. Advertise at least
  80% conditional detection only when both bounds are certified and the
  exact floor power at p=0.75 reaches 80% within the deployment cap.

The count rules are exact one-sided binomial test inversions, following the
principle behind [Clopper–Pearson confidence bounds](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm).
We test fixed bounds rather than estimating arbitrary endpoints. Each batch
size is fixed in advance; the six sizes are distinct experiments, not six
opportunities to stop a growing sample when it finally passes.

The 1% combined calibration-error allowance plus 4% deployment-test allowance
gives the registered conservative 5% budget for **any false calibration claim
OR a null deployment crossing**. This is a repeated-experiment guarantee,
not a posterior probability that a passing bound is true.

The threshold changed from 20 to 25 to reserve calibration risk. All
counterfactual comparisons below use that same threshold 25, so calibration
cost is not confused with the change in deployment error budget. Calibration
data are never recycled as fresh deployment evidence, and no old record is reset.

## Availability is a separate probability from detection power

For p0=0.50, p1=0.82 and deployment cap 96:

| Samples per source n | Deployment authorized | Both bounds certified / 80% report issued | End-to-end detection | Calibration spend | Expected total signal-side spend |
|---|---:|---:|---:|---:|---:|
| 0 | 0% | 0% | 0% | 0 | 0 |
| 8 | 0.39% | 0% | 0.39% | 1.92 | 1.930 |
| 32 | 1.00% | 0.014% | 1.00% | 7.68 | 7.705 |
| 128 | 6.63% | 1.39% | 6.61% | 30.72 | 30.888 |
| 512 | 37.85% | 33.39% | 37.74% | 122.88 | 123.838 |
| 2,048 | 97.27% | 97.27% | 96.99% | 491.52 | 493.981 |

The same fresh test, with bounds supplied for free, detects this signal with
probability 99.71% and uses an expected 21.08 deployment audits (spend 2.530).
That is a diagnostic counterfactual, not a cost-free competitor that learned
its own calibration. At n=2,048, calibration is over 99% of expected total
spend. No amortization across future tasks was assumed.

No claim is made that 2,048 is the minimum useful n: sizes between 512 and
2,048 were not searched. Nor is this symmetric allocation optimal. The
calibration sources have different gaps from their boundaries, so equal
sample counts need not be efficient. The frozen protocol exposes the cost;
it does not optimize it away.

## More calibration does not repair too few deployment audits

The exact conditional floor power and actual p=0.82 power under threshold 25:

| Deployment cap | Certified-floor power if p=0.75 | Actual conditional power if p=0.82 | Can issue an 80% floor-power report? |
|---|---:|---:|---|
| 24 | 40.53% | 70.18% | No |
| 64 | 77.70% | 97.97% | No |
| 96 | 87.88% | 99.71% | Yes, only after both calibration certificates |

The floor first reaches 80% at **68 audits**. Thus no calibration size in
this experiment can justify the requested 80% floor-power claim at the
primary 24-audit cap, or even at 64. More certainty about the model does not
create a longer test budget. That remains distinct from calibration failure.

The 87.88% floor at cap 96 is **conditional on an authorized test and a
signal source satisfying the floor**. It does not contradict the 37.74%
end-to-end result at n=512: many workflows abstain before any test begins.

## The boundary case is not cured by more samples

At n=2,048, the background certificate requires **at most 1,067 agreements**;
the signal certificate requires **at least 1,587 agreements**.

When the actual background is p0=0.55—exactly the valid ceiling—the chance of
certifying that ceiling is only **0.4495%**. At p0=0.50 it is 97.27%. The
deployment rule would be safe in both environments, but the empirical
certification procedure almost never authorizes the boundary environment.

This is not merely slow convergence in this procedure. Its exact rule
ensures `P_0.55(pass) <= 0.005` at every n by construction. The corresponding
signal-floor certificate at true p1=0.75 also passes with probability at most
0.005. A safety specification and a routinely certifiable operating point
need separation; treating them as interchangeable creates permanent
unavailability at the boundary.

This does not mean the boundary environment lacks signal, or that every
alternative calibration/deployment design is equally unavailable. It names
the cost of demanding this particular one-sided certificate before acting.

## Unsafe calibration sources remain visible

The registered out-of-bound cases are p0=0.60 and p1=0.70. They were not
discarded from the results. False-bound certification respects its 0.5%
budget; the two-bound failure respects 1%; the joint bad-event probability
respects the conservative 5% allocation throughout the grid.

There is an important selection caveat. At p0=0.60 and n=2,048, erroneous
authorization occurs with probability about **2.716e-13**. But **among those
rare authorizations**, null crossing by 96 has probability **9.95%**, because
the fresh background still violates the assumed 0.55 ceiling. Small
unconditional error does not imply small conditional error in every selected
subpopulation. The guarantee explicitly accounts for erroneous certificates;
it does not declare them impossible or turn a certificate into certainty.

Likewise, source stationarity and labels remain assumptions. The finite
record cannot prove that an unknown future task comes from the calibrated
source. Changing tasks, changing sensors, or shifting background rates needs
another scope argument, not just a retained JSON certificate.

## Verification and artifacts

[Implementation](../simulations/mrh_sa3h_paid_calibration.py),
[exact results](../simulations/mrh_sa3h_paid_calibration_results.json),
[tests](../simulations/test_mrh_sa3h_paid_calibration.py).

The artifact contains all **216 configuration records**, six calibration
cutoff/probability summaries, and **679 fresh-stream cap records**. Counts
and probability distributions are enumerated analytically, not sampled;
"paid" means charged in the synthetic cost model, not external expenditure.
Exact rational values accompany display decimals, including very small
failure probabilities and probabilities that round to one.

8,902 instrument assertions pass: exact binomial tails, rate-grid false-bound
coverage, cost and probability identities, independent short binary paths,
and SA-3G scaled-threshold regression. Four additional test groups pass,
including independently summed calibration count pairs, Boolean-event union
enumeration, and complete artifact replay. A second full execution is
byte-identical. No failed controls or post-result policy tuning.

- Source SHA-256: `425226f4d7e9017060cfc057bd10dc2d59217cd9ccdd4352443b274424efd3fe`.
- Artifact SHA-256: `468b017356d507039d50f9b1ce952093cb4b4914321bf0ac50b383c4fb0b331c`.

```sh
python3 -B simulations/mrh_sa3h_paid_calibration.py --controls-only
python3 -B simulations/mrh_sa3h_paid_calibration.py --output NEW_PATH.json
python3 -B -m unittest discover -s simulations -p test_mrh_sa3h_paid_calibration.py -v
```

The output option refuses overwrite. Earlier instruments and their artifacts
are unchanged. No real-agent integration or physics-ledger change.

## My perspective and the next useful question

The chain is now explicit: representation determines which signal can be
expressed; measurement determines which evidence can be acquired; an evidence
budget determines whether a decision can change; calibration determines
whether its safety and power claims are justified; **availability determines
whether that otherwise good procedure ever gets to run**.

The practical lead is not "recalibrate everything from scratch for every
decision." It is to ask whether a scoped calibration record can be reused
across genuinely matching tasks, while paying one calibration cost and
tracking the error budget across repeated deployments. That could amortize
the dominant cost, but only with explicit source continuity and lifetime
risk accounting—repeated fresh 4% tests are not one fleetwide 4% guarantee.
That is a next experiment, not a result claimed here.

Assumption still accepted: labeled IID sources match deployment. Practice
checked: exact fixed-size confidence claims do not license optional peeking
or evidence recycling. Operator objection upheld: this safe procedure is
mostly unavailable at small calibration budgets and expensive when available.
No substrate axiom was tested. SA-3H is complete.
