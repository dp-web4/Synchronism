# SA-3E results: reaching a gate and having power to reach it are different budgets

Codex, 2026-09-08. Exact synthetic control, not a new physics prediction or
a fitted explanation of every SA-3D trajectory.
[Registration](2026-09-08-codex-sa3e-evidence-budget-charter.md) committed and
pushed at `826005b9` before calculation. Threshold and all requested budget
curves were retained; no policy was changed.

**Eight audits can be enough to make reopening conceivable and still be far
from enough to make it probable. After evidence against a signal accumulates,
eight can become insufficient even in the best possible case.**

## What is exact here

The useful relation is supplied in advance: parity. A paid paired-channel
audit agrees with the target with probability 0.82 under parity and 0.5
under independent noise. Its likelihood ratio multiplies by 1.64 on agreement
and 0.36 on disagreement. Hold the evidence threshold at 20 and stop auditing
on crossing or budget exhaustion.

For four inherited evidence states and every cap from 0 through 32, rational
dynamic programming tracks all uncrossed states and absorbs crossing paths.
The reported probabilities are exact for this binary model, not estimates
from a seed cohort. It deliberately omits SA-3D's competing hypotheses,
initial exploration, and audit cadence; it is not a numerical replay of that
policy or a universal bound on its power.

## Reachability versus power

| Inherited evidence | Minimum possible audits | Cap for 50% power | Cap for 80% power | Cap for 90% power |
|---|---:|---:|---:|---:|
| 0.01 | 16 | 31 | Not reached by 32 | Not reached by 32 |
| 0.1 | 11 | 20 | Not reached by 32 | Not reached by 32 |
| 1 | 7 | 13 | 22 | 28 |
| 10 | 2 | 2 | 5 | 8 |

“Power” means probability of crossing by the cap when the supplied parity
alternative is true. The minimum possible column assumes every audit is
favorable. It is a feasibility check, not a sensible confidence in recovery.

For example, starting at one needs seven consecutive agreements to cross
at the first possible time. Eight audits offer only 24.93% power, despite
being above that seven-audit minimum. Starting at 0.1 makes an eight-audit
crossing impossible; 20 audits are needed just to exceed 50% power.

The full coarse view of the exact curves:

| Inherited evidence | Power by 8 audits | By 16 | By 24 | By 32 |
|---|---:|---:|---:|---:|
| 0.01 | 0% | 4.18% | 22.55% | 53.47% |
| 0.1 | 0% | 26.28% | 63.58% | 76.30% |
| 1 | 24.93% | 72.21% | 85.28% | 94.02% |
| 10 | 90.41% | 96.05% | 98.69% | 99.32% |

At a cadence of one opportunity every eight rounds, 22 audit opportunities
would require roughly 176 rounds, before considering late stopping. That is
a different time budget from the 64-round policy we tested. The calculation
does not authorize silently supplying those extra rounds.

## Cost remains a real constraint

Each audit costs 0.12 in the existing prediction-error utility units, not
currency. These are acquisition costs ONLY; no prediction-error term is
included in this control's cost table.

| Inherited evidence | Cap | Expected spend if parity | Expected spend if noise |
|---|---:|---:|---:|
| 0.1 | 8 | 0.960 | 0.960 |
| 0.1 | 32 | 2.625 | 3.835 |
| 1 | 8 | 0.930 | 0.959 |
| 1 | 32 | 1.744 | 3.771 |

When a real signal produces a crossing, measurements stop early. Noise
usually consumes nearly the whole cap. Buying more evidence therefore
improves detection power at substantial cost under the null. The fact that
an eight-audit budget cannot produce a certificate from state 0.1 does NOT
mean those measurements cannot improve their own current predictions; that
benefit is outside this power/cost control and was separately scored in SA-3D.

## The conditional error budget is not always 5%

From inherited evidence e, the conditional future crossing bound is e/20.
It is 0.05% at e=0.01, 0.5% at e=0.1, 5% at e=1, and 50% at e=10.
For example, from e=10 the exact null crossing probability by eight audits
is 33.98%. That does not violate a continuously maintained record's initial
5% bound: arriving at such a high-evidence state is itself uncommon under
the null. It WOULD be wrong to initialize a fresh record at ten and claim
that it has the same unconditional 5% guarantee as one initialized at one.

This distinction matters for continuity. Restarting the evidence calculation,
or selecting records based on their current evidence, changes the accounting.
Preserve the accumulated evidence and its error-budget scope; do not erase
unfavorable history just to make a gate easier to cross.

## What this changes in my proposed direction

The next design needs two checks, not just an evidence threshold:

1. **Feasibility:** can any allowed sequence of remaining measurements reach it?
2. **Power at cost:** under the declared useful alternatives, how likely is
   crossing within the allowed expense and time?

A policy that passes the first can still be practically inert. A policy
that improves the second can become too costly in noise. This is a concrete
design constraint for horizon management, not a newly universal horizon law.

I would carry a “reopening feasibility” field alongside a stop decision:
remaining experiments, inherited evidence, maximum attainable evidence,
declared detection-power assumptions, and the cost of acquiring it. A failed
feasibility check should say “this recovery contract cannot be met with the
remaining budget,” not “the outside contains no useful signal.”

This rung produces the calculation needed for that field, not a deployed
interface or an optimized acquisition policy. Unknown useful relationships,
changing regimes, composite noise, and real-agent performance remain open.

## Verification and artifacts

- [Exact instrument](../simulations/mrh_sa3e_evidence_budget.py), standard library.
- [Full record](../simulations/mrh_sa3e_evidence_budget_results.json): 264 exact
  curve points with rational and decimal probabilities, expected counts,
  spend, terminal evidence, and predeclared power-cap summaries.
- 2,147 assertions pass: exact mass conservation, bounded optional stopping
  including overshoot, conditional null bounds, monotonicity, reachability,
  and deterministic endpoints. All binary paths through cap eight are
  independently enumerated and agree with the dynamic program.
- Source compiles; fraction/display schemas and source hash match; full
  repeat is byte-identical. No controls failed or settings were changed.
- Source SHA-256:
  `8b9627af221e83e6b123aa41720eaf5aa05f3c86e36cb872854b7dff2c68f5c0`.
- Result SHA-256:
  `f9a504dc90b7876e555254bbc4073f9c0bf758713a93599a6bc06a64fbc932ac`.

```sh
python3 simulations/mrh_sa3e_evidence_budget.py --controls-only
python3 simulations/mrh_sa3e_evidence_budget.py
```

Output is stdout; `--output NEW_PATH` optionally creates a new result without
overwriting. This exact-control rung is complete; no new policy was evaluated.
