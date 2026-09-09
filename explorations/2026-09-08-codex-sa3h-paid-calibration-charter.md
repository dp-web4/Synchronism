# SA-3H registration: pay for the calibration certificate

Codex, 2026-09-08. dp approved continuing from SA-3G. Commit and push before
implementation/evaluation. This is a bounded exact experiment, not deployment.

## Question and access assumption

Can a finite, paid calibration record justify the null ceiling and power
floor previously supplied for free, and how much availability does it cost?

Assume access to TWO labeled IID calibration sources: background Bernoulli
agreement rate p0 and signal agreement rate p1. Source identity is known;
the rates are not given to the observer. Each sample costs 0.12, just like
a deployment audit. Both batches are bought in full before testing. No
reuse of their observations as deployment evidence. This access to labeled
sources is a strong supplied assumption, not a solved representation or
ground-truth problem. Calibration and deployment must be independent and
stationary within their respective source; no out-of-distribution claim.

This first step certifies predeclared bounds, rather than fitting arbitrary
new intervals or tuning a test to the resulting data.

## Frozen protocol

- Buy n background and n signal samples, n in {0,8,32,128,512,2048}.
- Fixed bounds: background ceiling u=0.55; signal floor l=0.75.
- Each one-sided calibration claim gets failure budget delta=0.005.
  From background count X, certify p0<=u only if
  `P[Binomial(n,u) <= X] <= delta`.
  From signal count Y, certify p1>=l only if
  `P[Binomial(n,l) >= Y] <= delta`.
  Empty batches certify neither bound. Use exact rational binomial tails;
  the rejection regions are monotone integer-count thresholds.
- If the null ceiling is not certified, abstain: no deployment audits, no
  detection, no safety-certified deployment. This is unavailability, not
  evidence that the environment has no signal.
- If the null ceiling is certified, run a FRESH independent audit record,
  evidence start=1, assumed alternative q=0.82, update factors q/u and
  (1-q)/(1-u). Gate threshold=25 (test error budget alpha=0.04). Stop at
  first crossing or m in {24,64,96} audits.
- The signal-floor certificate affects reporting only, not audit scheduling.
  Advertise >=80% conditional detection power only if both bounds are
  certified AND exact first-passage power at p=l within m is >=0.8.
  Otherwise say power uncertified; do not substitute q for a measured floor.
- Every (n,m) is a separate fixed-size experiment, not repeated looks at one
  growing calibration record. No "keep calibrating until it passes" without
  a new sequential error budget. No reset of an already-used evidence record.

Grid of true rates: p0 in {0.50,0.55,0.60}; p1 in {0.70,0.75,0.82,0.90}.
The p0=.60 and p1=.70 rows are explicit out-of-bound controls. Full grid:
6 sample sizes x 3 null rates x 4 signal rates x 3 deployment caps = **216
configuration records**. Enumerate calibration-count probabilities exactly,
then combine them with exact fresh-stream first passage. No Monte Carlo,
seed search, optimized allocation, threshold scan, or significance selection.

## Claims and their scopes

The calibration rules are one-sided exact-binomial test inversions, the same
principle underlying [Clopper–Pearson confidence limits](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm).
Each has probability <=0.005 of certifying a false bound. The two claims
jointly fail with probability at most 0.01 (union bound; independence is
available but not required for that inequality).

Given a valid null bound and a fresh stream from that background, the
likelihood product is a nonnegative supermartingale and threshold 25 gives
<=0.04 crossing probability. Calibration-error plus deployment-error
accounting therefore supports a conservative 0.05 joint budget for a false
calibration claim OR a null deployment crossing. Do not claim that 5% is the
sharp achievable bound, or spend unused conservatism after seeing outcomes.

These are repeated-experiment guarantees, NOT posterior probabilities that
a bound is true after seeing a passing certificate. Nor is unconditional
false-crossing control automatically conditional control among the selected
deployments: rare erroneous authorizations can concentrate there.

Power is conditional on deployment and an IID signal source covered by the
calibration claim. End-to-end detection additionally pays the probability
that calibration authorizes any deployment at all. The 80% statement does
not promise 80% end-to-end detection across calibration failures/abstentions.

## Required outputs

For each n, retain exact count acceptance thresholds and exact certification
probabilities at every tested p0/p1, including false-certification events.
For each m and source rate, retain exact crossing probability, expected
deployment audit count/spend, and the l=.75 power floor. Include the full
0..96 curves for independent replay, not just the chosen caps.

For all 216 configurations, report:

- probability of deployment authorization and abstention;
- probability of both calibration bounds being certified;
- probability of issuing the >=80% conditional-power report;
- end-to-end null crossing and signal detection, separately from conditional
  fresh-stream crossing/detection;
- probability of any false calibration claim, including the deliberately
  out-of-bound controls, and its union with null deployment crossing;
- fixed calibration cost 2n*0.12; expected deployment and total cost under
  null/signal; maximum allowance 2n+m. No silent amortization across tasks.

Compare against the same fresh-stream test with the bounds supplied for free
as a diagnostic counterfactual. Its bound is valid only on matching rate rows;
it is not a zero-cost achievable learned-calibration competitor. The primary
deployment cap remains 24; 64/96 expose the longer-horizon tradeoff. No
claim of optimal sample complexity from this conservative construction.

## Controls and stop

- Exact binomial mass conservation, endpoint counts, complementary tails,
  monotone acceptance regions, maximal/minimal cutoff boundary checks.
- Independent small-n enumeration of Bernoulli paths; check false-bound
  certification <=delta throughout a fixed p grid in [0,1] at steps 0.01.
- Independent small-horizon deployment path enumeration, threshold handling,
  stopped moments, cost/count bounds, and first-passage monotonicity.
- Regress against SA-3G by using its start=4/5 and threshold=20, which is
  pathwise equivalent to this start=1 and threshold=25. Do not mutate SA-3G.
- Verify probability factoring over independent calibration and deployment,
  joint calibration failure <=0.01 and joint bad-event budget <=0.05,
  truthful report conditions, cost identities, schemas, and source hashes.
- Compile, run tests, repeat the artifact byte-identically.

Assumption still accepted: correctly labeled, stationary calibration sources
represent future use. Standard-practice hazard: optional calibration peeking
or recycling training evidence as a fresh sequential test. Operator objection:
calibration and abstention may cost more than the recovery they enable. Report
that failure if it occurs; do not make samples free or optimize the grid.
No physics prediction, full-planner improvement, or real-agent integration.
Stop after the rung and its implications are documented and pushed.
