# SA-3G registration: which calibration error breaks which promise?

Codex, 2026-09-08. dp approved continuing with imperfect likelihood
calibration. Commit and push this protocol before implementation/evaluation.

## Question and scope

Separate three claims that a single word, "calibrated," can obscure:

1. Is threshold crossing rare under the declared null?
2. Is the stated detection probability valid under the named alternative?
3. Can the fixed evidence rule reach its threshold within the audit budget?

A normalized but wrong alternative need not invalidate a likelihood-ratio
test under a correct null. A wrong null can invalidate that guarantee.
Reachability is a pathwise property of the update rule, distinct from either
probabilistic claim. We will quantify the separation and the price of a
predeclared composite-null guard.

Use SA-3E's exact binary audit-stream reduction, **not** a new full-planner
benchmark: Z=1 means measured parity agrees with the subsequently revealed
target. A paired audit costs 0.12. In the inherited calibrated parity model
P(Z=1)=0.82; under independent fair noise it is 0.5. General Bernoulli p here
is an effective agreement rate, not a unique decomposition into sensor noise,
target noise, or another physical mechanism. No unpurchased information,
free calibration sample, new hypothesis, or real-agent integration.

## Frozen grid and evidence rules

- Starting evidence: 1/10 and 1.
- Threshold: 20; power target: 4/5.
- Assumed alternative agreement q: 13/20, 41/50, 9/10 (0.65, 0.82, 0.90).
- True IID agreement p: 1/2, 13/25, 11/20, 3/5, 13/20, 3/4, 41/50, 9/10
  (0.50, 0.52, 0.55, 0.60, 0.65, 0.75, 0.82, 0.90).
- Two fixed rules, indexed by null ceiling u:
  - nominal: u=1/2;
  - guarded: u=11/20.
- On agreement multiply by q/u; on disagreement by (1-q)/(1-u).
  Stop purchasing at the first evidence >=20, otherwise at the cap.
- Compute every cap 0 through 64. Primary comparison cap=24 (SA-3F's
  allowance); cap=8 is the old scarce-budget reference; cap=64 is a longer
  stress horizon, **not** permission to relabel it as a 24-audit improvement.

This is 96 exact curves, 6,240 cap records; no seed selection or Monte Carlo.
Report all combinations, not only settings that show an error or a benefit.
The guard and q values will not be fitted to the resulting curves.

The guarded rule tests a **different, broader null**: conditional agreement
probability at each purchased audit is <=0.55, given prior audit information.
It intentionally treats some weak positive association as null. Such an
association is not "independent noise" and detecting it with the nominal
rule is not a false alarm under the nominal rule's original point null.
Any claim of a 5% bound under that broader class must use the guarded rule.

Both starts are mathematical input states. The 1/10 start represents a
conditional continuation only when inherited from a valid continuous record
under the SAME rule/null. One cannot convert SA-3F's nominal evidence into
guarded evidence by changing a label, or reset evidence to recover power.
Starting at 1 gives a 5% bound; starting at 1/10 gives a conditional 0.5%
remaining crossing bound under the applicable null. These are not fleetwide
error budgets.

## Analytical obligations, before interpretation

For any assumed q>u, the one-step mean factor at true p is

`p*q/u + (1-p)*(1-q)/(1-u)`.

It is 1 at p=u and <=1 for p<=u. Thus the guarded product is a
nonnegative supermartingale for the conditional null class, even if p varies
with history; the nominal product remains a martingale under fair noise
even when q is a wrong alternative. The anytime crossing bound is start/20.
This is established sequential-testing machinery, not a new theorem:
[Ramdas et al., Game-theoretic statistics and safe anytime-valid inference](https://arxiv.org/abs/2210.01948).

The exact rule-level reachability envelope after m audits is
`start*(q/u)^m`, attained by all agreements. Changing the true p without
changing the rule does not change this envelope when both binary outcomes
remain possible; changing the rule to guard the null can change it.

Crossing is coordinatewise increasing in the audit outcomes (q>u), so among
IID alternatives p in [3/4,9/10], its lowest probability is attained at 3/4.
For every start/q/rule/cap report (a) nominal power evaluated at p=q and
(b) this interval-floor power, naming the interval explicitly. The interval
is supplied, not learned or statistically certified from a calibration set.
No detection certificate outside it. Failure of the floor to reach 80% does
not mean all members are below 80%; report the distinction.

## Outputs and controls

For each cap retain exact rational and display-decimal crossing probability,
survival probability, expected audit count/spend, expected stopped evidence,
best-case evidence, and whether the applicable null bound covers this p.
Retain first caps for 50%/80%/90% crossing (or null if absent by 64).
Separate the observer's assumed-q power from actual-p power and count/report
settings where a nominal 80% claim would be overstated under a tested p.
Include point-null and guarded-null ceilings in every report.

Controls:

- Independent binary-path enumeration through cap 8 checks first passage,
  expected cost, and stopped evidence, including deterministic endpoints.
- Exact probability conservation, monotonicity in cap and true p, budget
  limits, threshold handling, all-success envelope, and zero power before
  reachability.
- Exact one-step normalization/supermartingale checks. Under p=u, expected
  stopped evidence equals start; under p<=u it is <=start, with crossing
  bounded by start/20. Do NOT assert that bound under p>u.
- Small-horizon dynamic adversarial choices p in {0,u} check that the
  worst crossing probability over history-dependent null choices matches
  the IID p=u boundary (affinity permits endpoint maximization).
- SA-3E regression at q=.82, u=.5, p in {.5,.82}, both starts and caps 0..32.
- Match nominal single-component SA-3F power certificates at caps <=24.
- Compile, check schemas/source hashes, repeat the artifact byte-for-byte.

## Risks, interpretation, and stop

The analytical identities are expectations before evaluation; this is not a
blind test of whether supermartingales exist. The new output is the exact
miscalibration/power/cost map and a clear boundary on the observer's promises.
If mild deviations never breach the advertised bound within these horizons,
say so; invalidating the proof is not itself a measured bound violation.
If guarding loses most useful recovery, retain that cost rather than tune u.

Assumption not relaxed: a supplied uncertainty set meaningfully covers the
environment. Standard practice scrutinized: treating a correct likelihood
ratio calculation as a calibrated scientific claim. Operator objection:
synthetic calibration bounds are not available for free in a real agent.
No substrate dynamics or physics predictions are tested; ledger unchanged.
Stop after this rung. Learned calibration, alternative discovery, full-planner
transfer, and deployment need another design.
