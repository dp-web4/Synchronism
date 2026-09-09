# SA-3E registration: exact evidence-budget feasibility

Codex, 2026-09-08. Continuation requested by dp; register and push before
calculation. SA-3D's threshold is not being tuned or replaced.

## Question and deliberate simplification

How much distinguishing evidence can the fixed threshold require after an
observer has accumulated evidence against a signal? Separate absolute
reachability from useful detection probability. Compute exact probabilities,
not another Monte Carlo policy contest.

This is a simplified binary test, NOT a replay or fitted model of SA-3D's
multi-hypothesis states. The useful candidate is supplied: parity. For a
paired observation define agreement A = [Y equals H xor noisy_sensor].
The existing worlds imply independent A~Bernoulli(1/2) under noise and
A~Bernoulli(41/50) under parity. Their per-audit likelihood-ratio factors
are 41/25 on agreement and 9/25 on disagreement. An audit costs 0.12.

Freeze threshold at 20 and starting evidence at {1/100, 1/10, 1, 10}.
Interpret the starting value as the inherited state of a continuous evidence
record, not permission to initialize a fresh 5%-level test with wealth 10.
Evaluate EVERY audit cap M=0,...,32. Stop measuring at first crossing or
exhaustion of the cap. No outcome-based expansion beyond 32.

## Exact computation and reported metrics

Use rational dynamic programming over uncrossed states (number of audits,
number of agreements); keep crossed mass absorbing. For both null and
alternative, report first-crossing probability by every M, expected number
of audits before crossing/cap, and expected acquisition spend. Retain exact
fractions and floating-point display values.

Report the smallest cap in 0..32 achieving alternative crossing probability
at least 50%, 80%, and 90%; report unattained thresholds as null. Also report
the minimum number of all-agreement outcomes needed to reach 20. These are
predeclared summaries of the full curve, not selection of a new audit policy.

The conditional null bound from state e is e/20, NOT always 5%. A full
record initialized at one has the familiar 5% bound. Conditioning on prior
evidence changes the remaining risk budget; selecting/resetting records
requires separate accounting. No claim about unknown noise or discovered
representations is made.

## Independent controls

- Exact outcome probabilities normalize, and the null expectation of the
  one-step likelihood-ratio multiplier equals one.
- Total absorbed plus surviving probability equals one at every cap.
- Under the null, expected evidence at the bounded stopping time equals
  starting evidence exactly, including crossing overshoot.
- Crossing probability obeys e/20; power and expected audit count are
  nondecreasing in cap; expected count lies in [0,M].
- Brute-force all binary paths through cap 8 and compare crossing probability
  and expected audit count with the dynamic program, for each start/world.
- Zero power below the all-agreement reachability minimum; exact endpoint
  checks for always-agreement and always-disagreement generators.
- Compile, preserve source hash, and repeat output byte-for-byte.

The disclosed assumption is an idealized known useful channel and independent
audit outcomes. Expected spend here is acquisition cost ONLY: do not confuse
it with total prediction utility, which may benefit on the audit rounds.
Cadence, initial exploration, hypothesis mixtures, free-label availability,
and uncertain noise still matter in a full agent policy and are not optimized.

Prior art is ordinary likelihood-ratio sequential testing / absorbing-state
dynamic programming, with the anytime-valid context already cited in SA-3D.
No substrate axiom or physics prediction is tested. Stop after the exact
feasibility map and its interpretation; no audit-policy modification in this rung.
