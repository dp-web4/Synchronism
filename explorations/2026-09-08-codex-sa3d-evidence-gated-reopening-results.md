# SA-3D results: a valid gate can be too expensive in evidence to reopen

Codex, 2026-09-08. Synthetic decision experiment; no new physics or agent-
benchmark claim. [Protocol](2026-09-08-codex-sa3d-evidence-gated-reopening-charter.md)
was committed and pushed at `a323d705` before implementation/evaluation.
Threshold, schedule, seeds, models, and prices were not tuned after results.

**The evidence gate suppressed the specified false-reopening event, but
rescued none of the stopped informative runs in the new main cohort.**
The mathematical guarantee is useful and narrow. The combined recovery
policy, at this audit budget, is not a successful replacement for SA-3C.

## Intervention and the event actually controlled

The periodic policy buys both channels every eighth round if it would
otherwise stop. The gated version additionally requires E>=20 before an
unforced query after the base planner has first recommended stopping.
Initial exploration remains unrestricted. Audits remain paid and exempt.
Every purchased observation and later target updates E from the beginning;
there are no restarts, favorable-window selections, or hidden truth inputs.

E is the uniform mixture of cumulative likelihood ratios of the supplied
non-noise hypotheses against the specified independent-noise world. The
threshold is checked before buying the current observation. A final
recommendation is an unexecuted two-round-continuation diagnostic, not an
extra scored round.

## Larger noise cohort: fewer false reopenings, not zero cost

Each cell below uses 512 paired fresh seeds, 64 rounds each. A false reopening
here means ANY unforced non-stop action after the first stopping recommendation,
not merely an incorrect final recommendation. These are different events.

| Family | Periodic: ever reopened | Gated: ever reopened | Mean utility periodic → gated |
|---|---:|---:|---:|
| Restricted | 342/512 (66.8%) | 0/512 | 0.526588 → 0.524811 |
| Expanded | 383/512 (74.8%) | 0/512 | 0.534791 → 0.529054 |

Utility is prediction error plus acquisition price, lower better. Paired
gated-minus-periodic differences are −0.001777 ± 0.000816 (restricted) and
−0.005737 ± 0.000993 (expanded), where ± is one seed-level standard error,
not a confidence interval or a significance declaration.

Zero observed events does not establish zero population risk. The 5% bound
comes from the likelihood-ratio construction and its assumptions, not from
this Monte Carlo count. Gated threshold crossings still occurred in 22/512
restricted runs and 18/512 expanded runs; a crossing is not the same as an
actual post-stop query. Evidence can later decrease, or the observer may not
have stopped yet.

Gated observers still averaged 6.850 (restricted) and 6.779 (expanded) paid
audits per run. They ended with unnecessary autonomous queries in 2/512 and
1/512 runs respectively. All three cases had **never stopped**, so the gate
was deliberately inapplicable. This is a concrete boundary of the guarantee,
not an exception to it. The design does not control every unnecessary query.

## Fresh main cohort: the gate removes the recovery benefit

24 seeds per world/policy; each seed is new relative to SA-3B and SA-3C.
Mean total utility per round:

| Family | World | No audit | Periodic audit | Gated audit |
|---|---|---:|---:|---:|
| Restricted | Memory | 0.1301 | 0.1301 | 0.1301 |
| Restricted | Sensor | 0.3167 | 0.3026 | 0.3124 |
| Restricted | Noise | 0.5133 | 0.5389 | 0.5304 |
| Restricted | Parity | 0.5207 | 0.5378 | 0.5342 |
| Expanded | Memory | 0.1310 | 0.1310 | 0.1310 |
| Expanded | Sensor | 0.3099 | 0.2998 | 0.3081 |
| Expanded | Noise | 0.5205 | 0.5447 | 0.5351 |
| Expanded | Parity | 0.3858 | 0.3650 | 0.3808 |

Correct final autonomous channel choices:

| Family / informative world | No audit | Periodic | Gated |
|---|---:|---:|---:|
| Restricted / memory | 24/24 | 24/24 | 24/24 |
| Expanded / memory | 24/24 | 24/24 | 24/24 |
| Restricted / sensor | 19/24 | 22/24 | 19/24 |
| Expanded / sensor | 20/24 | 22/24 | 20/24 |
| Restricted / parity | 0/24 | 0/24 | 0/24 |
| Expanded / parity | 15/24 | 22/24 | 15/24 |

No gated main-cohort run reopens autonomously after stopping. Its slight
utility improvements over no audit in informative worlds can come from the
predictions made on paid audit rounds; they do not establish sustained
recovery. All memory runs already choose the right channel without stopping,
so this cohort supplies no memory-rescue test cases.

Expanded parity illustrates the price: gating increases utility by
0.015755 ± 0.010863 against periodic auditing and loses seven correct final
choices. Expanded sensor increases utility by 0.008307 ± 0.006377 and loses
two. These small main-cohort estimates are not universal effect sizes.
The restricted family still cannot express parity; neither safeguard fixes
that missing representation.

## Why the probability guarantee is valid—and limited

For an action chosen from past public history and any alternative k,

`L_k(t) = L_k(t-1) P_k(O_t,Y_t | action_t) / P_0(O_t,Y_t | action_t)`.

Under P0, conditional expectation sums the ratio against its denominator:
`E_0[L_k(t) | past] = L_k(t-1) sum_(o,y) P_k(o,y | action_t) = L_k(t-1)`.
The selected action is already fixed at this conditional step. A fixed
mixture retains nonnegativity, initial value one, and the martingale property.
Ville's inequality therefore bounds the probability of ever reaching 20 by
1/20. Every permitted post-stop unforced query requires a preceding crossing.
This applies standard anytime-valid inference; see
[Ramdas et al.](https://arxiv.org/abs/2210.01948), not a new MRH theorem.

The guarantee assumes the declared stationary noise world. It does not
certify closure, identify the right channel, bound losses under alternatives,
or apply unchanged to unknown/composite noise. It is per continuous evidence
record and per family. Repeated resets or multiple observers do not inherit
a shared 5% error budget: a joint claim needs an explicit allocation or other
valid joint construction. The 512 trials here estimate a per-run event rate;
they are not one fleet-wide 5% claim.

## Analytic diagnosis: certainty has an evidence budget

This is interpretation after observing the frozen results, not an additional
evaluated policy. In these likelihood tables a paired audit multiplies any
alternative's likelihood ratio by at most 1.8 (memory) or 1.64 (sensor/parity).
Consequently, after stopping at evidence E_s, even the best possible m audit
outcomes satisfy `E_future <= E_s * 1.8^m` until reopening. Non-audit stops
provide no distinguishing evidence and cannot increase E.

If `E_s < 20 / 1.8^8 = 0.181489`, even eight maximally favorable remaining
audits cannot reach the gate. With fewer remaining audits the necessary
starting evidence is higher. For example, E_s=0.1 needs at least ten audits
even under this optimistic bound, but the entire 64-round schedule contains
only eight. This bound does not claim that every failed run was in that
unreachable region; first-stop evidence snapshots were not registered as a
separate diagnostic and no count is asserted.

The mechanism is broader than an overly conservative constant. Evidence
against a signal accumulated before stopping must be overcome, and strict
error control can demand more measurements than the acquisition policy can
supply. A nominal reopening rule can therefore be practically—or in some
states mathematically—unreachable. Lowering the threshold afterward would
change the risk budget, not repair this registered result.

## What I would carry forward

SA-3A/B/C separated identifiability, representation, acquisition, and reopening.
SA-3D adds a constraint: **a confidence requirement must be priced in evidence,
not just attached to an acquisition policy.** A valid statistical safeguard
can have little recovery power within the available horizon.

An honest next design should report the null and controlled event, the scope
of its error budget, evidence already spent against alternatives, remaining
distinguishing experiments, and whether the threshold is even reachable.
Audit expenditure must remain visible; relabeling every query an "audit"
would evade the controlled event without making the system economical.

I would not ship this gate as a successful recovery policy on this evidence.
The next question is a registered feasibility-aware evidence budget (and its
utility cost), not a claim that calibration solved the stopping problem.
No real-agent transfer or private-repository access occurred in this rung.

## Reproduction and verification

- [Instrument](../simulations/mrh_sa3d_evidence_gate.py), using frozen SA-3B/C code.
- [Full result record](../simulations/mrh_sa3d_evidence_gate_results.json): 576
  main-cohort runs plus 2,048 null-cohort runs, 167,936 scored decisions total.
  Source hashes, seed metrics, audit/veto timing, posteriors, and evidence
  summaries are retained.
- 744 new checks plus 222 SA-3B and 1,028 SA-3C checks pass. They include
  exact rational likelihood/martingale controls and baseline regression;
  they are not 1,994 independent scientific findings.
- Source compiles; all 2,624 record schemas, source hashes, audit limits,
  posterior sums, and prior-evidence gate invariants pass. Full repeat output
  is byte-identical. No implementation control failed before evaluation.
- Source SHA-256:
  `88f00014fa83f10895b77fdcb32c0fb7578a014b42b43ed8ff80104fe05b12ae`.
- Result SHA-256:
  `76aa3024dd0be14e12ca85a47a4d911766defcb446b68fc5526ba4b01ad5fc6f`.

```sh
python3 simulations/mrh_sa3d_evidence_gate.py --controls-only
python3 simulations/mrh_sa3d_evidence_gate.py
```

Outputs go to stdout; `--output NEW_PATH` optionally creates a new artifact
and refuses overwrite. This rung stops here, including its negative result.
