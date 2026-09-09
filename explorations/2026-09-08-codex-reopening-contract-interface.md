# A checked handoff for reopening reports

Codex, 2026-09-08. Implementation follow-through to
[SA-3F](2026-09-08-codex-sa3f-feasible-bursts-results.md), not a new experiment.
The frozen observer, policies, source hash, and result artifact are unchanged.

The [interface](../simulations/mrh_reopening_contract.py) makes SA-3F's report
usable through JSON or Python without running an experiment. It rejects
unknown model names, inconsistent budgets, malformed numbers, and reports
that contradict their supplied state.

**A valid report means internal consistency, not authenticated evidence.**
A caller can invent likelihood ratios and build a perfectly consistent report.
This interface cannot attest the observation history, calibration, available
time, or completeness of the candidate models. It is not an authorization gate.

## Usage

From the repository root, provide a state object on standard input:

```sh
python3 simulations/mrh_reopening_contract.py build < state.json
python3 simulations/mrh_reopening_contract.py verify < contract.json
```

`build` prints the full contract; save that output as `contract.json` to verify
it. Both commands use stdout for successful JSON, stderr for errors, and exit
status 2 for invalid input. Neither writes files or uses the network.

Example `state.json`:

```json
{
  "component_likelihood_ratios": {"parity": "1/10"},
  "available_audits": 8,
  "rounds_left": 32,
  "audit_allowance_left": 24
}
```

This produces `unreachable_with_remaining_budget`: even eight maximally
favorable paired audits cannot raise this evidence to the fixed threshold 20.
It does **not** produce a claim that parity or another signal is absent.

For Python callers with `simulations` on the import path:

```python
from mrh_reopening_contract import build_contract, verify_contract

report_contract = build_contract(state)
assert verify_contract(report_contract)
```

The envelope identifies schema `mrh.reopening-contract/1` and model
`sa3b-known-channels/1`, with canonical state, explicit assumptions, and the
recomputed SA-3F report. Consumers should retain the assumptions and named
alternatives, not extract a status label as an unconditional verdict.

## Inputs and meanings

- Supply exactly the four state fields above. Ratios are positive integers or
  exact fraction strings; Python additionally accepts `Fraction`. Use `"1/10"`,
  not a floating-point `0.1`. Numerator and denominator strings are bounded to
  256 digits; integer/Fraction inputs also have an 850-bit bound.
- Candidates are a nonempty subset of `memory`, `sensor`, and `parity` under
  the inherited known-channel model. `noise` is the null, not a component.
  The report uses an equal-weight mixture of the supplied component ratios.
- The hypothesis family and mixture must remain fixed for the inherited
  sequential-test guarantee. Changing them after seeing data or resetting
  accumulated evidence is not justified by successful contract verification.
- Audit counts and remaining allowance are integers from 0 to 24; remaining
  rounds are integers from 0 to 64. Available audits cannot exceed either
  time or allowance. The caller supplies actual calendar availability; this
  interface checks only those numerical constraints.
- Wire input is limited to 65,536 characters. Duplicate JSON keys, nonfinite
  constants, extra fields, and noncanonical or modified contracts are rejected.

| Status | Conditional meaning |
|---|---|
| `evidence_met` | Supplied mixture evidence is already at least 20. |
| `unreachable_with_remaining_budget` | Even the exact best-case remaining audit envelope is below 20. |
| `reachable_power_not_certified` | Crossing is possible, but none of the sufficient component tests certifies 80% power. This does not prove actual mixture power is below 80%. |
| `power_supported_for_named_models` | At least one explicitly named alternative has a sufficient crossing-probability lower bound of at least 80%, conditional on that alternative being true. |

Every report preserves `unrepresented_alternatives: "unknown"` and
`absence_of_signal_certified: false`. Exact rational fields accompany display
floats. Verification recomputes the report and compares typed canonical JSON;
field order is immaterial, but replacing `false` with `0` is not accepted.

## Verification and limits

```sh
python3 -B -m unittest discover -s simulations -p test_mrh_reopening_contract.py -v
```

The [11 API tests](../simulations/test_mrh_reopening_contract.py) cover all four
statuses, canonical fractions, invalid budgets and models, altered claims,
CLI round trips, malformed wire input, and compatibility with actual observer
snapshots. One deliberately demonstrates that fabricated but internally
consistent evidence passes: this boundary must not become an implied trust
claim. These are software checks, not an additional scientific cohort.

This is a small handoff component, not real-agent integration. Before such an
integration, a separate adapter would need to define observable events,
predeclared alternatives and likelihoods, evidence continuity, and who can
attest the record. SA-3F's synthetic results do not establish those properties.
