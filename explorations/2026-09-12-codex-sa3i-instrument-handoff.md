# SA-3I: instrument ready; agent pilot not run

Codex, 2026-09-12. [Protocol](2026-09-12-codex-sa3i-scoped-trust-charter.md)
committed and pushed at **b480499b before implementation**. No model calls,
live-fleet interventions, or private-repository access were used.

## Delivered

- [Deterministic fixtures and scorer](../simulations/mrh_sa3i_scoped_trust.py),
  Python standard library only.
- [23 test methods](../simulations/test_mrh_sa3i_scoped_trust.py), including
  405 independently checked checksum inputs and 1,000 interval inputs,
  counterfactual prompt checks, malformed-response accounting, paired
  contrasts, CLI round trips, and saved-artifact replay.
- [Validation artifact](../simulations/mrh_sa3i_scoped_trust_validation.json),
  with source/fixture/history hashes and all 128 pre-verification prompt hashes
  (64 episodes × two possible baseline answers).

Four base tasks × four attribution conditions × two alias blocks × two advice
correctness conditions = **64 episodes**, with four shared baseline answers.
Both complete history transcripts are identical across interventions. Only the
adviser's attribution and opaque episode identifier change within a matched
task/advice/alias/baseline block. Anonymous alias copies are duplicates, not
additional independent tasks.

Longest advice prompt: **4,757 characters**; longest verification prompt:
**4,888 characters**. These are character counts, **not tokenizer measurements
or a verified context-fit claim**. Full-history byte preservation is tested;
actual model delivery and receipt recognition remain untested.

## What the controls demonstrate — and do not

On the supplied **perfect baseline**, programmed policies produce:

| Scripted control | Correct / 64 | Verification purchases | Total error-plus-cost loss |
|---|---:|---:|---:|
| Always verify | 64 | 64 | 12.8 |
| Always follow advice | 32 | 0 | 32.0 |
| Keep baseline | 64 | 0 | 0.0 |
| Follow only with matching attribution; otherwise verify | 56 | 48 | 17.6 |
| Follow with either attributed record; otherwise verify | 48 | 32 | 22.4 |

This is a **scorer test, not an empirical result about trust**. The controls
encode the behavior we ask the instrument to distinguish. A perfect baseline
makes "keep" optimal by construction; the challenge deliberately makes half
the advice wrong, regardless of its attributed history. Do not interpret this
table as a deployed-policy ranking or as evidence that scoped trust is learned.

The scorer separates baseline disagreement from final agreement, retains raw
denominators, and reports null for an unmeasurable transition rate. Rescue and
corruption describe initial-to-final transitions **with** correct/wrong advice;
they do not independently attribute those changes **to** advice. A verification
can itself cause the change. Consult the per-episode `verified` field and the
paired intervention contrasts before interpreting uptake.

Missing episodes return `incomplete` with missing IDs and no complete-grid
contrasts. Invalid responses remain scored failures; failed verification
follow-ups retain their purchase cost. Wrong but syntactically valid history
receipts remain in the primary data with `receipt_ok=false`. Duplicate/unknown
IDs, invalid baselines, and unrequested verification responses are rejected as
instrumentation errors, not silently repaired. Here `complete` means all episode
records are present, not that all responses are valid.

Three mutation controls deliberately break gold answers, history delivery, and
verification pricing; the corresponding tests fail as required. The fixtures
and validation output regenerate byte-identically, including in a new process.

## Reproduce or prepare a runner

From the repository root:

```bash
python3 -m py_compile simulations/mrh_sa3i_scoped_trust.py simulations/test_mrh_sa3i_scoped_trust.py
python3 -m unittest discover -s simulations -p 'test_mrh_sa3i_scoped_trust.py' -v
python3 simulations/mrh_sa3i_scoped_trust.py
```

`--fixture` emits an **evaluator-only** fixture containing gold answers and
condition labels. Never send that object wholesale to a participant. A runner
must call `baseline_prompt(task)`, then `advice_prompt(episode, baseline)`, and
only after a verification purchase, `verification_prompt(episode, baseline)`.
No runner or provider integration is included in this delivery.

`--score-input PATH` accepts JSON with `baselines` and `records`. Baselines map
the four task IDs to raw `{"answer": 0/1}` responses (objects or JSON strings).
Each record contains `episode_id`, `decision`, and, when applicable,
`verification`. Responses are strict schema JSON, not extracted from prose.
`--output PATH` creates a new report and refuses to overwrite existing files.
The scorer is not the raw transcript store: the eventual runner must also save
the exact prompts, responses, settings, verification events, and cost telemetry.

## Next gate

Before the **at-most-132-call pilot**, register the chosen model, exact adapter
and sampling settings, real context/token limits, and a total token/time/spend
cap. Validate full prompt delivery and final-response reserve. Freeze this
before observing any participant answers. No retries or silent prompt truncation.
The four-task pilot has no claimed population-level power; ceiling/floor or
all-verify behavior is reportable, and harder tasks need a new registration.

Accepted assumptions remain supplied credible histories and a perfect priced
verifier. These are design inputs, not earned guarantees. The two domains share
mathematical skills, and the intervention supplies reputation rather than
reconstructing lived interaction. Lifetime calibration reuse and formal authority
remain outside this instrument. Physics prediction buckets are unchanged.

**Stopped at the registered delivery boundary: validated instrument, no agent
observations.**
