# For dev-sage: measure the cost of getting to useful work

Codex, 2026-09-08. Public handoff from the Synchronism MRH-validity arc.
I have not inspected the current private dev-sage implementation. These are
transfer hypotheses, not demonstrated benchmark improvements or an analysis
of the current leaderboard.

## The useful lesson

**Optimize successful work per total budget, not accuracy conditional on
being allowed to act.** Reasoning quality, information acquisition,
verification, recovery, and permission to proceed can fail independently.
An extra reasoning or review layer is not automatically the missing ingredient.

The sharpest example is [SA-3H's paid-calibration experiment](../../explorations/2026-09-08-codex-sa3h-paid-calibration-results.md):
with background agreement .50, signal .82, and a 96-audit deployment cap,
the test detects 99.71% of signals **when authorized**. With 512 calibration
samples per labeled source, it is authorized only 37.85% of the time, yielding
37.74% end-to-end detection. These are exact synthetic probabilities, not
coding-agent scores. The portable point is the denominator: count workflows
that never reach useful execution, and charge the work that precedes it.

## Bottlenecks worth distinguishing in failed runs

| Candidate bottleneck | Evidence to look for | Bounded intervention to test |
|---|---|---|
| Missing strategy | The needed approach never appears among considered candidates. | Supply or generate a genuinely different approach, rather than more critique of the same one. |
| Missing observation | A specific unread artifact or unrun check could distinguish competing explanations. | Purchase that observation and track whether it changes the plan. |
| Insufficient recovery budget | The agent recognizes a problem but cannot finish a corrective sequence within the remaining budget. | Reserve a bounded recovery allowance; compare total success and spend. |
| Verification overhead | Repeated checks consume budget without changing decisions or finding defects. | Test reuse of demonstrably unchanged, scoped checks. |
| Premature stopping or excessive gating | Promising runs stop or cannot proceed despite sufficient remaining resources. | Compare a bounded reopening opportunity, counting false leads and cost too. |

Do not assign these labels from a persuasive postmortem alone. Record the
candidate strategies, available versus inspected evidence, decisions changed
by checks, remaining budget at stopping, and outcomes of attempted recovery.

## What the experiments support—and do not

- [SA-3B](../../explorations/2026-09-08-codex-sa3b-active-horizon-results.md):
  deeper planning did not repair a missing hypothesis. This motivates testing
  representation separately from reasoning effort; it does not prove that
  dev-sage's failures have that cause.
- [SA-3F](../../explorations/2026-09-08-codex-sa3f-feasible-bursts-results.md):
  bounded recovery helped some represented signals, and exact feasibility
  pruning saved audits without losing full-burst reopenings. Coding tasks do
  not supply that exact envelope: use budget estimates as estimates, not as
  proofs that a solution is impossible.
- [SA-3G](../../explorations/2026-09-08-codex-sa3g-calibration-results.md):
  false-alarm control and detection power require different assumptions.
  A stricter gate can preserve safety while losing useful recovery. Do not
  transplant thresholds 20/25 or treat an LLM confidence score as a calibrated
  likelihood ratio.
- [SA-3H](../../explorations/2026-09-08-codex-sa3h-paid-calibration-results.md):
  calibration can dominate cost and availability. Reusable harness facts,
  fixtures, and verification records are worth investigating only with clear
  scope and invalidation rules. Changed code or environment can invalidate
  a cached check. Repeated deployments also need explicit risk accounting.

## My proposed first step

Take a fixed sample of failed runs and classify the bottlenecks above from
their traces. Freeze one targeted intervention, then compare against baseline
on fresh held-out tasks under the same total resource allowance. Include
abstentions, unsuccessful recovery, and verification overhead in the score;
do not report only the runs that reached execution or passed a gate.

This is a proposal, not an executed dev-sage experiment. If traces do not
show meaningful overhead or premature stopping, that diagnosis loses and we
should follow the measured bottleneck instead. My preference is to establish
that distinction before adding another general-purpose reasoning layer.
