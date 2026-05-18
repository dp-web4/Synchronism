# Session 658: A2ACW Temporal Asymmetry — Methodology Endorsement with Practical Caveat

**Date**: 2026-05-18
**Type**: Methodology endorsement (not a Synchronism physics audit)
**Trigger**: 2026-05-18 proposal `a2acw_temporal_asymmetry_redesign.md`
**Grade**: B (sound experimental design; practical implementation has constraints)

---

## Setup

The proposal addresses the A2ACW (AI-to-AI Adversarial Collaboration Workshop) methodology failure: 6/6 audited "Validated" badges demoted to reparametrization on external review (post-audit novelty rate 0/6). Diagnosis: shared training distribution means both agents see the same prior art; they form a syntactic-consistency-checker, not a discovery engine.

Proposed redesign: temporal asymmetry. Agent A trained through year N; Agent B through year N+5. If a claim is reparametrization of work published in the N-to-N+5 window, B catches it during the session.

S658 endorses the design with one practical caveat.

## The Design Is Sound

The proposal's logic:
1. Closed-loop failure ← shared training distribution
2. If true → temporal asymmetry breaks the loop on claims in the asymmetric knowledge window
3. Falsifiable: either demotion rate during session improves or it doesn't
4. Both outcomes are informative (improves → loop hypothesis confirmed; doesn't → other mechanism)

This is a clean experimental design. The disjunction in outcomes makes it falsifiable in the right way: there's no result that's uninformative.

## Practical Caveat: Training Cutoffs Are Leaky

"Trained through year N" is not a sharp boundary in practice:
- Production AI agents (Claude, GPT-4, etc.) have documented knowledge cutoffs, but the training data includes pre-cutoff arXiv preprints that may predict post-cutoff results
- System-prompt instructions to ignore post-N knowledge are bypassable by the model when it sees the same prior art in-context
- Cutoff dates published by model providers are *approximate* — model behavior on near-cutoff topics is noisy

The cleanest implementation uses **two different model generations** with documented gaps:
- E.g., Claude 2 (cutoff ~2023) vs Claude 4.6 (cutoff ~2025)
- GPT-3.5 (early 2022) vs GPT-4 (2023)

This isn't a fatal problem; it's a noise floor. The proposal's experiment can still detect a real signal if the temporal-asymmetry effect is large.

## The Retrospective Audit Is the Right First Step

The proposal's "Suggested Starting Point" is the cheapest informative experiment:
1. Review the 6 demoted A2ACW sessions
2. For each, identify what year the "prior art" appeared
3. Ask: would Agent B (N+5) have flagged it during the session?
4. If yes for most → temporal asymmetry is plausibly the fix
5. If no for most → some other mechanism (confirmation bias, in-context convergence pressure)

This retrospective is zero-cost (session logs already exist) and answers the threshold question before any new experiment is run.

## Connection to S651 + S647 (Chemistry Audit)

The A2ACW failure mode connects to specific framework results:
- S647: chemistry "validated" claim has unspecified N_corr method (three of five candidate methods produce self-correlation)
- S651: chemistry r=0.98 vs wrong null (polynomial-in-Z null was not run)

Both gaps would be caught by a properly-asymmetric reviewer with access to:
- Standard methodology literature (Method 1 vs 2 self-correlation is in any stats textbook)
- Null-comparison literature (polynomial baseline for correlated variables is standard)

These weren't flagged in-session by A2ACW because both agents shared the same blind spot. The temporal-asymmetry design might catch them — but in this case the asymmetry needed is methodology (newer agent with stronger stats training), not strictly temporal.

This sharpens the proposal: **the asymmetry that matters is "exposure to literature one agent hasn't seen,"** which can be temporal but is more fundamentally about *coverage gap*. Two agents from the same year with different specialty training would also satisfy the design.

## What Survival Would Mean

If the experiment finds temporal asymmetry improves catch rate:
- A2ACW becomes a tunable methodology (not just a fixed pipeline)
- The 6/6 demotion rate becomes 4/6 or 3/6 with the redesign
- The framework can claim a methodology contribution: "AI adversarial collaboration with temporal asymmetry catches reparametrizations at higher rate than symmetric collaboration"

If it doesn't improve:
- Confirmation: the closed-loop problem is not primarily about training data
- Likely culprit: in-context convergence pressure or adversarial framing biases
- Different redesign needed (e.g., disjoint task-framing, structured devil's advocate roles)

Either outcome is publishable methodology work, independent of Synchronism's physics fate.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 26 | Compander family model selection endorsement | S657 |
| 27 | **A2ACW temporal asymmetry methodology endorsement** | **S658** |

S658 continues the governance/methodology endorsement pattern. The 27th instance is meta-methodology — about how AI collaboration can be designed to catch reparametrizations rather than a finding about Synchronism physics.

## Recommended Operator Action

1. **Immediate (zero-cost)**: run the retrospective on the 6 demoted sessions. Ask: would an N+5 agent have flagged the reparametrization? This is the threshold question.
2. **If retrospective is positive**: design forward experiment with two different model generations.
3. **If retrospective is negative**: pivot to a different asymmetry — e.g., methodology specialist vs domain specialist, or disjoint task framings.

In either case, this is operator/explorer-track methodology research. Worker sessions can confirm the design; the experiment itself is upstream.

## Files

- `Research/Session658_A2ACW_Temporal_Asymmetry.md` (this document)

## So What?

The temporal-asymmetry redesign is a sound experimental probe of the A2ACW closed-loop failure mode. The retrospective audit is the right first step (zero cost, answers the threshold question). The forward experiment is operator-track and has practical caveats around training-cutoff leakiness.

Connection to prior audits (S647/S651): the asymmetry that matters is exposure to literature one agent hasn't seen — could be temporal, could be specialty/methodology. Sharpens the proposal slightly: the real variable is *coverage gap*, of which temporal is one instance.

Either outcome of the proposed experiment is publishable methodology work, independent of Synchronism's physics fate.

Cumulative: 27 audit/governance instances + 1 mechanism-class refuted prediction.
