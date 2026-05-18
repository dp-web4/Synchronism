# Proposal: A2ACW Temporal Asymmetry Redesign

**Filed**: 2026-05-18
**Origin**: Pass 4 researcher in visitor log 2026-05-18
**Status**: Open

---

## The Problem

The current A2ACW (AI-to-AI Adversarial Collaboration Workshop) methodology has a known structural failure:
two AI agents sharing the same training distribution form a closed adversarial loop.
They function as a syntactic consistency checker — they can detect internal contradictions
but cannot generate out-of-distribution novelty.

Evidence: 6 of 6 audited "Validated" badges demoted to Reparametrization on external review.
Post-audit novelty rate: 0/6. This was diagnosed on the site as "in-distribution detection only."

The current proposed remedy — "evaluation by domain experts outside the training loop" — is a wish,
not a protocol. It requires external humans, is not scalable, and has no triggering criteria.

---

## The Proposed Experiment

**Key insight from the 2026-05-18 researcher visitor**: The shared-training-distribution problem
might be addressable within AI-only collaboration by introducing *temporal asymmetry* between agents.

### Design

- **Agent A**: Trained on physics literature through year N
- **Agent B**: Trained on physics literature through year N+5

Agent A proposes claims. Agent B challenges them with access to five more years of literature.
For any claim that was already known in year N, Agent B has no advantage. For any claim that
appeared in the literature between N and N+5, Agent B can flag it as prior art or contradict it.

This is a structural asymmetry that — if effective — would catch reparametrizations of results
published in the N-to-N+5 window without requiring human expert intervention.

### Prediction

If the closed-loop failure mode is *primarily* caused by shared training distribution:
- The temporal asymmetry design should increase the rate at which "Validated" claims
  are flagged as reparametrizations during the session (before human audit).
- The post-session demotion rate (the 6/6 figure) should decrease.

If the failure mode is *not* primarily caused by training distribution:
- No improvement should be observed, pointing to a different structural cause
  (e.g., confirmation bias in adversarial framing, in-context pressure to converge).

### Why This Matters

The A2ACW methodology is the framework's principal generative engine.
If the closed-loop problem can be partially fixed within the AI collaboration itself,
the 47:0 internal-claim:confirmed-prediction ratio becomes a tractable problem.
If not, the architecture requires a fundamentally different design (human-in-the-loop,
or reframing A2ACW as a consistency-checker-only rather than a discovery engine).

---

## Tractability

This is a methodologically tractable experiment:

1. Identify a physics domain with several claims that were published between 2020–2025
   (so a "year N = 2019" agent would miss them, a "year N+5 = 2024" agent would catch them).
2. Run standard A2ACW (both agents same training) — measure demotion rate by human expert.
3. Run temporal-asymmetry A2ACW (Agent A = 2019 cutoff, Agent B = 2024 cutoff) — measure demotion rate.
4. Compare.

The comparison is falsifiable: either the temporal asymmetry increases the in-session catch rate, or it doesn't.

---

## Connection to the Framework

This is not a physics question. It is a research methodology question about whether
the framework's own generative engine can be improved without abandoning AI collaboration.

The framework's most defensible published contribution is the A2ACW negative result:
*"AI adversarial collaboration at scale produces a 0% out-of-distribution discovery rate
(n=6), with self-diagnosed mechanism."*

If this proposal succeeds, the contribution upgrades from "here is why it fails"
to "here is how to partially fix it." If it fails, it strengthens the constraint.

Either way, the investigation is productive.

---

## Suggested Starting Point

- Review A2ACW session logs for the 6 demoted Validated claims.
- For each: identify what year the "prior art" appeared in the literature.
- Determine whether Agent B (N+5) would have flagged it during the session.
- This is a retrospective audit that can inform whether the temporal asymmetry design
  would have made a difference — before running any new experiments.
