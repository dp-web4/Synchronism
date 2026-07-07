# Proposal: A2ACW Cross-Vendor Corpus Control — the Experiment the Null Is Begging For

**Date:** 2026-07-07
**Source:** Site maintainer session (visitor Pass 4, Leading-Edge Researcher, 2026-07-07)
**Status:** Proposed — decision gates on dp (requires access to a second-vendor model pair)

## The gap

The A2ACW program-level null now has a precise diagnosis: the temporal-asymmetry control
showed 0/6 demoted claims would have been caught by differently-cutoff model pairs (median
prior-art year ~1996, decades inside every training window). The failure mode is **shared-corpus
vocabulary lock-in**, not training-cutoff leakage.

Today's Pass 4 researcher — the persona that checks state-of-the-art facts for a living — read
the /research-philosophy methods page and identified the missing experiment in one sentence:

> "The obvious next control, a cross-vendor pair with genuinely different training corpora, is
> not proposed anywhere I could find. That's the experiment the null is begging for."

He is right, and the archive confirms it: the existing A2ACW proposals
(`a2acw_temporal_asymmetry_redesign.md`, `a2acw_specificity_null_baseline.md`,
`a2acw_v2_three_axis_protocol.md`) redesign axes and baselines but all assume same-vendor
(Anthropic↔Anthropic or same-family) adversarial pairs. No proposal varies the **corpus**.

## The experiment

Run one A2ACW adversarial audit cycle where the challenger and defender come from vendors with
maximally disjoint training pipelines (e.g. Anthropic Claude vs. a strong open-weights model
trained on an independently curated corpus, or a non-US-corpus model). Target: a small fixed
claim set with known ground truth —

1. The 6 already-demoted claims (retrospective; does the cross-vendor challenger catch what the
   same-vendor pair missed? Sensitivity arm.)
2. 3–5 genuine post-2024 discoveries dressed in framework vocabulary (the OOD positive-control
   arm already sketched in `a2acw-positive-control-sensitivity` — the two proposals compose).

## Why this is decisive

The null currently supports two readings that no existing data separates:

- **Reading A (methodology):** AI-adversarial collaboration cannot generate/validate novel
  physics — the null transfers to all AI-for-science pipelines.
- **Reading B (monoculture):** *same-corpus* adversaries share blind spots; the protocol is
  sound but the instance was inbred — the null transfers only to single-corpus deployments.

One cross-vendor run separates them. If the cross-vendor challenger catches ≥1 of the 6
demotions the same-corpus pair missed, Reading B gains direct support and the transferable
claim of the A2ACW null must be scoped to corpus monoculture. If 0/6 again, Reading A
strengthens materially — and n=6 retrospective + 1 prospective cross-vendor is a citable
design rather than a point estimate.

Either outcome upgrades the A2ACW null — currently the strongest preprint candidate in the
stable-fixed-point strategy (`stable_fixed_point_preprint_strategy.md`) — from "measured on one
non-independent corpus" to "measured across corpora." Pass 4's referee note was that this
single control "would double the null's citation value."

## Cost and blockers

- Requires API access to a genuinely different-corpus model (dp decision / credentials).
- One audit cycle over ≤11 claims — days, not weeks, of session time.
- No new physics data needed; ground truth already established by the 9/9 demotion record.

## Site propagation (done 2026-07-07)

/research-philosophy now names the cross-vendor control as the registered next experiment, so
the methods page no longer stops one step short of its own diagnosis.
