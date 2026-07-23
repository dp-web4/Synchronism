# Proposal: Prior-Art Audit of the Locality No-Go Against Milgrom's Surface-Density Literature

**Date**: 2026-07-23
**Origin**: Maintainer WAKE, from site visitor Pass 4 (external-researcher persona, 2026-07-23)
**Status**: PROPOSED — literature execution, no new computation required
**Gates**: the preprint / citable-artifact decision (`stable_fixed_point_preprint_strategy.md`, gates on dp)

## The problem

The local-density no-go is the **one artifact the program currently stakes novelty on**
(/for-researchers presents it as the citable contribution; the preprint strategy lists it as
transferable null #1). Its provenance chain is already partially honest: the 2026-06-09 novelty
audit identified it as a quantified **instance of Milgrom's non-locality theorem**
(astro-ph/0510117), and the site credits "the core obstruction is *not* ours."

But today's external-researcher review identified a specific unexecuted step: **nobody has
searched Milgrom's surface-density formulation papers** (and the related σ-based MOND
literature) **for a prior *quantified* instance of the same statement.** The no-go's marketable
content is not the obstruction (Milgrom's) but the quantification: the ~1.7 dex cross-system
offset, the 10⁴–10⁶× cluster density mismatch, the ρ_crit ∝ V⁻² sign inversion, and the
switch/force-scale conflation diagnosis. If adjacent quantified statements exist in Milgrom's
surface-density discussions (e.g. the Σ† ≈ a₀/G transition surface density, MOND-Σ variants,
or the 2005 non-locality paper's own worked examples), the artifact shrinks from "quantified
corollary, first quantification ours" to "restatement with different units" — and the preprint
strategy should know which *before* dp commits effort.

The reviewer's framing cuts both ways, and honestly: "One paragraph situating the no-go against
Milgrom's σ-based papers **strengthens rather than weakens the claim if it survives**."

## Registered question (fixed before the search)

> Does any published statement in the Milgrom surface-density / modified-inertia literature
> quantify the failure of *local volumetric density* ρ(r) as the MOND organizing variable —
> specifically (a) a cross-system offset magnitude, (b) a cluster-scale density mismatch, or
> (c) a sign/scaling constraint equivalent to ρ_crit ∝ V⁻² — as opposed to stating the
> non-locality obstruction qualitatively?

**Verdict rule, pre-fixed both ways:**
- **Prior quantification found** → the no-go's novelty claim is demoted to "independent
  re-derivation"; /for-researchers reframes the artifact as a *synthesis* (the three-legged
  joint statement may still be novel as a package); preprint value drops accordingly and the
  strategy doc gets a dated note.
- **No prior quantification found** (after checking the named corpus below) → the no-go gains a
  "prior-art audited 2026-07-XX" line naming what was searched — which is precisely what a
  referee will demand — and the artifact's citability *strengthens*.

## Corpus to check (minimum)

1. Milgrom 1983a-c (the founding trilogy — surface-density scalings are already in 1983b).
2. Milgrom, "MOND as modified inertia," and the 2005 non-locality paper (astro-ph/0510117) —
   including its worked examples, not just the theorem.
3. Milgrom's Σ† / surface-density-threshold discussions (e.g. "The central surface density of
   dark halos" line of work, ~2009–2016; Brada & Milgrom on disc stability).
4. The MOND reviews (Famaey & McGaugh 2012 §6; Milgrom's Scholarpedia entry) for any
   density-vs-acceleration organizing-variable discussion.
5. Secondary: Banik & Zhao 2022 review; any explicit "why not volume density?" passages.

NotebookLM is the right tool if auth is restored (operator: re-login on CBP — expired since
07-22); otherwise WebFetch on arXiv abstracts + targeted PDF reads.

## Why this is a research-direction item, not a site fix

The site can only inherit whatever verdict this search produces. The decision it informs —
whether the program's single claimed citable artifact survives contact with its own prior-art
obligation — lives in the research repo and gates a dp-level strategy call. It is also the
program's own methodology applied to itself: the A2ACW lesson ("every audited 'contribution'
died on prior-art") says the prior on "our quantification is novel" should be LOW until the
walk is executed.

## Deliverable

`explorations/2026-07-XX-locality-nogo-prior-art-audit.md` with per-paper findings, plus the
one-paragraph situating text for /for-researchers (either framing, per the verdict rule).
