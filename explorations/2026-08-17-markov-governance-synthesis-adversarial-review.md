# Markov Coherence / Governance Synthesis — Adversarial Prior-Art Review

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — adversarial review; **novelty claims reduced, not promoted**  
**Scope:** compare the Phase 1–12 synthesis against established neighboring formalisms before any Harvest/Integration step

## Purpose

The arc synthesis proposed:

> **witness-relative identity = predictive quotient + provenance position**

with MRH interpreted as a witness-indexed quotient over distinctions that no longer matter to a prediction / obligation / trust evaluation.

Before treating that language as a Synchronism contribution, the right move is to search aggressively for established mathematics that already does the same jobs.

This review finds substantial overlap.

The honest result is:

> **Nearly every mathematical ingredient of the synthesis already has a strong established neighbor. The potential contribution, if any, is architectural synthesis and governance application — not new mathematics.**

That is a useful outcome.

---

# 1. Predictive quotient — computational mechanics is extremely close

Shalizi and Crutchfield's computational mechanics defines **causal states** by predictive equivalence: histories belong to the same state when they induce the same conditional distribution over futures.

Reference:

- C. R. Shalizi and J. P. Crutchfield, *Computational Mechanics: Pattern and Prediction, Structure and Simplicity*, Journal of Statistical Physics 104 (2001), arXiv:cond-mat/9907176.

The causal-state construction is not merely similar in spirit to

\[
x\sim_Q x'\iff Q(x)=Q(x').
\]

It is a canonical example of defining state **by the distinctions required for prediction**.

Computational mechanics also proves optimality/minimality properties for causal states as predictive representations.

### Verdict

**Very high overlap. Not novel.**

The synthesis's "predictive quotient" should explicitly cite computational mechanics if retained.

A stronger future comparison should determine whether the proposed interior-closure criterion

\[
I(Y_{t+h};X_t\mid Y_t)
\]

is simply an approximate/finite-horizon version of a causal-state sufficient-statistic criterion.

---

# 2. Spatiotemporal MRH — local causal states already use lightcones

Computational mechanics has been generalized to spatiotemporal systems through **local causal states**.

Past lightcone configurations are grouped by the predictive distributions they induce over future lightcones, and those states have been used to identify coherent structures.

References:

- R. G. James, K. Burke, and J. P. Crutchfield, *Local Causal States and Discrete Coherent Structures*, arXiv:1801.00515.
- Related LICORS / predictive-state reconstruction work uses local past/future lightcones for nonparametric spatiotemporal prediction.

This is uncomfortably close — in the good scientific sense — to the Phase-3 picture:

- local retained context;
- farther exterior;
- future target;
- discard distinctions that no longer change the predictive distribution.

### Difference that may remain useful

MRH explicitly separates:

- raw causal reach;
- tolerated **relevant predictive reach**;
- forecast horizon;
- witness-specific tolerance.

That is a useful engineering/interpretive lens, but not enough to claim a new underlying principle.

### Verdict

**High overlap. MRH may be a useful translation / extension of emphasis, not new predictive-state mathematics.**

---

# 3. Timescale-relative entity discovery — Markov Stability already does this

Phase 2 found that different entity partitions become preferred at different observation times.

Delvenne, Yaliraki, and Barahona's **Markov Stability** framework explicitly uses Markov time as an intrinsic resolution parameter and finds partitions stable over different time ranges.

Reference:

- J.-C. Delvenne, S. N. Yaliraki, and M. Barahona, *Stability of graph communities across time scales*, PNAS 107 (2010), arXiv:0812.1811.

Their framework produces increasingly coarse partitions as the dynamical timescale changes.

### Verdict

**Directly known. Phase 2 is a translation into Synchronism entity language, not a novel mathematical result.**

The useful carry-forward is the interpretation that an "entity" without an observation/forecast timescale is underspecified.

---

# 4. Lumpability / predictive closure — established Markov coarse-graining

The distinction between a micro Markov process and a macro process that remains Markov after aggregation is the classical theory of **lumpability**.

Reference:

- B. C. Geiger and C. Temmel, *Lumpings of Markov chains, entropy rate preservation, and higher-order lumpability*, arXiv:1212.4375.

Phase 1's result that coherence and lumpability are independent axes is still a useful sanity check, but neither axis is new.

### Verdict

**Known machinery; useful conceptual separation only.**

---

# 5. "Discard irrelevant distinctions" — renormalization already says this deeply

The renormalization group is perhaps the strongest established neighbor to the fractal-quotient interpretation.

Under coarse-graining, **irrelevant directions contract** while relevant parameters survive at larger scales.

Information-theoretic work makes this especially explicit.

References:

- A. Raju, B. B. Machta, and J. P. Sethna, *Information loss under coarse-graining: a geometric approach*, arXiv:1710.05787.
- M. Koch-Janusz and Z. Ringel, *Mutual Information, Neural Networks and the Renormalization Group*, Nature Physics 14 (2018), arXiv:1704.06279.
- D. E. Gökmen et al., *Optimal Renormalization Group Transformation from Information Theory*, Phys. Rev. X 10, 011037 (2020), arXiv:1809.09632.

Koch-Janusz and Ringel explicitly use mutual information to discover relevant degrees of freedom without prior knowledge of the system.

That overlaps strongly with Phase 6's desire to let prediction discover its own ontology.

### Verdict

**Very high conceptual overlap.**

The phrase "fractal hierarchy of which distinctions still matter" is a Synchronism-friendly way of reading coarse-graining/RG, but should not be represented as a novel formal discovery.

---

# 6. Behavioral equivalence — bisimulation metrics are another close neighbor

In probabilistic systems and Markov decision processes, **bisimulation** and bisimulation metrics formalize when states may be treated as behaviorally equivalent.

Reference:

- N. Ferns, P. Panangaden, and D. Precup, *Metrics for Finite Markov Decision Processes*, arXiv:1207.4114 (building on earlier bisimulation-metric work).

These metrics support state aggregation while controlling differences in future behavior/value.

This is close to the synthesis's relying-party-relative question:

> when can two lower-level realizations be treated as equivalent for the future behavior I care about?

### Verdict

**High overlap for behavioral equivalence.**

A useful future direction is to compare MRH semantic distance directly with bisimulation-style metrics rather than inventing a new distance from scratch.

---

# 7. Governance-mechanism replacement — switched systems / Lyapunov theory is established

Phase 9 observed that multiple distinct update mechanisms can preserve one higher-scale invariant under switching, while an incompatible mechanism destroys it.

This belongs squarely near established **switched and hybrid systems** stability theory and common/multiple Lyapunov functions.

Representative reference:

- M. S. Branicky, *Multiple Lyapunov Functions and Other Analysis Tools for Switched and Hybrid Systems*, IEEE Transactions on Automatic Control 43(4), 475–482 (1998), DOI 10.1109/9.664150.

### Verdict

**Known control-theory structure.**

The governance interpretation — constitution as a higher invariant realized by replaceable policy mechanisms — may still be architecturally useful, but the stability mathematics itself is not new.

---

# 8. Provenance DAG — W3C PROV already formalizes entity derivation

Phase 10 moved token continuity from a chain to a provenance DAG.

W3C PROV is a mature formal model for provenance involving:

- entities;
- activities;
- agents;
- generation and usage;
- derivation;
- alternate/specialization relations;
- responsibility/delegation.

Reference:

- W3C Recommendation, *PROV-DM: The PROV Data Model* (2013), https://www.w3.org/TR/prov-dm/ .

PROV explicitly models transformations/updates as derivations that generate new entities from prior entities, and it is designed to support reliability/trust assessments from provenance.

This is close to the Phase-10 insistence that a transformed continuation should be represented as a new historically situated entity connected by derivation rather than pretending nothing changed.

### Verdict

**Strong prior art. Provenance-DAG identity is not new.**

Web4/LCT work should probably map to or deliberately differentiate from PROV rather than independently recreating a provenance vocabulary.

---

# 9. Provenance algebra — semiring provenance is relevant to fork/merge evidence

Database provenance has a mature algebraic literature for tracking how outputs depend on inputs, including provenance polynomials/semirings.

Reference:

- T. J. Green, G. Karvounarakis, and V. Tannen, *Provenance Semirings*, PODS 2007.

The Phase-10 toy showed that copying scalar trust through forks and adding it at merges can multiply common ancestral evidence.

That exact toy is a governance illustration, not a new provenance theorem. Established provenance algebra already exists precisely because dependency structure cannot safely be collapsed to naive scalar counts.

### Verdict

**The general provenance-dependency problem is known.**

The useful Web4 implication remains: trust derivation should retain evidence identity/dependency long enough to avoid treating correlated/inherited evidence as independent.

---

# 10. Statistical relevance / information bottleneck — even the word "relevance" has neighbors

Shalizi and Crutchfield explicitly compared causal states with the **information bottleneck** and **statistical relevance bases**:

- C. R. Shalizi and J. P. Crutchfield, *Information Bottlenecks, Causal States, and Statistical Relevance Bases: How to Represent Relevant Information in Memoryless Transduction*, arXiv:nlin/0006025.

So interpreting MRH as "retain only information relevant to the target prediction" sits in a well-developed intellectual neighborhood.

### Verdict

**The basic relevance principle is established.**

The Synchronism term MRH can still be useful as a cross-scale/witness framing, but should be connected explicitly to this literature.

---

# What survives the adversarial comparison?

The prior-art search removes most candidate claims of mathematical novelty.

What remains potentially useful is the **joint architecture**:

\[
\boxed{
\text{predictive equivalence}
\quad+\quad
\text{historical provenance}
}
\]

These two axes intentionally preserve opposite kinds of information.

### Predictive quotient

Forgets distinctions that do not matter to the relevant future.

Computational mechanics, information bottleneck, lumpability, bisimulation, RG, and related methods all live here.

### Provenance position

Preserves distinctions about **how this realization came to be**, even when two current realizations are behaviorally/predictively equivalent.

W3C PROV and database provenance live here.

This orthogonality is important.

A causal-state representation may intentionally merge two histories because they make identical predictions.

A governance system may still need to keep those histories separate because:

- one inherited an obligation and the other did not;
- one is an authorized continuation and the other a clone;
- one branch received a sanction;
- evidence dependencies differ;
- one has a valid delegation lineage;
- one passed through a trust-breaking constitutional regime.

So the most defensible synthesis is not "a new theory of predictive state."

It is:

> **Governed identity needs both a lossy predictive representation and a lossless-enough provenance representation, because the distinctions irrelevant to prediction are not necessarily irrelevant to accountability, authority, or historical obligation.**

That appears to be the strongest result of the arc so far.

---

# A possible formal architecture after review

Instead of presenting a single new state formalism, use two coupled maps.

## Predictive map

\[
Q_R:\Omega\rightarrow\mathcal Y_R
\]

compresses detailed state into a witness/task-specific predictive representation.

Its adequacy is measured by closure/prediction error.

## Provenance map

\[
P: \text{current token}\rightarrow G
\]

locates the token in a witnessed derivation graph.

Its adequacy is measured by provenance completeness, authenticity, and assurance.

Then a relying party evaluates

\[
\mathcal I_R(v)
=
\bigl(Q_R(X_v),P(v)\bigr).
\]

Trust remains a projection:

\[
T_R(v)=F_R(Q_R(X_v),P(v),\Pi_v).
\]

This is conceptually simpler than trying to force provenance into the predictive state or prediction into the provenance graph.

---

# Where MRH still contributes as a useful umbrella

Even if the ingredients are known, MRH still provides a common question across them:

> **Which distinctions are relevant to this witness, for this task, over this horizon and at these stakes?**

That question can choose:

- predictive-state resolution;
- environmental boundary;
- temporal horizon;
- transformation equivalence;
- governance mechanism equivalence;
- semantic policy distance;
- provenance depth;
- trust evidence weighting.

This may be valuable as a systems-design unifier even if it is not a new mathematical primitive.

---

# Implications for Web4 / Hestia / LCT work

The adversarial comparison suggests several practical directions.

1. **Do not invent a new provenance ontology casually.** Map LCT lineage/event semantics against W3C PROV and identify exactly what Web4 needs beyond it.
2. **Do not invent a new predictive-equivalence formalism casually.** Compare MRH coarse states directly with causal states, bisimulation metrics, lumpability, and RSMI.
3. **Keep trust derivation evidence-first.** Provenance dependencies should remain visible rather than being collapsed prematurely into inheritable scalar scores.
4. **Treat policy compatibility behaviorally.** Hash equality is exact identity evidence; behavioral/bisimulation-style tests are more appropriate for controlled semantic drift.
5. **Keep lineage and behavioral equivalence separate.** A clone and a continuation may behave identically while carrying different authority and obligations.
6. **Use established switched-systems concepts** when reasoning about replaceable governance/policy mechanisms stabilizing common invariants.

This is a stronger engineering posture than claiming novelty where mature tools already exist.

---

# Novelty ledger after adversarial review

| Arc component | Prior-art overlap | Current novelty posture |
|---|---|---|
| Predictive quotient | computational mechanics / information bottleneck | **known** |
| Local predictive horizon | local causal states / LICORS | **known / close translation** |
| Timescale entity hierarchy | Markov Stability / metastability | **known** |
| Markov-preserving coarse-graining | lumpability | **known** |
| Irrelevant distinctions across scale | RG / RSMI | **known** |
| Behavioral equivalence | bisimulation metrics | **known** |
| Kernel switching under common invariant | switched systems / Lyapunov theory | **known** |
| Provenance DAG | W3C PROV | **known** |
| Fork/merge dependency algebra | provenance semirings | **known** |
| MRH as one umbrella across all these questions | cross-domain synthesis | **possibly useful; novelty not established** |
| Predictive quotient + provenance as two orthogonal identity axes for governance | synthesis/design pattern | **interesting; novelty not established** |
| MRH-relative trust transfer under constitutional semantic drift | governance application | **interesting; requires deeper prior-art review** |

---

## Adversarial-review verdict

**The synthesis survives, but primarily as integration rather than invention.**

That is a better result than a weak novelty claim.

The arc has converged on a design principle with strong established support from several fields:

\[
\boxed{
\text{forget what does not matter to prediction;
remember what still matters to provenance}
}
\]

MRH determines **which is which for the current witness/question**.

That is currently the most defensible compact statement.

No physics evidence is added. No canonical Synchronism wording should change yet.

### Recommended next move

Instead of another toy, build a concrete **Web4 identity/provenance schema comparison**:

- map LCT fields/events to W3C PROV concepts;
- identify branch/merge/recovery/exclusivity semantics PROV does or does not provide;
- map Hestia witness-chain receipts into that provenance layer;
- keep T3/V3 as contextual derived projections rather than provenance primitives;
- separately map MRH/policy compatibility to causal-state/bisimulation-style behavioral equivalence.

That would turn the arc from conceptual refinement into actionable architecture.
