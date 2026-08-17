# Markov Coherence / Governance Arc — Synthesis: MRH as Predictive Quotient + Provenance

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — synthesis / candidate formal refinement, **not canonical**  
**Scope:** integrates Phases 1–11; no physics claim; does not modify `FUNDAMENTALS.md`

## Starting point

The canonical Synchronism definition currently says:

> **MRH:** the boundary beyond which correlations become negligible for a given pattern's dynamics; the scale at which a pattern can treat everything else as averaged bulk behavior.

The exploration arc began by asking whether Markov machinery could make that idea operational.

It ended up doing more than that.

Across Phases 1–11, the same structure kept reappearing in physics-style toys and governance-style toys:

> **A higher-scale identity exists when many lower-scale distinctions can be discarded without materially damaging the prediction / obligation / trust question relevant to the witness.**

That suggests interpreting MRH not primarily as a geometric boundary but as a **witness-indexed quotient over irrelevant distinctions**.

This document compresses the arc into one provisional formal object.

---

# 1. Microstate, macrostate, and quotient

Let the detailed system state be

\[
X_t\in\Omega.
\]

Choose a coarse variable

\[
Y_t=Q(X_t),
\]

where

\[
Q:\Omega\rightarrow\mathcal Y.
\]

The map \(Q\) induces an exact equivalence relation

\[
x\sim_Q x'
\iff
Q(x)=Q(x').
\]

Thus a macrostate is literally an equivalence class of microstates:

\[
[x]_Q.
\]

The exploration's core question becomes:

> **When is the witness justified in using that quotient?**

The answer requires dynamics, not symbolism alone.

---

# 2. Interior predictive closure

A proposed coarse variable should discard micro-detail only when that detail contributes little additional predictive value for the macro-variable's future.

A natural criterion is

\[
\epsilon_{int}(Q;h)
=
I(Y_{t+h};X_t\mid Y_t).
\]

If

\[
\epsilon_{int}\le\epsilon^*_{int},
\]

then knowing the exact microstate adds little predictive information beyond the coarse state itself.

This is the **interior resolution horizon** discovered in Phase 5:

> below some resolution, additional implementation detail no longer matters to the macro-identity for the chosen forecast horizon.

So an entity has an informational inside as well as an outside.

---

# 3. Exterior predictive closure / MRH

Split environmental information into a retained boundary/context \(B_t\) and farther exterior \(E_t\).

Then define exterior predictive leakage

\[
\epsilon_{ext}(Q,B;h)
=
I(Y_{t+h};E_t\mid Y_t,B_t).
\]

An MRH boundary is adequate when

\[
\epsilon_{ext}\le\epsilon^*_{ext}.
\]

This formalizes the canonical phrase "beyond which correlations become negligible" as:

> **beyond which information adds less than the tolerated amount to the prediction relevant to this witness.**

Phase 3 demonstrated that this relevant predictive horizon can be strictly smaller than the raw causal horizon.

Something may be causally capable of influencing the target while being too weakly informative to justify explicit representation at the current MRH.

---

# 4. Temporal MRH / time-scale separation

Phases 1 and 2 established that entityhood is also temporal.

Let

- \(\tau_{internal}\) = time for lower-scale implementation detail to mix / be forgotten;
- \(\tau_{identity}\) = time for the macro-identity itself to be lost.

An entity-like scale exists over horizons satisfying approximately

\[
\boxed{
\tau_{internal}\ll h\ll\tau_{identity}
}
\]

while interior and exterior closure errors remain acceptably small.

This gives a precise interpretation of the earlier intuition that MRH is temporal/fractal.

At short \(h\), lower-level details remain relevant.

At intermediate \(h\), they can be quotiented out while a slower collective mode remains predictive.

At still longer \(h\), that collective identity may itself become implementation detail of another scale.

---

# 5. Coherence and emergence are separate axes

Phase 1 showed that a macro partition can be perfectly Markovian while weakly coherent.

So two properties must remain separate.

### Dynamic coherence

\[
\mathcal C
\sim
\frac{\tau_{identity}}{\tau_{internal}}.
\]

Large \(\mathcal C\): internal implementation changes many times before macro-identity is lost.

### Predictive closure

\[
\epsilon_{int}
=I(Y_{t+h};X_t\mid Y_t).
\]

Small \(\epsilon_{int}\): discarded micro-detail contributes little additional prediction.

An emergent entity-like variable needs both persistence and useful closure; neither alone is sufficient.

---

# 6. Relations can be entities

Phase 4 falsified the naive idea that dissonance simply means lack of coherence.

Strong positive and negative coupling can both create long-lived relational modes.

Thus the persistent variable may be neither constituent alone but a relation

\[
R(A,B).
\]

A relation becomes entity-like when it is itself a slow, predictively useful variable.

This is the first major generalization:

\[
\boxed{
\text{entities need not be aggregates; they can be persistent relations}
}
\]

---

# 7. Relations among relations / recursive quotients

Phase 5 then showed that higher-order relations can survive transformations that erase the identities of their lower-level implementation.

The recursive structure is

\[
X
\xrightarrow{Q_1}
Y^{(1)}
\xrightarrow{Q_2}
Y^{(2)}
\xrightarrow{Q_3}
\cdots
\]

where each higher map discards distinctions that have become fast / predictively irrelevant at the next scale.

When the hierarchy is genuinely nested, this gives quotient spaces

\[
\Omega
\rightarrow
\Omega/{\sim_1}
\rightarrow
(\Omega/{\sim_1})/{\sim_2}
\rightarrow\cdots
\]

The **fractal hierarchy is therefore a hierarchy of distinctions that still matter**, not necessarily a hierarchy of geometric object sizes.

This is the strongest formal reading of "fractal MRH" produced by the arc.

Important caveat: arbitrary task-dependent MRHs need not form a perfectly nested chain. Different witnesses/questions can induce overlapping or incompatible coarse-grainings. The nested quotient picture applies where genuine multiscale structure exists.

---

# 8. Identity as invariance under allowed transformation classes

Phases 5–8 successively replaced:

- microstate;
- components;
- particular relations;
- topology.

The higher identity survived as long as those transformations preserved the relevant macro-variable.

So introduce a transformation class

\[
\mathcal T_R
\]

for witness/task \(R\).

The class contains changes the witness is prepared to treat as **internal implementation change** rather than identity destruction.

Examples from the toys:

- state fluctuations;
- component turnover;
- relation turnover;
- topology rewiring.

The macro-variable \(Y\) is admissible only if it remains approximately invariant / predictively closed under the transformations being quotiented out.

This gives the evolving MRH notation

\[
\operatorname{MRH}
(Y,h,\epsilon,\mathcal T_R).
\]

---

# 9. Governance mechanisms can also be quotiented out

Phase 9 replaced the transition rule itself.

Different governance kernels

\[
K_1,K_2,K_3
\]

could all stabilize the same higher constitutional invariant.

Thus even the **dynamics implementing the entity** can become lower-MRH implementation detail.

Write provisionally

\[
K_a\sim_C K_b
\]

when mechanisms \(K_a,K_b\) belong to the same constitutionally admissible family for the relying party's question.

A common stochastic drift / common-Lyapunov-like condition is one possible operational test:

\[
\mathbb E[\Delta V_C\mid X,K]<0
\]

outside an acceptable region for every

\[
K\in\mathcal K_C.
\]

The exact criterion will be domain-specific.

The structural point is:

> **governance identity belongs at the level of preserved constraints, not necessarily at the level of one permanent enforcement algorithm.**

---

# 10. Predictive equivalence is not enough for historical identity

Phase 8 exposed an important distinction.

Two independently created systems may implement the same macro-invariant:

\[
Q(X_A)=Q(X_B).
\]

That establishes **type / relational equivalence**, not historical numerical identity.

So the quotient must be paired with provenance.

This yields the central synthesis:

\[
\boxed{
\text{witness-relative identity}
=
\text{predictive quotient}
+
\text{provenance position}
}
\]

The quotient answers:

> Which lower-level differences are irrelevant to what this entity is at this MRH?

The provenance position answers:

> Which historical continuation of that pattern is this?

Neither answers the other's question.

---

# 11. Provenance is a DAG, not a chain

Phase 10 replaced simple lineage with a directed provenance graph

\[
G=(V,E).
\]

A witnessed transition is an edge.

Reachability

\[
u\preceq v
\]

means \(v\) is a descendant of \(u\).

This is naturally a partial-order / ancestry relation.

It is not an identity equivalence relation because a fork

\[
A\rightarrow B,
\qquad
A\rightarrow C
\]

makes both \(B\) and \(C\) valid descendants without making

\[
B=C.
\]

Therefore separate:

- **type identity** — same relevant quotient / invariant class;
- **lineage relatedness** — connected by ancestry / common history;
- **token identity** — this particular current branch / node in provenance history.

---

# 12. Trust is a projection, not identity

Phase 10 also showed why inherited trust scores are dangerous under fork/merge.

Shared ancestral evidence can be counted repeatedly even though no new evidence has been created.

Therefore the primary object should be evidence/provenance, not an inheritable scalar trust value.

For relying party \(R\), write

\[
T_R(v)
=
F_R
\left(
\operatorname{UniqueEvidence}(\operatorname{Ancestors}(v)),
Q(v),
G_v
\right).
\]

The trust value is a **contextual projection** of:

- unique witnessed evidence;
- current relevant macro-properties;
- provenance/history;
- the relying party's own MRH and stakes.

This is consistent with Web4's broader claim that trust is contextual rather than absolute.

---

# 13. Continuous identity does not imply continuous trust transfer

Phase 11 added constitutional drift.

A descendant may have perfect witnessed lineage while its relevant behavior gradually moves far from the constitution a relying party originally admitted.

Thus trust transfer needs at least:

1. valid lineage;
2. current semantic compatibility;
3. relevant historical path conditions.

Provisionally:

\[
\operatorname{Transferable}_R(A\rightarrow B)
=
\operatorname{ValidLineage}(A,B)
\land
D_R(C_A,C_B)\le\delta_R
\land
H_R(path)\le\eta_R.
\]

This is structural notation, not a proposed universal formula.

The key result is:

\[
\boxed{
\text{identity continuity}\neq\text{trust continuity}
}
\]

Trust must be re-evaluated in the current witness context.

---

# 14. A provisional unified entity object

For a witness / task \(R\), a current entity token at provenance node \(v\) can be represented provisionally as

\[
\mathfrak E_R(v)
=
\Bigl(
Q_R,
Y_v,
h,
\epsilon_{int},
\epsilon_{ext},
\mathcal T_R,
G_v,
\Pi_v
\Bigr),
\]

where:

- \(Q_R\): coarse-graining / quotient map relevant to the witness;
- \(Y_v=Q_R(X_v)\): current macrostate / relational invariant;
- \(h\): forecast / obligation horizon;
- \(\epsilon_{int}\): tolerated lost prediction from discarded internal detail;
- \(\epsilon_{ext}\): tolerated lost prediction from discarded exterior context;
- \(\mathcal T_R\): transformations treated as identity-preserving implementation change;
- \(G_v\): provenance-DAG position / lineage history;
- \(\Pi_v\): witnessed evidence/provenance references.

Trust, authority, admissibility, and other evaluations are then functions of this object in a relying party's MRH rather than intrinsic scalar properties of the entity.

---

# 15. Minimal entity-like criterion

The arc suggests a compact candidate test for an entity-like macro-variable \(Y\) at witness scale \(R\):

### Time-scale separation

\[
\tau_{internal}\ll h\ll\tau_{identity}.
\]

### Interior closure

\[
I(Y_{t+h};X_t\mid Y_t)\le\epsilon_{int}.
\]

### Exterior closure

\[
I(Y_{t+h};E_t\mid Y_t,B_t)\le\epsilon_{ext}.
\]

### Transformation robustness

Allowed transformations

\[
T\in\mathcal T_R
\]

leave \(Y\)'s relevant predictive behavior inside tolerance.

### Token provenance

The current realization occupies a witnessed node/path in the provenance DAG if historical identity matters to the question.

This is not a theorem and not yet a canonical Synchronism definition.

It is the cleanest synthesis of what survived the arc.

---

# 16. The conceptual compression

The original intuitive picture was approximately:

\[
\text{small patterns form larger patterns, recursively}.
\]

The refined picture is:

\[
\boxed{
\text{fast distinctions}
\rightarrow
\text{slow predictive relations}
\rightarrow
\text{quotient out irrelevant implementation}
\rightarrow
\text{repeat}
}
\]

while provenance tracks **which historical token** realizes those patterns through time.

Or in words:

> **An emergent entity is a slow predictive equivalence class over many changing implementations; a historical entity is that class plus witnessed lineage.**

This captures component replacement, relational replacement, governance implementation replacement, and branching history without requiring literal material persistence.

---

# 17. Why this feels like a real refinement of MRH

The canonical MRH definition emphasizes the outer boundary beyond which correlation is negligible.

The arc suggests four simultaneous horizons:

### Exterior horizon

Which environment distinctions can be ignored?

### Interior horizon

Which implementation distinctions can be ignored?

### Temporal horizon

Which fast changes can be averaged away before the identity itself changes?

### Transformational horizon

Which classes of state/component/relation/mechanism change count as internal evolution rather than a new entity?

These are not four unrelated concepts.

They are all instances of one question:

> **Which distinctions remain relevant to the witness's current prediction / obligation / trust evaluation?**

That is a direct formal expansion of the word **Relevancy** in MRH.

---

# 18. Governance / Web4 consequence

The synthesis maps unusually well onto governance architecture.

A durable Web4 society should not be identified solely with:

- members;
- machines;
- vendors;
- topology;
- one policy engine;
- one policy hash forever;
- one scalar trust score.

Instead, a relying party can reason over:

1. the current governance invariants relevant to it;
2. the admissible transformation class at its assurance level;
3. the witnessed lineage connecting current state to prior trust anchors;
4. unique evidence receipts rather than inherited score claims;
5. semantic drift relative to the relying party's original admission state.

This turns "same governed entity" from a metaphysical assertion into an evidence-backed, MRH-relative evaluation.

---

# 19. What remains unresolved

This synthesis creates several research questions rather than closing them.

1. **Discovery:** Can useful quotient maps \(Q\) be learned robustly in large nonlinear systems without candidate observables?
2. **Approximate quotients:** What is the best formalism when closure is approximate and task-dependent rather than an exact partition?
3. **Non-nested MRHs:** How should incompatible/overlapping coarse-grainings across different witnesses be related?
4. **Transformation discovery:** Can \(\mathcal T_R\) itself be learned rather than specified?
5. **Semantic governance distance:** How should a relying party compare policy behavior without exhaustive enumeration?
6. **Fork/merge exclusivity:** What cryptographic or distributed-systems evidence is sufficient to distinguish continuation from split brain at different assurance levels?
7. **Trust dependence:** How should evidence correlations beyond exact receipt duplication be represented?
8. **Relations among transformations:** Can higher-order invariants exist not only over relations but over sequences/classes of transformations themselves?

The last question is especially interesting because it may be the next genuinely new layer rather than another restatement of the same quotient idea.

---

## Synthesis verdict

**The arc has reached a coherent formal candidate worth critique, but not canonization.**

What appears robust across the toys:

- coherence is time-scale separation, not mere persistence;
- emergence and predictive closure are distinct;
- entity scale is horizon-relative;
- relevant predictive reach differs from causal reach;
- relations and relations-among-relations can be entity-like slow variables;
- component, relation, topology, and even governance-mechanism replacement can be identity-preserving at a higher MRH;
- historical token identity requires provenance in addition to structural/predictive equivalence;
- provenance branches into a DAG;
- trust must derive from evidence/provenance rather than inherit naively through branches;
- continuous identity does not imply continuous trust transfer after semantic drift;
- MRH can be read as a witness-indexed quotient over distinctions that no longer matter.

No part of this changes the physics evidence ledger.

The appropriate next move is adversarial review of the synthesis — especially comparison against existing coarse-graining, bisimulation, causal-state, renormalization, switched-systems, provenance, and process-identity formalisms — before deciding whether any wording should migrate into the forum's Harvest/Integration path or into canonical fundamentals.
