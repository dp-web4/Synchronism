# Governance Arc — Phase 16: Evaluator Replacement Is Another Component-Turnover Problem

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — evaluator/proof-horizon formalism result  
**Code:** [`simulations/governance_phase16_evaluator_replacement.py`](../simulations/governance_phase16_evaluator_replacement.py)  
**Scope:** governance/evidence architecture toy; not a protocol standard, security proof, or physics claim

## Why this phase exists

Phase 15 ended with a recursion:

> once dependency declarations become executable contracts, the **evaluator itself** becomes part of the proof basis.

That means evaluator replacement is structurally the same kind of question that started this arc:

- components may be replaced while a larger identity persists;
- governance mechanisms may be replaced while a constitutional invariant persists;
- now evaluators may be replaced while a relying party hopes the old semantic decision still persists.

The naive question is:

> **Does the new evaluator give the same answers?**

That is necessary in some MRH, but it is not sufficient.

A second question is independent:

> **Does the new evaluator require the same evidence horizon?**

Those can disagree.

---

# 1. Three kinds of evaluator replacement

Use a baseline evaluator

\[
F_0(x)=0.6\,base+0.4\,witness
\]

with evidence dependency set

\[
D_0=\{base,witness\}.
\]

Compare three replacements.

## Exact refactor

\[
F_1(x)=\frac{3\,base+2\,witness}{5}.
\]

It is algebraically identical and consumes the same evidence.

## MRH-specific simplification

\[
F_2(x)=base.
\]

This is **not** globally equivalent to \(F_0\), but inside a restricted MRH in which

\[
witness=base
\]

it is exact while consuming only

\[
D_2=\{base\}.
\]

Thus a replacement can become valid only after a witness/task-specific quotient has made one lower-level distinction redundant.

## Hidden dependency expansion

\[
F_3(x)=F_0(x)+\text{rare conditional override}
\]

with dependency set

\[
D_3=\{base,witness,mode,override\}.
\]

On the currently observed region where `mode=0`, it is behaviorally identical to \(F_0\).

But it has widened the evidence horizon and contains a latent semantic branch outside that observed region.

---

# 2. Result — behavioral equivalence and proof-horizon equivalence are independent

The exact exhaustive toy gives:

| replacement | correlated MRH | observed MRH | all states | evidence relation |
|---|---|---|---|---|
| exact refactor \(F_1\) | exact | exact | exact | same dependencies |
| simplification \(F_2\) | exact | differs by up to 0.4 | differs by up to 0.4 | dependencies contract |
| hidden expansion \(F_3\) | exact | exact | differs by up to 0.15 | dependencies expand |

So the categories are genuinely different.

### Drop-in compatible

\[
F_a\equiv_R F_b
\quad\text{and}\quad
D_b=D_a.
\]

### Compatible with proof-horizon contraction

\[
F_a\equiv_R F_b
\quad\text{and}\quad
D_b\subset D_a.
\]

The new evaluator preserves the relevant result while requiring fewer distinctions.

### Behaviorally compatible but proof-horizon expanding

\[
F_a\equiv_R F_b
\quad\text{but}\quad
D_b\supset D_a.
\]

The current answers may match perfectly, yet an old proof bundle is no longer complete for the new evaluator.

### Semantic change

\[
F_a\not\equiv_R F_b.
\]

Old reliance must be reconsidered for that MRH.

This gives a compact result:

\[
\boxed{
\text{same answer}\not\Rightarrow\text{same proof basis}
}
\]

---

# 3. Old proof bundles cannot simply be relabeled

Suppose an old proof bundle contains exactly

\[
B_0=\{base,witness\}.
\]

Then \(F_1\) can be recomputed from it.

Inside the correlated MRH, \(F_2\) can also be recomputed from it because its new dependency set is a subset.

But \(F_3\) cannot legitimately be recomputed from the old bundle because it now needs

\[
mode,override.
\]

Even though its current observed answers are identical to the old evaluator, the old evidence basis is incomplete under the new plan.

Therefore:

> **proof transfer across evaluator replacement should mean re-derivation under the new evaluator, not renaming an old result as though the plan had never changed.**

This mirrors the earlier lineage result:

- evidence may migrate;
- provenance may establish authorized continuity;
- the new claim should still be derived under the new plan.

---

# 4. Evaluator equivalence is MRH-relative

The simplification \(F_2\) is the cleanest example.

Globally it is wrong relative to \(F_0\).

Inside the restricted state set

\[
\Omega_R=\{x:witness=base\},
\]

it is exact.

So evaluator compatibility should not be written merely as

\[
F_a\equiv F_b.
\]

It is more honestly something like

\[
\boxed{
F_a\equiv_{R,\Omega_R,\epsilon}F_b
}
\]

meaning that the two evaluators are behaviorally equivalent to tolerance \(\epsilon\) over the states/questions relevant to relying party \(R\).

That is directly analogous to the earlier result that entity identity itself is timescale/MRH-relative.

---

# 5. A proof-horizon contraction is a genuine higher-level gain

The \(F_2\) case is especially interesting.

If a higher-level invariant makes `witness` perfectly determined by `base` inside the accepted domain, retaining both variables is redundant for this decision.

The replacement evaluator can quotient out one distinction:

\[
D_0=\{base,witness\}
\rightarrow
D_2=\{base\}.
\]

This is the evaluator-level version of component replacement / relational coarse-graining:

> a distinction that mattered at the lower scale becomes implementation detail at the higher one.

So a better evaluator is not necessarily one with more inputs or more complexity.

It may be one that proves the same relevant proposition while depending on **less** evidence.

That gives a useful engineering objective:

\[
\boxed{
\text{preserve required behavior while minimizing justified dependency footprint}
}
\]

subject to assurance and safety constraints.

This is closely related to sufficient-statistic / minimal-representation ideas; no mathematical novelty is claimed.

---

# 6. Current-corpus equivalence is weak evidence

The hidden-expansion evaluator \(F_3\) is exactly identical to \(F_0\) on every state in the currently observed `mode=0` region.

A regression suite containing only that region would pronounce it perfect.

Outside the region, it differs on three of the sixteen exhaustive states, by as much as 0.15.

If the rare `mode=1` branch occurs with probability \(p_m\), and the other toy variables are uniform, the actual divergence probability is

\[
p_{diverge}=0.375\,p_m.
\]

Random testing requires approximately:

| rare-mode probability | divergence probability | trials for 95% detection | trials for 99% detection |
|---:|---:|---:|---:|
| 10% | 3.75% | 79 | 121 |
| 1% | 0.375% | 798 | 1,226 |
| 0.1% | 0.0375% | 7,988 | 12,279 |
| 0.01% | 0.00375% | 79,885 | 122,803 |

This restates the Phase-15 rare-branch lesson at the evaluator-replacement level.

A replacement cannot earn strong equivalence merely by matching a finite observed corpus.

---

# 7. Evaluator replacement has three orthogonal axes

The arc now suggests separating at least three questions.

## Semantic compatibility

Does the new evaluator preserve the relevant decision relation?

\[
d_R(F_a,F_b)\le\epsilon_R.
\]

## Evidence compatibility

Can the new evaluator operate from the old proof horizon?

\[
D_b\subseteq D(B_a).
\]

## Historical/authority continuity

Was the evaluator replacement itself an authorized, witnessed governance transition?

This is provenance, not behavior.

A byte-identical evaluator installed through an unauthorized path can be behaviorally equivalent while lacking governance continuity.

A legitimate upgrade can be historically continuous while materially changing semantics.

Again the axes must not be collapsed.

---

# 8. Candidate replacement relation

For relying party / task \(R\), old evaluator \(F_a\), new evaluator \(F_b\), and old proof bundle \(B_a\), define provisionally:

\[
\operatorname{CompatibleReplacement}_R(F_a\to F_b)
\]

only when:

\[
\operatorname{AuthorizedLineage}(F_a,F_b)
\]

and

\[
d_R(F_a,F_b)\le\epsilon_R
\]

and the new required evidence is available/current:

\[
D_b\subseteq D(B_a)\cup D(B_{delta}).
\]

The new proof is then **re-derived**:

\[
B_b=F_b(B_a\cup B_{delta}),
\]

with provenance linking the new result to the old evidence / prior proof where appropriate.

This is not a proposed Web4 normative formula. It is the structural decomposition that seems necessary.

---

# 9. The evaluator also changes the temporal MRH

Phase 14 tied semantic freshness to dependency sets.

Therefore changing

\[
D_a\rightarrow D_b
\]

changes the invalidation horizon even when the current output does not change.

A dependency-expanding evaluator:

- needs more evidence;
- can be invalidated by more event classes;
- produces shorter-lived proof bundles, all else equal.

A dependency-contracting evaluator can do the opposite.

So evaluator selection affects not only computational cost and answers, but the **temporal durability of evidence**.

That is a direct connection between the parallel Hestia trust-model work and the Synchronism MRH arc.

---

# 10. Governance implication

A governed system should probably not represent an evaluator upgrade merely as:

```text
version 7 -> version 8
```

A relying party may need evidence about:

- old and new evaluator/plan commitments;
- authorized upgrade lineage;
- declared/enforced evidence contracts;
- semantic compatibility basis and its MRH/domain;
- dependency expansion or contraction;
- whether existing proof bundles were reusable as evidence;
- whether a new proof was actually re-derived;
- what state basis the new result commits to.

This makes evaluator evolution itself a first-class witnessed governance act.

---

## Phase-16 verdict

**PASS — evaluator replacement closes the component-replacement recursion at a higher fractal level.**

What survived:

1. Behavioral equivalence and evidence-horizon equivalence are independent.
2. An evaluator can be equivalent only inside a particular MRH/domain.
3. A replacement can preserve behavior while expanding the proof horizon; old proof bundles then become incomplete even before any answer changes.
4. A replacement can preserve behavior while contracting the proof horizon; this is a legitimate higher-level compression gain.
5. Proofs should be re-derived under a new evaluator rather than inherited as if the plan were unchanged.
6. Evaluator lineage/authorization is orthogonal to both semantics and evidence sufficiency.
7. Changing evaluator dependencies changes semantic-freshness / invalidation horizons.

What did not happen:

- no claim that finite-domain equality proves general evaluator equivalence;
- no universal semantic-distance metric;
- no claim that fewer dependencies are always safer or better;
- no protocol standard or physics result.

---

# Where the arc points next

The next question follows immediately:

> **What is the continuity object when evaluators themselves form a governed lineage and old proofs must be renewed across upgrades?**

The likely answer is not "trust inheritance."

It looks more like:

\[
\boxed{
\text{evidence persists; plans evolve; claims are re-derived.}
}
\]

That would make a proof/claim lineage distinct from both entity lineage and evaluator lineage.

The next phase should model those three DAGs together and ask when a claim can be renewed without returning to the full authoritative history.
