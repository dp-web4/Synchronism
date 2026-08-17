# Markov Coherence Arc — Phase 5: Relations Among Relations as Slow Invariants

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 5 toy-model result  
**Code:** [`simulations/markov_phase5_relational_hierarchy.py`](../simulations/markov_phase5_relational_hierarchy.py)  
**Scope:** formalism probe only; no physics claim

## Why this phase changed

Phase 4 showed that relations themselves can be metastable objects. A stable anti-aligned relation can be just as coherent as an aligned one, so the next step should not be ordinary object aggregation.

The new question is:

> **Can a dynamical search discover a hierarchy of relations-among-relations without being told what those relational variables are?**

If yes, the fractal hierarchy can be expressed more generally as a hierarchy of slow predictive invariants, rather than only nested spatial aggregates.

---

## Construction

Use eight binary microvariables

\[
s_1,\ldots,s_8\in\{-1,+1\}.
\]

For construction only, group them into four base pairs:

\[
(1,2),(3,4),(5,6),(7,8).
\]

Define the following relational variables conceptually:

\[
r_1=s_1s_2,\quad
r_2=s_3s_4,\quad
r_3=s_5s_6,\quad
r_4=s_7s_8.
\]

Then define relations among those relations:

\[
q_1=r_1r_2=s_1s_2s_3s_4,
\]

\[
q_2=r_3r_4=s_5s_6s_7s_8,
\]

and finally

\[
Q=q_1q_2=\prod_{i=1}^{8}s_i.
\]

The stochastic dynamics contains several classes of state transformations:

1. **fast pair-gauge moves:** flip both members of one base pair; these erase pair orientation while preserving its relation \(r_i\);
2. **intermediate within-quartet moves:** alter the constituent pair relations while preserving \(q_1\) or \(q_2\);
3. **slower across-quartet moves:** alter the quartet relations while preserving \(Q\);
4. **rare single flips:** break all relational levels containing that variable.

The hierarchy is therefore deliberately present in the transformation structure.

But the analysis is **not supplied the hierarchy**.

Instead it searches every nonconstant parity observable of the eight variables:

\[
\prod_{i\in S}s_i,
\qquad
\emptyset\neq S\subseteq\{1,\ldots,8\}.
\]

There are

\[
2^8-1=255
\]

such observables.

For this translation-invariant random walk on the binary state space, every parity observable is an exact eigenfunction of the transition operator. Its eigenvalue directly measures one-step persistence; larger positive eigenvalues correspond to slower-decaying relational modes.

Thus the search is simply:

> among all 255 candidate microstate and relational observables, which ones are dynamically slowest?

---

# Result — the blind spectral ranking reconstructs the relational hierarchy

The slowest modes were:

| rank / class | observable | eigenvalue \(\lambda\) | correlation time \(\tau=-1/\ln|\lambda|\) |
|---|---|---:|---:|
| 1 | \(Q=s_1s_2s_3s_4s_5s_6s_7s_8\) | 0.9000 | 9.491 |
| 2–3 | \(q_1=s_1s_2s_3s_4\), \(q_2=s_5s_6s_7s_8\) | 0.6300 | 2.164 |
| 4–7 | \(r_1,r_2,r_3,r_4\) | 0.5750 | 1.807 |
| next | individual \(s_i\) | 0.5475 | 1.660 |

So the dynamical search recovered, in order:

\[
\boxed{
\text{relation among quartet-relations}
\;>\;
\text{quartet-relations}
\;>\;
\text{pair-relations}
\;>\;
\text{individual states}
}
\]

without being told those relational variables beforehand.

This is the strongest formalism result in the arc so far — not because the mathematics is new, but because it forces a refinement of the conceptual picture.

---

# The key correction: higher-scale identity can survive loss of lower-scale identity

Consider a fast pair-gauge move that flips both members of a pair:

\[
(s_1,s_2)\rightarrow(-s_1,-s_2).
\]

Both constituent microstates change completely, but

\[
r_1=s_1s_2
\]

is unchanged.

At the next level, transformations can change \(r_1\) and \(r_2\) while preserving

\[
q_1=r_1r_2.
\]

And still larger transformations can alter the quartet relations while preserving

\[
Q=q_1q_2.
\]

So persistence does not require persistence of constituents.

This gives a cleaner definition of an emergent variable:

> **An emergent identity is a slow invariant under a class of faster lower-scale transformations.**

The important object is therefore not necessarily a collection of things. It may be an equivalence class over many changing lower-scale realizations.

---

# A possible mathematical reading of "fractal emergence"

The previous picture was approximately:

\[
\text{small entities}
\rightarrow
\text{larger entities}
\rightarrow
\text{still larger entities}.
\]

Phase 5 suggests replacing that with:

\[
\boxed{
\text{fast degrees of freedom}
\rightarrow
\text{slow invariants of their transformations}
\rightarrow
\text{slow invariants among those invariants}
\rightarrow\cdots
}
\]

This is much more general.

A molecule need not be interpreted only as a bag of atoms. A biological identity need not be the persistence of particular molecules. A social structure need not preserve particular people. The higher-scale identity can reside in **relations that remain invariant while lower-scale constituents turn over**.

Those examples are analogies, not claims that this toy captures their actual dynamics.

---

# Entityhood becomes transformation-relative

This also refines the Phase-2 result that entityhood is timescale-relative.

A candidate entity is now relative to at least three things:

\[
E=(\text{observable},\text{transformation class},\text{timescale}).
\]

An observable is entity-like when:

1. many faster transformations change its micro-realization;
2. those transformations leave the observable approximately invariant;
3. the observable remains predictively useful over the chosen horizon;
4. slower transformations eventually change or destroy it.

In spectral language, the hierarchy is encoded by separated slow modes.

In Synchronism language, the hierarchy may be encoded by **patterns of Intent relations that recur despite continuous lower-scale redistribution**.

The second phrasing is interpretation; the first is the tested mathematical structure.

---

# MRH consequence

MRH now needs to answer not just "how far away does relevant information extend?" but also:

> **Which lower-scale distinctions can be discarded because they lie inside the equivalence class defining the current entity?**

That suggests two complementary boundaries:

### Exterior MRH

How much environmental information must be retained for prediction?

### Interior resolution horizon

How much microstate information must be retained before further detail becomes predictively irrelevant to the macro-identity?

For a macro-observable \(Y=f(X)\), one can ask for a retained internal representation \(Z(X)\) such that

\[
I(Y_{t+h};X_t\mid Z_t)\le\epsilon_{int}.
\]

This is formally parallel to the exterior MRH criterion:

\[
I(Y_{t+h};E_t\mid B_t)\le\epsilon_{ext}.
\]

Thus an emergent entity may occupy an **informational shell**:

- excessive internal detail is irrelevant below it;
- sufficiently distant external detail is irrelevant beyond it;
- the predictive identity lives between those two discarded domains.

This is provisional, but it gives MRH a natural partner on the inside.

---

# Relations among relations are not merely higher-order correlations

A warning is important here.

Any product such as

\[
Q=q_1q_2
\]

can be written algebraically whether or not it matters dynamically.

What makes a higher-order relation entity-like is not that it can be defined, but that the dynamics selects it as a **slow, persistent, predictively useful mode**.

So the criterion is not:

> "Can I write a relation among relations?"

It is:

> "Does that relation survive transformations that erase its lower-level implementation?"

That keeps the hierarchy empirical/model-dependent instead of allowing arbitrary symbolic nesting.

---

# What is known versus what is Synchronism-specific

The mathematics here is standard spectral analysis of Markov operators on a deliberately simple binary group. Slow eigenmodes, invariants, conserved quantities, and multiscale coarse-graining are established ideas.

Therefore this phase makes **no novelty claim**.

What may be useful in Synchronism is the synthesis:

1. entityhood as time-scale-separated recurrence;
2. emergence as predictive closure;
3. MRH as tolerated exterior predictive loss;
4. interaction as MRH overlap plus relational mode structure;
5. recursive emergence as slow relations among changing lower-scale relations.

Whether that synthesis is genuinely useful outside toy models remains open.

---

## Phase-5 verdict

**PASS — relational recursion survived a blind dynamical search.**

What changed:

1. A higher-scale entity need not preserve lower-scale constituent identities.
2. Persistent relations can outrank constituent states in dynamical lifetime.
3. Relations among relations can themselves be the slowest modes.
4. A promising operational definition of emergence is **slow invariance under faster transformations**.
5. The fractal hierarchy may be better represented by nested invariants / equivalence classes than by nested objects.
6. MRH may have a complementary **interior resolution** notion: irrelevant micro-detail below, irrelevant environmental detail beyond.

What did not happen:

- no new mathematics;
- no physics prediction;
- no evidence that parity observables are universal;
- no demonstration yet that arbitrary nonlinear dynamics yield similarly clean recursive relational levels.

---

# Next question

The current toy makes the hierarchy especially easy because parity observables are exact eigenfunctions. The next test should remove that privilege.

Construct a nonlinear stochastic system where the useful relational variables are **not known analytically**, then ask a data-driven method to discover slow observables from trajectories alone.

The sharper question is:

> **If the analyst is given only raw microstate time series, can a relational hierarchy emerge from prediction itself — without supplying candidate relations?**

If that works, the result would move from "the correct relational observable is spectrally slow" toward the stronger claim that **predictive dynamics can discover the relational ontology it needs.**
