# Markov Coherence Arc — Phase 6: Predictive Dynamics Can Recover a Hidden Relational Ontology

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 6 toy-model result  
**Code:** [`simulations/markov_phase6_data_driven_relational_discovery.py`](../simulations/markov_phase6_data_driven_relational_discovery.py)  
**Scope:** formalism probe only; no physics claim

## Question

Phase 5 recovered a recursive hierarchy of relations among relations, but the analysis still benefited from an unusually convenient mathematical fact: parity observables were exact eigenfunctions of the chosen process.

The stronger question is:

> **If the learner is given only raw microstate transitions, with no candidate relational observables supplied, can the slow predictive structure itself recover the hidden hierarchy?**

This matters for the component-replacement problem. If larger identity can survive replacement of lower-level components, then the useful ontology should ideally be discoverable from persistence/prediction rather than inserted by an analyst.

---

## Construction

Use the same eight binary microvariables

\[
s_1,\ldots,s_8\in\{-1,+1\},
\]

with a deliberately hidden hierarchy of transition channels resembling Phase 5:

- frequent pair-orientation replacements;
- slower changes of relations within quartets;
- still slower changes across quartets;
- rare single-variable changes.

But two changes remove the exact analytic privilege.

### 1. Add state-dependent nonlinear transition rates

A weak random microstate energy is drawn once from random local fields and pair couplings. Proposed moves are accepted with a Metropolis probability

\[
\min\left(1,e^{-\beta\Delta E}\right).
\]

Transition probabilities therefore depend on the current microstate. The clean translation symmetry of Phase 5 is broken, and the intended relational variables are no longer exact eigenfunctions.

### 2. Give the learner only trajectory counts

A long trajectory is generated from the 256-state Markov process. The learner receives only observed microstate-to-microstate transition counts.

It then:

1. constructs a reversible empirical transition operator from the counts;
2. diagonalizes that operator;
3. ranks the learned nonconstant eigenmodes by persistence.

No pair definitions, quartet definitions, parity functions, or relational features are passed to the eigensolver.

Only **after** learning are the discovered modes compared against the hidden construction variables as diagnostics:

\[
Q=\prod_{i=1}^{8}s_i,
\]

\[
q_1=s_1s_2s_3s_4,\qquad q_2=s_5s_6s_7s_8,
\]

and

\[
r_1=s_1s_2,\ldots,r_4=s_7s_8.
\]

---

# Result — the leading empirical modes recover the hidden relations

Independent validation of the committed construction with 300,000 sampled transitions produced approximately:

| learned mode | eigenvalue | correlation time | strongest hidden diagnostic |
|---:|---:|---:|---|
| 1 | 0.9641 | 27.38 | \(Q\), \(|\rho|\approx0.9998\) |
| 2 | 0.9028 | 9.78 | \(q_1\), \(|\rho|\approx0.9949\) |
| 3 | 0.8890 | 8.50 | \(q_2\), \(|\rho|\approx0.9946\) |
| 4–5 | ~0.775 | ~3.9 | strong mixtures dominated by base relations \(r_3,r_4\) |
| later slow modes | ~0.71–0.77 | ~2.9–3.8 | remaining base-relation subspace / nuisance mixtures |

The committed script uses 700,000 transitions for a more stable estimate, but the qualitative ordering is already clear in the smaller independent run.

The striking part is not that the hidden variables are correlated with learned modes; the process was built to contain that hierarchy. The important methodological point is:

> **The learner was never told what a relation was.**

It was given only raw state-transition statistics, and the slow predictive eigenspace independently reconstructed the high-level relational variables.

---

# The hierarchy survives loss of exact algebraic form

In Phase 5, one could reasonably object that parity variables were privileged by the model's exact symmetry.

Here the random state-dependent energy mixes and distorts the slow modes. Yet the top modes remain extremely close to the hidden relational hierarchy:

\[
Q\rightarrow q_1,q_2\rightarrow r_i\text{-dominated subspace}\rightarrow\text{faster microstructure}.
\]

So the hierarchy is not merely an artifact of supplying a clever coordinate system to the analysis.

At least in this finite toy, the predictive dynamics itself contains enough information to infer a coordinate system close to the one in which the hierarchy is simple.

That is a stronger result than Phase 5.

---

# Component replacement becomes a quotient-space idea

The original motivating intuition was:

> a larger entity can remain "the same" while its components are replaced.

The Markov arc now gives this a more precise interpretation.

Suppose many microstates \(x\) correspond to the same slow macro-observable \(Y(x)\). Then replacement, turnover, or internal rearrangement can move the system through the equivalence class

\[
[x]_Y=\{x':Y(x')=Y(x)\}
\]

without changing the macro-identity.

The entity is therefore not identified with one microstate or one inventory of components. It is associated with a **quotient of microstate space** under distinctions that the slow predictive dynamics treats as irrelevant.

In plain language:

> **Identity is what remains invariant across the transformations that count as internal replacement.**

This seems substantially closer to the intended meaning of recurrence in Synchronism than literal constituent persistence.

---

# But invariance alone is still insufficient

A caution from earlier phases remains essential.

One can construct infinitely many invariants or coarse labels. What selects an entity-like identity is not merely that some function remains unchanged.

The stronger candidate conditions are now:

1. **dynamical persistence** — the macro-observable changes slowly;
2. **internal replacement tolerance** — many faster micro-transformations leave it approximately unchanged;
3. **predictive closure** — discarded micro-detail contributes little additional predictive power;
4. **external closure / MRH** — sufficiently distant exterior detail contributes little additional predictive power;
5. **discoverability** — the slow variable can in principle be inferred from raw dynamics rather than requiring an arbitrary semantic label.

The fifth condition is new in this phase and may be especially useful as a guard against analyst-imposed ontology.

---

# Relations among relations become learned coordinates

Phase 5 phrased the hierarchy as

\[
\text{fast degrees of freedom}
\rightarrow
\text{slow invariants}
\rightarrow
\text{slow invariants among those invariants}.
\]

Phase 6 adds another layer:

\[
\boxed{
\text{raw dynamics}
\rightarrow
\text{learned slow coordinates}
\rightarrow
\text{relational hierarchy in those coordinates}
}
\]

This suggests that an ontology need not be fundamental or externally imposed. A useful ontology can be an **emergent coordinate system selected by predictive persistence**.

That is perhaps the most interesting conceptual result of this phase.

The microvariables remain perfectly real as degrees of freedom. But they are not necessarily the most useful variables for describing identity at longer horizons.

The learned slow coordinates answer a different question:

> Which distinctions continue to matter after the fast transformations have done their work?

---

# A possible refinement of MRH

MRH has so far evolved from a spatial intuition into a target-, horizon-, tolerance-, and task-relative predictive boundary.

The component-replacement result suggests an even more general formulation:

> **An MRH is a relevance boundary in information space, not necessarily a boundary in physical space.**

There are at least two directions of irrelevance:

- outward: external variables whose additional predictive value is negligible;
- inward: microstate distinctions whose additional predictive value is negligible.

A macro-identity may therefore be specified by the minimal predictive representation that sits between these discarded domains.

If \(Z_t\) is a learned representation of the current state sufficient to predict target \(Y_{t+h}\), then one wants approximately

\[
I(Y_{t+h};X_t\mid Z_t)\le\epsilon_{int}
\]

and

\[
I(Y_{t+h};E_t\mid Z_t)\le\epsilon_{ext}.
\]

This is still a research scaffold, not a finalized Synchronism definition.

---

# What this does and does not establish

The mathematics remains standard spectral / transfer-operator reasoning on a synthetic Markov process. The process was deliberately constructed with a hidden relational hierarchy.

Therefore this result does **not** show that arbitrary natural systems must possess such hierarchies, and it does not provide a novel physics prediction.

What it does establish in the toy is narrower and useful:

> A recursive relational identity can remain discoverable from raw dynamics even after the exact analytic variables defining that identity are hidden from the learner and distorted by state-dependent nonlinear transition rates.

That clears one methodological obstacle for the broader idea.

---

## Phase-6 verdict

**PASS — the relational hierarchy survived data-driven discovery from raw microstate transitions.**

What changed:

1. The learner no longer needs candidate relational variables supplied in advance.
2. The dominant learned slow mode recovers the highest-order relation almost perfectly.
3. The next modes recover the intermediate relational level.
4. Lower relational modes appear as a slower subspace even when exact symmetry is broken.
5. Component replacement can be formalized as motion within an equivalence class that preserves a slow predictive identity.
6. A useful emergent ontology may be a coordinate system selected by predictive persistence rather than an inventory of persistent components.
7. MRH increasingly looks like a boundary in **relevance/information space**, with both internal and external directions.

Still absent:

- a natural-system demonstration;
- a nonlinear continuous-state example;
- a proof that the discovered representation is unique;
- a principled rule for selecting how many slow coordinates constitute one entity;
- new physics.

---

# What now looks like the sharper question

The next problem is no longer simply "can we discover the slow variables?"

It is:

> **When does a collection of learned slow variables deserve to be treated as one entity rather than merely several correlated observables?**

That points directly at identity through component turnover.

A strong next test would create a system in which:

- individual components enter and leave;
- no fixed subset of components persists;
- a relational organization survives;
- the learned macrostate remains predictive across the turnover events.

If a data-driven representation tracks the organization continuously while its constituent identities are explicitly replaced, that would be a much closer operational analogue of the biological/social/organizational identity question that motivated this branch.
