# Markov Coherence Arc — Phase 8: Identity as a Relational Attractor Under Membership and Topology Turnover

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 8 toy-model result  
**Code:** [`simulations/markov_phase8_topology_turnover.py`](../simulations/markov_phase8_topology_turnover.py)  
**Scope:** formalism / governance bridge; no physics claim

## Why this phase changed

Phase 7 allowed every component to be replaced while a larger relational identity persisted, but it preserved a fixed role topology. That left an obvious loophole:

> perhaps the larger identity was really stored in the slots.

So Phase 8 removes another layer of persistence.

The stronger question is:

> **Can a recognizable macro-organization persist after both its original members and its original relations are gone?**

This matters beyond the physics-language arc. It is directly relevant to governance systems, institutions, software societies, hubs, fleets, and other entities whose participants and interaction graph can both change while the larger organization is expected to remain continuous.

---

## Construction

Use 12 current members connected by a cycle graph.

Each member carries a binary state. The organization is not defined by a particular member, slot label, or edge. Instead it is defined relationally:

> neighboring members tend to occupy opposite states.

Equivalently, the macro-observable is the fraction of current edges that satisfy that anti-aligned relationship:

\[
S(G,x)
=\frac{1}{|E|}
\sum_{(u,v)\in E}
\mathbf 1[x_u\neq x_v].
\]

The process contains three event classes.

### 1. Relational repair

Most events are local updates. A member looks only at its current neighbors and probabilistically adjusts toward the current relational organization.

### 2. Membership turnover

A current member is removed and replaced by a brand-new identity with a random initial state.

The new member does not inherit its predecessor's state.

### 3. Topology turnover

A 2-opt graph move changes current adjacency while keeping the population connected as one even cycle.

Thus not only the members but also the particular pairwise relations are continually replaced.

For every episode we explicitly track:

- all original member identities;
- all original edges.

The episode is evaluated only once **both sets have vanished**.

At that point there is no original component and no original relation left to carry the initial organization directly.

---

# Result

Using 1,000 episodes per condition:

| replacement/event | rewiring/event | mean time until all original members + edges are gone | mean current edge satisfaction | fraction with at least 10/12 relations satisfied |
|---:|---:|---:|---:|---:|
| 0.02 | 0.001 | ~1867 steps | 0.850 | 0.829 |
| 0.03 | 0.003 | ~1260 steps | 0.825 | 0.749 |
| 0.05 | 0.003 | ~749 steps | 0.800 | 0.685 |

For comparison, a uniformly random assignment of binary states on the same 12-cycle has

\[
\mathbb E[S]=0.5,
\]

and the exact probability of satisfying at least 10 of the 12 relations by chance is only

\[
P(S\ge10/12)\approx0.0327.
\]

So in the moderate-turnover regimes, the relational organization remains strongly present after complete turnover of both the original members and the original edges.

This is not an immortal identity. Raising turnover eventually erodes it. The point is that continuity is not reducible to material or relational inventory continuity.

---

# The conceptual update

Phase 7 suggested:

> identity is continuously maintained by a relation that recruits changing components into an equivalence class of implementations.

Phase 8 strengthens that:

> **Identity can persist as a relational attractor even when both the components and the specific relations implementing it are replaced.**

The persistent object is therefore not:

- the member set;
- the edge set;
- the state of any fixed component;
- or any particular pairwise relation.

It is a **constraint on the evolving relation field**.

This is closer to a process ontology than an object ontology.

---

# A hierarchy of what can turn over

The arc now distinguishes at least four increasingly strong continuity tests:

1. **state turnover** — microstate changes while components remain;
2. **component turnover** — members change while roles/relations remain;
3. **relation turnover** — particular connections change while an organizational relation remains;
4. **organizational turnover** — if the governing relational attractor itself changes, the macro-identity is no longer the same one at that MRH.

So the phrase "same entity" is incomplete unless the witness specifies which transformation classes are considered identity-preserving.

A more explicit object is

\[
E=(Y,\mathcal T, h,\epsilon),
\]

where

- \(Y\) is the macro-observable / relational organization;
- \(\mathcal T\) is the class of transformations treated as internal replacement;
- \(h\) is the temporal horizon;
- \(\epsilon\) is tolerated predictive/identity loss.

This dovetails with the evolving MRH notation:

\[
\operatorname{MRH}(Y,h,\epsilon,\mathcal T).
\]

MRH is therefore not merely spatial, and not merely temporal. It is indexed by the **transformations the witness agrees to quotient out**.

---

# Fractal MRH as quotienting

This gives a potentially clean mathematical interpretation of the old intuition that MRH is fractal.

At one scale, a particular component replacement is a major event.

At a higher scale, the same event is quotiented out as implementation detail:

\[
X\sim_{\mathcal T}X'
\]

if the transformation from microstate \(X\) to \(X'\) preserves the macro-variable relevant to the current witness and horizon.

A higher MRH therefore sees an equivalence class

\[
[X]_{\mathcal T}
\]

where a lower MRH sees many distinct states.

The fractal hierarchy is not just a hierarchy of sizes. It is a hierarchy of **which distinctions still matter**.

That may be the cleanest formal expression of MRH produced by this arc so far.

---

# A necessary refinement: invariance is not enough

There is a governance-relevant problem with defining identity only as a relational invariant.

Suppose an entirely separate system is built tomorrow with the same macro-organization.

It may realize the same invariant, but it is not automatically the **same historical entity**.

So two concepts must be separated:

### Type identity

The systems instantiate the same relational organization / equivalence class.

### Token continuity

One current entity is the historical continuation of a particular earlier entity.

Token continuity requires more than structural similarity. It needs a **causal/provenance path** through the sequence of transformations.

This suggests:

\[
\boxed{
\text{persistent token identity}
=
\text{relational invariance}
+
\text{witnessed lineage}
}
\]

That second term is crucial. Otherwise a perfect clone can silently claim continuity merely by matching the invariant.

---

# Why this matters for Web4 / Hestia / hub governance

Hestia currently frames governance as the basis for hub admissibility: identity, composed law, witnessed action history, and trust derived from that record give the hub something checkable rather than merely asserted.

The Markov arc suggests a deeper way to interpret that architecture.

A governed society should probably **not** be identified with:

- its current agent roster;
- one model vendor;
- one machine;
- a fixed delegation graph;
- or one static policy implementation.

Those are all candidate lower-scale components or relations that may legitimately turn over.

The more durable governance identity is closer to:

1. a set of slow constitutional / policy invariants;
2. a witnessed chain showing how each allowed transformation followed from the previous state;
3. current evidence that the transformed society still satisfies those invariants.

That maps naturally onto the distinction above:

- **governance invariant** supplies type continuity;
- **witness/provenance chain** supplies token continuity.

This may matter especially for future federation, fleet evolution, agent replacement, delegation changes, policy-entity replacement, and hub migration.

The question for a relying hub then becomes more precise than "is this the same society?"

It becomes:

> **Which transformations occurred, which governance invariants survived them, and is there a witnessed causal chain connecting the present state to the state I previously trusted?**

That is a much stronger basis for continuity than a persistent identifier alone.

---

# Governance MRH

This also suggests a governance-specific MRH.

For a relying party \(R\), action class \(A\), horizon \(h\), and risk tolerance \(\epsilon\):

\[
\operatorname{MRH}_{gov}(R,A,h,\epsilon)
\]

could denote the minimum governance context required before additional interior/exterior history changes the relying party's decision by less than \(\epsilon\).

At low stakes, the relevant MRH might contain only:

- current identity key;
- current role/delegation;
- recent witnessed conduct.

At higher stakes, the MRH might expand to include:

- longer provenance history;
- policy lineage;
- peer witness consistency;
- device/constellation evidence;
- prior governance transitions;
- federation ancestry.

Thus "how much history must be carried forward?" becomes a relevance question, not a universal retention rule.

This is directly consonant with the broader Web4 principle that trust should scale to context and stakes rather than collapse to one global score.

---

# Important caveat

This toy still preserves one architectural class: the graph remains a connected even cycle.

So Phase 8 has not shown that identity survives arbitrary topology changes. The current macro-identity is partly the attractor induced by that architecture class.

That is not a defect in the result; it identifies the next boundary.

The next test should allow the **architecture class itself** to mutate:

- cycle to modular graph;
- modular graph to another topology;
- changing population size;
- changing local update rules;
- perhaps changing governance mechanism while preserving some higher invariant.

Then the hard question becomes:

> **How much mechanism can change before a higher-level relation no longer has legitimate continuity?**

That is much closer to constitutional evolution than component replacement.

---

## Phase-8 verdict

**PASS — relational organization survived complete component and edge turnover in a nontrivial parameter regime.**

What changed:

1. Fixed component identity is not required for macro-identity.
2. Fixed pairwise relations are not required either.
3. The persistent object can be a constraint / attractor over an evolving relation field.
4. MRH can be interpreted as a hierarchy of distinctions quotiented out by the witness.
5. Structural invariance alone is insufficient for historical identity.
6. **Relational invariance + witnessed lineage** is a stronger candidate for token continuity.
7. This maps naturally onto governance systems where roster, delegation graph, agents, and implementation can change while constitutional identity persists.

Still absent:

- arbitrary topology change;
- changing organization size;
- changing update/governance mechanism;
- a theorem specifying when lineage + invariance is sufficient;
- any physics claim.

# Next question

Phase 9 should attack the architecture itself.

Instead of asking whether identity survives changing components inside one architecture class, allow several governance/dynamical mechanisms to implement the same higher relation.

Then ask:

> **Can a macro-identity migrate between implementations without a discontinuity in its predictive / governance invariants?**

That would be the analogue of constitutional continuity across a change of machinery — and likely the most direct bridge from this arc into Web4/Hestia/hub governance.
