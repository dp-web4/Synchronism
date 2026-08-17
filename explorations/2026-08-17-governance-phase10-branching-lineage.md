# Governance Arc — Phase 10: Branching Lineage, Token Identity, and Why Trust Must Carry Provenance

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — lineage/provenance continuation of the Markov coherence arc  
**Code:** [`simulations/governance_phase10_lineage_dag.py`](../simulations/governance_phase10_lineage_dag.py)  
**Scope:** formalism / governance bridge; no physics claim

## Why the mathematical object changed

Phases 1–9 progressively moved the persistent object upward:

- state;
- metastable entity;
- relation;
- relation among relations;
- relational attractor;
- constitutional invariant shared across multiple governance mechanisms.

Phase 9 still treated historical continuity as a chain:

\[
E_0\rightarrow E_1\rightarrow E_2\rightarrow\cdots
\]

Real governed entities do not remain chains.

They fork, restore, federate, migrate, split, and merge.

So the next question is:

> **What becomes of token identity when witnessed lineage branches?**

At this point a Markov chain is no longer the natural primary representation. The relevant object is a **directed provenance graph**.

---

# The first correction: lineage is not an equivalence relation

Suppose

\[
A\rightarrow B,
\qquad
A\rightarrow C.
\]

Both \(B\) and \(C\) can possess perfectly valid witnessed causal lineage from \(A\).

But that does not imply

\[
B=C.
\]

So "shares valid lineage with" cannot mean "is numerically the same token as."

The natural structures are different:

### Constitutional/type equivalence

\[
B\sim_{constitution}C
\]

can be an equivalence relation: two current systems instantiate the same relevant governance invariants.

### Ancestry

\[
A\preceq B
\]

means there is a valid witnessed path from \(A\) to \(B\).

This is naturally a **partial order / reachability relation** on a provenance DAG, not an equivalence relation.

### Token identity

A current token is one historically situated branch of that DAG.

This means the word "same" must be split at least three ways:

1. same **type / constitutional class**;
2. same **lineage family / common ancestry**;
3. same **current historical token**.

Those are not interchangeable.

---

# Forks create descendants, not duplicate numerical identity

At a witnessed fork

\[
A\rightarrow\{B,C\},
\]

both descendants inherit legitimate ancestry from \(A\).

A useful rule is:

> **A fork preserves ancestry but ends exclusive token continuation.**

After the branch event, \(B\) and \(C\) should be treated as distinct current tokens even if:

- they begin with byte-identical state;
- they implement the same constitution;
- they inherit the same keys through some migration mechanism;
- they initially contain the same trust history.

Otherwise a copy operation would create two simultaneous instances each entitled to claim singular numerical identity.

This is the governance analogue of distributed-system split brain.

---

# Clone versus continuation

Phase 8 already separated structural similarity from historical continuity.

The DAG formulation makes that exact.

### Continuation

There is a witnessed authorized edge

\[
A\rightarrow B.
\]

### Clone

A system \(C\) independently reproduces the same state / constitution but has no valid incoming lineage edge from \(A\).

It may satisfy

\[
C\sim_{constitution}A,
\]

but not

\[
A\preceq C.
\]

Thus:

\[
\boxed{\text{state equality does not imply token continuity}}
\]

and neither does constitutional equivalence.

---

# Restore-from-backup is a fork unless exclusivity is established

This becomes particularly important for recovery.

Suppose a backup of \(A\) is restored as \(B\).

If the original \(A\) is provably gone and the recovery protocol witnesses the transition, it may be reasonable for relying parties to treat \(B\) as the exclusive continuation.

But if \(A\) might still be active, then restoration creates

\[
A\rightarrow B
\]

while \(A\) itself also continues.

That is a fork.

So safe numerical-identity transfer requires something like **fencing / exclusivity evidence**, not merely possession of old state or keys.

This suggests an additional continuity ingredient:

\[
\boxed{
\text{exclusive continuation}
=
\text{valid lineage}
+
\text{no concurrent sibling continuation}
}
\]

at whatever assurance level the relying party requires.

A persistent identifier alone cannot establish that.

---

# Merge creates a new descendant with multiple parents

Now suppose two branches later merge:

\[
B,C\rightarrow D.
\]

Then \(D\) has valid ancestry from both branches.

But saying

\[
D=B=C
\]

would erase the fact that \(B\) and \(C\) had distinct histories after the fork.

A cleaner interpretation is:

> **A merge creates a new token whose provenance has multiple parents.**

It may inherit constitutional continuity from both and evidence from both, but the merge event itself is a first-class historical transformation.

The resulting identity graph is therefore a DAG, not a chain and not a tree.

---

# The trust problem: naive inheritance manufactures evidence

Branching creates a concrete failure mode for trust systems.

Assume a root society has 100 unique witnessed evidence receipts.

At every generation:

1. the society forks in two;
2. both children inherit the full ancestral history;
3. each branch adds 10 genuinely new receipts;
4. the branches merge.

Compare two accounting methods.

### Wrong: inherit scalar score, then add parent scores at merge

### Safer: carry evidence provenance and derive trust from the union of unique receipts

The toy produces:

| fork/merge generation | naive inherited score | unique evidence receipts | artificial inflation |
|---:|---:|---:|---:|
| 0 | 100 | 100 | 1.00x |
| 1 | 220 | 120 | 1.83x |
| 2 | 460 | 140 | 3.29x |
| 3 | 940 | 160 | 5.88x |
| 4 | 1,900 | 180 | 10.56x |
| 5 | 3,820 | 200 | 19.10x |
| 6 | 7,660 | 220 | 34.82x |
| 7 | 15,340 | 240 | 63.92x |
| 8 | 30,700 | 260 | 118.08x |

Only 160 genuinely new receipts were created across eight generations.

Yet naive scalar inheritance turns the original evidence into an apparent score over **118 times** larger than the actual unique evidence base.

No fraud is required. The inflation arises purely from counting common ancestry multiple times.

This gives a strong governance rule:

\[
\boxed{
\text{trust evidence must merge by provenance-aware union, not scalar addition}
}
\]

A nonlinear trust function may still be applied afterward, but deduplication / dependence accounting must happen before treating inherited evidence as independent support.

---

# This strongly favors derived trust over stored trust

The result is particularly relevant to Web4/Hestia architecture.

If a trust score is stored as an inheritable scalar, fork/merge histories make its semantic meaning unstable.

If trust is instead **derived at read time from witnessed evidence with receipts/provenance**, then a relying party can:

- deduplicate common ancestry;
- distinguish inherited from branch-local evidence;
- weight evidence by age/context/assurance;
- detect conflicting branch histories;
- recompute under its own MRH and risk tolerance.

So the provenance DAG should be treated as primary data; trust scores are projections of that data for a particular evaluation context.

This matches the broader Web4 principle that trust is contextual rather than absolute.

---

# Trust does not conserve by scalar arithmetic

A merge should therefore not compute

\[
T_D=T_B+T_C.
\]

Even if \(T\) were a perfectly meaningful scalar within each branch, the two values are statistically dependent because they share ancestry.

A safer abstraction is

\[
T_D(\mathrm{MRH})
=
F_{\mathrm{MRH}}
\left(
\operatorname{UniqueEvidence}(B\cup C)
\right),
\]

where \(F\) is the relying party's contextual trust derivation.

The critical operation is the provenance-aware union before evaluation.

---

# Branching also makes MRH-relative identity explicit

After a fork, two descendants may be:

- the same **constitutional type** at one MRH;
- members of the same **lineage family** at another;
- distinct **historical tokens** at a finer governance MRH.

Thus "are these the same entity?" remains incomplete without the witness scale and question.

A hub asking whether two nodes implement compatible governance may legitimately quotient out the branch distinction.

A hub deciding whether a prior delegation, payment obligation, vote, sanction, or authorization transfers to both branches may absolutely not quotient it out.

So identity itself has a Markov Relevancy Horizon:

\[
\operatorname{MRH}_{identity}
(
\text{question},
\text{history},
\text{risk},
\text{tolerance}
).
\]

---

# A provisional event vocabulary for governed identity

The DAG suggests treating transformations as explicit typed events rather than trying to force all history into one "same/not-same" bit.

### Continuation

One parent, one child, compatible invariant, witnessed authorized transition.

### Migration

Continuation in which implementation/location changes but the governance class is preserved.

### Amendment

One lineage continues but the constitutional class itself changes under an authorized higher-order rule.

### Fork

One parent produces multiple simultaneously valid descendants.

### Clone

Structural/type equivalence without valid causal lineage.

### Merge

A new descendant has multiple parent lineages.

### Recovery

A descendant restores earlier state; whether it is continuation or fork depends on exclusivity/fencing evidence.

### Extinction

A token has no valid continuing descendant.

This vocabulary is likely more useful to Web4 than trying to assign one universal metaphysical identity answer.

---

# LCT implication

A Linked Context Token / identity object should therefore be thought of less like a permanent name attached to one immutable object and more like a **provenance-addressable evolving token**.

At minimum, a governance-aware lineage representation needs to express:

- parent token(s);
- transformation/event type;
- governing rule authorizing the transition;
- invariant/constitution version before and after;
- evidence/witnesses for the transition;
- whether the event establishes exclusive continuation or a branch;
- descendant/fork relationships where visible;
- unique evidence references rather than only inherited aggregate scores.

This is a conceptual implication, not a claim about current LCT implementation completeness.

---

## Phase-10 verdict

**PASS — branching forces identity and trust to become provenance-DAG concepts.**

What changed:

1. Historical ancestry is a directed reachability relation, not an equivalence relation.
2. Constitutional/type equivalence, common lineage, and current token identity must be separated.
3. A fork preserves ancestry but creates distinct current tokens.
4. A merge creates a descendant with multiple parents rather than making the parents retroactively identical.
5. Recovery from backup needs exclusivity/fencing evidence to avoid accidental split-brain identity.
6. Scalar trust inheritance is unsafe under fork/merge because shared ancestral evidence can be multiplied without creating new evidence.
7. Trust should be derived from provenance-aware unique evidence in the relying party's MRH.

What did not happen:

- no universal trust formula is proposed;
- no claim that all historical identity questions have one objective answer;
- no claim that a DAG alone proves exclusivity or prevents Sybil attacks;
- no physics result.

---

# The next problem: constitutional drift and Ship of Theseus across amendments

Fork/merge history exposes one more unresolved question.

A society can maintain continuous lineage while slowly changing its constitution:

\[
C_0\rightarrow C_1\rightarrow C_2\rightarrow\cdots
\]

Each amendment may be locally acceptable, yet after enough small changes the endpoint may no longer resemble the constitution that an old relying party originally trusted.

This is a governance Ship-of-Theseus problem.

So the next question is:

> **How much accumulated constitutional drift can occur before trust should stop transferring automatically, even though lineage remains perfectly continuous?**

That appears to require a path-length / semantic-drift notion on the provenance DAG, indexed again by the relying party's MRH and stakes.
