# Governance Arc — Phase 18: A Generic Relevance Contract

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — architectural synthesis / Hestia-Web4 bridge  
**Code:** [`simulations/governance_phase18_relevance_contract.py`](../simulations/governance_phase18_relevance_contract.py)  
**Scope:** candidate common evidence-basis abstraction; not a protocol standard, security proof, or physics claim

## Why this phase exists

Phases 13–17 repeatedly produced the same metadata under different names:

- what proposition is being evaluated;
- which subject/token the proposition is about;
- which evaluator/plan produced the answer;
- which authoritative state the answer is based on;
- what evidence/dependencies the evaluator may consume;
- whether the supplied evidence is complete **for this claim**;
- which intervening changes invalidate the result;
- what assurance supports the evidence;
- which entity/evaluator/prior-claim lineage the result descends from;
- how to escalate/re-expand the evidence horizon when the compressed basis is insufficient.

If continuity, trust, policy, authority, occupancy, and assurance each invent a separate version of this structure, the system will accumulate subtly incompatible meanings of words like `complete`, `fresh`, `basis`, `current`, and `derived`.

So the question is:

> **Is there one small meta-contract that can describe the relevance boundary of many different Web4/Hestia claims without replacing their domain-specific proof formats?**

The answer appears to be yes, with an important restraint:

> **use a common envelope/interface, not one giant universal proof object.**

---

# 1. Existing architecture already contains fragments of this contract

This is not a greenfield abstraction.

## Hestia `ReadBasis`

Current Hestia dashboard code carries:

- `mode` — e.g. `windowed-projection` versus `full-traversal`;
- `window` — bounded rows considered;
- `complete` — whether the answer rests on the whole chain;
- `note` — what the relying party should know before accepting the compressed reading.

Its default is deliberately fail-closed: an unstated basis is treated as compressed/incomplete, never implicitly complete.

That directly corresponds to **coverage/completeness metadata**.

## Hestia trust derivation dependency inventory

The current trust path declares:

- `DERIVATION_EVENT_TYPES`;
- `DERIVATION_KEYS`.

The chain reader projects to those model-declared evidence needs.

That directly corresponds to an **evaluator dependency contract**.

Phase 15 already noted the remaining hardening opportunity: make the evaluator incapable of undeclared reads rather than merely keeping the declaration beside the code.

## Hestia evidence/occupancy typing

The governance work distinguishes evidence class and occupancy basis specifically to avoid presenting a weaker basis as though it were a stronger one.

That corresponds to **typed assurance/basis semantics** rather than one scalar confidence number.

## Hub AssuranceReceipt

The hub already has a portable signed AssuranceReceipt mechanism intended to let an external relying party verify A2 evidence without trusting Hestia itself.

That corresponds to **cryptographically portable basis evidence**.

## LCT lineage and the witness chain

LCT lineage already records token evolution, while the witnessed chain provides the authoritative historical record from which derived views can be projected.

That corresponds to **subject lineage + archival provenance**.

So Phase 18 is mostly a refactoring insight:

> several parallel pieces are implementing aspects of one higher-level relation between claim, evidence, plan, state, and relying party.

---

# 2. Candidate abstract shape

Call the meta-object a **Relevance Contract** for now.

For one derived proposition, write conceptually:

\[
\mathcal R
=
(Q,S,F,B,D,C,A,L,X),
\]

where:

- \(Q\): claim/query being evaluated;
- \(S\): subject/token/context;
- \(F\): evaluator/plan identity;
- \(B\): authoritative basis state/commitment;
- \(D\): dependency contract;
- \(C\): coverage/completeness statement;
- \(A\): assurance evidence/type;
- \(L\): relevant lineage references;
- \(X\): escalation/re-expansion rule.

This tuple is descriptive, not a normative wire schema.

A useful implementation may expose these through a trait/interface while preserving very different domain-specific receipts underneath.

---

# 3. Dependencies need temporal semantics

Phase 17 showed that a dependency name is insufficient.

At minimum, dependencies may be:

### Snapshot-sensitive

The claim needs a current authenticated value/commitment.

Example shape:

```text
current role occupancy
current law hash
current registry membership
```

### Path-sensitive

The claim needs relevant historical events since some anchor.

Example shape:

```text
no sanction since admission
valid lineage since anchor
no trust-breaking constitutional excursion
```

So the dependency contract needs something conceptually like:

\[
D=D_{snapshot}\cup D_{path}(anchor).
\]

This distinction determines both:

- what evidence must be fetched;
- what future events invalidate the claim.

---

# 4. Fetch horizon and invalidation horizon should be one contract

This is one of the strongest recurring results of the arc.

The same dependency declaration should drive both directions in time.

## Backward-looking question

What evidence from history/current state must be read?

\[
G_R=\operatorname{Project}(G,D_R).
\]

## Forward-looking question

Which future changes make the existing proof stale?

\[
\operatorname{Invalidate}(B,e)
\iff
\operatorname{Touches}(e,D_R).
\]

If fetch dependencies and invalidation dependencies are maintained as separate hand-written lists, they can drift.

So:

\[
\boxed{
\text{one dependency contract should determine both evidence acquisition and semantic freshness.}
}

That is the most concrete software consequence of the MRH formalization so far.

---

# 5. Completeness must always name its question

`complete` by itself is dangerous.

The arc has repeatedly needed two different statements:

### Complete for claim

\[
C_{claim}=true
\]

means the supplied basis contains everything the declared evaluator requires to evaluate this proposition under its stated contract.

### Complete history

\[
C_{history}=true
\]

means the entire authoritative historical record is present.

A small proof bundle may legitimately satisfy the first and not the second.

Therefore the generic abstraction should avoid a bare universal `complete=true` where possible.

Prefer something like:

```text
complete_for: exclusive-continuation@plan-hash
full_history_complete: false
```

or equivalent typed semantics.

---

# 6. Composition is where the contract earns its keep

A consequential action often requires multiple propositions:

\[
Q_{action}
=Q_{continuity}
\land
Q_{policy}
\land
Q_{trust}
\land
Q_{authority}.
\]

Let each proposition carry relevance contract \(\mathcal R_i\).

The action-level dependency horizon is then approximately the union:

\[
D_{action}=\bigcup_i D_i.
\]

But composition requires more than set union.

Each component must be:

- complete for its stated claim;
- semantically current or validly carry-forwardable to the action state;
- about compatible subjects/lineages;
- produced under a recognized plan;
- supported by assurance the relying party accepts for that sub-claim.

The toy deliberately keeps assurance as a vector:

```text
[A2 continuity, A2 policy, A1 trust, A2 authority]
```

rather than averaging it into a synthetic `A1.75`.

That preserves the exact weak surface.

---

# 7. Toy composition result

The toy gives four proof contracts at different historical positions.

A later identity event invalidates the older policy and trust claims because their dependency contracts include identity continuity.

A newer continuity proof and authority proof remain current.

Before renewal, composition fails with explicit reasons:

```text
policy: stale via identity@15
trust: stale via identity@15
```

After issuing **new** policy/trust claims at the current basis, composition succeeds.

Then an intentionally cheap trust reading with

```text
complete_for_claim = false
```

is composed with the valid continuity proof.

The composition fails explicitly as:

```text
trust-cheap: incomplete-for-claim
```

rather than silently treating the cheap projection as a full trust basis.

That is exactly the evidence-substitution failure Hestia's `ReadBasis`, evidence classes, and occupancy basis are already trying to prevent in separate places.

---

# 8. The result is not necessarily deny

A basis failure should not automatically be confused with a substantive negative verdict.

For example:

- `policy compatible = false` is a claim result;
- `policy proof stale` means the system lacks a current basis for deciding that claim;
- `policy proof incomplete` means the supplied projection cannot establish it;
- `policy evaluator unavailable` is an infrastructure condition.

Those should not collapse to the same boolean.

This fits the broader Hestia three-verdict / escalation posture.

At the meta-level, a relevance contract can support statuses like:

```text
basis-current
basis-needs-refresh
basis-needs-expansion
basis-unverifiable
```

while leaving the domain evaluator to determine the substantive claim.

That keeps **epistemic failure** separate from **negative evidence**.

---

# 9. Escalation is MRH expansion

Earlier Hestia work treats escalation as finding a sufficiently permitted resolver rather than automatically sending everything to a human.

The relevance-contract analogue is:

> when the current compressed basis is insufficient, identify the least expansion needed to make the question decidable at the required assurance.

Examples:

### Windowed trust projection insufficient

Expand:

```text
2,000-row projection -> full relevant chain traversal
```

### Exclusive-continuation proof uncertain

Expand:

```text
local lineage proof -> authoritative registry/fencing proof
```

### New evaluator adds path-sensitive dependency

Expand:

```text
old proof + delta -> historical backfill to admission anchor
```

### A2 evidence insufficient for stakes

Expand:

```text
A2 receipt -> stronger assurance ceremony / evidence source
```

So the escalation function \(X\) is not merely an error message.

It describes how to **grow the MRH until the evidence basis is sufficient for this relying party's question**.

This is a satisfying connection between MRH and practical governance.

---

# 10. A common interface is safer than a universal receipt

There is a tempting but probably wrong implementation move:

> define one enormous `EvidenceReceipt` containing every possible field for every subsystem.

That would couple unrelated domains and create a rapidly ossifying schema.

A better architectural direction is a small common interface, conceptually:

```text
RelevanceBasis {
    claim_kind()
    subject_ref()
    evaluator_plan_ref()
    basis_state_ref()
    dependency_contract()
    coverage()
    assurance_basis()
    lineage_refs()
    escalation_path()
}
```

Then:

- `ReadBasis` can implement/project part of it;
- an AssuranceReceipt can expose its assurance/state basis through it;
- a trust proof can expose its evaluator/dependencies;
- a continuity proof can expose lineage/registry basis;
- an occupancy proof can expose its occupancy/evidence basis.

The existing domain-specific objects remain authoritative for their own semantics.

The common layer enables generic composition, freshness checking, UI disclosure, and escalation.

---

# 11. The contract itself must be provenance-bearing

Phase 15's warning applies recursively.

A relevance contract cannot simply be an unsigned explanation attached to a result.

Load-bearing parts must be bound to the result/proof where appropriate:

- evaluator/plan hash;
- dependency-contract version/hash;
- basis state/root/position;
- claim kind/subject;
- coverage statement;
- assurance references.

Otherwise an intermediary could pair a strong result with a conveniently weaker description of how it was obtained.

So:

\[
\boxed{
\text{the relevance boundary is part of the evidence, not commentary about the evidence.}
}

---

# 12. MRH as a contract, not a radius

This arc began with MRH described mostly as a contextual boundary around an entity.

The accumulated formalization now supports a more general operational reading:

\[
\boxed{
\operatorname{MRH}_R(Q)
=\text{the smallest justified contract of distinctions that must be retained to evaluate }Q\text{ for }R.
}

"Smallest" remains aspirational; exact minimization may be difficult.

The important word is **justified**.

The contract says:

- what distinctions were retained;
- why they are sufficient;
- what was omitted;
- what future change matters;
- when the boundary must expand.

Spatial radius, graph depth, temporal window, evidence class, policy dimension, and lineage depth become different coordinates of the same general idea rather than competing definitions of MRH.

---

## Phase-18 verdict

**PASS — the parallel governance work appears to contain multiple partial implementations of one relevance-basis abstraction.**

What survived:

1. A small common relevance-contract interface could unify basis/completeness semantics without replacing domain-specific proofs.
2. Evidence fetch and semantic invalidation should derive from the same dependency contract.
3. Dependencies need temporal semantics such as current-state versus path-sensitive.
4. `complete` must be claim-relative, not an unqualified boolean assertion.
5. Composition should preserve typed assurance rather than averaging heterogeneous evidence strength.
6. Stale/incomplete/unverifiable basis is epistemically different from a substantive negative verdict.
7. Escalation can be interpreted as deliberate MRH expansion until the evidence basis becomes sufficient.
8. The relevance contract itself must be bound/provenanced as part of the evidence.

What did not happen:

- no proposal to replace `ReadBasis`, AssuranceReceipt, LCT lineage, or domain receipts with one mega-schema;
- no universal minimization algorithm for MRH;
- no normative Hestia/Web4 API change;
- no physics result.

---

# The convergence point

The parallel strands now line up unusually well:

\[
\text{Synchronism MRH}
\longleftrightarrow
\text{predictive/relevance quotient}
\longleftrightarrow
\text{proof dependency horizon}
\longleftrightarrow
\text{Hestia evidence basis}
\longleftrightarrow
\text{Web4 relying-party context}.
\]

The common operation is not "compress everything."

It is:

\[
\boxed{
\text{retain exactly the distinctions you can justify as relevant;
make the boundary inspectable;
expand it when the question or stakes change.}
}

That is probably mature enough now to bring back into the original Synchronism MRH wording as a **candidate refinement**, while keeping the existing canonical text untouched until reviewed.
