# Governance Arc — Phase 17: Claims Have Their Own Lineage

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — claim/evaluator/entity lineage synthesis  
**Code:** [`simulations/governance_phase17_claim_renewal.py`](../simulations/governance_phase17_claim_renewal.py)  
**Scope:** governance/evidence architecture toy; not a protocol standard, complexity proof, or physics claim

## The recursion after Phase 16

Phase 16 separated three questions around evaluator replacement:

- did the entity remain a valid historical continuation?
- did the evaluator remain an authorized/compatible continuation?
- is the old evidence horizon sufficient for the new evaluator?

One more distinction is necessary.

A derived claim such as

> "this society is an exclusive continuation,"  
> "this policy is compatible," or  
> "this trust projection is 0.78"

also has a history.

If the entity evolves and the evaluator evolves, the new claim should not pretend to be numerically the same historical artifact as the old claim merely because it inherited some evidence.

So there are at least three DAGs:

\[
G_E = \text{entity lineage},
\]

\[
G_F = \text{evaluator/plan lineage},
\]

\[
G_C = \text{claim lineage}.
\]

They interact, but they answer different questions.

---

# 1. A claim should be renewed, not inherited

Suppose claim \(C_0\) was derived at state basis \(s_0\) by evaluator \(F_0\):

\[
C_0=F_0(B_0;s_0),
\]

where \(B_0\) is a complete-for-claim proof bundle.

Later:

- the entity is now at descendant \(E_1\);
- the evaluator is now authorized plan \(F_1\);
- the authoritative state is \(s_1\).

The correct operation is not

\[
C_0\mapsto C_0\text{ with a new timestamp}.
\]

It is

\[
\boxed{
C_1=F_1(B_1;s_1)
}
\]

with provenance recording that \(C_1\) was renewed using evidence/proofs descended from or overlapping \(B_0\).

So:

> **evidence can be reused; claims should be re-derived.**

This is the claim-level analogue of "trust is derived, not asserted."

---

# 2. Incremental renewal can avoid full history

The synthetic history contains 100,000 events.

A relying party's anchor is at position 10,000.

The old claim was derived at position 90,000.

The current head is 100,000.

The old evaluator's path-sensitive dependencies are:

\[
D_0=\{identity,policy\}.
\]

The old complete-for-claim proof contains 476 relevant historical objects from anchor to old basis.

Between old basis and current head, only 53 additional `identity`/`policy` events occur.

Thus same-plan renewal can reuse the 476-object old proof and fetch only the 53 relevant delta objects.

The equivalent full rebuild from anchor to current head would inspect/retain 529 relevant objects.

The new material is therefore about 10.0% of the corresponding rebuild in this toy.

The number is illustrative; the structural point is:

\[
\boxed{
\text{renewal cost can scale with relevant delta rather than total history.}
}

This is standard incremental-maintenance territory, not new computation theory.

---

# 3. Dependency contraction makes renewal cheaper

Suppose the replacement evaluator legitimately needs only the `identity` path in the current MRH.

Then the relevant delta contains only six objects.

The older proof may contain more material than the new evaluator needs, but that does not hurt completeness.

The new evaluator can project a narrower evidence horizon from the retained proof and current delta.

This is the proof-renewal version of Phase 16's dependency contraction:

\[
D_1\subset D_0.
\]

A higher-level quotient can make old distinctions unnecessary.

---

# 4. Dependency expansion splits into two very different cases

This phase exposed a distinction that Phase 16 did not yet need.

A newly introduced dependency can be:

1. **current-state dependent**, or
2. **path dependent**.

Those have radically different renewal costs.

## Current-state expansion

Suppose evaluator \(F_1\) newly requires a current authenticated trust-state value.

If the proposition genuinely depends only on the current committed state, the old proof does not need historical trust events to become complete.

The renewal can use:

- the 53 old-dependency delta objects;
- one current authenticated state/commitment proof.

The toy assigns that snapshot proof a schematic cost of 21 objects, for 74 new objects total.

This is still local to the current head.

## Path-sensitive expansion

Now suppose evaluator \(F_1\) newly requires:

> whether any relevant sanction events occurred since the relying party's original admission anchor.

That fact is historical.

The old proof never retained sanctions because \(F_0\) did not depend on them.

So \(C_1\) cannot be renewed merely from \(B_0\) plus post-\(s_0\) delta.

The missing historical dependency must be backfilled from the anchor to the old basis.

In the toy:

- 153 sanction objects occur in the omitted historical span;
- 13 occur in the new delta;
- together with the 53 old-dependency delta objects, 219 new objects must be fetched.

This is about 31.5% of the corresponding expanded-path full rebuild.

The exact percentage is unimportant.

The important distinction is:

\[
\boxed{
\text{new dependency}\neq\text{new current-state dependency}.
}

A newly added **path dependency** can reopen history that an older MRH safely discarded for the old question.

---

# 5. The scaling difference is structural

Hold the post-proof delta fixed at 10,000 events.

With the toy rates, same-plan renewal expects about 60 new relevant old-dependency events regardless of how old the admission anchor is.

But a newly introduced path-sensitive sanction dependency scales with the omitted historical span:

| old history span | same-plan new material | expanded path new material |
|---:|---:|---:|
| 10,000 | 60 | 100 |
| 100,000 | 60 | 280 |
| 1,000,000 | 60 | 2,080 |
| 10,000,000 | 60 | 20,080 |

So there are two different asymptotics:

### Stable dependency contract

\[
\text{renewal cost}\sim O(\text{relevant delta}).
\]

### Newly introduced historical dependency

\[
\text{renewal cost}\sim O(\text{previously omitted relevant history}).
\]

Again, this is not a complexity theorem for Web4. It is the architectural consequence of what the proof horizon retained.

---

# 6. MRH compression is reversible only if the authoritative record remains available

This produces an important clarification.

When a proof horizon quotients out irrelevant history, the history is irrelevant **for that claim at that time**.

It is not necessarily safe to destroy the authoritative record.

A future evaluator or relying party can ask a broader question whose dependency set includes distinctions that old proof bundles discarded.

Then the system needs to expand the MRH again and retrieve those distinctions from the authoritative provenance source.

Thus:

\[
\boxed{
\text{proof compression may be lossy locally while archival provenance remains lossless enough globally.}
}

This aligns with the existing Hestia posture that the witness chain remains authoritative even when cheap projected readings are used for ordinary decisions.

---

# 7. Three lineages, one renewal event

A useful renewal event therefore links three independent histories.

## Entity lineage

\[
E_0\leadsto E_1
\]

Answers:

> which governed historical token is this?

## Evaluator lineage

\[
F_0\leadsto F_1
\]

Answers:

> which authorized evaluation plan is now in force?

## Claim lineage

\[
C_0\leadsto C_1
\]

Answers:

> which previous claim/evidence basis helped produce this newly derived claim?

The new claim should commit to all three where relevant.

A conceptual claim receipt might contain:

```json
{
  "claim_id": "...",
  "subject_entity": "E1",
  "entity_lineage_anchor": "E0",
  "evaluator_plan": "F1",
  "prior_claim": "C0",
  "basis_state": "state-root/current-head",
  "evidence_basis": "...",
  "delta_basis": "...",
  "complete_for_claim": true
}
```

This is descriptive architecture, not a proposed wire format.

---

# 8. A renewed claim is not necessarily a continuation of the old conclusion

Suppose \(C_0\) said:

> trusted for action X.

After new evidence or evaluator evolution, \(C_1\) may say:

> not trusted for action X.

The claim lineage remains perfectly valid.

So claim continuity means provenance continuity of the **evaluation process**, not preservation of the verdict.

This is subtle but important:

\[
\boxed{
\text{claim-lineage continuity}\neq\text{claim-value continuity}.
}

That parallels the earlier distinction:

- entity lineage can remain continuous while constitution drifts;
- evaluator lineage can remain continuous while semantics change;
- claim lineage can remain continuous while the conclusion flips.

Provenance records legitimate transformation. It does not promise sameness of result.

---

# 9. Proof reuse needs a dependency-coverage test

Let old proof \(B_0\) cover dependency contract \(D_0\).

Let new evaluator require \(D_1\).

Partition new dependencies into:

\[
D_1=D_{current}\cup D_{path}.
\]

A renewal can remain local if:

- retained old evidence covers every needed historical/path dependency up to old basis;
- current/delta evidence covers the interval from old basis to new basis;
- newly introduced current-state dependencies can be proven at the current state;
- any newly introduced path dependency is backfilled to the required anchor.

So the useful question is not merely:

\[
D_1\subseteq D_0?
\]

but:

> **Does the union of retained proof + authenticated delta + available backfill cover the temporal semantics of every dependency required by the new plan?**

That is a more complete formulation of proof-horizon compatibility.

---

# 10. Hestia/Web4 implication

The pattern suggests a practical separation for future governance evidence.

A derived trust/policy/continuity answer can be cached or carried as a proof bundle, but it should name:

- the entity/token basis;
- the evaluator/plan basis;
- the authoritative state basis;
- dependency classes and whether each is snapshot- or path-sensitive;
- the historical anchor for path-sensitive claims;
- completeness for the stated claim;
- the prior claim/proof from which incremental renewal occurred, if any.

Then evaluator upgrades can tell the system whether:

- old proof bundles remain sufficient;
- a small delta is enough;
- a current authenticated snapshot must be added;
- historical backfill is required;
- or full re-derivation/escalation is unavoidable.

This is considerably more precise than a universal TTL or "recompute everything" policy.

---

## Phase-17 verdict

**PASS — entity, evaluator, and claim continuity are distinct provenance objects.**

What survived:

1. Claims should be newly derived under the current evaluator and state, not inherited by relabeling old results.
2. Complete old proof bundles can make renewal incremental when dependency semantics remain compatible.
3. Dependency contraction can reduce the future proof horizon.
4. New current-state dependencies can often be satisfied locally at the current head.
5. New path-sensitive dependencies can force historical backfill across distinctions the old MRH had legitimately omitted.
6. Therefore claim-specific compression should not imply destruction of the authoritative provenance source.
7. Claim-lineage continuity does not imply claim-value continuity.

What did not happen:

- no new incremental-computation mathematics;
- no claim that all historical queries can be efficiently backfilled;
- no universal event taxonomy or dependency semantics;
- no protocol standard or physics result.

---

# A compact convergence

The arc now has a recurring operation at several scales:

\[
\boxed{
\text{retain an invariant / sufficient proof,}
\quad
\text{discard lower-MRH distinctions,}
\quad
\text{re-expand when a new question makes them relevant again.}
}
\]

That is true for:

- microstate versus entity;
- components versus relation;
- implementation versus constitution;
- full provenance versus claim proof;
- old evaluator versus replacement evaluator;
- old claim versus renewed claim.

The next question is whether this can be represented as one generic **relevance contract** shared by Hestia/Web4 surfaces, rather than each subsystem inventing its own notion of basis, completeness, dependencies, anchor, and escalation.
