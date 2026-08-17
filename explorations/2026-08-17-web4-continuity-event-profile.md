# Web4 Continuity Event Profile — Proof-Carrying Identity Transitions

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — candidate profile / architecture exploration; **not a normative Web4 proposal**  
**Depends on:** `2026-08-17-web4-lct-prov-continuity-architecture.md`

## Purpose

The prior mapping found that W3C PROV / ORG already provide most generic lineage semantics.

The remaining Web4 problem is narrower:

> **How can a governed identity transition carry enough cryptographic, authorization, witness, and assurance evidence for a relying party to evaluate continuity without treating lineage as automatic trust or authority inheritance?**

This document sketches a candidate event profile around that problem.

The design rule is:

\[
\boxed{
\text{lineage evidence}
\neq
\text{authority inheritance}
\neq
\text{trust inheritance}
}
\]

Each must be established separately.

---

# 1. Continuity is an event, not merely an edge

The current LCT lineage tuple

```json
{
  "parent": "lct:web4:...",
  "reason": "rotation",
  "ts": "..."
}
```

is useful for traversal.

But a governed transition has more structure than an edge.

Represent it conceptually as a first-class Activity:

\[
E_{parents}
\xrightarrow{A_{transition}}
E_{results}.
\]

The Activity can then carry:

- authorization;
- governing law / plan;
- roles;
- state commitments;
- witness evidence;
- assurance;
- continuity classification.

The LCT lineage edge becomes a pointer into this richer event graph.

---

# 2. Candidate envelope

Illustrative only:

```json
{
  "@context": [
    "https://www.w3.org/ns/prov",
    "https://web4.io/contexts/continuity.jsonld"
  ],
  "@type": ["prov:Activity", "web4:ContinuityEvent"],

  "event_id": "evt:web4:...",
  "event_type": "rotation|migration|fork|merge|recovery|amendment|retirement",
  "ts": "2026-08-17T00:00:00Z",

  "parents": ["lct:web4:..."],
  "results": ["lct:web4:..."],

  "authorization": {
    "society": "lct:web4:society:...",
    "actor": "lct:web4:...",
    "role": "lct:web4:role:...",
    "law_ref": "policy:web4:...",
    "law_hash": "sha256:...",
    "decision_ref": "evt:web4:decision:..."
  },

  "continuity": {
    "mode": "exclusive|branching|merge|unknown",
    "scope": "lct:web4:society:...",
    "assurance_profile": "A1|A2|...",
    "fencing_evidence": ["receipt:web4:..."],
    "registry_commitment": "sha256:..."
  },

  "pre_state": {
    "lct_hashes": ["sha256:..."],
    "law_hash": "sha256:...",
    "witness_head": "sha256:..."
  },

  "post_state": {
    "lct_hashes": ["sha256:..."],
    "law_hash": "sha256:...",
    "witness_head": "sha256:..."
  },

  "evidence_refs": ["receipt:web4:..."],
  "witnesses": ["lct:web4:..."],
  "sig": "..."
}
```

The exact field names are not important yet.

The structural separation is.

---

# 3. Event identity should commit to the event, not name it arbitrarily

A continuity event is valuable only if its references cannot be changed after signing.

A reasonable pattern is:

\[
\text{event\_id}
=H(\operatorname{canonical}(event\setminus signatures)).
\]

Signatures / witness attestations then cover the event digest.

The committed material should include at minimum:

- event type;
- parent set;
- result set;
- pre/post state commitments;
- governing-law reference/hash;
- continuity mode and scope;
- evidence references;
- timestamp / ordering material.

This is ordinary cryptographic hygiene, not a novel mechanism.

---

# 4. The authorization rule: the **pre-state law** must authorize the transition

A governance transition must not be able to justify itself using only the law it creates.

For ordinary transitions:

\[
\text{Authorize}(A_{transition})
\leftarrow
C_{pre}.
\]

For constitutional amendment:

\[
C_{old}
\xrightarrow{\text{authorized under }C_{old}}
C_{new}.
\]

This sounds obvious, but making it explicit prevents a subtle circularity:

> a new constitution cannot retroactively be the sole authority that made its own adoption valid.

Genesis is the special case because there is no prior internal constitution; it must name the external/founding authority/evidence under which genesis is being claimed.

So `authorization.law_hash` is normally a commitment to the **pre-event** governing plan.

The post-event law belongs separately in `post_state`.

---

# 5. Event profiles

## 5.1 Genesis

Topology:

\[
\varnothing\rightarrow A.
\]

PROV reading:

- issuance is Activity;
- LCT is generated Entity;
- society / founders are associated Agents;
- founding procedure is the Plan.

Minimum evidence questions:

- what cryptographic binding was created?
- who issued it?
- what witnesses observed the issuance?
- under what founding authority / genesis procedure?
- what registry/state commitment includes the new token?

Continuity mode: not applicable / `genesis`.

---

## 5.2 Rotation

Topology:

\[
A\rightarrow B.
\]

Typical semantic relation:

\[
B\;\texttt{prov:wasRevisionOf}\;A.
\]

Web4-specific evidence:

- authorization under parent/pre-state law;
- proof of control / approved key replacement;
- declared overlap window;
- relationship / delegation migration evidence where applicable;
- parent invalidation / retirement after cutover;
- fencing evidence if `B` claims **exclusive** continuation.

Important:

During the overlap window, both tokens may legitimately be active.

So the continuity state can change over the event lifecycle:

\[
\text{overlap}
\rightarrow
\text{exclusive after fence/retirement}.
\]

An event or follow-up finalization receipt should make that transition inspectable.

---

## 5.3 Migration

Migration changes location / machine / implementation context.

Two cases should not be conflated.

### Same LCT remains bound

The migration is provenance about the hosting realization; no new LCT need be invented merely because a machine changed.

### Binding / key identity changes

Then the operation resembles rotation:

\[
A\rightarrow B
\]

with migration-specific evidence.

For exclusive migration, the old runtime must be fenced strongly enough for the relying party's assurance requirement.

Otherwise the safer classification is `unknown` or `branching`.

---

## 5.4 Fork

Topology:

\[
A\rightarrow B,C,\ldots
\]

Continuity mode:

```text
branching
```

Each child can possess valid ancestry from \(A\).

No child thereby acquires numerical identity with its siblings.

The fork event should explicitly commit to the **known result set** if the society controls the fork atomically.

But a relying party must not infer that the listed children are globally exhaustive unless the event carries evidence for that stronger claim.

A fork therefore preserves ancestry, not singular identity.

---

## 5.5 Merge

Topology:

\[
A,B,\ldots\rightarrow C.
\]

Continuity mode:

```text
merge
```

The result is a new token with multiple provenance parents.

The merge does **not** imply:

\[
T_C=T_A+T_B
\]

for trust, nor

\[
Authority_C=Authority_A\cup Authority_B
\]

for permissions.

Common ancestral evidence must be deduplicated before contextual trust derivation.

Authorities / delegations need explicit transition semantics under the law authorizing the merge.

Thus the event should reference separate evidence / decisions for:

- provenance merge;
- authority mapping;
- obligation transfer;
- trust recomputation.

Lineage alone does none of those automatically.

---

## 5.6 Recovery

Recovery is the most dangerous continuity event because it often looks identical to exclusive continuation while actually creating a fork.

Suppose backup state \(S_A\) is restored as \(B\).

The provenance statement

\[
B\;\texttt{wasDerivedFrom}\;S_A
\]

is straightforward.

The identity question is not.

If old runtime/token \(A\) may still operate, then the safe classification is:

```text
branching or unknown
```

not `exclusive`.

Only when the relying party accepts evidence that the prior live continuation was fenced / invalidated should the restored token be presented as an exclusive continuation.

Recovery therefore forces the central rule:

\[
\boxed{
\text{state restoration}\neq\text{exclusive identity restoration}.
}
\]

---

## 5.7 Amendment

A constitutional / policy amendment may leave the current LCT token unchanged.

The provenance event still matters:

\[
C_{old}ightarrow C_{new}.
\]

Represent the law versions as Entities / Plans and the amendment as an Activity.

Evidence should include:

- old law hash;
- new law hash;
- authority / role;
- decision / vote / escalation evidence as required by old law;
- witnesses;
- optional behavior/probe receipts.

The event proves **valid amendment lineage**.

It does not prove that every prior relying party should transfer trust automatically.

---

## 5.8 Retirement / revocation

A token may be invalidated because it is:

- compromised;
- superseded;
- expired;
- voluntarily retired;
- sanctioned.

Map this to `prov:wasInvalidatedBy` / an invalidation Activity.

If there is a successor, link that separately through derivation.

Do not overload revocation to mean historical erasure; provenance should remain inspectable.

---

# 6. Exclusive continuation is a **negative-knowledge** claim

Positive provenance is easy:

> here is a witnessed path from \(A\) to \(B\).

Exclusive continuity is stronger:

> there is no other legitimate live continuation of \(A\) relevant to this authority context.

That requires evidence of absence.

In an open distributed system, a global proof of nonexistence is generally unavailable.

So exclusivity must be **scoped**.

A meaningful claim is closer to:

> within society / registry / authority scope \(S\), under commitment \(R_t\), token \(B\) is the sole active recognized continuation of \(A\).

Thus a continuity receipt should carry something like:

```json
{
  "mode": "exclusive",
  "scope": "lct:web4:society:...",
  "registry_commitment": "...",
  "fencing_evidence": ["..."],
  "assurance_profile": "A2"
}
```

The relying party decides whether that scoped evidence is sufficient.

There is no honest universal boolean `exclusive=true` detached from scope and assurance.

---

# 7. A registry / accumulator could make exclusivity evidence stronger

The negative-knowledge problem suggests one concrete mechanism direction.

If a society maintains an authoritative append-only registry of active lineage state, it can commit to a root:

\[
R_t=H(\text{registry state at }t).
\]

Then an exclusive-continuation claim can provide evidence such as:

- parent \(A\) is retired / fenced in committed state \(R_t\);
- child \(B\) is active and names \(A\) as parent;
- no other **recognized active child of \(A\)** exists in that committed registry state.

Depending on registry design, the third item could use:

- an authenticated map proof;
- an accumulator / non-membership proof;
- a signed registry query receipt;
- quorum attestations over the registry state.

This still proves only:

> no competing recognized child exists **inside the committed registry scope**.

It cannot prove that no unauthorized clone exists anywhere in the universe.

That limitation is appropriate and should be explicit.

It is also directly MRH-compatible: exclusivity is always evaluated inside a relevant authority horizon.

---

# 8. `continuation mode` is not a trust verdict

Suppose a token presents:

```text
mode = exclusive
assurance = A2
```

That does not mean "trust this token."

It means the relying party has evidence for one narrow proposition:

> at A2 assurance, within the named scope, this token is claimed/proven to be the sole recognized continuation of the parent.

The relying party still evaluates:

- current policy;
- behavior;
- witness quality;
- trust evidence;
- stakes;
- its own admission anchor.

Keep the semantics narrow.

---

# 9. Authority cannot be inferred from ancestry

This is especially important after forks.

If parent \(A\) has delegation \(D\), a fork

\[
A\rightarrow B,C
\]

does not automatically imply:

\[
D(B)=D(C)=D(A).
\]

Doing that could double authority by copying identity history.

Likewise a merge does not automatically sum parent authorities.

Authority transfer must be a separately authorized governance event or explicit term of the continuity activity.

A useful verifier rule is:

> **Never derive capability/delegation inheritance solely from `wasDerivedFrom`.**

Lineage is evidence of descent, not permission.

---

# 10. Obligations also need explicit transfer semantics

The same applies to obligations:

- debts;
- sanctions;
- votes;
- contracts;
- pending appeals;
- service commitments;
- liability.

After a fork, does each child inherit the whole obligation, a partition, or neither?

There is no universal answer.

That is a law / contract question.

So continuity events should allow explicit references to obligation-transition decisions, but should not encode one default metaphysical rule.

Again:

\[
\text{ancestry}\neq\text{obligation inheritance}.
\]

---

# 11. Trust inheritance is even more dangerous

The prior fork/merge toy showed why.

If common ancestral evidence is copied into both branches and then scalar scores are summed on merge, the same evidence is counted repeatedly.

Therefore continuity events should carry **evidence identity / provenance references**, not inherited scalar trust quantities.

Trust is derived after deduplication / dependency handling:

\[
T_R(C)
=F_R\left(
\operatorname{UniqueEvidence}(Ancestors(C)),
\text{current state},
\text{stakes}
\right).
\]

This is a relying-party projection, not lineage metadata.

---

# 12. A verifier pipeline

A continuity verifier can remain modular.

## Stage 1 — structural validity

Check:

- event schema;
- parent/result cardinality consistent with event type;
- canonical IDs / hashes;
- timestamp / ordering structure.

## Stage 2 — cryptographic validity

Check:

- signatures;
- binding proofs;
- witness receipts;
- registry commitments.

## Stage 3 — governance authorization

Check:

- actor / role;
- **pre-state** law / plan;
- decision evidence;
- quorum / escalation requirements as declared by that law.

## Stage 4 — provenance consistency

Check:

- parents exist / are resolvable;
- derivation ordering is sensible;
- invalidated entities are not silently reused contrary to the claimed event semantics;
- fork / merge graph remains explicit.

## Stage 5 — continuity evidence

Classify:

- exclusive;
- branching;
- merge;
- unknown.

For exclusivity, validate scope + fencing / registry evidence at the claimed assurance level.

## Stage 6 — relying-party evaluation

Separately evaluate:

- semantic policy compatibility;
- evidence quality;
- trust projection;
- obligation / authority transfer;
- stakes.

Only Stage 6 decides what the relying party will actually do.

This keeps protocol evidence and trust interpretation cleanly separated.

---

# 13. Event-specific structural rules

A preliminary profile could enforce simple cardinalities.

| event | parents | results | default continuity shape |
|---|---:|---:|---|
| genesis | 0 | ≥1 | genesis |
| rotation | 1 | 1 | exclusive **only after fence/finalization evidence** |
| migration | 1 | 0 or 1 | same-token migration or exclusive/unknown successor |
| fork | ≥1 | ≥2 | branching |
| merge | ≥2 | ≥1 | merge |
| recovery | ≥1 + backup evidence | ≥1 | unknown until fencing established |
| amendment | policy parent(s) | policy result(s) | token identity may remain continuous |
| retirement | ≥1 | 0 | invalidation / no successor asserted by this event |

These are profile conveniences, not universal metaphysical rules.

---

# 14. The special problem of event finality

A continuity transition should not be considered complete merely because intent was recorded.

For example, rotation may involve:

1. intent / authorization;
2. new token creation;
3. relationship migration;
4. old-token fencing;
5. final witnessed completion.

If step 4 fails, the intended exclusive rotation may actually leave a fork.

So high-consequence continuity operations should have an explicit finality model.

Conceptually:

\[
\text{INTENT}
\rightarrow
\text{COMMIT}
\rightarrow
\text{FINALIZED}
\]

with the final state determining the continuity classification actually supportable by evidence.

This matches a lesson already appearing independently in Hestia governance work: a state-changing governance act should not report success until the terminal witness / durable consequence exists.

---

# 15. Recovery is the best adversarial test case

Any candidate profile that handles easy rotation but fails recovery is not done.

A robust recovery scenario should test at least:

1. old node crashes;
2. backup is restored on a new node;
3. old node later reappears;
4. both possess valid historical keys/state;
5. registry/network partitions may delay fencing evidence;
6. relying parties have different assurance requirements.

Desired behavior:

- lineage to the same ancestor remains visible;
- restoration never silently implies singular continuity;
- concurrent instances become a fork / ambiguity;
- exclusive status is earned only after scoped fencing evidence exists;
- prior trust is re-evaluated rather than duplicated automatically.

This is the highest-value future simulation / implementation test.

---

# 16. A subtle MRH result: numerical identity is authority-scoped

The negative-knowledge problem adds another refinement to the MRH arc.

At a low-assurance / local MRH, a relying party may accept:

> "this key-valid descendant is probably the continuation I care about."

At a higher-stakes MRH, it may require:

> "the authoritative society registry shows no competing live sibling and the old instance is fenced at A2+ assurance."

So even **exclusive token identity** is not one context-free observable.

A better statement is:

\[
\operatorname{ExclusiveContinuation}_R(A,B)
\]

where \(R\) includes:

- authority scope;
- assurance profile;
- registry horizon;
- evidence freshness;
- stakes.

This does not make identity arbitrary.

It makes the evidence threshold and scope explicit.

---

# 17. Architectural compression

The continuity-event profile can now be compressed into four independent questions.

### Descent

\[
\text{Did }B\text{ come from }A?
\]

PROV / lineage evidence.

### Exclusivity

\[
\text{Is }B\text{ the only recognized live continuation of }A\text{ in scope }S?
\]

Fencing / registry / assurance evidence.

### Compatibility

\[
\text{Does }B\text{ still satisfy what relying party }R\text{ cared about?}
\]

MRH-relative behavioral / policy evaluation.

### Reliance

\[
\text{Given descent, exclusivity, compatibility, evidence, and stakes, what should }R\text{ permit?}
\]

Trust / authority / obligation decision.

These should not be collapsed into one `same_identity` boolean.

---

## Profile verdict

**Useful enough to carry forward; not ready for normative standardization.**

What the profile adds beyond ordinary PROV semantics:

1. governed authorization under a named pre-state law;
2. cryptographically committed event evidence;
3. explicit continuity classification;
4. scope + assurance for exclusivity claims;
5. fencing / registry evidence for recovery and migration;
6. strict separation of lineage from authority, obligation, and trust inheritance;
7. finality semantics for multi-step continuity operations;
8. an MRH-relative interpretation of numerical identity evidence.

The most important new sentence is probably:

> **Exclusive continuation is a scoped, assurance-backed negative-knowledge claim — not something that follows from a parent pointer.**

---

# Next move

Stress-test the profile with a small split-brain recovery model.

The toy should include:

- authoritative registry state;
- network partition;
- old-node resurrection;
- restore from backup;
- delayed fencing;
- two relying parties with different assurance thresholds.

Measure / classify the intervals during which:

- descent is known;
- exclusivity is false or unknown;
- both branches remain key-valid;
- the registry eventually converges on one recognized continuation;
- different relying parties legitimately make different continuity decisions.

That would turn the most important architectural gap into an operational example.
