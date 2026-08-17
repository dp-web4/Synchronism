# Web4 / LCT Continuity Architecture — Mapping to W3C PROV and ORG

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — architecture mapping / design exploration; **not a Web4 standard change**  
**Scope:** carries the Markov coherence / governance arc into Web4 identity and hub governance

## Why this document exists

The preceding arc converged on a provisional identity architecture:

\[
\text{witness-relative identity}
=
\text{predictive quotient}
+
\text{provenance position}.
\]

An adversarial prior-art review then found that the two sides have strong established neighbors:

- predictive equivalence / coarse-graining: computational mechanics, lumpability, bisimulation, information bottleneck, RG, Markov Stability;
- historical derivation: W3C PROV and related provenance formalisms.

So the next responsible step is not to invent another lineage ontology.

It is to ask:

> **What does the current Web4 LCT model already express, what can standard PROV/ORG express directly, and what Web4-specific evidence semantics remain genuinely necessary?**

This document performs that mapping.

---

# 1. Current LCT position

The current Web4 LCT core specification already contains the right broad pieces:

- a cryptographically bound `lct_id` and `subject`;
- a dynamic MRH of bound, paired, and witnessing relationships;
- policy / capabilities / constraints;
- T3/V3 tensors;
- attestations;
- an optional `lineage` array;
- revocation state;
- explicit rotation semantics.

The current lineage entry is intentionally small:

```json
{
  "parent": "lct:web4:mb32:previous...",
  "reason": "genesis|rotation|fork|upgrade",
  "ts": "2025-09-01T00:00:00Z"
}
```

Rotation creates a new LCT, points its lineage at the parent, permits a temporary overlap window, migrates relationships, and then retires the parent.

This is already a provenance skeleton.

What it does **not** yet fully express is the stronger governance problem exposed by the arc:

- multiple-parent merge semantics;
- recovery / restore semantics;
- exclusive continuation versus intentional branching;
- fencing evidence when exclusivity is claimed;
- the authority / law that authorized a transformation;
- policy state before and after the transition;
- semantic compatibility with a relying party's old admission state;
- assurance level of the continuity claim;
- provenance-aware evidence inheritance through forks and merges.

Those should not automatically become more ad-hoc fields inside `lineage` until existing provenance standards are accounted for.

---

# 2. W3C PROV already supplies most of the historical graph grammar

Primary references:

- W3C PROV-DM: https://www.w3.org/TR/prov-dm/
- W3C PROV-O: https://www.w3.org/TR/prov-o/
- W3C PROV Constraints: https://www.w3.org/TR/prov-constraints/
- W3C Organization Ontology: https://www.w3.org/TR/vocab-org/

PROV's core distinction is useful here:

### `prov:Entity`

A thing with a fixed aspect for the purpose of the provenance account.

### `prov:Activity`

Something that occurs over time and acts upon / transforms entities.

### `prov:Agent`

Something bearing responsibility for an activity or entity.

PROV then connects them with relations such as:

- `prov:wasGeneratedBy`;
- `prov:used`;
- `prov:wasDerivedFrom`;
- `prov:wasRevisionOf`;
- `prov:wasInvalidatedBy`;
- `prov:wasAttributedTo`;
- `prov:wasAssociatedWith`;
- qualified associations carrying a `prov:Plan` and `prov:Role`.

It also has:

- `prov:specializationOf`;
- `prov:alternateOf`;
- generation and invalidation events;
- explicit event-ordering constraints.

This already gives us the graph grammar needed for most LCT lifecycle transitions.

---

# 3. W3C ORG is even closer to governed-society continuity than expected

The W3C Organization Ontology extends PROV for organizational history.

Its `org:ChangeEvent` is a subclass of `prov:Activity` and is explicitly intended for events such as:

- mergers;
- major restructuring;
- other changes after which a resulting organization may warrant a distinct identity / URI.

It permits:

- one or more `org:originalOrganization` inputs;
- one or more `org:resultingOrganization` outputs;
- derivation of the resulting organization from the originals through `prov:wasDerivedFrom`.

This means the generic semantics of:

\[
A\rightarrow B,C
\]

and

\[
B,C\rightarrow D
\]

are already well-covered as organizational change / provenance patterns.

So **fork and merge are not ontology gaps**.

The gap is whether the transition is:

- cryptographically bound;
- authorized by the correct governing law;
- witnessed;
- complete enough to support the relying party's continuity decision;
- exclusive when it claims to be exclusive.

Those are Web4 evidence questions, not missing PROV graph primitives.

---

# 4. Direct LCT → PROV / ORG mapping

| Web4 / LCT concept | PROV / ORG analogue | Notes |
|---|---|---|
| LCT version / token state | `prov:Entity` | Treat each historically fixed LCT realization as an Entity |
| human / AI / society actor | `prov:Agent` | Society can also be `org:Organization` |
| non-agentive Web4 subject | `prov:Entity` | Not every Web4 entity type should be forced into PROV Agent |
| LCT issuance | `prov:Activity` | Generates the LCT entity |
| birth society | `prov:Agent` / `org:Organization` | Associated with issuance activity |
| birth witnesses | agents participating in attestation activities | Their signed attestations remain Web4 evidence |
| binding proof | generated evidence entity | Cryptographic semantics remain Web4-specific |
| rotation | `prov:wasRevisionOf` + derivation activity | New token derived from prior token |
| upgrade | `prov:wasRevisionOf` or typed derivation | Depends whether semantics are revision-like |
| fork | one source Entity → multiple derived Entities | Can be one change activity or several explicit derivations |
| merge | multiple source Entities → one derived Entity | Directly representable in PROV / ORG |
| recovery from backup | derivation from backup + prior lineage state | Requires Web4-specific exclusivity/fencing semantics |
| revocation / retirement | `prov:wasInvalidatedBy` | PROV invalidation says this account considers the entity no longer available after the event |
| policy / constitution version | `prov:Plan` | A plan can be attached to an Agent's qualified association with an Activity |
| actor role in transition | `prov:Role` | Useful for sovereign, witness, arbiter, operator, etc. |
| witnessed transformation | qualified Activity + evidence Entities | Web4 signatures / witness receipts give verifiability |
| membership / role structure | W3C ORG + Web4 role LCTs | Current-state relationship semantics are not primarily provenance |
| current MRH graph | Web4-specific current context graph | Do **not** replace it with PROV; provenance and current relevance answer different questions |
| T3/V3 computation | Activity using evidence entities and plan/version | Result tensor/projection is generated output, not primitive truth |

This strongly suggests **reuse rather than replacement**.

---

# 5. Keep current-state context and historical provenance separate

One architectural mistake would be to make PROV replace the LCT's MRH.

It should not.

They answer different questions.

### MRH / current context

> What relationships and distinctions matter to this entity / witness **now**, for the current interaction?

### Provenance

> By what historical activities and derivations did this present realization come to exist?

A witness may legitimately prune an old relationship from its operational MRH while still needing that relationship preserved in historical provenance.

Conversely, a provenance ancestor may be historically important while currently irrelevant to a low-stakes interaction.

That is exactly the predictive-quotient / provenance orthogonality discovered by the arc.

---

# 6. The genuinely Web4-specific layer: proof-carrying continuity

PROV can say:

\[
B\;\texttt{prov:wasDerivedFrom}\;A.
\]

It does not, by itself, make that assertion cryptographically trustworthy.

Nor does it answer:

- who was authorized to create the derivation;
- under which law;
- which witnesses observed it;
- what assurance level supports the claim;
- whether the old token was actually fenced off;
- whether a hidden concurrent sibling still exists;
- whether the resulting policy remains acceptable to a particular relying party.

That is where Web4 has something distinct to add.

The useful target is therefore:

\[
\boxed{
\text{PROV semantics}
+
\text{Web4 verifiable evidence}
+
\text{MRH-relative relying-party evaluation}
}
\]

—not a new replacement provenance ontology.

---

# 7. Candidate primitive: a witnessed Continuity Event

Rather than expanding each LCT into an ever-growing history document, represent an identity transition as an immutable signed / witnessed event.

Conceptually:

```json
{
  "@type": ["prov:Activity", "web4:ContinuityEvent"],
  "event_id": "evt:web4:...",
  "event_type": "rotation|migration|fork|merge|recovery|amendment|revocation",

  "parents": ["lct:web4:old-a..."],
  "results": ["lct:web4:new-b..."],

  "authorized_by": {
    "society": "lct:web4:society:...",
    "authority": "lct:web4:role:...",
    "law_hash": "sha256:...",
    "decision_ref": "evt:web4:decision:..."
  },

  "continuation_claim": "exclusive|branching|merge|unknown",
  "fencing_evidence": ["receipt:web4:..."],

  "pre_state": {
    "policy_hash": "sha256:...",
    "assurance_profile": "A2"
  },
  "post_state": {
    "policy_hash": "sha256:...",
    "assurance_profile": "A2"
  },

  "witnesses": ["lct:web4:..."],
  "evidence_refs": ["receipt:web4:..."],
  "ts": "2026-08-17T00:00:00Z",
  "sig": "..."
}
```

This is intentionally **illustrative**, not a proposed normative wire schema.

The important design choices are structural:

1. transition is a first-class event / activity;
2. parents and results are arrays, so fork and merge are natural;
3. authorization is evidence attached to the activity;
4. law / plan version is named;
5. exclusivity is a claim with supporting evidence, not inferred from topology;
6. policy / assurance before and after are visible;
7. trust is not embedded as an inherited scalar.

---

# 8. Why `exclusive continuation` is the real semantic gap

PROV can represent an entity being invalidated and a successor being generated.

But a relying party still needs to know whether the claim

> "there is only one legitimate live continuation"

is sufficiently supported.

This is a distributed-systems **fencing** problem, not merely a provenance relation.

A recovery from backup illustrates it.

Suppose old token \(A\) is restored as \(B\).

If the old instance is demonstrably fenced / invalidated, then a relying party may accept:

\[
A\rightarrow B
\]

as exclusive continuation.

If the old instance may still be operating, the correct history is effectively a fork:

\[
A\rightarrow A',B.
\]

The byte-identical restored state and possession of old keys are not enough to distinguish those cases.

So a Web4 continuity event should distinguish:

- `exclusive`;
- `branching`;
- `merge`;
- `unknown`.

And an `exclusive` claim should carry **fencing / exclusivity evidence at a named assurance level**.

No universal assurance threshold belongs in the standard; the relying party decides whether the evidence is sufficient for the stakes.

---

# 9. Do not use `prov:alternateOf` as a shortcut for clones

PROV-O describes `alternateOf` as two entities presenting aspects of the same thing, with examples such as alternate serializations / backup copies.

That is **not automatically the same as an independently created clone**.

A byte-identical new Web4 society with no authorized lineage from the old one can be:

- behaviorally equivalent;
- constitutionally equivalent;
- perhaps a copy of public state;

without being the same historical token or necessarily an alternate representation of the same underlying Web4 presence.

So clone semantics should remain conservative:

\[
\text{structural equality}
\not\Rightarrow
\text{historical continuity}.
\]

Use PROV derivation / alternate relations only when their semantics genuinely apply.

---

# 10. Law maps naturally to `prov:Plan`

A particularly clean mapping appears for governance.

PROV defines `prov:Plan` as an entity describing intended actions / steps, and a qualified `prov:Association` can connect:

- an Activity;
- the responsible Agent;
- the Agent's Role;
- the Plan followed.

That maps closely onto a witnessed Web4 governance action:

\[
\text{activity}
+
\text{authority}
+
\text{role}
+
\text{law version}.
\]

For example, a policy amendment could be represented semantically as:

- amendment = `prov:Activity` / `web4:AmendmentEvent`;
- society/operator/council = `prov:Agent`;
- constitutional role = `prov:Role`;
- signed law version = `prov:Plan`;
- prior policy = used Entity;
- new policy = generated Entity / revision;
- decision / signatures / chain hashes = Web4 evidence.

This is more interoperable than inventing a parallel model for activity/agent/role/law.

---

# 11. Trust tensors should be outputs with provenance, not inherited properties

The fork/merge phase showed that copying scalar trust down branches and adding it at merges can multiply shared ancestral evidence without creating any new evidence.

PROV's activity/entity model suggests the cleaner representation:

\[
\text{evidence entities}
\xrightarrow{\text{evaluation activity + plan}}
\text{trust projection entity}.
\]

For a relying party \(R\):

\[
T_R(v)
=
F_R(
\text{unique evidence},
\text{current state},
\text{provenance},
\text{stakes}
).
\]

A T3/V3 result can therefore carry provenance such as:

- evidence inputs actually consumed;
- derivation-law / model version;
- computation time;
- computing agent / oracle;
- assurance profile;
- resulting tensor.

This matches Hestia's current evidence-first direction much better than treating the scalar as hereditary identity state.

---

# 12. Policy hashes are necessary but not sufficient for trust transfer

A policy hash answers exact equality:

\[
H(C_a)=H(C_b)?
\]

It does not answer semantic compatibility.

The constitutional-drift phase showed that:

- many tiny legitimate amendments can accumulate into material behavioral drift;
- two equal-size textual changes can differ greatly in relevance to a relying party;
- endpoint similarity can hide a historically risky excursion.

So a future relying party needs a second kind of evidence beyond hashes.

A promising non-normative concept is a **Policy Behavior Receipt**.

It might record:

```json
{
  "policy_hash": "sha256:...",
  "probe_set_hash": "sha256:...",
  "probe_results_hash": "sha256:...",
  "engine_version": "...",
  "assurance_profile": "...",
  "witness_refs": ["..."],
  "ts": "..."
}
```

The receipt does **not** say whether the new policy is trustworthy.

It gives a relying party reproducible behavioral evidence from which it can evaluate semantic drift under its own MRH / risk weighting.

This keeps Web4's core principle intact:

> evidence is inspectable; interpretation remains contextual.

---

# 13. Candidate relying-party primitive: Admission Anchor

The continuity arc also suggests that some state belongs to the **relying party**, not the entity.

When hub \(R\) admits entity \(A\), it can locally retain an admission anchor:

\[
\mathcal A_R(A)
=
(
\text{LCT node},
\text{law hash},
\text{assurance profile},
\text{evidence head},
\text{relevant probe set},
\text{risk policy}
).
\]

Later, descendant \(B\) does not simply announce "I am still trusted."

It presents evidence allowing \(R\) to answer:

1. Is there valid witnessed derivation from my anchor to this token?
2. Did the path fork, merge, recover, or amend?
3. If exclusive continuation is claimed, what fencing evidence supports it?
4. What governance plan / law authorized each material transition?
5. Has relevant policy behavior drifted outside *my* tolerance?
6. Does the evidence remain sufficient at the current stakes?

Different relying parties can anchor at different historical points and legitimately reach different answers.

That is MRH applied directly to governance continuity.

---

# 14. Proposed four-layer architecture

The comparison now suggests a clean separation of concerns.

## Layer 1 — identity / current context

**Web4 LCT + MRH**

- cryptographic presence;
- current relationships;
- current capabilities / policy references;
- current contextual state.

## Layer 2 — historical semantics

**W3C PROV + W3C ORG-compatible event graph**

- generation;
- derivation;
- revision;
- invalidation;
- agents / roles / plans;
- organizational change;
- fork / merge topology.

## Layer 3 — verifiable evidence

**Web4 witnessing / receipts / assurance**

- signatures;
- witness chain references;
- authority proofs;
- law hashes;
- fencing evidence;
- assurance profiles;
- evidence dependency identities.

## Layer 4 — contextual evaluation

**MRH / relying-party policy / T3-V3 projection**

- semantic policy relevance;
- trust derivation;
- evidence weighting;
- risk / stakes;
- acceptance / re-admission decisions.

This prevents one layer from being overloaded to answer questions belonging to another.

---

# 15. What Web4 probably should *not* invent

Based on this comparison, avoid independently inventing:

- generic Entity / Activity / Agent provenance classes;
- generic derivation relations;
- generic revision / generation / invalidation vocabulary;
- generic organizational merger / restructuring graph semantics;
- a new definition of fork merely to express one parent with many children;
- a new definition of merge merely to express many parents with one child.

Those have mature prior art.

Web4 can instead specialize / profile them where cryptographic and governance semantics are needed.

---

# 16. What Web4 may genuinely need beyond PROV / ORG

The comparison leaves a much smaller and more defensible list.

### 16.1 Cryptographically bound provenance claims

PROV says what a derivation statement *means*; Web4 needs portable evidence that a specific accountable society/agent actually made and witnessed it.

### 16.2 Governed authorization of transitions

The transition should name the authority, role, law/plan version, and decision evidence that made it admissible.

### 16.3 Assurance-tagged exclusivity / fencing

PROV invalidation is a provenance statement. Web4 needs evidence supporting the stronger operational claim that no competing live continuation remains, at an explicitly named assurance profile.

### 16.4 Evidence dependency preservation

Fork / merge must retain enough receipt identity / dependency information to avoid treating inherited common evidence as independent new evidence.

### 16.5 MRH-relative semantic compatibility

Two governance mechanisms or constitutions may be exact-different but behaviorally equivalent for one relying party and materially different for another.

### 16.6 Relying-party admission anchors

Continuity of *trust* must be evaluated relative to the state on which that party actually relied, not a universal global genesis state.

These are the high-value Web4-specific problems.

---

# 17. A concrete implication for the current LCT `lineage` field

The current optional `lineage` array is not wrong.

It is simply carrying more semantic weight than a three-field parent/reason/time tuple can eventually support.

A good evolution path would be:

1. keep the lightweight `lineage` pointers for easy traversal / compatibility;
2. allow each lineage edge to reference a first-class continuity event / receipt;
3. model event semantics in a PROV/ORG-compatible way;
4. keep Web4-specific authorization, witnessing, assurance, and fencing evidence in the event profile;
5. derive trust separately from evidence rather than encoding trust inheritance in lineage.

Conceptually:

```json
"lineage": [
  {
    "parent": "lct:web4:...",
    "reason": "rotation",
    "event_ref": "evt:web4:continuity:...",
    "ts": "..."
  }
]
```

Again: illustrative only.

---

# 18. One current specification tension worth revisiting later

The LCT spec's strongest design principle says Web4 provides inspectable evidence and does not prescribe universal trust thresholds.

Elsewhere, birth-certificate language still says things such as:

- "Trust Bootstrap: Society's trust transfers to new entity";
- birth-certificate LCT initial trust is "High (inherited)";
- self-issued LCTs have "low trust until witnessed."

Those phrases may be intended as shorthand, but after the governance arc they are semantically stronger than necessary.

A cleaner future phrasing might be:

> a society-issued birth certificate provides **stronger / richer inherited evidence context**, which a relying party may choose to weight highly.

That preserves the design principle that trust itself is not a protocol-declared scalar inheritance.

This is a candidate audit note, not a proposed immediate spec mutation.

---

# 19. Relation to the MRH refinement

This mapping also clarifies why provenance must not itself be aggressively coarse-grained.

At an operational MRH, a witness may safely quotient out:

- component identity;
- implementation mechanism;
- irrelevant environmental state;
- minor behavioral differences.

But provenance distinctions can remain important even when predictive distinctions vanish.

Two current states may satisfy

\[
Q_R(X_A)=Q_R(X_B)
\]

while one is an authorized descendant and the other is a clone.

Thus the combined architecture becomes:

\[
\boxed{
\underbrace{Q_R(X)}_{\text{forget irrelevant predictive distinctions}}
\quad+
\underbrace{P(v)}_{\text{retain relevant historical distinctions}}
}
\]

and the relying party chooses, through MRH, which portions of each are relevant to the current decision.

---

# 20. Architectural verdict

**PASS — the prior-art comparison shrinks the Web4 problem in a useful way.**

What is already supplied by mature standards:

1. provenance entities, activities, and responsible agents;
2. derivation / revision / generation / invalidation;
3. roles and plans associated with activities;
4. organizational change events;
5. multiple inputs / outputs supporting fork and merge histories;
6. temporal consistency constraints.

What remains distinctly valuable for Web4:

1. cryptographic binding of identity / provenance claims;
2. witnessed, governed authorization of transformations;
3. portable evidence receipts;
4. assurance-tagged exclusivity / fencing evidence;
5. provenance-aware trust evidence through branch/merge;
6. MRH-relative behavioral / constitutional compatibility;
7. relying-party-local admission anchors and trust re-evaluation.

The design principle can therefore be compressed to:

\[
\boxed{
\text{reuse provenance semantics; add verifiable governance evidence.}
}
\]

That is a much tighter claim than "Web4 needs a new provenance ontology."

---

# Next move

The next useful step is concrete rather than conceptual:

1. take the current LCT lifecycle (`genesis`, `rotation`, `fork`, `upgrade`, revocation);
2. add the missing governance events (`migration`, `merge`, `recovery`, `amendment`, extinction);
3. write each as a small PROV/ORG-compatible continuity-event profile;
4. define exactly what evidence a relying party can verify for each event;
5. especially stress-test `recovery` and `exclusive continuation`, since those are where ordinary provenance semantics are least sufficient.

That can produce a candidate **Web4 continuity-event profile** without changing the LCT core specification yet.
