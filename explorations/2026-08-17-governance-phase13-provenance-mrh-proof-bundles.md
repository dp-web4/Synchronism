# Governance Arc — Phase 13: Provenance Has an MRH Too

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — provenance / governance formalism result  
**Code:** [`simulations/governance_phase13_provenance_mrh.py`](../simulations/governance_phase13_provenance_mrh.py)  
**Scope:** proof-bundle architecture toy; not new provenance mathematics, not a protocol standard, no physics claim

## Question

Phase 12 ended with an obvious scaling problem.

A long-lived Web4 identity can accumulate a very large provenance DAG:

- rotations;
- migrations;
- forks and merges;
- policy amendments;
- witness receipts;
- trust evidence;
- recoveries;
- unrelated historical activity.

Yet a relying party making one current decision usually does not need all of it.

So ask:

> **What part of the provenance graph is relevant to one current continuity / authority / trust question?**

This is the MRH operation applied to history itself.

---

# Prior art check first

This is **not** a new mathematical idea.

Established provenance work already includes closely related concepts:

- **why-provenance / witnesses** — subsets of input sufficient to explain a query result;
- **dependency provenance** — provenance as dependency/program-slicing analysis;
- **query-specific provenance graphs** — compute only the part of provenance relevant to a why/why-not question;
- **PROV graph abstraction** — hide or group parts of a PROV graph while preserving validity / justified dependencies.

Representative references:

- James Cheney, Amal Ahmed, Umut A. Acar, *Provenance as Dependency Analysis*, arXiv:0708.2173.
- Paolo Missier et al., *ProvAbs: model, policy, and tooling for abstracting PROV graphs*, arXiv:1406.1998; later validity-preserving PROV abstraction work in *Future Generation Computer Systems* 111 (2020).
- Seokki Lee, Sven Köhler, Bertram Ludäscher, Boris Glavic, *Efficiently Computing Provenance Graphs for Queries with Negation*, arXiv:1701.05699.

So the correct novelty posture is:

> **query-specific provenance slicing is known; MRH supplies a cross-domain interpretation and Web4 supplies specific claims/evidence types.**

---

# 1. Full provenance versus claim-specific proof bundle

Let the complete provenance/evidence graph be

\[
G=(V,E).
\]

Let relying party \(R\) ask claim/question

\[
q_R.
\]

A proof bundle is a subgraph / evidence set

\[
G_{q_R}\subseteq G
\]

that is sufficient for the verifier to evaluate that claim.

In the strongest exact case:

\[
V_R(q_R,G_{q_R})
=
V_R(q_R,G).
\]

For an approximate projection such as trust, one may instead require a declared tolerance:

\[
D\left(
V_R(q_R,G_{q_R}),
V_R(q_R,G)
\right)
\le\epsilon_R.
\]

This document deliberately avoids calling the bundle *minimal*.

Finding mathematically minimal provenance/witness sets can be difficult, and the simple toy below computes only a **sufficient backward dependency closure**.

That is enough for the architectural point.

---

# 2. Different claims require different history

The same present token can support several questions.

## Positive descent

> Is current token \(B\) a valid descendant of my admission anchor \(A\)?

Relevant history:

- anchor;
- lineage path \(A\leadsto B\);
- transition events;
- authorization / law references needed to validate them;
- required witnesses.

Unrelated branches and events are irrelevant.

## Exclusive continuation

> Is \(B\) the sole recognized live continuation of \(A\) in authority scope \(S\)?

Needs everything above **plus**:

- fencing / retirement evidence;
- current authoritative registry commitment;
- non-membership / no-sibling evidence in that scope;
- assurance and freshness metadata.

So exclusivity has a strictly larger provenance MRH than positive descent.

## Policy compatibility

> Does the current descendant still satisfy the policy behavior I relied on?

Needs:

- relevant continuity path;
- policy amendments touching the relying party's relevant dimensions;
- behavior/probe receipts needed to evaluate those dimensions.

It need not include every amendment to unrelated policy surfaces.

## Trust projection

> What trust evidence does my current derivation actually consume?

Needs:

- unique evidence receipts actually used;
- dependency identities for deduplication;
- derivation model / law version;
- declared read basis;
- enough lineage/context to interpret those receipts.

It does not automatically need the entire chain.

## Full audit

> Tell me everything that happened.

Here the provenance MRH intentionally expands to the entire graph.

This is the control that prevents “MRH” from becoming a blanket excuse to discard inconvenient history.

---

# 3. Toy construction

The simulation constructs a synthetic provenance graph containing:

- one 80-transition valid lineage from a relying party's admission anchor to the current token;
- authorization, law, and two witness objects per transition;
- eight policy-relevance probe/check objects;
- 60 unique trust-evidence receipts;
- an authenticated-registry/fencing proof for exclusivity;
- 10,000 unrelated provenance event groups by default.

It then creates five claims:

1. descent;
2. exclusive continuation;
3. policy compatibility;
4. trust projection;
5. full audit.

For each claim it computes the simple backward dependency closure.

Again:

> **the closure is sufficient by construction, not asserted to be globally minimal.**

---

# 4. Baseline result

With 10,000 unrelated event groups and a one-million-entry registry, the synthetic full graph contains:

\[
32,596
\]

proof/provenance objects.

The claim-specific closures are:

| claim | proof objects | fraction of full graph |
|---|---:|---:|
| descent | 491 | 1.506% |
| exclusive continuation | 515 | 1.580% |
| policy compatibility | 508 | 1.558% |
| trust projection | 554 | 1.700% |
| full audit | 32,596 | 100% |

The percentages are not empirical Web4 performance measurements.

They merely demonstrate the logical point:

> a claim can depend on a small historical cut even when the complete provenance record is much larger.

---

# 5. Irrelevant history can grow without enlarging a targeted proof

Sweep the number of unrelated synthetic event groups:

| unrelated event groups | full graph | descent | exclusive | compatibility | trust |
|---:|---:|---:|---:|---:|---:|
| 0 | 596 | 491 | 515 | 508 | 554 |
| 100 | 916 | 491 | 515 | 508 | 554 |
| 1,000 | 3,796 | 491 | 515 | 508 | 554 |
| 10,000 | 32,596 | 491 | 515 | 508 | 554 |
| 50,000 | 160,596 | 491 | 515 | 508 | 554 |

That constancy is deliberately built into the dependency structure.

It is not a discovery about real governance systems.

But it gives the intended MRH principle operational form:

\[
\boxed{
\text{history may grow while the evidence relevant to one present claim remains bounded.}
}
\]

The full audit correctly grows with history.

---

# 6. Negative evidence changes the scaling problem

Exclusive continuation contains the difficult claim:

> there is no other recognized active sibling in registry scope \(S\).

Naively, one might imagine transmitting the entire active registry to prove that.

An authenticated dictionary / Merkle-like registry can instead commit to the whole state with one root and provide a logarithmic membership/non-membership path.

The toy therefore gives the exclusivity proof

\[
\lceil\log_2 N\rceil
\]

registry-path proof objects for \(N\) registry entries.

Result:

| registry entries | path objects | exclusive proof objects |
|---:|---:|---:|
| 10 | 4 | 499 |
| 100 | 7 | 502 |
| 1,000 | 10 | 505 |
| 10,000 | 14 | 509 |
| 1,000,000 | 20 | 515 |
| 1,000,000,000 | 30 | 525 |

This is a schematic authenticated-map model, **not a cryptographic proof or prescribed Web4 accumulator**.

The architectural lesson is more important:

> **a cryptographic commitment can let a witness quotient out a large body of irrelevant registry detail while retaining verifiable evidence for the exact fact it needs.**

That is very close to MRH in engineering form.

---

# 7. Proof-carrying quotient

The earlier Synchronism formulation was:

\[
X\sim_R X'
\]

when the distinction between \(X\) and \(X'\) is irrelevant to witness/task \(R\).

For provenance, we can say something analogous without discarding the authoritative record itself.

The complete history remains intact:

\[
G.
\]

But the relying party receives a claim-specific projection:

\[
G
\xrightarrow{Q_{q_R}}
G_{q_R}.
\]

The critical difference from lossy summarization is that the projection can carry **proofs / commitments** tying it back to the authoritative state.

So a useful phrase is:

> **proof-carrying quotient** — a reduced relevant view plus verifiable commitments showing what larger state it is a view of.

This is descriptive language, not a claim of a novel cryptographic primitive.

---

# 8. Hestia is already doing a primitive version of this

A current Hestia trust-derivation change provides a concrete architectural parallel.

The trust code now declares:

- `DERIVATION_EVENT_TYPES` — exactly the event families the model reads;
- `DERIVATION_KEYS` — exactly the event-data keys the model reads.

The reader filters/prunes the witness chain to those declared evidence needs rather than retaining every field of every row.

Just as important, `DashboardSnapshot` carries a `ReadBasis` describing:

- whether the read is `windowed-projection` or `full-traversal`;
- the window when bounded;
- whether the basis is complete;
- a note telling the relying party when it may need to escalate.

The default is deliberately fail-closed: an unstated basis is treated as compressed / incomplete rather than silently implied complete.

That is already extremely close to the architecture being derived here:

\[
\boxed{
\text{model-declared evidence needs}
+
\text{explicit basis/completeness metadata}.
}
\]

So this phase is not asking Hestia to reverse direction.

It suggests generalizing a pattern already present in trust derivation.

---

# 9. Generalize `ReadBasis` into a proof-basis concept

A future continuity / trust / policy proof could carry something conceptually like:

```json
{
  "basis": {
    "claim": "exclusive-continuation",
    "mode": "dependency-slice",
    "anchor": "lct:web4:...",
    "current": "lct:web4:...",
    "scope": "lct:web4:society:...",
    "complete_for_claim": true,
    "full_history_complete": false,
    "commitments": ["registry-root:...", "chain-head:..."],
    "omission_policy": "...",
    "generated_at": "..."
  }
}
```

Important distinction:

### `complete_for_claim`

The producer claims that the bundle includes everything its declared verification algorithm requires for this specific proposition.

### `full_history_complete`

The bundle contains the entire provenance history.

Those are very different properties.

A 500-object proof bundle can be complete for a descent claim while obviously not being a complete 160,000-object historical audit.

This vocabulary could prevent a common evidence substitution error:

> presenting a compressed answer as though it were a complete historical one.

---

# 10. Provenance MRH is question-relative, not just depth-relative

The current LCT MRH specification contains a `horizon_depth` expressed in relationship hops.

The arc suggests that provenance relevance cannot be captured by hop count alone.

Two facts 80 transitions back may differ:

- one is the relying party's admission anchor and is load-bearing;
- another is an unrelated witness event and irrelevant.

Likewise a very recent event may be irrelevant to a trust query if the current model does not read its type or fields.

So a provenance MRH is better indexed by:

\[
(
\text{claim},
\text{anchor},
\text{stakes},
\text{model},
\epsilon,
\text{assurance}
)
\]

than by chronological depth alone.

---

# 11. Exact claims and approximate projections must be separated

Some claims should usually be verified exactly relative to their declared scope:

- valid signature;
- valid lineage path;
- required quorum present;
- registry non-membership proof valid;
- pre-state law authorized the event.

Others are naturally contextual projections:

- trust score;
- behavioral compatibility;
- risk;
- relevance.

For the first class, omitted required evidence should fail verification.

For the second, compression may be acceptable provided the basis/tolerance is declared.

Do not use one generic `confidence` number to blur the difference.

---

# 12. Proof bundles should be generated by the evaluator's declared dependencies

The safest architecture is not:

> "someone decided these fields looked relevant."

It is closer to the seam Hestia has already begun building:

> **the evaluation model declares the event types / fields / proof classes it consumes, and the evidence reader derives its slice from that declaration.**

Then a test can assert:

- every dependency the model reads is present;
- fields the model does not read may be pruned;
- adding a newly read field without updating the declaration fails tests;
- the proof basis states whether history coverage is complete or bounded.

That turns relevance from informal judgment into an inspectable contract between evaluator and evidence layer.

---

# 13. Different proof bundles can share one immutable provenance source

This avoids a false choice:

### Option A

Always transmit / process the entire sacred chain.

### Option B

Throw old history away for efficiency.

A better structure is:

\[
\boxed{
\text{retain authoritative full provenance; derive purpose-specific verifiable views.}
}
\]

The full history remains available for:

- high-stakes escalation;
- audit;
- model change;
- dispute;
- forensic analysis.

Cheap routine decisions operate over declared proof slices.

That is precisely how MRH should behave: relevance limits the current computation, not reality.

---

# 14. Another fractal step

The arc began by asking how a physical entity can ignore distant environmental detail.

Now the same pattern appears in its own history.

At one scale:

\[
\text{microstate detail}
\rightarrow
\text{macro predictive state}.
\]

At another:

\[
\text{full environmental state}
\rightarrow
\text{relevant boundary context}.
\]

Now:

\[
\text{full historical provenance}
\rightarrow
\text{relevant proof subgraph}.
\]

All three are versions of:

> **retain the distinctions required for this prediction / obligation / verification; aggregate or omit the rest while preserving a route to escalation.**

This is making the word **Relevancy** in MRH increasingly literal.

---

# 15. The important asymmetry: the quotient view is disposable, provenance is not

Predictive coarse-graining can often throw information away permanently for the local computation.

Governance provenance should be more conservative.

A proof bundle is a **view**, not the authoritative history.

Why?

Because tomorrow's question can have a different MRH.

An event irrelevant to today's trust projection may become decisive in:

- a later dispute;
- a sanctions review;
- a fork/merge analysis;
- a changed trust model;
- a different relying party's risk context.

So:

\[
Q_{q_1}(G)\neq Q_{q_2}(G)
\]

and neither should destroy \(G\).

That gives governance a strong version of the underlying MRH discipline:

> **relevance is a property of the witness/question, not a declaration that omitted reality never mattered.**

---

## Phase-13 verdict

**PASS as architecture synthesis; no novelty claim.**

What survived:

1. Provenance itself has a question-relative MRH.
2. A full provenance DAG and a sufficient proof bundle are different objects.
3. Positive descent, exclusive continuation, policy compatibility, trust projection, and full audit require different proof subgraphs.
4. Negative/exclusivity evidence can use authenticated commitments to avoid transmitting irrelevant registry bulk.
5. A proof bundle should declare its basis and whether it is complete **for the claim**, separately from whether it contains full history.
6. Evaluation models should declare their evidence dependencies so slices are inspectable and testable rather than hand-selected.
7. The authoritative provenance source should remain intact; MRH produces disposable/query-specific views.
8. Hestia's current `DERIVATION_EVENT_TYPES` / `DERIVATION_KEYS` + `ReadBasis` already instantiate an early form of this pattern.

What did not happen:

- no new provenance-slicing mathematics;
- no proof that the synthetic closure is minimal;
- no prescribed Merkle/accumulator structure;
- no claim that the toy object counts predict real bandwidth or runtime;
- no protocol mutation;
- no physics result.

The compact result is:

\[
\boxed{
\text{keep the full history; carry only the proof horizon required by the present question.}
}
\]

---

# Next move

The remaining architectural question is now less about representation and more about **composition**.

A real relying-party decision may simultaneously require:

- continuity proof;
- exclusivity proof;
- policy-compatibility proof;
- trust evidence;
- authority / delegation proof.

Each has its own relevant slice.

The next question is:

> **When several MRH-specific proof bundles are composed for one action, how do we prevent duplicated evidence, inconsistent anchors, incompatible freshness windows, or one weak proof basis silently weakening the whole decision?**

That points toward a typed **evidence manifest** / proof-bundle composition rule:

\[
B_{action}
=
B_{continuity}
\cup
B_{policy}
\cup
B_{trust}
\cup
B_{authority},
\]

with explicit shared anchors, assurance levels, freshness, completeness, and dependency identity.

That is where I would continue.
