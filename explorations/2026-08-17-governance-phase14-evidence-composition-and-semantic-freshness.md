# Governance Arc — Phase 14: Evidence Composition and Semantic Freshness

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — evidence-composition formalism result  
**Code:** [`simulations/governance_phase14_evidence_composition.py`](../simulations/governance_phase14_evidence_composition.py)  
**Scope:** governance/evidence architecture toy; not a protocol standard, security proof, or physics claim

## Question

Phase 13 ended with several claim-specific proof bundles:

\[
B_{continuity},
B_{policy},
B_{trust},
B_{authority},\ldots
\]

A real action often needs several at once:

\[
B_{action}
=
B_{continuity}
\cup
B_{policy}
\cup
B_{trust}
\cup
B_{authority}.
\]

That creates a new problem.

Even if every bundle was valid when produced, can they safely be composed **now**?

The first instinct is to add timestamps and TTLs.

That turns out to be too weak.

---

# 1. Freshness is not merely age

Suppose a policy proof was generated one second ago.

If the policy changed 100 ms later, the proof may already be stale for a policy-sensitive action.

Now suppose a continuity proof was generated an hour ago.

If the intervening million events were all unrelated telemetry and no continuity state changed, the proof may still describe the relevant continuity fact perfectly well.

So:

\[
\boxed{
\text{recent}\neq\text{current}
}
\]

and

\[
\boxed{
\text{old}\neq\text{stale}.
}
\]

Chronological age matters for some cryptographic / operational reasons, but it is not sufficient to determine semantic freshness.

---

# 2. Dependency-relative validity

Let proof bundle \(B\) declare the event classes / state dimensions on which its claim depends:

\[
D(B).
\]

Let the bundle be generated at event position \(t_B\).

Then a simple semantic-currentness condition is:

\[
\boxed{
\operatorname{Current}(B,t)
\iff
\nexists e\in(t_B,t]
:\operatorname{type}(e)\in D(B)
}
\]

or, more generally, no intervening state change intersects the dependency set the verifier says can alter the proposition.

This is standard dependency / cache-invalidation thinking applied to governance evidence.

The MRH interpretation is useful:

> **events outside the proof's dependency horizon are temporally irrelevant to that claim.**

---

# 3. Toy event distribution

The simulation uses five event classes:

| event class | probability per event |
|---|---:|
| unrelated / noise | 0.90 |
| trust-relevant | 0.04 |
| policy-relevant | 0.03 |
| authority-relevant | 0.02 |
| continuity-relevant | 0.01 |

Illustrative proof dependencies are:

\[
D(B_{continuity})
=\{continuity\},
\]

\[
D(B_{trust})
=\{trust,continuity\},
\]

\[
D(B_{policy})
=\{policy,continuity\},
\]

\[
D(B_{authority})
=\{authority,policy,continuity\}.
\]

The composed action manifest depends on the union:

\[
D(B_{action})
=
\bigcup_i D(B_i)
=
\{trust,policy,authority,continuity\}.
\]

Nothing about those particular dependencies is claimed universal.

The point is to study composition.

---

# 4. Validity horizon follows relevant-event rate

If events are independent with relevant-event probability \(p_B\), then the number of events until invalidation is geometric with mean

\[
\mathbb E[\tau_B]=\frac1{p_B}.
\]

For the toy:

| proof | invalidation hazard | exact mean validity horizon | simulated mean | median | 95th percentile |
|---|---:|---:|---:|---:|---:|
| continuity | 0.01 | 100.000 | 100.302 | 70 | 299 |
| trust | 0.05 | 20.000 | 20.062 | 14 | 59 |
| policy | 0.04 | 25.000 | 24.918 | 17 | 74 |
| authority | 0.06 | 16.667 | 16.630 | 12 | 49 |
| composed action | 0.10 | 10.000 | 10.012 | 7 | 29 |

The exact values are toy parameters.

The structural result is:

\[
\boxed{
\text{composing claims usually expands the dependency set and shortens the common validity horizon.}
}
\]

The action bundle is not as temporally durable as its most durable component.

---

# 5. Composition creates a shared relevance horizon

This is not exactly a conventional "weakest link" rule.

Suppose:

- continuity changes rarely;
- policy changes more often;
- trust evidence changes still differently.

The composed action remains supported only while **every required proposition** remains current.

So its invalidating event set is the union:

\[
D_{action}
=
D_1\cup D_2\cup\cdots.
\]

If event classes overlap across proofs, the hazards are not simply additive, but the same logical rule remains:

> any change that invalidates any load-bearing claim can invalidate the composed decision basis.

This gives an action a **shared temporal MRH** over the state distinctions relevant to all required proofs.

---

# 6. Global state-root invalidation is safe but often unnecessarily expensive

One simple design is:

> any change to the society state produces a new root, therefore every proof tied to the old root must be regenerated.

In the toy, that means every event invalidates every proof:

\[
p=1,
\qquad
\mathbb E[\tau]=1\text{ event}.
\]

Dependency-aware invalidation lets the composed action survive unrelated events:

\[
\mathbb E[\tau_{action}]=10\text{ events}
\]

under the illustrative event mix, while the continuity proof alone survives on average 100.

The efficiency gain is not free.

To do this safely, the proof must declare enough about its dependencies that a verifier can establish that intervening events were noninterfering for that proposition.

That declaration becomes part of the evidence basis.

---

# 7. A short trace shows why TTL alone is insufficient

Generate all four component proofs at event 0.

Then observe:

```text
1 noise
2 noise
3 policy
4 noise
5 trust
6 authority
```

At event 6, every proof is exactly six events old.

Yet under the declared dependencies:

| proof | current? | reason |
|---|---|---|
| continuity | yes | no continuity event occurred |
| trust | no | trust event at 5 |
| policy | no | policy event at 3 |
| authority | no | policy/authority events occurred |
| composed action | no | several required component proofs are stale |

So a TTL rule such as "proofs younger than 10 events are fresh" would accept stale evidence.

A rule such as "every event invalidates every proof" would reject the still-current continuity proof.

Both flatten relevance.

---

# 8. Irrelevant activity should not age semantic evidence

The toy then inserts only unrelated/noise events after proof generation:

| unrelated events | continuity current | trust current | action current |
|---:|---|---|---|
| 0 | yes | yes | yes |
| 10 | yes | yes | yes |
| 100 | yes | yes | yes |
| 1,000 | yes | yes | yes |
| 100,000 | yes | yes | yes |

That is deliberate by construction.

It expresses the intended principle:

\[
\boxed{
\text{irrelevant change should not automatically destroy relevant evidence.}
}
\]

This is the temporal form of the same component-replacement idea that started the arc.

---

# 9. But semantic freshness does not replace absolute-time validity

A proof can remain semantically unchanged while becoming unusable for other reasons:

- signing key expires or is revoked;
- assurance receipt exceeds a relying party's permitted age;
- challenge / nonce freshness window closes;
- external evidence source has a maximum accepted staleness;
- protocol replay limits apply.

So a practical validity predicate is more like:

\[
\operatorname{Valid}_R(B,t)
=
\operatorname{CryptoValid}(B,t)
\land
\operatorname{SemanticallyCurrent}(B,t)
\land
\operatorname{FreshEnough}_R(B,t).
\]

MRH/relevance-aware invalidation complements wall-clock freshness; it does not abolish it.

---

# 10. The Frankenproof problem

Composition creates another failure mode.

A relying party can be shown:

- a continuity proof from state epoch \(s_1\);
- a trust proof from \(s_2\);
- a policy proof from \(s_3\);
- an authority proof from \(s_4\).

Each can be individually valid in isolation.

But there may never have been one system state in which **all four propositions held simultaneously**.

For example:

1. old policy was restrictive when the high trust score was derived;
2. policy later became permissive;
3. an old trust receipt is combined with the new policy proof;
4. authority changed again;
5. the bundle presents the best surviving fact from each epoch.

Call this a **Frankenproof** descriptively: composition of individually genuine evidence into a collectively misleading state claim.

No fraud in the individual receipts is required.

---

# 11. Composition therefore needs state compatibility, not just signature validity

A typed evidence manifest should let the verifier establish one of two things.

### Common-state basis

All component proofs commit to the same sufficiently current society / chain / registry state.

or

### Valid carry-forward

A proof originated at an older state, but every intervening event is demonstrably irrelevant to the dependencies that proof declares.

Formally, for proof generated at state \(s_i\) and action evaluated at current state \(s_*\):

\[
B_i(s_i)
\rightsquigarrow
B_i(s_*)
\]

only if the transition path

\[
s_i\leadsto s_*
\]

contains no event invalidating \(D(B_i)\).

This is stronger than checking timestamp age.

---

# 12. Candidate evidence-manifest fields

A claim-specific bundle might eventually expose something like:

```json
{
  "claim_type": "policy-compatibility",
  "subject": "lct:web4:...",
  "admission_anchor": "lct:web4:...",
  "authority_scope": "lct:web4:society:...",

  "basis_state": {
    "chain_head": "sha256:...",
    "registry_root": "sha256:...",
    "position": 123456
  },

  "dependencies": {
    "event_types": ["law_amended", "policy_changed"],
    "fields": ["law_hash", "decision", "role_lct"]
  },

  "basis": {
    "complete_for_claim": true,
    "full_history_complete": false
  },

  "assurance_profile": "A2",
  "generated_at": "...",
  "not_after": "...",
  "evidence_refs": ["..."]
}
```

Illustrative only.

The core idea is that **dependency declarations travel with the proof**.

---

# 13. Action manifests should compose requirements, not scores

Suppose an action policy requires:

```text
continuity: exclusive-supported, >= A2
policy: compatible, >= A2
authority: valid, >= A2
trust: derivation available, >= A1
```

The composition rule should check each proposition separately.

A strong trust proof cannot upgrade weak continuity evidence.

A strong continuity proof cannot compensate for expired authority.

So:

\[
\boxed{
\text{evidence assurance is typed by claim; it should not collapse into one average confidence score.}
}
\]

The same follows for completeness.

`complete_for_claim=false` on one required exact proposition should not be averaged away by several complete bundles.

---

# 14. Evidence identity must survive composition

Proof bundles can share receipts.

For example, one witnessed event might support:

- continuity;
- trust derivation;
- authority history.

When manifests are unioned, the receipt should remain one evidence object referenced by several claims.

Do not count it three times merely because three proof bundles cite it.

This is the same fork/merge lesson from Phase 10 at the evidence-manifest level:

\[
\boxed{
\text{compose evidence by identity-aware union, not scalar addition.}
}

---

# 15. Hestia's declared derivation dependencies point toward this architecture

Hestia already declares the event types and event-data keys used by its trust derivation and projects the witness-chain read accordingly.

Phase 13 interpreted that as an early claim-specific provenance MRH.

Phase 14 suggests the next logical property:

> a derived trust/evidence result can also declare **which changes would invalidate that result**.

If the declared evidence dependency set and the invalidation dependency set are generated from the same model boundary, then:

- evidence fetching;
- cache invalidation;
- proof composition;
- freshness explanation

all derive from one inspectable contract.

That is architecturally attractive because divergence between those four lists would otherwise be dangerous.

---

# 16. A proof bundle has two horizons

The previous work focused on **evidence breadth**:

> which provenance objects are required?

Phase 14 adds **evidence lifetime**:

> which subsequent transformations can occur before the proof must be refreshed?

So for claim \(q\), one can think provisionally of:

\[
B_q
=
(
G_q,
D_q
),
\]

where:

- \(G_q\) is the relevant historical/evidence slice;
- \(D_q\) is the class of future changes that invalidate it.

This is a nice duality:

### Backward MRH

Which past distinctions must I retain to justify this claim?

### Forward MRH

Which future changes remain irrelevant before I must recompute it?

That is another way of making MRH explicitly temporal.

---

# 17. Relations among proofs

The arc has now repeated its earlier relational move at a new scale.

Initially:

- entities were treated individually;
- then relations among entities became the persistent object;
- then relations among relations.

Now:

- proof bundles are individually valid objects;
- but the **relations among proof bundles** determine whether a composed action basis is coherent.

The relevant relations include:

- same subject;
- compatible anchor;
- compatible authority scope;
- compatible state commitment;
- no invalidating intervening events;
- shared evidence identity;
- assurance requirements.

So validity of the whole is not a sum of validity of the parts.

It is itself a relational property.

That is precisely the pattern that started this exploration.

---

## Phase-14 verdict

**PASS as architecture sharpening; no novelty claim.**

What changed:

1. Evidence age and evidence semantic freshness are different.
2. A proof remains current across changes outside its declared dependency set.
3. Proof composition expands the dependency set and usually shortens the joint validity horizon.
4. Individually valid proofs can form an invalid/misleading composition if their state bases are incompatible.
5. Carry-forward from an older proof state requires evidence that intervening events were irrelevant to that proof's dependencies.
6. Assurance/completeness must remain typed by claim rather than averaged into one confidence value.
7. Shared evidence must be deduplicated by identity during composition.
8. A proof bundle has a backward provenance horizon and a forward invalidation horizon.
9. Relations among proof bundles become the higher-scale object determining action-level coherence.

What did not happen:

- no claim that the toy event distribution represents Hestia/Web4 deployment rates;
- no replacement for TTLs, revocation checks, or cryptographic expiry;
- no normative manifest schema;
- no security proof against adversarial dependency declarations;
- no physics result.

The compact result is:

\[
\boxed{
\text{proof freshness is determined by relevant change, not elapsed time alone.}
}

---

# Next move

The obvious adversarial question is now:

> **Who is allowed to declare the dependency set?**

If a proof producer can omit a load-bearing dependency, it can make stale evidence appear permanently current.

So the next phase should attack the dependency declaration itself.

Possible controls:

- derive dependencies mechanically from the evaluator / policy implementation where possible;
- bind dependency declarations to evaluator/model version;
- test mutation sensitivity: change each declared/undeclared state class and verify whether the claim changes;
- compare dynamic declarations with static dependency analysis;
- let relying parties demand broader conservative dependency classes at higher assurance;
- treat an undeclared basis as weakest / incomplete, as Hestia already does for `ReadBasis`.

That is the next load-bearing seam.
