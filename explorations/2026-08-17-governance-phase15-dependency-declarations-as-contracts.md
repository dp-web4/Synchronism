# Governance Arc — Phase 15: Dependency Declarations Must Be Contracts, Not Claims

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — adversarial evidence-basis result  
**Code:** [`simulations/governance_phase15_dependency_contract.py`](../simulations/governance_phase15_dependency_contract.py)  
**Scope:** evaluator/evidence architecture toy; no protocol standard, security proof, or physics claim

## The problem after Phase 14

Phase 14 proposed dependency-relative semantic freshness:

\[
\operatorname{Current}(B,t)
\iff
\nexists e\in(t_B,t]
:\operatorname{type}(e)\in D(B),
\]

where \(D(B)\) is the set of changes that can invalidate proof bundle \(B\).

That creates an immediate attack surface:

> **Who gets to say what belongs in \(D(B)\)?**

If the proof producer simply declares a conveniently small dependency set, then relevant changes can occur without invalidating the proof.

The MRH becomes an assertion rather than evidence.

So the next principle is:

\[
\boxed{
\text{the evaluator must not be trusted to merely describe its own dependencies.}
}
\]

---

# 1. A toy hidden dependency

Use a deterministic evaluator:

```python
score = base_score + 0.1 * witness_quality

if mode == "rare":
    score += rare_override
```

Suppose the declared dependency set is mistakenly:

```text
base_score
witness_quality
mode
```

and omits:

```text
rare_override
```

On ordinary execution, the declaration appears correct.

On the rare branch, it is wrong.

This is the entire problem in miniature.

---

# 2. Runtime tracing on one ordinary execution is insufficient

Run the evaluator with:

```text
mode = normal
```

Observed reads:

```text
base_score
mode
witness_quality
```

The hidden dependency is invisible because that branch did not execute.

Run with:

```text
mode = rare
```

and the read set becomes:

```text
base_score
mode
rare_override
witness_quality
```

So:

\[
\boxed{
\text{observed dependencies in one execution}
\neq
\text{all possible dependencies of the evaluator}.
}
\]

This is standard control-flow / program-analysis territory, not a novel result.

---

# 3. Random mutation testing has a rare-branch problem

Suppose an undeclared dependency matters only in fraction \(p\) of test contexts.

After \(n\) independent random tests, probability of exercising it at least once is:

\[
P_{detect}
=
1-(1-p)^n.
\]

The toy gives:

| hidden-dependency activation rate | 10 tests | 100 tests | 1,000 tests | 10,000 tests |
|---:|---:|---:|---:|---:|
| 1% | 9.56% | 63.40% | 99.996% | ~100% |
| 0.1% | 1.00% | 9.52% | 63.23% | 99.995% |
| 0.01% | 0.10% | 1.00% | 9.52% | 63.21% |

For probability \(\gamma\) of exercising the hidden branch at least once:

\[
n
\ge
\frac{\ln(1-\gamma)}{\ln(1-p)}.
\]

Required trials:

| activation rate | 95% detection | 99% detection |
|---:|---:|---:|
| 10% | 29 | 44 |
| 1% | 299 | 459 |
| 0.1% | 2,995 | 4,603 |
| 0.01% | 29,956 | 46,050 |

So mutation testing is valuable, but absence of a discovered dependency is not strong evidence of absence when the dependency activates rarely.

---

# 4. Make the declaration executable

A stronger design is to place an evidence-access boundary around the evaluator.

The evaluator declares:

\[
D_{declared}.
\]

The evidence view then permits reads **only** from that set.

Conceptually:

```text
read declared field      -> allowed
read undeclared field    -> hard failure
```

In the toy:

### Normal path with incomplete declaration

Passes, because the evaluator reads only declared fields.

### Rare path with incomplete declaration

Fails closed:

```text
evaluator attempted undeclared read: rare_override
```

### Rare path after correcting declaration

Passes.

This changes the dependency manifest from documentation into an enforceable interface contract.

---

# 5. Why this is stronger

Without access enforcement:

\[
D_{declared}
\subsetneq
D_{actual}
\]

can silently produce a wrong but plausible answer.

With access enforcement, the same mismatch produces:

\[
\text{evaluation failure}.
\]

That is preferable for a load-bearing governance computation.

The evaluator cannot accidentally consume evidence that the basis omitted.

This does **not** prove the declaration is complete over all possible paths.

But it guarantees an important local property:

> **if an execution succeeds, every evidence field it actually consumed came through the declared interface.**

That is much stronger than a comment or independently maintained list.

---

# 6. The declaration should be generated from the evaluator boundary where possible

Several assurance methods can reinforce one another.

### Static over-approximation

Analyze all possible reads / dependencies in the deterministic evaluator implementation.

Strength:
- can see unexecuted branches.

Weakness:
- indirect access, dynamic dispatch, generated code, reflection, external calls, and opaque models can make exact analysis difficult.

### Runtime access tracing

Record what one actual execution consumed.

Strength:
- exact for that execution.

Weakness:
- under-approximates unexecuted branches.

### Mutation / sensitivity testing

Perturb candidate inputs and observe output changes.

Strength:
- catches semantic dependencies the implementation inventory may have misunderstood.

Weakness:
- incomplete, especially for rare paths / interactions.

### Access enforcement

Disallow undeclared reads.

Strength:
- converts declaration mismatch into visible failure rather than silent wrong output.

Weakness:
- still requires the declaration to be conservatively broad enough for all legitimate execution paths.

A high-assurance posture combines these rather than pretending one is sufficient.

---

# 7. Control dependencies matter too

Suppose `rare_override` itself is not read on the normal path.

The normal result still depends on the fact that:

```text
mode != rare
```

So the discriminator `mode` is load-bearing.

If `mode` changes, the proof must be recomputed because the evaluator may enter a different dependency region.

This means a useful dependency system must account not merely for direct data reads but also for **control dependencies**.

In general:

\[
D(B)
\supseteq
D_{data}
\cup
D_{control}.
\]

This is another reason simple field grepping is not a complete foundation.

---

# 8. Hidden evaluator state is also a dependency

The evidence fields are only one side.

The result can also depend on:

- evaluator code version;
- model weights / model version;
- prompt / system policy;
- configuration constants;
- external dictionary / ontology version;
- trust-law RDF;
- random seed / nondeterministic sampling regime;
- external service responses.

So a reproducible proof basis needs to commit not just to input evidence but to the **evaluation plan**.

Conceptually:

\[
Result
=F_{	heta}(Evidence),
\]

where \(\theta\) is itself provenance/evidence.

A result computed under evaluator \(F_{\theta_1}\) should not silently be presented as though it came from \(F_{\theta_2}\).

This maps naturally to the earlier PROV idea of a `prov:Plan` attached to the evaluation Activity.

---

# 9. Evidence basis and invalidation basis should derive from the same contract

Phase 13 asked:

> which historical evidence does this model need to read?

Phase 14 asked:

> which future changes invalidate the resulting proof?

Those are two sides of one dependency specification.

If evaluator \(F\) declares evidence dependency classes

\[
D_F,
\]

then the architecture should try to derive both:

### Fetch projection

\[
G_F=\operatorname{Project}(G,D_F)
\]

### Invalidation rule

\[
\operatorname{Invalidate}(B_F,e)
\iff
\operatorname{Touches}(e,D_F).
\]

Maintaining unrelated hand-written lists for fetch and invalidation would create exactly the drift this exercise is trying to eliminate.

---

# 10. This directly resonates with a real Hestia failure mode that was already caught

A recent Hestia trust-derivation change documents a concrete version of this problem.

An earlier inventory believed derivation used eleven fields because it counted only direct `.get(...)` calls.

It missed reads routed through the `entry_str(...)` helper.

A full walk found **28** event-data keys instead.

The commit explicitly notes that a projection omitting one of those fields would not necessarily fail — it could derive a wrong trust number that looked valid.

That is almost exactly the Phase-15 toy, except real.

The current code responds by placing `DERIVATION_KEYS` and `DERIVATION_EVENT_TYPES` beside the derivation code as load-bearing declarations and pruning through that list.

Phase 15 suggests the logical next hardening:

> eventually make the TrustModel/evaluator boundary enforce those declared reads so declaration and consumption cannot silently diverge.

This is architectural extrapolation, not a claim that Hestia currently has that enforcement.

---

# 11. A candidate evaluator manifest

A future evaluator / trust-law artifact could conceptually expose:

```json
{
  "evaluator_id": "sha256:...",
  "evaluator_version": "...",
  "plan_hash": "sha256:...",

  "input_contract": {
    "event_types": ["..."],
    "fields": ["..."],
    "external_sources": ["..."]
  },

  "invalidation_contract": {
    "event_types": ["..."],
    "ontology_versions": ["..."],
    "policy_surfaces": ["..."]
  },

  "analysis_basis": {
    "static": "...",
    "runtime_trace": "...",
    "mutation_suite": "..."
  }
}
```

Illustrative only.

The important idea is that the model version and dependency contract are cryptographically bound together.

---

# 12. Assurance should govern how conservative the dependency horizon is

A low-stakes relying party might accept:

- runtime-traced dependencies;
- a bounded evidence window;
- coarse event-class invalidation.

A higher-stakes party might require:

- evaluator code / plan hash;
- conservative static dependency over-approximation;
- access enforcement;
- mutation/sensitivity coverage;
- full chain escalation where dependency completeness is uncertain.

So the dependency horizon itself has assurance.

Again, Web4 should expose the basis; the relying party decides what is sufficient for the act.

---

# 13. The MRH cannot be self-serving

This is a more general philosophical correction to the arc.

MRH says a witness can discard distinctions that are irrelevant to its current question.

But the witness cannot simply declare inconvenient distinctions irrelevant and call the problem solved.

Relevance needs an external criterion:

- predictive sufficiency;
- verifier equivalence;
- evaluator dependency;
- task performance;
- tolerated error;
- independent evidence.

So:

\[
\boxed{
\text{MRH is justified compression, not permission to ignore.}
}
\]

That applies equally to physics modeling and governance evidence.

---

## Phase-15 verdict

**PASS — and it adds a necessary epistemic constraint to the MRH architecture.**

What survived:

1. Dependency-relative freshness is only as trustworthy as the dependency declaration.
2. Single-run tracing under-approximates dependencies on unexecuted branches.
3. Random mutation testing can require tens of thousands of trials for rare dependencies and cannot prove completeness by itself.
4. An evidence-access boundary can make undeclared reads fail closed instead of silently producing wrong results.
5. Static analysis, runtime tracing, mutation testing, and access enforcement provide complementary assurance.
6. Control dependencies and evaluator/model/config versions belong in the dependency basis.
7. Evidence fetching and proof invalidation should derive from the same dependency contract.
8. The evaluator version / plan and dependency manifest should be cryptographically bound together.
9. Hestia's recent 11-versus-28-field correction is a concrete example of why this matters.
10. MRH requires a testable justification for discarded distinctions; it cannot be self-declared convenience.

What did not happen:

- no claim that exact static dependency analysis is always possible;
- no claim that access tracing alone proves semantic sufficiency;
- no claim that random testing proves absence of hidden dependencies;
- no normative evaluator-manifest schema;
- no security proof;
- no physics result.

The compact result is:

\[
\boxed{
\text{a relevance boundary becomes trustworthy when crossing it is mechanically detectable.}
}
\]

---

# Natural next move

The next question is now evaluator evolution.

If a trust/policy/continuity evaluator changes from

\[
F_{\theta_1}
\rightarrow
F_{\theta_2},
\]

then:

- its evidence dependencies can change;
- its invalidation horizon can change;
- old cached proof bundles may no longer be valid under the new evaluator;
- the same historical evidence can produce a different projection.

That means the **evaluator itself behaves like a replaceable component** from the earlier emergence arc.

The obvious next question is:

> **When can two evaluator versions be treated as equivalent for a relying party's purpose, and when does an evaluator upgrade require recomputing historical trust / continuity projections?**

That loops the governance architecture directly back to Phase 9's mechanism-equivalence problem — now with concrete proof/evidence consequences.
