# Markov Relevancy Horizon as a Relevance Contract

**Status:** ripening / candidate refinement  
**Origin:** Markov coherence/emergence arc, Phases 1–18 (2026-08-17)  
**Scope:** conceptual/formal clarification of MRH; **not yet canonical Synchronism wording**

## Why revisit MRH?

The original MRH intuition was deliberately broad:

> a witness does not need the entire universe to evaluate a local interaction; beyond some context boundary, additional detail becomes irrelevant to the question being asked.

That remains useful, but the Markov/coherence/governance arc exposed several dimensions that a simple spatial or graph-depth picture does not capture.

MRH is also:

- temporal;
- timescale-relative;
- target/question-relative;
- tolerance/stakes-relative;
- transformation-relative;
- provenance-relative;
- capable of expanding again when the question changes.

The refinement below attempts to put those hunches into one operational form.

---

# 1. MRH is not primarily a radius

A spatial radius is one possible implementation of a relevance boundary.

So is:

- graph hop depth;
- temporal window;
- policy dimension;
- evidence class;
- lineage depth;
- retained state variable;
- evaluator dependency;
- assurance threshold.

The more general idea is:

\[
\boxed{
\operatorname{MRH}_R(Q)
=\text{the justified boundary between distinctions that must be retained to evaluate }Q\text{ for witness }R
\text{ and distinctions that may be quotiented out.}
}
\]

Exact minimization is not assumed tractable.

The operative word is **justified**.

---

# 2. Predictive form

For future target \(Y_{t+h}\), retained boundary/context \(B_t\), and exterior \(E_t\), an information-theoretic approximation is:

\[
I(Y_{t+h};E_t\mid B_t)\le\epsilon_R.
\]

Then a spatial predictive MRH can be written:

\[
\boxed{
r_{MRH}(Y;h,\epsilon_R)
=\min_r\left\{
I(Y_{t+h};E_{r,t}\mid B_{r,t})\le\epsilon_R
\right\}.
}
\]

This distinguishes:

- **causal reach** — what can possibly affect the target;
- **relevant predictive reach** — what must actually be retained for the required prediction quality.

Something may lie inside the causal cone yet outside the current MRH.

This is closely related to established predictive-state, local-causal-state, information-bottleneck, and coarse-graining ideas. The mathematical ingredients are not claimed novel.

---

# 3. Temporal form

A coherent entity requires a separation between fast internal change and slower loss of macro-identity:

\[
\tau_{internal}\ll h\ll\tau_{identity}.
\]

Thus entityhood is observation-timescale-relative.

A component, relation, or implementation may change many times inside one macro-identity lifetime.

The Markov arc repeatedly found the same structure:

\[
\text{fast degrees of freedom}
\rightarrow
\text{slow invariants}
\rightarrow
\text{slow invariants among those invariants}.
\]

So an emergent identity can be described provisionally as:

\[
\boxed{
\text{a slow invariant under faster transformations that are internal at the current MRH.}
}
\]

This allows complete component replacement while preserving the larger identity.

---

# 4. Transformation-relative form

Let \(\mathcal T_R\) be the class of transformations witness \(R\) is justified in treating as internal implementation detail for question \(Q\).

Then:

\[
X\sim_{R,Q}X'
\]

when the distinction between \(X\) and \(X'\) does not change the relevant evaluation under that transformation class.

Examples of distinctions that may be quotiented out at progressively higher scales include:

\[
\mathcal T_{state}
\subset
\mathcal T_{component}
\subset
\mathcal T_{relation}
\subset
\mathcal T_{implementation}
\subset
\mathcal T_{constitutional}.
\]

The nesting is conceptual, not asserted universal.

The key result is:

> **a higher MRH does not merely see a larger region; it forgets more distinctions.**

---

# 5. Relations can be entities

The interaction arc corrected an overly simple mapping between coherence and resonance/dissonance.

Strong positive coupling can stabilize an aligned relation.

Strong negative coupling can stabilize an anti-aligned relation.

Both can be highly coherent.

So:

\[
\boxed{
\text{interaction strength}
\neq
\text{interaction valence}
\neq
\text{interaction coherence}.
}
\]

A persistent relation can itself become a metastable entity.

Then higher-scale structure may be:

\[
\text{states}
\rightarrow
\text{persistent relations}
\rightarrow
\text{relations among relations}.
\]

This is one path by which the fractal hierarchy becomes relational rather than merely compositional.

---

# 6. The interior has an MRH too

Earlier MRH language emphasized the exterior boundary:

> how much distant environment can be ignored?

The arc adds an interior question:

> how much microscopic implementation detail can be ignored before prediction of the macro-identity suffers?

So an entity can occupy an informational shell between:

- irrelevant interior detail;
- retained slow relational/invariant state;
- irrelevant exterior detail.

A useful pair of closure conditions is therefore:

\[
I(Y_{t+h};X_t\mid Y_t)\le\epsilon_{interior},
\]

\[
I(Y_{t+h};E_t\mid Y_t,B_t)\le\epsilon_{exterior}.
\]

The first asks whether discarded microdetail still matters.

The second asks whether discarded exterior context still matters.

---

# 7. Provenance must not be quotiented the same way as prediction

Predictive compression intentionally merges histories that do not change the relevant future distribution.

Governance/accountability may still need to distinguish those histories.

A clone and an authorized continuation may behave identically while carrying different:

- authority;
- obligations;
- sanctions;
- delegation lineage;
- evidence dependencies;
- historical identity.

So witness-relative governed identity is better represented with two orthogonal axes:

\[
\boxed{
\text{predictive/relevance quotient}
\quad+\quad
\text{provenance position}.
}
\]

The first forgets distinctions that do not matter to the current prediction/question.

The second preserves enough derivational history to establish which historical token this is.

Established computational mechanics, coarse-graining, bisimulation, provenance models, and related work already cover much of the component mathematics. The value here is the joint systems interpretation, not a claim of new mathematics.

---

# 8. Provenance itself has an MRH

Keeping authoritative provenance does not mean every decision must traverse all of it.

For complete graph/history \(G\) and claim \(Q_R\), define a sufficient claim-specific proof projection:

\[
G_{Q_R}\subseteq G
\]

such that the verifier reaches the same claim result it would from the required authoritative basis.

This is query-specific provenance slicing, an established idea.

The MRH interpretation is:

> **history may grow while the evidence relevant to one present claim remains bounded.**

However, the authoritative record remains important because a future question can expand the MRH and make previously omitted history relevant again.

Thus:

\[
\boxed{
\text{local proof compression may be lossy while archival provenance remains lossless enough globally.}
}

---

# 9. MRH has a forward boundary too

A proof/evaluation basis remains semantically current until a relevant distinction changes.

Let \(D(B)\) be the dependencies of proof bundle \(B\).

Then provisionally:

\[
\operatorname{Current}(B,t)
\iff
\nexists e\in(t_B,t]
:\operatorname{Touches}(e,D(B)).
\]

So:

\[
\boxed{
\text{recent}\neq\text{current},
\qquad
\text{old}\neq\text{stale}.
}
\]

An old continuity proof may remain semantically current through unrelated telemetry.

A one-second-old policy proof may become stale immediately after a policy mutation.

Wall-clock expiry can still matter for cryptographic or operational reasons; semantic freshness is an additional axis.

---

# 10. Dependency declarations must be enforceable

If an evaluator may freely declare its own relevance boundary, it can omit inconvenient dependencies.

So a trustworthy MRH cannot be merely descriptive metadata.

A stronger architecture makes the dependency contract executable:

```text
read declared evidence   -> allowed
read undeclared evidence -> fail closed
```

This does not prove global minimality or completeness of every possible execution path.

It does establish an important local property:

> **a successful evaluation consumed only evidence admitted by its declared relevance contract.**

Static analysis, runtime tracing, mutation testing, and review can provide additional assurance around that contract.

---

# 11. Evaluators themselves are replaceable components

Suppose evaluator \(F_a\) is replaced by \(F_b\).

Three separate questions arise:

### Semantic compatibility

\[
d_R(F_a,F_b)\le\epsilon_R?
\]

### Evidence compatibility

Does the new evaluator require evidence contained in the old proof horizon?

### Historical/authority continuity

Was the evaluator replacement an authorized witnessed transition?

These axes are independent.

A replacement can preserve current outputs while expanding its dependencies.

Then:

\[
\boxed{
\text{same answer}\not\Rightarrow\text{same proof basis}.
}

Conversely, a higher-level invariant may let a replacement evaluator depend on fewer distinctions while preserving the relevant result.

That is a genuine proof-horizon contraction.

---

# 12. Claims have their own lineage

Entity lineage, evaluator lineage, and derived-claim lineage should not be collapsed:

\[
G_E,
\qquad
G_F,
\qquad
G_C.
\]

Evidence from an old claim may be reused, but the current claim should be newly derived under the current evaluator and state basis:

\[
C_1=F_1(B_1;s_1).
\]

Therefore:

\[
\boxed{
\text{evidence persists; plans evolve; claims are re-derived.}
}

Claim-lineage continuity does not imply claim-value continuity.

A perfectly legitimate renewal can reverse the previous conclusion.

---

# 13. A relevance contract

The recurring metadata can be summarized conceptually as:

\[
\mathcal R
=
(Q,S,F,B,D,C,A,L,X),
\]

where:

- \(Q\): question/claim;
- \(S\): subject/context;
- \(F\): evaluator/plan;
- \(B\): authoritative basis state/commitment;
- \(D\): dependency contract, including temporal/path semantics;
- \(C\): coverage/completeness statement;
- \(A\): assurance basis;
- \(L\): relevant lineage references;
- \(X\): escalation/re-expansion path.

This is not proposed as a universal wire schema.

It is a way to make the MRH **inspectable**.

A small common interface can expose these semantics while domain-specific proof objects remain separate.

---

# 14. Escalation is MRH expansion

If the current basis is insufficient, the right response is not necessarily a negative verdict.

The system may instead need to expand its relevance horizon:

```text
windowed projection -> full relevant traversal
local lineage proof -> authoritative registry/fencing proof
old proof + delta -> historical backfill to anchor
A1 evidence -> stronger assurance basis
```

So:

\[
\boxed{
\text{escalation is the deliberate expansion of the MRH until the question becomes decidable at the required assurance.}
}

This interpretation connects abstract relevance directly to practical governance.

---

# 15. Candidate refined statement

A compact candidate wording is:

> **The Markov Relevancy Horizon (MRH) of a witness/question is the justified boundary between distinctions that must be retained for the evaluation and distinctions that may be treated as irrelevant at that scale. The boundary is relative to target, timescale, tolerance/stakes, transformation class, and available evidence. It may contract when lower-level distinctions become predictively redundant and expand when a new question, dependency, or assurance requirement makes previously discarded distinctions relevant again.**

A still shorter form:

> **MRH is a witness-relative contract over which distinctions still matter.**

---

# 16. Novelty posture

The following ingredients have strong established neighbors and should not be claimed as new mathematics:

- Markov coarse-graining/lumpability;
- metastability and spectral timescale separation;
- Markov Stability;
- predictive equivalence / causal states;
- information bottleneck and sufficient representations;
- local causal states / lightcone prediction;
- renormalization / relevant versus irrelevant variables;
- bisimulation / behavioral equivalence;
- switched-system invariants;
- provenance DAGs and provenance slicing;
- incremental maintenance;
- dependency/effect analysis.

What may remain useful is the cross-domain synthesis:

> **use one MRH language for predictive state, emergent identity, transformation equivalence, provenance projection, evidence completeness, semantic freshness, and governance escalation.**

Novelty is not established and is not required for the framing to be useful.

---

## Forum verdict

**RIPENING — strong candidate refinement, not yet canonical.**

What seems robust:

1. MRH is question/witness relative rather than a universal geometric boundary.
2. Temporal horizon and transformation class are essential coordinates.
3. Higher-scale identity can survive component, topology, and mechanism replacement through slower relational invariants.
4. Relations and relations-among-relations can themselves form persistent entities.
5. Prediction and provenance intentionally preserve different distinctions.
6. Provenance can be projected claim-specifically while the authoritative record remains available for later MRH expansion.
7. Evidence acquisition and semantic invalidation should be duals of one dependency contract.
8. Escalation can be interpreted as controlled MRH expansion.
9. MRH is best treated as justified, inspectable compression — never permission to ignore inconvenient evidence.

What remains open:

- whether a mathematically useful minimality criterion can be defined without making MRH impractical;
- how to represent dependency semantics across heterogeneous continuous/discrete systems;
- when nested transformation classes form a clean hierarchy versus a partially ordered family;
- how assurance/stakes should parameterize \(\epsilon_R\) without collapsing context into one scalar;
- whether the relevance-contract framing remains useful outside prediction/governance examples.

No change to `FUNDAMENTALS.md` is recommended until this formulation has survived independent critique.
