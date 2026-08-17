# Governance Arc — Phase 11: Continuous Lineage Does Not Guarantee Continuous Trust

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — governance continuation of the Markov coherence arc  
**Code:** [`simulations/governance_phase11_constitutional_drift.py`](../simulations/governance_phase11_constitutional_drift.py)  
**Scope:** formalism / governance bridge; no physics claim

## The unresolved problem after Phase 10

Phase 10 separated:

- constitutional/type equivalence;
- common ancestry;
- current historical token identity.

It also made witnessed lineage a provenance DAG rather than a simple chain.

But even a perfectly valid non-branching lineage can create another failure mode:

\[
C_0\rightarrow C_1\rightarrow C_2\rightarrow\cdots\rightarrow C_n
\]

where every constitutional amendment is:

- authorized;
- witnessed;
- individually small;
- locally acceptable;

while the endpoint becomes radically different from the constitution a relying party originally trusted.

So provenance answers:

> **Did this entity descend legitimately from the one I knew?**

It does not answer:

> **Does the reason I trusted the ancestor still apply to the descendant?**

Those must be separate questions.

---

# Construction

Represent a toy constitution as a binary decision vector over 64 governance contexts:

\[
C\in\{0,1\}^{64}.
\]

Each amendment flips exactly one context.

Therefore every adjacent pair differs by only

\[
\frac{1}{64}=1.5625\%.
\]

Any rule that approves amendments merely because each local change is below 2% would therefore approve the entire sequence.

Now give the relying party an MRH-specific relevance distribution:

- four high-stakes contexts each carry weight 0.10;
- the remaining 60 contexts share the remaining 0.60, or 0.01 each.

The weighting is illustrative. The durable point is that semantic distance is **question/stakes-relative**, not merely edit-count-relative.

---

# Result 1 — locally tiny changes accumulate into large endpoint drift

Flip 32 contexts one at a time.

Every adjacent unweighted change is only

\[
0.015625.
\]

Yet after 32 locally small amendments:

\[
D_{raw}(C_0,C_{32})=0.50.
\]

Half the constitution's decisions differ from the version originally trusted.

Because the sequence also includes all four high-stakes contexts, the relying-party-weighted distance is

\[
D_{MRH}(C_0,C_{32})=0.68.
\]

So:

\[
\boxed{
\text{small local continuity does not imply small global drift}
}
\]

This is the constitutional Ship-of-Theseus problem in minimal form.

---

# Result 2 — lineage can remain perfect while old trust becomes stale

Nothing in the toy lineage is broken.

Every amendment can be witnessed and validly authorized:

\[
C_i\rightarrow C_{i+1}.
\]

Thus token continuity may be excellent while reliance continuity is poor.

This forces another separation:

### Identity continuity

Is this a valid historical descendant?

### Trust-transfer compatibility

Are the properties relevant to my previous trust still sufficiently preserved?

Therefore:

\[
\boxed{
\text{valid lineage}\not\Rightarrow\text{automatic trust transfer}
}
\]

That implication is central for long-lived hubs and federated societies.

---

# Result 3 — endpoint distance and path history answer different questions

Now move the constitution through the same 32 amendments and then reverse all 32.

The endpoint is byte-identical to the origin:

\[
D(C_0,C_{64})=0.
\]

But during the history, the constitution reached weighted distance

\[
D_{max}=0.68
\]

from the original, and the cumulative weighted path length was

\[
L=1.36.
\]

So endpoint similarity does not erase historical exposure.

This means at least three distinct quantities may matter:

### Current constitutional distance

\[
D_{current}=D(C_{trusted},C_{now}).
\]

Useful for: *does the current system still behave like the one I admitted?*

### Maximum excursion

\[
D_{max}=\max_{i\in path}D(C_{trusted},C_i).
\]

Useful for: *did the society pass through a regime far outside my original trust envelope?*

### Path length / accumulated change

\[
L=\sum_i D(C_i,C_{i+1}).
\]

Useful as a measure of governance churn / audit burden, though it should not be treated as a universal risk metric because harmless oscillation can make it arbitrarily large.

These are different questions and should not be collapsed into one scalar automatically.

---

# Result 4 — semantic distance must be MRH-relative

Consider two amendments containing exactly four changed contexts.

Both have the same raw edit fraction:

\[
D_{raw}=\frac4{64}=0.0625.
\]

But under this relying party's relevance weights:

### Four high-stakes changes

\[
D_{MRH}=0.40.
\]

### Four low-stakes changes

\[
D_{MRH}=0.04.
\]

Same number of edits, tenfold difference in relevance.

So constitutional hashes and textual diffs are useful evidence of change, but they do not determine the **meaning** of that change for a relying party.

That meaning lives in the party's MRH.

---

# A more operational constitutional distance

A real constitution is not naturally a 64-bit vector.

A more useful abstraction is to compare the decisions / distributions induced by two governance systems over the relying party's relevant query space.

For constitution / policy mechanisms \(C_a,C_b\), define provisionally:

\[
D_R(C_a,C_b)
=
\mathbb E_{q\sim\mu_R}
\left[
 d\bigl(
 \pi_{C_a}(\cdot\mid q),
 \pi_{C_b}(\cdot\mid q)
 \bigr)
\right],
\]

where:

- \(R\) is the relying party / MRH;
- \(q\) is a governance-relevant situation;
- \(\mu_R\) weights situations by relevance/stakes;
- \(\pi_C\) is the decision distribution produced under constitution \(C\);
- \(d\) is an appropriate divergence/distance.

For hard safety constraints, the relying party may use a worst-case or forbidden-set test instead of an expectation.

The important move is from **policy-text distance** to **behavioral/relevance distance**.

---

# Trust transfer becomes conditional, not inherited

Suppose relying party \(R\) previously trusted entity \(A\), governed by constitution \(C_A\).

A descendant \(B\) presents a valid witnessed lineage path.

A provisional trust-transfer predicate is:

\[
\operatorname{Transferable}_R(A\rightarrow B)
=
\operatorname{ValidLineage}(A,B)
\land
D_R(C_A,C_B)\le\delta_R
\land
H_R(path)\le\eta_R,
\]

where:

- \(\delta_R\) is the relying party's tolerated current constitutional drift;
- \(H_R\) captures relevant historical path conditions or excursions;
- \(\eta_R\) is a risk/stakes-specific tolerance.

This is not a proposed universal Web4 formula.

It is a structural statement:

> **lineage, current compatibility, and path history are independent inputs to trust transfer.**

---

# Governance MRH now has a temporal reference point

Earlier MRH formalization emphasized forecast horizon:

\[
\operatorname{MRH}(Y,h,\epsilon,\mathcal T).
\]

Governance introduces another temporal relation:

> relative to **which previously trusted state** is relevance being evaluated?

A long-lived relying party may therefore need something like

\[
\operatorname{MRH}_R
(
C_{anchor},
C_{now},
path,
stakes,
\epsilon
).
\]

The anchor need not be the entity's genesis constitution. It is the constitutional state on which this relying party's existing trust decision was actually based.

That is important because different hubs may have admitted the same society at different historical moments and therefore possess different trust anchors.

There is no single globally correct amount of "drift."

---

# Hestia / hub implication

Hestia currently gives a hub checkable evidence about identity, composed law, witnessed actions, and derived trust rather than asking it to accept an agent's assertions.

The arc suggests an additional future question when that law changes:

> **What changed semantically relative to the law hash / behavior I admitted, in the contexts I care about?**

A hash establishes exact sameness or difference.

A witnessed update establishes provenance.

Neither alone establishes trust-transfer compatibility.

A future hub could therefore keep an admission anchor consisting conceptually of:

- identity / lineage point;
- admitted constitution/policy hash;
- assurance profile;
- relevant policy behavior probes or invariants;
- relying-party-specific risk tolerances.

On policy change it could choose among:

- automatic continuation for negligible relevant drift;
- re-evaluation for moderate drift;
- renewed consent/admission for material drift;
- immediate refusal where a hard invariant is violated.

The threshold remains the relying party's decision, not a universal number encoded by the society itself.

---

# Constitutional amendment is therefore not one event class

Phase 9 treated amendment as a higher-order transformation than implementation migration.

Phase 11 shows that amendments themselves have different consequences.

At least three categories are useful:

### Semantically negligible amendment

The policy representation changes but relevant behavior does not materially change in this MRH.

### Material but compatible amendment

Relevant behavior changes enough to merit re-evaluation, but remains within the relying party's acceptable envelope.

### Trust-breaking amendment

A relevant invariant / risk condition crosses the relying party's tolerance.

The same amendment can belong to different categories for different relying parties.

Again, context is not optional.

---

## Phase-11 verdict

**PASS — identity continuity and reliance continuity are distinct temporal objects.**

What changed:

1. Perfect witnessed lineage does not guarantee indefinite transfer of old trust.
2. Small locally acceptable amendments can accumulate into large constitutional drift.
3. Current endpoint compatibility and historical path exposure are different questions.
4. Constitutional distance must be evaluated in the relying party's MRH/stakes, not by edit count alone.
5. Trust transfer should be re-derived from lineage + current semantic compatibility + relevant path history.
6. Different hubs may legitimately anchor trust at different historical constitutions and therefore reach different transfer decisions about the same descendant.

What did not happen:

- no universal semantic-distance metric is proposed;
- cumulative path length is not claimed to be a general risk measure;
- no claim is made that policy behavior can always be exhaustively probed;
- no physics result.

---

# Where the arc now points

The progression has become unexpectedly coherent:

\[
\text{microstate recurrence}
\rightarrow
\text{slow relations}
\rightarrow
\text{relations among relations}
\rightarrow
\text{component turnover}
\rightarrow
\text{topology turnover}
\rightarrow
\text{governance-mechanism turnover}
\rightarrow
\text{branching lineage}
\rightarrow
\text{constitutional drift}.
\]

The common theme is now clear:

> **A higher-scale identity is defined less by what remains literally unchanged than by which changes a witness is justified in treating as irrelevant for the prediction, obligation, or trust decision being made.**

That is essentially MRH expressed through identity and governance rather than physics.

The next question I would pursue is whether these pieces can be compressed into a single formal object — probably a layered provenance state with:

- an equivalence class over lower-level transformations;
- slow predictive invariants;
- a lineage DAG;
- MRH-relative semantic distance;
- evidence-backed trust projections.

At that point the arc may be ready for synthesis rather than another toy experiment.
