# Markov Coherence Arc — Phase 1: Metastability, Internal Dynamics, and Lumpability

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 1 toy-model result  
**Code:** [`simulations/markov_coherence_toy.py`](../simulations/markov_coherence_toy.py)  
**Scope:** formalism probe only; no physics claim

## Question

Can the intuitive Synchronism phrase **"a coherent entity is a stable recurring pattern"** be sharpened using Markov-chain mathematics without immediately collapsing into either:

1. mere persistence, where a frozen state looks maximally coherent; or
2. a renaming of standard metastability?

The Phase-1 candidate is a **time-scale separation**:

\[
\mathcal C(C)=\frac{\tau_{\rm exit}(C)}{\tau_{\rm internal}(C)}.
\]

Here:

- \(\tau_{\rm exit}\) measures how long the process remains within a candidate macro-pattern \(C\);
- \(\tau_{\rm internal}\) measures how rapidly microstates *within* \(C\) relax/rearrange.

The intuition is simple: an entity should be able to change internally many times while remaining recognizably the same macro-pattern.

---

## Toy construction

Use six microstates arranged as two 3-state clusters:

\[
A=\{0,1,2\},\qquad B=\{3,4,5\}.
\]

The internal transition matrix for each dynamic cluster is

\[
W=
\begin{pmatrix}
0.62&0.23&0.15\\
0.18&0.64&0.18\\
0.20&0.20&0.60
\end{pmatrix}.
\]

At each step, probability \(e\) leaks to the other cluster; otherwise the process follows \(W\). In the symmetric arm every microstate has the same leakage \(e\), so the partition \(A/B\) is **exactly strongly lumpable**. A perturbed arm makes the state leakage rates differ by ±25%, breaking exact lumpability while preserving the same qualitative two-cluster structure.

This construction is intentionally transparent. Phase 1 is a sanity check, not a search for novelty.

---

## Metrics

For the restricted substochastic matrix \(Q\) on cluster \(A\), let \(\rho\) be its dominant eigenvalue. We use

\[
\tau_{\rm exit}=\frac{1}{1-\rho}
\]

as a metastable residence-time proxy.

Condition \(Q\) on remaining inside \(A\), producing a row-normalized internal chain. If \(|\lambda_2|\) is its second-largest eigenvalue magnitude, use

\[
\tau_{\rm internal}=\frac{1}{1-|\lambda_2|}
\]

as an internal **relaxation-time proxy**. This is not claimed to be a rigorous total-variation mixing time; later phases can replace it with stricter measures.

Strong-lumpability error is measured as the largest difference, among states in the same source block, in aggregate transition probability to a target block.

---

## Result 1 — increasing leakage destroys time-scale separation

The internal relaxation proxy is approximately constant:

\[
\tau_{\rm internal}\approx1.756\ \text{steps}.
\]

For the exactly lumpable arm:

| leakage \(e\) | \(\tau_{exit}\) | \(\mathcal C=\tau_{exit}/\tau_{internal}\) |
|---:|---:|---:|
| 0.005 | 200.000 | 113.884 |
| 0.010 | 100.000 | 56.942 |
| 0.020 | 50.000 | 28.471 |
| 0.050 | 20.000 | 11.388 |
| 0.100 | 10.000 | 5.694 |
| 0.200 | 5.000 | 2.847 |
| 0.300 | 3.333 | 1.898 |

So the candidate measure behaves as intended in the transparent case: as the boundary becomes more permeable, internal microdynamics have fewer opportunities to cycle before macro-identity is lost.

This is **not a discovery**; it is the required sanity check.

---

## Result 2 — persistence alone fails the entity intuition

Replace \(W\) by the identity matrix while keeping leakage \(e=0.02\).

The cluster still has

\[
\tau_{exit}=50,
\]

so a persistence-only score would call it just as stable as the dynamic cluster at the same leakage.

But internally it never mixes:

\[
|\lambda_2|=1,\qquad \tau_{internal}=\infty.
\]

Therefore the proposed dynamic-coherence ratio gives

\[
\mathcal C=0.
\]

This is the most useful Phase-1 result conceptually.

> **Persistence is necessary for entity-like coherence, but not sufficient. A coherent entity must preserve macro-identity while continuing to transform internally.**

That maps more closely to Synchronism's existing entity definition — recurrence across tick sequences — than a static attractor does.

---

## Result 3 — coherence and lumpability are different axes

In the symmetric arm, the two-cluster partition is exactly lumpable for **every** tested leakage value, even while \(\mathcal C\) collapses from ~114 to ~1.9.

So:

> **A macrostate can be perfectly Markovian after coarse-graining without being strongly coherent.**

Conversely, in the perturbed arm the lumpability error rises with leakage:

| leakage \(e\) | lumpability error |
|---:|---:|
| 0.005 | 0.0025 |
| 0.010 | 0.0050 |
| 0.020 | 0.0100 |
| 0.050 | 0.0250 |
| 0.100 | 0.0500 |
| 0.200 | 0.1000 |
| 0.300 | 0.1500 |

while a metastable cluster can still exist.

This separates two notions that were easy to blur in the initial intuition:

- **coherence / entity persistence:** time-scale separation;
- **emergent predictive closure:** quality of the macrostate coarse-graining.

They may interact, but they are not the same mathematical property.

---

## What I would change in the initial proposal after Phase 1

The original phrase "entity = metastable pattern possessing an effective Markov blanket" was too compressed.

A better provisional decomposition is:

### Entity-like coherence

\[
\boxed{\text{fast internal dynamics} \ll \text{slow identity loss}}
\]

### Emergent level

\[
\boxed{\text{discarded micro-detail adds little predictive value at the macro scale}}
\]

### MRH

\[
\boxed{\text{distant exterior state adds little predictive value beyond the local boundary}}
\]

These are three separate conditional/time-scale questions. Their **intersection**, rather than any one alone, may be the interesting Synchronism object.

---

## A stronger candidate picture

Let \(Y_t\) denote the macro-identity of a candidate entity, \(X_t\) its microstate, and \(E_t\) state outside a proposed MRH boundary.

The object of interest may be a scale where all three hold simultaneously:

1. **metastability**
   \[
   \tau_{exit}\gg\tau_{internal};
   \]

2. **micro-detail closure**
   \[
   I(Y_{t+1};X_t\mid Y_t)\ll1;
   \]

3. **exterior closure**
   \[
   I(Y_{t+1};E_t\mid Y_t,B_t)\ll1.
   \]

That gives a very concrete interpretation of a fractal entity:

> A scale becomes entity-like when the process forgets its exact internal implementation and most of its distant environment **faster than it forgets its own macro-identity**.

That sentence is now the main hypothesis I would carry into Phase 2/3.

---

## Next exploration

### Phase 2 — do not hand-pick the macrostate

The current toy tells us what happens after the analyst declares \(A/B\) to be the entities. The next nontrivial step is to remove that privilege.

Construct a chain with several plausible partitions and ask whether a criterion combining:

- metastable residence time;
- internal relaxation;
- approximate lumpability / predictive closure;
- compression benefit

can **discover** the useful macrostate without being told where it is.

The key danger is circularity: if the scoring function simply rewards the structure used to construct the chain, nothing has been learned.

### Phase 3 — MRH radius

After that, move to a spatially coupled process and measure

\[
I(Y_{t+1};E_{>r,t}\mid Y_t,B_{\le r,t})
\]

as a function of radius \(r\). If there is a knee where exterior predictive contribution collapses, MRH gains an operational meaning as a **predictive-information horizon**.

That is the part of the arc I currently find most promising.

---

## Phase-1 verdict

**PASS as a formalism sanity check; no novelty claim.**

What survived:

- mere persistence is insufficient;
- time-scale separation captures the dynamic notion of recurrent identity better;
- lumpability/predictive closure is distinct from coherence;
- MRH should remain a third, exterior-dependence concept;
- the interesting object may be the coincidence of all three at a particular scale.

What did **not** happen:

- no new Markov mathematics;
- no evidence for the Synchronism substrate;
- no physics prediction;
- no reason yet to elevate the candidate score into the canonical framework.

The toy did its job: it made the next question sharper.
