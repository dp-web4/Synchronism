# Markov Coherence Arc — Phase 2: Timescale-Dependent Entity Discovery

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 2 toy-model result  
**Code:** [`simulations/markov_phase2_partition_discovery.py`](../simulations/markov_phase2_partition_discovery.py)  
**Scope:** formalism probe only; no physics claim

## Question

Phase 1 assumed the analyst already knew the candidate entity boundary. Phase 2 removes that privilege:

> Given a hierarchical Markov process and many possible partitions, can a dynamical criterion recover useful macro-entities without being told where they are?

A second question emerged during construction:

> Is there one privileged entity decomposition, or does the appropriate entity scale depend on the observation timescale?

The second question turned out to be the more interesting one.

## Construction

An 8-state Markov chain was built with two nested dynamical scales:

- four strongly coupled 2-state pairs: `(0,1),(2,3),(4,5),(6,7)`;
- nested inside two more weakly coupled 4-state groups: `(0,1,2,3),(4,5,6,7)`;
- with weak leakage between the two 4-state groups.

A small fixed perturbation was added to every transition probability so the intended partitions were not exactly symmetric or perfectly lumpable.

All set partitions of the 8 states into 2–4 blocks, excluding singleton blocks, were enumerated: **714 candidate partitions**.

For each partition we measured dynamic coherence, approximate strong-lumpability error, and **Markov Stability** across observation times.

## Result 1 — the hierarchy is recovered without supplying the partition

The highest-Stability partition among all 714 candidates was:

| observation time | discovered best partition |
|---:|---|
| 1 | `(01)(23)(45)(67)` |
| 2 | `(0123)(4567)` |
| 3 | `(0123)(4567)` |
| 5 | `(0123)(4567)` |
| 10 | `(0123)(4567)` |
| 20 | `(0123)(4567)` |
| 50 | `(0123)(4567)` |

Thus the algorithm recovers both intentionally embedded levels — but **at different dynamical timescales**.

At one step, the 2-state pairs are the most meaningful persistent communities. Once the observation window extends beyond the pair-scale mixing/leakage time, the two larger 4-state structures become the more stable description.

This is known behavior of Markov Stability in multiscale systems; the Synchronism-relevant consequence is the interpretation.

## Result 2 — the larger-scale entities are substantially more dynamically coherent

For the four 2-state blocks, the geometric-mean coherence is

\[
\mathcal C_{pairs}\approx3.72.
\]

For the two 4-state blocks,

\[
\mathcal C_{macro}\approx16.92.
\]

The larger structures retain identity for many more internal relaxation cycles than the smaller pairs.

That does **not** make the pair level unreal. It means the pair level is the useful entity description only over a shorter temporal MRH.

## Result 3 — approximate closure is good at both levels

Approximate strong-lumpability errors were

\[
\epsilon_{lump,pairs}\approx0.0124,
\qquad
\epsilon_{lump,macro}\approx0.0169.
\]

Both partitions therefore provide reasonably autonomous coarse descriptions despite the perturbation.

This reinforces the Phase-1 separation:

- coherence is a time-scale-separation property;
- predictive closure is a coarse-graining property;
- neither uniquely selects a universal entity scale.

Observation timescale matters.

# The conceptual update

The initial arc language implicitly treated an "entity" as though a single correct boundary existed and then repeated fractally.

Phase 2 suggests a better statement:

> **Entityhood is indexed by observation timescale. A subset may be a coherent predictive entity over one temporal MRH and merely an internal degree of freedom of a larger entity over another.**

This sharpens MRH: it should not be only spatial/contextual. It likely needs an explicit **forecast/observation timescale**.

A provisional notation is

\[
\operatorname{MRH}(Y;h,\epsilon),
\]

where \(Y\) is the pattern/entity, \(h\) is the prediction or observation horizon, and \(\epsilon\) is tolerated lost predictive information.

This avoids a category mistake: asking for "the MRH of an entity" without specifying *over what duration the entity must remain predictively closed*.

## The fractal hierarchy may therefore be spectral

A more precise hierarchy than simple spatial nesting is

\[
\text{fast modes}\rightarrow\text{intermediate metastable modes}\rightarrow\text{slow modes}.
\]

Each spectral gap creates a range of times over which some degrees of freedom can be forgotten while a slower collective identity remains predictable.

That gives a potentially rigorous interpretation of fractal emergence:

> **A new scale appears when fast internal modes can be averaged away while slower collective modes remain dynamically distinguishable.**

This is standard multiscale-dynamics mathematics. The interesting Synchronism contribution, if any, would be connecting it consistently to MRH, recurrence-based entityhood, and interaction classes — then seeing whether that combined language transfers usefully.

## Phase-2 verdict

**PASS as a useful sharpening; still no novelty claim.**

What changed:

1. The method recovered nested entity partitions without being told them.
2. The preferred partition depended strongly on observation time.
3. Entityhood therefore appears naturally **timescale-relative**.
4. MRH should probably be parameterized by a forecast/observation horizon.
5. A fractal hierarchy may be most naturally understood through spectral time-scale separation rather than spatial nesting alone.

What did not happen: no new mathematics, no Synchronism-specific result, no physics evidence, and no reason to canonize a particular coherence score yet.

# Next: Phase 3 changes slightly

The original Phase-3 proposal asked for a spatial MRH radius using one-step prediction. That is too easy for any local-update model: the one-step causal neighborhood is known by construction.

The more meaningful question is:

> **How does the minimum predictive boundary grow with forecast horizon?**

So Phase 3 should measure

\[
r_{MRH}(h,\epsilon)
\]

rather than just \(r_{MRH}(\epsilon)\).

If the process has finite propagation speed, this should generate a causal-cone-like growth law. If correlations decay or coherent structures screen information, the effective MRH may grow more slowly than the raw causal cone.

That distinction — **causal reach versus relevant predictive reach** — is exactly what the "R" in MRH ought to capture.
