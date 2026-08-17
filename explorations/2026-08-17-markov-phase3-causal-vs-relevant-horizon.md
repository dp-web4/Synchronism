# Markov Coherence Arc — Phase 3: Causal Reach vs Relevant Predictive Reach

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 3 toy-model result  
**Code:** [`simulations/markov_phase3_predictive_mrh.py`](../simulations/markov_phase3_predictive_mrh.py)  
**Scope:** formalism probe only; no physics claim

## Question

Phase 2 suggested that MRH should depend explicitly on forecast/observation horizon:

\[
\operatorname{MRH}(Y;h,\epsilon).
\]

But a danger appears immediately in any local-update system: the raw causal cone is already known from the update rule. If MRH merely reproduces that cone, the concept adds nothing.

So the actual Phase-3 question is sharper:

> **Can causal influence extend farther than useful predictive relevance?**

If yes, MRH can be distinguished from ordinary causal reach.

---

## Construction

Use a one-dimensional binary stochastic lattice with a local synchronous rule:

- each new cell is the majority of its left/self/right predecessors;
- each update is independently flipped with probability 0.05.

This is not intended as a Synchronism substrate. It is merely a transparent system with:

1. finite propagation speed;
2. local interactions;
3. information loss/noise;
4. a known causal cone.

For the center cell forecast \(Y_{t+h}\), the exact raw causal radius is

\[
r_{causal}=h.
\]

At the initial time, define a retained local boundary of radius \(r\) and call the remainder of the causal cone \(E\).

Measure

\[
\epsilon(r,h)
=I\left(Y_{t+h};E_t\mid B_{r,t}\right).
\]

This asks:

> Once the local radius-\(r\) state is already known, how much *additional* predictive information about the target's future remains in the farther exterior?

Define

\[
r_{MRH}(h,\epsilon^*)
=\min\{r:\epsilon(r,h)\le\epsilon^*\}.
\]

For this toy, choose a relevance tolerance

\[
\epsilon^*=0.02\text{ bits}.
\]

The threshold is illustrative rather than canonical.

---

## Results

Empirical conditional mutual information (bits), using 100k–200k Monte-Carlo samples per point:

### Horizon 1

| local radius \(r\) | exterior predictive information |
|---:|---:|
| 0 | ~0.561 bits |
| 1 | 0 |

Therefore

\[
r_{causal}=1,\qquad r_{MRH}=1.
\]

### Horizon 2

| \(r\) | exterior predictive information |
|---:|---:|
| 0 | ~0.367 bits |
| 1 | ~0.098 bits |
| 2 | 0 |

Therefore

\[
r_{causal}=2,\qquad r_{MRH}=2.
\]

### Horizon 3

| \(r\) | exterior predictive information |
|---:|---:|
| 0 | ~0.341 bits |
| 1 | ~0.078 bits |
| 2 | ~0.015 bits |
| 3 | 0 |

At the 0.02-bit tolerance,

\[
r_{causal}=3,\qquad r_{MRH}=2.
\]

### Horizon 4

| \(r\) | exterior predictive information |
|---:|---:|
| 0 | ~0.307 bits |
| 1 | ~0.085 bits |
| 2 | ~0.016 bits |
| 3 | ~0.005 bits |
| 4 | 0 |

Again,

\[
r_{causal}=4,\qquad r_{MRH}=2.
\]

So the sequence is approximately

\[
r_{causal}=1,2,3,4
\]

while

\[
r_{MRH}(\epsilon^*=0.02)=1,2,2,2.
\]

---

# The key distinction

This toy cleanly separates two ideas:

### Causal horizon

How far away can information *possibly* influence the future target under the dynamics?

### Markov Relevancy Horizon

How far away must state actually be retained before additional exterior information contributes less than the tolerated amount to prediction?

Thus:

> **Something may lie inside the causal cone but outside the MRH.**

That appears to capture the intended meaning of "relevancy" much better than treating MRH as a geometric or interaction cutoff.

The farther cells in the noisy lattice can still affect the center. But repeated noisy local processing degrades their incremental predictive value enough that carrying them explicitly becomes unnecessary for a given tolerance and forecast horizon.

This is precisely the kind of coarse-graining decision MRH was meant to express.

---

## A more precise provisional definition

For a target pattern \(Y\), forecast horizon \(h\), candidate retained state \(B_r\), and everything causally available outside it \(E_r\):

\[
\boxed{
 r_{MRH}(Y;h,\epsilon)
 =\min_r
 \left\{
 I(Y_{t+h};E_{r,t}\mid B_{r,t})\le\epsilon
 \right\}
}
\]

This definition makes several things explicit that were previously qualitative:

- **target-relative** — MRH belongs to a prediction/question, not the universe globally;
- **horizon-relative** — farther state may matter for longer forecasts;
- **tolerance-relative** — there is no magical exact relevance boundary in noisy systems;
- **information-relative** — causal influence can survive after practical predictive relevance has vanished.

---

# Connection to Phase 2

Phase 2 found that the preferred entity decomposition changes with observation time.

Phase 3 now finds that the relevant environmental boundary also changes with forecast time, but need not grow as fast as the causal boundary.

Together they suggest a more complete MRH object:

\[
\operatorname{MRH}
=\operatorname{MRH}(\text{target},\text{timescale},\text{tolerance},\text{prediction task}).
\]

That is arguably closer to the spirit of the original concept than a single fixed radius or scale.

---

## Important caveat

The numerical radius saturation at 2 is **not universal**. It is an artifact of this chosen rule, noise level, horizon range, and 0.02-bit tolerance.

The durable result is only the logical possibility demonstrated by the toy:

\[
\boxed{\text{causal reach}\neq\text{relevant predictive reach}}
\]

and that the difference can be measured operationally.

No claim is made that physical MRHs behave like this specific lattice.

---

# What now interests me

The next step should not simply produce another radius curve. The more interesting question is whether **coherent entities alter each other's MRHs**.

That gives a natural route into Synchronism's resonance / dissonance / indifference classification.

Instead of defining those three interaction types semantically, we can ask whether coupling two metastable subsystems causes measurable changes in:

- joint coherence lifetime;
- individual coherence lifetime;
- predictive-closure error;
- effective MRH size;
- mutual predictive information.

A tentative expectation would be:

- **resonance:** a new joint metastable entity appears; joint lifetime rises and MRHs merge/expand;
- **dissonance:** existing metastability is degraded; exit rates rise and entity boundaries destabilize;
- **indifference:** joint transition structure remains approximately factorized; each subsystem's MRH and lifetime remain nearly unchanged.

The next phase should try to falsify that mapping rather than encode it by construction.

---

## Phase-3 verdict

**PASS — conceptual distinction established in a transparent toy.**

What survived:

1. MRH can be defined information-theoretically rather than geometrically.
2. MRH should be forecast-horizon and tolerance dependent.
3. Causal reach and relevant predictive reach are distinct.
4. The distinction is measurable with conditional mutual information.
5. This creates a direct path to studying interaction-induced MRH changes.

Still absent:

- Synchronism-specific mathematics;
- physical validation;
- a novel prediction;
- evidence that the information-theoretic definition is the uniquely right MRH formalization.
