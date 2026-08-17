# Markov Coherence Arc — Phase 4: Interaction Classes Are About Relational Structure, Not Stability Alone

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 4 toy-model result  
**Code:** [`simulations/markov_phase4_coupling_classes.py`](../simulations/markov_phase4_coupling_classes.py)  
**Scope:** formalism probe only; no physics claim

## Question

Phase 3 suggested using interaction-induced changes in coherence/MRH to operationalize Synchronism's three interaction classes:

- resonance,
- dissonance,
- indifference.

The first instinct was:

- resonance → more coherent joint state,
- dissonance → less coherent / destabilized state,
- indifference → unchanged independent states.

Phase 4 shows that this is too simple.

## Construction

Use two persistent binary subsystems \(A,B\in\{-1,+1\}\).

Each subsystem has strong self-memory and a tunable coupling \(J\) to the other subsystem's current state.

- \(J>0\): favors alignment,
- \(J<0\): favors anti-alignment,
- \(J=0\): exact factorization.

The update probabilities are logistic functions of self-memory plus coupling. The model is deliberately symmetric and transparent.

Measures:

1. stationary correlation \(\langle AB\rangle\),
2. mutual information \(I(A;B)\),
3. observational transfer information \(I(A_{t+1};B_t\mid A_t)\),
4. intervention-like one-step causal sensitivity of \(A\)'s transition law to \(B\),
5. residence time of the aligned relation \(AB=+1\),
6. residence time of the anti-aligned relation \(AB=-1\).

## Result 1 — indifference is cleanly identifiable

At \(J=0\):

- stationary correlation \(\approx 0\),
- mutual information \(\approx 0\),
- transfer information \(\approx 0\),
- causal sensitivity \(=0\).

So in this toy, **indifference corresponds naturally to approximate factorization**:

\[
P(A',B'\mid A,B)\approx P(A'\mid A)P(B'\mid B),
\]

plus negligible cross-predictive information.

That is stronger and more operational than merely saying "the systems do not interact much."

## Result 2 — positive and negative coupling both create strong joint coherence

For strong positive \(J\), the stationary distribution concentrates on aligned states:

\[
(++), (--),
\]

with \(\langle AB\rangle\to +1\).

For strong negative \(J\), it concentrates just as strongly on anti-aligned states:

\[
(+-), (-+),
\]

with \(\langle AB\rangle\to -1\).

In both cases:

- mutual information approaches 1 bit,
- cross-predictive dependence becomes strong,
- the favored relational state becomes long-lived.

Thus:

> **Anti-alignment can be every bit as coherent as alignment.**

So coherence magnitude alone cannot distinguish resonance from dissonance.

This is the main correction from Phase 4.

## Result 3 — stability and interaction valence are orthogonal

The toy exposes at least two separate axes:

### Axis A — coupling / relational dependence

How strongly does each subsystem's next-state distribution depend on the other?

Possible measures:

\[
I(A_{t+1};B_t\mid A_t)
\]

or an intervention-based transition sensitivity.

### Axis B — relational geometry / sign

What structured relationship is stabilized?

For binary variables, the simplest order parameter is

\[
R=\langle AB\rangle.
\]

- \(R>0\): aligned relation,
- \(R<0\): anti-aligned relation,
- \(R\approx0\): no persistent sign relation.

A third axis is then needed for **joint persistence** — how long the relationship itself survives.

So an interaction is not adequately described by one scalar "coherence."

A more faithful object is something like

\[
\mathcal I(A,B)
=
(\text{coupling strength},\text{relation structure},\text{joint lifetime}).
\]

## Consequence for Synchronism terminology

The naive mapping

\[
\text{resonance}=\text{coherence},\qquad
\text{dissonance}=\text{decoherence}
\]

should be rejected.

A dissonant interaction can generate an extremely stable **anti-coherent / oppositional relation** while remaining highly organized.

That actually fits many ordinary physical examples better: destructive interference is not random noise; it is a very precise phase relationship.

So a better provisional translation may be:

### Indifference

Cross-dynamics approximately factorize:

\[
I(A_{t+1};B_t\mid A_t)\approx0,
\qquad
I(B_{t+1};A_t\mid B_t)\approx0.
\]

The systems may coexist but neither substantially enters the other's predictive state.

### Resonance

Coupling produces a persistent **constructive relational mode** — a joint macrostate whose relationship reinforces a shared pattern or increases a relevant joint lifetime.

### Dissonance

Coupling produces a persistent **oppositional relational mode** or systematically suppresses a target mode. The relation may itself be highly coherent.

This keeps "dissonance" from being confused with randomness.

## The deeper update: coherence may live on relations, not only entities

Phases 1–3 treated coherence primarily as a property of a state-set / entity.

Phase 4 suggests another layer:

> **Relations themselves can be metastable objects.**

The persistent thing need not be \(A\) or \(B\) separately. It may be

\[
R(A,B)=\text{aligned},
\]

or

\[
R(A,B)=\text{anti-aligned}.
\]

That means an emergent entity can be a **stable relation among entities**, not merely a larger spatial aggregate.

This is potentially important for Synchronism's fractal picture. At a higher MRH, what becomes the new "state variable" may be the relationship class itself.

Examples in generic terms:

\[
(A,B)\longrightarrow R_{AB}
\]

and then

\[
(R_{AB},R_{CD},\ldots)\longrightarrow R_{higher}.
\]

The hierarchy is therefore not necessarily "small objects combine into bigger objects." It can be:

\[
\text{states}\rightarrow\text{persistent relations}\rightarrow\text{relations among relations}.
\]

That is a stronger and more general notion of emergence.

## Where MRH enters

If \(B_t\) materially improves prediction of \(A_{t+1}\), then \(B\) lies inside \(A\)'s predictive MRH for that horizon and tolerance.

As coupling falls to zero, this contribution disappears and the systems leave each other's effective MRHs.

Thus interaction can be described as **MRH deformation**:

\[
J\uparrow
\quad\Rightarrow\quad
\text{cross-predictive relevance rises}
\quad\Rightarrow\quad
\text{MRHs merge or overlap more strongly}.
\]

But the *sign/structure* of the resulting relation is a separate question from whether the MRHs overlap.

That gives a clean separation:

- **MRH overlap:** are the systems dynamically relevant to each other?
- **relation structure:** what mode does that coupling stabilize?
- **relation lifetime:** how coherent is that mode?

## Important caveat

This binary spin-like model makes "alignment" and "anti-alignment" obvious by construction. It therefore cannot establish that those exact order parameters generalize.

The durable lesson is structural:

\[
\boxed{\text{interaction strength}\neq\text{interaction valence}\neq\text{interaction coherence}}
\]

Any future mapping of resonance/dissonance/indifference should preserve those distinctions.

## Phase-4 verdict

**PASS — and it falsified the first naive mapping.**

What survived:

1. Indifference maps naturally to approximate dynamical factorization / negligible cross-predictive information.
2. Resonance and dissonance cannot be distinguished by coherence magnitude alone.
3. Oppositional relations can be highly coherent.
4. Relations themselves can be metastable/emergent state variables.
5. Interaction can be described partly as MRH overlap/deformation, while relation structure remains a separate axis.

What failed:

- "resonance increases coherence, dissonance decreases coherence" as a general rule.

## Next question

The natural Phase 5 is now different from the original plan.

Instead of recursively coarse-graining **states only**, recursively coarse-grain **persistent relations**.

Can a hierarchy be discovered in which:

\[
\text{microstates}
\rightarrow
\text{pair relations}
\rightarrow
\text{relations among pair-relations}
\rightarrow\cdots
\]

and does the predictive closure improve at some relational level?

If so, Synchronism's "fractal entities" might admit a much more precise reading:

> **Emergent scales are stable predictive relations among lower-scale degrees of freedom, not necessarily aggregates of them.**

That is where I would go next.
