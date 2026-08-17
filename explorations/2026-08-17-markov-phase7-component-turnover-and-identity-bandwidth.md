# Markov Coherence Arc — Phase 7: Component Turnover and Identity Bandwidth

**Date:** 2026-08-17  
**Status:** `[ACTIVE-MRH]` — Phase 7 toy-model result  
**Code:** [`simulations/markov_phase7_component_turnover.py`](../simulations/markov_phase7_component_turnover.py)  
**Scope:** formalism probe only; no physics claim

## Question

Phases 5–6 showed that a slow relational variable can survive transformations that erase lower-scale implementation details, even when the learner is not given the relational vocabulary in advance.

The stronger question is the one motivating this phase:

> **Can a learned larger-scale identity remain continuous while every original component is physically replaced?**

And if yes:

> **Is there a measurable replacement rate beyond which that identity can no longer remain coherent?**

This turns component replacement from a philosophical example into a timescale-separation problem.

---

## Construction

Use eight role positions arranged only through a relational constraint. Each role is occupied by a replaceable component carrying a binary state

\[
x_i\in\{-1,+1\}.
\]

The role topology favors an alternating collective mode. Define fixed role signs

\[
r_i=(+1,-1,+1,-1,\ldots),
\]

so the hidden collective order parameter is

\[
M=\frac1N\sum_i r_i x_i.
\]

The learner is **not** given \(M\), the role signs, or component identities.

Two event classes occur:

1. **Relational repair/update:** one occupant updates from the current collective relation.
2. **Component turnover:** one current occupant is removed, replaced by a brand-new component identity, and initialized randomly.

Thus the implementation can change independently of the larger relational organization.

A weak fixed local-field perturbation breaks perfect symmetry, preventing the learned mode from being merely an artifact of an exactly symmetric toy.

The analysis has two arms:

### Raw-state learning arm

The learner sees only raw 8-bit transitions, estimates an empirical reversible Markov operator, and extracts its slowest nonconstant mode.

Only afterward is that learned mode compared with the hidden staggered relational order parameter.

### Explicit ancestry arm

A separate simulation tracks component IDs solely to determine when all original occupants have been replaced.

Component identity is therefore used as a diagnostic, never as an input to the learned macrostate.

---

# Three relevant timescales

This experiment naturally exposes three distinct temporal scales.

### Repair time

\[
\tau_{repair}
\]

How quickly local disruption introduced by turnover is absorbed back into the relational mode.

### Ancestry-loss time

\[
\tau_{ancestry}
\]

How long until every component present at the start has been replaced at least once.

For uniform random replacement of one of \(N\) roles at turnover probability \(q\), the mean full-turnover time is the coupon-collector result

\[
\tau_{ancestry}=\frac{N H_N}{q},
\]

where \(H_N\) is the \(N\)-th harmonic number.

For \(N=8\),

\[
N H_N\approx21.7429.
\]

### Macro-identity time

\[
\tau_{identity}
\]

The correlation time of the slow learned macro-mode,

\[
\tau_{identity}=-\frac{1}{\ln\lambda_2},
\]

where \(\lambda_2\) is the largest nontrivial eigenvalue.

The interesting regime is

\[
\boxed{
\tau_{repair}\ll\tau_{ancestry}<\tau_{identity}
}
\]

because then local replacement is repaired quickly, the entire original component population disappears, and yet the larger identity remains the slower variable.

---

# Result 1 — the macro-identity outlives its complete original membership

At the default turnover probability

\[
q=0.08,
\]

the exact transition operator gives approximately

\[
\lambda_2=0.998659,
\]

and therefore

\[
\tau_{identity}\approx744.96\text{ microsteps}.
\]

The expected time until all eight original occupants have disappeared is

\[
\tau_{ancestry}\approx271.79\text{ microsteps}.
\]

So

\[
\boxed{
R_{turn}
=\frac{\tau_{identity}}{\tau_{ancestry}}
\approx2.74
}
\]

The larger dynamical identity persists for nearly three complete component-population replacement times.

This is the operational version of the component-replacement intuition:

> **The identity is not stored in the continued presence of its original components.**

It is stored in the slower relational dynamics that repeatedly recruits current components into the same macro-pattern.

---

# Result 2 — identity persistence has a replacement bandwidth

Sweep the turnover probability \(q\):

| turnover \(q\) | \(\tau_{identity}\) | \(\tau_{ancestry}\) | \(R_{turn}=\tau_{identity}/\tau_{ancestry}\) |
|---:|---:|---:|---:|
| 0.02 | 5951.2 | 1087.1 | 5.47 |
| 0.04 | 2301.3 | 543.6 | 4.23 |
| 0.06 | 1212.0 | 362.4 | 3.34 |
| 0.08 | 745.0 | 271.8 | 2.74 |
| 0.10 | 503.3 | 217.4 | 2.31 |
| 0.12 | 362.5 | 181.2 | 2.00 |
| 0.16 | 213.9 | 135.9 | 1.57 |
| 0.20 | 141.4 | 108.7 | 1.30 |
| 0.25 | 93.3 | 87.0 | 1.07 |
| 0.30 | 66.5 | 72.5 | 0.92 |
| 0.40 | 39.1 | 54.4 | 0.72 |

The crossing occurs between roughly

\[
q=0.25\text{ and }0.30.
\]

Below that region, the macro-identity typically has a longer characteristic lifetime than the full original ancestry.

Above it, component replacement outruns the relational mode's ability to remain the slow variable.

That suggests a useful quantity:

\[
\boxed{
R_{turn}=\frac{\tau_{identity}}{\tau_{ancestry}}
}
\]

which may be read as a **turnover resilience ratio**.

- \(R_{turn}>1\): identity outlives complete original membership.
- \(R_{turn}\approx1\): ancestry turnover and identity decay occur on comparable scales.
- \(R_{turn}<1\): implementation turnover is faster than persistent macro-identity.

The precise threshold is toy-specific. The structural result is not.

---

# Identity is maintained, not possessed

This phase suggests a conceptual correction that may be more important than the numerical result.

The phrase

> "the entity has an identity"

sounds as if identity were a property stored somewhere in the entity's components.

The dynamics suggest a better formulation:

> **Identity is continuously maintained by a relation that recruits changing components into an equivalence class of implementations.**

The difference is substantial.

At time \(t_0\), implementation \(X^{(0)}\) realizes macrostate \(Y\).

Later,

\[
X^{(0)}\rightarrow X^{(1)}\rightarrow X^{(2)}\rightarrow\cdots
\]

with every original component eventually absent, while

\[
f(X^{(0)})
\approx f(X^{(1)})
\approx f(X^{(2)})
=Y.
\]

The continuity belongs to the equivalence relation \(f\), not to material ancestry.

---

# This makes the temporal/fractal MRH more explicit

Previous phases established:

1. entityhood depends on observation/forecast timescale;
2. exterior MRH depends on forecast horizon and tolerated predictive loss;
3. lower-scale distinctions can become irrelevant internally;
4. relations among relations can become the slow variables.

Component turnover now adds another requirement:

> A distinction that matters at one timescale can become merely implementation detail at a slower scale.

For the occupant-level witness, component identity is highly relevant: this particular component exists, changes, and leaves.

For the slower relational witness, *which particular component currently occupies the role* can become irrelevant so long as replacement happens within the repair bandwidth.

Thus MRH appears naturally **fractal and temporal**:

\[
\boxed{
\operatorname{MRH}
=\operatorname{MRH}(Y,h,\epsilon,\mathcal T)
}
\]

where

- \(Y\) is the target macro-pattern,
- \(h\) is the forecast/observation horizon,
- \(\epsilon\) is tolerated predictive information loss,
- \(\mathcal T\) is the class/rate of transformations being treated as internal implementation change.

The same event can therefore be:

- identity-destroying at one MRH;
- ordinary component turnover at a higher/slower MRH.

That is not inconsistency. It is scale-relative relevance.

---

# A two-sided temporal shell

Phase 5 proposed an informational shell bounded by irrelevant internal detail and irrelevant external detail.

Phase 7 makes that shell explicitly temporal.

A coherent macro-entity requires something like

\[
\boxed{
\tau_{internal}\ll h\ll\tau_{identity}
}
\]

for an observation horizon \(h\):

- processes much faster than \(h\) are averaged into implementation detail;
- processes much slower than \(h\) define the persistent identity visible at that scale.

Component turnover can itself sit on either side of the boundary depending on its rate.

So the "fractal" hierarchy need not be defined primarily by size. It can be defined by nested bands of dynamical relevance:

\[
\text{fast replacement}
\rightarrow
\text{relational repair}
\rightarrow
\text{macro persistence}
\rightarrow
\text{slower relation among macro-relations}
\rightarrow\cdots
\]

This is perhaps the cleanest formal expression yet of the earlier intuition that MRH is both scale-dependent and temporal.

---

# Ship of Theseus, but with a measurable criterion

The classic replacement puzzle asks when an object remains "the same" after its parts are replaced.

This formalism does not solve semantic identity in general.

It does replace the binary question with a dynamical one:

> **Does the candidate macro-variable remain a slower, predictively closed mode than the process replacing its implementation?**

If yes, then at that chosen MRH the replacements are internal transformations rather than identity-destroying events.

If no, the distinction can no longer be discarded.

That makes identity explicitly relative to prediction scale rather than metaphysical sameness.

---

# Phase-7 verdict

**PASS — component replacement can be cleanly separated from macro-identity loss by timescale.**

What survived:

1. A learned relational identity can outlive complete replacement of its original component population.
2. Identity persistence is a process, not a conserved material inventory.
3. Turnover has a measurable bandwidth: above a sufficient replacement rate, macro-identity no longer outruns ancestry loss.
4. The ratio \(R_{turn}=\tau_{identity}/\tau_{ancestry}\) is a useful provisional diagnostic.
5. MRH is naturally temporal as well as informational/spatial.
6. The same lower-scale change can be identity-destroying in one MRH and irrelevant implementation turnover in another.

What did not happen:

- no claim that \(R_{turn}\) is a universal entity metric;
- no new mathematics;
- no physics validation;
- no demonstration yet that the hierarchy survives open-ended role creation/deletion rather than fixed role topology;
- no proof that learned identities are unique.

---

# Next question — identity without fixed slots

This toy still grants the higher-scale entity a hidden scaffold: eight persistent **roles** remain even though their occupants turn over.

That leaves an important loophole.

Perhaps the "identity" is actually stored in the fixed role topology.

The next phase should therefore remove fixed slots as well.

Let components enter, leave, split, merge, or change adjacency so that neither component IDs nor a fixed coordinate system survives. Then ask whether a slower relational variable can still be learned.

That would test a much stronger claim:

> **Can identity persist when both the components and the coordinate system used to name them are transient?**

If so, the natural object is no longer an invariant over a fixed collection of variables. It is an invariant over a changing relational graph.

That would move the arc from relations-among-relations toward **relations among transformations of relational structure**.
