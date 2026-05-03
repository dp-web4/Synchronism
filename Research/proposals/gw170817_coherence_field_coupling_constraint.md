# Proposal: GW170817 Constraint — Does the Coherence Field Modify Gravitational Wave Propagation?

**Filed:** 2026-05-03  
**Source:** Site maintainer session — Pass 4 researcher feedback surfaced a specific observational gap  
**Priority:** HIGH — potential immediate falsification or scope clarification

---

## The Constraint

GW170817 (2017 binary neutron star merger with optical counterpart AT2017gfo) produced the bound:

```
|c_GW − c| / c < 10⁻¹⁵
```

This single measurement ruled out essentially all "dark sector with derivative coupling" modifications of gravity at galactic scales. Specifically killed: TeVeS, most Horndeski-class theories, Galileon models, and a broad swath of scalar-tensor parameter space where the scalar field couples derivatively to gravity.

---

## The Research Question

Does the Synchronism coherence field C(ρ) modify gravitational wave propagation speed?

This question has three possible answers, each with different implications:

### Case 1: Yes, the field couples derivatively to the metric

If the coherence field ∂_μC or C itself multiplies the kinetic term for gravitons, then the effective GW speed would differ from c. The modification would be:

```
c_GW² = c² × (1 + α_C × C(ρ))
```

for some coupling coefficient α_C. GW170817 would then require |α_C × C(ρ)| < 10⁻¹⁵ at galactic densities. Since C(ρ) ~ 0.01–0.99 at galactic densities, this effectively kills α_C entirely — the coupling is zero to 15 decimal places. That means the coherence field cannot affect GW propagation, which raises the question of what dynamical channel it uses to modify galactic rotation curves.

**Implication:** If Case 1, Synchronism is likely already falsified by GW170817 in the same way TeVeS was.

### Case 2: The field couples non-derivatively (potential coupling only)

If C(ρ) enters only via potential terms (not derivatives of the metric), graviton kinetic terms are unchanged and c_GW = c exactly. This is consistent with GW170817.

**But:** This raises a new question. If the coherence field doesn't couple to the gravitational kinetic term, how does it modify gravitational dynamics at galactic scales? Standard GR dynamics arise from the Einstein-Hilbert kinetic term. A potential-only coupling would contribute like a position-dependent cosmological constant, which is not what the galaxy rotation ansatz appears to be doing.

**Implication:** Survives GW170817 but needs an explicit dynamical mechanism. The current framework gives no field equations for C(ρ) — it's specified, not derived. This is the same kinematic-layer problem already identified (Sessions: lorentz_invariance_gap_kinematic_layer.md).

### Case 3: The framework is not a field theory — it's a parameterization

If Synchronism makes no claim that C(ρ) is a dynamical field (no Lagrangian, no action, no propagating modes), then GW170817 is simply not applicable. C(ρ) would be a phenomenological ansatz about density-dependent coherence, with no prediction about GW propagation.

**Implication:** Survives GW170817 trivially, but only because it makes no claims at that level. This is the honest current state (no action principle, no equation of motion). The cost is: GR violations at galactic scales need a different mechanism — possibly a non-geometric, information-theoretic one that doesn't propagate at all.

---

## Why This Matters

1. **Immediate credibility issue:** Researchers in modified gravity will ask about GW170817 within minutes of reading the site. Currently the site says nothing. The silence reads as ignorance, not as a deliberate scope restriction.

2. **Positive framing opportunity:** Case 3 is not a failure — it's a precise scope statement. If Synchronism is a parameterization (not a field theory), it cannot be falsified by GW170817, and that's a feature, not a bug. But the framework needs to say this explicitly.

3. **Connection to the kinematic-layer gap:** The GW170817 question is the same question as "what is the equation of motion for C(ρ)?" If there is no EOM, there is no propagating mode, and GW constraints don't apply. If there IS an EOM (a necessary step toward a real theory), it will have to be a Lorentz-covariant one — which connects directly to `lorentz_invariance_gap_kinematic_layer.md`.

---

## Proposed Research Tasks

1. **Classify the current framework status** on the field/parameterization axis. Is C(ρ) a dynamical field or a static ansatz? The answer determines which observational constraints apply.

2. **If field**: derive the kinetic term for C(ρ) from a Lagrangian consistent with Lorentz invariance and check whether GW speed is modified.

3. **If parameterization**: write a scope statement clarifying that GW constraints (TeVeS, Horndeski) apply to theories with action principles, and Synchronism currently lacks one. This is honest and distinguishes the framework from the alternatives it might otherwise be confused with.

4. **Write a short page `/gw170817-constraint`** that walks through the reasoning above — even if the conclusion is "this constraint doesn't apply to a parameterization" — because researchers WILL ask and deserves an answer.

---

## Related Gaps

- `lorentz_invariance_gap_kinematic_layer.md` — GW speed is one face of the Lorentz-invariance question
- The Born rule / N_corr / dual-C trilogy — three faces of the missing kinematic layer
- GW170817 is a fourth concrete constraint that the kinematic layer, once specified, must survive

---

## Precedent in the Literature

The DHOST (Degenerate Higher-Order Scalar-Tensor) theories were specifically engineered after GW170817 to preserve c_GW = c while still modifying gravity. They do this by choosing Lagrangian coefficients that decouple the tensor kinetic term from the scalar field. This class of theories survived GW170817 while TeVeS didn't.

If Synchronism ever develops an action principle, DHOST-type decoupling is the template for surviving GW170817. But this requires committing to a specific Lagrangian structure — which is exactly what the kinematic-layer proposals have been advocating.
