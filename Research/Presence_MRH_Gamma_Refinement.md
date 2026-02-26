# Presence / MRH / γ Refinement

**Date**: 2026-02-26
**Source**: Nova forum contributions (synchronism-site/forum/nova/)
**Status**: Adopted for site; research integration in progress

---

## Summary

Three coherent refinements tighten the interpretation of the Synchronism equation's parameters without changing the mathematics:

```
C(ρ) = tanh(γ · log(ρ/ρ_crit + 1))
```

### 1. ρ: From Density to Presence

**Previous**: ρ = density (g/cm³, kg/m³)
**Refined**: ρ = presence — a scalar representation of compatible structural elements available within a given MRH, sufficient to support emergent coherence.

Presence is the scalar projection of a multidimensional compatibility vector. Density is one form of presence, but presence also encompasses temperature, energy levels, catalytic surfaces, field gradients, and lower-fractal scaffolding.

Formally: ρ = f(compatibility vector)

Operational constraints:
1. Quantifiability — must map to measurable observables
2. Domain transparency — transformation from observables → presence must be explicit
3. MRH dependence — presence is defined relative to a specific MRH
4. Falsifiability — incorrect presence definition must produce failed predictions

### 2. MRH: Tighter Operational Definition

**Previous**: "The boundary beyond which interactions become statistically irrelevant"
**Refined**: "The minimal set of interacting degrees of freedom whose state transitions materially influence the coherence evolution of a defined system"

Key properties: Minimal, Markovian, Relevancy-bound, Scale-relative.

Operational criteria:
- **Predictive sufficiency**: Removing any element inside MRH degrades coherence prediction
- **Predictive closure**: Adding elements outside MRH does not materially improve prediction

If inclusion of additional DOF changes predicted coherence behavior, the MRH was incorrectly specified.

### 3. γ: MRH Coupling Density

**Previous**: γ = 2/√N_corr (coupling strength)
**Refined**: γ also encodes MRH structural interaction efficiency

γ ∝ λ · K_MRH / D_MRH

Where:
- λ = interaction strength scaling factor
- K_MRH = connectivity (interaction density between elements)
- D_MRH = dimensionality (effective degrees of freedom)

This gives the N_corr formula a structural interpretation: it's not just "how many particles correlate" but "how efficiently the interaction network converts presence into coherence."

---

## Compatibility with Existing Work

### MRH_COMPLEXITY_FORMALIZATION.md

The existing formalization defines MRH as H = (ΔR, ΔT, ΔC) — spatial, temporal, and complexity extent. Nova's "minimal set of interacting DOF" maps directly to the ΔC (complexity) dimension, while adding operational criteria (predictive sufficiency/closure) that constrain how the boundary is drawn. These are complementary, not contradictory.

### GAMMA_UNIFICATION.md

γ = 2/√N_corr gives the formula. Nova's γ ∝ λ·K/D gives the structural reason WHY the formula works. The factor of 2 comes from phase-space contraction; N_corr captures the connectivity-to-dimensionality ratio of the correlated ensemble.

### Coupling-Coherence Experiment

The experiment already partially validated the presence reframe by treating ρ as coupling probability (rate of compatible information sharing) rather than physical density. The Hill function fit (n ≈ 3.1) suggests cooperative binding dynamics — consistent with presence as compatibility-dependent rather than quantity-dependent.

### Compatibility Lens Insight

The earlier observation that "ρ should be density of COMPATIBLE elements" (Research/Compatibility_Lens_Insight.md) is now formalized by Nova's presence definition. The compatibility lens was the intuition; presence is the formalization.

---

## Open Questions

### Can N_corr be expressed as K_MRH / D_MRH?

If N_corr ∝ K/D, then:
- γ = 2/√N_corr = 2/√(K/D) = 2·√(D/K)
- This would connect both γ formulas into a single expression
- Testable: compute K and D for systems where N_corr is known, check if K/D ∝ N_corr

### Does presence-as-compatibility predict p_crit better than density alone?

The coupling-coherence experiment's derived p_crit attempt used information-theoretic arguments. With the presence reframe, p_crit should be derivable from the compatibility structure of the MRH, not just the raw information content.

### Fractal implications

If γ is MRH-dependent and MRH is scale-relative, then γ at one fractal level is structurally constrained by the MRH at that level. This means the γ = 2 classical limit (N_corr = 1, single uncorrelated entities) corresponds to the maximally simple MRH — one entity, full coupling, no dimensional dilution.

---

## Source Documents

- `synchronism-site/forum/nova/Refining-From-Density-to-Presence.md`
- `synchronism-site/forum/nova/Refining-Markov-Relevancy-Horizon.md`
- `synchronism-site/forum/nova/linking-MRH-to-gamma.md`
