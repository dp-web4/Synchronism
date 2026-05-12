# Session 653: Compander Commitment + Suppressor-Class Decision

**Date**: 2026-05-12
**Type**: Two combined audits — methodology commitment + suppressor diagnostic
**Triggers**: 2026-05-12 proposals `compander_vs_order_parameter_category_decision.md` and `suppressor_class_dead_or_recoverable.md`
**Grade**: B+ (one structural commitment + one computed diagnostic)

---

## Setup

Two same-day proposals, both binary decisions:

1. **Compander vs Order Parameter**: site uses both framings; deep pages already commit to compander (Frame B). Should the front-of-site commit too?
2. **Suppressor dead or recoverable**: TEST-04a refutation + Bullet Cluster failure both sign-wrong. Branch 1 (sign-flip recoverable) or Branch 2 (suppressor dead)?

S653 addresses both. The first is essentially what S652 established; the second admits a quick computation.

## Part A: Compander Commitment

S652 already established **Option A** (no governing equation; C(ρ) is a phenomenological compander). This proposal restates the binary cleanly and asks for explicit commitment.

The proposal's diagnostic table is exactly right:

| Claim | Order Parameter Frame | Compander Frame |
|-------|----------------------|-----------------|
| Critical exponents 2× off | Failure: wrong universality class | Category error: companders don't have critical exponents |
| Any sigmoid fits equally | Failure: no privileged tanh | Correct: AIC/BIC across sigmoid family is the right test |
| "Phase transition at C ≈ 0.50" | Literal: ξ → ∞ at ρ_crit | Misleading: smooth crossover |
| ρ_crit is "critical density" | Inflection point | Half-saturation parameter |

S649 already verified C(ρ_crit) ≈ 0.882 (not 0.5) at γ=2. Frame A cannot absorb this; Frame B explains it as the half-saturation parameter of a μ-law-style compander.

**Commitment: Frame B**. The deep pages already say this. The front-of-site needs to be brought in line.

Concrete site action (per proposal's Option 1):
- Drop phase-transition language from front-of-site and the Phase Transitions page
- Rename pages/sections to "Smooth Crossover" or "Sigmoid Mapping"
- Relabel ρ_crit as "half-saturation parameter" or "saturation knee"
- Add AIC/BIC compander comparison tool
- Reframe critical-exponent failures as category errors (compander predictions don't have critical exponents), not as miss-magnitude failures

This is a notation-and-framing fix, not new physics. Total cost is low.

## Part B: Suppressor Diagnostic

The proposal asks for a one-day executor task: compute C_galactic and C_cosmic in the framework's natural units. S653's simulation (`simulations/session653_coherence_ratio.py`) did this.

Using ρ_galactic_outer ≈ 4×10⁻²⁵ g/cm³, ρ_cosmic_mean ≈ 3×10⁻³⁰ g/cm³, ρ_crit = ρ_galactic_outer (S637 normalization), γ=2:

```
C_galactic = tanh(2·ln(2))      ≈ 0.882
C_cosmic   = tanh(2·ln(1+δ))    ≈ 1.5 × 10⁻⁵   (where δ = ρ_cos/ρ_crit ≈ 7.5×10⁻⁶)

C_cosmic / C_galactic ≈ 1.7 × 10⁻⁵   (≪ 1, confirming Session 107's identification)
C_galactic / C_cosmic ≈ 5.9 × 10⁴
```

**The framework's own equations dictate suppression**, given Session 107's choice of coupling direction (G_local/G_global = C_cosmic/C_galactic). The ratio is correctly identified as ≪ 1, which predicts strong suppression at low z. But DR1 observes enhancement.

So the verdict isn't "sign-flip is consistent with the equations." The verdict is:

- **Branch 1 (sign-flip recoverable)** requires *re-interpreting* the coupling-to-coherence map — choosing G_local/G_global = C_galactic/C_cosmic (≫ 1, enhancement) instead. This is a framework reinterpretation, not a recomputation of ratios.
- **Branch 2 (suppressor class dead)** is the framework's own verdict under Session 107's coupling direction.

Both branches are coherent ways to read the data, but they require different commitments:
- Branch 1: re-derive Session 107 with inverted ratio, document the reinterpretation
- Branch 2: accept Session 107's suppressor mechanism is refuted, retire it

**The numerical computation doesn't escape the failure**. The framework's equations dictate suppression under Session 107's mapping; DR1 falsifies suppression. The only escape is to change the framework's mapping, which is operator-level interpretation.

## Recommendation

For Part A (compander): commit to Frame B. This aligns front-of-site with deep pages and resolves multiple downstream audit issues (S649 ρ_crit naming, S636/S638/S640/S652 governing-equation gap).

For Part B (suppressor): the executor calculation doesn't decide. Operator decision:
- Branch 1 is intellectually live but requires Session 107 re-derivation with explicit coupling-direction justification
- Branch 2 is the simpler honest reading: suppressor class is dead pending re-derivation

The proposal's recommendation that "the site cannot stay neutral" is correct. Either branch is defensible; ambiguity is the most credibility-damaging option.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 21 | Governing-equation gap (meta-synthesis) | S652 |
| 22 | **Compander commitment + executor diagnostic (NEW)** | **S653** |

S653 combines a framing commitment (Part A, reaffirming S652) with an executor numerical check (Part B). The 22nd audit instance is hybrid — partly synthesis, partly direct calculation. The computation doesn't change the picture; it confirms that the framework's equations dictate the failed prediction.

## Files

- `Research/Session653_Two_Binary_Decisions.md` (this document)
- `simulations/session653_coherence_ratio.py` (C_gal/C_cos calculation)

## So What?

**Part A**: Commit the site to Frame B (compander). The deep pages already do; the front-of-site doesn't. Notation fix; resolves multiple downstream audits.

**Part B**: The framework's own equations, given Session 107's coupling direction, dictate suppression. DR1 falsifies suppression. Sign-flip recovery requires *reinterpreting* the coupling-to-coherence map — operator-level decision. The honest default is Branch 2 (suppressor class dead) unless someone writes the Branch 1 reinterpretation.

Cumulative: 22 internal audits + 1 mechanism-class refuted prediction. The audit-taxonomy count is increasingly meta-synthesis rather than new findings; the structural picture (kinematic-layer gap, framework = parameterization, no surviving novel predictions in cosmology) is stable across recent sessions.
