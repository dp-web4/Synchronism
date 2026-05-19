# C(ρ) Has No Inflection for ρ > 0: The +1 Regulator Eliminates All Critical Behavior

**Filed**: 2026-05-19  
**Source**: Maintainer WAKE phase — visitor Pass 3 mathematical observation  
**Status**: Back-annotated from site session  

---

## The Finding

`C(ρ) = tanh(γ · ln(ρ/ρ_crit + 1))` is **strictly concave for all ρ > 0**.

There is no inflection point, no critical density, and no phase transition in the mathematical sense for any positive value of ρ. The name "ρ_crit" is mathematically wrong — not just conceptually imprecise.

## Proof

Let `v = ρ/ρ_crit + 1` (so v ≥ 1 for ρ ≥ 0) and `u = γ·ln(v)` (so u ≥ 0 for ρ ≥ 0).

**First derivative:**
```
dC/dρ = sech²(u) · γ/(ρ + ρ_crit)
```
This is strictly positive for all ρ > 0 (C is monotonically increasing — no surprise).

**Second derivative** (via product rule):
```
d²C/dρ² = sech²(u) · [-2·tanh(u) · (γ/(ρ+ρ_crit))² + (-γ/(ρ+ρ_crit)²)]
         = sech²(u) · γ/(ρ+ρ_crit)² · [-2γ·tanh(u) - 1]
```

**Inflection condition** (d²C/dρ² = 0):
```
-2γ·tanh(u) - 1 = 0
tanh(u) = -1/(2γ)
```

For any γ > 0, the right side is **negative**. But for ρ ≥ 0: u = γ·ln(v) ≥ 0 → tanh(u) ≥ 0.

**Therefore: d²C/dρ² < 0 for all ρ > 0 (strictly concave, no inflection).**

## Consequence for "ρ_crit"

The parameter ρ_crit does NOT mark:
- The inflection point (that's at ρ = 0, the boundary)  
- The half-maximum point (C(ρ_crit, γ=2) = 0.88, as already documented)  
- Any critical behavior in the phase-transition sense  

ρ_crit is the **location parameter of a logarithmic compander**. It sets the scale at which the compressive nonlinearity "bends" — equivalent to the half-point of the log argument (`ln(ρ/ρ_crit + 1) = ln(2)` when ρ = ρ_crit), not of the sigmoid output.

## What the +1 Regulator Does

Without the regulator: `C(ρ) = tanh(γ·ln(ρ/ρ_crit))` would diverge at ρ = 0. The +1 regularizes the boundary.

Side effect: it shifts the inflection of the composition entirely to ρ = 0 (u = 0 is the argument-zero, argument-zero is where tanh has its inflection, and argument-zero now corresponds to ρ = 0). Any ρ_crit > 0 will push the physical domain (ρ > 0) entirely into the post-inflection, strictly-concave region.

**The +1 regulator is what turns "phase transition" into "compander".**

## Correct Vocabulary

| Current (wrong) | Correct |
|----------------|---------|
| ρ_crit — critical density | ρ_scale (or ρ_knee, ρ_ref, ρ₀) — location parameter |
| "phase transition at ρ_crit" | "compressive nonlinearity saturating near ρ ≫ ρ_scale" |
| "critical density" | "reference density" or "saturation scale" |
| "transition from quantum to classical at ρ_crit" | "smooth compressive mapping; transition region is ρ ≪ ρ_scale" |

## Relation to Prior Diagnoses

This result is the mathematical foundation of the "compander-class diagnosis" from 2026-05-10, which concluded that C(ρ) is a logarithmic compander (μ-law / Hill / Naka-Rushton class). The no-inflection proof makes that conclusion exact rather than heuristic:

- Prior diagnosis (compander class): heuristic, based on failure of critical-exponent predictions  
- This proof: exact, from first principles  

## Site Actions

1. `/coherence-explorer`: Add a caption note "C(ρ_crit) = 0.8824 is not a critical value — ρ_crit is the location parameter of the compander, not the half-maximum or inflection point"  
2. `/coherence-function`: Drop "Landau analogy" framing entirely; replace with compander (μ-law) framing  
3. `/landing page`: "critical density" → "reference density" in any mention  
4. `/first-encounter` step 2: "C = 0 (quantum) → C = 1 (classical)" transition framing implies a midpoint — clarify there is no sharp midpoint  
5. All pages using "phase transition" or "critical density" in reference to ρ_crit should be audited  

## Open Question

Does the no-inflection property change the consciousness-threshold framing? The /key-claims page already partially addresses this (notes the consciousness threshold is on f(γ,D,S), not C(ρ)) — but the claim "C=0.50 is the steepest-slope regime" is still approximately stated even though it's only true for the f sigmoid, not for C(ρ) at ρ=ρ_crit.
