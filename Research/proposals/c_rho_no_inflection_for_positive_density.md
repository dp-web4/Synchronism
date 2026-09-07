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

---

## CORRECTION (2026-09-07): the result is correct and **coordinate-dependent**, and it has been used outside its coordinate

**Source**: maintainer track, from the synchronism-site visitor researcher persona, 2026-09-07.
Back-annotated because this proposal's conclusion was being cited on the site as a general
"C(ρ) has no privileged value," which is false.

The proof above is right: in **linear ρ**, C is strictly concave for ρ > 0 and has no inflection.
But the framework never works in linear ρ. Every plot, the log-density argument itself, the
Coherence Explorer axis, and the whole "spans 80 orders of magnitude" framing are in **log ρ** —
and in that coordinate the inflection exists.

With `u = ρ/ρ_crit`, maximise `dC/d(ln ρ) = u·dC/du` for `C = tanh(γ ln(1+u))`:

```
ln f  = ln u + ln γ + ln(1 - C²) - ln(1+u)
d/du  = 1/u - 1/(1+u) - 2C·γ/(1+u) = 0
      = 1/(u(1+u)) = 2Cγ/(1+u)
```

**Inflection condition in log-density:  C* = 1 / (2 γ u*)**

- At **γ = 1/2** (the value SPARC selects at 0.489 and DESI DR2 at 0.487): C = u/(u+2), so the
  condition reads u/(u+2) = 1/u ⟹ u² − u − 2 = 0 ⟹ **u* = 2, C* = 0.500 exactly.**
- At γ = 2: u* ≈ 0.416, **C* ≈ 0.601**.

So at the framework's own empirically preferred γ, **C = 0.50 is exactly the point of maximum
sensitivity of coherence to log-density** — the one value in [0,1) that is dynamically distinguished.

### What this does and does not overturn

- **Does not overturn:** "ρ_crit is not a critical point." That stands — this is a saturation knee,
  there is no self-consistency loop, no free energy, no critical exponents, and the Critical
  Exponents failure is unaffected. An inflection in a monotone sigmoid is not a phase transition.
- **Does not overturn:** the C ≈ 0.50 consciousness-threshold demotion. That demotion is correct,
  but it rests on **circularity** — the eight "independent methods" inherit one unvalidated
  calibration and none of them measures C — not on geometry.
- **Does overturn:** the argument "dC/dρ is maximised at ρ = 0, therefore C = 0.50 is not
  dynamically privileged." That argument is false in the coordinate the framework actually uses,
  and it was carried on the site's `/consciousness-demo` page until 2026-09-07.

### The lesson this instance carries

This is an **over-refutation**, the same class as the a₀ "8σ", the ΔBIC = +184 effective-N inflation,
and the +17.95σ vs 8.7σ Cassini gap: a demotion argued with more force than the mathematics supports.
The direction of the error is *against* the framework, which is why it survived — the audit machinery
is tuned to catch overclaiming and does not symmetrically catch over-refuting. A reader who checks the
derivative finds the program refuting itself with false algebra, and that costs exactly what
overclaiming costs. **Corrections must be error-checked in the same direction as claims.**
