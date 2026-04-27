# Session 636: C(ρ) Is Not a Mean-Field Order Parameter — Diagnosis Sharpens

**Date**: 2026-04-27
**Type**: Theoretical audit responding to 2026-04-27 back-annotation
**Grade**: A- (clean structural argument, narrow scope)

---

## Trigger

`Research/proposals/coherence_function_meanfield_diagnosis.md` (2026-04-27). Pass 3 graduate-student visitor diagnosed three documented C(ρ) failures (β ~ 2× off, melting points 53% off, T_c 6.5× off) as one failure: "C(ρ) is a mean-field order parameter and inherits mean-field's failure modes."

The proposed remedy: Ginzburg-Landau expansion + Wilson-Fisher RG corrections to recover correct critical exponents.

The proposal also flagged a null result as productive: "If the Ginzburg-Landau expansion cannot be defined for C(ρ) (no underlying Hamiltonian)... that is the honest diagnosis."

This session pursues the structural question first: **is C(ρ) actually mean-field in the rigorous sense?**

## What "Mean-Field" Requires

A Landau / mean-field order parameter satisfies a **self-consistency equation**: the order parameter m appears on both sides of its defining relation.

Ising mean-field: m = tanh(β J z m + β h). The argument β J z m depends on m itself. Solving this transcendental equation is what produces the phase transition at T_c — m = 0 above T_c (one solution), m ≠ 0 below T_c (three solutions, two stable).

Without self-consistency, there's no phase transition mechanism. tanh of an external function is just a smooth crossover.

## What C(ρ) Is

The Synchronism coherence function:

```
C(ρ) = tanh(γ · log(ρ/ρ_crit + 1))
```

The argument of tanh is `γ · log(ρ/ρ_crit + 1)`. This is a function of:
- ρ — external density (input)
- γ — fixed parameter
- ρ_crit — fixed parameter

**C does not appear in the argument.** There is no self-consistency.

C(ρ) is tanh(f(ρ)) where f is an external function of ρ. It is **not** the solution to a self-consistent Landau equation. It is a parametrized S-curve.

## Consequences

If C(ρ) is not self-consistent, then:

1. **No phase transition in the rigorous sense.** The "transition at ρ_crit" is the inflection point of an analytic tanh, not a thermodynamic singularity. (Confirmed in S633: Taylor expansion at ρ_crit is regular polynomial; no critical exponent in the phase-transition sense.)

2. **No underlying Landau free energy.** Self-consistent order parameters arise as ∂F/∂m = 0. Without self-consistency, there's no F to extremize.

3. **No Ginzburg-Landau expansion is defined.** GL expansion presupposes a free-energy functional F[m, ∇m]. C(ρ) doesn't have one.

4. **No universality class.** Universality classes characterize critical-point fluctuations of order parameters. C(ρ) has no critical point in this sense.

5. **No Wilson-Fisher RG corrections to apply.** WF corrects mean-field exponents below upper critical dimension. C(ρ) doesn't have mean-field exponents — it doesn't have any critical exponents.

The proposal's null-result option holds: **C(ρ) has no Landau expansion because it's not a Landau order parameter.**

## What About the "β = 0.5" Prediction?

The proposal says "β_predicted ≈ 0.5" for the framework. The actual archive prediction is `β · γ = 0.5` (simulations/chemistry/critical_exponent_test.py, Chemistry Session #29).

If γ = 1 is enforced: β = 0.5 (mean-field exponent). The "2× off" claim presumes this.

But γ is itself fitted in the framework — γ = 2/√N_corr, with N_corr extracted per-system. The simulation derives γ from observed β to keep the formula β·γ = 0.5 satisfied:

> From `critical_exponent_test.py` line 122: *"Let's not assume γ - instead derive it FROM β to test consistency"* → γ_predicted = 1/(2β).

This is a fit, not a prediction. For 3D Ising β = 0.326, the formula gives γ = 1.53. For mean-field β = 0.5, γ = 1.0. For 2D Ising β = 0.125, γ = 4.0.

The framework's "critical-exponent prediction" β·γ = 0.5 has γ as a degree of freedom set per-system. **It's a one-parameter family with no test against observation.** The "factor of 2" mismatch is what you get if you also enforce γ = 1; the framework doesn't enforce that.

## What the 53% Melting Errors Actually Are

`Research/discoveries/gamma-universal-boundary.md` lists:
- Melting points: 53% average error
- Critical exponents: 2× off (the simulation-level claim above)
- T_c: not directly addressed in this honest-limitations note

These are **correlation failures**: the framework correlates γ with material properties via empirical formulas, and the formulas don't fit melting points well.

This is different from "mean-field's failure on first-order transitions." Mean-field for first-order transitions IS unreliable, but Synchronism isn't doing mean-field for melting — it's doing γ-correlation. Correlation failure ≠ mean-field failure.

The proposal's unification ("three failures are one failure: mean-field") is appealing but doesn't match the structure. The three failures are:
1. β × γ = 0.5 is a non-prediction (γ is fit) → no critical-exponent failure to fix
2. Melting points: empirical correlation failure (γ vs T_m), not a mean-field failure
3. T_c: empirical correlation failure (η vs T_c, refuted at YBCO), not a mean-field failure

The unifying diagnosis is sharper but different: **none of these are predictions in the rigorous sense.** They are empirical correlations between fitted parameters and observed quantities. Their failures are correlation failures, not theory failures, because there's no theory making the predictions.

## What This Adds

**Seven site-archive audits, all same failure mode:**

| Claim | Source | Failure |
|-------|--------|---------|
| BTFR n≈2.2 (S631) | #48 | "not rigorous" / refuted |
| α² in A (S631) | #66 | α=1.0 fiducial |
| 500 Mpc (S632) | #4 | dimensionally inconsistent |
| 80 orders (S633) | site only | range vs smoothness |
| 47 contributions (S634) | #582 says 30 | 57% overcount |
| /galaxy-rotation badge (S635) | scorecard | depends on refuted DM + uncomputed ΔBIC |
| Mean-field diagnosis of C(ρ) (S636) | this session | C(ρ) is not mean-field (no self-consistency); diagnosis is too generous |

S636 is meta on the proposal: the visitor's diagnosis ("C is mean-field") is *closer* to the truth than the site's framing ("C is universal"), but still overstates C(ρ)'s structural status. The honest disposition is "C(ρ) is a phenomenological S-curve correlated with system properties." That's weaker than mean-field, weaker than universal, and matches what the simulations actually do.

## Recommended Site-Side Actions (Operator)

1. `/chemistry-limitations` should not present C(ρ)'s failures as "mean-field failures fixable by RG." The structural claim "mean-field" is itself unsupported. Honest framing: **C(ρ) is an empirical S-curve. Its correlations with material properties are imperfect (53% melting error etc.); these are correlation imperfections, not theory failures, because no underlying theory predicts the correlations.**
2. Drop any "universality class" framing for C(ρ). It doesn't have one because it isn't an order parameter.
3. The proposal's research questions about Ginzburg-Landau expansion and upper critical dimension don't have answers because the premises don't hold. State this as the disposition.

## What I Did NOT Do

- I did not attempt the Ginzburg-Landau expansion the proposal requested. The structural argument shows it's not defined; running through the formality would be performative.
- I did not search for whether ANY reformulation of Synchronism makes C(ρ) self-consistent. That's a constructive theoretical question — if the operator wants to pursue making C(ρ) into a real order parameter, that's a separate session with a different shape (theory development, not audit).

## So What?

The mean-field critique was a step toward honesty but didn't go far enough. The framework's C(ρ) isn't mean-field; it's a fit function. The "failures" aren't mean-field failures; they're empirical-correlation failures. The proposal's null-result option (no Landau expansion exists) is the correct disposition.

Same pattern as S631–635: archive's structure is consistent with the audit; the public framing overshoots. The proposal's framing also slightly overshoots (treating C as "mean-field" rather than "S-curve"), but the proposal's null-result option already anticipated the correct disposition.

## Files

- `Research/Session636_CRho_Not_Mean_Field.md` (this document)
- Updates to SESSION_FOCUS.md
