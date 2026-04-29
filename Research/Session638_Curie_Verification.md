# Session 638: Verification of the Curie-Paramagnet Reduction of C(ρ)

**Date**: 2026-04-29
**Type**: Verification of site-explorer-track analysis (8th Site-Archive-Audit instance)
**Trigger**: Same-day proposals `coherence_function_landau_reduction_question.md` and `coherence_function_curie_paramagnet_reduction.md`
**Grade**: B+ (verification, not new analysis — but tightens the cosmology-quantum tracks)

---

## Setup

The morning of 2026-04-29 produced two same-day proposals:

1. **Question** (06:11): Does C(ρ) reduce to Landau theory near criticality? If yes, the chemistry corpus is validating Landau (1937), not Synchronism specifically.
2. **Answer** (08:13): No — it reduces to *less than* Landau. C(ρ) is the equilibrium response of a single binary variable in an external log-density field — a Curie paramagnet (non-interacting), without a critical point, broken symmetry, or universality class.

The answer proposal contains a complete derivation. S638 verifies it independently via computer algebra.

## What the Proposal Claims

For C(ρ) = tanh(γ · log(ρ/ρ_crit + 1)) with h(γ, ρ) = γ · log(ρ/ρ_crit + 1):

1. **Free energy structure**:
   ```
   F(C, ρ) = ((1+C)/2) ln(1+C) + ((1−C)/2) ln(1−C) − h·C
   ```
   Equilibrium ∂F/∂C = 0 gives C = tanh(h).

2. **Taylor expansion**:
   ```
   F = (1/2) C² + (1/12) C⁴ + (1/30) C⁶ + (1/56) C⁸ + (1/90) C¹⁰ + ... − h·C
   ```
   Coefficient of C^(2n) is **1/[2n(2n−1)]**.

3. **All coefficients positive** → F is convex around C=0 → no critical point.

4. **C ≥ 0 only** (h ≥ 0 since log(x+1) ≥ 0 for x ≥ 0) → no Z₂ symmetry.

5. **ρ_crit is not a critical density** (C(ρ_crit) > 0 for any γ > 0).

## Verification Results

Code: `simulations/session638_curie_verification.py`

| Claim | Method | Result |
|-------|--------|--------|
| 1. Equilibrium gives C = tanh(h) | sympy `solve` of ∂F/∂C = 0 | ✓ Returns `[tanh(h)]` |
| 2. Coefficients = 1/[2n(2n−1)] | sympy `series` to C¹², compare to expected | ✓ All 6 orders (n=1..6) match exactly: 1/2, 1/12, 1/30, 1/56, 1/90, 1/132 |
| 3. All positive | inspection of expansion | ✓ Every quadratic-and-higher coefficient is +1/[2n(2n−1)] > 0 |
| 4. C never negative | numerical evaluation across γ ∈ {0.1, 1, 2, 5} and ρ/ρ_crit ∈ [10⁻⁴, 10⁴] | ✓ Confirmed C ≥ 0 always |
| 5. C(ρ_crit) ≠ 0 | direct evaluation | ✓ At γ=2: C(ρ_crit) = 0.882; at γ=1: 0.600; at γ=0.5: 0.333 |

The math holds.

## Two Script Issues (Reporting, Not Substance)

I am noting these to keep the verification honest:

- **Sympy didn't simplify the substitution check for Claim 1.** It returned `[tanh(h)]` from `solve` (definitive), but `simplify` on the substituted expression left a residual involving `log(1±tanh(h))`. The identity `(1/2) log((1+C)/(1−C)) = arctanh(C)` requires explicit rewrite to collapse. The math is correct (the solver confirms it); sympy just didn't recognize the identity automatically.
- **My "smoothness" numerical check used a fixed-magnitude jump threshold.** For γ ≥ 1 the curve is *steep* (not discontinuous), and adjacent samples can differ by more than 0.01. The curve is smooth (analytic, in fact); my metric was inappropriate. The "false" verdicts are spurious.

Neither issue affects the verification's substantive conclusion.

## Interpretation

The Curie reduction is correct. C(ρ) is the equilibrium response of a single binary variable in an external field h(ρ) = γ·log(ρ/ρ_crit + 1):

```
F[binary] = − T · S[binary] − h · C
         = (T/2)[(1+C) ln(1+C) + (1−C) ln(1−C)] − h · C
```

Setting T=1 recovers the framework's tanh form. This is **MaxEnt over a two-state system in an external log-density field** — a textbook Curie paramagnet, with no interactions and no collective behavior.

This is structurally weaker than:
- **Mean-field Landau**: has a critical point (where the C² coefficient changes sign), broken Z₂ symmetry, defines a universality class with predicted exponents (β=1/2, etc.). C(ρ) has none of these.
- **Ising mean field with self-consistency**: m = tanh(βJzm) — the field depends on ⟨m⟩ itself. This is what the framework's notation suggests but the math doesn't implement. Dropping the self-consistency strips out the critical-point structure.

## What This Sharpens

S636 found C(ρ) is "not mean-field" via the self-consistency check (the argument depends on external ρ only, not on ⟨C⟩). S638 confirms the proposal-track answer that the system is even simpler than that — it is the *equilibrium of a single non-interacting two-state variable*, the Curie response. The 53% melting-point error, 6.5× YBCO T_c error, and 0/7 fractal-coherence-bridge are exactly what a non-interacting Curie response would produce when applied to interacting phase transitions: the model lacks the structure those phenomena require.

## Audit-Channel Taxonomy Now 8 Modes

| # | Type | Session |
|---|------|---------|
| 1 | Quantitative refutation + mislabeling | S631 |
| 2 | Dimensional inconsistency | S632 |
| 3 | Structural overclaim | S633 |
| 4 | Count discrepancy | S634 |
| 5 | Domain-level badge overclaim | S635 |
| 6 | Category error (mean-field) | S636 |
| 7 | Derivation succeeds but predicts undetectable signal | S637 |
| 8 | **External-track derivation independently verified (NEW)** | **S638** |

The 8th mode is qualitatively different: rather than the worker session producing the analysis, the site explorer track produced a complete analysis and the worker session verified it. The verification track is now operational.

## What Survives, Honestly Stated (from the proposal, with my agreement)

After this verification, the framework's substantive content is:

1. The choice of effective field h = γ · log(ρ/ρ_crit + 1) — log-density coupling is a phenomenological choice, not derived
2. The CLT-flavored γ = 2/√N_corr scaling — open whether derivable from a microscopic model
3. The cross-scale claim — true by construction (any sigmoid + log-density field will fit), weakly informative
4. TEST-04, TEST-07, TEST-02 predictions — survive as effective-field-theory predictions, not as derivations from collective coherence

## Site Action (Operator Queue)

Per the proposal's recommendations, with this verification:
- Rename ρ_crit → ρ_∗ (or annotate prominently as "field-zero offset, not critical density")
- Relabel "phase boundary" plots as "saturation curves"
- Add a `/coherence-function-as-paramagnet` derivation page
- Sharpen "89% chemistry validated" caveat to "monotonicity-across-corpus correlation universality"
- The Ising mean-field justification for tanh has already been retracted on parameter-derivations; this verification gives the precise alternative model

## Files

- `Research/Session638_Curie_Verification.md` (this document)
- `simulations/session638_curie_verification.py` (sympy + numpy verification)

## So What?

The Curie identification is the precise structural diagnosis the framework was missing. It locates the framework: a phenomenological saturation response with no microscopic basis in collective coherence. This is a clean and honest position — not destruction, but bounding. Combined with S629–S637, the framework's predictive content is now fully characterized: MOND in the testable cosmology regime (S637), Curie-paramagnet response in the chemistry/condensed-matter regime (S638). Neither contributes a discriminating experimental test against existing physics.

The path forward suggested by the proposal — adding a self-consistency loop where h depends on ⟨C⟩ — would re-introduce critical-point structure and make the framework a genuine phase-transition theory. Whether that path is worth pursuing is an operator-level decision; if pursued, it would be a constructive arc rather than a stress-test arc.
