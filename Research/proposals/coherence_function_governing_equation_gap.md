# Proposal: The Governing Equation Gap — What Equation Does C(ρ) Solve?

**Filed from**: synchronism-site maintainer session 2026-05-11  
**Origin**: Four-persona visitor review, Pass 3 and Pass 4 both independently surfaced this question  
**Status**: Open research gap — no current treatment in archive or site

---

## The Gap

The site correctly states on `/coherence-function`:

> "The tanh shape is *motivated* by mean-field theory, not derived from it… Any sigmoid satisfying the four constraints (logistic, erf, arctan, Hill) would have been an equally valid choice."

And crucially:

> "C(ρ) is not a self-consistency equation — ρ goes in, C comes out, with no feedback loop."

This is the correct description of C(ρ) as a **forward map**. The problem is that the motivation invokes mean-field theory, but the Ising/Landau mean-field tanh is specifically the solution to the self-consistency equation:

```
m = tanh(βJzm)
```

The tanh appears in Landau theory not as a postulate but as the **fixed-point solution** of this equation. It carries with it symmetry breaking, a diverging susceptibility at T_c, and critical exponents. C(ρ) shares the functional form but none of this physical content, because there is no fixed-point equation.

## The Question

**What field equation, if any, is C(ρ) the solution of?**

Three possible answers:

### Option A: There is no governing equation (C(ρ) is purely phenomenological)
C(ρ) is a postulated forward map from density to coherence — a **phenomenological compander** chosen for its saturation properties. In this case:
- The mean-field motivation is decorative (shares a form, not a derivation)
- The tanh has no privileged status over Hill / logistic / erf (as the site admits)
- The correct framing is: "We chose the tanh-log compander from the family of sigmoidal squashing functions because it has the right asymptotic behavior at both ends"
- This is a valid position — μ-law audio companding, Naka-Rushton photoreceptor response, and Hill enzyme kinetics are all honest phenomenological companders with no governing equation

### Option B: C(ρ) is the mean-field solution to a coherence self-consistency equation (not yet derived)
The conjecture: there exists an equation of the form

```
C = F[C, ρ, γ]
```

whose fixed point is C(ρ) = tanh(γ ln(ρ/ρ_crit + 1)). If this equation can be written down, it would:
- Give the tanh a physical justification (not just motivation)
- Imply specific critical exponents at the fixed point (which can be tested)
- Provide a natural definition of ρ_crit as the critical point of the self-consistency loop
- Make the parameter γ derivable rather than ansatz

**Challenge**: The "+1" log-regulator asymmetrizes the sigmoid in a way that does not arise naturally from symmetric mean-field theory. Any governing equation would need to account for why ρ = 0 maps to C = 0 (not C = -1). This may require a non-equilibrium formulation.

### Option C: C(ρ) is the solution to a dynamic equation at steady state
The framework has been criticized for lacking time evolution. If C(ρ) is the **steady-state** solution to a dynamic equation:

```
dC/dt = -∂V(C)/∂C + η(ρ, γ)
```

where V(C) is a double-well potential and η(ρ, γ) is a density-dependent driving term, then:
- C(ρ) would be the fixed point at dC/dt = 0
- Critical behavior would emerge when the potential barrier height → 0
- Time evolution would describe how systems equilibrate to C(ρ)
- This would be the kinematic layer that multiple visitor passes and explorer sessions have identified as missing

## Connection to Known Diagnoses

This gap is related to but distinct from previously filed proposals:

- **`coherence_function_landau_reduction_question.md`** (2026-04-29): Asks if C(ρ) reduces to Landau theory — answered "no" (compander, not order parameter). The governing-equation gap is the upstream question: not whether it *reduces to* Landau, but whether it *has* a Landau-like origin equation at all.

- **`dual_coherence_functions_kinematic_bifurcation.md`** (2026-05-07): Identifies the dual-C problem (C(ρ) vs C(γ,D,S)). The governing equation question applies to both branches.

- **`rho_crit_asymmetry_saturation_knee.md`** (2026-05-08): The "+1" regulator is precisely what prevents ρ_crit from being a true critical density. Any governing equation must explain why this asymmetry exists.

## Research Program

1. **Test Option A explicitly**: Can a polynomial-in-density fit achieve equal or better fit to any dataset where C(ρ) is applied? If yes, the forward-map interpretation is confirmed. (The chemistry null-model comparison would answer this for chemistry; similar null comparison needed for galaxy rotation.)

2. **Attempt Option B construction**: Write down the simplest possible mean-field self-consistency equation whose fixed point approximates tanh(γ ln(ρ/ρ_crit + 1)). If the "+1" regulator cannot be accommodated naturally, this option is effectively closed.

3. **Assess Option C feasibility**: What is the minimal dynamic equation whose steady state reproduces C(ρ)? Does it predict measurable relaxation timescales? Are those timescales consistent with known decoherence timescales in any test system?

## Why This Matters

The governing equation gap determines **whether C(ρ) is foundational or phenomenological**. If foundational (Option B or C), the framework has a dynamical content that can generate novel predictions about time evolution and critical phenomena. If phenomenological (Option A), the framework's predictive power is bounded by the compander family it belongs to — and the honest framing is "a well-chosen squashing function with physical motivation for its parameters," not "the coherence equation."

The site's /honest-assessment and /coherence-function pages have already moved most of the way to the Option A framing. The remaining step is to make this explicit at the front of site: **not "motivated by mean-field theory" (which implies shared physics) but "shares the functional form of mean-field solutions" (which is the accurate description of a compander).**

---

*Visitor-sourced research gap. Pass 3 grad student: "If C(ρ) is not a self-consistency equation, why is the tanh form privileged over Hill / logistic / erf?" Pass 4 researcher: "An open-loop forward map from density to coherence doesn't predict dynamics; it labels regimes. What predicts time evolution?"*
