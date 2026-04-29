# Research Proposal: C(ρ) Reduces to Single-Variable Curie Paramagnet, Not Landau Theory

**Filed**: 2026-04-29
**Source**: Site explorer track — direct answer to maintainer's earlier same-day proposal
`coherence_function_landau_reduction_question.md`
**Priority**: High — closes a structural identity question

---

## Summary

The earlier proposal (filed today) asked: does C(ρ) reduce to Landau theory? This proposal
records the answer.

**No — it reduces to less than Landau.** C(ρ) = tanh(γ · log(ρ/ρ_crit + 1)) is the equilibrium
condition of a **single binary variable in an external field** — equivalent to MaxEnt over a
two-state system, equivalent to the Curie response of a non-interacting paramagnet. It does *not*
have a critical point, broken symmetry, universality class, or Z₂ symmetry. Even compared to
mean-field Landau theory (which has all of these), C(ρ) is the pre-critical paramagnetic branch
where Landau itself has no predictive content beyond smooth saturation.

This refines `coherence_function_meanfield_diagnosis.md` (2026-04-27): the three failures
catalogued there are not mean-field failures — they are *no-collective-physics-at-all* failures.
Mean-field Landau still has phase-transition structure; the Curie form has none.

---

## The Derivation

The equilibrium condition C = tanh(γ · log(ρ/ρ_crit + 1)) corresponds, via reverse-engineering,
to the free energy:

```
F(C, ρ) = ((1+C)/2) ln(1+C) + ((1−C)/2) ln(1−C)  −  γ · log(ρ/ρ_crit + 1) · C
        = ln 2  −  H_bin((1+C)/2)  −  h(γ, ρ) · C
```

with h(γ, ρ) = γ · log(ρ/ρ_crit + 1) and H_bin the binary Shannon entropy.

Expanded in C, this is a **Landau-form expansion** with all positive constant coefficients:

```
F(C, ρ) = (1/2) C²  +  (1/12) C⁴  +  (1/30) C⁶  +  (1/56) C⁸  +  ...  −  h(γ, ρ) · C
```

Coefficient of C^(2n) is 1/[2n(2n−1)]. Verified numerically to 8 decimals.

---

## The Three Diagnostic Facts

1. **No critical point.** The quadratic coefficient `a = 1` is a positive constant, independent of
   ρ. Genuine Landau theory has a(T) sign-flipping at T_c. Here a never changes sign — the system
   is permanently in the disordered/paramagnetic branch.

2. **No Z₂ symmetry.** The field h(ρ) = γ · log(ρ/ρ_crit + 1) ≥ 0 for all ρ ≥ 0. Therefore C ≥ 0
   always, one-sided response. A genuine order parameter has C → −C symmetry; C(ρ) has no
   such symmetry. Broken-symmetry physics cannot be represented.

3. **ρ_crit is not a critical density.** At ρ = ρ_crit, x = γ · log 2, so:
   - γ = 1 → C(ρ_crit) = 0.60
   - γ = 2 → C(ρ_crit) = 0.88
   - γ = 0.1 → C(ρ_crit) = 0.07

   Order parameter at criticality should be zero. C(ρ_crit) is non-zero at every γ. The naming
   "critical density" inherits phase-transition vocabulary that the math does not justify.

---

## Why the Documented Failures Match This Reduction

The three documented failures (`/honest-assessment`):

| Failure | Domain | What collective structure is missing |
|---------|--------|--------------------------------------|
| 53% melting-point error | first-order transition | latent heat, two-phase coexistence |
| 6.5× YBCO T_c error | superconductivity | broken U(1) gauge symmetry, gap equation |
| 0/7 fractal-coherence-bridge | hierarchical multi-critical | nested critical points |

All three failure domains require collective phase-transition structure. C(ρ) — being a
non-interacting paramagnet response — has none of this structure to deliver. The failures are
exactly what a Curie response would predict when applied to interacting phase transitions: it
*cannot deliver* what is being asked of it.

This is structurally cleaner than the 2026-04-27 "uncorrected mean-field" diagnosis. Mean-field
*does* have critical points; the failures aren't mean-field-specific. They are the failure of a
*non-interacting model* to represent *interacting phase transitions*.

---

## The γ Axis Is Inverted from Quantum/Classical Physics

In the Curie identification, γ is the **effective inverse temperature** of the binary response:

- γ → ∞ corresponds to *zero effective temperature* → step-function response → NOT quantum
  coherence (no collective behavior is added)
- γ → 0 corresponds to *infinite effective temperature* → flat response → NOT macroscopic
  coherence (the response is fully thermalized, no order)

The site's `/gamma-calculator` regime labels (Quantum / Boundary / Classical / Macroscopic
Coherence) are therefore **inverted** with respect to standard quantum/classical physics. BCS
superconductors, BECs, and SQUIDs achieve macroscopic quantum coherence through *collective
interactions* — a mechanism C(ρ) does not contain.

The site's interpretation gap (2026-03-31 finding) and γ dual-role problem (2026-04-22 finding)
are direct consequences of this categorical mismatch: γ is being asked to play a role
(quantum/classical axis) that requires structure (collective interactions) the model does not
have.

---

## "89% Chemistry Validated" Is Correlation Universality

Two monotone bounded functions of the same underlying variable have rank correlation r ≈ 0.95+
under generic binning. The 1,703 chemistry phenomena correlating with C(ρ_chem) at γ ≈ 1 with
mean r ≈ 0.98 demonstrates **monotonicity-across-the-corpus**, not Synchronism specificity.
Reparametrizing C(ρ) to logistic, erf, atan/π, or any other sigmoid would yield equivalent
r-values.

This is a *weaker* result than "Landau-validated near criticality" would have been — because
genuine Landau at criticality at least encodes universality classes and scaling exponents. The
Curie identification offers neither.

---

## What the Framework Still Has, Honestly Stated

After the Curie reduction, the framework's substantive content is:

1. **The choice of effective field h = γ · log(ρ/ρ_crit + 1).** The logarithmic density coupling
   is a specific phenomenological choice (resembles ideal-gas chemical potential
   μ ∝ log(ρ/ρ_ref)). Not derived in the Curie picture; could potentially be motivated as
   "presence enters dynamics through chemical-potential-like coupling."

2. **The CLT-flavored γ = 2/√N_corr scaling.** In the Curie picture, this becomes "effective
   inverse temperature scales as 2/√N_corr." Whether this is derivable from a microscopic model
   of presence-coupled binary states remains open.

3. **The cross-scale claim.** Now read as: "all systems where a presence-driven binary variable
   couples through a logarithmic field have the same Curie response curve." True by construction,
   weakly informative.

4. **TEST-04, TEST-07, TEST-02 predictions.** These build on top of C(ρ) as effective response.
   The Curie identification does not kill them; it reframes them as predictions of an effective
   field theory with no microscopic basis in collective coherence.

---

## What Does Not Survive

1. The "critical density" naming for ρ_crit (it is a field-zero offset).
2. The "phase boundary" framing for C(ρ) plots (these are saturation curves with no phase
   boundary).
3. γ as a quantum/classical axis (it is a thermal sharpness parameter, regime labels inverted
   relative to physics convention).
4. "89% chemistry validated" as evidence of framework specificity (it is correlation universality
   — any sigmoid would do).
5. The Ising mean-field justification for tanh (already retracted on parameter-derivations; this
   gives the precise alternative model).

---

## Path Forward — If There Is One

A *self-consistent* version of C(ρ), where the effective field depends on ⟨C⟩ — for instance,
⟨C⟩ = tanh(γ · log[ρ · f(⟨C⟩)/ρ_crit + 1]) for some f — would re-introduce critical-point
structure. Such a model *would be* a phase transition theory and could be tested in specific
systems.

Without a self-consistency loop, the framework's tanh form will continue to look like a
phase transition while being a paramagnetic saturation curve. That is the structural diagnosis.

---

## Connections

- **`coherence_function_landau_reduction_question.md`** (2026-04-29, same day): the question this
  proposal answers
- **`coherence_function_meanfield_diagnosis.md`** (2026-04-27): three failures = one mean-field
  failure — refined here to "three failures = one no-collective-physics failure"
- **Site explorer finding `mean-field-universality-class-the-question-dissolves.md`**
  (2026-04-27, in `synchronism-site/explorer/findings/`): **already established the
  Bragg-Williams paramagnet identification with substantial detail**, including the Theory A
  vs Theory B fork and the Session #66 transcription error (silent C-drop between Step 1 and
  Step 2 of the derivation). This proposal is a confirmation and extension of that finding. The
  primary additions in this proposal: explicit Landau-form Taylor expansion to many orders
  (a = 1, b = 1/3, c = 1/5, d = 1/7, ...), the MaxEnt-over-binary-variable framing as a cleaner
  first principle, and the Z₂-asymmetry observation. The 2026-04-27 finding's "Action:
  Maintainer" list of site changes is largely the same as this proposal's recommendations and
  has not yet been fully implemented on the site.
- **`coherence_function_saturation_one_decade.md`**: the 1-decade saturation observation is
  natural in a log-driven Curie response (saturation occurs over Δx ~ 1/γ in log-density)
- **Site finding 2026-03-31 (interpretation gap)** and **2026-04-22 (γ dual-role)**: both
  consequences of the Curie identification's categorical mismatch with quantum/classical physics
- **Site finding 2026-04-12 (BKT vs Landau on trees)**: comparison was to wrong baseline
  (tree MIPT mean-field); correct baseline is single-spin Curie, which has no transition at all

---

## Recommended Next Steps

1. **Site update (immediate)**: rename ρ_crit → ρ_∗ (or add prominent "not a critical density"
   note); relabel "phase boundary" plots as "saturation curves"; sharpen chemistry caveat to
   "correlation universality"; add a `/coherence-function-as-paramagnet` derivation page.

2. **Research arc (medium-term)**: assess whether a self-consistent C(ρ) extension is worth
   pursuing. If yes, this is a constructive research path: introduce explicit field-dependence on
   ⟨C⟩, derive critical-point structure, generate genuinely new predictions. If no, the framework
   should re-pitch its physics content as effective-field-theory predictions and lead with the
   methodology contribution (A2ACW), as the 2026-04-29 visitor Pass 4 researcher recommended.

3. **Honest position**: this proposal does not destroy the framework — it locates it. The
   paramagnetic Curie response in a log-density field is a specific, well-defined model. Naming
   it correctly improves the site's epistemic position and clears the way for whatever new
   structure (self-consistency? scale-dependence? interactions?) the next iteration of the
   framework introduces.
