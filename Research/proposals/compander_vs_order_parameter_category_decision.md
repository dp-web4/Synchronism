# Proposal: C(ρ) Must Commit to a Category — Compander or Order Parameter

**Filed**: 2026-05-12  
**Track**: Maintainer WAKE phase  
**Trigger**: Four-persona visitor review; Grad Student and Researcher both independently identified the equivocation as load-bearing

---

## The Ambiguity

The site currently uses two incompatible framings for C(ρ):

**Frame A: Order parameter (Landau mean-field)**
- Language: "phase transition," "critical density," "critical exponents," "mean-field universality"
- Expectation: C obeys RG scaling laws, β/ν/γ exponents match universality class
- Commitment: the framework must specify which universality class and why

**Frame B: Phenomenological compander (sigmoid map)**
- Language: "motivated, not derived," "any sigmoid satisfying four constraints is equally valid," "tanh is a functional choice"
- Expectation: AIC/BIC comparison across compander family (tanh, logistic, erf, Hill, Naka–Rushton, Kubo) is the correct diagnostic
- Commitment: the framework does not inherit mean-field machinery; failures on critical exponents are category errors, not calibration misses

Both frames appear on the same pages. The deep theory pages (Coherence Function, Parameter Derivations) consistently admit Frame B. The front-of-site pitch and the Phase Transitions page use Frame A language.

---

## Why This Is a Research Decision, Not an Editorial One

The category determines what counts as failure and what counts as confirmation:

| Claim | If Order Parameter | If Compander |
|-------|-------------------|--------------|
| Critical exponents ~2× off | Failure: wrong universality class | Category error: companders don't have critical exponents |
| Any sigmoid fits equally | Failure: no privileged tanh | Correct: AIC/BIC test is what matters |
| "Phase transition at C≈0.50" | Literal: needs ξ→∞ at ρ_crit | Misleading: smooth crossover, no divergence |
| ρ_crit is "critical density" | Meaningful: inflection = critical point | Misnomer: it's the half-saturation parameter |

The deep pages have *already committed* to Frame B — they explicitly state no self-consistency loop, no Z₂ symmetry, no Hamiltonian. The critical-exponent failures are correctly diagnosed as "expected for mean-field applied outside its scope." That IS Frame B.

The front-of-site needs to commit to the same frame.

---

## The Specific Inconsistency

At γ=2, C(ρ_crit) = tanh(2·ln(2)) ≈ 0.8824, **not 0.5**.

This means ρ_crit is *not* the inflection point of the sigmoid — it is the saturation knee. In Frame A (order parameter), a critical density should be the inflection point of the order parameter curve, where susceptibility diverges. The "+1" regulator in the log shifts the sigmoid so ρ_crit is nowhere near the inflection. This asymmetry is documented on the site (Coherence Explorer displays "C at ρ_crit = 0.8824") but not interpreted.

- Frame A cannot absorb this: ρ_crit is not the critical point
- Frame B explains it naturally: the half-saturation parameter of a μ-law compander is not the inflection of the sigmoid

---

## Three Resolution Options

**Option 1: Commit to Frame B (compander)**
- Drop all phase-transition language from the front-of-site and Phase Transitions page
- Rename the page "Smooth Crossover" or "Sigmoid Mapping"
- Relabel "ρ_crit" as "half-saturation parameter" or "saturation knee"
- Add AIC/BIC compander comparison tool
- Reframe critical-exponent failures as category errors, not misses
- This is the most intellectually honest option given what the deep pages already admit

**Option 2: Commit to Frame A (order parameter)**
- Derive a self-consistency loop for C(ρ) (what determines ρ from C(ρ)?)
- Specify the symmetry class (Z₂? U(1)? SU(2)?)
- Derive critical exponents from first principles
- Predict the universality class and compare to data
- Accept that the current 2× exponent miss is a quantitative failure to be fixed, not a category error to excuse
- This requires substantially more mathematical development than the framework currently has

**Option 3: Hybrid (scope restriction)**
- C(ρ) is a compander when applied across 80 orders of magnitude
- C(ρ) behaves as an order parameter *only* at the specific calibration regime for which a mean-field description is valid (well above upper critical dimension)
- Phase-transition language restricted to that regime with explicit scope markers elsewhere
- Requires specifying which systems satisfy the scope restriction and which don't

---

## Recommendation

Option 1 (commit to Frame B) is both the most honest given existing deep-page content and the most productive research direction. The compander-class framing opens the AIC/BIC comparison tool, which is the single most useful diagnostic the site currently lacks. It also resolves the ρ_crit naming issue cleanly.

Option 2 would be more ambitious but would require a Hamiltonian and symmetry analysis not currently present. The visitor researcher and grad student both independently identified the missing self-consistency loop as the central gap.

---

## Consequence for the Site

The Phase Transitions page needs to either:
(a) Explicitly state "C(ρ) is a sigmoid compander, not an order parameter — the phase-transition analogy is motivational only. Critical exponent comparisons are category errors." OR  
(b) Present a derivation of C(ρ) as an order parameter, including a self-consistency loop, symmetry, and Hamiltonian.

Currently it does neither — it shows the critical-exponent miss and attributes it to "applying outside its intended scope," which is Frame B language without committing to Frame B.

---

## Back-Annotation Note

This proposal was triggered by the site's visitor review identifying the compander/order-parameter equivocation as load-bearing for four consecutive personas. The diagnosis is in the site's own deep pages; the proposal is to propagate that commitment to the research core documentation.
