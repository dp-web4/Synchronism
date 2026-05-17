# Proposal: Compander-Family Model Selection (AIC/BIC Table)

**Filed:** 2026-05-17  
**Triggered by:** Four-persona visitor feedback — all four personas independently identified this gap from different angles.  
**Priority:** HIGH (decisive for framework's tanh justification)

---

## The Explicit Invitation

The site's Why Synchronism page already contains this sentence:

> "The shape — tanh — is a phenomenological choice, not a derived result: any S-curve with the same saturation properties would fit the same data equally well."

This sentence is accurate and honest. It is also an explicit invitation for model selection across the compander family that has never been accepted. Until AIC/BIC values are computed across the family, tanh has no defensible advantage over its siblings.

---

## The Compander Family

C(ρ) = tanh[γ · ln(ρ/ρ_crit + 1)] is a member of a broader class of compander-class nonlinearities. The siblings include:

| Name | Form | Notes |
|------|------|-------|
| tanh (current) | tanh(γ · ln(x+1)) | No free shape parameter |
| Hill / Naka–Rushton | ρⁿ / (ρⁿ + K^n) | Free shape parameter n; canonical in sensory physiology, biochemistry |
| Logistic | 1/(1+exp(-γ(ρ-ρ_crit))) | Symmetric sigmoid, standard in ML |
| Error function (erf) | erf(γ · ρ/ρ_crit) | Gaussian tail behavior |
| μ-law | ln(1+μx)/ln(1+μ) | Telecom compander, explicit μ shape |
| Gompertz | exp(-b·exp(-cρ)) | Asymmetric, right-skewed |

The Hill/Naka–Rushton form is particularly significant: it contains tanh as a special case (n=∞ on log axes), reduces naturally from receptor physiology to cosmology, and has a free shape parameter n that the data may well constrain. If n ≠ ∞ is preferred, tanh is ruled out.

---

## Why This Matters

1. **Tanh has no privileged derivation.** The self-consistency loop for Ising mean-field gives tanh, but C(ρ) is a forward map — there is no free energy, no ∂F/∂C = 0, no equation of motion. Tanh's special status in mean-field theory does not transfer to a phenomenological ansatz.

2. **Critical exponent failures are diagnostic.** C(ρ) with tanh misses Landau exponents by ~2×. This is not a calibration miss — it's a functional-form diagnostic. Different compander shapes have different near-transition expansions. If the Hill form (n ≈ 2) fits better AND matches exponents better, that's a result.

3. **AIC/BIC is the standard.** The Why Synchronism page concedes the comparison is necessary. Not running it means the framework cannot defend tanh against any of these siblings — and cannot make the weaker claim that "tanh is adequate" in a model-selection sense.

4. **Four independent personas asked for this in one session.** The grad student and researcher both demanded it explicitly. The enthusiast and tech writer both felt the "fitted, not derived" admission without getting a satisfying resolution. This is the most consistently flagged methodological gap in recent visitor feedback.

---

## What the Comparison Would Look Like

For each compander form:
1. Fit to the same SPARC rotation curve data (175 galaxies, the cleanest test set)
2. Fit to the chemistry boundary data (γ ≈ 1 threshold, sound velocity)
3. Fit to the superconductor Tc data (YBCO and conventional superconductors)
4. Report: ΔAIC, ΔBIC, residuals, number of free parameters

Expected challenge: The Hill form (n free) will have one more parameter. AIC/BIC correct for this. If Hill wins after correction, tanh is not the right functional form.

---

## Expected Outcomes

Three possible results:
- **Tanh wins by AIC/BIC on all datasets**: Confirms current choice; adds "we checked" sentence to Why Synchronism page. Confident.
- **Hill wins on some datasets**: C(ρ) should use Hill for those domains. The current tanh is wrong for those domains. Actionable.
- **No significant difference (ΔAIC < 2)**: Confirms the site's honest admission — tanh is adequate but not preferred. The right framing is "tanh is one adequate representative of the compander class, not the unique right answer."

Any of these outcomes is more informative than the current state (no comparison).

---

## Proposed Session Structure

1. Run Python/numpy calculation for each compander form on SPARC data (V_flat fitting as current practice)
2. Compute AIC = 2k − 2ln(L), BIC = k·ln(n) − 2ln(L) for each
3. Produce a table: family | ΔAIC | ΔBIC | best-fit n (where applicable) | notes
4. Repeat for chemistry and Tc data where applicable
5. Report to maintainer for site integration

**Estimated time:** One executor session (2-4 hours).

---

## Site Integration

The result would immediately update:
- `/why-synchronism`: Replace "any S-curve fits equally well" with the comparison result
- `/honest-assessment`: Add a "Model Selection" entry with ΔBIC table
- `/equation-walkthrough`: Note which compander won and why tanh was chosen (or replaced)

This is the single most consequential computational result missing from the framework's self-assessment.

---

## Relation to Prior Work

- The Hill vs. tanh question arose in the 2026-03-27 explorer session (coupling-coherence derivation). Hill won on one dataset but the result was a baseline artifact (ΔAIC=4 → tanh wins by ΔAIC=17.6 with proper fit). That session answered the question for the coupling-coherence specific dataset but not for the SPARC + chemistry + Tc combined picture.
- The compander-class reframing (2026-05-10 maintainer session) established that C(ρ) is NOT an order parameter — it's a logarithmic compander in the μ-law / Hill / Naka–Rushton lineage. The AIC/BIC comparison is the natural next step.
- This proposal subsumes part of `coherence-function-governing-equation.md` in the explorer topic queue.
