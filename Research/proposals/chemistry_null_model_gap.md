# Proposal: Chemistry Correlation Claims Require Null Model Comparison

**Date:** 2026-05-10
**Origin:** Maintainer session, triggered by Pass 3 (Graduate Physics Student) visitor feedback
**Status:** Open — Action Required Before Publishing Chemistry Results

---

## The Gap

The site's chemistry "validation" presents sound velocity (r=0.982), electronegativity (r=0.979), atomic volume (r=0.956), and ~89% of 1,703 phenomena as "Validated" against the γ ≈ 1 boundary. The relevant null model has never been computed.

**The fundamental problem:** Sound velocity, electronegativity, and atomic volume are themselves **near-monotonic functions of atomic structure**. Any smooth monotonic function of atomic number Z (or of mean atomic mass, or of electron count) will achieve r → 1 on the same dataset — by construction, not by physics.

The comparison being made is: *does the Synchronism coherence parameter C(ρ,γ) predict r = 0.98?*

The question that should be asked is: *does C(ρ,γ) predict r = 0.98 significantly better than any other smooth monotonic function fit through the same data?*

These are categorically different questions.

---

## The Null Model

To determine whether r=0.98 is evidence for Synchronism specifically (rather than for "density is related to chemical properties in a monotonic way"), the analysis needs:

1. **Polynomial null**: Fit a degree-2 or degree-3 polynomial in Z to the same 1,703 phenomena. What r does it achieve?
2. **MOND null**: Fit a MOND-type formula (a function of g/a₀) to the same phenomena. What r does it achieve?
3. **Landau null**: Fit a generic tanh function with two free parameters (not the Synchronism-specific C(ρ,γ)). What r does it achieve?

If any of these achieves r ≥ 0.98, then r=0.98 for Synchronism is not evidence of framework specificity — it is evidence that density and chemical properties are monotonically related (which was already known before Synchronism was proposed).

---

## What This Means for Site Claims

If the null comparison shows:
- **C(ρ,γ) significantly beats the null** → The "Validated" badge may be defensible, but needs the null comparison documented
- **C(ρ,γ) ties the null** → Chemistry correlations are Reparametrization: Synchronism captures the density-monotonicity that any sigmoidal function would capture
- **C(ρ,γ) beats the null by a small margin** → Chemistry is a Reparametrization of Landau-class behavior; boundary cases only differ from the null

Current best-estimate: **the null comparison will show tie or marginal win**, because:
- The framework was calibrated to chemistry data (i.e., the boundary values were fitted to make correlations high)
- The quantities with highest r (sound velocity, electronegativity) are known to be among the most monotonic with Z
- The template-bias concern (sessions 134-2660 used similar analysis patterns) suggests the sampling may preferentially find density-monotonic phenomena

---

## Proposed Resolution

**Path A: Compute the null (preferred)**
Run the polynomial/MOND/Landau null comparisons on the publicly available data. If Synchronism wins meaningfully (Δr > 0.05 against the best null), the "Validated" badges stand with this comparison documented. If not, downgrade to "Reparametrization — consistent with density-monotonicity" everywhere.

**Path B: Pre-emptive honesty (interim)**
While the null computation is pending, downgrade chemistry badges from "Validated" to "Reparametrization or Validated — pending null model comparison." This is what the site currently needs: the existing caveat language about template bias is necessary but not sufficient; it does not address the monotonicity null.

**Path C: Restrict claims to specific cases**
Identify the phenomena where the correlation is NOT monotonic in Z — cases where chemistry reverses the expected trend and Synchronism tracks the reversal. Those cases would be genuine evidence. There are likely to be few of these.

---

## Why This Matters

If the chemistry "validation" is vacuous (comparing against the wrong null), then:
1. The 89% validation rate is not evidence for the framework — it is evidence that density is a determinant of chemical properties, which is not novel
2. The "Validated" badges on r=0.98+ correlations are misleading to readers who understand what r=0 means as the null (it is not the relevant null)
3. A grad student reviewer would catch this immediately and conclude the chemistry evidence is junk — which damages the site's credibility on everything else

The site's honesty standards require that we apply the correct null model before claiming validation.

---

## Related Work

- Site `/honest-assessment` already acknowledges: "Until the rate is decomposed into independent prospective predictions, '89% validated' is better read as '89% consistent with the γ ≈ 1 boundary'"
- Site `/chemistry-limitations` documents the failure cases (melting points 53% error, critical exponents 2× off)
- Memory entry: "Validated badges keep collapsing to Reparametrization (4-of-4 audited)"
- Explorer topic `chemistry-gamma-circularity.md` is open with related concern about γ = 2/√N_corr self-correlation

This proposal is distinct from the circularity concern: the circularity question is whether γ itself is circularly defined. This null-model question is whether r=0.98 is meaningful given the baseline monotonicity.

---

## Recommended Next Steps

1. **Executor track**: Run polynomial (degree 2 and 3 in Z) regression on the same 1,703 phenomena. Report the resulting r values. Δr = r(Sync) - r(polynomial) is the only meaningful figure.
2. **Site (interim)**: Downgrade the "89% Validated" badge on `/honest-assessment` to "Reparametrization — pending null model comparison."
3. **Site (interim)**: Add one sentence to the chemistry section of honest-assessment and gamma-boundary: "The relevant null is not r=0 — it is r(polynomial in Z). Until that comparison is run, r=0.98 is evidence of density-monotonicity, not specifically of Synchronism."
