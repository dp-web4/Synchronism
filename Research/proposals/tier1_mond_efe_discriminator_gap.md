# Proposal: Tier-1 Tests Are Not Discriminators — MOND+EFE Covers Them All

**Date**: 2026-05-13
**Source**: Visitor Pass 4 (researcher-level) — site review
**Status**: Open — requires numerical comparison to MOND+EFE predictions

---

## The Gap

After TEST-03 (galaxy rotation scatter, R²=0.14 < 0.20 kill) and TEST-04a (cosmological growth, sign-reversed) both failed, the site's remaining "Active Discriminating Tests" are:

- **TEST-01**: SPARC environment — σ_int dependence on large-scale density
- **TEST-02**: Wide binary density — anomaly amplitude correlates with local matter density
- **TEST-05**: RAR environment partition — scatter partitions by environment

These three are all **environment-dependent predictions**. The site claims they are "Synchronism-specific" and "directly distinguish from MOND." This claim is **not supported** by a numerical comparison to the alternative it invokes.

## The Problem

MOND's **External Field Effect (EFE)** has predicted environment-dependent gravitational dynamics since Bekenstein & Milgrom (1984). AQUAL and QUMOND (Milgrom 1983, Bekenstein-Milgrom 1984; see also Famaey & McGaugh 2012 reviews) explicitly make a system's internal dynamics dependent on the ambient gravitational field. Wide binaries in the presence of an external field behave differently than those in isolation — exactly what TEST-02 measures. RAR scatter partitioned by environment is a prediction of AQUAL/QUMOND EFE dynamics — exactly what TEST-05 measures.

Pittordis & Sutherland (2023) and Banik et al. (2024) have explicitly tested EFE on Gaia DR3 wide binaries. The site's own /galaxy-rotation page does not contain the strings "EFE," "External Field Effect," "AQUAL," "QUMOND," "Banik," or "Pittordis." The discriminator claim is being made without engaging the primary alternative it claims to supersede.

## What Is Needed

For TEST-01, TEST-02, and TEST-05 to survive as discriminating Tier-1 tests, the following must be computed and documented:

1. **MOND+EFE numerical prediction for the same observable.** Not a qualitative statement — a quantitative functional form. For TEST-02 (wide binary anomaly vs local density ρ_ext), MOND+EFE predicts a specific ξ(a_ext/a₀) function where a_ext = a_ext(ρ_ext). What does Synchronism predict, and at what ρ_ext does it diverge numerically from MOND+EFE?

2. **For TEST-05 (RAR environment):** AQUAL/QUMOND predict σ_int(RAR) depends on the external gravitational field. What is the predicted slope? What is Synchronism's predicted slope? Do they differ?

3. **For TEST-01 (SPARC environment):** Synchronism predicts density-dependent scatter via C(ρ_env). MOND+EFE predicts scatter via a_ext/a₀. These could be degenerate (if C(ρ_env) ∝ a_ext/a₀), numerically distinct (if the functional forms diverge), or orthogonal (if they partition the scatter differently by galaxy type). Which is it?

## Diagnostic Branches

**Branch A — Degenerate predictions:** Synchronism's environment-dependent predictions are identical in functional form to MOND+EFE on all three tests. If so, TEST-01, TEST-02, and TEST-05 cannot discriminate between the frameworks. The honest conclusion: zero remaining Tier-1 discriminators after TEST-03 and TEST-04a.

**Branch B — Numerically distinct predictions:** Synchronism predicts a different dependence on ρ_ext (or a_ext) than MOND+EFE, and the predictions diverge measurably at some parameter value accessible to Gaia DR4 or SPARC. This would constitute genuine discriminating Tier-1 content.

**Branch C — Different partitioning variable:** Synchronism's variable is local mass density ρ (measured from baryonic maps), while MOND's EFE variable is external acceleration a_ext (measured from Newtonian g of surrounding mass). For isolated galaxies these are correlated but not identical — the distinction may be empirically accessible.

## Recommended Action

1. Compute the MOND+EFE prediction for the TEST-02 wide binary observable (anomaly amplitude vs ρ_ext or a_ext) using QUMOND equations and the Banik+2024 fitting framework.

2. Compute Synchronism's prediction for the same observable using C(ρ_ext) and the current parameterization.

3. If Branch A: add explicit "MOND-degenerate" labels to TEST-01, TEST-02, TEST-05 on the site. Revise the "Active Discriminating Tests" section to show zero active discriminators.

4. If Branch B or C: document the quantitative divergence point and specify the dataset, redshift, and column required to test it.

## Why This Is Research-Critical

The Synchronism program's scientific standing rests on having tests that can distinguish it from existing theories. If all remaining active tests are MOND+EFE-degenerate, the framework has produced:

- Zero predictions that a positive experimental result would uniquely confirm
- One prediction where a negative result (DESI TEST-04a) uniquely refuted

That asymmetry — refutable but not confirmable with existing test designs — is the current honest state of the program. Documenting it is the advance; pretending otherwise is the liability.

The site's own Pass 4 reviewer noted: "The site has zero published discriminators that would in principle distinguish a positive result from MOND+EFE doing exactly what Milgrom always said it would."

---

## Related Proposals

- `session107_disfavored_by_desi_dr1.md` — TEST-04a refutation chain
- `test04a_mechanism_class_sign_failure.md` — mechanism-class diagnosis
- `rar_sigma_int_environment_slope_derivation.md` — TEST-05 slope derivation attempt
