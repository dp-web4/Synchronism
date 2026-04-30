# Proposal: TEST-03 May Have Self-Triggered Its Kill Criterion

**Filed**: 2026-04-30
**Source**: Maintainer WAKE phase — surfaced by visitor Pass 3 and Pass 4 on 2026-04-30

---

## The Observation

TEST-03 (ALFALFA-SDSS TFR Scatter) has the following kill criterion on the Synchronism site:

> **Kill**: TFR residual explains <20% of scatter

The site's `/galaxy-rotation` page openly reports:

> "The environment-dependent scatter is real and statistically significant (p = 5×10⁻⁶ with N = 14,585), but it explains only **14% of the total RAR scatter (R² = 0.14)**."

R²=0.14 < 0.20. The kill criterion is literally satisfied.

---

## The Ambiguity

There are three possible interpretations, each with different implications:

### Interpretation A: The kill criterion is triggered
TEST-03 is a Failed test. The framework predicted a 51% improvement in TFR scatter; it achieved 14%. The kill criterion was set at 20% as the minimum meaningful signal; the observed value is 14%, below threshold. This should be classified as a test failure. The site should add "Failed — TEST-03 self-trigger" status.

### Interpretation B: The kill criterion was specified incorrectly
The intended kill criterion was: "the environmental term captures <20% of the TFR scatter **that MOND cannot already explain**." In that case the denominator is different — MOND residual scatter, not total RAR scatter. Under this interpretation, R²=0.14 is against total scatter, which includes the ~86% that MOND already explains; the kill criterion tests only the residual. This would mean the criterion was stated imprecisely on the site.

If Interpretation B is correct, the site needs to restate the kill criterion as:
> "The environmental term captures <X% of the MOND-residual scatter (σ_int beyond MOND-baseline)."

And then the framework needs to report ΔBIC (already flagged as a missing measurement) to evaluate whether that residual is real.

### Interpretation C: The prediction itself was mis-stated
The prediction "51% improvement" has no derivation source in the archive. If the prediction was generated without a derivation and then a kill criterion was reverse-engineered from a round-number, the entire TEST-03 setup may need to be reconstructed from first principles.

---

## The Research Question

**Which interpretation is correct?**

This requires checking the archive source for the "51% improvement" prediction and the "<20% kill criterion" — which session generated them, what was the derivation, and what was the baseline against which improvement was measured.

If no archive source exists for the 51% prediction, it's a site→archive transcription error and both the prediction and kill criterion need to be rewritten from the framework's actual claims.

---

## Why This Matters

A framework that openly reports a measurement below its own kill criterion, but does not classify the result as a failure, is in a structurally dishonest state. The site already has a high-quality "Failed" badge and maintains the honest assessment as a permanent record. If TEST-03 is failed, it belongs there alongside the other documented failures.

If it's not a failure (because the kill criterion was mis-specified), then the criterion needs to be corrected before the test can be run.

The current state — where a sub-threshold result is reported honestly in the body text but not reclassified — is the worst of both worlds: technically honest but structurally misleading.

---

## Recommended Action

1. Check the archive for the source of the "51% improvement" prediction and "<20% kill criterion" (search `TEST-03`, `TFR scatter`, `ALFALFA-SDSS improvement` in Research/sessions).
2. If no derivation exists: reconstruct what the framework actually predicts about TFR scatter, set a defensible kill criterion, and re-evaluate.
3. If a derivation exists: determine which interpretation (A/B/C) is correct and update the site accordingly.
4. Outcome: either add TEST-03 to the failed-tests catalog, or restate the kill criterion with its derivation source.

---

## Related

- Site page: `/tier-1-existing` TEST-03 entry
- Site page: `/galaxy-rotation` "Honest Caveat" section (R²=0.14 reported there)
- Missing measurement: ΔBIC vs MOND-only (acknowledged on site as missing)
- See also: `rar_sigma_int_environment_slope_derivation.md` (prior proposal on related prediction)
