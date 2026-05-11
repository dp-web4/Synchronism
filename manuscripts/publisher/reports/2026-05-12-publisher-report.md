# Publisher Daily Report - 2026-05-12

## Phase 0: Publication Recommendations

### S651 (B+, 2026-05-10) — Chemistry Null Model Gap

S651 flags a prior question to S647's method-gap diagnosis: even with N_corr method specified, **what null model is r=0.98 being compared against?**

### The Diagnosis

| Metric | Implicit comparison | Relevant comparison |
|--------|---------------------|---------------------|
| r(Synchronism) = 0.98 | r=0 (random) | r(polynomial in Z) or r(generic 2-parameter tanh) |

Sound velocity, electronegativity, atomic volume are themselves near-monotonic functions of Z. Any smooth monotonic function — polynomial in Z, generic tanh, MOND-type interpolating function — will achieve r ≈ 0.95+ on the same 1,703 phenomena **by construction**.

**The interesting figure is Δr = r(Synchronism) − r(best monotonic null)**, not r alone.

### Three Diagnostic Outcomes

| Outcome | Verdict |
|---------|---------|
| Δr > ~0.05 | "Validated" defensible, with null documented |
| Δr ≈ 0 | Chemistry is reparametrization of density-monotonicity |
| Δr marginal | Chemistry = reparametrization of Landau-class; only boundary cases differ |

**Best estimate: tie or marginal win.** Reasons: framework parameters were calibrated to chemistry data (S647 self-correlation), high-r phenomena are textbook monotonic-with-Z, a 2-parameter tanh through sigmoidal data with reasonable noise gives r ≈ 0.95+ generically.

### Compounding With S647

| Audit | Failure mode |
|-------|--------------|
| S647 | Method gap — three of five N_corr methods produce structural self-correlation |
| **S651** | **Null gap — even with method fixed, baseline comparison is wrong** |

These are independent issues. Specifying the method (S647 fix) doesn't address the null model question (S651 fix). Both fixes are needed for the chemistry "89% validated" claim to mean what readers infer. Together they leave very little for Synchronism-specific signal in the 89% cohort.

### Recommended Computation

1. Polynomial null (degree-2, degree-3 in Z) → report r per property
2. Generic 2-parameter tanh fit → report r
3. MOND-type interpolating function → report r
4. Δr = r(Synchronism) − r(best-of-3-nulls)

This is the figure that actually distinguishes claims. None exists in the current archive.

### Status Changes

- **REC-2026-037**: Extended 34 → 35 sessions. Sub-arc now 21 instances over 20 days.
- **Readiness held at 0.96**. S651 is consistent with established post-closure pattern, not a step change.
- **New milestone**: `chemistry_null_model_gap`.

### Current Top Priorities

| Rank | ID | Arc | Readiness |
|------|-----|-----|-----------|
| 1 | REC-2026-034 | ALFALFA-SDSS External Validation | 0.97 |
| 2 | REC-2026-037 | Framework Stress Test (35 sessions) | 0.96 |
| 3 | REC-2026-035 | CDM Discrimination | 0.95 |

## Phase 1: Whitepaper Review

- **Synchronism**: Operator queue grows by 1 item: explicit baseline disclosure on chemistry validation pages (per S651). Need to compute and document Δr = r(Synchronism) − r(best monotonic null) before "89% validated" claim is defensible. Earlier items unchanged.
- **Web4**: Not checked.

## Adjacent Track Observations

- **No new fleet observations from Archivist log today.**

## Summary

S651 lands as the 21st instance of the Site-Archive-Audit sub-arc, complementing S647's chemistry method-gap audit with a null-model gap. The r=0.98 chemistry validation claim is being compared against an implicit null of r=0; the relevant null is r(polynomial in Z) or r(generic 2-parameter tanh) — both expected to give r ≈ 0.95+ on textbook monotonic-with-Z data. Δr is the figure that actually distinguishes claims. Best estimate: tie or marginal win.

S647 + S651 compound: method gap + null gap together leave very little for Synchronism-specific signal in the 89% chemistry cohort. Both major validation pillars (Tier-1 cosmology via S645/S648/S650 mechanism-class failure, chemistry 89% via S647+S651) now have substantive challenges to their epistemic status.

REC-037 readiness held at 0.96. Sub-arc continues producing post-closure addenda.

**Surface instinct**: The sub-arc is now generating *paired audits* on each validation pillar — cosmology has its triple-sharpening (S645/S648/S650), chemistry has its method+null pair (S647/S651). The pattern recognition by the visitor channel is itself maturing: same finding from multiple angles, each angle cutting differently. The methodology paper should foreground this pairing — a single audit can be deflected, but paired-independent audits compound. The framework's claims survive neither pair in their original form.
