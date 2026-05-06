# Session 644: ρ_crit — Calibration Consistency, Not Independent Prediction

**Date**: 2026-05-06
**Type**: Site-Archive-Audit (14th instance, post-arc-closure)
**Trigger**: 2026-05-05 proposal `rho_crit_derivation_calibration_vs_prediction.md`
**Grade**: B+ (structural clarification with concrete next-step recommendation)

---

## Setup

The Pass 4 visitor flagged: ρ_crit = A · V_flat² takes V_flat as **input**. The "5% agreement" of theoretical A (0.0294) with empirical A (0.028) is calibration consistency — not predictive success. To convert this into a real test, the framework needs to predict V_flat (or ρ_crit) from independent observables.

S644 confirms the diagnosis from the archive and recommends Path C (independent β_J measurement from velocity dispersion).

## The Structural Issue

From Session #66 (and confirmed by S631 audit):

```
ρ_crit = A · V_flat²
A = 4π / (β_J² · G · R₀²) ≈ 0.0294
```

with β_J ≈ 1.1 calibrated from SPARC (Jeans-length / half-light radius), R₀ = 8 kpc chosen reference, V_flat measured per galaxy.

The "validation" comes from comparing:
- **Theoretical A**: computed from G, β_J=1, R₀=8 kpc → 0.0294
- **Empirical A**: fitted from observed (ρ, V_flat) pairs → 0.028

These match to 5%. But the inputs to "theoretical A" (β_J, R₀) are themselves calibrated to SPARC data. The agreement is internal consistency between two parameterizations of the same data, not an independent prediction.

**No closed predictive loop exists**: the formula consumes V_flat (as input) and produces ρ_crit (as output). Nothing in the chain predicts V_flat from independent observables, then uses that predicted V_flat to test ρ_crit.

## Why This Is the Same Pattern as Prior Audits

The proposal's diagnosis maps onto the kinematic-layer gap (S641, S642):
- ρ_crit's "derivation" assumes V_flat is given (treats it as kinematic input)
- The kinematic substrate (what determines V_flat from underlying physics) is unspecified
- β_J is fitted, not derived; R₀ is chosen, not derived

The Jeans criterion is dimensionally correct — it gives the right *form* of the relationship between density and velocity scale. But the *constant* (A ≈ 0.028) requires three calibrated inputs (β_J, R₀, and the implicit choice that V_flat is the right velocity scale to use, not σ or some other dispersion measure).

## What Would Make This a Prediction

The proposal's three paths:

**Path A (BTFR-based)**: Use stellar mass + gas fraction → predict V_flat via Tully-Fisher → predict ρ_crit → check density profiles. Partially overlaps TEST-09 (BTFR audit, S631). Does not really independent — TFR is itself calibrated.

**Path B (cosmological)**: Try to derive ρ_crit from H₀ and z_form alone. Dimensional analysis shows V_flat² will reappear via the only velocity scale available. Equivalent in structure to the MOND coincidence a₀ ≈ cH₀/(2π), which the framework correctly classifies as dimensional analysis, not derivation. Likely cannot escape circularity without a kinematic layer.

**Path C (velocity dispersion)**: β_J is the ratio λ_Jeans / R_half. λ_Jeans depends on density and velocity dispersion σ. R_half is photometric. If σ is measured *independently* of V_flat (not derived from rotation), then β_J becomes an independent prediction. Then ρ_crit = V_flat² / (G·β_J²·R_half²) becomes a real test: does it produce A = 0.028 ± 0.003 across SPARC galaxies *without* circularity?

**Path C is the genuine independent-prediction route** if SPARC σ data is available. Cost is $0; novelty is high. This is the framework's first opportunity to convert a calibration into a prediction at the galactic scale.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 13 | Regime-label inversion under operational formula | S643 |
| 14 | **Calibration consistency presented as independent prediction (NEW)** | **S644** |

S639 found two metrics under one label. S643 found two formulas under one label. S644 finds **one formula with one observable on both sides** — the input gets dressed as if it were prediction. Same audit-channel mechanism: a presentation strengthens what the math actually delivers.

## Connection to S631 (Original α² Audit)

S631 established that α in A = 4π/(α² G R₀²) is a fiducial parameter (α = 1.0), not the fine-structure constant. S644 extends: even with α correctly read as fiducial, the formula still requires V_flat as input, so it doesn't independently predict ρ_crit. The α relabeling fix (recommended in S631) is necessary but not sufficient — the structural calibration-vs-prediction issue remains after that fix.

## Recommended Site Action

**Short-term**:
- Relabel `/parameter-derivations` ρ_crit row from "Validated | Jeans Criterion | 5% Agreement" to "Semi-derived | Jeans Criterion | Calibration Consistent" or "Internally Consistent | Calibration to SPARC."
- Add explicit note: "V_flat is an input to the formula. The 5% agreement between theoretical and empirical A is consistency between two parameterizations of SPARC data, not an independent prediction."

**Medium-term (research direction)**:
- Implement Path C using SPARC velocity dispersion data where available. This is the cleanest path to an independent test of A.
- A Path C result would either:
  - Confirm A = 0.028 ± 0.003 from σ-derived β_J → first genuine prediction at this scale.
  - Find different A from σ-derived β_J → reveals which velocity scale the framework actually uses, and clarifies whether σ or V_flat is the kinematic input.

**Long-term**:
- The kinematic-layer question (S641, S642): until N_corr and β_J have first-principles derivations, ρ_crit remains a calibrated parameter. This is consistent with the framework's overall structural status as a parameterization without an action principle (S642).

## Files

- `Research/Session644_RhoCrit_Calibration_Vs_Prediction.md` (this document)
- No simulation needed — definitional clarification

## So What?

ρ_crit's "5% agreement" is real — but it's agreement between two parameterizations of the same data, not an independent prediction. The site's "Validated" badge overstates the claim. Path C (use σ, not V_flat, to get β_J) is the first available route to converting calibration into prediction. Whether it succeeds is an empirical question that SPARC data can answer at zero cost.

This is the 14th audit instance and follows the same pattern: site claims more than archive supports, the gap can be located precisely, and a corrective action is concrete (relabel, plus an optional research path). The kinematic-layer gap continues to be the single root that all 14 audits trace to.

One proposal remains pending: session107 DESI DR1 (2026-05-05).
