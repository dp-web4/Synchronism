# Session #560: Grand Synthesis XXVI — The Extension Arc

**Date**: 2026-02-07
**Status**: Synthesis (no tests)

## Overview

Grand Synthesis XXVI covers the Extension Arc (Sessions #556-559), the fourteenth arc of the research program. This arc pushed beyond the galaxy-level model's limits: characterizing the within-galaxy scatter structure, investigating the 4 physical outliers, and testing whether a radial gradient model could improve the correction.

## The Extension Arc (Sessions #556-559)

### Session #556: Inner RC Anatomy
- Within-galaxy scatter = 77% noise + 23% structured excess
- Inner scatter 4.5× outer, but noise identical (1.03×) — excess is physical
- AR(1) autocorrelation: lag-1 r=0.77, 4.2× run length
- c_V predicts offset gradient: r=-0.440 (p=1.2e-6)
- Only 2.3% of within-galaxy scatter pooled-predictable
- Outer radii at noise floor (scatter/noise = 1.2)
- Grade: A

### Session #558: Outlier Deep Dive
- All 4 outliers are gas-rich dwarfs (f_gas 84-98th percentile)
- Implied M/L range 0.33-0.98 (all physically reasonable)
- Distance fixes of 23-45% would suffice for each
- 3/4 positive residuals, not significantly clustered (16th percentile)
- Each unique among nearest neighbors (2-4σ)
- NGC2915 most extreme (c_V=0.29, BCD morphology)
- Grade: A-

### Session #559: Offset Gradient Model
- Gradient model reduces total RAR scatter 5.6% (LOO)
- BUT worsens outer radii 70% (0.045 → 0.076 dex)
- Gradient genuinely predictable (LOO R²=0.428, c_V r=-0.462)
- Gradient orthogonal to offset (ΔLOO=-0.0006 as 7th variable)
- AIC/BIC favor gradient, but it destroys the outer noise floor
- Tradeoff (better inner, worse outer), not pure improvement
- Physics: M/L gradient within galaxies
- Grade: A

## Arc Synthesis: The Model's Boundaries

The Extension Arc answers three questions about the model's limits:

### Q1: What creates the within-galaxy scatter?
**Answer**: 77% measurement noise, 23% structured physical excess. The structured component is spatially correlated (AR(1) with r=0.77) and concentrated at inner radii, where non-circular motions, mass model decomposition errors, and bar/spiral effects dominate. It is NOT predictable from galaxy-level properties (2.3% pooled R²).

### Q2: Why does the model fail for 4 galaxies?
**Answer**: All 4 are gas-rich dwarfs where measurement systematics are largest. Each could be explained by plausible distance errors (23-45%), M/L errors (factor 0.7-2.0), or inclination errors (5-15°). They are at the extreme tail of the f_gas distribution, where the linear gas correction may saturate. They are not a cluster but occupy the "extreme dwarf" corner.

### Q3: Can we go beyond a constant offset?
**Answer**: Yes, but at a cost. The gradient model (LOO R²=0.428) improves total scatter by 5.6% and inner radii by 10.9%, but worsens outer radii by 70%. The gradient contains real physics (within-galaxy M/L gradients) but is orthogonal to the offset — it doesn't improve galaxy-level predictions at all. The choice between constant and gradient models depends on whether you value total scatter (gradient) or outer-radius precision (constant).

## The Complete Arc Structure

| # | Arc | Sessions | Key Finding |
|---|-----|----------|-------------|
| 1 | Research | 433-438 | Initial exploration |
| 2 | Foundation | 441-448 | Core variables identified |
| 3 | Variable Selection | 449-454 | 6-var model built |
| 4 | Robustness | 455-460 | Model validated |
| 5 | Derivation | 525-527 | Coefficients MOND-derivable |
| 6 | Interpretation | 526-530 | Offset = M/L correction |
| 7 | Rehabilitation | 528-532 | N_corr sign fixed |
| 8 | Perspective | 533-538 | Boost/offset architecture |
| 9 | Resolution | 539-541 | No missed variables |
| 10 | Diagnostic | 542-545 | Bulk-driven, MOND>CDM |
| 11 | Information | 547-550 | logL irreplaceable, distance needed |
| 12 | Validation | 551-553 | Errors random, universal predictor |
| 13 | Characterization | 554-557 | Model fingerprint |
| **14** | **Extension** | **556-559** | **Beyond galaxy-level limits** |

## The Research Program: State of Knowledge

After 160 sessions and 14 arcs, the research program has established:

**The Model**: offset = -3.379 + 1.897×logV - 0.548×logL - 0.218×c_V - 0.451×f_gas + 0.147×logV×c_V + 0.181×logL×f_gas

**What It Is**: A MOND M/L correction. All 6 signs predicted by MOND. The V-L ratio [4.01, 4.27] brackets MOND's prediction of 4.0. The model corrects luminosity → mass, not M/L directly.

**How Well It Works**: LOO R² = 0.938 [0.893, 0.965]. 67% of galaxies within 10%. Outer RAR at 1.14× noise.

**Where It Fails**: 4 gas-rich dwarfs (3.1% of sample), all with plausible measurement explanations. Inner radii have 4.5× structural excess that no galaxy-level model can fix.

**What It Can't Do**: Predict within-galaxy structure (2.3% pooled R²). The gradient model is a tradeoff, not an improvement. The model has extracted all galaxy-level information from SPARC.

## Open Directions

1. **Cross-dataset validation**: Apply to LITTLE THINGS, THINGS, or other HI surveys
2. **Theoretical derivation**: Complete MOND derivation of all 6 coefficients
3. **Cosmological context**: Connect the model to halo abundance matching
4. **Machine learning exploration**: Whether nonlinear methods can extract more from the existing features (Session #495 says no, but with more data...)
5. **The gradient question**: Whether radius-resolved M/L fitting could recover the gradient model's predictions from first principles

## Files Created

- `Research/Session560_Grand_Synthesis_XXVI.md`: This document

---

*Session #560: Grand Synthesis XXVI*
*Fourteen arcs complete*
*Grand Total: 1645/1645 verified*

**Key finding: Extension Arc (Sessions #556-559) pushes beyond galaxy-level limits. Inner scatter is physical (4.5× outer, noise 1.03×). Outliers are gas-rich dwarfs with plausible errors. Gradient model improves total 5.6% but worsens outer 70% — tradeoff not improvement. Model has extracted all galaxy-level information from SPARC. Fourteen arcs complete.**
