# Session #427: Radially-Resolved RAR Correction

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Sessions 417 and 422 showed the R_eff effect amplifies outward and c_V predicts the inner offset. This session builds a radially-resolved correction model that predicts the RAR residual at each radius, not just a constant per galaxy.

## Central Result: Galaxy-Level Correction Dominates; Radial Resolution Adds Modestly

| Model | Point-level RMS | CV RMSE | Improvement |
|-------|----------------|---------|-------------|
| Standard RAR | 0.236 | 0.236 | — |
| V+R+c_V (constant per galaxy) | 0.162 | 0.166 | **30%** |
| V+R+c_V + radial resolved | 0.159 | 0.161 | +2% additional |
| Full model (V+R+L+c_V + radial + interactions) | 0.137 | — | In-sample only |

The galaxy-level correction provides the bulk of the improvement. Adding radial resolution yields only 2.2% additional improvement in cross-validation.

## Key Findings

### 1. Variance Decomposition (Test 8)

| Source | Variance share |
|--------|---------------|
| Between-galaxy | **63%** |
| Within-galaxy | 37% |

Nearly two-thirds of the point-level RAR scatter in the MOND regime comes from systematic galaxy-to-galaxy differences — exactly what V+R+c_V captures.

### 2. c_V Predicts Radial Slope (Test 4)

Each galaxy's offset has a radial slope (d_offset/d_log(r/R_eff)). c_V is the strongest predictor:
- r(slope, c_V) = **-0.59**
- r(slope, V) = +0.16, r(slope, R) = -0.11, r(slope, L) = -0.11
- V+R+c_V → slope: R² = 0.42

The c_V coefficient for slope is -1.12: high-c_V galaxies have **negative slopes** (offset decreases outward), while low-c_V galaxies have positive slopes (offset increases outward). This is physically consistent with Session 422's finding that c_V predicts inner offset.

### 3. The Full Radially-Resolved Model (Test 5)

```
resid = -3.89 + 1.80×logV - 0.27×logR - 0.26×logL + 0.79×c_V
        + 0.48×log(r/R_eff) + 0.15×log(r/R_eff)×logR - 0.66×log(r/R_eff)×c_V
```

The c_V × log(r/R_eff) interaction (-0.66) is the most important: c_V's positive effect on offset reverses at large radii. This means:
- Inner regions: c_V → positive offset (concentrated mass → more acceleration)
- Outer regions: c_V effect disappears (consistent with Session 422: r = -0.003)

### 4. Point-Level Improvement (Tests 2, 7)

| Step | RMS (dex) | Reduction |
|------|-----------|-----------|
| Standard RAR | 0.236 | — |
| + V+R+c_V (galaxy-level) | 0.162 | 31% |
| + V+R+L+c_V (galaxy-level) | 0.150 | 37% |
| + Full radial model | 0.137 | 42% |

## Why Radial Resolution Adds Little in CV

The 37% within-galaxy variance is largely **measurement noise** — variations from point to point within a single galaxy's rotation curve that cannot be predicted from galaxy-level properties. The true systematic radial structure (c_V × r/R_eff interaction) is small compared to this noise.

This is consistent with Session 417's finding that inter-galaxy scatter dominates (inter/intra ratio = 1.70).

## Physical Interpretation

The result confirms that the RAR offset is primarily a **galaxy-level** property, not a radially-varying one. The outward amplification (Session 417) is real but reflects the fact that inner regions (where the RAR is more Newtonian and accurate) dilute the galaxy-level MOND-regime signal, not that the offset itself changes dramatically with radius.

The c_V × radial interaction reveals that concentrated-mass galaxies (high c_V) have a distinct radial profile: positive offset in the inner region that diminishes outward. Extended, diffuse galaxies (low c_V) show the opposite: their negative offset grows outward.

## Grade: B+

A useful and informative session. The variance decomposition (63/37) is an important quantitative finding. The c_V → slope correlation (r = -0.59) deepens our understanding. However, the main practical conclusion (radial resolution adds only 2.2% in CV) is a negative result — the galaxy-level model is sufficient. This is valuable to know but less exciting than a positive discovery.

## Files Created

- `simulations/session427_radial_correction.py`: 8 tests
- `Research/Session427_Radial_Correction.md`: This document

---

*Session #427 verified: 8/8 tests passed*
*Grand Total: 805/805 verified*

**Key finding: 63% of point-level RAR scatter is between-galaxy (captured by V+R+c_V). Radial resolution adds only 2.2% in CV. c_V strongly predicts the radial slope (r = -0.59): concentrated galaxies have decreasing offset outward, diffuse galaxies increasing. Full radial model RMS = 0.137 but mostly in-sample improvement. Galaxy-level correction is sufficient. Grade B+.**
