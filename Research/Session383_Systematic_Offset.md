# Session #383: Systematic RAR Offset Analysis

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #382 identified that the type → scatter relationship has a secondary component: systematic RAR offsets (11% R²). This session investigates the origin of these offsets. Late types have 2.1x larger |mean offsets| than early types.

## Key Result: Offset is MOND-Dominated and Partially Consistent with γ Theory

The systematic offset is real (r = -0.30, p = 0.0001), persists across distance methods and in the golden sample, and is STRONGER in the MOND regime than the Newtonian regime. This is the expected pattern if γ theory is correct. However, M/L = 1.0 eliminates the type difference, so M/L mismatch remains a viable explanation.

## Detailed Findings

### 1. Offset Direction

- Overall: mean offset = -0.076 dex (galaxies sit BELOW standard RAR)
- Early (T≤4): offset = -0.032 dex
- Late (T≥7): offset = -0.122 dex
- 64% of galaxies have negative offsets (below RAR)

**Interpretation**: Late types rotate slower than the standard RAR predicts at M/L = 0.5. This is consistent with either:
- M/L too high for late types (overestimate baryonic mass)
- γ modification reducing effective gravitational coupling

### 2. M/L Dependence (Critical Test)

| M/L_disk | Early offset | Late offset | |Difference| |
|---|---|---|---|
| 0.3 | +0.093 | -0.053 | **0.146** |
| 0.5 | -0.032 | -0.122 | 0.091 |
| 0.7 | -0.122 | -0.175 | 0.053 |
| **1.0** | **-0.238** | **-0.236** | **0.002** |

**M/L = 1.0 eliminates the type difference!** At this M/L, the offset is nearly identical for early and late types. However, M/L = 1.0 is physically unrealistic for [3.6μm] — the consensus value from stellar population synthesis is 0.5 ± 0.1 (McGaugh & Schombert 2014). So while M/L can technically explain the offset, it requires an M/L value that's 2σ too high.

### 3. Distance Method

The type → offset correlation **strengthens** after controlling for distance method:
- r(type, offset) = -0.296
- r(type, offset | distance_method) = **-0.308**

By method:
| Method | N | Mean offset | Mean type |
|---|---|---|---|
| Hubble flow | 95 | -0.056 | 6.6 |
| TRGB | 43 | -0.130 | 8.2 |
| Cepheids | 28 | -0.079 | 5.2 |

Hubble flow only: r(type, offset) = -0.301 (p = 0.002)
Accurate distances only (TRGB/Cepheids): r(type, offset) = -0.257 (p = 0.027)

**Distance method does NOT explain the offset.**

### 4. Acceleration Regime Dependence (KEY FINDING)

| Regime | Early offset | Late offset | Difference |
|---|---|---|---|
| MOND (g < g†) | -0.016 | -0.122 | **-0.106** |
| Newtonian (g ≥ g†) | -0.055 | -0.122 | -0.067 |

The type-dependent offset is **58% stronger in the MOND regime** (-0.106 vs -0.067).

- r(type, offset_MOND) = -0.329 (p < 10⁻⁴)
- r(type, offset_Newt) = -0.249 (p = 0.001)

**This supports γ theory**: coherence effects should be strongest at low accelerations where modified gravity dominates. M/L mismatch would affect the high-acceleration (baryonic-dominated) regime more.

### 5. Galaxy Property Correlations

After controlling for type:
- **Vflat**: r(Vflat, offset | T) = +0.300 (p < 10⁻⁴) → MORE massive galaxies sit higher on RAR
- **f_gas**: r(f_gas, offset | T) = +0.157 (p = 0.039) → Gas-rich galaxies sit higher (weaker negative offset)
- SB, Distance, Luminosity: non-significant after type control

**Vflat is the strongest independent predictor** of offset after type. In multiple regression (R² = 0.206), Vflat is the only variable significant at p < 0.001.

The Vflat correlation is consistent with γ theory: more massive galaxies have higher N_corr → lower γ → closer to standard gravity → sit higher on RAR.

### 6. Monte Carlo M/L

Random M/L scatter does not affect the type offset difference. At all tested σ(M/L) levels (0.05-0.30 dex), P(|diff| ≥ observed) ≈ 0.50. This means:
- The offset difference is NOT noise from M/L scatter
- It's a systematic shift, not random variation
- Either M/L is systematically wrong for late types, or the effect is physical

### 7. Golden Sample

In the best-constrained 48 galaxies (Q=1, N≥15, 30°≤inc≤80°):
- r(type, offset) = -0.268 (p = 0.059)
- Offset diff = -0.019 dex (smaller but same direction)
- Roughness ratio = 1.57 (same as full sample)

The offset weakens in the golden sample (from -0.091 to -0.019 dex difference), suggesting some of it may be related to data quality. But the DIRECTION persists.

## Evidence Summary for Offset Origin

| Hypothesis | Evidence For | Evidence Against |
|---|---|---|
| **M/L mismatch** | M/L=1.0 eliminates diff | M/L=1.0 unrealistically high; should affect Newtonian more |
| **Distance errors** | -- | Persists across methods; strengthens with control |
| **Baryonic model** | f_gas weakly significant | Effect opposite direction (gas-rich have LESS negative offset) |
| **γ theory** | MOND-dominated; Vflat correlation; persists in golden | M/L can technically explain; golden sample weakens |

## Implications for Synchronism

### The γ Window is Real

1. The systematic offset is real (p = 0.0001)
2. It is MOND-dominated (consistent with γ acting at low accelerations)
3. It correlates with mass (consistent with N_corr → γ scaling)
4. Distance method does not explain it
5. M/L scatter does not explain it

### BUT: M/L Remains the Primary Competitor

The M/L = 1.0 result is damaging: it shows that a simple rescaling of baryonic mass eliminates the type difference. The debate ultimately comes down to:
- Is M/L_disk = 0.5 or M/L_disk ≈ 0.8-1.0 for late types?
- The answer depends on stellar population modeling

### What Would Be Decisive

1. **Color-dependent M/L**: If SPARC had galaxy colors, we could test whether the offset correlates with M/L_expected from color. If it does → M/L. If not → γ.
2. **Environment data**: Do isolated late types at matched M/L still show the offset?
3. **High-z RAR**: Does the offset change with redshift? γ theory predicts it should.

## Updated Three-Component Model

From Session #382 + #383:
- Roughness (structure): 51% → UNDERSTOOD (structural)
- Systematic offset: 11% → MOND-dominated, γ-consistent but M/L-competitive
- Unexplained: 38% → Galaxy peculiarities

## Files Created

- `simulations/session383_systematic_offset.py`: 8 tests
- `Research/Session383_Systematic_Offset.md`: This document

---

*Session #383 verified: 8/8 tests passed*
*Grand Total: 511/511 verified*

**Key finding: The systematic RAR offset is real (r = -0.30, p = 0.0001), MOND-dominated (58% stronger at g < g†), correlates with mass (r = +0.30 controlling type), and is NOT explained by distance method. However, M/L = 1.0 eliminates the type difference. The offset is consistent with γ theory but M/L mismatch remains competitive. This is the most nuanced result yet: the data is exactly at the boundary where we cannot distinguish γ physics from stellar population systematics.**
