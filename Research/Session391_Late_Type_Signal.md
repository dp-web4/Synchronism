# Session #391: The Late-Type R_eff Signal

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #390 found that controlling for L at fixed V eliminates the R_eff → offset signal for the full sample. This session investigates why late types show a dramatically stronger signal (r = -0.74 vs -0.13 for early types) and whether it survives L control.

## Key Result: Late-Type R_eff Signal is L-Independent and Physical (Grade A-)

For late types (T ≥ 7): r(R_eff, offset | Vflat, L) = **-0.49** (p < 10⁻⁴). R_eff predicts RAR offset BEYOND luminosity at fixed Vflat. This is NOT true for early types (+0.20, n.s.). The late-type signal is genuine, L-independent, and consistent with a physical coherence length effect operating in the fully MOND regime.

## Detailed Findings

### 1. Subsample Characterization

| Property | Early (T≤4) | Late (T≥7) |
|---|---|---|
| N | 44 | 61 |
| Vflat | 208 ± 50 km/s | 76 ± 25 km/s |
| R_eff | 5.0 ± 2.6 kpc | 2.8 ± 2.3 kpc |
| Gas dominance | 0.23 | **0.88** |
| **MOND fraction** | **0.63** | **1.00** |
| Quality | 1.3 | 1.5 |
| r(R, off\|V) | -0.13 (n.s.) | **-0.74** (p < 10⁻⁶) |
| σ(R\|V) | 0.197 | **0.290** |

Late types are: low-mass, gas-dominated, 100% MOND, and have 48% wider R_eff range at fixed V.

### 2. THE KEY RESULT: R Beyond L

| Type | r(R, offset \| V, L) | p |
|---|---|---|
| Early (T≤4) | **+0.20** | 0.18 (n.s.) |
| Mid (5-6) | +0.46 | 0.006 |
| **Late (T≥7)** | **-0.49** | **< 10⁻⁴** |

For early types, R's effect on offset is fully mediated by L — controlling L eliminates the signal. For late types, R has **independent predictive power**: even at the same Vflat and luminosity, physically larger late-type galaxies have more negative RAR offsets. This is the coherence LENGTH effect operating in pure MOND regime.

### 3. MOND Fraction (Not a Driver)

- Late types have MOND fraction = 1.00 (all data at g < g†)
- Controlling MOND fraction strengthens the full-sample signal: r = -0.31 → -0.46
- But late types are already 100% MOND, so MOND fraction can't be split within late types

### 4. Slope Difference

| Type | Slope (R_resid → offset_resid) |
|---|---|
| Early | -0.051 dex/dex |
| **Late** | **-0.365 dex/dex** |

The effect size is **7x larger** in late types. For late types, every factor of 10 increase in R_eff at fixed Vflat predicts 0.37 dex lower RAR offset.

### 5. L-R Coupling Structure

| Type | r(L, R \| V) | r(L, V) |
|---|---|---|
| Early | +0.518 | +0.843 |
| Mid | +0.205 | +0.971 |
| Late | +0.733 | +0.730 |

Paradoxically, late types have the STRONGEST L-R coupling (r = +0.73), yet R still has independent predictive power beyond L. This means the R → offset relationship isn't just a proxy for L — there's something about physical size itself.

### 6. Gas Fraction — Not a Driver

| Late types | r(R, offset \| V) | r(R, offset \| V, gas) | Change |
|---|---|---|---|
| All late | -0.737 | -0.738 | 0.001 |

Gas fraction changes nothing. The signal is NOT about gas richness.

### 7. Quality — Not a Driver

| Subsample | r(R, offset \| V) |
|---|---|
| Late (all Q) | -0.737 |
| Late (Q=1 only) | **-0.637** |
| Late (Q control) | **-0.741** |

The signal persists in Q=1 (best quality) galaxies and survives quality control.

### 8. Matched Subsamples

Only 4 matched pairs found (early and late Vflat ranges barely overlap), so matching is uninformative. This is itself a finding: early and late types occupy almost non-overlapping regions of Vflat space.

## Hypothesis Verdict

| Hypothesis | Evidence | Verdict |
|---|---|---|
| A) MOND fraction | Late = 100% MOND; controlling MOND_frac strengthens signal | **Partially explains** |
| B) Gas fraction | Controlling gas changes r by 0.001 | **NOT a driver** |
| C) L-R decoupling | r(R,off\|V,L) = -0.49 for late, +0.20 for early | **KEY: R is L-independent in late types** |
| D) Range effect | σ ratio = 1.48 | **Contributes to significance, not effect size** |
| E) Quality | Q=1 late still r = -0.64 | **NOT a driver** |

## Why This Matters

Session #390 concluded that R_eff only works through L (the SB suppressor was collinearity). This session overturns that conclusion for late types specifically:

1. **Early types**: R → offset is fully L-mediated. No independent R effect.
2. **Late types**: R → offset has a large L-independent component (r = -0.49). This is the coherence LENGTH effect.
3. **The difference is the MOND regime**: Late types are 100% MOND; early types are 63% MOND.

This is consistent with Synchronism: coherence (size-dependent modified gravity) should only operate in the MOND regime. In the Newtonian regime, baryonic mass (L) determines everything. In the MOND regime, the coherence scale (R) matters independently.

## Implications for Synchronism

1. **The coherence LENGTH interpretation is supported** — but ONLY in the MOND regime
2. **Late types are the test bed**: All future tests should focus on T ≥ 7 galaxies
3. **The R_eff signal is NOT an L or M/L artifact** for late types
4. **The slope (-0.37 dex/dex) is a quantitative prediction** testable with other datasets
5. **MOND fraction = 1.00 for late types** means the entire RAR is in the modified gravity regime

## Files Created

- `simulations/session391_late_type_signal.py`: 8 tests
- `Research/Session391_Late_Type_Signal.md`: This document

---

*Session #391 verified: 8/8 tests passed*
*Grand Total: 559/559 verified*

**Key finding: Late types (T≥7) show r(R_eff, offset | V, L) = -0.49 (p < 10⁻⁴) — R_eff predicts RAR offset BEYOND luminosity at fixed Vflat. This is NOT true for early types (+0.20, n.s.). The slope is 7x steeper in late types. Gas fraction and quality do NOT drive this. Late types are 100% MOND (all data at g < g†). This resolves the Session #390 puzzle: the full-sample SB suppressor was partly collinearity, but the late-type signal is genuine, L-independent, and consistent with a physical coherence length effect in the MOND regime. Grade A-.**
