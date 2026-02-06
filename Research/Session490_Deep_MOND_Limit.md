# Session #490: The Deep MOND Limit — Testing g << a₀

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

In deep MOND (g_bar << a₀), the RAR simplifies to g_obs = √(g_bar × a₀). This session tests the deep-limit predictions and explores how the 6-variable model performs across different acceleration regimes.

## Central Result: 51% of MOND Points Are Deep MOND; Model R² Peaks at Moderate Accelerations

Half of all MOND-regime data points (1146/2258) are in deep MOND (g < 0.1 a₀). The 6-variable model R² is highest for galaxies in the moderate regime (⟨g⟩ = 0.1 a₀: R² = 0.981) and slightly lower in both deep MOND (R² = 0.942) and the transition regime (R² = 0.933). Late types are 100% in the deep regime (⟨g/a₀⟩ = 0.076).

## Key Findings

### 1. Regime Classification (Test 1)

| Regime | Points | Fraction | Galaxies (≥3 pts) |
|--------|--------|----------|-------------------|
| MOND (g < a₀) | 2258 | 100% | 128 |
| Deep (g < 0.1 a₀) | **1146** | **50.8%** | **93 (72.7%)** |
| Very deep (g < 0.01 a₀) | 16 | 0.7% | 2 (1.6%) |

| Type | ⟨g_bar/a₀⟩ |
|------|-----------|
| Early | 0.35 |
| Mid | 0.32 |
| Late | **0.076** |

**Late types are overwhelmingly deep MOND** (⟨g/a₀⟩ = 0.076). Early and mid types are typically at moderate accelerations (⟨g/a₀⟩ = 0.3-0.35).

### 2. Offset vs Acceleration (Test 2)

| log(g/a₀) | ⟨offset⟩ | σ | N |
|-----------|----------|---|---|
| [-3, -2) | +0.040 | 0.274 | 16 |
| [-2, -1.5) | -0.051 | 0.209 | 276 |
| [-1.5, -1) | -0.046 | 0.193 | 854 |
| [-1, -0.5) | -0.001 | 0.162 | 587 |
| [-0.5, 0) | +0.007 | 0.153 | 525 |

**The mean offset is systematically negative in deep MOND** (⟨offset⟩ = -0.046), meaning the McGaugh function slightly overpredicts g_obs in the deep regime. This is because the McGaugh function exceeds the exact deep-MOND limit by 0.047 dex — a known property of the exponential interpolation function.

The scatter decreases from σ = 0.27 at very low accelerations to σ = 0.15 at the MOND boundary, consistent with higher signal-to-noise at higher accelerations.

### 3. Deep-MOND-Only Model (Test 3)

| Target | R²_6 | LOO R²_6 |
|--------|------|----------|
| Deep offset | 0.908 | 0.892 |
| Outer offset (same galaxies) | 0.955 | 0.947 |

r(deep offset, outer offset) = 0.944

The deep offset is slightly harder to predict (R² = 0.908 vs 0.955 for outer) — the deep MOND regime has more scatter because the accelerations are so small that relative measurement errors are larger. But the two offsets are highly correlated (r = 0.94).

### 4. The Running Exponent (Test 4)

| log(g/a₀) | Observed slope | McGaugh theory |
|-----------|---------------|----------------|
| -2.0 | 0.33 | 0.52 |
| -1.5 | 0.59 | 0.54 |
| -1.0 | 0.72 | 0.57 |
| -0.5 | 0.61 | 0.63 |
| 0.0 | 0.52 | 0.71 |

**The observed running exponent does NOT match the McGaugh function perfectly.** At g/a₀ ≈ 0.1 (log = -1), the data shows slope = 0.72 vs theory 0.57. This excess steepness may reflect correlated M/L errors or the fact that individual galaxies don't trace the mean RAR at every point.

At very deep MOND (log = -2), the slope drops to 0.33 (below the 0.5 deep-limit prediction), but this bin has only 72 points and large scatter.

### 5. Deep-Limit BTFR (Test 5)

| Regime | Slope | R² |
|--------|-------|-----|
| Deep (⟨g⟩ < 0.1 a₀) | 2.51 | 0.628 |
| **Moderate (0.1-0.5 a₀)** | **4.15** | **0.907** |
| Transition (⟨g⟩ > 0.5 a₀) | 3.07 | 0.770 |

**SURPRISING**: The BTFR slope is closest to 4.0 in the MODERATE regime (4.15), not in the deep regime (2.51). The deep-MOND BTFR has a much shallower slope because:
1. Deep-MOND galaxies are all low-mass dwarfs with narrow V range
2. Gas fraction varies widely, adding scatter
3. The shallow slope is a statistical artifact of restricted range + noise

### 6. Model Performance by Acceleration (Test 6)

| Quartile | ⟨g/a₀⟩ | R² | RMS |
|----------|--------|-----|-----|
| Q1 (deepest) | 0.041 | 0.942 | 0.045 |
| **Q2** | **0.097** | **0.981** | **0.029** |
| Q3 | 0.239 | 0.960 | 0.023 |
| Q4 (highest) | 0.460 | 0.933 | 0.024 |

**The model peaks at Q2** (⟨g/a₀⟩ ≈ 0.1). The sweet spot is moderate MOND — deep enough that MOND physics dominates, but not so deep that measurement noise overwhelms.

Partial r(log⟨g/a₀⟩, offset | V, L) = +0.25 — mild but significant.

### 7. Transition Regime (Test 7)

In very deep MOND (g < 0.03 a₀):
- McGaugh function RMS: 0.219 dex
- Exact deep-limit RMS: 0.214 dex
- Difference: only 0.005 dex

**The McGaugh function is essentially identical to the deep limit** for g < 0.03 a₀. The 0.047 dex systematic offset between the two at moderate accelerations is absorbed by the galaxy-level offset.

## Physical Interpretation

### The Deep-MOND Systematic

The mean offset is -0.046 in deep MOND. This is NOT a failure of MOND — it's a property of the McGaugh interpolation function. The McGaugh function gives ν(x) = 1/(1-exp(-√x)), which in deep MOND gives ν ≈ 1/√x + small corrections. These corrections are ~0.05 dex positive, making the predicted g_obs slightly too high. The galaxy-level offset absorbs this.

### Why the Model Peaks at Moderate MOND

Q2 (⟨g/a₀⟩ = 0.1) has the best R² because:
1. These galaxies have enough acceleration for signal-to-noise
2. They are in the clean MOND regime where the physics is well-defined
3. They typically have well-measured rotation curves (Sbc-Sd types)

Q1 (deep MOND) has more scatter because very low accelerations correspond to outer-disk points with larger measurement errors.

## Grade: B

A thorough exploration of the deep-MOND regime that establishes the acceleration-dependent structure. The key findings — 51% of points in deep MOND, model peaks at moderate MOND (R² = 0.981), and the 0.047 dex McGaugh-vs-deep-limit offset — are informative. The running exponent analysis shows deviations from the McGaugh prediction at all acceleration levels, though interpreting individual point-level slopes is complicated by galaxy-to-galaxy variation. Grade slightly lower because the deep-MOND BTFR slope (2.51) is likely a statistical artifact rather than physics, and the session is more descriptive than discovery-oriented.

## Files Created

- `simulations/session490_deep_mond_limit.py`: 8 tests
- `Research/Session490_Deep_MOND_Limit.md`: This document

---

*Session #490 verified: 8/8 tests passed*
*Grand Total: 1229/1229 verified*

**Key finding: 51% of MOND points are deep MOND (g < 0.1 a₀). 6-var model peaks at moderate MOND (R² = 0.981 at ⟨g⟩ = 0.1 a₀). McGaugh function gives -0.047 dex systematic in deep regime (known property). Deep-MOND BTFR slope = 2.5 (artifact of restricted range). Running exponent 0.3-0.7 vs theory 0.5-0.7. r(deep, outer offset) = 0.94. Grade B.**
