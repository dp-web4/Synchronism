# Session #466: The Baryonic Tully-Fisher Relation — A Deep Dive

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The Baryonic Tully-Fisher Relation (BTFR: M_bar ∝ V⁴) is the tightest galaxy scaling relation known. In MOND, it follows directly from the deep-MOND limit. This session examines the BTFR in our SPARC sample and its connection to the RAR offset and 5-variable model.

## Central Result: BTFR Residuals Anti-Correlate With RAR Offset at r = -0.83

The BTFR residual (excess mass at fixed V) is strongly anti-correlated with the RAR offset (excess g_obs at fixed g_bar): r = -0.83 controlling for V. This means galaxies that are "too massive" for their velocity are also "too faint" gravitationally — consistent with M/L overestimation being the dominant source of scatter in both relations.

## Key Findings

### 1. BTFR Slope and Normalization (Test 1)

| Method | Slope | Intercept | RMS (dex) |
|--------|-------|-----------|-----------|
| OLS (M\|V) | 3.604 | 2.418 | 0.300 |
| Orthogonal | 3.995 | 1.615 | 0.315 |
| Inverse OLS | 4.022 | 1.559 | — |
| MOND prediction | 4.000 | 1.798 | — |

The OLS slope (3.60) is shallower than MOND predicts (4.00) due to Malmquist-type bias — low-mass galaxies with upward scatter in M are preferentially included. The **orthogonal regression slope (3.99) matches MOND perfectly**. The inverse OLS (4.02) also matches. This is a textbook demonstration of regression method dependence.

R²(BTFR) = 0.896, r(logV, logM) = 0.947.

### 2. BTFR Intrinsic Scatter (Test 2)

| Component | Value (dex) |
|-----------|-------------|
| Observed scatter | 0.300 |
| Velocity error contribution (4×0.04) | 0.160 |
| M/L uncertainty contribution | 0.150 |
| **Intrinsic scatter** | **0.205** |

The intrinsic BTFR scatter (0.21 dex) is larger than the RAR offset scatter (0.155 dex). This makes sense: the BTFR collapses all radial information into two numbers (M, V), while the RAR offset captures the full radial structure.

### 3. Gas-Rich vs Gas-Poor BTFR (Test 3)

| Subset | N | Slope | RMS (dex) |
|--------|---|-------|-----------|
| Gas-rich (f_gas > 0.5) | 28 | 3.636 | 0.243 |
| Gas-poor (f_gas ≤ 0.5) | 100 | 3.406 | 0.302 |

**Gas-rich galaxies have 20% less scatter (0.243 vs 0.302 dex)**, confirming that M/L uncertainty is a dominant scatter source — gas-rich galaxies have less M/L uncertainty because gas mass is directly measured. Both slopes are below 4.0 (OLS bias), but gas-rich is closer.

### 4. BTFR Residual vs RAR Offset (Test 4)

| Correlation | Value |
|-------------|-------|
| r(BTFR_resid, offset) | **-0.753** |
| r(BTFR_resid, offset \| V) | **-0.831** |
| r(BTFR_resid, offset \| V, L) | -0.576 |
| ΔR² adding BTFR_resid to 5-var | +0.004 |

The strong anti-correlation (-0.83 at fixed V) reveals the underlying physics: both the BTFR residual and RAR offset are driven by M/L. A galaxy with over-estimated M/L has too much inferred baryonic mass (positive BTFR residual) AND too little gravitational acceleration relative to that mass (negative RAR offset). They anti-correlate because they measure the same error from opposite sides.

**The BTFR residual adds nothing to the 5-variable model (ΔR² = 0.004)** — the model already captures this information through logV and logL.

### 5. BTFR Curvature (Test 5)

| Model | RMS | ΔBIC |
|-------|-----|------|
| Linear | 0.300 | 0 |
| Quadratic | 0.274 | **-18.4** |

**The BTFR has significant curvature (ΔBIC = -18).** The quadratic coefficient is c₂ = +1.64, meaning the BTFR steepens at high velocities. This is expected in MOND: the transition from deep MOND (slope 4) to Newtonian (slope undefined — no BTFR) creates curvature at intermediate accelerations.

The scatter also varies with velocity: lowest quartile (dwarfs) has σ = 0.43 dex, while upper quartiles have σ ≈ 0.20-0.23 dex.

### 6. a₀ From the BTFR Normalization (Test 6)

| Method | a₀ (×10⁻¹⁰ m/s²) |
|--------|-------------------|
| BTFR (slope fixed at 4) | **1.865 ± 0.121** |
| MOND standard | 1.200 |
| cH₀/(2π) | 1.042 |

**The BTFR-derived a₀ = 1.87 is 5.5σ from MOND and 6.8σ from cH₀/(2π).** This is much higher than from point-level RAR fitting (a₀ ≈ 1.04). The discrepancy arises because:

1. The BTFR uses M_bar = M_star + M_gas with M/L_disk = 0.5 (fixed assumption)
2. The gas mass estimate uses our simplified f_gas proxy (V²_gas / V²_bar at flat region)
3. The BTFR normalization is sensitive to the mean M/L of the sample

This confirms Session 459's finding: a₀ and M/L are degenerate. A higher M/L shifts the BTFR normalization up, requiring a higher a₀.

### 7. BTFR by Hubble Type (Test 7)

| Type | N | OLS Slope | RMS (dex) | ⟨f_gas⟩ |
|------|---|-----------|-----------|---------|
| S0-Sa | 12 | 2.87 | 0.217 | 0.18 |
| Sab-Sb | 26 | 2.81 | 0.183 | 0.11 |
| Sbc-Sc | 30 | **4.23** | **0.161** | 0.21 |
| Scd-Sm | 19 | 3.94 | 0.204 | 0.39 |
| Im-BCD | 41 | 2.22 | 0.339 | 0.50 |

**Sbc-Sc galaxies have the tightest BTFR (0.161 dex) and the slope closest to 4 (4.23).** This is the "Goldilocks" type: massive enough for well-measured rotation curves, but late enough that M/L variation is modest. Im-BCD galaxies have by far the highest scatter (0.339 dex) — their irregular dynamics and small masses make V_flat poorly defined.

With slope fixed at 4, the type dependence is dramatic: Im-BCD scatter is 0.46 dex (nearly 3× the Sbc-Sc value of 0.16 dex). Late types also show a systematic offset: Scd-Sm galaxies are 0.16 dex below the mean BTFR.

## Physical Interpretation

### The BTFR-RAR Connection

The BTFR and RAR are two projections of the same underlying gravitational physics:

- **BTFR**: Integrates all mass and velocity into two numbers (M_bar, V_flat)
- **RAR**: Preserves radial structure (g_bar(r), g_obs(r) at each radius)
- **5-variable model**: Adds structural information (c_V, f_gas) to bridge the two

The -0.83 correlation between BTFR residual and RAR offset shows they share a common driver: M/L estimation error. When M/L is overestimated, M_bar rises (positive BTFR residual) while the predicted g_RAR rises above g_obs (negative offset).

### Why a₀(BTFR) ≠ a₀(RAR)?

The BTFR gives a₀ = 1.87 while the point-level RAR gives a₀ = 1.04. This 80% discrepancy comes from how M_bar is computed:
- **Point-level**: g_bar uses the mass model at each radius
- **BTFR**: M_bar = M/L × L + M_gas, summing all mass at once

The difference reflects the radial M/L structure: the effective M/L varies with radius, and the global average M/L (used in BTFR) differs from the local M/L weights (used in RAR). This is another manifestation of the M/L-a₀ degeneracy.

### The BTFR Curvature

The significant quadratic term (ΔBIC = -18) means the BTFR is not a perfect power law. MOND predicts this: the exact BTFR is only a power law in the deep MOND limit (g << a₀). At intermediate accelerations (g ~ a₀), the interpolation function introduces curvature. The positive c₂ (+1.64) means the slope steepens at high V — massive galaxies follow a steeper BTFR than dwarfs.

## Grade: B+

A thorough characterization of the BTFR that reveals several interesting findings: the -0.83 anti-correlation with RAR offset (confirming M/L as the dominant scatter driver), significant curvature (ΔBIC = -18), the orthogonal slope matching MOND perfectly (3.99 vs 4.00), and a₀(BTFR) = 1.87 revealing the M/L-a₀ degeneracy from a new angle. The gas-rich subsample has 20% less scatter as expected. However, the BTFR adds nothing new to the 5-variable model (ΔR² = 0.004), confirming it contains a subset of the information already captured.

## Files Created

- `simulations/session466_btfr_deep_dive.py`: 8 tests
- `Research/Session466_BTFR_Deep_Dive.md`: This document

---

*Session #466 verified: 8/8 tests passed*
*Grand Total: 1061/1061 verified*

**Key finding: The BTFR orthogonal slope (3.995) matches MOND perfectly (4.000). BTFR residuals anti-correlate with RAR offset at r = -0.83 (controlling V), confirming M/L drives scatter in both relations. The BTFR has significant curvature (ΔBIC = -18). Gas-rich scatter is 20% lower. a₀(BTFR) = 1.87 ± 0.12 — far from MOND (1.2) or cH₀/(2π) (1.04) — revealing the M/L-a₀ degeneracy from the BTFR side. BTFR residuals add nothing to the 5-var model (ΔR² = 0.004). Grade B+.**
