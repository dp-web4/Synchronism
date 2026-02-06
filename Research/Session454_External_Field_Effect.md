# Session #454: The MOND External Field Effect in SPARC Residuals

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

In MOND, the internal dynamics of a galaxy are affected by the external gravitational field (EFE). This session tests whether the ~3% residual scatter correlates with environment proxies.

## Central Result: No Evidence for EFE in SPARC Residuals

| Environment proxy | r(X, best residual) | ΔR² if added |
|-------------------|---------------------|-------------|
| log Distance | -0.021 | — |
| Outer RC slope | -0.059 | +0.0006 |
| V_decline (V_last/V_max) | +0.030 | +0.0003 |
| Neighbor count (±3 Mpc) | +0.046 | — |
| EFE proxy (nearby × declining) | -0.044 | +0.0003 |
| Inclination | -0.126 | +0.0024 |

**All correlations are weak** (|r| < 0.13). No proxy achieves significant ΔR². The F-test for adding outer_slope to the 5-variable model gives p = 0.46 (not significant).

## Key Findings

### 1. No Distance Dependence (Test 1)
- r(log D, best residual) = -0.021 — essentially zero
- Nearby galaxies (D < 5 Mpc): mean residual +0.018
- Distant galaxies (D > 40 Mpc): mean residual +0.001
- No systematic trend across distance bins

### 2. Outer RC Shape (Test 2)
- 22/128 galaxies have declining RCs (V_decline < 0.9)
- r(outer_slope, residual) = -0.059 (raw), -0.067 (partial)
- The outer RC slope is already partially captured by c_V (r = -0.39)

### 3. Density Proxy (Test 3)
- Isolated (0-5 neighbors): mean residual -0.001
- Dense (>15 neighbors): mean residual +0.001
- No density dependence in residuals

### 4. Deep MOND Regime (Test 5)
- 33 galaxies in deep MOND (f_MOND > 0.8, logV < 1.9)
- Even in this regime (where EFE should be strongest): r(log D, residual) = +0.10 (weak)
- Nearby deep MOND: residual -0.015 vs distant: -0.003 — not significant

### 5. RC Flatness as 6th Variable (Test 6)
- Adding outer_slope: ΔR² = +0.0006 (F = 0.56, p = 0.46)
- Adding V_decline: ΔR² = +0.0003
- **Neither improves the model**

## Physical Interpretation

The MOND EFE is a real theoretical prediction, but it cannot be detected in SPARC for several reasons:

1. **SPARC lacks sky coordinates**: We can only use distance as an environment proxy, which is very crude
2. **The sample is heterogeneous**: SPARC was not designed for environment studies
3. **The EFE is subtle**: Predicted effects are ~5-10% on the outer RC, which is below our noise floor
4. **c_V already captures some of the signal**: r(outer_slope, c_V) = -0.39

A proper EFE test requires: (a) a volume-complete survey, (b) known external field values, (c) galaxies with extended rotation curves well into the MOND regime.

## Implications for the ~3% Residual

The ~3% residual scatter is NOT explained by any observable in our dataset:
- Not distance/environment (this session)
- Not inclination, quality, or measurement noise (Session 449)
- Not any interaction or nonlinear term (Session 449, 451)

The residual appears to be **genuine Gaussian noise** from: stellar population diversity, non-equilibrium dynamics, or unmodeled physics below the detection threshold.

## Grade: B

A clean negative result that confirms the ~3% residual is unstructured. The EFE analysis was limited by SPARC's lack of sky coordinates. The result is methodologically sound but doesn't advance understanding beyond confirming the noise floor.

## Files Created

- `simulations/session454_external_field_effect.py`: 8 tests
- `Research/Session454_External_Field_Effect.md`: This document

---

*Session #454 verified: 8/8 tests passed*
*Grand Total: 981/981 verified*

**Key finding: No evidence for MOND External Field Effect in SPARC residuals. All environment proxies have |r| < 0.13 with the 5-variable model residual. RC flatness adds ΔR²=0.0006 (p=0.46). The ~3% residual scatter is unstructured — truly random with respect to all available observables. The 5-variable model is essentially saturated. Grade B.**
