# Session #463: The 5-Variable Model With Improved Interpolation

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session 461 found an improved interpolation function (α=0.458, a₀=1.276). Does the 5-variable model benefit from using it?

## Central Result: The Standard RAR Is Actually Better for the 5-Variable Model

| Metric | Standard (α=0.5, a₀=1.2) | Improved (α=0.458, a₀=1.276) |
|--------|--------------------------|-------------------------------|
| R² | **0.872** | 0.858 |
| RMS | **0.0556** | 0.0579 |
| LOO | **0.0589** | 0.0612 |
| Point-level RMS | **0.136** | 0.139 |

The improved interpolation function makes the model slightly **worse**. This is because the 5-variable model's coefficients were effectively absorbing the interpolation function error through the V and V×c_V terms.

## Key Findings

1. **Offsets are highly correlated**: r(standard, improved) = 0.995. The improved RAR shifts all offsets by ~0.024 dex, preserving the galaxy-to-galaxy structure.

2. **Variance budget shifts**: V explains less variance with the improved RAR (12.5% vs 17.8%), while L explains more (49.3% vs 44.4%). The c_V contribution is nearly identical (13.5% vs 13.1%).

3. **Residuals are nearly identical**: r(standard residual, improved residual) = 0.992. The model's explanatory power comes from galaxy-level physics, not interpolation function correction.

4. **Most galaxies unchanged**: 109/128 galaxies have |Δresid| < 0.01 dex. Only 5 improved, 14 worsened.

5. **Coefficients are stable**: logV (+2.77 → +2.77), logL (-0.49 → -0.50), c_V (+2.29 → +2.32), V×c_V (-0.92 → -0.93). Only f_gas shifts notably (-0.18 → -0.14).

## Physical Interpretation

The 5-variable model and the standard RAR form a **co-adapted system**. The logV coefficient (+2.77) partially compensates for the interpolation function's bias in the deep MOND regime. When we use a "better" interpolation function, this compensation becomes over-correction, slightly worsening the model.

This confirms that the 5-variable model is robust: it captures galaxy-level physics (M/L variation, phantom DM, gas fraction) regardless of the interpolation function choice. The point-level and galaxy-level analyses are complementary but largely independent.

## Grade: B

A clean comparison that produces the expected null result. The 5-variable model is robust to interpolation function changes because it operates at the galaxy level, where interpolation function errors average out. Not graded higher because the result, while confirming robustness, doesn't advance understanding.

## Files Created

- `simulations/session463_improved_rar_model.py`: 8 tests
- `Research/Session463_Improved_RAR_Model.md`: This document

---

*Session #463 verified: 8/8 tests passed*
*Grand Total: 1037/1037 verified*

**Key finding: The improved RAR (α=0.458) slightly worsens the 5-variable model (R²: 0.872→0.858) because the model's coefficients already compensate for standard RAR imperfections. Offsets are r=0.995 correlated — the galaxy-to-galaxy structure is preserved. The 5-variable model is robust to interpolation function choice. Grade B.**
