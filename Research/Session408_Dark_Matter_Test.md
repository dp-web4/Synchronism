# Session #408: Dark Matter Halo Test — Can NFW Scatter Produce r=-0.74?

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Tests whether standard ΛCDM dark matter halos can explain the R_eff → RAR offset correlation. In ΛCDM, galaxies with different sizes at fixed rotation speed may have different dark matter halo concentrations, which would create size-dependent RAR offsets.

## Central Result: ΛCDM Direction Correct, Strength Probably Insufficient

| Aspect | ΛCDM prediction | Observed | Match? |
|--------|----------------|----------|--------|
| Direction | r < 0 (less conc. → larger, lower g_obs) | r = -0.74 | YES |
| Magnitude | σ = 0.05-0.15 dex | σ = 0.14 dex | MARGINAL |
| Correlation strength | r ≈ -0.1 to -0.3 | r = -0.74 | **NO (2.5-7× too strong)** |
| DM mediation | ~100% | ~15-18% | **NO** |

## Key Findings

### 1. DM Fraction Does NOT Fully Mediate R_eff Effect

| Model | r | Mediation |
|-------|---|-----------|
| r(R_eff, offset \| V) | -0.74 | baseline |
| r(R_eff, offset \| V, DM_frac) | -0.60 | 18% |
| r(R_eff, offset \| V, V_ratio) | -0.47 | 36% |
| r(R_eff, offset \| V, DM_frac, V_ratio) | -0.63 | 15% |

**CAVEAT**: These DM indicators (DM_frac, V_ratio) involve V_obs, which also enters the offset through g_obs. The mediation test is partially circular. True independent DM indicators would require simulation data.

### 2. The Correlation Strength Problem

Standard ΛCDM with abundance matching predicts:
- R_eff ∝ R_vir × λ (halo spin parameter)
- Halo concentration c is weakly anti-correlated with λ: r(λ, c) ≈ -0.1 to -0.3
- Expected r(R_eff, c | V) ≈ -0.1 to -0.3
- Observed r(R_eff, offset | V) = -0.74

The observed correlation is **2.5-7× stronger** than expected from standard abundance matching. This is a significant tension with ΛCDM predictions.

### 3. What Would Resolve This?

If hydrodynamic simulations (EAGLE, IllustrisTNG, FIRE) produce r(R_eff, RAR_offset | V_flat) ≈ -0.7 for late-type galaxies, then ΛCDM could explain our result through complex baryonic physics. If they produce r ≈ -0.2, then ΛCDM fails.

**This is the single most important comparison that could falsify or validate our interpretation.**

## Grade: B+

The session provides a thorough ΛCDM comparison but is limited by our inability to access simulation data. The key finding (correlation too strong for standard abundance matching) is meaningful but depends on theoretical expectations that could be revised.

## Files Created

- `simulations/session408_dark_matter_test.py`: 8 tests
- `Research/Session408_Dark_Matter_Test.md`: This document

---

*Session #408 verified: 8/8 tests passed*
*Grand Total: 669/669 verified*

**Key finding: ΛCDM can explain the DIRECTION of the R_eff effect (less concentrated halos → larger R_eff → lower g_obs) and marginally the MAGNITUDE (σ=0.14 dex vs expected 0.05-0.15). But the CORRELATION STRENGTH (r=-0.74) appears 2.5-7× stronger than expected from standard abundance matching (r≈-0.1 to -0.3). DM fraction mediates only 18% of the effect. Definitive test requires comparison with hydrodynamic simulations. Grade B+.**
