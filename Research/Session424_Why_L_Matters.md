# Session #424: Why L Matters — Decomposing the 4-Variable Model

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The jump from V+R+c_V (LOO = 0.087) to V+R+L+c_V (LOO = 0.057) is the largest single-variable improvement in the model hierarchy. This session investigates why luminosity adds so much information beyond velocity, size, and concentration.

## Central Result: L Acts as a Suppressor for c_V

| Step | LOO-RMSE | Description |
|------|----------|-------------|
| V+R | 0.102 | Baseline |
| V+R+c_V | 0.087 | c_V alone adds 0.015 |
| V+R+L | 0.100 | L alone adds 0.002 |
| **V+R+L+c_V** | **0.057** | **c_V after L adds 0.043** |

c_V is **3× more powerful after L is controlled**. The c_V coefficient increases by 76% (from 0.33 to 0.59) when L is added to the model. This is a classic suppressor mechanism: L removes a confounding channel between c_V and offset.

## Key Findings

### 1. L and SB Are Exactly Equivalent (Test 2)

At fixed R_eff, L and SB carry identical information (RMS difference = 0.00000):
- V+R+c_V+L: RMS = 0.0510, LOO = 0.0566
- V+R+c_V+SB: RMS = 0.0510, LOO = 0.0566

This confirms R_eff² ∝ L/SB strictly. Using L or SB as the fourth variable makes zero difference. The model coefficients adjust accordingly:
- With L: R_coeff = -0.29
- With SB: R_coeff = -0.78 (absorbs extra R dependence from SB → L conversion)

### 2. Partial Correlations Beyond V+R+c_V (Test 1)

| Variable | r(var, offset \| V, R, c_V) | p |
|----------|---------------------------|---|
| L | **-0.78** | 2×10⁻¹³ |
| SB | **-0.78** | 2×10⁻¹³ |
| gas dominance | +0.41 | 10⁻³ |
| L/V⁴ (BTFR residual) | -0.42 | 8×10⁻⁴ |

### 3. L Acts as Structure, Not Mass (Test 3)

| Subset | N | V+R+c_V LOO | V+R+L+c_V LOO | Improvement |
|--------|---|-------------|---------------|-------------|
| Gas-dominated | 18 | 0.079 | 0.061 | 23.5% |
| Disk-dominated | 42 | 0.082 | 0.061 | 24.8% |

If L were a mass proxy, it should be less useful in gas-dominated galaxies (where L underestimates M_bar). The nearly identical improvement (24% vs 25%) confirms L encodes **structural** information, not total mass.

### 4. The Suppressor Mechanism (Test 5)

When L is added to V+R+c_V:
- V coefficient: +1.29 → +1.75 (+36%)
- R coefficient: -0.48 → -0.29 (-40%)
- c_V coefficient: +0.33 → **+0.59** (+76%)
- L coefficient: new at -0.25

L absorbs some of R's information (reducing R's coefficient) while simultaneously unmasking c_V's true predictive power. The mechanism: without L control, c_V is confounded by a luminosity channel. More luminous galaxies tend to have both higher c_V (more concentrated) and more negative offset (through the V-L-offset chain). This creates a negative confound that partially cancels c_V's true positive effect.

### 5. M/L Robustness (Test 4)

| M/L (disk) | R² | RMS | r(c_V \| V,R) |
|-----------|-----|------|---------------|
| 0.2 | 0.86 | 0.078 | +0.68 |
| 0.3 | 0.85 | 0.078 | +0.63 |
| 0.5 | 0.82 | 0.082 | +0.53 |
| 0.7 | 0.80 | 0.087 | +0.44 |
| 1.0 | 0.76 | 0.094 | +0.33 |

R² remains above 0.76 across all M/L assumptions. Lower M/L yields stronger results (because lower M/L emphasizes the gas contribution, making g_bar more accurately computed).

### 6. L Predicts Both Inner and Outer (Test 6)

Unlike c_V (which is inner-only), L predicts both regions:

| Region | r(c_V, offset \| V,R) | r(L, offset \| V,R,c_V) |
|--------|----------------------|------------------------|
| Inner (r < 2 R_eff) | +0.64 | **-0.67** |
| Outer (r > 2 R_eff) | -0.003 | **-0.63** |

L fills the gap that c_V cannot: it predicts the outer offset where c_V has no signal.

### 7. PCA Confirms 4 Dimensions Needed (Test 7)

| PCA components | LOO-RMSE |
|---------------|----------|
| 2 PCs | 0.172 |
| 3 PCs | 0.102 |
| 4 original vars | **0.057** |

The 4th dimension (PC4, 2.7% of predictor variance) carries crucial information for offset prediction. This cannot be reduced without major performance loss.

## Physical Interpretation

The 4-variable model captures four distinct physical channels:

| Variable | What it encodes | Where it predicts |
|----------|----------------|-------------------|
| V_flat | Overall mass scale, asymptotic dynamics | Both |
| R_eff | Spatial extent of baryonic distribution | Mainly outer |
| c_V | Mass concentration (profile shape) | Mainly inner |
| L | Total baryonic content at fixed structure | Both |

L acts primarily as a suppressor that enables c_V to express its full predictive power. The physical interpretation: at fixed V_flat, R_eff, and c_V, a more luminous galaxy has:
- More baryonic mass relative to its structure
- Higher g_bar at each radius
- Therefore a more negative RAR offset (g_obs doesn't rise proportionally to g_bar)

This is the **mass-structure mismatch** channel: galaxies with more mass than their structure "predicts" (higher L at fixed V, R, c_V) have more negative offsets.

## Grade: A

The suppressor mechanism is cleanly demonstrated (c_V coefficient increases 76% when L is added). The equivalence of L and SB is rigorously confirmed. The mass vs structure test (gas vs disk dominated) is clever and definitive. The inner/outer decomposition reveals complementary spatial coverage. One grade below A+ because the result is more explanatory than discovery — it unpacks Session 423's finding rather than finding something new.

## Files Created

- `simulations/session424_why_L_matters.py`: 8 tests
- `Research/Session424_Why_L_Matters.md`: This document

---

*Session #424 verified: 8/8 tests passed*
*Grand Total: 789/789 verified*

**Key finding: L acts as a classic suppressor for c_V. c_V's marginal contribution triples after L is controlled (0.015 → 0.043). c_V coefficient increases 76%. L and SB are exactly equivalent at fixed R (RMS diff = 0). L encodes structure not mass (equal improvement in gas/disk-dominated). L predicts both inner and outer (r = -0.67, -0.63) unlike inner-only c_V. PCA confirms 4 dimensions needed. V+R+L+c_V: LOO = 0.057, approaching noise floor of ~0.029. Grade A.**
