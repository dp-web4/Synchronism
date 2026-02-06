# Session #381: Environment vs Structure Disentanglement

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The Gas Fraction Control Arc (Sessions #376-378) established that morphology → RAR scatter is real (p = 5×10⁻⁶). But the fundamental ambiguity remained: is the effect driven by **environment** (Synchronism prediction: sparse environment → lower N_corr → higher γ → more scatter) or by **intrinsic structure** (late-type galaxies have messier kinematics)?

This session directly attacks that ambiguity using internal SPARC proxies and rotation curve shape analysis.

## Key Result: STRUCTURE FAVORED (5-1 evidence score)

The evidence strongly favors the **structural interpretation**: late-type galaxies have more RAR scatter because they have rougher rotation curves, not because they're in sparser environments.

## Detailed Findings

### 1. Composite Environment Proxy

Constructed an isolation index from four observables:
- Gas fraction (high = more isolated, gas not stripped)
- HI-to-stellar mass ratio (high = more isolated)
- Surface brightness (low = more isolated, LSBs preferentially in voids)
- Vflat (low = less massive, less likely in dense groups)

**Critical result**: The isolation proxy has r = +0.841 with Hubble type, meaning it's essentially a morphology proxy, not an independent environment measure.

| Partial Correlation | r | p | Favors |
|---|---|---|---|
| r(isolation, scatter \| type) | -0.035 | 0.65 | Structure |
| r(type, scatter \| isolation) | +0.155 | 0.04 | Structure |

After controlling for type, isolation has ZERO predictive power for scatter. But type retains marginal significance after controlling for isolation.

### 2. Within-Type Scatter Analysis

**THE critical test**: If environment matters, "isolated-like" late-types should have more scatter than "group-like" late-types.

| Type Bin | N | r(isolation, scatter) | Isolated/Dense ratio |
|---|---|---|---|
| Early (T≤3) | 28 | +0.152 (n.s.) | 1.170 |
| Mid (T=4-6) | 50 | -0.023 (n.s.) | 0.958 |
| Late-spiral (T=7-8) | 26 | +0.204 (n.s.) | 1.104 |
| Irregular (T≥9) | 67 | -0.044 (n.s.) | 0.970 |

**After removing type means**: r(isolation, residual scatter) = +0.003 (p = 0.97)

Within late types (T≥7, N=93): r = -0.001, isolated/dense ratio = 0.990

**Verdict**: Isolation predicts NOTHING within any morphological type. This is strong evidence against the environment hypothesis using these proxies.

### 3. HI Mass Ratio

With corrected MHI loading (was missing in initial run):

- Overall: r(log HI/stellar, scatter) = +0.166 (p = 0.028)
- **After controlling type**: r = -0.028 (p = 0.72) → Signal VANISHES
- **Within late types**: r = -0.018 (p = 0.86)
- **Within early types**: r = +0.072 (p = 0.63)
- Type-matched Mann-Whitney: z = +0.56, p = 0.58

The HI ratio-scatter correlation is entirely driven by the type-HI confound (r(HI ratio, type) = +0.796).

### 4. Surface Brightness Within Types

- Overall: r(log SB, scatter) = -0.240 (p = 0.001)
- **After controlling type**: r = -0.092 (p = 0.23) → Weakens substantially
- **After controlling roughness**: r = -0.021 (p = 0.78) → VANISHES completely

Within-type SB-scatter correlations are weak and non-significant in all bins. Only late types show a marginal LSB/HSB ratio of 1.189.

### 5. Gas Richness Extremes (Late Types Only)

Among 93 late-type galaxies:
- Gas-poor: σ = 0.121
- Gas-mid: σ = 0.123
- Gas-rich: σ = 0.109
- Ratio (rich/poor): **0.900** - gas-rich late types have LESS scatter!
- Mann-Whitney: p = 0.51 (non-significant)

This is the **opposite direction** from the environment prediction. If gas-rich means more isolated, and isolation should increase scatter (via γ), we'd expect gas-rich to have MORE scatter.

### 6. Rotation Curve Shape Analysis (KEY FINDING)

| Correlation | r | p | R² |
|---|---|---|---|
| Roughness → scatter | +0.711 | <10⁻¹⁰ | **50.6%** |
| \|Mean offset\| → scatter | +0.393 | <10⁻⁶ | 15.4% |
| Outer/inner ratio → scatter | -0.232 | 0.002 | 5.4% |

**Roughness explains half the scatter variance!** This is the single strongest predictor found in any session.

**After controlling roughness**:
- r(type, scatter | roughness) = +0.039 (p = 0.61) → Type is NON-SIGNIFICANT

**After controlling type**:
- r(roughness, scatter | type) = +0.692 (p < 10⁻¹⁰) → Roughness remains very strong

**Interpretation**: The type → scatter relationship is almost entirely mediated by roughness. Late-type galaxies have more scatter because they have rougher rotation curves (1.575x rougher than early types), not because of their environment.

### 7. Residual Pattern Analysis

- Offset fraction (systematic shift): 0.509 (51%)
- Fluctuation fraction (point-to-point): 0.491 (49%)

Late types have slightly more systematic offsets (55.1% vs 49.2%), marginally significant (p = 0.036).

Lag-1 autocorrelation: high for both types (early: 0.632, late: 0.590), difference non-significant (p = 0.27). Both types show smooth, correlated residuals - this is expected from systematic RAR misfit.

## Evidence Summary

| Test | r | p | Favors |
|---|---|---|---|
| Isolation → scatter \| type | -0.035 | 0.65 | **Structure** |
| Type → scatter \| isolation | +0.155 | 0.04 | **Structure** |
| Roughness → scatter | +0.711 | <10⁻¹⁰ | **Structure** |
| Type → scatter \| roughness | +0.039 | 0.61 | **Structure** |
| \|Mean offset\| → scatter | +0.393 | <10⁻⁶ | Environment |
| Within-type isolation | +0.003 | 0.97 | **Structure** |

**Score: Structure 5, Environment 1, Neutral 0**

## Implications for Synchronism

### The Hard Truth

The NP2 signal (morphology → RAR scatter) that we've carefully validated through Sessions #376-378 appears to be primarily a **structural** effect, not an environmental one. The rotation curve roughness of late-type galaxies - driven by non-circular motions, warps, bars, and kinematic irregularities - accounts for essentially all of the type → scatter correlation.

### What This Does NOT Mean

1. It does NOT disprove Synchronism's γ theory
2. It does NOT mean environment has no effect on the RAR
3. It does NOT mean NP2 is wrong (the scatter difference is real)
4. It means the **interpretation** needs updating: the observed scatter difference is not evidence for environment-dependent γ

### What This DOES Mean

1. NP2 should be **downgraded from "STRONGLY SUPPORTED" to "OBSERVED BUT NOT EVIDENCE FOR ENVIRONMENT"**
2. The critical test for environment-dependent γ remains **untested** - it requires actual environment data, not morphological proxies
3. The roughness finding is important: it identifies a specific structural mechanism for the scatter

### Remaining Hope for Environment Effect

The |mean offset| correlation (r = 0.393, the only environment-supporting evidence) suggests that systematic RAR offsets (not noise) contribute 15% of scatter variance. These offsets COULD reflect different effective γ values. But this is speculative without environment data.

### Key Caveat

The isolation index constructed from SPARC observables has r = 0.841 with Hubble type. This means it's NOT an independent environment measure - it's essentially morphology by another name. **We cannot distinguish environment from structure using internal SPARC data alone.** External environment catalogs (Karachentsev UNGC, Kourkchi-Tully groups, Chae et al. external field estimates) are essential.

## Lessons Learned

1. **Internal proxies cannot separate correlated effects**: When isolation correlates r = 0.84 with morphology, partial correlation is powerless to distinguish them
2. **Roughness is the hidden variable**: Rotation curve point-to-point variation explains 50% of scatter variance - this is a major finding
3. **Gas richness goes the WRONG way**: Gas-rich late types have LESS scatter (ratio 0.900), opposite to the environment prediction
4. **Be honest about negative results**: This session challenges the Synchronism interpretation of NP2, and that's valuable
5. **The type → scatter effect is MEDIATED by roughness**: r(type, scatter | roughness) = 0.039 means type has no independent effect once roughness is controlled

## Files Created

- `simulations/session381_environment_disentangle.py`: 8 tests
- `Research/Session381_Environment_Disentangle.md`: This document

---

*Session #381 verified: 8/8 tests passed*
*Grand Total: 495/495 verified*

**Key finding: The morphology → RAR scatter relationship is primarily STRUCTURAL, not environmental. Rotation curve roughness (point-to-point variation) explains 50.6% of scatter variance and completely mediates the type effect (r(type, scatter | roughness) = 0.039, n.s.). Within-type isolation proxies show zero signal (r = 0.003). Gas-rich late types have LESS scatter (ratio 0.90), opposite to the environment prediction. NP2 remains observationally valid but should be reinterpreted: the scatter difference reflects intrinsic kinematic complexity, not evidence for environment-dependent γ. External environment data (Karachentsev UNGC, Chae et al. external fields) is required for a definitive test.**
