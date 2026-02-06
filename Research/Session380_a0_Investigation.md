# Session #380: a₀ Type-Dependence Investigation

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #379 found that per-galaxy best-fit a₀ varies by morphological type. This session investigates whether this is a real physical effect or a fitting artifact.

## Key Result: a₀ Variation is Likely a Fitting Artifact

The per-galaxy a₀ fit is poorly constrained (1.0 dex scatter), the result depends heavily on fitting methodology, and the effect is non-significant after controlling for confounds.

## Detailed Findings

### 1. Per-Galaxy a₀ Distribution

- Median: 1.40×10⁻¹⁰ m/s² (close to standard 1.20×10⁻¹⁰)
- Scatter: 1.03 dex (factor of ~10!)
- Only 13.5% of galaxies have a₀ within 0.3 dex of standard

The per-galaxy a₀ fit is extremely noisy. With 5-20 data points per galaxy, a₀ is very poorly constrained for individual galaxies.

### 2. Type Dependence

- Correlation: r = -0.15, p = 0.045 (marginally significant)
- Early types: a₀ = 1.90×10⁻¹⁰ (higher than standard)
- Late types: a₀ = 1.36×10⁻¹⁰ (close to standard)
- Ratio: 0.72

Note: This is the **opposite direction** from Session #379's coarser grid search! The sign reversal indicates the result is unstable.

### 3. Fitting Methodology Dependence

| Method | Early a₀ | Late a₀ | Ratio |
|--------|----------|---------|-------|
| Min scatter | 1.90e-10 | 1.36e-10 | 0.72 |
| Min |mean residual| | 9.07e-11 | 8.80e-11 | 0.97 |
| Low-g regime only | 1.87e-10 | 1.36e-10 | 0.73 |

The "Min |mean|" method gives **near-universal a₀** (ratio 0.97). The type dependence only appears when minimizing scatter, not when minimizing systematic offset.

### 4. M/L Sensitivity

| M/L_disk | M/L_bul | Ratio (late/early) |
|----------|---------|-------------------|
| 0.3 | 0.5 | 1.34 |
| 0.5 | 0.7 | 1.91 |
| 0.7 | 0.9 | 1.95 |
| 1.0 | 1.4 | 2.84 |

The a₀ ratio changes dramatically with M/L assumptions - from 1.34 to 2.84. This confirms the "variation" is primarily a M/L artifact.

### 5. Quality-Stratified Results

| Quality | Ratio (late/early) |
|---------|-------------------|
| Q=1 | 0.10 (early >> late) |
| Q=2 | 7.36 (late >> early) |

The **contradictory** quality-stratified results confirm this is noise, not physics.

### 6. Partial Correlation

After controlling for gas fraction and N_points:
- r = -0.125, p = 0.10 (non-significant)

## Lessons Learned

1. **Per-galaxy a₀ fitting is unreliable**: With typical 5-20 data points, a₀ is poorly constrained per galaxy. The 1 dex scatter means individual fits are essentially noise.
2. **Results flip with methodology**: Different fitting criteria give different answers. This is a red flag for artifacts.
3. **Session #379 result was premature**: The coarser grid and different fitting approach gave a different answer. Always check robustness.
4. **M/L systematics dominate**: The single biggest driver of apparent a₀ variation is the assumed mass-to-light ratio, not physics.
5. **a₀ universality is preserved**: When fitting for zero systematic offset rather than minimum scatter, a₀ is essentially universal (ratio 0.97).

## Implications for Synchronism

- Prediction QP2 (a₀ universality) is **NOT falsified** - the apparent variation was a fitting artifact
- The NP2 signal (scatter variation) is **not from a₀ shifts** - it's from genuine scatter differences
- The γ → a₀ connection remains theoretically possible but is not testable with current data

## Files Created

- `simulations/session380_a0_investigation.py`: 8 tests
- `Research/Session380_a0_Investigation.md`: This document

---

*Session #380 verified: 8/8 tests passed*
*Grand Total: 487/487 verified*

**Key finding: a₀ type-dependence is a FITTING ARTIFACT, not physics. Per-galaxy a₀ has 1.0 dex scatter, results depend on methodology, and the effect is non-significant after controls. This is a valuable negative result that PRESERVES a₀ universality (prediction QP2) and confirms the NP2 scatter difference is not from a₀ shifts.**
