# Session #459: Joint M/L and a₀ Fit — Breaking the Degeneracy

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session 458 found that SPARC's point-level best-fit a₀ ≈ 1.04 × 10⁻¹⁰ at fixed M/L = 0.5. But M/L and a₀ are degenerate: higher M/L → higher g_bar → higher preferred a₀. This session explores the joint (M/L, a₀) parameter space to determine whether the a₀ = cH₀/(2π) preference is genuine or an artifact of the M/L assumption.

## Central Result: Strong M/L-a₀ Degeneracy — Cannot Distinguish Theories

The data constrain a **line** in (M/L, a₀) space, not a point:

**M/L_disk ≈ -0.315 × a₀(×10⁻¹⁰) + 0.861**

Both (M/L=0.50, a₀=1.04) and (M/L=0.48, a₀=1.20) are nearly equally good fits. The RMS difference is only 0.8%.

## Key Findings

### 1. 2D Grid: Point-Level RMS (Test 1)

| M/L_disk | Best a₀ (×10⁻¹⁰) | RMS (dex) |
|----------|-------------------|-----------|
| 0.35 | 1.50 | 0.183 |
| 0.40 | 1.30 | 0.180 |
| 0.45 | 1.15 | 0.178 |
| **0.50** | **1.05** | **0.177** |
| **0.55** | **0.95** | **0.177** |
| 0.60 | 0.85 | 0.177 |
| 0.65 | 0.80 | 0.177 |

**Global minimum: M/L = 0.55, a₀ = 0.95 × 10⁻¹⁰.** But the RMS varies by only 0.001 dex across the M/L = 0.45-0.65 range — the minimum is extremely shallow.

### 2. Galaxy-Level Offset Scatter (Test 2)

Galaxy-level scatter is minimized at M/L = 0.70, a₀ = 1.50 × 10⁻¹⁰ — higher than the point-level optimum. Galaxy-level and point-level analyses weight different regimes, leading to different preferred parameters.

### 3. Profile Likelihood (Test 7)

Marginalizing over M/L at each a₀:

| a₀ (×10⁻¹⁰) | Best M/L | RMS | Note |
|-------------|----------|------|------|
| 0.70 | 0.675 | 0.1779 | |
| 0.79 | 0.625 | 0.1769 | |
| **0.89** | **0.575** | **0.1766** | **Profile minimum** |
| 0.98 | 0.550 | 0.1767 | |
| 1.04 | ~0.525 | 0.1772 | Planck cH₀/(2π) |
| 1.17 | 0.475 | 0.1779 | |
| 1.20 | ~0.475 | 0.1780 | MOND |

**Profile likelihood minimum: a₀ = 0.89 × 10⁻¹⁰ at M/L = 0.575.** This is lower than both MOND and cH₀/(2π). However, the profile is extremely flat — the difference between the minimum and a₀ = 1.2 is only 0.8%.

### 4. Gas-Dominated Galaxies (Test 4)

For f_gas > 0.5 galaxies (N=28), M/L should matter least:
- Best-fit: M/L = 0.75, a₀ = 0.70 × 10⁻¹⁰
- RMS varies 25% with M/L — **still M/L-dependent**, contradicting the expectation
- This surprising M/L sensitivity for gas-dominated galaxies may reflect: (a) gas-dominated galaxies still have stellar disks contributing g_bar, or (b) the gas fraction is computed at the flat region, not everywhere along the RC

### 5. Three-Parameter Fit (Test 5)

Best 3D: M/L_disk = 0.60, M/L_bul = 0.80, a₀ = 0.80 × 10⁻¹⁰
- Improvement over default: only **0.8%** in RMS
- The default (M/L_disk=0.5, M/L_bul=0.7, a₀=1.2) is nearly optimal

### 6. Subsample Dependence (Test 6)

**With default M/L (0.5, 0.7):**

| Subsample | Best a₀ (×10⁻¹⁰) |
|-----------|-------------------|
| Late (T≥7) | 0.89 |
| Early (T<7) | 1.19 |
| Gap: 0.30 | |

**With best-fit M/L (0.6, 0.8):**

| Subsample | Best a₀ (×10⁻¹⁰) |
|-----------|-------------------|
| Late (T≥7) | 0.77 |
| Early (T<7) | 0.89 |
| Gap: 0.12 | |

The subsample gap is reduced from 0.30 to 0.12 — the joint fit helps but doesn't eliminate the dependence. Late types consistently prefer lower a₀ than early types.

### 7. The M/L-a₀ Degeneracy Line (Test 7)

```
M/L_disk ≈ -0.315 × a₀(×10⁻¹⁰) + 0.861
```

This line passes through:
- (a₀=1.20, M/L=0.48) — near standard MOND
- (a₀=1.04, M/L=0.53) — near cH₀/(2π) with Planck H₀
- (a₀=0.89, M/L=0.58) — profile likelihood minimum

All three points give essentially the same RMS (within 0.8%).

## Physical Interpretation

### The Degeneracy Is Fundamental

M/L and a₀ enter the RAR through g_bar: g_bar ∝ M/L × L / R ∝ M/L. Higher M/L shifts all g_bar values to the right on the RAR plot, which is equivalent to shifting a₀ to the right. The two parameters are nearly perfectly degenerate at the ~1% level.

### Why Session 458's Result Was Misleading

Session 458 found a₀_best = 1.04 × 10⁻¹⁰ at fixed M/L = 0.5. This is correct but incomplete: the result says "M/L = 0.5 prefers a₀ ≈ 1.04" — it does NOT say "SPARC proves a₀ = cH₀/(2π)." The agreement with cH₀/(2π) is a consequence of the specific M/L choice.

### What Would Break the Degeneracy

1. **Independent M/L from stellar population synthesis**: If SPS models firmly predict M/L = 0.5 ± 0.05 for disk galaxies, then a₀ = 1.04 ± 0.10 would follow
2. **Purely gas-dominated galaxies**: DDO-class galaxies with f_gas > 0.9 at all radii would eliminate M/L dependence — but SPARC has very few such galaxies
3. **Joint fit with informative priors**: Bayesian analysis with M/L priors from photometry could break the degeneracy

## Grade: B+

A necessary and technically sound analysis that correctly identifies the M/L-a₀ degeneracy as the main limitation of the a₀ constraint. The profile likelihood analysis is well-executed. However, the result is somewhat disappointing — we cannot distinguish a₀ = 1.2 from a₀ = 1.04 with SPARC data alone when M/L is free. The gas-dominated analysis is limited by the small subsample size and residual M/L sensitivity.

## Files Created

- `simulations/session459_joint_ml_a0.py`: 8 tests
- `Research/Session459_Joint_ML_a0.md`: This document

---

*Session #459 verified: 8/8 tests passed*
*Grand Total: 1013/1013 verified*

**Key finding: M/L and a₀ are strongly degenerate in SPARC: M/L ≈ -0.315 × a₀(×10⁻¹⁰) + 0.861. Profile likelihood minimum at a₀ = 0.89, M/L = 0.575, but the surface is flat (RMS varies only 0.8% from a₀=0.7 to 1.5). The Session 458 result (best-fit a₀ ≈ 1.04 at M/L=0.5) reflects the M/L choice, not a genuine cH₀/(2π) preference. Subsample gap reduced from 0.30 to 0.12 with joint fit. Cannot break the degeneracy without independent M/L or purely gas-dominated galaxies. Grade B+.**
