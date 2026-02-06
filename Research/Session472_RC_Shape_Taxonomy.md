# Session #472: Rotation Curve Shape Taxonomy

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Rotation curves come in distinctive shapes: rising, flat, declining, peaked. This session quantifies these shapes, classifies them, and tests whether shape parameters add predictive power beyond c_V.

## Central Result: 66% of RCs Are Flat, and Shape = f(mass)

Two-thirds of rotation curves are flat (|outer slope| < 5%). Rising RCs are found in dwarfs (⟨logV⟩ = 1.93, ⟨T⟩ = 7.3) and declining RCs in massive galaxies (⟨logV⟩ = 2.22, ⟨T⟩ = 4.7). No shape parameter adds to the 5-variable model — c_V already captures the relevant shape information.

## Key Findings

### 1. Shape Classification (Test 2)

| Class | N | % | ⟨logV⟩ | ⟨T⟩ | ⟨c_V⟩ |
|-------|---|---|--------|-----|-------|
| Rising | 22 | 17% | 1.93 | 7.3 | 0.64 |
| Flat | 85 | 66% | 2.04 | 6.5 | 0.85 |
| Declining | 21 | 16% | 2.22 | 4.7 | 1.00 |
| Peaked | 14 | 11% | — | — | — |

The classification follows the mass sequence perfectly: rising RCs in slow-rotating dwarfs, flat RCs in intermediate spirals, declining RCs in fast-rotating early-type disks.

### 2. The Universal Rotation Curve (Test 5)

Mean V/V_flat at r/R_eff by mass class:

| r/R_eff | Dwarfs (30-70) | Low mass (70-130) | Massive (130-200) | Giant (200-350) |
|---------|---------------|-------------------|-------------------|-----------------|
| 0.5 | 0.42 | 0.57 | 0.73 | 0.95 |
| 1.0 | 0.63 | 0.77 | 0.96 | 1.03 |
| 2.0 | 0.84 | 0.93 | 1.03 | 1.06 |
| 3.0 | 0.89 | 0.97 | 1.03 | 1.04 |
| 5.0 | 0.91 | 0.99 | 1.04 | 1.02 |

**Giant galaxies reach V_flat by 0.5 R_eff; dwarfs don't reach it until 5 R_eff.** This is the mass-dependent rotation curve shape: massive galaxies have concentrated mass distributions that produce quickly-rising, peaked-and-declining RCs. Dwarfs have diffuse mass distributions producing slowly-rising RCs that only flatten in the deep MOND regime.

### 3. Self-Similarity (Test 6)

| r/R_eff | CV(V_obs) | CV(V/V_flat) | Reduction |
|---------|-----------|-------------|-----------|
| 1 | 0.68 | 0.25 | 63% |
| 3 | 0.54 | 0.11 | 79% |
| 5 | 0.47 | 0.09 | 80% |

Scaling by V_flat reduces scatter by 63-80%, confirming partial self-similarity. The remaining scatter (CV = 0.09-0.25) reflects genuine shape diversity — the mass-dependent inner structure that c_V quantifies.

### 4. Shape Residuals (Test 7)

Shape residuals at different radii are strongly correlated (r = 0.60-0.90), indicating a single "shape mode." A galaxy that's above-average V/V_flat at r = R_eff is also above-average at r = 5 R_eff. This single mode is c_V — the velocity concentration already in the 5-variable model.

No shape parameter improves the 5-var model by more than ΔR² = 0.006 (outer slope), and none improve LOO significantly.

## Physical Interpretation

The rotation curve shape encodes two things:
1. **Mass concentration** (bulge-to-disk ratio → c_V)
2. **Mass** (total dynamical mass → V_flat)

Both are already in the 5-variable model. Additional shape parameters (outer slope, peak ratio, roughness) are collinear with (V, c_V) and add no independent information.

## Grade: B

A competent characterization of RC shapes that confirms c_V captures the essential shape information. The universal rotation curve by mass bin (Test 5) and self-similarity analysis (Test 6) are informative. However, the main conclusion is negative (no shape parameter improves the model), and the analysis doesn't reveal any surprising new physics.

## Files Created

- `simulations/session472_rc_shape_taxonomy.py`: 8 tests
- `Research/Session472_RC_Shape_Taxonomy.md`: This document

---

*Session #472 verified: 8/8 tests passed*
*Grand Total: 1101/1101 verified*

**Key finding: 66% of RCs are flat, 17% rising (dwarfs), 16% declining (massive). Giant galaxies reach V_flat by 0.5 R_eff; dwarfs not until 5 R_eff. V_flat scaling reduces scatter by 63-80%. Shape residuals at different radii are correlated (single mode = c_V). No shape parameter adds to the 5-var model. Grade B.**
