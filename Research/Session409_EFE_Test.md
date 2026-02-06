# Session #409: MOND External Field Effect (EFE) Test

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Tests whether the MOND External Field Effect can explain the R_eff → RAR offset correlation. In MOND, an external gravitational field suppresses the internal MOND boost, which could create size-dependent offsets if larger galaxies are systematically more isolated.

## Central Result: EFE Cannot Explain the R_eff Effect

**Mediation by ALL environment indicators: 2.3%**

| Environment control | r(R_eff, offset | V, env) | Mediation |
|--------------------|-----------------------|-----------|
| Distance | -0.707 | 4.1% |
| g_ext/g_int | -0.707 | 4.1% |
| g_ext (neighbors) | -0.731 | 0.9% |
| All combined | -0.721 | **2.3%** |

### EFE Sensitivity Test

| Subsample | N | r(R_eff, offset | V) |
|-----------|---|---------------------|
| Low g_int (most EFE-susceptible) | 30 | -0.736 |
| High g_int (least susceptible) | 31 | -0.736 |

**Identical** — the effect does not depend on EFE susceptibility.

## Why EFE Fails to Explain the Effect

1. **Too small**: EFE produces ~0.1 dex effects; our offset range is ~0.3 dex
2. **No mediation**: Environment indicators explain only 2.3% of R_eff signal
3. **No sensitivity**: Low and high g_int galaxies show identical R_eff correlations
4. **g_ext/g_int too small**: Mean g_ext/g_int ≈ 0.002 — EFE is negligible for these galaxies

## Grade: B+

Clean test but limited by the lack of proper environment data. The SPARC sample consists primarily of field galaxies where EFE is inherently weak, so a definitive EFE test would require galaxies in denser environments.

## Files Created

- `simulations/session409_efe_test.py`: 8 tests
- `Research/Session409_EFE_Test.md`: This document

---

*Session #409 verified: 8/8 tests passed*
*Grand Total: 677/677 verified*

**Key finding: MOND EFE mediates only 2.3% of the R_eff → RAR offset correlation. Low-g_int and high-g_int galaxies show identical R_eff effects (both r=-0.74). EFE is too weak for SPARC field galaxies (mean g_ext/g_int = 0.002) and too small in magnitude (~0.1 dex vs ~0.3 dex). Grade B+.**
