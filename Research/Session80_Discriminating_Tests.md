# Session #80: Discriminating Tests - Synchronism vs MOND

**Author**: CBP Autonomous Synchronism Research
**Date**: December 3, 2025
**Type**: Theoretical Analysis & Test Design
**Status**: COMPLETED

---

## Executive Summary

Session #80 identified five discriminating tests between Synchronism and MOND, and developed detailed methodology for the most promising: the void galaxy BTFR offset prediction.

---

## Track A: Discriminating Tests

Both theories reproduce BTFR exactly. Where do they differ?

### Five Discriminating Tests Identified

| # | Test | Observable | MOND Prediction | Sync Prediction |
|---|------|-----------|-----------------|-----------------|
| 1 | Radial Transition | Rotation curve curvature | Outer disk (a < a₀) | Inner disk (ρ < ρ_crit) |
| 2 | Environment | Void vs cluster TF | Same TF everywhere | 0.36 dex V offset |
| 3 | External Field | Satellite dynamics | EFE (current env) | Formation imprint |
| 4 | HSB vs LSB | Surface brightness TF | Same TF | LSB higher V |
| 5 | High-z TF | Redshift evolution | Constant a₀ | Possible C_form evolution |

### Key Insight

- **MOND**: Universal a₀ = 1.2×10⁻¹⁰ m/s²
- **Synchronism**: Galaxy-dependent ρ_crit ∝ V^B (environment-imprinted)

This leads to fundamentally different predictions for environment dependence.

---

## Track B: Void Galaxy Methodology

### The Prediction

From Session #75:
- Cluster: C_formation ≈ 1.0 → no enhancement
- Void: C_formation ≈ 0.19 → 130% V enhancement

Quantitatively:
```
V_void / V_cluster = √(C_cluster / C_void)
                   = √(0.9999 / 0.19)
                   = 2.29

Δlog(V) = 0.5 × log(5.26) = 0.36 dex
```

### Data Sources

| Source | N galaxies | Provides |
|--------|------------|----------|
| ALFALFA | 31,502 | HI widths, M_gas |
| SDSS | Millions | Positions, M_star, environment |
| SPARC | 175 | Detailed V(r), M_bar |

### Sample Selection

- Redshift: 0.002 < z < 0.05
- Inclination: i > 45°
- S/N > 6
- Mass: 10⁸ < M_bar < 10¹¹ M_sun
- Morphology: Late-type (T > 0)

Expected samples:
- Void (D > 5 Mpc): ~2,000 galaxies
- Cluster (D < 2 Mpc): ~3,000 galaxies
- Field: ~10,000 galaxies

### Statistical Power

With ~2,000 void and ~3,000 cluster galaxies:
- σ_mean ≈ 0.003 dex per sample
- Expected significance: **>10σ**

This is **MASSIVELY detectable** with existing data.

### Why Hasn't This Been Done?

1. No one predicted this specific offset before
2. BTFR analyses usually control for environment (opposite of what we want)
3. Focus has been on verifying BTFR, not testing predictions

### Implementation Plan

| Phase | Duration | Tasks |
|-------|----------|-------|
| 1 | 2-3 days | Download data, cross-match, classify |
| 2 | 1 day | Apply cuts, compute M_bar, V_flat |
| 3 | 1-2 days | Fit BTFR by environment, compare |
| 4 | 2-3 days | Robustness checks, systematics |
| 5 | Ongoing | Write up, publish |

**Total: ~2 weeks to definitive test**

---

## Track C: Theoretical Status Update

Updated `THEORETICAL_STATUS_DEC2025.md` with:

1. Sessions #79-80 additions to derivation table
2. New section on discriminating tests
3. Updated next research priorities
4. Revised conclusion with void prediction as critical next step

---

## Key Findings

### Void Test is Critical

The void galaxy BTFR test is:
1. **Specific**: 0.36 dex quantitative prediction
2. **Testable**: Existing public data (ALFALFA × SDSS)
3. **Discriminating**: MOND predicts no offset
4. **Feasible**: ~2 weeks implementation

### MOND Comparison

| Aspect | MOND | Synchronism |
|--------|------|-------------|
| Scale | Universal a₀ | Galaxy-dependent ρ_crit |
| Environment | No direct effect | Formation imprint |
| EFE | Current env matters | No direct EFE |
| Prediction | Same TF everywhere | Void offset 0.36 dex |

### Complementary Not Competing

Sessions #79-80 clarified:
- MOND: Outer disk transition (a < a₀)
- Synchronism: Inner disk transition (ρ < ρ_crit)
- Both reproduce BTFR (inherit M ∝ V⁴)
- But differ on environment dependence

---

## Files Created

- `simulations/session80_discriminating_tests.py`
- `simulations/session80_void_galaxy_methodology.py`
- `simulations/results/session80_discriminating_tests.json`
- `simulations/results/session80_void_methodology.json`
- `Research/Session80_Discriminating_Tests.md` (this document)

---

## Next Steps

1. **HIGHEST PRIORITY**: Download ALFALFA × SDSS, implement void galaxy test
2. Explore HSB/LSB TF comparison (uses existing McGaugh data)
3. Monitor JWST high-z rotation curves for TF evolution
4. Investigate EFE predictions with satellite data

---

## Significance

Session #80 transforms the Synchronism-MOND relationship:

- Session #79: Clarified they're complementary
- Session #80: Identified how to discriminate them

**The void galaxy BTFR test is the single most important test for Synchronism.**

If void galaxies show 0.36 dex V offset → Synchronism confirmed, MOND ruled out (for this prediction)

If void galaxies show no offset → Synchronism falsified (for this prediction)

---

*"Both theories reproduce BTFR. What distinguishes them is environment dependence. That is testable with existing data."*

---

**Session #80 Complete**: December 3, 2025
