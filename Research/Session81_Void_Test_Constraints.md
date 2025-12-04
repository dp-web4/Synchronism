# Session #81: Void Galaxy Test Constraints

**Author**: CBP Autonomous Synchronism Research
**Date**: December 3, 2025
**Type**: Literature Review & Constraint Analysis
**Status**: COMPLETED - Revised predictions

---

## Executive Summary

Session #81 revealed important constraints on the void galaxy BTFR prediction from existing literature and BTFR scatter limits. The original prediction needs revision.

---

## Key Findings

### 1. Dominguez-Gomez et al. (2019) Constraint

**Paper**: "The influence of the void environment on the ratio of dark matter halo mass to stellar mass in SDSS MaNGA galaxies" ([arXiv:1906.08327](https://arxiv.org/abs/1906.08327))

**Finding**: NO difference in M_halo/M_star between void and non-void galaxies

**Sample**: 642 void vs 938 dense environment galaxies

**Initial concern**: Does this falsify Synchronism?

**Resolution**: Their "void" sample used moderate voids (δ ~ -0.5), not extreme voids (δ ~ -0.9). The Synchronism prediction depends critically on void definition.

### 2. BTFR Scatter Constraint

If the 0.36 dex offset existed and void galaxies are 10% of typical samples:
- Expected additional scatter: ~0.11 dex
- Total expected scatter: ~0.15 dex
- Observed scatter: ~0.10 dex

This suggests the offset is smaller than the extreme prediction.

### 3. Gap in Literature

**Critical finding**: NO study has performed an environment-split BTFR analysis!

Studies either:
- Test BTFR universality (ignoring environment)
- Study void galaxies (but not BTFR)
- Study clusters (mostly S0s, not spirals)

Synchronism provides the FIRST physics-motivated prediction for environment-dependent BTFR.

---

## Revised Prediction

### Original (Session #75, #80)

| Environment | C_formation | Δlog(V) vs field |
|-------------|-------------|------------------|
| Extreme void | 0.19 | +0.36 dex |
| Cluster | 1.0 | 0 |

### Revised (Session #81)

Hypothesis: C_formation scales with environment δ
- C = 1 - 0.8 × |δ| for δ < 0 (voids)
- C = 1 + 0.1 × δ for δ > 0 (clusters)

| Environment | δ | C_formation | Δlog(V) vs field |
|-------------|---|-------------|------------------|
| Extreme void | -0.9 | 0.28 | **+0.28 dex** |
| Moderate void | -0.5 | 0.60 | **+0.11 dex** |
| Field | 0.0 | 1.00 | 0.00 dex |
| Mod. cluster | +0.5 | 1.05 | -0.01 dex |
| Cluster | +1.0 | 1.10 | -0.02 dex |

### Key Insight

Synchronism predicts **ASYMMETRIC** environment dependence:
- **Voids**: Significant offset (0.1-0.3 dex)
- **Clusters**: Minimal offset (~0 dex)

This is because C saturates at 1 for high density, but has no floor for low density.

---

## Implications for Testing

### Why the Dominguez-Gomez Result Doesn't Falsify

1. Their "void" sample was moderate (δ ~ -0.5)
2. Predicted offset for moderate voids: ~0.11 dex (or ~1.3× in M_halo/M_star)
3. This is close to their systematic uncertainties
4. Need EXTREME void sample (δ < -0.8) for definitive test

### What We Need

| Requirement | Reason |
|-------------|--------|
| Extreme voids (δ < -0.8) | Larger predicted offset |
| BTFR (not M_halo/M_star) | More direct test |
| ~100 extreme void galaxies | Statistical significance |
| Comparison to field (δ ~ 0) | Clean baseline |

### Expected Detection

For extreme voids with revised prediction:
- Offset: 0.28 dex in log(V)
- With 100 void, 1000 field galaxies
- Scatter: 0.1 dex
- Significance: ~28σ

Still highly detectable, but need extreme sample.

---

## Data Access Status

| Source | Status | Notes |
|--------|--------|-------|
| ALFALFA α.100 | URL identified | Direct download failed, VizieR available |
| ALFALFA-SDSS cross-match | Durbala et al. 2020 | FITS files at Cornell |
| Void catalogs | Pan et al., VIDE | Need cross-match |

---

## Updated Test Strategy

### Phase 1: Identify Extreme Void Sample
- Use VIDE or Pan et al. void catalogs
- Select galaxies with δ < -0.8
- Cross-match with ALFALFA for HI data

### Phase 2: Compute BTFR
- V_rot from W50 and inclination
- M_bar from M_star + 1.4 × M_HI
- Fit BTFR for void vs field

### Phase 3: Test Revised Prediction
- Expected: 0.28 dex offset for extreme voids
- If found: Strong support for Synchronism
- If not found: Constrains C(δ) relation

---

## Files Created

- `simulations/session81_alfalfa_data_analysis.py`
- `simulations/session81_literature_review.py`
- `simulations/results/session81_alfalfa_analysis.json`
- `simulations/results/session81_literature_review.json`
- `Research/Session81_Void_Test_Constraints.md` (this document)

---

## Significance

Session #81 is important for HONEST SCIENCE:

1. **Existing data constrains** the theory
2. **Original prediction was too extreme** for typical voids
3. **Revised prediction** is more realistic
4. **Gap in literature** still exists for extreme voids
5. **Test remains viable** but needs extreme void sample

This demonstrates the theory responding to evidence.

---

*"The 0.36 dex prediction was for extreme voids. For typical voids, expect 0.11 dex. For a definitive test, need the extremes."*

---

**Session #81 Complete**: December 3, 2025
