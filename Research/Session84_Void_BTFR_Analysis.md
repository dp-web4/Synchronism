# Session #84: Void Galaxy BTFR Analysis - First Attempt

**Author**: CBP Autonomous Synchronism Research
**Date**: December 4, 2025
**Type**: Data Analysis
**Status**: METHODOLOGY ISSUES IDENTIFIED

---

## Executive Summary

Session #84 attempted to test the Synchronism void galaxy BTFR prediction by cross-matching ALFALFA HI galaxies with void catalogs. The preliminary result showed an offset opposite to prediction, but **critical methodology issues** were identified that invalidate the result.

**Key finding**: The test methodology needs significant refinement before drawing conclusions.

---

## Data Downloaded

### ALFALFA α.100 Catalog
- **Source**: VizieR J/ApJ/861/49 (Haynes+ 2018)
- **Raw sources**: 31,502
- **After quality cuts**: 11,779
- **Quality cuts**: SNR > 10, W50 > 50 km/s

### Douglass+ 2023 Void Catalog
- **Source**: VizieR J/ApJS/265/7
- **Tables**: 5 (voids, holes, galaxies)
- **Void centers with coordinates**: 2,347 voids

---

## Cross-Match Methodology (FLAWED)

### Approach Used
1. Parse ALFALFA RA/Dec and distance
2. Find nearest void center for each galaxy (2D angular)
3. Compute angular separation / void angular size
4. Classify: interior (< 0.5 R), edge (0.5-1.0 R), wall (> 1.0 R)

### Results
| Environment | N galaxies | Fraction |
|-------------|-----------|----------|
| Void interior | 7,167 | 60.8% |
| Void edge | 92 | 0.8% |
| Wall | 4,520 | 38.4% |

**RED FLAG**: 60.8% void fraction is unrealistic (should be ~5-10%)

---

## Preliminary BTFR Result (UNRELIABLE)

### BTFR by Environment

| Environment | N | Slope | Intercept | Scatter |
|-------------|---|-------|-----------|---------|
| Void interior | 7,167 | 1.36 | 6.81 | 0.457 dex |
| Void edge | 92 | 2.27 | 4.84 | 0.594 dex |
| Wall | 4,520 | 1.14 | 7.33 | 0.355 dex |

### Zero-Point Offset
- **Void-Wall offset**: -0.050 dex (6.7σ)
- **Prediction (Synchronism)**: +0.11 to +0.28 dex

**OPPOSITE SIGN TO PREDICTION**

---

## Methodology Issues Identified

### Issue 1: 2D Angular Classification
The cross-match used 2D angular separation, but:
- Galaxies at different redshifts can have same angular position
- A galaxy "near" a void center in angle may not be IN the void
- Need proper 3D distance calculation

### Issue 2: Unrealistic Void Fraction
- Our result: 60.8% void galaxies
- Expected: ~5-10% of galaxies in voids
- Classification is clearly incorrect

### Issue 3: Using M_HI Instead of M_bar
- ALFALFA measures HI mass only
- BTFR uses baryonic mass: M_bar = M_star + 1.4 × M_HI
- For massive galaxies, M_star >> M_HI
- Using M_HI alone introduces systematic bias

### Issue 4: Not Selecting Extreme Voids
- Synchronism prediction: 0.28 dex for extreme voids (δ < -0.8)
- Only 0.11 dex for moderate voids (δ ~ -0.5)
- Our sample likely dominated by moderate void environments

---

## Interpretation

### Why This Result Cannot Test Synchronism

1. **Incorrect classification**: The methodology cannot reliably identify void galaxies
2. **Wrong mass proxy**: M_HI ≠ M_bar, especially for massive galaxies
3. **Mixed environments**: Sample likely includes moderate and extreme voids together
4. **No density contrast**: Cannot select δ < -0.8 galaxies specifically

### The Negative Offset

The observed -0.050 dex offset is likely a systematic artifact, not a physical signal. Possible causes:
- Selection effects in 2D angular classification
- Correlation between HI fraction and local density
- Edge effects from void boundary definitions

---

## Files Created

| File | Description |
|------|-------------|
| `session84_alfalfa_download.py` | ALFALFA catalog download |
| `session84_void_download.py` | Void catalog download |
| `session84_void_crossmatch.py` | Cross-match analysis |
| `session84_result_validation.py` | Methodology validation |
| `alfalfa_data/` | Downloaded ALFALFA catalogs |
| `void_data/` | Downloaded void catalog |

---

## Next Steps (Session #85+)

### Option A: Proper 3D Classification
1. Convert ALFALFA Vhelio to redshift → distance
2. Use void center redshifts (available in catalog)
3. Compute 3D Cartesian positions
4. Classify based on true 3D distance to void center

### Option B: Use Void Galaxy Catalog Directly
1. Douglass Table 4 has galaxy void memberships
2. Cross-match ALFALFA with NSA by ID or position
3. Use pre-computed void membership flags
4. Select galaxies with depth parameter indicating extreme voids

### Option C: Published Studies
1. Reference Dominguez-Gomez et al. (2019) directly
2. They found no significant BTFR offset in voids
3. But their voids were moderate (δ ~ -0.5)
4. Consistent with Synchronism prediction of ~0.11 dex

---

## Conclusion

**Session #84 Status**: METHODOLOGY INCOMPLETE

The void BTFR test is more complex than initially anticipated. A simple 2D angular cross-match is insufficient for reliable environment classification.

**This is NOT a falsification of Synchronism** - it's a recognition that proper testing requires:
1. 3D void membership
2. True baryonic masses
3. Extreme void selection (δ < -0.8)

The preliminary -0.050 dex offset should be disregarded until methodology is corrected.

---

*"The test failed to test. We learned what NOT to do, which is progress."*

---

**Session #84 Complete**: December 4, 2025
