# Session #172: High-Precision Subset and 3D Density Analysis

**Date**: December 23, 2025
**Status**: Critical methodological insights
**Arc**: Real Data Application

---

## Executive Summary

Session #172 followed the Session #171 recommendation to analyze the high-precision CF4 subset (SNIa, TRGB, Cepheids) and investigated an unexpected finding: 3D density classification gives **opposite results** to angular density classification.

### Key Findings

1. **High-precision subset**: 958 galaxies with 5-8% distance errors
2. **3D vs Angular discrepancy**: Correlation only r = 0.054
3. **3D density in distance bins**: Persistently shows Synchronism direction (ratio 2.5-5.0)
4. **Angular density in distance bins**: Persistently shows opposite (ratio 0.7-1.1)

---

## 1. High-Precision Subset Analysis

### Sample Characteristics

| Method | N galaxies | Mean σ_DM (mag) | Rel. distance error |
|--------|-----------|-----------------|---------------------|
| SNe Ia | 920 | 0.165 | 7.6% |
| TRGB | 32 | 0.106 | 4.9% |
| Cepheids | 52 | 0.083 | 3.8% |
| **Combined** | **958** | - | **~6%** |

For comparison:
- Tully-Fisher: 20.6% error
- Fundamental Plane: 24.4% error

### Velocity-Environment Test

Using angular overdensity:

| Environment | N | Mean |v| (km/s) |
|-------------|---|---------------------|
| Void (Q1) | 223 | 517 ± 36 |
| Overdense (Q4) | 202 | 609 ± 42 |
| **Ratio** | - | **0.849** |

Result: **Opposite to Synchronism** (predicts ratio > 1)

### Distance-Binned Analysis

| Distance | Ratio | Direction |
|----------|-------|-----------|
| 10-50 Mpc | 0.67 | Opposite |
| 30-70 Mpc | 0.86 | Opposite |
| **50-90 Mpc** | **1.05** | **Synchronism** |
| 100-140 Mpc | 0.87 | Opposite |

Only the 50-90 Mpc bin shows Synchronism direction.

---

## 2. Group Peculiar Velocities (table4.dat)

Analyzed 26,077 groups with published peculiar velocities.

### Angular Density Classification

| Environment | Mean |Vpec| (km/s) |
|-------------|----------------------|
| Void (Q1) | 580 ± 8 |
| Overdense (Q4) | 1047 ± 14 |
| **Ratio** | **0.553** |

Result: **Opposite to Synchronism**

---

## 3. Critical Finding: 3D vs Angular Density

### The Discrepancy

Using supergalactic coordinates (SGX, SGY, SGZ):

| Metric | Void/Overdense Ratio |
|--------|---------------------|
| Angular density | 0.55 |
| **3D density** | **4.23** |

The 3D density gives a dramatically **opposite** result!

### Investigation

```
Correlation (angular δ vs 3D δ): r = 0.054
```

The two density metrics are essentially **uncorrelated**.

### Key Insight: Distance Bias

| Classification | Mean distance (Mpc) |
|----------------|---------------------|
| Angular void | 109 |
| Angular overdense | 177 |
| **3D void** | **214** |
| **3D overdense** | **70** |

The 3D "voids" are preferentially at **larger distances**!

This would create a false Synchronism signal because:
- |v| correlates with distance (r = 0.40)
- Larger distance → larger velocity errors → larger apparent |v|

### BUT: Distance-Controlled Tests

When controlling for distance in narrow bins:

**3D Density (distance-matched):**
| Distance | Ratio | Direction |
|----------|-------|-----------|
| 30-60 Mpc | 2.49 | Synchronism |
| 50-80 Mpc | 3.52 | Synchronism |
| 80-110 Mpc | 4.92 | Synchronism |

**Angular Density (distance-matched):**
| Distance | Ratio | Direction |
|----------|-------|-----------|
| 30-60 Mpc | 0.85 | Opposite |
| 50-80 Mpc | 0.69 | Opposite |
| 80-110 Mpc | 1.09 | Synchronism |
| 120-150 Mpc | 0.95 | Opposite |

---

## 4. Interpretation

### The Puzzle

Two density metrics give **opposite** conclusions:

1. **Angular density** (2D projection): Void galaxies have LOWER |v|
2. **3D density** (supergalactic coords): Void galaxies have HIGHER |v|

Both persist when controlling for distance.

### Possible Explanations

1. **Angular density misclassifies true voids**
   - Projects along line of sight
   - Mixes foreground/background
   - May misidentify fingers-of-god as voids

2. **3D density captures real structure**
   - Uses physical coordinates
   - But may be affected by redshift-space distortions

3. **Neither captures "true" voids**
   - Need proper void catalogs from independent surveys
   - Void centers, not just low-density regions

4. **True signal is small**
   - Expected ~15-35% enhancement
   - Easily swamped by selection effects

---

## 5. Recommendations

### Immediate

1. **Obtain published void catalogs**
   - Pan et al. 2012 (SDSS DR7 voids)
   - Sutter et al. 2012 (SDSS void catalog)
   - Cross-match CF4 with void-centric coordinates

2. **Forward model selection effects**
   - Simulate Malmquist bias
   - Predict expected velocity-density correlation from errors alone

3. **Alternative tests**
   - Galaxy infall into clusters (independent of peculiar velocity errors)
   - Void outflow profiles
   - Bulk flow analysis

### Longer-term

4. **Wait for DESI peculiar velocities**
   - Smaller errors
   - Larger samples
   - Better void identification

---

## 6. Files Created

| File | Description |
|------|-------------|
| `session172_high_precision_subset.py` | High-precision subset analysis |
| `session172b_group_velocities.py` | Group peculiar velocity analysis |
| `session172c_3d_density_investigation.py` | 3D vs angular density investigation |
| `session172_high_precision_analysis.png` | Visualization |
| `session172b_group_velocities.png` | Group velocity figures |
| `session172c_density_investigation.png` | Density comparison figures |

---

## 7. Key Takeaway

**The Synchronism signal in CF4 data remains ambiguous due to methodology sensitivity.**

- Angular density: suggests opposite to prediction
- 3D density: suggests Synchronism direction
- Neither is definitive due to selection effects

The path forward requires:
1. Independent void catalogs
2. Forward modeling of biases
3. Complementary observational tests

---

*Session completed: December 23, 2025*
