# Session #174: Forward Modeling Selection Effects

**Date**: December 24, 2025
**Status**: Critical methodological advance
**Arc**: Real Data Application

---

## Executive Summary

Session #174 developed forward models to quantify selection effects in CF4 peculiar velocity data. The key finding is **surprising and important**:

### Main Result

**CF4 shows LESS velocity-environment correlation than pure selection effects predict.**

This means:
- After accounting for distance-dependent errors
- CF4 peculiar velocities are **LOWER** in low-density regions
- This is **OPPOSITE** to the Synchronism prediction

---

## 1. The Selection Effect Problem

### Distance-Dependent Errors

CF4 peculiar velocities are dominated by distance errors:

| Statistic | Value |
|-----------|-------|
| Mean σ_v | 1758 km/s |
| Mean |v_pec| | 626 km/s |
| Signal-to-noise | 0.43 |
| Correlation (σ_v, distance) | r = 0.82 |

**Key insight**: Velocity errors scale with distance, creating spurious correlations.

### Why Forward Modeling?

Previous sessions showed conflicting results depending on methodology:
- Session #172: 3D vs angular density gave opposite conclusions
- Session #173: Void catalog results were methodology-sensitive

Forward modeling answers: **What signal do we EXPECT from selection effects alone?**

---

## 2. Forward Model Design

### Null Hypothesis

1. **True velocities**: σ_v = 300 km/s (ΛCDM), NO environment dependence
2. **Distance errors**: Match CF4 distribution (σ_dm = 0.46 mag)
3. **Positions**: Same as CF4 (same environment classification)

### Test Procedure

1. Generate random true velocities (no correlation with environment)
2. Add realistic distance errors
3. Compute observed peculiar velocities
4. Measure velocity-environment ratio
5. Compare to CF4 observed ratio

---

## 3. Results by Density Metric

### 3D Density Classification

| Quantity | Null Model | CF4 Observed |
|----------|------------|--------------|
| Void mean distance | 163 Mpc | 154 Mpc |
| Overdense mean distance | 59 Mpc | 50 Mpc |
| **Distance ratio** | **2.76** | **3.08** |
| Void <|v|> | 2252 km/s | 842 km/s |
| Overdense <|v|> | 799 km/s | 336 km/s |
| **Velocity ratio** | **2.82 ± 0.06** | **2.50** |
| Z-score | - | **-11.07** |

**Finding**: CF4 shows LESS ratio than null (2.50 vs 2.82).

### Angular Density Classification (More Reliable)

| Quantity | Null Model | CF4 Observed |
|----------|------------|--------------|
| Void mean distance | 94 Mpc | 94 Mpc |
| Overdense mean distance | 117 Mpc | 117 Mpc |
| **Distance ratio** | **0.80** | **0.80** |
| Void <|v|> | - | 521 km/s |
| Overdense <|v|> | - | 689 km/s |
| **Velocity ratio** | **0.93 ± 0.02** | **0.756** |
| Z-score | - | **-10.56** |

**Finding**: CF4 shows significantly lower ratio than null (0.756 vs 0.93).

---

## 4. Interpretation

### Best-Fit True Enhancement

To match CF4 observations, we need:

| Density Metric | Best-Fit Enhancement |
|----------------|----------------------|
| 3D density | 0.70 (30% lower in voids) |
| Angular density | 0.60 (40% lower in voids) |

### What This Means

**CF4 galaxies in low-density regions have TRUE peculiar velocities that are 30-40% LOWER than those in high-density regions.**

This is the **OPPOSITE** of the Synchronism prediction (15-35% higher in voids).

### Possible Explanations

1. **Physical suppression**: Low-density regions genuinely have lower velocity dispersion
   - Consistent with fewer nearby mass sources
   - Standard ΛCDM prediction: voids have lower velocity dispersion

2. **Synchronism ruled out**: The coherence function prediction is wrong
   - G_eff enhancement in voids does not produce velocity enhancement

3. **CF4 bias corrections**: Published v_pec may already include corrections
   - Would reduce spurious distance-error correlations
   - Need to check CF4 methodology

4. **Wrong environment proxy**: Angular/3D density may not capture true voids
   - Need real void catalogs from independent surveys

---

## 5. Key Methodological Insights

### Distance Bias Comparison

| Metric | Void/Overdense Distance Ratio | Reliability |
|--------|-------------------------------|-------------|
| 3D density | 3.08 | Low (severe bias) |
| Angular density | 0.80 | Higher (minimal bias) |

Angular density is more reliable because it uses only sky positions, not distances.

### Forward Model Validation

The null model correctly predicts:
- |v| vs distance correlation (r ≈ 0.4)
- Expected velocity ratio from errors alone
- Distribution shape of observed velocities

---

## 6. Implications for Synchronism

### The Challenge

Synchronism predicts:
```
v_void/v_overdense ≈ 1.15-1.35
```

CF4 shows:
```
v_void/v_overdense ≈ 0.76 (angular density)
```

After correcting for selection effects:
```
True enhancement ≈ 0.60 (40% LOWER in voids)
```

### Possible Reconciliations

1. **Scale mismatch**: Synchronism effect may only appear at different scales
2. **Density definition**: Angular density may not capture coherence transition regions
3. **Velocity field vs individual**: CF4 uses processed velocities, not raw measurements
4. **Sample bias**: CF4 is dominated by TF/FP galaxies with specific selection

### What Would Support Synchronism

- Void/overdense ratio > 1.0 with angular density
- OR ratio > null prediction for any reliable environment metric
- Cluster dispersion test showed promising results (ratio = 3.28)

---

## 7. Files Created

| File | Description |
|------|-------------|
| `session174_forward_model.py` | Initial null hypothesis model |
| `session174b_null_model_refinement.py` | Refined with CF4 positions |
| `session174c_angular_density_test.py` | Angular density analysis |
| `session174_forward_model.png` | Visualization |
| `session174b_null_refinement.png` | Refined model figures |
| `session174c_angular_test.png` | Angular density figures |

---

## 8. Cumulative Progress (Sessions #169-174)

| Session | Key Finding |
|---------|-------------|
| #169 | Mock data: 8.4σ detection possible |
| #170 | Velocity-quartile: >> 10σ significance |
| #171 | Distance errors dominate |
| #172 | 3D vs angular discrepancy |
| #173 | Methodology-sensitive; cluster promising |
| **#174** | **Forward model: CF4 is OPPOSITE to Synchronism** |

### Overall Assessment

The forward modeling shows that after accounting for selection effects:
- CF4 does NOT support Synchronism peculiar velocity predictions
- Low-density regions have LOWER velocities, not higher
- The cluster dispersion test remains the most promising approach

---

## 9. Next Steps

### Immediate

1. **Check CF4 v_pec definition**: Are bias corrections already applied?
2. **Obtain high-quality velocity sample**: SN Ia only
3. **Use proper void catalogs**: Pan+2012, Sutter+2012

### Alternative Tests

4. **Cluster infall profiles**: Less affected by distance errors
5. **Bulk flow measurements**: Different systematic effects
6. **Void outflow profiles**: Direct test of gravity modification

### Theoretical

7. **Refine Synchronism prediction**: Account for velocity field smoothing
8. **Consider alternative signatures**: Beyond peculiar velocities

---

*Session completed: December 24, 2025*
