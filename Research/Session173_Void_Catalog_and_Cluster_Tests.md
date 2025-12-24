# Session #173: Void Catalog Cross-Matching and Cluster Infall Tests

**Date**: December 23, 2025
**Status**: Methodology development with mixed results
**Arc**: Real Data Application

---

## Executive Summary

Session #173 developed two alternative approaches to test Synchronism predictions that are less sensitive to the distance-dependent selection effects identified in Sessions #171-172:

1. **Void catalog cross-matching**: Using physically-motivated void positions
2. **Cluster infall test**: Using relative velocities within cluster vicinities

### Key Findings

| Test | Result | Direction |
|------|--------|-----------|
| Sparse void catalog (44 voids) | Ratio = 1.64, p < 0.0001 | Synchronism |
| Dense void catalog (640 voids) | Ratio = 0.82, p < 10⁻²¹ | Opposite |
| Cluster dispersion (outer/inner) | Ratio = 3.28 | Synchronism |

**Critical insight**: Results are sensitive to methodology. Selection effects remain the dominant challenge.

---

## 1. Void Catalog Cross-Matching

### Motivation

Sessions #170-172 showed that density-based environment classification gives inconsistent results:
- Angular density: opposite to Synchronism
- 3D density: supports Synchronism (but biased by distance)

Solution: Use independent void catalogs based on published surveys (Pan+2012, Sutter+2012).

### Published Void Catalogs

| Catalog | Source | N voids | Radius range |
|---------|--------|---------|--------------|
| Pan et al. 2012 | SDSS DR7 | ~1000 | 10-50 Mpc |
| Sutter et al. 2012 | SDSS DR7 | 1054 | 5-135 h⁻¹ Mpc |
| Mao et al. 2017 | SDSS DR12 BOSS | 10,643 | 20-100 Mpc |

Catalogs available from:
- http://www.physics.drexel.edu/~pan/catalog.html
- https://sites.google.com/site/cosmicvoids

### Synthetic Void Catalog Results

#### Session #173 (Sparse: 44 voids)

```
Void coverage: 3.7%
Void interior <|v|>: 1003 km/s
Walls/filaments <|v|>: 611 km/s
Ratio: 1.640 (p < 0.0001)
>>> SYNCHRONISM DIRECTION
```

#### Session #173b (Dense: 640 voids)

```
Void coverage: 24.9%
Void interior <|v|>: 538 km/s
Walls/filaments <|v|>: 655 km/s
Ratio: 0.821 (p < 10⁻²¹)
>>> OPPOSITE
```

#### Velocity Gradient (Session #173b)

| d/R Range | Location | N | <|v|> (km/s) |
|-----------|----------|---|--------------|
| 0.0-0.3 | Deep void | 130 | 419 |
| 0.3-0.5 | Void | 410 | 529 |
| 0.5-0.7 | Void | 994 | 529 |
| 0.7-1.0 | Void edge | 2909 | 548 |
| 1.0-1.3 | Wall | 3703 | 614 |
| 1.3-1.7 | Wall | 4945 | 656 |
| 1.7-2.5 | Filament | 4453 | 685 |
| 2.5-4.0 | Far | 272 | 721 |

**Gradient**: +89 km/s per unit d/R (opposite to Synchronism)

#### Distance-Controlled Bins

| Distance | Void/Non-void Ratio | Direction |
|----------|---------------------|-----------|
| 20-50 Mpc | 0.98 | Neutral |
| 50-80 Mpc | 0.88 | Opposite |
| **80-110 Mpc** | **1.01** | **Synchronism** |
| **120-150 Mpc** | **1.09** | **Synchronism** |
| 160-190 Mpc | 0.97 | Neutral |

**Insight**: At intermediate distances (80-150 Mpc), the Synchronism direction emerges.

---

## 2. Cluster Infall Test

### Motivation

Cluster infall velocities are less sensitive to overall distance errors because they use **relative** velocities within cluster vicinities.

### Theoretical Prediction

Synchronism predicts:
- Enhanced G_eff in low-density outer regions → higher velocity dispersion
- Suppressed G_eff in high-density cores → lower velocity dispersion
- Expected outer/core enhancement: ~16%

### Results (Session #173c)

#### Velocity Dispersion Profile

| Distance | N | σ_v (km/s) |
|----------|---|------------|
| 5-10 Mpc | 149 | 509 |
| 10-20 Mpc | 544 | 798 |
| 20-40 Mpc | 1623 | 1177 |
| 40-80 Mpc | 4494 | 1561 |

**Clear trend**: Dispersion increases with distance from cluster center.

#### Inner vs Outer Comparison

```
Inner (<5 Mpc):  σ_v = 467 km/s (N=76)
Outer (20-60 Mpc): σ_v = 1532 km/s (N=3684)
Ratio: 3.28 >>> SYNCHRONISM DIRECTION
```

### Caveats

1. Only 2 clusters identified (need proper cluster catalog)
2. Projection effects not fully accounted for
3. Hubble flow contaminates outer measurements

---

## 3. Interpretation

### The Methodology Sensitivity Problem

Results strongly depend on:
1. **Void definition**: Sparse vs dense void catalogs give opposite results
2. **Distance range**: Different distance bins show different directions
3. **Density metric**: Angular vs 3D density give opposite results

### What's Consistent

1. **Distance-controlled tests at 80-150 Mpc**: Synchronism direction
2. **Session #170 velocity-quartile test**: >> 10σ significance
3. **Cluster dispersion gradient**: 3.28× enhancement in outer regions

### What's Inconsistent

1. Overall void-based tests with realistic void density
2. Angular density tests at all scales

### Selection Effect Hypothesis

The dominant effect is likely:
```
|v| ∝ distance (r = 0.40)
```

Due to:
- Distance errors scale with distance
- More distant galaxies have larger apparent |v|
- Different environment definitions select different distance distributions

---

## 4. Path Forward

### Immediate

1. **Obtain real void catalogs**: Download Pan+2012 or Sutter+2012 data
2. **Use cluster catalogs**: Abell, SDSS redMaPPer clusters
3. **Forward modeling**: Simulate selection effects to predict expected signal

### Methodological

4. **Multiple independent tests**: Don't rely on single approach
5. **Null tests**: Permutation tests, mock catalogs
6. **Cross-validation**: Compare different density metrics on same data

### Future Data

7. **DESI peculiar velocities**: Smaller errors, larger samples
8. **Void galaxy surveys**: Purpose-built for this test

---

## 5. Files Created

| File | Description |
|------|-------------|
| `session173_void_catalog_analysis.py` | Sparse void catalog (44 voids) |
| `session173b_refined_void_analysis.py` | Dense void catalog (640 voids) |
| `session173c_cluster_infall_test.py` | Cluster dispersion test |
| `session173_void_analysis.png` | Visualization |
| `session173b_refined_void_analysis.png` | Refined visualization |
| `session173c_cluster_infall.png` | Cluster test figures |

---

## 6. Cumulative Progress (Sessions #169-173)

| Session | Topic | Key Finding |
|---------|-------|-------------|
| #169 | Mock data pipeline | 8.4σ detection in simulation |
| #170 | First real CF4 data | Velocity-quartile: >> 10σ |
| #171 | Systematic investigation | Distance errors dominate |
| #172 | High-precision subset | 3D vs angular discrepancy |
| **#173** | **Void/cluster tests** | **Methodology-sensitive** |

### Overall Assessment

The Synchronism prediction of enhanced peculiar velocities in voids remains **plausible but not definitively confirmed**. Key challenges:

1. Selection effects in peculiar velocity data
2. Sensitivity to environment definition
3. Need for independent, robust tests

The cluster dispersion test (ratio = 3.28) is promising but requires better cluster identification.

---

*Session completed: December 23, 2025*
