# Session #446: c_V Tracks in the Acceleration Plane

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

All previous analysis used galaxy-level offsets. This session tests whether c_V creates visible, distinct tracks in the point-level g_obs vs g_bar acceleration plane.

## Central Result: c_V Creates Distinct Tracks, Strongest in Inner Regions

| Region | Track separation (high - low c_V) |
|--------|----------------------------------|
| All MOND | **+0.078 dex** (p < 0.001) |
| Inner (r < R_eff) | **+0.227 dex** |
| Outer (r > 2R_eff) | +0.009 dex |

**The c_V effect is almost entirely an inner phenomenon.** The track separation is 25× larger in inner regions than outer regions. This is physically intuitive: mass concentration affects where the acceleration peaks (inner regions), while the outer asymptote is set by the total mass.

## Key Findings

### 1. Binned RAR by c_V Tercile (Test 1)

The mean separation (high - low c_V) across g_bar bins: **+0.11 dex**. The separation increases toward lower g_bar (deeper MOND regime), reaching +0.43 dex at log(g_bar) ≈ -9.6. At high g_bar (Newtonian regime), the separation is smaller.

### 2. Modified RAR Formula (Test 3)

A simple c_V-dependent RAR:
```
g_obs ≈ g_RAR × 10^(-0.121 + 0.105 × c_V)
```

The c_V-only correction improves point-level scatter by only 2.3% — because c_V alone captures only the geometric component (~5% of variance). The full V+L+c_V correction gives 21.5%.

The effective a₀ varies from 0.83×a₀ (low c_V) to 0.91×a₀ (high c_V) — a ~10% range.

### 3. Inner vs Outer Divergence (Test 5)

| Region | Low c_V mean | High c_V mean | Separation |
|--------|-------------|--------------|-----------|
| Inner (r < R_eff) | -0.199 | +0.028 | **+0.227** |
| Outer (r > 2R_eff) | -0.027 | -0.019 | +0.009 |

The inner separation is 25× the outer separation. This is the most direct evidence that c_V captures an inner mass distribution effect, not a total mass effect.

### 4. Combined c_V + BTFR Tracks (Test 6)

The 2×2 split (c_V × BTFR residual) creates four distinct RAR tracks:

| Quadrant | Mean residual |
|----------|-------------|
| High c_V, Below BTFR | **+0.098** |
| Low c_V, Below BTFR | -0.016 |
| High c_V, Above BTFR | -0.055 |
| Low c_V, Above BTFR | **-0.177** |

The extreme quadrants (high c_V + below BTFR vs low c_V + above BTFR) span 0.275 dex — over a factor of 1.9 in g_obs at fixed g_bar.

### 5. Galaxy-Level = Point-Level (Test 7)

| Approach | RMS (log g) | Improvement |
|----------|-----------|-------------|
| Standard RAR | 0.182 | — |
| Galaxy-level V+L+c_V | 0.143 | 21.5% |
| Point-level V+L+c_V fit | 0.141 | 22.3% |

The coefficients are nearly identical (c_V: 0.44 galaxy vs 0.41 point-level). The galaxy-level approach misses only 0.8% of the point-level fit — confirming the correction is a galaxy-level property, not a point-level one.

## Physical Interpretation

The acceleration plane visualization confirms the physical picture: c_V captures how concentrated the mass distribution is, which affects where the galaxy sits on the RAR primarily through its inner acceleration profile. In the outer parts (where the flat rotation curve dominates), all galaxies converge toward the same RAR regardless of c_V.

This is consistent with the MOND interpretation: the algebraic RAR (g_obs = g_bar/(1 - e^(-√(g_bar/a₀)))) assumes each point is independent, but in reality the inner acceleration is influenced by the global mass distribution. A concentrated mass (high c_V) creates a steeper inner potential, pushing g_obs above the algebraic prediction.

## Grade: B+

A clean visualization of the c_V effect in the acceleration plane. The inner vs outer separation (25×) is the highlight — it provides the most direct physical interpretation of what c_V captures. The 2×2 quadrant analysis shows the combined M/L + geometry effect spans nearly 0.3 dex. The galaxy vs point-level comparison confirms the correction's nature.

## Files Created

- `simulations/session446_acceleration_tracks.py`: 8 tests
- `Research/Session446_Acceleration_Tracks.md`: This document

---

*Session #446 verified: 8/8 tests passed*
*Grand Total: 933/933 verified*

**Key finding: c_V creates distinct tracks in the acceleration plane (separation +0.078 dex, p<0.001). The effect is 25× stronger in inner regions (0.23 dex) than outer (0.009 dex) — it's an inner mass concentration effect. Galaxy-level and point-level approaches give equivalent results (21.5% vs 22.3%). Combined c_V × BTFR tracks span 0.28 dex. Grade B+.**
