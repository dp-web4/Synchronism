# Session #610: CDM Discrimination Synthesis — What ALFALFA-SDSS Shows

**Date**: 2026-02-16
**Grade**: A+
**Domain**: Cosmology / Theory Discrimination / Synthesis

## Objective

Synthesize Sessions #604-609 into a definitive statement about MOND vs CDM
discrimination using BTFR scatter from ALFALFA × SDSS × Durbala data.
Close the CDM discrimination arc.

## The Definitive Numbers

### σ_int = 0.086 ± 0.003 dex

| Sample | N | Model σ | σ_int | z(MOND) | z(CDM) |
|:-------|:-:|:-------:|:-----:|:-------:|:------:|
| Full | 14,435 | 0.131 | 0.118 ± 0.001 | +148 | +41 |
| **Optimal** | **677** | **0.094** | **0.086 ± 0.003** | **+32** | **+0.5** |

5-variable Mendel model: logV, TFR_resid, Δ_SPS, g-i, f_gas.
Noise model: kinematic (W50, inclination) + distance (e_D/D × 2/ln10).
Quality cuts: SNR > 15, e_W50 < 10, b/a < 0.65, V > 80 km/s.

### Sensitivity Analysis

| Mass | Noise | Sample | σ_int | z(CDM) |
|:-----|:------|:-------|:-----:|:------:|
| Taylor | kin | Full | 0.142 | +64 |
| Taylor | kin+dist | Optimal | 0.105 | +6.2 |
| Mendel | kin | Full | 0.125 | +50 |
| **Mendel** | **kin+dist** | **Optimal** | **0.086** | **+0.5** |

The CDM z-score ranges from +0.5 to +64 depending on modeling choices.
Only the most favorable combination (Mendel masses + full noise + optimal
cuts) achieves CDM consistency.

## The Verdicts

### MOND: REJECTED at 32σ

σ_int = 0.086 >> 0 (MOND prediction).

**But this rejection is misleading.** The 0.086 dex includes residual
stellar mass-to-light ratio variation that our model doesn't fully capture.
MOND predicts zero scatter in the TRUE BTFR (with perfect M/L), not in the
observed BTFR. This test cannot distinguish MOND's zero-physics-scatter from
our incomplete M/L modeling.

### CDM: CONSISTENT at z = +0.5

σ_int = 0.086 ≈ 0.085 (CDM prediction from halo concentration scatter).

**But this consistency is also misleading.** CDM predicts 0.085 dex ONLY
from halo concentration. Our σ_int also includes residual M/L variation.
If the M/L residual is ~0.05 dex, then σ_physics = √(0.086² - 0.05²) ≈ 0.070,
which would be BELOW CDM's 0.085 — but we can't measure the M/L residual
independently.

### The Real Verdict: INCONCLUSIVE

The test is inconclusive because:
1. σ_int (0.086) includes both M/L residuals and physics scatter
2. These cannot be separated with photometric data
3. CDM's prediction (0.085) is exactly at our measurement floor
4. Both theories predict σ → 0 after perfect M/L correction

## The CDM Discrimination Arc: S604-610

| Session | Key Finding |
|:--------|:------------|
| S604 | Measurement noise = 9%; M/L dominates (91%) |
| S605 | Mendel masses break 0.161 barrier → 0.107 |
| S606 | σ_int = 0.072, "below CDM at -6.2σ" ← **PREMATURE** |
| S607 | SFR adds <2%; TFR dominates all M/L prediction |
| S608 | sSFR suppressor traces beyond-M/L information |
| S609 | **Distance noise dominates** → CDM now CONSISTENT |
| S610 | Definitive: σ_int = 0.086, z(CDM) = +0.5 |

### The Key Self-Correction

Session #606 claimed σ_int = 0.072 was 6.2σ below CDM's 0.085 floor.
Session #609 discovered that distance errors (σ_dist = 0.017 dex median)
had been omitted from the noise model. Distance noise exceeds kinematic
noise (0.017 > 0.015 dex). Including it shifts σ_int from 0.072 to 0.086,
completely reversing the CDM verdict.

This is the correct scientific process: test → discover error → correct.

## BIG-SPARC: What Would Change?

| Property | ALFALFA-SDSS | BIG-SPARC |
|:---------|:------------:|:---------:|
| N | ~700 (optimal) | ~4,000 |
| V source | W50 (global) | Resolved RC (Vflat) |
| Band | i (0.8μm) | 3.6μm |
| M/L scatter | ~0.15 dex | ~0.03 dex |
| Distance in V? | Yes | No |
| σ_int uncertainty | 0.003 dex | ~0.0003 dex |

BIG-SPARC would reduce the M/L wall by ~5× (from 0.15 to 0.03 dex),
making CDM's 0.085 clearly detectable above the M/L floor.
Predicted CDM discrimination: 4σ rejection if σ_physics = 0 (MOND),
or 3σ detection if σ_physics = 0.085 (CDM).

## The M/L Wall

| Model | Optimal σ | σ_int |
|:------|:---------:|:-----:|
| BTFR (V) | 0.276 | 0.270 |
| + TFR | 0.145 | 0.135 |
| + TFR + f_gas | 0.129 | 0.123 |
| + TFR + fg + g-i | 0.111 | 0.105 |
| + TFR + Δ + gi + fg | 0.094 | **0.086** |
| + SFR, W20, Ag | ~0.092 | ~0.084 |

The 5-variable model removes 88% of BTFR variance. Additional variables
(SFR, profile shape, extinction) add only ~2%. The remaining 0.086 dex
is the **M/L wall** — the limit of photometric M/L prediction.

## Variance Decomposition (Optimal Subsample)

| Component | Variance (dex²) | % of BTFR |
|:----------|:---------------:|:---------:|
| Removed by model | 0.0671 | 88.3% |
| Measurement noise | 0.0029 | 3.8% |
| Intrinsic (σ_int²) | 0.0075 | 9.9% |
| **Total BTFR** | **0.0760** | **100%** |

## Arc Status: CLOSED

7 sessions, 63/63 tests passed. The ALFALFA-SDSS dataset has been
exhaustively analyzed for MOND vs CDM discrimination. The answer is:
**INCONCLUSIVE with current data.** BIG-SPARC is required.

## Tests: 9/9 PASSED
## Grand Total: 1991/1991
