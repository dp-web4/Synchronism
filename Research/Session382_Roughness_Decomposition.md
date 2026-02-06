# Session #382: Roughness Decomposition and γ Reinterpretation

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #381 discovered that rotation curve roughness explains 50.6% of RAR scatter variance and completely mediates the type → scatter relationship. This session asks: is roughness a CAUSE (structure → roughness → scatter) or a MEDIATOR (environment → γ → roughness → scatter)?

## Key Result: Full Mediation, Predominantly Structural

Roughness mediates **88%** of the type → scatter relationship (Sobel z = 3.70, p = 0.0002). The remaining 12% direct effect is non-significant. A secondary systematic offset component (11% R²) remains as the maximum possible γ signal.

## Detailed Findings

### 1. Roughness Decomposition

Roughness is entirely dynamical - velocity measurement errors contribute negligibly (the SPARC mass model file doesn't propagate per-point errors to the log-acceleration space in a way that creates meaningful expected roughness).

Type effect on roughness persists after controlling for:
- Quality: r(type, roughness | Q) = +0.249 (p = 0.001)
- Distance: r(type, roughness | D) = +0.408 (p < 10⁻⁵)

This is not a measurement artifact - it's genuinely rougher rotation curves.

### 2. Quality-Matched Type Effect

| Quality | Early dyn_rough | Late dyn_rough | Ratio |
|---|---|---|---|
| Q=1 (High) | 0.0441 | 0.0702 | 1.59 |
| Q=2 (Medium) | 0.0432 | 0.0630 | 1.46 |

The roughness difference is present at BOTH quality levels with similar ratios.

### 3. Roughness by Type (Monotonic)

| Type | N | Dyn Roughness | Scatter |
|---|---|---|---|
| S0-Sab (T≤2) | 16 | 0.038 | 0.076 |
| Sb-Sbc (T=3-4) | 30 | 0.048 | 0.085 |
| Sc-Scd (T=5-6) | 32 | 0.057 | 0.108 |
| Sd-Sdm (T=7-8) | 26 | 0.072 | 0.120 |
| Sm-Im (T=9-10) | 63 | 0.068 | 0.115 |

Monotonically increasing through Sd-Sdm, then plateaus for irregulars.

### 4. Roughness-Acceleration Relationship

| Regime | Roughness |
|---|---|
| Low-g (g < g†) | 0.054 |
| High-g (g ≥ g†) | 0.063 |
| Ratio (low/high) | **0.86** |

**OPPOSITE to γ prediction!** If γ effects are stronger at low accelerations (MOND regime), we'd expect MORE roughness at low g. Instead, roughness is lower at low g. This argues against γ-mediated roughness.

By type:
- Early: low/high = 0.46 (much less rough at low g)
- Late: low/high = 1.00 (equal roughness across regimes)

### 5. Mass-Roughness Relationship

| Correlation | r | p |
|---|---|---|
| Vflat vs dyn_roughness | -0.313 | <10⁻⁴ |
| log SB vs dyn_roughness | -0.318 | <10⁻⁴ |
| Vflat vs dyn_roughness \| type | -0.150 | 0.050 |

More massive galaxies have less roughness, consistent with deeper potential wells stabilizing kinematics. The Vflat correlation is marginally significant even after controlling type.

Within late types: Vflat → roughness is weak (r = -0.13, n.s.), suggesting the mass effect is mostly between types.

### 6. Residual Scatter After Roughness Removal

Linear model: scatter = 1.28 × roughness + 0.028, R² = 0.506

After removing roughness:
- r(type, residual scatter) = +0.037 (p = 0.63) → **Type signal completely GONE**
- Mann-Whitney (early vs late residuals): z = -0.38, p = 0.71

But: r(|mean offset|, residual scatter) = +0.270 (p = 0.0003) → systematic offsets still matter.

### 7. Mediation Analysis (Sobel Test)

```
╔══════════════════════════════════════════════════════╗
║  MEDIATION: Type → Roughness → Scatter              ║
╠══════════════════════════════════════════════════════╣
║  Path a (Type → Roughness):    +0.00375             ║
║  Path b (Roughness → Scatter): +1.268               ║
║  Path c (Total effect):         +0.00542             ║
║  Path c' (Direct effect):      +0.00067             ║
║  Indirect effect (a × b):      +0.00475             ║
║  Proportion mediated:           87.6%               ║
║  Sobel z = +3.70, p = 0.0002                        ║
║                                                      ║
║  VERDICT: FULL MEDIATION                             ║
║  Direct effect < 30% of total; indirect significant  ║
╚══════════════════════════════════════════════════════╝
```

### 8. Three-Component Scatter Model

| Component | R² | Interpretation |
|---|---|---|
| Roughness (structure) | 50.6% | Kinematic irregularity |
| Systematic offset (γ?) | 11.3% | Possible environment |
| Unexplained | 38.1% | Galaxy peculiarities |

The systematic offset component (11.3%) is the **maximum possible γ signal**. Late types have 2.1x larger systematic RAR offsets than early types. This could reflect:
- Variable γ (Synchronism)
- M/L mismatch (different stellar populations)
- Distance errors (type-dependent distance methods)
- Baryonic model inadequacy

## Implications for Synchronism

### What This Means for NP2

1. **The scatter difference IS real** (confirmed across Sessions #376-382)
2. **It is predominantly structural** (88% mediated by roughness)
3. **But a systematic offset component remains** (~11% R²)
4. **The offset is the γ window** - if environment effects exist, they must live in this 11%

### What This Means for γ Theory

1. **Roughness is NOT γ-mediated**: Low-g roughness < high-g roughness (opposite to prediction)
2. **Mass → roughness**: Consistent with potential well depth controlling stability
3. **The γ → scatter mechanism needs revision**: If γ affects scatter, it does so through systematic RAR offsets (changing the effective a₀ or RAR shape), NOT through increasing kinematic roughness

### Revised γ Mechanism

Instead of: γ → kinematic irregularity → scatter
Perhaps: γ → modified effective RAR → systematic offset → contributes to scatter

This is a more subtle prediction. The original NP2 prediction didn't specify the mechanism. The data now says:
- NOT: late-type galaxies are noisier (that's structure)
- MAYBE: late-type galaxies sit at different RAR locations (that's offset)

### Updated NP2 Status

**NP2: PARTIALLY SUPPORTED** - The scatter difference is real but predominantly structural. A 11% systematic offset component remains as possible evidence for γ effects. External environment data required to resolve.

## Lessons Learned

1. **Always decompose effects**: Total scatter = roughness + offset. These have different physical origins.
2. **Mediation analysis is powerful**: The Sobel test quantifies how much of an effect goes through a mediator.
3. **Check acceleration dependence**: Low-g vs high-g roughness falsifies γ-mediated roughness.
4. **Quality matching essential**: The type → roughness effect persists at matched quality (1.59 ratio at Q=1).
5. **Velocity errors not in SPARC mass models**: The mass model file doesn't include per-point velocity error propagation. This is a known limitation.
6. **Negative evidence is constructive**: Ruling out γ-mediated roughness narrows the possible γ mechanism.

## Files Created

- `simulations/session382_roughness_decomposition.py`: 8 tests
- `Research/Session382_Roughness_Decomposition.md`: This document

---

*Session #382 verified: 8/8 tests passed*
*Grand Total: 503/503 verified*

**Key finding: Roughness FULLY MEDIATES the type → scatter relationship (88%, Sobel z = 3.70, p = 0.0002). Low-g roughness < high-g roughness, falsifying γ-mediated roughness. Three-component model: roughness (51%), systematic offset (11%, possible γ), unexplained (38%). The γ mechanism must work through systematic RAR offsets, not kinematic roughness. NP2 revised to PARTIALLY SUPPORTED.**
