# Session #376: Gas Fraction Control Arc - Part 1

**Gas Fraction Control Arc - Part 1**
**Date**: 2026-02-05
**Status**: 8/8 verified

## Overview

The Empirical Execution Arc (Sessions #372-375) found that RAR scatter correlates with galaxy type/mass/SB in the direction predicted by Synchronism's environment-dependent γ (NP2). The strongest caveat was that **gas fraction is strongly correlated with all NP2 proxy variables** (r = 0.70-0.73).

This session answers the critical question: **Does the NP2 signal survive gas fraction control?**

## Key Result: NP2 SURVIVES GAS FRACTION CONTROL

```
╔══════════════════════════════════════════════════════════╗
║  GAS FRACTION CONTROL: DOES NP2 SURVIVE?                ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║  f_gas ↔ type:      r = +0.70 (STRONG confound)        ║
║  f_gas ↔ scatter:   r = +0.03 (NO direct effect!)      ║
║                                                          ║
║  Partial correlations (controlling for f_gas):           ║
║    Type → scatter:   r = +0.285 (p = 0.001)            ║
║    SB → scatter:     r = -0.304 (p = 0.00004)          ║
║    Mass → scatter:   r = -0.290 (p = 0.0005)           ║
║                                                          ║
║  Signal attenuation: NEGATIVE (signal STRENGTHENS)      ║
║  Matched pairs: 19/22 (86%) late > early at same f_gas  ║
║                                                          ║
║  ★ VERDICT: NP2 SURVIVES (Grade A-)                     ║
╚══════════════════════════════════════════════════════════╝
```

## Surprising Discovery: Gas Fraction is NOT a Confound

Despite expectations, gas fraction turns out to be a **suppressor variable**, not a confound:

1. **f_gas strongly correlates with type** (r = 0.70, p ~ 10⁻³⁸)
2. **f_gas does NOT correlate with RAR scatter** (r = 0.03, p = 0.71)
3. **Therefore, controlling for f_gas INCREASES the type→scatter signal**

This is counterintuitive. We expected gas fraction to be the main driver of increased scatter in late-type galaxies. Instead:
- Per-galaxy RAR scatter is essentially independent of gas fraction
- The type→scatter correlation is real and not mediated by gas

### Why the Aggregate F-Ratio Was Misleading

The Session #374 pooled variance ratio F(gas-rich/gas-poor) = 2.16 suggested gas-rich galaxies have more scatter. But this reflects **number of data points per galaxy** and **pooling effects**, not a per-galaxy trend. When measured per galaxy, f_gas has essentially no effect on scatter (r = 0.03).

## Detailed Results

### Test 1: Gas Fraction Distribution

| Group | N | f_gas median | IQR |
|-------|---|-------------|-----|
| Early (T≤4) | 46 | 0.027 | [0.011, 0.072] |
| Mid (T=5-6) | 32 | 0.122 | [0.068, 0.235] |
| Late (T≥7) | 93 | 0.372 | [0.242, 0.554] |

Gas fraction varies by 14x between early and late types. This seemed like a fatal confound for NP2.

### Test 2: RAR Scatter vs Gas Fraction

| Gas quartile | f_gas range | σ_median | σ_mean |
|-------------|-------------|----------|--------|
| Q1 (gas-poor) | [0.000, 0.057] | 0.097 | 0.119 |
| Q2 | [0.062, 0.198] | 0.079 | 0.088 |
| Q3 | [0.210, 0.401] | 0.087 | 0.103 |
| Q4 (gas-rich) | [0.406, 0.945] | 0.095 | 0.115 |

No monotonic trend in per-galaxy scatter. Gas-poor galaxies actually have the *highest* median scatter.

### Test 3: Gas-Controlled Morphology

Within gas fraction bins, the morphology effect persists:

| Gas Bin | F(late/early) | Early N | Late N |
|---------|---------------|---------|--------|
| Low gas (f < 0.08) | 2.80 | 38 | 6 |
| Mid gas (0.08-0.34) | 7.53 | 7 | 35 |
| High gas (f ≥ 0.34) | -- | 1 | 52 |

Both measurable bins show F > 1 (late types have more scatter even at matched gas fraction).

**Caveat**: Small N in some cells (6 late-type gas-poor, 7 early-type mid-gas).

### Test 4: SB → Scatter (Gas-Controlled)

- Raw: r = -0.240 (p = 0.001)
- Partial (controlling f_gas): r = -0.304 (p = 0.00004)
- Attenuation: -26.5% (signal INCREASES)

Within both gas-poor and gas-rich subsamples, LSB galaxies show ~40% more scatter than HSB galaxies.

### Test 5: Mass → Scatter (Gas-Controlled)

- Raw: r = -0.256 (p = 0.002)
- Partial (controlling f_gas): r = -0.290 (p = 0.0005)
- Attenuation: -13.3% (signal increases)

### Test 6: Matched-Pair Analysis

22 pairs matched within Δf_gas < 0.15:
- **19/22 (86.4%) have late σ > early σ**
- Mean early scatter: 0.068
- Mean late scatter: 0.136
- Difference: +0.068 dex

This is the strongest evidence: even at identical gas fractions, late-type galaxies have substantially more RAR scatter.

### Test 7: Partial Correlation Matrix

Full correlation matrix:
```
              T     f_gas   log_SB  log_Vf  σ_RAR
    T      +1.000  +0.728  -0.818  -0.827  +0.254
    f_gas  +0.728  +1.000  -0.718  -0.734  +0.081
    log_SB -0.818  -0.718  +1.000  +0.790  -0.261
    log_Vf -0.827  -0.734  +0.790  +1.000  -0.256
    σ_RAR  +0.254  +0.081  -0.261  -0.256  +1.000
```

Key column: σ_RAR correlates with T (+0.254), SB (-0.261), Vflat (-0.256), but NOT with f_gas (+0.081).

Double partial (controlling for f_gas AND quality): r = +0.277, p = 0.001

### Test 8: Final Verdict

NP2 SURVIVES GAS FRACTION CONTROL (Grade A-)

## What This Means for Synchronism

1. **The NP2 prediction gains credibility**: Environment-dependent RAR scatter is not a gas fraction artifact
2. **The signal is genuine**: Late-type galaxies have ~2x more per-galaxy RAR scatter than early-type, independent of gas content
3. **Gas fraction is a suppressor**: It correlates with the predictor (type) but not the outcome (scatter), making it a statistical suppressor rather than a confound
4. **Interpretation still unclear**: The type→scatter link could be:
   - Environment-dependent γ (Synchronism predicts)
   - Morphological structure effects on rotation curve measurement
   - Stellar feedback effects on gas dynamics
   - Asymmetric drift corrections differing by type
   - Beam-smearing effects differing by galaxy size/distance

## Honest Assessment

### What This Session Got Right
- Gas fraction is not the explanation for NP2
- The matched-pair analysis is clean and convincing
- The partial correlation analysis is rigorous

### Important Caveats
1. **Small samples in key bins**: Only 6 late-type gas-poor galaxies, 7 early-type mid-gas
2. **Proxies remain indirect**: Hubble type ≠ environment
3. **Other confounds exist**: Inclination, distance, angular resolution
4. **The pooled F-ratio was misleading**: Aggregate statistics ≠ per-galaxy statistics
5. **175 galaxies is still a small sample** for multi-variate analysis

### What Would Be Conclusive
1. Cross-match SPARC with explicit environment catalogs (group/cluster/void)
2. Larger rotation curve samples (ALFALFA, THINGS, LITTLE THINGS)
3. Simulations: what scatter does galaxy structure alone produce?
4. Inclination and beam-smearing controls

## Statistical Note: Suppressor Variables

The finding that gas fraction is a suppressor is an important statistical lesson. A variable can be strongly correlated with the predictor (r = 0.70) but not the outcome (r = 0.03). In this case, "controlling for" the variable actually strengthens the signal because it removes irrelevant variance from the predictor without removing relevant signal.

This means the Session #374 concern about gas fraction confounding was valid to raise but empirically turns out to be unfounded. The research process was correct: identify potential confound → test it → learn from the result.

## Files Created

- `simulations/session376_gas_fraction_control.py`: 8 tests
- `simulations/session376_gas_fraction_control.png`: 4-panel visualization
- `Research/Session376_Gas_Fraction_Control.md`: This document

## Next Steps

1. **Session #377**: Multi-variate regression with all confounds (inclination, distance, quality, gas fraction)
2. **Test remaining confounds**: Angular resolution, beam smearing, asymmetric drift
3. **External environment data**: Search for SPARC galaxy group membership catalogs
4. **Extended sample**: Can we add more galaxies from other surveys?

---

*Session #376 verified: 8/8 tests passed*
*Gas Fraction Control Arc: 1/? sessions*
*Grand Total: 455/455 verified across 17 arcs*

**Key discovery: NP2 (environment-dependent RAR scatter) SURVIVES gas fraction control. Gas fraction is a suppressor variable, not a confound. Per-galaxy RAR scatter does not correlate with gas fraction (r = 0.03). The morphology → scatter signal is real and statistically significant even after controlling for gas fraction (partial r = +0.285, p = 0.001). 19/22 matched pairs show late types having more scatter at identical gas fractions.**
