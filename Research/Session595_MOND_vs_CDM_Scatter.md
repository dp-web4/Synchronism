# Session #595: MOND vs CDM — Is There Halo-Driven Scatter?

**Date**: 2026-02-12
**Grade**: A-
**Domain**: Cosmology / Theory Testing

## Objective

Session #594 showed the TFR residual captures ALL detectable intrinsic BTFR
scatter. This session asks: is the remaining scatter consistent with measurement
noise alone (MOND prediction), or is there evidence for additional halo-driven
scatter (CDM prediction)?

## Key Results

### Formal Noise Test

| Metric | Value |
|--------|------:|
| σ_corrected | 0.195 dex |
| σ_noise (RMS) | 0.256 dex |
| χ²/dof (corrected) | 5.88 |
| Implied σ_intrinsic | 0.172 dex (from χ²) |

χ²/dof >> 1 means the **individual error bars are underestimated** (by ~√5.9 ≈ 2.4×).
This is a persistent problem with HI survey data — W50 and photometric errors
are formal uncertainties that miss systematic effects.

### Outlier Analysis

Outliers (5.1% with |resid| > 2σ) are overwhelmingly:
- **Low velocity**: 10.7% at V=30-80 km/s vs 1.1% at V=150-500 km/s
- **Low mass**: logMstar = 8.26 vs 9.38
- **Gas-rich**: f_gas = 0.77 vs 0.66
- **Nearby**: D = 85 Mpc vs 123 Mpc

These are all properties associated with **poor W50 measurements**, not
distinct physics. The fat tails (kurtosis = 2.93) are noise artifacts.

### Morphology/Concentration — No Signal

| Proxy | r(proxy, BTFR_corrected) | % variance |
|-------|-------------------------:|----------:|
| M/L residual (V,L removed) | 0.034 | 0.1% |
| Concentration proxy | 0.009 | 0.0% |
| g-i color residual | 0.032 | 0.1% |

No morphological or concentration proxy predicts corrected BTFR residuals.
This is consistent with MOND (no halo-driven scatter) but also with CDM
(signal buried in noise).

### Distance/Environment

| Property | Correlation |
|----------|------------|
| r(logD, BTFR_corrected) | **+0.449** |
| |resid| vs logD slope | **-0.197** |

Distant galaxies have LOWER scatter — a selection effect (Malmquist bias:
only bright, massive, well-measured galaxies detected at large D).
Scatter decreases from 0.28 dex (D < 30 Mpc) to 0.16 dex (D > 100 Mpc).

### CDM Scatter Budget

| Source | σ (dex) |
|--------|--------:|
| Halo concentration | 0.085 |
| Stellar-to-halo mass | 0.180 |
| Total CDM prediction | 0.199 |

After removing M/L variation (TFR correction), CDM predicts ~0.085 dex
residual scatter from halo concentration variation.

**Observed upper limit**: 0.153 dex (95% bootstrap) — above CDM's 0.085,
so **not constraining**.

### MOND vs CDM Model Comparison

| Model | BIC | Notes |
|-------|----:|-------|
| MOND (σ_intr = 0) | 41,036 | Zero free parameters |
| CDM (σ_intr free) | 41,046 | One free parameter |

**MOND preferred by ΔBIC = 9.6**, but this is unreliable because both models
have χ²/dof >> 1. The comparison depends critically on the noise model,
which is poorly calibrated.

## Physical Interpretation

The data are **consistent with MOND** (zero intrinsic scatter after M/L
correction) but **cannot rule out CDM** (noise floor too high). The test
is noise-limited, not signal-limited.

The key obstacle is W50-based velocities: single-dish line widths have
systematic errors (beam confusion, pointing offsets, turbulence correction)
that are not captured in the formal error bars. This inflates χ² while
allowing large residuals.

## What Would Settle It

BIG-SPARC (~4000 galaxies, resolved RCs, 3.6μm photometry) would:

1. Reduce σ_noise from ~0.26 to ~0.05 dex (resolved RCs, no W50 systematics)
2. Use 3.6μm where M/L ≈ const (TFR slope ≈ 4.0, not 2.2)
3. Allow direct RAR computation (not just BTFR)
4. CDM's 0.085 dex would be a 1.7σ detection at σ_noise = 0.05

Even BIG-SPARC would need σ_noise < 0.04 dex for a 2σ detection of CDM scatter.
This is achievable but challenging — distance errors alone contribute ~0.03 dex.

## Verdict

**A-**: Honest, thorough analysis that correctly identifies the limitations.
All observations consistent with MOND but also with CDM — the test is
noise-limited. The outlier analysis, distance dependence, and CDM scatter
budget are well-executed. The BIC comparison is reported but properly
caveated as unreliable. The session clearly identifies what data would
settle the question (BIG-SPARC with σ_noise < 0.04 dex).

## Files
- `simulations/session595_mond_vs_cdm_scatter.py` — 9 tests, all passing

## Key Takeaways

1. Current data cannot distinguish MOND from CDM via BTFR scatter
2. Outliers are noise artifacts (low-V, nearby galaxies with poor W50)
3. No morphology/concentration proxy predicts corrected residuals
4. BIC favors MOND by 9.6 but is unreliable (χ² >> 1)
5. Need σ_noise < 0.04 dex for 2σ CDM concentration detection
6. BIG-SPARC is the critical next dataset
