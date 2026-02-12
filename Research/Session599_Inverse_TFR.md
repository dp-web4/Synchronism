# Session #599: Inverse TFR — Predicting Rotation Velocity from Photometry

**Date**: 2026-02-12
**Grade**: A-
**Domain**: Cosmology / Distance Estimation

## Objective

Test whether the M/L correction that gives 51% improvement in the forward
BTFR also improves the INVERSE TFR (predicting velocity from photometry),
which is the form used for Tully-Fisher distance estimation.

## Key Results

### The Forward-Inverse Asymmetry

| Direction | Correction | σ (dex) | Improvement |
|-----------|-----------|--------:|----------:|
| Forward (logMbar vs logV) | TFR residual | 0.195 | **51.4%** |
| Inverse (logV vs logMbar) | Color | 0.161 | **1.2%** |

**The 51% forward improvement collapses to ~1-2% in the inverse direction.**

This asymmetry is the central finding. It arises because:
- Forward: V is known → TFR residual = direct M/L measurement → massive correction
- Inverse: V is unknown → only photometric proxies (color, f_gas) → weak correction

### Inverse Velocity Prediction Models

| Model | σ_LOO (logV) | ΔV/V |
|-------|------------:|---------:|
| logMbar (SPS) | 0.163 | 45.7% |
| logMbar + color | 0.161 | 45.0% |
| logL_i only | 0.171 | 48.3% |
| logL_i + color | 0.169 | 47.5% |
| logL_i + color + f_gas | 0.164 | 45.9% |
| logL_i + color + f_gas + logMHI | 0.160 | 44.5% |

The best photometric model (L+color+f_gas+MHI) gives 2.2% improvement over
simple Mbar — negligible for practical purposes.

### Forward vs Inverse TFR Slopes

| Relation | Slope | 1/Slope |
|----------|------:|--------:|
| Forward TFR (logL vs logV) | 2.18 | 0.46 |
| Inverse TFR (logV vs logL) | 0.23 (=1/4.31) | 4.31 |
| Bisector | 3.07 | 0.33 |

The product of forward and inverse slopes = r² = 0.507. The inverse TFR slope
(4.31) is much steeper than 1/forward (4.58) due to regression dilution.

### Velocity Precision by Mass Range

| V range (km/s) | N | σ(logV) | ΔV/V |
|----------------|----:|-------:|-----:|
| 20-50 | 1,636 | 0.201 | 58.9% |
| 50-100 | 4,502 | 0.138 | 37.4% |
| 100-150 | 4,453 | 0.097 | 25.1% |
| **150-250** | **3,413** | **0.090** | **23.0%** |
| 250-500 | 432 | 0.104 | 27.0% |

Distance precision improves from ~59% (dwarfs) to ~23% (spirals with
V=150-250 km/s). The sweet spot is L* to massive galaxies.

### Malmquist Bias

| Property | Forward | Inverse |
|----------|--------:|--------:|
| Skewness | +0.030 | -1.147 |
| Mean shift (near→far) | +0.834 | +0.020 |
| σ ratio (far/near) | 0.932 | 1.034 |

The forward TFR has massive Malmquist bias (mean shifts by 0.83 dex from
near to far galaxies). The inverse TFR is nearly immune (0.02 dex shift),
confirming its use in distance estimation.

### Peculiar Velocity Precision

At 100 Mpc (H₀=75 km/s/Mpc):
- σ(V_pec) ≈ 2800 km/s — far too large for individual galaxies
- Typical peculiar velocities: 300-500 km/s
- S/N per galaxy: ~0.1

**Conclusion**: ALFALFA W50 data cannot measure individual peculiar velocities.
Only statistical ensembles (100+ galaxies per bin) or resolved rotation curves
can achieve useful precision.

## Physical Interpretation

The forward-inverse asymmetry reveals a fundamental truth about the MOND
offset predictor: **it is a FORWARD tool, not an INVERSE tool**.

In the forward direction, observing V directly measures the gravitational
field, and the TFR residual tells you how much luminosity deviates from the
M/L=const expectation. This is a strong signal (r=0.87).

In the inverse direction, you have photometric quantities (L, color, f_gas)
but no dynamical information. Color constrains M/L weakly (r≈0.5) because:
1. Color-M/L has large intrinsic scatter (~0.2 dex)
2. Color doesn't capture SFH diversity at fixed mass
3. W50 errors dominate the velocity scatter (~0.16 dex)

The Mbar model slightly outperforms L+color because SPS mass estimation
already encodes color→M/L information. Adding color on top of SPS masses
is partially redundant.

## Verdict

**A-**: Clean analysis revealing the fundamental forward-inverse asymmetry
(51% vs 2% improvement). The V-dependent precision profile is well-mapped
(23% at V=150-250, 59% at V<50). The Malmquist bias quantification
confirms inverse TFR's advantage. The negative result (color barely helps)
is informative — it explains why TF distance surveys don't routinely use
M/L corrections.

## Files
- `simulations/session599_inverse_tfr.py` — 9 tests, all passing

## Key Takeaways
1. M/L correction: 51% forward improvement → only 2% inverse improvement
2. V is the information carrier — without V, M/L correction is weak
3. Best inverse precision: σ(logV) = 0.090 at V=150-250 km/s (ΔD/D=23%)
4. W50 data cannot measure individual peculiar velocities (S/N~0.1)
5. Mbar slightly outperforms L+color in inverse direction (SPS encodes color)
6. Forward TFR has 0.83 dex Malmquist bias; inverse has only 0.02 dex
