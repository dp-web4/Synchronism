# Session #604: BTFR Error Budget — Forward Modeling the Noise Floor

**Date**: 2026-02-15
**Grade**: A
**Domain**: Cosmology / Error Analysis

## Objective

Forward-model the BTFR scatter from known measurement errors to determine
how much is noise vs intrinsic physics. Previous sessions estimated noise
statistically (S594: "52% noise") but never propagated actual per-galaxy
error bars through the BTFR formula.

## Method

Separate error sources into three categories:
1. **Kinematic** (pure measurement): e_W50, e_b/a → affect V_rot
2. **HI mass** (measurement): e_logMHI → affects M_gas
3. **SPS model** (systematic): e_logMsT → affects M_star

Propagate both analytically (partial derivatives) and via Monte Carlo (500 trials).
Use maximum likelihood to extract intrinsic scatter with per-galaxy errors.

## Key Results

### Measurement Noise is Only 9% of BTFR Variance

| Error Source | σ (dex) | % of BTFR variance |
|:-------------|:-------:|:-------------------:|
| W50 (line width) | 0.107 | 7.1% |
| b/a (inclination) | 0.026 | 0.4% |
| logMHI (HI mass) | 0.046 | 1.3% |
| **Measurement total** | **0.119** | **8.8%** |
| logMstar (SPS model) | 0.046 | 1.3% |
| **All errors** | **0.131** | **10.1%** |

W50 dominates kinematic noise (95% of kinematic variance). Analytic and
Monte Carlo estimates agree to 0.9%.

### ML Intrinsic Scatter = 0.375 dex

Using per-galaxy measurement noise (kin+HI only):
- σ_int = 0.3753 ± 0.0024 dex (captures M/L + physics)
- χ²/dof = 1.021 (well-calibrated)
- σ_int is 87% of observed BTFR variance

Including SPS model uncertainty:
- σ_int = 0.3723 ± 0.0024 dex (physics only)
- M/L contribution to scatter: 0.048 dex (1.4% of variance)

**Critical finding**: M/L scatter measured through e_logMsT is tiny (0.048 dex).
The real M/L variation is much larger — it's the *astrophysical* variation in
stellar mass-to-light ratio, not the *measurement* uncertainty in SPS fitting.

### TFR Correction Removes Both Noise and Intrinsic

| | σ_total | σ_meas | σ_int | meas% |
|:--|:-------:|:------:|:-----:|:-----:|
| BTFR | 0.401 | 0.119 | 0.375 | 8.8% |
| TFR-corrected | 0.195 | 0.051 | 0.182 | 6.7% |
| Reduction | 51.4% | 57.5% | 51.4% | — |

TFR correction reduces noise by 58% and intrinsic scatter by 51% — roughly
proportional. The 58% noise reduction occurs because V-L correlation partially
cancels kinematic errors: when V is overestimated, TFR predicts L should be
higher, which moves the BTFR prediction in the right direction.

### Noise Floor Achievable with Quality Cuts

| Quality cut | N | σ_meas (dex) | < CDM 0.085? |
|:------------|:-:|:------------:|:------------:|
| Full sample | 14,435 | 0.118 | No |
| e_W50 < 10 | 8,084 | 0.070 | YES |
| e_W50 < 5 | 2,465 | 0.069 | YES |
| Optimal (SNR>15, eW50<10, b/a<0.65, V>80) | 677 | 0.050 | YES |

**Best achievable noise floor: 0.050 dex** with the optimal subsample.
This is 0.6× the CDM concentration scatter prediction (0.085 dex), meaning
CDM-level signals are in principle detectable.

### But TFR Doesn't Remove Enough M/L Variation

In the optimal subsample after TFR correction:
- σ_int = 0.161 ± 0.005 dex
- This is 16.1σ above CDM's 0.085 dex prediction
- And 34σ above MOND's zero prediction

The residual 0.161 dex is NOT noise (noise floor is 0.046 dex). It's
**remaining M/L variation not captured by TFR**. The TFR removes 51% of
intrinsic scatter but significant M/L diversity persists at fixed V and L.

### Error Budget Varies with Velocity

| V range | N | σ_obs | σ_meas | noise% |
|:--------|:-:|:-----:|:------:|:------:|
| 20-50 | 1,635 | 0.699 | 0.158 | 5.1% |
| 50-80 | 2,535 | 0.470 | 0.134 | 8.1% |
| 80-120 | 3,992 | 0.339 | 0.111 | 10.8% |
| 120-180 | 4,096 | 0.275 | 0.102 | 13.7% |
| 180-350 | 2,126 | 0.261 | 0.103 | 15.5% |

Noise fraction *increases* with velocity (5% → 16%) because intrinsic scatter
decreases faster than noise at high V. But even at high V, 85% of scatter
is intrinsic.

### Reconciliation with S594

S594 reported "52% noise / 48% intrinsic". This session finds only 9% measurement
noise. The discrepancy: S594 defined "noise" as the TFR-corrected scatter
(0.195 dex), which includes residual M/L variation beyond the TFR model.
S604 separates *measurement noise* (W50+b/a+MHI errors) from *astrophysical
scatter* (real M/L variation). The "52% noise" from S594 is really
"52% = 9% measurement + 43% M/L variation removed by TFR".

### Error-Weighting and Non-Gaussianity

After dividing by per-galaxy errors, Student-t df improves from 12.3 (raw) to
12.4 (weighted) — negligible improvement. S602 found df=5.15 for raw residuals.
The discrepancy (12.4 vs 5.15) occurs because S602 used σ_btfr for all galaxies
while S604 uses per-galaxy σ_int, which absorbs more of the tail structure.

3σ outlier rate: 1.32% (vs 0.27% Gaussian) — 4.9× excess, confirming S602's
heavy-tail finding. Same AGC 251924 appears in the worst-χ² list.

## Synthesis

### The Nine-Word Summary
**Measurement noise is 9%; M/L dominates; TFR captures half.**

### What We Learned
1. BTFR scatter is 91% astrophysical, 9% measurement noise
2. W50 errors dominate kinematic noise (95% of kinematic variance)
3. Quality cuts can reduce noise floor to 0.050 dex (below CDM 0.085)
4. BUT residual intrinsic scatter after TFR (0.161 dex) still >> 0.085
5. TFR removes ~half of M/L variation, leaving the other half
6. MOND vs CDM discrimination requires removing M/L to < 0.085 dex residual
7. This means BIG-SPARC (resolved RCs) OR a better M/L predictor

### Implications for Future Work
The noise floor (0.050 dex with cuts) is already BELOW the CDM prediction.
The bottleneck is not measurement noise — it's **residual M/L variation**.
If we could predict M/L to better than 0.161 dex (compared to current 0.195→0.161
after TFR), CDM discrimination becomes possible even with W50 data.

This reframes the problem: not "we need better data" but "we need a better
M/L predictor". Colors fail (S594: 0% gain beyond TFR). What else could work?
- Multi-band TFR (S598: g-band adds 0%)
- SFR or metallicity from SDSS spectroscopy
- Structural parameters (Sérsic index, effective radius)
- Machine learning on full SDSS photometry

## Tests: 9/9 PASSED
## Grand Total: 1937/1937
