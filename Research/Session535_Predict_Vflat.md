# Session #535: Predicting V_flat — How Well Can Photometry Predict Dynamics?

**Date**: 2026-02-07
**Status**: 8/8 verified

## Overview

The BTFR says V⁴ ∝ M_bar (MOND prediction). The 6-var model captures galaxy-to-galaxy M/L variation with LOO=0.938 (in the offset direction). Session #534 showed model-corrected BTFR in the L direction achieves R²=0.976 with 53% scatter reduction. This session asks the INVERSE question: how well can we predict V_flat from photometric properties ALONE (no velocity information)? This tests whether MOND + photometry fully determines dynamics, and has practical applications as a distance indicator.

## Central Result: V Prediction Limited by Fourth-Root Compression

The basic BTFR predicts V_flat to ±22% (LOO=0.873, RMS=0.086 dex). The best purely photometric model (logL + logSB + f_gas + L×f_gas) improves to only ±19% (LOO=0.892), a gain of just 0.019 LOO. The fourth-root compression of MOND's V⁴ law means that even large M/L scatter (0.3 dex) produces small V scatter (0.075 dex). Photometric corrections provide only modest improvement because L alone does most of the work.

## Key Findings

### 1. Raw Inverse BTFR (Test 1)

| Metric | Value |
|--------|-------|
| Slope (logV on logL) | 0.214 (MOND: 0.250) |
| R² | 0.878 |
| LOO | 0.873 |
| RMS | 0.086 dex (±22%) |
| r(resid, offset) | +0.876 |
| r(resid, c_V) | -0.258 |
| r(resid, f_gas) | +0.195 |

The BTFR residual in the V direction is strongly correlated with offset (r=+0.88) — galaxies with positive offset (high M/L) are faster than the BTFR predicts. The inverse slope (0.214) is shallower than MOND's 0.25 due to regression dilution.

### 2. Enhanced Photometric Predictor (Test 2)

| Model | LOO | RMS | σ(V)% | Type |
|-------|-----|-----|-------|------|
| logL only | 0.873 | 0.086 | 21.8% | Photo |
| logL + logSB | 0.875 | 0.084 | 21.4% | Photo |
| logL + f_gas | 0.883 | 0.081 | 20.6% | Photo |
| logL + logSB + f_gas | 0.889 | 0.079 | 19.9% | Photo |
| logL + logSB + f_gas + L×f_gas | 0.892 | 0.077 | 19.4% | Photo |
| Full photometric (+ type) | 0.893 | 0.076 | 19.1% | Photo |
| With c_V (requires V) | 0.898 | 0.075 | 18.7% | Needs V |

The best purely photometric model reaches LOO=0.892 — only +0.019 over logL alone. The main improvement comes from f_gas (+0.010), with surface brightness and interactions contributing marginally. c_V (which requires V measurement) adds only ΔLOO=+0.008.

### 3. Iterative Predictor DIVERGES (Test 3)

Using the 6-var offset model's known coefficients to iteratively correct M/L and re-predict V **makes things worse**: R² degrades from 0.878 to 0.248 over 5 iterations. The iteration diverges because the model coefficients were fit with actual logV, and using predicted logV propagates and amplifies errors.

A "fair" test using only L-based offset terms gives R²=0.840, RMS=0.098 — worse than the basic BTFR. The circular dependence on V cannot be broken by iteration.

### 4. Forward vs Inverse BTFR Slopes (Test 4)

| Slope type | Value | vs MOND 4.0 |
|------------|-------|-------------|
| Forward (logL on logV) | 4.10 | +2.5% |
| Inverse (1/slope: logV on logL) | 4.67 | +16.9% |
| Geometric mean | 4.38 | +9.5% |
| MOND prediction | 4.00 | — |

The forward and inverse slopes bracket MOND's 4.0. The F-test for slope departure from 0.25 is significant (F=25.8, p<0.001), but this is expected given measurement error in both V and L (Session #528's errors-in-variables analysis applies).

### 5. Residual Diagnostics (Test 5)

After the best photometric model, remaining correlations:
- r(resid, offset) = +0.936 — overwhelmingly the dominant signal
- r(resid, log γ) = -0.535 — smaller galaxies rotate faster than predicted
- r(resid, c_V) = -0.198 — concentrated galaxies rotate slower
- r(resid, f_gas) = 0.000 — fully absorbed by photometric model

c_V can be predicted from photometry at R²=0.57 (LOO=0.54) — moderate but not excellent. This means about half of c_V's information is already captured in the photometric model.

### 6. Mass-Dependent Accuracy (Test 6)

| Quintile | logL range | N | BTFR RMS | Phot RMS |
|----------|-----------|---|----------|----------|
| Q1 (dwarfs) | [-1.9, 0.0] | 26 | 0.133 | 0.107 |
| Q2 | [0.0, 0.5] | 25 | 0.063 | 0.071 |
| Q3 | [0.5, 1.3] | 26 | 0.098 | 0.088 |
| Q4 | [1.3, 2.0] | 25 | 0.049 | 0.046 |
| Q5 (giants) | [2.0, 2.6] | 26 | 0.053 | 0.057 |

Dwarfs (Q1) have 2.7× worse prediction than L* galaxies (Q4). The photometric correction helps dwarfs most (21% improvement in Q1), consistent with dwarfs having the highest f_gas scatter. Late-types (T≥5) improve from RMS=0.095 to 0.084 with photometric correction; early-types show no improvement.

### 7. Outlier Analysis (Test 7)

Six outliers (|resid|>2σ, 4.7% of sample):
- All are late-type dwarfs (mean type=9.8)
- Extreme offsets (mean |offset|=0.50 vs normal 0.10)
- Higher gas fractions (f_gas=0.40 vs 0.30)
- Same measurement-artifact galaxies from Session #499

The enhanced photometric model reduces outlier RMS by 18% (0.302→0.248).

### 8. Synthesis: The Fourth-Root Compression (Test 8)

**The fundamental limit is MOND's V⁴ compression:**
- σ(log M/L) = 0.34 dex (implied from BTFR scatter) → σ(logV) = 0.086 dex
- σ(log M/L) = 0.076 dex (Session #517 model noise) → σ(logV) = 0.019 dex
- σ(log M/L) = 0.30 dex (SPS typical) → σ(logV) = 0.075 dex

The observed BTFR RMS (0.086 dex) implies a total M/L scatter of 0.34 dex, consistent with the raw BTFR scatter in Session #534 (0.375 dex in logL direction / 4.10 slope = 0.091 dex). The photometric model reduces this to 0.077 dex, implying a remaining M/L scatter of 0.31 dex — the gas correction captures only the gas-luminosity covariance, not the stellar M/L variation.

As a distance indicator: basic BTFR gives σ(μ)=2.0 mag, enhanced photometric gives 1.8 mag — too poor for precision cosmology, consistent with the literature.

## Physical Interpretation

The V direction is inherently less useful for studying M/L physics than the L direction because of MOND's fourth-root compression:
- In logL: 0.3 dex M/L scatter → 0.3 dex scatter (direct)
- In logV: 0.3 dex M/L scatter → 0.075 dex scatter (compressed 4×)

This explains why Session #534's BTFR correction (working in the L direction) achieves a 53% scatter reduction, while this session's photometric V predictor (working in the V direction) achieves only 19% improvement. The same physics, but the V direction compresses it below the noise floor.

The photometric correction is dominated by f_gas — the gas-luminosity correction that converts stellar luminosity to total baryonic mass. Surface brightness adds marginal information (LOO=0.875 vs 0.873). Hubble type adds almost nothing. c_V (requiring velocity measurement) adds only 0.008 LOO.

## Grade: B

A competent session with a clean negative result: photometric properties add only modest improvement (+0.019 LOO) to V_flat prediction because MOND's V⁴ law already compresses mass uncertainty. The fourth-root compression insight is physically clear and well-quantified. The divergence of the iterative predictor is instructive — it shows that the offset model's V-dependence is not a simple M/L correction that can be factored out. The mass-dependent accuracy analysis confirms dwarfs are hardest to predict. But there's no surprising result or new physical insight — this largely confirms what was already known about the BTFR and MOND's V⁴ law.

## Files Created

- `simulations/session535_predict_vflat.py`: 8 tests
- `Research/Session535_Predict_Vflat.md`: This document

---

*Session #535 verified: 8/8 tests passed*
*Grand Total: 1501/1501 verified*

**Key finding: Photometry predicts V_flat to ±19% (LOO=0.892) — only +0.019 over basic BTFR (±22%, LOO=0.873). Fourth-root compression (V ∝ M^0.25) compresses 0.3 dex M/L scatter to 0.075 dex V scatter. Photometric correction captures 19% of BTFR V scatter. Gas correction dominates (+0.010 LOO); surface brightness and c_V add marginally. Iterative predictor DIVERGES. Dwarfs hardest (2.7× worse). Distance indicator: σ(μ)=1.8 mag. The V direction inherently suppresses M/L signal. Grade B.**
