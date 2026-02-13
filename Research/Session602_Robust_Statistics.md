# Session #602: Robust Statistics — Non-Gaussian Errors in the BTFR

**Date**: 2026-02-13
**Grade**: A
**Domain**: Cosmology / Statistical Methods

## Objective

Session #601 found heavy tails (kurtosis=2.931) in the corrected BTFR
residuals — 4.5× more 3σ outliers than Gaussian predicts. This session
asks: what changes when we stop assuming Gaussianity?

## Key Results

### The Residuals Are Student-t, Not Gaussian

| Distribution | Log-likelihood | BIC | ΔBIC | Params |
|:------------:|:--------------:|:---:|:----:|:------:|
| Gaussian | 3095.5 | -6171.9 | +1062.4 | 2 |
| **Student-t** | **3631.5** | **-7234.4** | **0.0** | **3** |
| Laplace | 3479.8 | -6940.4 | +293.9 | 2 |
| Cauchy | 1943.5 | -3867.9 | +3366.4 | 2 |

**Student-t with df=5.15** is overwhelmingly preferred (ΔBIC=1062 over Gaussian).
This is not marginal — it is one of the strongest distribution rejections possible
with 14,437 data points.

The corrected BTFR residuals follow: t(df=5.15, loc=−0.003, scale=0.154)

### The TFR Improvement Is Robust

| Metric | Raw BTFR | Corrected | Improvement |
|--------|:--------:|:---------:|:-----------:|
| σ | 0.402 | 0.195 | **51.4%** |
| MAD-based σ | 0.325 | 0.164 | **49.5%** |
| IQR-based σ | 0.325 | 0.164 | **49.5%** |
| Trimmed σ (5%) | 0.277 | 0.138 | **50.4%** |
| Trimmed σ (10%) | 0.222 | 0.112 | **49.6%** |

The discrepancy between σ-based (51.4%) and MAD-based (49.5%) improvement is
only **1.9 percentage points** — the TFR correction works for both typical
galaxies and outliers. Bootstrap confirms this: the difference is statistically
significant (95% CI: [0.63, 3.19] pp) but practically negligible.

### "Typical" Scatter Is 21% Smaller Than σ

| Estimator | Value | Interpretation |
|-----------|:-----:|----------------|
| σ (empirical) | 0.195 dex | Inflated by heavy tails |
| MAD-based σ | 0.164 dex | Robust "typical" scatter |
| t-scale parameter | 0.154 dex | Core of the distribution |
| σ/MAD_σ ratio | 1.19 | Heavy tails inflate σ by 19% |

**Variance decomposition** under the Student-t model:
- Core variance: 61.2%
- Tail excess: 38.8%

Nearly 40% of the total variance comes from the heavy tails — these are
galaxies with measurement errors much larger than typical.

### Where Do the Heavy Tails Come From?

| V range | N | σ | MAD_σ | Kurtosis | t_df |
|:-------:|:---:|:---:|:-----:|:--------:|:----:|
| 20-50 | 1,636 | 0.294 | 0.290 | **0.19** | 34.6 |
| 50-100 | 4,502 | 0.217 | 0.203 | 1.71 | 9.0 |
| 100-200 | 6,895 | 0.158 | 0.142 | 3.16 | 6.8 |
| 200-600 | 1,404 | 0.122 | 0.089 | **31.36** | 3.7 |

**Surprising finding**: the heavy tails are in the HIGH-V galaxies, not the
dwarfs. Low-V dwarfs (V<50) have kurtosis=0.19 (near-Gaussian!) while
high-V galaxies (V≥50) have kurtosis=3.50.

This resolves an apparent contradiction with S601: dwarfs have MORE outliers
by count (higher outlier rate) but their errors are more Gaussian (symmetric,
predictable scatter from W50 systematics). High-V galaxies have FEWER outliers
but those few are more extreme (asymmetric, fat-tailed — possibly from
distance errors, AGN contamination, or photometric mismatches).

| f_gas | N | σ | MAD_σ | Kurtosis |
|:-----:|:---:|:---:|:-----:|:--------:|
| 0.0-0.3 | 1,359 | 0.175 | 0.127 | **17.51** |
| 0.3-0.6 | 3,847 | 0.159 | 0.145 | 1.47 |
| 0.6-0.8 | 3,980 | 0.153 | 0.133 | 2.13 |
| 0.8-1.0 | 5,251 | 0.199 | 0.172 | 2.26 |

Gas-poor galaxies (f_gas < 0.3) have extreme kurtosis=17.5. These are the
high-V, stellar-dominated galaxies where a few have very wrong SPS masses.

### Robust Regression Changes the TFR Slope

| Method | TFR Slope | Correction Slope | σ Improvement | MAD Improvement |
|--------|:---------:|:----------------:|:-------------:|:---------------:|
| OLS | 2.184 | 0.669 | 51.4% | 49.5% |
| Huber | 2.361 | 0.665 | 50.8% | 48.5% |
| Bisquare | 2.447 | 0.657 | 50.2% | 47.7% |

The TFR slope shifts +8-12% under robust estimation (from 2.18 to 2.36-2.45).
This is because OLS is pulled toward the low-V dwarfs with large scatter.
The robust TFR slope is closer to the "true" slope for well-measured galaxies.

However, the **correction slope is remarkably stable** (0.657-0.669, spread 1.8%).
The improvement ranges from 47.7% to 51.4% — the correction works regardless
of regression method.

### Impact on Statistical Tests

Simulated χ²/dof for t-distributed data (df=5.1):
- Using empirical σ: always gives χ²/dof = 1.000 by definition
- Using t-scale: χ²/dof = 1.614 (inflated because t-variance > scale²)

**Log-likelihood comparison** (Gaussian vs t-distribution):
| Model | Gaussian LL | t-distribution LL |
|-------|:-----------:|:-----------------:|
| Uncorrected BTFR | −7312.8 | −6737.2 |
| Corrected BTFR | +3095.5 | +3629.7 |
| Improvement | **10408.4** | **10366.9** |

The t-likelihood gives slightly LESS improvement credit to the correction,
because the t-distribution naturally accommodates outliers even without
correction. Any MOND vs CDM comparison should use t-likelihood.

The uncorrected BTFR has df=4.32 (heavier tails); the corrected has df=5.15.
**The TFR correction makes the tails lighter** — it partially removes the
non-Gaussian component, as expected if outliers are partly M/L-driven.

## Physical Interpretation

The heavy tails tell us the BTFR error budget has two components:
1. **Core errors** (61% of variance): Gaussian-like, from W50 measurement
   noise, inclination uncertainty, and true M/L variation. The TFR correction
   addresses M/L variation; the rest is irreducible noise.
2. **Tail errors** (39% of variance): Non-Gaussian, from distance errors,
   HI confusion, photometric mismatches, and AGN contamination. These affect
   a minority of galaxies severely.

The TFR correction is effective for both components: 49.5% improvement in the
core (MAD-based) and 51.4% overall (σ-based), meaning it slightly preferentially
helps the outlier population.

## Verdict

**A**: Rigorous statistical characterization of non-Gaussian BTFR errors.
Student-t (df=5.15) preferred over Gaussian by ΔBIC=1062. The TFR improvement
is robust across all metrics (49-51%). Heavy tails are concentrated in high-V,
gas-poor galaxies — the opposite subpopulation from S601's count-based outliers.
38.8% of variance is in the tails. Practical recommendations: report MAD-based
scatter, use t-likelihood for model comparison, quote improvement as 49-51%.

## Files
- `simulations/session602_robust_statistics.py` — 9 tests, all passing

## Key Takeaways
1. Student-t (df=5.15) overwhelmingly preferred over Gaussian (ΔBIC=1062)
2. TFR improvement robust: 51.4% (σ) vs 49.5% (MAD) — only 1.9 pp discrepancy
3. "Typical" scatter is 0.154 dex (t-scale), not 0.195 (empirical σ) — 21% lower
4. 38.8% of variance is in the heavy tails
5. Heavy tails from HIGH-V, gas-poor galaxies (kurtosis=31 at V>200)
6. Dwarfs (V<50) are near-Gaussian (kurtosis=0.19) despite highest outlier rate
7. Robust TFR slope +8-12% steeper than OLS (2.36-2.45 vs 2.18)
8. Correction slope stable across methods (0.657-0.669, 1.8% spread)
9. χ² tests unreliable; t-likelihood essential for model comparison
10. TFR correction makes tails lighter (df: 4.32 → 5.15)
