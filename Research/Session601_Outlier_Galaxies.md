# Session #601: Outlier Galaxies — Where Does the Predictor Fail?

**Date**: 2026-02-13
**Grade**: A-
**Domain**: Cosmology / ALFALFA-SDSS Outlier Analysis

## Objective

Identify and characterize 3σ outliers in the TFR-corrected BTFR across 14,437
ALFALFA-SDSS galaxies. Which galaxies does the M/L correction fail for, and
what makes them special?

## Key Results

### Outlier Census

| Threshold | N outliers | Rate | Gaussian expected |
|:---------:|:----------:|:----:|:-----------------:|
| >2σ | 734 | 5.08% | 4.55% |
| >3σ | 177 | 1.23% | 0.27% |
| >4σ | 36 | 0.25% | 0.01% |
| >5σ | 12 | 0.08% | ~0.00% |

The 3σ outlier rate is **4.5× the Gaussian expectation** (1.23% vs 0.27%).
The distribution has heavy tails: kurtosis = 2.931, skewness = +0.287.

### Two Outlier Classes

| Property | Over-massive (93) | Under-massive (84) | Normal (14,260) |
|----------|:-----------------:|:-------------------:|:---------------:|
| logV | 1.853 | 1.541 | 2.046 |
| logMbar | 9.900 | 7.936 | 9.995 |
| f_gas | 0.981 | 0.621 | 0.717 |
| g-i | 0.440 | 0.480 | 0.570 |
| Distance | 133.4 Mpc | 16.0 Mpc | 118.1 Mpc |
| W50 | 110 km/s | 50.5 km/s | 182 km/s |

**Over-massive outliers**: Gas-dominated (f_gas=0.981), moderate distance, blue.
These galaxies have more baryonic mass than expected at their V and L.
→ Possible: HI confusion, distance overestimate, true high M_gas.

**Under-massive outliers**: Nearby (16 Mpc), very low V and W50, low-mass.
→ Dominated by W50 systematics in nearby dwarfs (W50 ≠ V_flat).

### Velocity-Dependent Outlier Rate

| V range (km/s) | N | N outliers | Rate | Mean |residual| |
|:---------------:|:---:|:----------:|:----:|:---------------:|
| 20-50 | 1,636 | 89 | 5.4% | 0.232 |
| 50-80 | 2,535 | 46 | 1.8% | 0.181 |
| 80-120 | 3,992 | 22 | 0.6% | 0.145 |
| 120-180 | 4,096 | 14 | 0.3% | 0.116 |
| 180-300 | 2,052 | 6 | 0.3% | 0.088 |
| 300-600 | 126 | 0 | 0.0% | 0.092 |

The outlier rate drops 18× from low-V to high-V (5.4% → 0.3%). Dwarfs
dominate the outlier population: 89/177 (50%) are at V < 50 km/s.

### Distance-Dependent Outlier Rate

| D range (Mpc) | N | N outliers | Rate |
|:--------------:|:---:|:----------:|:----:|
| 5-25 | 395 | 84 | **21.3%** |
| 25-50 | 744 | 9 | 1.2% |
| 50-100 | 4,039 | 15 | 0.4% |
| 100-175 | 6,955 | 50 | 0.7% |
| 175-250 | 2,302 | 19 | 0.8% |

Nearby galaxies (D < 25 Mpc) have an extreme 21.3% outlier rate — these are
the nearby dwarfs where W50 is a poor proxy for V_flat.

### TFR-M/L Consistency

| Sample | r(TFR_resid, Bell M/L) | Median |M/L deviation| |
|--------|:----------------------:|:-----------------------:|
| All | 0.302 | 0.131 dex |
| Normal | 0.301 | 0.131 dex |
| Outliers | 0.300 | 0.168 dex (1.28×) |

Outliers show only 28% larger M/L deviation from the TFR→M/L relation than
normal galaxies. The TFR residual and Bell color give consistent M/L
estimates even for most outliers — they are NOT M/L anomalies but rather
measurement-error anomalies.

### Inclination and SNR Effects

| b/a range | N | Outlier rate | σ_corr |
|:---------:|:---:|:-----------:|:------:|
| 0.20-0.30 | 1,787 | 0.5% | 0.168 |
| 0.30-0.40 | 2,011 | 1.2% | 0.190 |
| 0.70-0.85 | 3,534 | 1.7% | 0.200 |

Nearly face-on galaxies (b/a > 0.7) have 3.4× the outlier rate of edge-on
galaxies (1.7% vs 0.5%). The inclination correction sin(i) diverges as
b/a → 1, amplifying W50 errors.

| SNR range | N | Outlier rate | σ_corr |
|:---------:|:---:|:-----------:|:------:|
| 6.5-10 | 8,192 | 0.9% | 0.186 |
| 10-20 | 5,062 | 1.4% | 0.196 |
| 20-50 | 1,108 | 2.3% | 0.233 |
| 50-100 | 70 | 5.7% | 0.282 |

Counter-intuitively, **higher SNR = more outliers**. This is a selection
effect: high-SNR detections are nearby gas-rich dwarfs (the same population
driving the distance and velocity trends).

### TFR Correction: What It Fixes and Creates

| Category | N |
|----------|:-:|
| Raw BTFR outliers | 192 |
| Corrected BTFR outliers | 177 |
| Persistent (outlier in both) | 44 |
| Fixed by TFR correction | 148 |
| Created by TFR correction | 133 |

The TFR correction **fixed 148 raw outliers** but **created 133 new ones**.
The 44 persistent outliers are measurement-error driven (median D = 12.5 Mpc).
The 148 "fixed" outliers were at high f_gas (0.537) and far distances — these
had M/L-driven BTFR scatter that the TFR residual successfully removed.

### Most Extreme Outliers

| Rank | AGC | σ | V_rot | logMbar | g-i | f_gas | D (Mpc) | Interpretation |
|:----:|:---:|:-:|:-----:|:-------:|:---:|:-----:|:-------:|----------------|
| 1 | 251924 | +9.1 | 273 | 12.46 | 3.80 | 0.01 | 163.8 | Extreme color (photometry error?) |
| 2 | 174681 | +7.8 | 82 | 12.02 | 3.15 | 0.00 | 133.8 | Extreme color, stellar-dominated |
| 3 | 189093 | +7.1 | 130 | 11.50 | 3.35 | 0.04 | 122.4 | Extreme color |
| 4 | 221639 | +6.7 | 42 | 9.19 | -0.03 | 1.00 | 15.1 | Nearby dwarf, pure gas |
| 5 | 323359 | +6.6 | 65 | 9.85 | -0.46 | 1.00 | 106.1 | Extreme blue, pure gas |

The top 3 outliers all have g-i > 3.0 — extreme colors far outside the
normal range (median 0.57). These are almost certainly photometric errors
or mismatches, not real astrophysics. The next tier (ranks 4-10) are
gas-dominated dwarfs (f_gas > 0.99) where the HI mass dominates the
baryonic budget.

### Follow-Up Targets

17 galaxies meet the criteria: outlier AND SNR > 20 AND 50 < V < 250 km/s.
These are the most promising targets for resolved rotation curve observations
that could reveal whether the predictor truly fails or the W50 measurement
is at fault.

## Physical Interpretation

The outlier population tells a clear story:

1. **Most outliers are measurement artifacts, not physics.** The 89 low-V
   dwarfs (V < 50 km/s) are dominated by W50 systematics — W50 is a poor
   proxy for V_flat in slowly rotating galaxies with turbulent gas motions.

2. **Nearby galaxies are most affected.** The 21.3% outlier rate at D < 25 Mpc
   reflects ALFALFA's sensitivity to nearby dwarfs that would be undetected
   at greater distances.

3. **The TFR correction is genuinely effective.** It fixed 148 outliers whose
   scatter was M/L-driven (removed by the V-L residual correction). The 133
   "new" outliers it created were hidden beneath the larger raw scatter.

4. **Heavy tails are real.** Kurtosis = 2.931 means the error distribution
   is strongly non-Gaussian. This matters for any statistical test assuming
   Gaussian errors (e.g., χ² tests in S595).

5. **The model does NOT fail for specific galaxy types.** It fails where the
   DATA fails — poor W50, extreme inclination, photometric errors.

## Verdict

**A-**: Clean identification and characterization of the outlier population.
The key finding is that outliers are measurement-driven, not physics-driven:
dominated by low-V dwarfs with W50 systematics, nearby galaxies, and
photometric errors. The TFR correction effectively removes M/L-driven
outliers but cannot fix measurement errors. Heavy tails (kurtosis = 2.931)
confirm non-Gaussian errors throughout the sample. 17 follow-up targets
identified for resolved observations.

## Files
- `simulations/session601_outlier_galaxies.py` — 9 tests, all passing

## Key Takeaways
1. 177 3σ outliers (1.23%) — 4.5× Gaussian expectation
2. Dominated by low-V dwarfs (89/177) — W50 systematics, not physics
3. 21.3% outlier rate at D < 25 Mpc (nearby dwarfs)
4. Over-massive outliers: gas-dominated (f_gas=0.981), distant
5. Under-massive outliers: nearby (16 Mpc), very low V
6. TFR correction fixed 148 raw outliers, created 133 new ones (revealed)
7. Top 3 most extreme outliers have g-i > 3.0 — photometric errors
8. Heavy tails (kurtosis=2.931) — errors are non-Gaussian
9. 17 promising follow-up targets for resolved rotation curves
