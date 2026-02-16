# Session #608: The Suppressor Effect — Why sSFR Flips Sign

**Date**: 2026-02-16
**Grade**: A
**Domain**: Cosmology / Statistical Methodology

## Objective

Session #607 discovered that sSFR22 (22μm specific SFR) is NOT significant
raw (r=-0.012, p=0.21) but highly significant in the full model (t=+15.7)
with a POSITIVE coefficient — a sign flip. This session investigates the
suppressor mechanism and its implications for CDM discrimination.

## Key Results

### The Sign Flip is Real and Strong

| Context | r(sSFR, logMbar) | p |
|:--------|:----------------:|:-:|
| Raw | -0.177 | ~0 |
| Partial (V, TFR, g-i, f_gas controlled) | **+0.194** | 1.7e-99 |

Raw: high sSFR → low M/L → less stellar mass → negative Mbar residual.
After controlling M/L proxies (g-i, f_gas): the M/L channel is removed,
and what remains has the OPPOSITE sign.

### Not a Strict Classical Suppressor

| Criterion | Value |
|:----------|:-----:|
| r²(sSFR, logMbar) | 0.0314 |
| ΔR² when added to model | 0.0013 |
| Ratio (ΔR²/r²) | 0.04 |

The formal suppressor criterion (ΔR² > r²) is NOT met. sSFR is a
"semi-suppressor" — it contributes some direct prediction through the M/L
channel AND cleans other predictors. The sign flip is real but the mechanism
is partial overlap removal, not pure suppression.

### Δ_SPS is the Most Affected Predictor

| Predictor | r(X, sSFR) | Coefficient Δ |
|:----------|:----------:|:-------------:|
| Δ_SPS | -0.510 | **+33.5%** |
| g-i | -0.442 | +1.6% |
| f_gas | +0.318 | -0.6% |
| logV | -0.222 | +0.4% |
| TFR | -0.082 | -0.1% |

sSFR changes the Δ_SPS coefficient by 33.5% — far more than any other
predictor. This makes physical sense: Δ_SPS (Taylor-Mendel mass difference)
encodes SPS method disagreements that correlate strongly with sSFR
(r = -0.510). Adding sSFR cleans Δ_SPS of the part that's sSFR-redundant.

### Physical Mechanism

sSFR is 80.5% independent of g-i color. After removing color dependence,
sSFR still correlates with BTFR residuals (r = +0.083, p = 3e-19).

The color-independent sSFR likely encodes:
1. **Dust-embedded star formation** (22μm traces warm dust, invisible in optical)
2. **Recent gas accretion history** (more recent → higher sSFR at fixed color)
3. **In CDM**: halo assembly history → concentration → BTFR residual

### M/L vs CDM: Population Splits

| Population | σ(5-var) | σ(6-var) | sSFR gain |
|:-----------|:--------:|:--------:|:---------:|
| Gas-rich | 0.094 | 0.092 | 2.3% |
| **Gas-poor** | **0.023** | **0.022** | **3.7%** |
| Low-V | 0.139 | 0.133 | 4.3% |
| High-V | 0.072 | 0.071 | 1.6% |

Gas-poor galaxies benefit more from sSFR (+3.7% vs +2.3%), consistent with
sSFR primarily correcting **stellar** M/L rather than CDM halo scatter.
Low-V galaxies also benefit more (4.3% vs 1.6%).

### Velocity-Dependent σ_int

| V range | σ_int(5-var) | σ_int(6-var) | Δ% |
|:--------|:------------:|:------------:|:--:|
| 50-80 | 0.157 | 0.150 | +4.2% |
| 80-120 | 0.115 | 0.111 | +3.4% |
| 120-180 | 0.077 | 0.074 | +3.4% |
| 180-350 | 0.073 | 0.074 | **-2.4%** |

sSFR reduces scatter at V < 180 km/s but INCREASES it slightly at V > 180.
The V-dependence (Spearman r of |resid| vs logV) decreases from -0.291 to
-0.257 — a 12% reduction. This is consistent with sSFR partially absorbing
the CDM mass-dependent concentration scatter.

### Residual Structure (6-var Mendel)

- **Kurtosis**: 10.66 (very heavy tails, worse than S602's raw 5.15)
- **Skewness**: +1.93 (strong positive skew)
- **3σ outlier rate**: 1.51% (5.6× Gaussian)
- **Distance correlation**: r = +0.176, still the main systematic
- **logV²**: adds 1.1% (mild non-linearity)
- **logV × sSFR**: adds 1.9% (interaction between mass and star formation)

### Progressive σ_int (Mendel masses)

| Model | σ_int | z(CDM) |
|:------|:-----:|:------:|
| BTFR (V only) | 0.407 | +119 |
| + TFR | 0.156 | +65.5 |
| + TFR + g-i | 0.156 | +65.5 |
| + TFR + g-i + f_gas | 0.129 | +48.5 |
| + TFR + Δ + g-i + fg | 0.111 | +32.7 |
| + TFR + Δ + gi + fg + sSFR | 0.107 | +28.3 |

g-i adds exactly 0% beyond TFR (confirming S594). f_gas adds 17%.
Δ_SPS adds 14%. sSFR adds 4%. All of these reduce σ_int from the M/L
channel, not from CDM physics.

## Synthesis

### The Key Insight

The suppressor effect proves sSFR provides information BEYOND M/L
(partial r = +0.194 after controlling all M/L proxies). This is expected
in CDM (sSFR ↔ halo assembly) but harder to explain in MOND (no halo
property to correlate with).

However, the population split (gas-poor benefits more) suggests the
primary mechanism is still **M/L correction**, not CDM absorption.
sSFR likely captures dust-embedded star formation that g-i and f_gas miss,
improving the stellar M/L estimate.

### Revised CDM Interpretation

The model going below CDM's 0.085 floor (S606-607) should be interpreted as:
1. **Model variables correlate with halo concentration** — expected in CDM
2. **The correlation is INDIRECT** (through M/L, not direct halo measurement)
3. **σ_int < 0.085 does NOT refute CDM** — it means variables are not halo-independent
4. **A clean CDM test requires halo-independent variables** — none are available

### What This Means for MOND

In MOND, sSFR should only reduce scatter through M/L (no halos).
Finding a significant sSFR signal AFTER controlling M/L proxies (partial
r = +0.194) is consistent with MOND IF the signal comes from improved M/L
prediction (e.g., dust-embedded SF). The population split favors this
interpretation. MOND is NOT refuted by the suppressor effect.

## Tests: 9/9 PASSED
## Grand Total: 1973/1973
