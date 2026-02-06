# Session #398: Quantitative Calibration Arc Synthesis

**Date**: 2026-02-06
**Status**: Synthesis (no tests)

## Arc Overview: Sessions #395-397

This arc tested whether Synchronism's specific quantitative prediction (γ = 2/√N_corr) matches the data. The answer: the qualitative structure is strongly confirmed but the specific formula is wrong, and the correction should be LOCAL, not global.

## Arc Timeline

| Session | Topic | Grade | Key Finding |
|---------|-------|-------|-------------|
| #395 | Functional Form | B | All forms equivalent at fixed V+L; cannot discriminate; γ=2/√N_corr amplitude 5-10× too large |
| #396 | Recalibration | B+ | Slope is 32% of theory; α ≈ 0.26 (not 0.5); radial dependence discovered; multiplicative confirmed |
| #397 | Radial Profile | **A-** | LOCAL N_corr(r) far superior (r=0.78 vs 0.59); **30% scatter reduction**; MOND-specific; measurement error confound |

## The Discovery Narrative

1. **Session #395 (Testing the Formula)**: We tested γ = 2/√N_corr directly against the data. ALL functional forms (1/√N_corr, 1/N_corr, R_eff, R/V²) give identical results at fixed V+L because they reduce to monotonic functions of R_eff. The amplitude of 2/√N_corr is dramatically too large, predicting 0.1-0.8 dex offsets where observations show 0.05-0.15 dex.

2. **Session #396 (Recalibrating)**: Fitting the actual amplitude gives ε ≈ -0.63 (effective, from linear model), compared to theory's 2.0. The slope is 32% of predicted. Best-fit power law index α = 0.26 (not 0.5). But a crucial discovery: the correction varies WITHIN galaxies (r = +0.24 with radius, p = 10⁻¹⁴).

3. **Session #397 (Going Local)**: Following the radial finding, we test a LOCAL N_corr model: N_corr(r) = V(r)²/(r × a₀). This is dramatically better than global N_corr (r = 0.78 vs 0.59 with RAR residual). It reduces per-galaxy RAR scatter by **30%** — the first practical improvement. But measurement uncertainty (e_vobs/V_obs) is correlated with radius, partially confounding the result.

## Quantitative Results

### What the Theory Predicted vs What We Found

| Quantity | Theory | Data | Agreement |
|----------|--------|------|-----------|
| Direction of size effect | Negative (larger → below RAR) | Negative | ✓ |
| Correction is multiplicative | Yes | Yes | ✓ |
| MOND-specific | Yes | Yes | ✓ |
| ε (normalization) | 2.000 | ~0.6 | ✗ (3× off) |
| α (power law index) | 0.500 | ~0.26 | ✗ |
| Global N_corr sufficient | Yes | No — LOCAL is far superior | ✗ |
| Scatter reduction | Yes (predicted) | **30% with local model** | ✓ (direction) |

### The Local vs Global Revolution

| Model | r with RAR residual | RMS | Per-galaxy scatter |
|-------|--------------------|----|-------------------|
| Global N_corr | 0.593 | 0.184 | 0.114 dex (0% improvement) |
| **Local N_corr(r)** | **0.779** | **0.143** | **0.079 dex (30% improvement)** |
| No model | — | 0.235 | 0.114 dex |

The local model improves over global by:
- r: 0.593 → 0.779 (+31%)
- RMS: 0.184 → 0.143 (-22%)
- Per-galaxy scatter: 0.114 → 0.079 (-30%)

### The Measurement Error Confound

The radial trend is partially confounded with measurement uncertainty:
- r(r/R_eff, Δresid) raw = +0.243
- r(r/R_eff, Δresid | e_vobs/V_obs) = +0.091

This doesn't eliminate the signal (still p = 0.005) but reduces it by 63%. The true physical signal may be 40-60% of the observed local N_corr improvement.

## Implications for Synchronism Theory

### What Must Change
1. **γ = 2/√N_corr is quantitatively wrong** — the coefficient "2" must be ~0.3-0.6
2. **The power law may not be α = 0.5** — data prefer α ≈ 0.26 or logarithmic
3. **N_corr should be LOCAL, not global** — the coherence scale depends on the local environment within the galaxy

### What This Suggests for the Theory
The local model N_corr(r) = V(r)²/(r × a₀) can be rewritten as:

N_corr(r) = g_obs(r) / a₀ × (r/V(r)²) × V(r)² / r
           = V(r)² / (r × a₀)

At each radius, this is the ratio of the local centripetal acceleration to a₀. Where this is large (high V, small r), the system is in the Newtonian regime. Where it's small (low V, large r), the system is deep in the MOND regime.

The success of the LOCAL model suggests that the coherence correction operates POINT-BY-POINT, not as a galaxy-wide property. This is more consistent with a modified gravity interpretation (where the modification depends on the local gravitational field) than with a global coherence interpretation (where the entire galaxy acts as one quantum system).

### Revised Prediction
Instead of:
> g_obs = g_RAR × (1 + 2/√N_corr)  [GLOBAL, one correction per galaxy]

The data suggest:
> g_obs(r) = g_RAR(r) × (1 + ε/√N_corr(r))  [LOCAL, varies with r]

where N_corr(r) = V(r)²/(r × a₀) and ε ≈ 0.3-0.6 (to be precisely calibrated after removing measurement error confound).

## Updated Master Scorecard

| ID | Prediction | Status | Key Evidence | Grade |
|----|-----------|--------|-------------|-------|
| NP1 | a₀ = cH₀/(2π) | SUPPORTED | 94% of MOND | B+ |
| NP2 | Morphology → scatter | STRONGLY SUPPORTED | p = 5×10⁻⁶ | A- |
| NP6 | N_corr → offset | SUPPORTED (LOCAL better) | r = 0.78 (local) | A- |
| NP7 | R_eff MOND-dominated | SUPPORTED | 3x at g < g† | A- |
| NP8 | R_eff M/L-independent | SUPPORTED | 7/7 tests | A- |
| NP9 | R_eff L-independent (late) | SUPPORTED | r = -0.49 | A- |
| NP10 | Dynamical radius confirms | SUPPORTED | R_max r = -0.47 | A |
| **NP11** | **γ = 2/√N_corr quantitative** | **NOT CONFIRMED** | **amplitude 3× too large** | **C** |
| **NP12** | **Local N_corr superior** | **SUPPORTED** | **30% scatter reduction** | **A-** |
| **NP13** | **Radial trend MOND-specific** | **SUPPORTED** | **late only, absent early** | **A-** |

## Arc Grade: B+

The arc is genuinely productive: it falsifies the specific γ formula, discovers the superiority of local N_corr, and achieves the first scatter reduction. But the measurement uncertainty confound prevents full credit.

## Arc Statistics

- Sessions: 3 (#395-397) + synthesis (#398)
- Tests: 24/24 verified
- Grand Total: 599/599 verified

---

*Grand Total: 599/599 verified (24 from this arc)*

**Arc summary: The specific formula γ = 2/√N_corr is NOT confirmed (amplitude 3× too large, power law index wrong). But the qualitative prediction (size matters in MOND, multiplicative, directionally correct) IS confirmed. The breakthrough discovery: LOCAL N_corr(r) = V(r)²/(r × a₀) is far superior to global N_corr (r = 0.78 vs 0.59), achieving a 30% per-galaxy scatter reduction — the first practical improvement to the RAR. The radial trend is specific to late types (MOND regime) but partially confounded with measurement uncertainty.**
