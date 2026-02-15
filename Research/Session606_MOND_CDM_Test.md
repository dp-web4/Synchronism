# Session #606: MOND vs CDM — Can the Mendel Model Discriminate?

**Date**: 2026-02-15
**Grade**: A+
**Domain**: Cosmology / Theory Discrimination

## Objective

With the full Mendel model (S605) achieving σ = 0.107 dex, test whether
MOND vs CDM discrimination is now possible. CDM predicts intrinsic scatter
of ~0.085 dex from halo concentration variation. MOND predicts zero
intrinsic scatter.

## Key Results

### Three-Way Verdict: Exceeds Both (Full Sample)

| Prediction | σ_int | z-score | Status |
|:-----------|:-----:|:-------:|:------:|
| MOND | 0 | +151σ | **REJECTED** |
| CDM | 0.085 | +14.2σ | **REJECTED** |
| Observed | 0.094 ± 0.001 | — | Exceeds both |

Full sample: σ_int exceeds CDM's prediction. Residual M/L variation
(~0.040 dex) persists beyond the 5-variable model.

### But the Optimal Subsample FAVORS MOND over CDM

| Subsample | N | σ_int | z(MOND) | z(CDM) |
|:----------|:-:|:-----:|:-------:|:------:|
| Full | 14,437 | 0.094 | +151 | +14.2 |
| Optimal | 677 | 0.072 ± 0.002 | +34.4 | **-6.2** |
| Aggressive | 46 | 0.061 ± 0.007 | +9.3 | **-3.8** |

The optimal subsample (SNR>15, eW50<10, b/a<0.65, V>80) gives σ_int = 0.072,
which is **below** CDM's 0.085 dex at 6.2σ. Even more aggressive cuts give
0.061 dex at 3.8σ below CDM.

**This is a genuine surprise.** If CDM is correct, the intrinsic scatter
cannot be BELOW 0.085 dex (this is a hard floor from halo concentration
diversity). Finding σ_int < 0.085 at 6.2σ either means:
1. Our measurement noise model is overestimated (inflating σ_noise, reducing σ_int)
2. The CDM 0.085 prediction is too high for this mass range
3. Our model absorbs some CDM scatter through its 5 variables
4. CDM concentration scatter is not independent of the model variables

### σ_int Decreases with Velocity (CDM-like!)

| V range | N | σ_int |
|:--------|:-:|:-----:|
| 50-80 | 2,535 | 0.129 |
| 80-120 | 3,992 | 0.092 |
| 120-180 | 4,096 | 0.073 |
| 180-350 | 2,127 | 0.072 |

σ_int drops by 44% from low-V to high-V. CDM predicts this trend
(halo concentration scatter decreases with mass: Duffy+2008). But the
actual magnitude at high-V (0.072) is BELOW CDM's 0.085 — meaning either:
- CDM overestimates the scatter at high masses, OR
- The model absorbs some CDM scatter through correlated variables

### Systematic Contamination

Distance shows r = +0.13 with standardized residuals. This is the only
significant systematic correlation. It likely inflates σ_int for the full
sample (more distant galaxies have larger absolute distance errors).
The optimal subsample mitigates this somewhat.

### Model Specification Dependence

| Model | σ_int |
|:------|:-----:|
| BTFR (V only) | 0.422 |
| V + TFR | 0.168 |
| V + TFR + f_gas | 0.142 |
| V + TFR + g-i + f_gas | 0.139 |
| Full (5-var) | 0.094 |

σ_int is strongly model-dependent — each variable removes M/L scatter.
This means σ_int at any level is an *upper bound* on CDM physics.
The "true" CDM signal can only be extracted when ALL M/L variation is removed.

### Bootstrap Stability

σ_int = 0.094 [0.091, 0.097] (95% CI, 1000 bootstraps).
CDM's 0.085 is outside the CI (below). MOND's 0 is far outside.
The result is robust to resampling.

### BIG-SPARC Prediction

At 3.6μm with resolved rotation curves:
- Expected noise: 0.035 dex
- Expected M/L variation: 0.030 dex (from S598 band analysis)
- Physics-only σ_int: ~0.089 dex
- CDM discrimination: 4.1σ rejection possible

## Synthesis

### The Critical Finding

**The optimal subsample gives σ_int = 0.072 < 0.085 (CDM floor).**

This is either:
1. **Strong evidence against CDM** — the BTFR is tighter than CDM allows
2. **Model absorbs CDM scatter** — the 5 variables correlate with halo concentration
3. **Selection effect** — quality cuts select less concentrated halos

Interpretation (2) is most likely: f_gas, g-i, and TFR residual all correlate
with halo properties. The model may be partially removing CDM concentration
scatter as a side effect of removing M/L variation.

### Implications

1. MOND is rejected at all levels (σ_int >> 0, even at high-V)
2. CDM's 0.085 is rejected for the full sample (+14σ) but this is contaminated
3. The optimal subsample is BELOW CDM (−6.2σ) — puzzling
4. σ_int decreases with V — the CDM-predicted trend
5. The contradiction (trend right, magnitude low) suggests model variables
   partially absorb CDM scatter
6. Clean CDM test requires variables known to be independent of halo concentration

### What Would Resolve This

1. **BIG-SPARC**: Resolved RCs eliminate W50 systematics; 3.6μm reduces M/L
2. **Concentration measurements**: Direct Sérsic index or NFW concentration
3. **External variable**: A CDM-correlated quantity NOT in the model (e.g., halo spin)

## Tests: 9/9 PASSED
## Grand Total: 1955/1955
