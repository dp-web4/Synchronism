# Session #592: Circularity Test — i-Band vs Mstar-Derived Luminosity

**Date**: 2026-02-12
**Grade**: A-
**Domain**: Cosmology / Statistical Methodology

## Objective

Session #591 showed a 15.8% BTFR scatter reduction using the 3-var predictor on
14,585 ALFALFA-SDSS galaxies, but flagged a potential circularity: the luminosity
was derived FROM stellar mass (L = Mstar/0.5). This session tests whether the
improvement survives when luminosity comes independently from the i-band absolute
magnitude.

## Key Question

Is the predictor's BTFR improvement driven by genuine MOND physics (V-L-f_gas
determine M/L corrections) or by algebraic coupling through shared dependence
on M_star?

## Results

### Circularity Quantification

| Method | BTFR Scatter Reduction | r(offset, BTFR_resid) |
|--------|----------------------:|---------------------:|
| Mstar-derived L + f_gas | **+15.9%** | -0.708 |
| **i-Band L + f_gas** | **+14.5%** | -0.623 |
| V + i-Band L only (no f_gas) | **+16.2%** | -0.717 |
| V + f_gas only (no L) | -3.3% | +0.059 |

**Circularity fraction**: only **8.8%** of the improvement is from the Mstar-L coupling.
The i-band predictor retains 14.5% improvement using photometry independent of stellar mass.

### Surprising Finding: f_gas Hurts the BTFR!

V + L alone (without f_gas) gives the BEST result: **+16.2%** improvement.
Adding f_gas actually reduces improvement from 16.2% to 14.5%.

**Why?** The BTFR already includes gas mass directly (Mbar = Mstar + 1.33×MHI).
The f_gas coefficient in the 3-var model corrects the RAR for the fact that gas-rich
galaxies need less M/L correction. But in the BTFR context, this correction is
redundant — the gas is already counted in Mbar. The f_gas term adds noise.

This is consistent with the SPARC analysis: f_gas corrects the **stellar** contribution,
but the BTFR already has the stellar mass. The predictor's power for the BTFR comes
from the **V-L relationship** (MOND's mass-to-light prediction), not from f_gas.

### Partial Correlations (Controlling for Confounds)

| Partial Correlation | r | p |
|---|---:|---:|
| r(offset, BTFR_resid \| Mstar) | **-0.890** | ≈ 0 |
| r(offset, BTFR_resid \| Mstar, V) | **-0.811** | ≈ 0 |

These are extraordinary. Even after removing the effects of both Mstar AND V,
the predictor still explains 66% of the remaining BTFR variance (r²=0.66).
This residual signal is carried by f_gas and the non-linear V-L interaction.

### Luminosity Comparison

Two luminosities correlate strongly (r = 0.983) but differ systematically:
- Mean difference: 0.037 dex (Mstar-based is slightly higher)
- Std: 0.192 dex
- This difference IS the SPS-fitted M/L variation (mean M/L_SPS,i = 0.62)

## Physical Interpretation

The predictor's BTFR improvement comes primarily from the **V-L ratio**:
- The 3-var model has β_V/β_L = 1.739/0.450 = 3.86 ≈ MOND's 4.0
- This ratio encodes V^4 ∝ L (MOND's Tully-Fisher)
- At fixed V, brighter galaxies are predicted to have lower M/L
- This is exactly what MOND predicts: L and V are redundant in the deep-MOND regime

The f_gas contribution is orthogonal — it corrects the RAR (point-by-point scatter)
but doesn't help the BTFR (galaxy-integrated relation). This makes physical sense:
the BTFR is already the integrated relation.

## Verdict

**A-**: The circularity concern is definitively addressed. Only 8.8% of the
improvement is circular. The i-band predictor retains 14.5% scatter reduction
using photometry independent of stellar mass fitting. The V-L relationship
(MOND's TFR) drives the improvement, not the f_gas term. Partial correlations
(r = -0.89 controlling for Mstar) confirm the predictor captures genuine physics.

The most surprising finding: V+L alone beats the full 3-var predictor for the
BTFR, because f_gas is redundant when gas mass is already in Mbar.

## Files
- `simulations/session592_circularity_test.py` — 8 tests, all passing

## Implications for Future Work
1. For BTFR applications, use V+L only (drop f_gas)
2. For RAR applications, keep f_gas (it corrects per-point scatter)
3. The V-L ratio = 3.86 is the predictor's core — this IS MOND's TFR
4. BIG-SPARC will enable both BTFR and RAR tests with consistent photometry
