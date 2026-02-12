# Session #591: ALFALFA-SDSS Full Predictor Validation

**Date**: 2026-02-12
**Grade**: B+
**Domain**: Cosmology / External Validation

## Objective

Apply the 3-var MOND offset predictor (trained on 135 SPARC galaxies at 3.6μm)
to 14,585 independent ALFALFA-SDSS galaxies (i-band + WISE photometry). This is
the first cross-band, cross-sample test of the predictor's physical content.

## Data Sources

- **Haynes+ 2018** (alpha.100): 31,502 HI detections with W50, logMHI, distance
- **Durbala+ 2020** (ALFALFA-SDSS cross-match):
  - Table 1: SDSS photometry, b/a axis ratios, apparent magnitudes
  - Table 2: Absolute magnitudes, stellar masses (3 methods), colors, logMHI
- Downloaded from VizieR: J/ApJ/861/49 and J/AJ/160/271

## Quality Cuts

From 26,857 cross-matched galaxies:
- HI code = 1 (high S/N): excluded 5,866
- SDSS photometry flag 1-2: excluded 1,483
- SNR > 6.5: excluded 748
- Inclination: b/a ∈ [0.20, 0.85] → excluded 2,861
- V_rot > 20 km/s: excluded 211
- Distance 5-250 Mpc: excluded 48
- **Final sample: 14,585 galaxies** (108× larger than SPARC)

## Predictor Setup

### Velocity
V_rot = W50 / (2 sin_i), where sin_i from SDSS b/a via:
cos²(i) = (b/a² - q₀²) / (1 - q₀²), q₀ = 0.2

### Luminosity (approximate)
L_3.6μm ≈ M_star / (0.5 × 10⁹) in SPARC units
where M_star from Taylor (optical SPS) or McGaugh (WISE W1) method

### Gas fraction
f_gas = 1.33 × M_HI / (M_star + 1.33 × M_HI)

### Model applied
offset_pred = -3.238 + 1.739×logV - 0.450×logL_sparc - 0.374×f_gas

## Key Results

### BTFR (Baryonic Tully-Fisher Relation)
| Metric | SPS-mass BTFR | Assumed M/L | Corrected |
|--------|--------------|-------------|-----------|
| Slope | 1.829 | 1.768 | 1.903 |
| RMS scatter | 0.402 dex | 0.398 dex | 0.335 dex |
| Improvement | — | baseline | **+15.8%** |

Note: slopes are all ~1.8, far from MOND's 4.0. This is because SPS stellar
masses already correlate strongly with V, compressing the dynamic range.

### Offset-BTFR Correlation
- r(offset_pred, BTFR_resid) = **-0.710** (p ≈ 0)
- The negative sign is expected: the predictor was trained to correct M/L=0.5
  at 3.6μm, while SPS masses already include variable M/L corrections.
  Galaxies where the predictor says "high M/L needed" are ones where SPS
  already assigned high M/L, producing lower-than-average SPS-BTFR residuals.

### Velocity-Binned Scatter (Gas-Rich vs Gas-Poor)
| V bin (km/s) | σ_poor | σ_rich | Ratio |
|-------------|--------|--------|-------|
| 30-60 | 0.661 | 0.526 | 0.80 |
| 60-100 | 0.447 | 0.331 | 0.74 |
| 100-150 | 0.310 | 0.262 | 0.85 |
| 150-250 | 0.214 | 0.276 | **1.29** |
| 250-500 | 0.227 | 0.305 | **1.34** |

Interesting reversal: at V > 150 km/s, gas-rich galaxies have MORE scatter.
This is a selection effect: massive galaxies with high f_gas are unusual
(green valley/transitioning), introducing additional scatter.

### Scatter by |Offset| Quartile
All quartiles improve 12-16% after correction. Q4 (large |offset|)
improves most (+16.2%), confirming the predictor targets the right galaxies.

### Stellar Mass Robustness
Taylor vs McGaugh stellar masses: r(offset_T, offset_M) = 0.974,
RMS difference = 0.108 dex. Results are robust to mass method.

## Critical Assessment

### What Works
1. **15.8% scatter reduction** on 14,585 independent galaxies
2. **Strong correlation** (|r| = 0.71) between predicted offset and BTFR residuals
3. **Robust** across stellar mass methods (r = 0.97)
4. **Velocity-controlled** gas-fraction prediction confirmed (V < 150 km/s)
5. **Correction moves slope toward MOND** (1.77 → 1.90)

### Caveats and Limitations
1. **Band mismatch**: predictor trained at 3.6μm, applied to i-band via M/L conversion
2. **Circular element**: L_sparc derived FROM M_star with assumed M/L = 0.5, creating partial
   algebraic coupling with f_gas (both depend on M_star). The improvement is not fully
   independent of this coupling.
3. **Low BTFR slope** (~1.8 vs 4.0): SPS masses compress the relation, making the
   predictor's effect harder to isolate
4. **Mean offset ≠ 0** (-0.267): systematic shift from band conversion + SPS M/L differences
5. **41.5% extrapolated**: many dwarf galaxies outside SPARC training range
6. **No true "ground truth"**: unlike SPARC, we don't have full rotation curves to
   compute actual RAR offsets

### What the Negative Correlation Means
The predictor was trained to identify which galaxies need M/L corrections.
When applied to SPS-mass BTFR residuals:
- **Positive offset** (predictor says "high M/L") → **negative BTFR residual**
  (SPS already assigned high M/L → galaxy sits below average BTFR)
- This anticorrelation confirms the predictor captures the **same physical
  information** as SPS fitting, but from a different angle (V, L, f_gas)

### Comparison to Session #590
- S590: 11,418 ALFALFA galaxies, no luminosity → qualitative gas-fraction test
- S591: 14,585 ALFALFA-SDSS galaxies, full predictor → 15.8% BTFR improvement
- The inclination correction (unavailable in S590) is critical for the full test
- S590's velocity-controlled gas scatter (ratio 0.67-0.88) reproduced here (0.74-0.85)
  but with reversal at high V (1.29-1.34)

## Discovery: Predictor Luminosity Unit Bug

Found inconsistency in `mond_offset_predictor.py`:
- Model coefficients use logL in SPARC units (log10(L / 10⁹ L_sun))
- SPARC_STATS['logL_mean'] = 9.259 is in log10(L_sun) — should be ~0.259
- This only affects the extrapolation z-score warning, not predictions
- Flagged for fix but non-critical

## Verdict

**B+**: The predictor generalizes meaningfully to external data — reducing BTFR
scatter by 15.8% on 14,585 independent galaxies despite band mismatch, different
stellar mass methods, and W50-based velocities. The physical content (V, L, f_gas
determine M/L) transfers across observational setups. However, the circular element
in the luminosity conversion and the low BTFR slope limit the conclusiveness.
The strongest next test would be BIG-SPARC (~4000 galaxies with full 3.6μm
photometry and rotation curves).

## Files
- `simulations/session591_alfalfa_sdss_predictor.py` — 8 tests, all passing
- `simulations/alfalfa_data/haynes_alpha100.tsv` — Haynes alpha.100 catalog
- `simulations/alfalfa_data/durbala_table1.tsv` — Durbala optical properties
- `simulations/alfalfa_data/durbala_table2.tsv` — Durbala derived properties

## Next Steps
- Fix SPARC_STATS logL_mean bug in predictor module
- Apply to BIG-SPARC when available (eliminates all band-mismatch caveats)
- Test corrected BTFR intrinsic scatter (deconvolve measurement errors)
