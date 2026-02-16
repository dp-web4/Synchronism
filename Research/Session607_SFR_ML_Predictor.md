# Session #607: SFR as M/L Predictor — Does Star Formation Break the Barrier?

**Date**: 2026-02-16
**Grade**: A
**Domain**: Cosmology / M/L Prediction

## Objective

Session #606 achieved σ_int = 0.094 (full) and 0.072 (optimal subsample)
using 5 variables. The Durbala+2020 catalog contains additional columns that
were never downloaded: three independent SFR estimates, a third stellar mass,
internal extinction corrections, plus W20 (profile width) from Haynes.

Can specific star formation rate (sSFR) reduce BTFR scatter beyond TFR + g-i
+ f_gas + Δ_SPS?

## New Data Downloaded

| Column | Source | N (quality) | Coverage |
|:-------|:-------|:-----------:|:--------:|
| logSFR22 | 22μm unWISE | 11,648 | 80.7% |
| logSFRN | GALEX NUV | 7,200 | 49.9% |
| logSFRG | GSWLC-2 | 6,895 | 47.8% |
| logMsG | GSWLC-2 | 6,895 | 47.8% |
| W20 | Haynes α.100 | 14,426 | 99.9% |
| Ag, Ai | Durbala internal ext. | 14,435 | 100% |

## Key Results

### sSFR22 Does NOT Correlate with BTFR Residuals

| SFR indicator | r(sSFR, BTFR_resid) | N | p |
|:-------------|:--------------------:|:-:|:-:|
| 22μm (WISE) | **-0.012** | 11,648 | **0.21** |
| NUV (GALEX) | -0.207 | 7,200 | 3e-70 |
| GSWLC-2 | -0.192 | 6,895 | 5e-58 |

**The 22μm SFR is NOT significantly correlated with BTFR residuals** (p=0.21).
This is surprising — 22μm traces obscured star formation via warm dust heated
by young stars. It should predict M/L.

NUV and GSWLC SFRs DO correlate strongly (r ≈ -0.2), consistent with the
expectation that high sSFR → low M/L → negative BTFR residual.

The discrepancy occurs because 22μm SFR and optical M/L are nearly
orthogonal: 22μm traces embedded star formation that optical colors already
miss. The BTFR uses optical stellar masses (from Taylor/Mendel), so
dust-obscured SFR doesn't affect the measured M/L.

### All SFR Indicators Add ~0% Beyond TFR

| Model | σ (dex) | Gain beyond TFR |
|:------|:-------:|:---------------:|
| BTFR | 0.3898 | — |
| + TFR | 0.1848 | — |
| + TFR + sSFR22 | 0.1838 | 0.6% |
| + TFR + g-i | 0.1843 | 0.3% |
| + TFR + f_gas | 0.1645 | 11.0% |
| + TFR + sSFR + f_gas | 0.1645 | 11.0% |

**sSFR adds 0.6% beyond TFR — effectively zero.** This confirms S594's finding
that TFR absorbs all M/L information from color and SFR. The only variable that
consistently adds beyond TFR is f_gas (11%).

### All Three SFRs Are Equivalent

On the common 3-SFR subsample (N=5,502), all three add ~0% beyond TFR.
Inter-correlations: r(22μm, NUV) = 0.735, r(NUV, GSWLC) = 0.910.

### W20/W50 Profile Shape: 0% Gain

Despite r(W20/W50, BTFR_resid) = +0.153 (highly significant), W20/W50 adds
exactly 0% beyond TFR. The TFR residual already captures the profile shape
information.

W20/W50 median = 1.119, strongly anti-correlated with logV (r = -0.691).
This confirms that low-mass galaxies have broader, more Gaussian HI profiles
(rising RCs) while massive galaxies have flat-topped profiles (flat RCs).

### Third Mass (GSWLC-2) is Worst

| Method | BTFR σ | TFR-corrected σ |
|:------:|:------:|:---------------:|
| Taylor | 0.329 | **0.144** |
| Mendel | 0.354 | **0.122** |
| GSWLC | 0.360 | 0.151 |

GSWLC masses give the worst raw AND worst TFR-corrected BTFR.
Mendel remains best post-correction by a wide margin.

### Kitchen Sink: 8-Variable Model (Taylor masses)

| Variable | β | t-stat |
|:---------|:-:|:------:|
| logV | +1.707 | +108.7 |
| TFR | +0.768 | +112.3 |
| Δ_SPS | -0.105 | -10.5 |
| g-i | +0.661 | +72.1 |
| f_gas | +1.069 | +100.2 |
| **sSFR22** | **+0.042** | **+15.7** |
| W20/W50 | +0.026 | +2.9 |
| **Ag** | **+0.152** | **+14.6** |

σ = 0.133 dex. Gain beyond 5-var: only 2.1%.

**Paradox**: sSFR22 is NOT significant raw (r=-0.012, p=0.21) but IS
significant at 15.7σ in the kitchen-sink model! This means sSFR22 provides
information that is ONLY useful in combination with the other 5 variables —
a suppressor effect. The 22μm SFR captures dust-embedded star formation
orthogonal to the optical variables.

Ag (dust extinction) also significant at 14.6σ — internal dust independently
predicts M/L.

### CDM Discrimination: sSFR Deepens the Puzzle

| Model | σ_int | z(CDM) |
|:------|:-----:|:------:|
| 5-var (Mendel, full) | 0.086 ± 0.001 | +1.1 |
| **6-var (+sSFR22)** | **0.081 ± 0.001** | **-7.1** |
| 5-var (optimal, N=600) | 0.064 ± 0.002 | -10.2 |
| **6-var (optimal)** | **0.061 ± 0.002** | **-12.2** |

Adding sSFR22 reduces σ_int by 5.8% and pushes the FULL sample below CDM's
0.085 floor at 7.1σ. The optimal subsample reaches 0.061, now 12.2σ below CDM.

**This strengthens the S606 finding.** Either:
1. CDM's concentration scatter is overestimated
2. The model absorbs CDM scatter through correlated variables
3. Some systematic reduces σ_int below reality

### sSFR Is Partially Independent of Color

r(logsSFR22, g-i) = -0.411. Moderate anti-correlation (blue = high sSFR),
but 83% of sSFR variance is independent of color.

Partial r(sSFR, Mbar | V, TFR, g-i, f_gas) = +0.065 (p = 1.3e-6).
The partial sign is POSITIVE (opposite to raw) — another suppressor signature.

## Synthesis

### The Five-Word Summary
**SFR adds little; TFR dominates.**

### What We Learned
1. sSFR22 (22μm) is NOT correlated with raw BTFR residuals (r=-0.012, p=0.21)
2. NUV and GSWLC sSFR DO correlate (r ≈ -0.2) but add 0% beyond TFR
3. W20/W50 profile shape adds 0% beyond TFR
4. GSWLC stellar masses are worst of three SPS methods
5. Kitchen-sink (8-var) gains only 2.1% beyond 5-var
6. sSFR22 IS significant at 15.7σ in the full model — suppressor effect
7. Ag (extinction) significant at 14.6σ — dust matters
8. Adding sSFR pushes σ_int to 0.081, now 7.1σ BELOW CDM for full sample
9. Optimal subsample: σ_int = 0.061, 12.2σ below CDM
10. TFR residual remains the dominant M/L predictor by far

### Implications

The TFR is an extraordinary information compressor. It captures M/L information
from color, SFR, profile shape, and extinction — all in a single V-L residual.
No additional variable provides more than 2% improvement beyond TFR + f_gas.

The deepening CDM puzzle (now 12.2σ below the 0.085 floor) increasingly
suggests interpretation (2) from S606: the model variables correlate with
halo concentration, partially absorbing CDM scatter. sSFR likely correlates
with halo formation time (which determines concentration), explaining why
adding it pushes σ_int further below CDM.

This is important: **the model may be REMOVING CDM scatter rather than failing
to detect it.** If so, CDM is not rejected — but the BTFR test is compromised
by variable selection. A clean CDM test requires variables known to be
independent of halo concentration.

## Tests: 9/9 PASSED
## Grand Total: 1964/1964
