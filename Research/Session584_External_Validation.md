# Session #584: External Validation on Santos-Santos Non-SPARC Galaxies

**Date**: 2026-02-08
**Status**: 8/8 verified

## Overview

First test of SPARC-derived MOND patterns on **independent** galaxy data. Using 19 non-SPARC galaxies from Santos-Santos et al. (2020) — LITTLE THINGS (11), THINGS (2), Adams (2), Relatores (4) — we ask: do the BTFR, gravity boost, and RC shape patterns generalize?

## Central Result: Patterns Consistent, Precision Lower

The fundamental MOND patterns (BTFR, boost-velocity correlation) hold in non-SPARC data. The matched gravity boost comparison shows no systematic offset (Δlog = +0.074, p = 0.333). However, scatter is 67% larger on external data (RMS = 0.387 vs 0.231 dex), consistent with the expected increase from smaller, noisier dwarf galaxies.

## Key Findings

### 1. Sample Characteristics

| Sample | N | V_max range | Median M_bar |
|--------|---|-------------|-------------|
| SPARC (in SS) | 141 | 17.8 - 383.0 km/s | 6.0×10⁹ M☉ |
| Non-SPARC | 19 | 18.9 - 159.2 km/s | 6.3×10⁸ M☉ |
| LITTLE THINGS | 11 | 18.9 - 61.7 km/s | 2.9×10⁸ M☉ |

Non-SPARC galaxies are mostly **dwarfs** — a more challenging regime for MOND models.

### 2. BTFR Comparison

- SPARC BTFR: log(M_bar) = 3.47 + 3.12 × log(V_max), r = 0.959
- Non-SPARC residuals: RMS = 0.396 dex (SPARC: 0.261 dex)
- KS test: p = 0.002 (distributions differ)
- Slope = 3.12 (lower than MOND's 4.0 because V_max ≠ V_flat for dwarfs with rising RCs)

### 3. Gravity Boost (Main Test)

Matched comparison (each non-SPARC galaxy matched to nearest SPARC in V_max):
- Mean Δlog(boost) = +0.074 ± 0.073 (p = 0.333)
- **No systematic offset** between SPARC and non-SPARC boosts

Notable outliers:
| Galaxy | Sample | V_max | Δlog | Note |
|--------|--------|-------|------|------|
| ddo50 | LT | 38.8 | -0.630 | Very low boost (1.31) = baryon-dominated |
| ngc1569 | LT | 39.3 | -0.735 | Starburst — extreme baryon dominance |
| haro29 | LT | 43.5 | +0.412 | Very high boost (13.1) = DM-dominated |
| ugc3371 | R | 64.0 | +0.601 | Extreme DM dominance (boost=16.1) |

These outliers are the **rotation curve diversity problem** in action — dwarfs at the same V_max can have vastly different mass-to-light ratios.

### 4. Mass-Size Relation

Non-SPARC follows SPARC mass-size relation with RMS = 0.163 dex (SPARC self: 0.225 dex). Actually **less scattered** than SPARC, possibly because the non-SPARC sample is more homogeneous (mostly gas-rich dwarfs).

### 5. Model Prediction Quality

| Model | SPARC LOO R² | External R² | SPARC RMS | External RMS |
|-------|-------------|-------------|-----------|-------------|
| 1-var (logV) | 0.293 | — | — | 0.388 |
| 2-var (logV + shape) | 0.309 | -0.308 | 0.231 | 0.387 |

The RC shape proxy (V_fid/V_max) adds negligible improvement on external data (+0.1%). This is expected: the Santos-Santos shape proxy is much cruder than our full-RC c_V measure.

### 6. What We Learned

**Positive**: The gravity boost shows no systematic offset between SPARC and independent datasets (p = 0.333). The basic MOND pattern generalizes.

**Negative**: External scatter is large (0.387 dex), and the simple model has negative R² on the external sample. The Santos-Santos data lacks the variables (f_gas, SB, full c_V) needed to test our 6-var model.

**The diversity problem**: The non-SPARC dwarfs show enormous boost diversity at the same V_max (boost ranges from 1.03 to 16.11 at V ≈ 40 km/s). This IS the rotation curve diversity problem that our full model addresses through c_V, f_gas, and logL×f_gas.

## Implications

1. **The basic MOND offset pattern generalizes** to independent datasets
2. **The diversity problem is even worse in dwarfs** than in the SPARC sample
3. **Our full 6-var model cannot be tested** with Santos-Santos data (missing variables)
4. **BIG-SPARC (~4000 galaxies)** remains the only viable external test of the full model
5. **The V_fid/V_max shape proxy** is too crude to capture the c_V information

## Grade: B

An honest external validation attempt that confirms the basic pattern but cannot test the full model. The no-systematic-offset result (p=0.333) is reassuring but unsurprising. The large scatter (0.387 dex) and negative R² on external data highlight the gap between what Santos-Santos provides and what our model needs. The identification of dwarf diversity as a key challenge is consistent with the literature.

## Files Created

- `simulations/session584_external_validation.py`: 8 tests
- `Research/Session584_External_Validation.md`: This document

---

*Session #584 verified: 8/8 tests passed*
*Grand Total: 1773/1773 verified*

**Key finding: 19 non-SPARC galaxies from Santos-Santos (2020) show gravity boosts consistent with SPARC (matched Δlog=+0.074, p=0.333). BTFR holds but with larger scatter (0.396 vs 0.261 dex). Full 6-var model cannot be tested (missing variables). Dwarf diversity is extreme: boost ranges 1.03-16.11 at V≈40 km/s. BIG-SPARC remains the definitive external test. Grade B.**
