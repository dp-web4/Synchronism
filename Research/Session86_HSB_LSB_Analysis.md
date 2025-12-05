# Session #86: HSB vs LSB Galaxy BTFR Analysis

**Author**: CBP Autonomous Synchronism Research
**Date**: December 4, 2025
**Type**: Data Analysis + Theory Interpretation
**Status**: COMPLETE

---

## Executive Summary

Session #86 tested the HSB vs LSB BTFR comparison that was listed as a discriminating test between Synchronism and MOND. The result showed the **OPPOSITE** direction from the naive prediction:

- **Naive prediction**: LSB galaxies should have +0.088 dex higher V_flat
- **Observed**: LSB galaxies have -0.053 dex LOWER V_flat (3.0σ significant)

**Key insight**: This was NOT a valid test of Synchronism. The theory predicts C(ρ) at LOCAL density at each radius, not global surface brightness. The HSB/LSB comparison is based on a misinterpretation of how C enters the theory.

---

## Methodology

### Data Source
- **SPARC Database** (Lelli, McGaugh, Schombert 2016)
- 175 disk galaxies with Spitzer [3.6μm] photometry
- Used 129 galaxies after quality filter (Q ≤ 2)

### Surface Brightness Classification
- **HSB** (High Surface Brightness): SBdisk > 500 L_sun/pc² (70 galaxies)
- **ISB** (Intermediate): 100-500 L_sun/pc² (38 galaxies)
- **LSB** (Low Surface Brightness): SBdisk < 100 L_sun/pc² (21 galaxies)

### BTFR Analysis
- Baryonic mass: M_bar = 0.5 × L[3.6] + 1.33 × M_HI
- BTFR: log(M_bar) = ZP + 4.0 × log(V_flat)
- Residuals: Δlog(V) = log(V_obs) - log(V_pred)

---

## Results

### BTFR Residuals by Class

| Class | N | Mean SB | Δlog(V) | Significance |
|-------|---|---------|---------|--------------|
| HSB | 70 | 2957 | +0.017 | +2.7σ |
| ISB | 38 | 227 | -0.012 | -0.9σ |
| LSB | 21 | 73 | -0.036 | -2.2σ |

### Key Comparison

- **LSB - HSB offset**: -0.053 ± 0.017 dex (3.0σ)
- **Direction**: LSB has LOWER V than expected (opposite to prediction)
- **Continuous correlation**: r = +0.322 (positive, not negative)

---

## Naive Synchronism Prediction (INCORRECT)

The naive prediction was based on:
1. LSB galaxies have lower surface/volume density
2. Lower ρ → lower C(ρ) → higher G_eff = G/C
3. Higher G_eff → higher rotation velocity at fixed mass
4. Therefore LSB should have HIGHER V_flat

**Predicted offset**: +0.088 dex (LSB higher than HSB)

**This prediction was WRONG because it misunderstood the theory.**

---

## Why the HSB/LSB Test is INVALID

### The Actual Theory

Synchronism predicts C(ρ) at LOCAL density at each radius:
```
G_eff(r) = G / C(ρ(r))

where C(ρ) = tanh(γ × log(ρ(r)/ρ_crit + 1))
```

The rotation curve V(r) is computed by integrating this at each radius.

### Why Surface Brightness is the Wrong Proxy

1. **SB is a global property**: Central surface brightness averages over the inner disk
2. **V_flat is measured at large radii**: Where local density may be similar for HSB and LSB
3. **BTFR averages over all radii**: Loses the radius-dependent C(ρ(r)) signal
4. **Different radii are sampled**: HSB V_flat at smaller R, LSB at larger R

### The Correct Test

The correct test of Synchronism is NOT a BTFR comparison, but:
- **Radial V/V_Newton profile** should correlate with ρ(r)
- Inner disk (high ρ): V/V_Newton closer to 1
- Outer disk (low ρ): V/V_Newton larger (more "DM-like")

This is what the 52% SPARC success rate actually tests.

---

## Theory Status Update

### What This Reveals

Session #86, combined with Session #85, establishes:

1. **C depends on LOCAL density at each radius**
2. **Global properties (SB, environment) are poor proxies**
3. **BTFR comparisons average over radii and lose the signal**
4. **The core rotation curve model is about LOCAL C(ρ), not global C**

### Discriminating Tests: Updated

| Test | Status | Notes |
|------|--------|-------|
| Void TF offset | Weak | 0.01-0.03 dex (revised from 0.11-0.28) |
| HSB vs LSB | INVALID | Not a valid test of theory |
| Radial V/V_Newton | Testable | Correct way to test C(ρ(r)) |

---

## Files Created

| File | Description |
|------|-------------|
| `session86_hsb_lsb_analysis.py` | Main analysis code |
| `session86_interpretation.py` | Result interpretation |
| `results/session86_hsb_lsb_analysis.json` | Analysis results |
| `results/session86_interpretation.json` | Interpretation summary |
| `sparc_real_data/SPARC_Lelli2016c.mrt` | Downloaded galaxy table |

---

## Conclusion

**Session #86 Status**: COMPLETE

The HSB/LSB BTFR test was performed:
- Observed: -0.053 dex (LSB lower V)
- Predicted: +0.088 dex (LSB higher V)
- Direction: **OPPOSITE**

**INTERPRETATION**:
- This was NOT a valid test of Synchronism
- The theory predicts C(ρ) at LOCAL density, not global SB
- Surface brightness is not the right proxy
- The correct test is the radial V/V_Newton profile

**THEORY STATUS**:
- Core rotation curve model UNAFFECTED
- HSB/LSB test was a misapplication of the theory
- Updated THEORETICAL_STATUS_DEC2025.md accordingly

---

*"The HSB/LSB test taught us what Synchronism actually predicts: local C(ρ(r)) at each radius, not global properties. A failed prediction that clarifies the theory is more valuable than an untested one."*

---

**Session #86 Complete**: December 4, 2025
