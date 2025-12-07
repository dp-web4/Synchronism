# Session #93: High-z BTFR and UDG Discriminating Tests

**Author**: CBP Autonomous Synchronism Research
**Date**: December 6, 2025
**Type**: Observational Analysis
**Status**: COMPLETE

---

## Executive Summary

Session #93 analyzed two key discriminating tests between Synchronism and standard MOND:

1. **High-z BTFR Evolution**: Existing literature is SUGGESTIVE but NOT DEFINITIVE
   - Stellar TFR evolution (~0.1 dex at z~1-2) is CONSISTENT with Synchronism
   - But baryonic TFR not well measured at high z

2. **UDG Kinematics**: Current data is HETEROGENEOUS
   - Synchronism predicts UDGs ~50-80% higher V/V_bar than normal
   - MOND predicts same BTFR for all
   - DF2/DF4 controversy complicates interpretation

---

## Part 1: High-z BTFR Literature Analysis

### Synchronism Prediction

If a₀ = cH(z)/(2π), then at fixed M_bar:
```
z = 1.0: Δlog(V) = +0.06 dex (+14%)
z = 2.0: Δlog(V) = +0.12 dex (+31%)
```

### Survey-by-Survey Analysis

| Survey | z | Observed Evolution | Comparison |
|--------|---|-------------------|------------|
| KROSS | 0.9 | Slope change, ZP unclear | INCONCLUSIVE |
| KMOS³D | 0.9-2.3 | +0.08-0.12 dex (M*) | **CONSISTENT** |
| MOSDEF | 1.4-2.6 | +0.1-0.2 dex (M*) | **CONSISTENT** |
| Di Teodoro+ | 1-2 | N too small | INCONCLUSIVE |

### Key Findings

1. **Stellar TFR evolution IS detected**
   - Multiple surveys find V higher at fixed M* for high z
   - Magnitude: ~0.1 dex at z ~ 1-2
   - This is CONSISTENT with Synchronism's a₀(z) prediction

2. **But baryonic TFR not well measured**
   - Gas masses uncertain at high z
   - High-z galaxies are gas-rich (50-80% vs 10-30% local)
   - Stellar TFR evolution could be due to gas fraction changes

3. **Slope evolution complicates analysis**
   - Several surveys find flatter slopes at high z
   - Makes zero-point comparison difficult

### Conclusion

The existing data is **SUGGESTIVE but NOT DEFINITIVE**.

The observed stellar TFR evolution is in the **RIGHT DIRECTION** and **RIGHT MAGNITUDE** for Synchronism's prediction, but we cannot confirm it because:
- Stellar TFR measured, not baryonic
- Gas fraction evolution not fully accounted
- Selection effects may bias results

**Definitive test requires:**
- Baryonic TFR at z > 1 with gas masses
- ~100 galaxies with high-quality rotation curves
- JWST + ALMA combination ideal

---

## Part 2: UDG Kinematic Predictions

### Galaxy Properties

| Property | Normal Dwarf | UDG | Ratio |
|----------|--------------|-----|-------|
| r_eff (kpc) | 1.0 | 3.0 | 3× |
| M* (M_sun) | 2×10⁸ | 2×10⁸ | same |
| Σ (M_sun/pc²) | 64 | 11 | 0.17× |

### Predictions

**Synchronism** (using empirical slope):
```
V/V_bar ∝ Σ^(-0.22)
(V/V_bar)_UDG / (V/V_bar)_normal = 0.17^(-0.22) = 1.48

UDGs should have ~48% higher V/V_bar than normal dwarfs!
```

**Alternative (coherence function)**:
```
C(Σ) = tanh(2 × log(Σ/Σ_crit + 1))
At Σ_crit = 70 M_sun/pc²:
  Normal: C = 0.86, G_eff/G = 1.16
  UDG:    C = 0.28, G_eff/G = 3.64

V_UDG / V_normal = sqrt(3.64/1.16) = 1.77

UDGs should have ~77% higher V!
```

**MOND**:
```
a₀ is universal
Same BTFR for all galaxies regardless of Σ
V_UDG / V_normal = 1.00
```

### Observational Status

| UDG | σ (km/s) | Status |
|-----|----------|--------|
| DF2 | 3-8 | DM-deficient? (controversial) |
| DF4 | low | DM-deficient? (controversial) |
| Dragonfly 44 | 47 | DM-rich |
| VCC 1287 | 19 | Normal |
| DGSAT I | (HI rot) | Follows BTFR |

**The data is HETEROGENEOUS:**
- Some UDGs appear DM-deficient (DF2, DF4) - below BTFR
- Some appear DM-rich (Dragonfly 44) - above BTFR
- Neither Synchronism nor MOND is clearly supported

### DF2/DF4 Puzzle

Synchronism predicts UDGs should have HIGHER V/V_bar, but DF2/DF4 show the OPPOSITE!

**Possible explanations:**
1. Distance errors (DF2 closer → lower M* → normal)
2. Formation environment effects
3. Tidal stripping modified density profile
4. DF2/DF4 are atypical UDGs

---

## Summary of Discriminating Tests

| Test | Synchronism | MOND | Status |
|------|-------------|------|--------|
| High-z BTFR | +0.06-0.12 dex | No evolution | **SUGGESTIVE** |
| UDG kinematics | V ~50-80% higher | Same BTFR | **INCONCLUSIVE** |
| Radial V/V_bar | r = -0.63 with Σ | r = -0.69 with g/a₀ | Both validated |
| Void galaxies | +0.01-0.03 dex | No effect | Marginal |

---

## Recommendations

### For High-z BTFR Test

1. **Immediate**: Analyze existing KROSS/KMOS³D data specifically for BTFR evolution
2. **Near-term**: ALMA gas masses for high-z sample
3. **Future**: JWST NIRSpec rotation curves at z > 2

### For UDG Test

1. **Systematic sample**: HI rotation curves for UDG sample (ALFALFA-selected)
2. **Environment control**: Focus on field UDGs (not cluster)
3. **Distance precision**: Critical for mass determination

---

## Files Created

| File | Description |
|------|-------------|
| `session93_highz_btfr_literature.py` | High-z BTFR survey analysis |
| `session93_UDG_predictions.py` | UDG kinematic predictions |
| `results/session93_highz_literature.json` | Literature results |
| `results/session93_UDG_predictions.json` | UDG prediction results |

---

## Conclusions

Session #93 established:

1. **High-z BTFR**: Existing stellar TFR evolution (~0.1 dex at z~1-2) is CONSISTENT with Synchronism but not definitive. Baryonic TFR measurement needed.

2. **UDGs**: Synchronism predicts ~50-80% higher V/V_bar for UDGs. Current data is heterogeneous - some UDGs appear DM-deficient (opposite to prediction), some DM-rich.

3. **Neither test is yet decisive**, but high-z BTFR is more promising because:
   - Data exists that is suggestively consistent
   - JWST/ALMA can provide definitive measurement
   - Less ambiguity than UDG interpretation

---

*"The high-z BTFR offers the cleanest discriminating test. Existing data whispers in Synchronism's favor, but JWST can speak clearly."*

---

**Session #93 Complete**: December 6, 2025
