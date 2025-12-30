# Session #201: Precision a₀ Measurement Analysis

**Date**: December 30, 2025
**Machine**: CBP
**Status**: MAJOR INSIGHT DISCOVERED

---

## Executive Summary

Session #201 analyzed the precision measurement of the critical acceleration a₀ as a test distinguishing Synchronism from MOND.

**Key Discovery**: The difference between Synchronism and MOND goes beyond the 12% a₀ value difference. The **interpolating functions themselves differ qualitatively** in the deep MOND regime.

---

## The Two a₀ Values

| Framework | a₀ | Source |
|-----------|-----|--------|
| Synchronism | 1.05 × 10⁻¹⁰ m/s² | Derived: c H₀ Ω_m^φ |
| MOND | 1.2 × 10⁻¹⁰ m/s² | Empirical fit to BTFR |

**Ratio**: 0.874 (12.5% difference)

---

## Interpolating Function Comparison

### The Functions

**Synchronism**:
```
C(a) = Ω_m + (1 - Ω_m) × (a/a₀)^(1/φ) / [1 + (a/a₀)^(1/φ)]
G_eff/G = 1/C(a)
```

**MOND (standard)**:
```
ν(y) = 1/2 + √(1/4 + 1/y)  where y = a_N/a₀
a = a_N × ν(y)
```

### Enhancement Factor Comparison

| a/a₀ | Sync G_eff/G | MOND ν | Ratio |
|------|--------------|--------|-------|
| 0.01 | 2.84 | 10.5 | 0.27 |
| 0.10 | 2.23 | 3.70 | 0.60 |
| 0.50 | 1.71 | 2.00 | 0.85 |
| 1.00 | 1.52 | 1.62 | 0.94 |
| 2.00 | 1.37 | 1.37 | 1.00 |
| 10.0 | 1.15 | 1.09 | 1.06 |

### Deep MOND Limit (CRITICAL!)

**Synchronism**: G_eff/G → 1/Ω_m ≈ **3.17** (bounded)

**MOND**: ν → √(a₀/a_N) → **∞** as a_N → 0 (unbounded)

This is a **qualitative difference**, not just quantitative!

---

## Implications for Testing

### For Typical Galaxies (a ~ a₀)

At a/a₀ ~ 1, both frameworks give similar enhancement (~1.5-2×).
The 12% a₀ difference translates to ~3% velocity difference.
This is **marginally detectable** with current data (~5-10% errors).

### For Deep MOND Systems (a << a₀)

This is where the frameworks diverge dramatically:
- MOND predicts arbitrarily large enhancement
- Synchronism caps enhancement at 3.17×

**Test systems**:
- Ultra-faint dwarfs (UFDs)
- Outer rotation curves of gas-rich dwarfs
- Extended HI disks

**Challenge**: These systems also have the most indifferent mass (in Synchronism interpretation), complicating the comparison.

---

## Rotation Curve Analysis

### Example: NGC 1560-like Galaxy

| Model | V at 5 kpc | Enhancement |
|-------|------------|-------------|
| Newtonian | 25 km/s | 1.0× |
| Synchronism | 38 km/s | 1.5× |
| MOND | 59 km/s | 2.4× |

The **21 km/s difference** between Synchronism and MOND is much larger than the ~5-10 km/s measurement errors.

### Interpretation

If observations show V ~ 60 km/s → Favors MOND
If observations show V ~ 40 km/s → Favors Synchronism (+ indifferent mass)

**Current MOND fits work well** with observed rotation curves. This suggests:
1. MOND interpolation is approximately correct, OR
2. Synchronism + significant indifferent mass mimics MOND

---

## The Indifferent Mass Complication

In deep MOND systems, Synchronism predicts:
- G_eff enhancement: max 3.17×
- Indifferent mass: must make up the difference

For UFDs with M_dyn/M_* ~ 100-1000:
- G_eff contributes: ~3×
- Indifferent mass: ~30-300×

This is **consistent** with early-universe accretion of indifferent patterns.

But it means **we cannot distinguish pure G_eff from G_eff + indifferent mass** in these systems without additional constraints.

---

## Proposed Tests

### 1. Deep MOND Rotation Curve Shape

**Prediction difference**:
- MOND: V continues to rise (or stay flat) as V ∝ M^0.25
- Synchronism: V flattens at lower value, then rises only with mass

**Test**: Extended HI rotation curves of isolated gas-rich dwarfs

### 2. Baryonic Tully-Fisher Relation Scatter

**If Synchronism is correct** with universal a₀ = 1.05×10⁻¹⁰:
- BTFR should have small scatter
- But normalization shifted from MOND expectation

**If MOND is correct** with empirical a₀ = 1.2×10⁻¹⁰:
- BTFR fits well (by construction)
- Synchronism would show systematic residuals

### 3. Mass Discrepancy-Acceleration Relation (MDAR)

**MDAR**: Plot M_dyn/M_baryon vs acceleration

**Synchronism prediction**:
```
M_dyn/M_baryon = (1 + f_indiff) × G_eff/G
               = (1 + f_indiff) / C(a)
```

**MOND prediction**:
```
M_dyn/M_baryon = ν(a_N/a₀)
```

At low acceleration, these diverge:
- Synchronism: bounded enhancement × (1 + f_indiff)
- MOND: unbounded enhancement

---

## Falsification Criteria

### Synchronism Falsified If:

1. Deep MOND systems show enhancement > 3.17× WITHOUT additional mass
2. BTFR fits with Synchronism a₀ give systematically worse χ²
3. UFDs show dynamics consistent with MOND ν-function, not C(a)

### MOND Falsified If:

1. Deep MOND enhancement shows bounded behavior
2. Correlation between inferred "dark matter" and system age/environment
3. Lensing shows less mass than dynamics predict (opposite of M_dyn > M_lens)

---

## Conclusions

### Key Findings

1. **12% a₀ difference is secondary**
   - The interpolating functions differ more importantly
   - Deep MOND limit is qualitatively different

2. **Bounded vs Unbounded Enhancement**
   - Synchronism: max G_eff/G = 3.17
   - MOND: enhancement can be arbitrarily large
   - This is the cleaner test

3. **Indifferent Mass Complicates Testing**
   - Synchronism invokes indifferent mass to make up the difference
   - Need independent constraints on indifferent mass fraction

4. **Current Data Likely Insufficient**
   - Both frameworks can fit typical galaxies
   - Deep MOND systems are rare and uncertain
   - Statistical analysis of large samples needed

### Next Steps

1. Implement full MCMC fitter for both models
2. Apply to SPARC-like sample
3. Focus on deep MOND systems specifically
4. Seek independent constraints on indifferent mass

---

## Files Created

- `simulations/session201_a0_precision.py` - a₀ sensitivity analysis
- `simulations/session201_rotation_curve_fitter.py` - Rotation curve models

---

*Session #201: The Synchronism vs MOND distinction is not just about a₀, but about the bounded vs unbounded nature of the enhancement factor in the deep MOND limit.*
