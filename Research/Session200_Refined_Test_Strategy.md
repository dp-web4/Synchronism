# Session #200: Refined Test Strategy for Synchronism

**Date**: December 30, 2025
**Machine**: CBP
**Status**: STRATEGY DOCUMENT

---

## Context

Sessions #199-200 investigated cluster mass measurements as a test of Synchronism. Key findings:

1. **M_dyn/M_lens ≈ 1.1** observed, but Synchronism predicts G_eff/G ≈ 1.9
2. **Velocity anisotropy** (β ~ 0.3-0.4) explains the discrepancy
3. **Caustic mass method** is less anisotropy-dependent but still shows ~1.2-1.3
4. **Abell 520** is consistent with framework (complex geometry, not anomalous physics)

The challenge: Both ΛCDM and Synchronism can fit the data by invoking different "corrections."

---

## The Core Problem

| What we predict | What we observe | The gap |
|-----------------|-----------------|---------|
| M_dyn/M_lens = 1.9 | M_dyn/M_lens = 1.1 | Factor ~1.7 |

**ΛCDM says**: Ratio should be 1.0, deviations are "systematics"
**Synchronism says**: Ratio should be 1.9, reduced by anisotropy to 1.1

Both interpretations are self-consistent. We need a DIFFERENTIAL test.

---

## Refined Test Strategy

### Test 1: RADIAL TREND (Cleanest)

**The Prediction**:
- Inner cluster (high a): M_dyn/M_lens ~ 1.3-1.5
- Outer cluster (low a): M_dyn/M_lens ~ 2.0-2.5
- Trend should be monotonic increase with radius

**Why it works**:
- Independent of absolute calibrations
- If anisotropy β increases with radius, it PARTIALLY cancels G_eff increase
- But residual trend should remain
- ΛCDM predicts flat ratio = 1.0

**Data needed**:
- Same clusters with M_dyn AND M_lens at multiple radii
- Weak lensing profiles (available from HSC, DES, KiDS)
- Velocity dispersion profiles (rarer)
- Caustic mass profiles (available from CIRS, HeCS)

**Key papers to examine**:
- Umetsu et al. (2016) - Stacked lensing profiles
- Rines et al. (2016) - Caustic profiles
- CLASH program - Multi-radius data

### Test 2: ANISOTROPY-CORRECTED MASSES

**The Prediction**:
- After proper β(r) correction, M_dyn,corrected/M_lens = G_eff/G = 1/C(a)
- No residual "systematic" after correction

**Challenge**:
- Proper anisotropy requires orbit modeling
- Very few clusters have this done properly
- Model-dependent

**Data needed**:
- Schwarzschild or Made-to-Measure orbit models
- Phase-space reconstructions
- Individual cluster analyses, not stacks

**Key papers**:
- Mamon et al. - MAMPOSSt method
- Binney & Mamon (1982) - Jeans modeling fundamentals
- van der Marel et al. - Orbit superposition

### Test 3: M_caustic vs M_σ WITHIN CLUSTERS

**The Prediction**:
- M_caustic/M_σ = 1/(1-β)
- For β = 0.3-0.4: ratio = 1.4-1.7
- Observed: ~1.1-1.3

**Interpretation**:
- If β is lower than expected (~0.1-0.2), ratio would be ~1.1-1.2 ✓
- OR both methods have similar calibration biases

**Data needed**:
- Same clusters with both caustic AND dispersion masses
- Consistent apertures
- HeCS-SZ has this (Rines et al.)

### Test 4: DEEP MOND REGIME (Ultra-faint Dwarfs)

**The Prediction**:
- For a << a₀: G_eff/G → 1/Ω_m ≈ 3.17
- Ultra-faint dwarfs should show maximal enhancement
- BTFR slope → 0.25 in this limit

**Advantage**:
- Galaxies, not clusters (different systematics)
- Deep MOND regime (no ambiguity about a/a₀)
- Large sample from Gaia + spectroscopy

**Challenge**:
- These are dispersion-supported, not rotation-supported
- Need proper dynamical modeling
- Dark matter halos in ΛCDM are also maximally dominant here

**Data**:
- SPARC database (some UFDs)
- Gaia + spectroscopic surveys
- Simon & Geha (2007) and subsequent work

### Test 5: SATELLITE TIMING ARGUMENTS

**The Prediction**:
- Timing mass of Local Group uses G
- If G_eff applies, timing mass is overestimated
- Or equivalently, satellites move faster than expected

**Current status**:
- M31-MW timing gives M_total ~ 5×10^12 M_sun
- This is in rough agreement with other estimates
- No obvious anomaly

**Challenge**:
- Timing argument has large uncertainties
- Proper motion errors dominate
- Hard to distinguish G from G_eff

### Test 6: EXTERNAL FIELD EFFECT (EFE)

**MOND-specific Prediction**:
- External gravitational field affects internal dynamics
- Galaxies in clusters should differ from isolated galaxies
- Synchronism: Does C(a) depend on external a_ext?

**Current formulation**:
- Our C(a) uses local acceleration only
- Could be extended: C = C(a_local + a_external)?
- This would be a refinement

**Data**:
- Compare satellite galaxy kinematics (in cluster field)
- To field galaxy kinematics (isolated)
- Some MOND studies have examined this

---

## Priority Ranking

| Test | Feasibility | Discriminating Power | Priority |
|------|-------------|---------------------|----------|
| Radial trend | Medium | High | **#1** |
| β-corrected M | Low | High | #2 |
| M_caustic/M_σ | High | Medium | #3 |
| Ultra-faint dwarfs | Medium | Medium | #4 |
| Timing arguments | Low | Low | #5 |
| EFE | Low | High | Future |

---

## Immediate Next Steps

### Session #201-202: Precision a₀ Measurement

From Session #198 roadmap, this is parallel track. Use ultra-faint dwarfs to distinguish:
- Synchronism: a₀ = 1.05 × 10⁻¹⁰ m/s²
- MOND empirical: a₀ = 1.2 × 10⁻¹⁰ m/s²

12% difference is testable with precision data.

### Session #203-205: Radial Trend Compilation

Focus on collecting radial M_dyn/M_lens data:
1. Identify clusters with BOTH lensing AND dynamical radial profiles
2. Compile from literature
3. Test for radial trend

### Session #206+: Theoretical Refinements

If observational tests inconclusive:
1. Consider refinements to C(a) formula
2. Explore external field effects
3. Investigate indifferent pattern properties

---

## What Would Falsify Synchronism?

1. **Radial trend is FLAT**
   - M_dyn/M_lens ≈ 1.0 at ALL radii
   - After proper anisotropy correction

2. **a₀ clearly wrong**
   - Precision measurements give a₀ = 1.2 × 10⁻¹⁰
   - Not 1.05 × 10⁻¹⁰

3. **Indifferent mass properties wrong**
   - Self-interaction cross-section > 1 cm²/g
   - Would contradict Bullet Cluster

4. **Ultra-faint dwarfs show wrong G_eff**
   - Deep MOND should have G_eff/G → 3.17
   - If clearly different, formula is wrong

---

## What Would Validate Synchronism?

1. **Radial trend matches G_eff/G profile**
   - Increasing ratio with radius
   - Matches 1/C(a) prediction

2. **a₀ = 1.05 × 10⁻¹⁰ measured**
   - Precision BTFR or UFD analysis
   - Distinguishes from MOND empirical value

3. **Anisotropy-corrected M_dyn/M_lens = G_eff/G**
   - After proper β modeling
   - Removes the "systematic" explanation

4. **Consistent across scales**
   - Dwarfs, spirals, groups, clusters
   - All show same G_eff pattern

---

## Conclusion

The key challenge is that both ΛCDM and Synchronism can fit current data by adjusting "nuisance parameters" (anisotropy, systematics).

The path forward is:
1. **Differential tests** (radial trends, not absolute values)
2. **Precision measurements** (a₀ value)
3. **Independent methods** (caustic vs dispersion vs lensing)

Session #200 established the strategy. Sessions #201+ will execute.

---

*"The test of a first-rate intelligence is the ability to hold two opposed ideas in mind at the same time and still retain the ability to function." - F. Scott Fitzgerald*

*Applied to physics: Both ΛCDM and Synchronism currently function. The data will eventually distinguish them.*
