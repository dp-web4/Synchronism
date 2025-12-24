# Session #177 Brief: Void Galaxy BTFR Test

**Prepared**: 2025-12-24
**Priority**: HIGH - Unanimous recommendation from three-model peer review
**Type**: Falsification Test

---

## Objective

Test whether galaxies in different environments (void vs cluster) follow the same or different Baryonic Tully-Fisher Relation (BTFR).

**Synchronism prediction**: Void galaxies should show systematically higher V/V_bar (or lie on a different BTFR intercept) than cluster galaxies of the same baryonic mass.

**MOND prediction**: All galaxies lie on the same universal BTFR regardless of environment (a₀ is universal).

**ΛCDM prediction**: Scatter increases with environment, but no clean systematic offset based purely on local density.

---

## Why This Test Matters

From the three-model peer review (2025-12-24):

> **Nova (GPT-5.2)**: "This becomes hard to 'explain away' with baryonic modeling and is closer to your own 'directly measures coherence' claim."

> **Grok-4**: "A clean 2× factor would exceed [ΛCDM predictions]—your prediction is stronger, but only if voids are truly underdense."

> **Gemini-3-pro**: "If galaxies in deep voids and galaxies in rich clusters fall on the exact same BTFR, Synchronism is falsified. MOND would survive."

This is the **most discriminating test** identified by all three external reviewers.

---

## Data Sources

### Primary: SPARC Database
- 175 galaxies with high-quality rotation curves
- V_flat measurements available
- Baryonic mass (gas + stellar) computed
- Location: `/Synchronism/data/sparc/` or download from http://astroweb.cwru.edu/SPARC/

### Environment Classification Options

1. **2MRS Density Field** (preferred)
   - 2MASS Redshift Survey reconstruction
   - Provides δ = (ρ - ρ_mean)/ρ_mean at galaxy positions
   - Need to cross-match SPARC galaxies by coordinates

2. **SDSS DR7 Environment Catalog**
   - Tempel et al. void catalog
   - Group/cluster catalogs available
   - Coverage limited to SDSS footprint

3. **Cosmicflows-4 Density Reconstruction**
   - Already used in Sessions #169-175
   - Smoothed density field available
   - May have selection effects we identified

### Cross-matching
- SPARC provides RA/Dec for all galaxies
- Match to environment catalogs within ~1 arcmin tolerance
- Record environment metric (δ or void/field/cluster classification)

---

## Methodology

### Step 1: Data Preparation
```python
# Load SPARC data
# Key columns: Name, Distance, Vflat, Mstar, Mgas, Mbary
# Compute V_bar = sqrt(G * M_bary / R) at appropriate radius
```

### Step 2: Environment Cross-match
```python
# For each SPARC galaxy:
#   - Get RA, Dec, Distance
#   - Query environment catalog
#   - Assign δ (density contrast) or classification
#
# Define bins:
#   - Void: δ < -0.5 (or explicit void catalog membership)
#   - Field: -0.5 < δ < 1.0
#   - Cluster: δ > 5.0 (or group/cluster membership)
```

### Step 3: BTFR Analysis
```python
# Standard BTFR: log(V_flat) = a * log(M_bary) + b
#
# Fit separately for:
#   - All galaxies (baseline)
#   - Void galaxies only
#   - Field galaxies only
#   - Cluster galaxies only
#
# Compare:
#   - Intercept (b) differences
#   - Slope (a) consistency
#   - Scatter in each bin
```

### Step 4: Statistical Tests
```python
# 1. F-test: Does environment-split model fit better than single BTFR?
# 2. Residual analysis: Plot BTFR residuals vs environment
# 3. Bootstrap: Uncertainty on intercept differences
# 4. KS test: Are void and cluster residual distributions different?
```

---

## Specific Predictions

From Session #176 and the coherence function:

| Environment | ρ/ρ_cosmic | C(ρ) | M_dyn/M_bary | V/V_bar |
|-------------|------------|------|--------------|---------|
| Cluster core | 10,000 | ~1.0 | 1.00 | 1.00 |
| Field | 1.0 | 0.65 | 1.54 | 1.24 |
| Void interior | 0.2 | 0.49 | 2.04 | 1.43 |

**Quantitative prediction**:
- Void galaxies should show V_flat ~15-20% higher than cluster galaxies at fixed M_bary
- This corresponds to BTFR intercept shift of ~0.06-0.08 dex

---

## Success Criteria

### Synchronism SUPPORTED if:
- Void galaxies systematically above cluster galaxies on BTFR (>2σ)
- Intercept difference Δb > 0.05 dex in direction predicted
- Residual correlation with environment (r > 0.3)

### Synchronism FALSIFIED if:
- All galaxies lie on same BTFR regardless of environment (<0.02 dex difference)
- No significant residual-environment correlation
- This would support universal MOND over density-dependent Synchronism

### Inconclusive if:
- Sample size too small after environment cuts
- Environment classification unreliable
- Systematics dominate (morphology correlates with environment)

---

## Controls Required

The reviewers emphasized controlling for confounders:

1. **Morphology**: Void galaxies tend to be later-type. Match by Hubble type or restrict to late-types only.

2. **Surface brightness**: LSB galaxies preferentially in voids. Match by central surface brightness or include as covariate.

3. **Stellar mass estimation**: Different SFH in voids could bias M_star. Use gas-dominated galaxies where M_gas >> M_star.

4. **Distance**: More distant galaxies have larger V_flat errors. Weight by distance uncertainty or restrict to nearby sample.

---

## Files to Create

1. `simulations/session177_void_btfr_test.py` - Main analysis script
2. `simulations/session177_environment_crossmatch.py` - Data preparation
3. `data/sparc_environment_catalog.csv` - Cross-matched dataset
4. `simulations/results/session177_btfr_analysis.json` - Results
5. `Research/Session177_Void_Galaxy_BTFR_Test.md` - Full session notes

---

## Expected Duration

- Environment cross-matching: 1-2 hours
- BTFR analysis: 1-2 hours
- Statistical tests and controls: 1-2 hours
- Documentation: 1 hour

**Total**: ~4-6 hours for thorough analysis

---

## References

- McGaugh et al. (2016) - SPARC database paper
- Lelli et al. (2017) - BTFR with SPARC
- Tempel et al. (2014) - SDSS void catalog
- Lavaux & Hudson (2011) - 2MRS density field

---

## Note on Outcome

This test has real falsification potential. If the result is negative (no environment dependence), we must honestly report that Synchronism's density-dependent prediction is not supported, and MOND's universal a₀ is preferred.

The goal is truth, not confirmation.

---

*Brief prepared by Claude Opus 4.5 following three-model peer review (Nova GPT-5.2, Grok-4, Gemini-3-pro)*
