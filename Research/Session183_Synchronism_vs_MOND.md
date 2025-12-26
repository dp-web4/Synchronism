# Session #183: Synchronism vs MOND Discrimination

**Date**: December 26, 2025
**Machine**: CBP
**Status**: ✓ COMPLETE - Identified discriminating tests

---

## Executive Summary

Both Synchronism and MOND explain TDG mass discrepancies. This session identifies the **fundamental difference** between them and proposes **discriminating observational tests**.

**Key Finding**: The Chae et al. 2020 "EFE detection" may actually support Synchronism, since external field strength correlates with environment density.

---

## Fundamental Difference

| Theory | Mass Discrepancy Depends On |
|--------|----------------------------|
| **MOND** | Acceleration (a < a₀) |
| **Synchronism** | Density (ρ < ρ_t) |

This leads to qualitatively different predictions:

### MOND Predictions
- Same discrepancy at all densities for fixed acceleration
- External Field Effect (EFE): External gravity suppresses internal dynamics
- Contours in (ρ, a) plane are **horizontal**

### Synchronism Predictions
- Same discrepancy at all accelerations for fixed density
- No EFE: Only local density matters
- Contours in (ρ, a) plane are **vertical**

---

## Discriminating Tests

### Test A: Same Acceleration, Different Density

Compare TDGs vs isolated dwarfs at same V_rot:

| Object | V_rot | ρ/ρ_t | MOND M_dyn/M_bary | Sync M_dyn/M_bary |
|--------|-------|-------|-------------------|-------------------|
| TDG in tidal stream | 30 km/s | 0.1 | ~2.0 | **2.23** |
| Isolated dwarf | 30 km/s | 0.5 | ~2.0 | **1.71** |
| Cluster dwarf | 30 km/s | 2.0 | ~2.0 | **1.37** |

**Prediction**:
- MOND: All have same M_dyn/M_bary (acceleration-only)
- Synchronism: TDGs > isolated > cluster (density-dependent)

**Observable**: Match TDGs to isolated dwarfs by V_rot, compare M_dyn/M_bary

---

### Test B: BTFR Scatter Correlations

The Baryonic Tully-Fisher Relation scatter should correlate with:

| Theory | BTFR Scatter Correlation with Environment |
|--------|------------------------------------------|
| MOND | None (scatter is measurement error) |
| Synchronism | **Negative** (low ρ → high V residual) |

**Observable**:
1. Measure BTFR for SPARC galaxies
2. Estimate environment density from galaxy catalogs
3. Correlate BTFR residuals with environment

Simulation shows: r = -0.84 for Synchronism, r ≈ 0 for MOND

---

### Test C: External Field Effect

The EFE is a **MOND-specific prediction** with no analogue in Synchronism.

| Theory | Dwarf in Cluster Outskirts |
|--------|---------------------------|
| MOND | Reduced discrepancy (EFE suppression) |
| Synchronism | Increased discrepancy (low local ρ) |

**Critical Issue**: External field and local density are **correlated** in the cosmic web!

The Chae et al. 2020 "EFE detection" (4σ) may actually be measuring:
- MOND: g_ext effect on internal dynamics
- Synchronism: ρ_local effect on G_eff

**Needed Test**: Find systems where g_ext and ρ_local are **uncorrelated**
- Possible candidates: Galaxies at same distance from cluster but different local overdensities
- Tidal debris at different stages of expansion

---

### Test D: Density-Acceleration Plane

Map M_dyn/M_bary across the (log ρ, log a) plane:

- MOND: Horizontal contours
- Synchronism: Vertical contours

**Observable**: Compile sample of dwarfs with measured (ρ, a, M_dyn/M_bary)

---

## The Chae et al. 2020 Reinterpretation

### Original Claim
- EFE detected at 4σ in SPARC galaxies
- Rotation curves decline in outer parts when external field is strong
- Interpreted as MOND's EFE (violation of Strong Equivalence Principle)

### Synchronism Interpretation
The detected signal could be:

1. **External field → high density environment**
   - Galaxies in clusters/groups have both high g_ext AND high ρ
   - High ρ → G_eff ≈ G → Newtonian rotation curves
   - Appears as "EFE" but is actually density effect

2. **Observable correlation**
   - g_ext ∝ ρ_environment (from gravitational attraction)
   - Effect attributed to EFE may be density-dependent G_eff

### Discriminating Prediction

If Chae et al. detection is:
- **MOND EFE**: Effect should correlate with g_ext at fixed ρ
- **Synchronism ρ-effect**: Effect should correlate with ρ at fixed g_ext

**Proposed Test**: Re-analyze SPARC data with independent density estimates
1. Use galaxy catalog to estimate local ρ for each SPARC galaxy
2. Estimate g_ext independently
3. Partial correlation analysis: Does rotation curve decline correlate with ρ or g_ext?

---

## Existing Evidence Assessment

| Observation | MOND Interpretation | Synchronism Interpretation |
|-------------|--------------------|-----------------------------|
| TDG mass excess | Low a → MOND boost | Low ρ → G_eff enhancement |
| Chae 2020 EFE | g_ext suppresses boost | ρ_high → G_eff ≈ G |
| Fornax dwarfs disturbed | EFE reduces protection | Low ρ in tidal environment |
| BTFR tightness | Fundamental law | Average ρ determines average G_eff |

**Current Status**: Both theories can explain existing observations
**Needed**: Tests that disentangle ρ from a

---

## Proposed Observational Program

### Priority 1: TDG vs Isolated Dwarf Comparison
- Match TDGs to isolated dwarfs by V_rot
- Compare M_dyn/M_bary
- If TDGs systematically higher → Synchronism

### Priority 2: BTFR Environment Correlation
- Compute BTFR residuals for SPARC
- Correlate with environment density
- If significant correlation → Synchronism

### Priority 3: SPARC Reanalysis
- Independent ρ estimates for all galaxies
- Partial correlation of rotation curve shape with (ρ, g_ext)
- Determine which drives the "EFE" signal

### Priority 4: Cluster Edge Dwarfs
- Find dwarfs at same cluster-centric distance but different local ρ
- Compare velocity dispersions
- MOND: Same (same g_ext)
- Synchronism: Different (different ρ)

---

## Mathematical Summary

### MOND
```
a_eff = a_N / μ(a_N/a₀)

μ(x) = x / √(1 + x²)  [simple interpolation]

M_dyn/M_bary = 1/μ(a/a₀)
```

External Field Effect:
- If g_ext > a₀: Internal dynamics become Newtonian
- If g_int < g_ext < a₀: Intermediate suppression

### Synchronism
```
C(ρ) = Ω_m + (1 - Ω_m) × (ρ/ρ_t)^(1/φ) / [1 + (ρ/ρ_t)^(1/φ)]

G_eff = G / C(ρ)

M_dyn/M_bary = G_eff / G = 1/C(ρ)
```

No External Field Effect - only local density matters.

---

## Conclusions

1. **Synchronism and MOND make different predictions** despite both explaining TDGs

2. **The key difference is ρ-dependence vs a-dependence**

3. **Chae et al. 2020 EFE detection is ambiguous** - could be Synchronism density effect

4. **Discriminating tests are feasible** with existing data (SPARC + environment catalogs)

5. **TDG vs isolated dwarf comparison** is the cleanest test

---

## Files Created

- `simulations/session183_synchronism_vs_mond.py`
- `simulations/session183_sync_vs_mond.png`
- `Research/Session183_Synchronism_vs_MOND.md`

---

## References

1. [Chae et al. 2020](https://iopscience.iop.org/article/10.3847/1538-4357/abbb96) - "Testing the Strong Equivalence Principle: Detection of the External Field Effect in Rotationally Supported Galaxies"

2. [Triton Station Blog](https://tritonstation.com/2020/12/18/statistical-detection-of-the-external-field-effect-from-large-scale-structure/) - McGaugh on EFE detection

3. [MOND Scholarpedia](http://www.scholarpedia.org/article/The_MOND_paradigm_of_modified_dynamics) - EFE formalism

---

*Session #183: Identified discriminating tests between Synchronism and MOND*
