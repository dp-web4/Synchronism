# Session #176: Cluster Velocity Dispersion Profiles

**Date**: December 24, 2025
**Focus**: Theoretical Framework for Cluster Dynamics
**Follow-up to**: Session #175 (Real Data Application Arc Conclusion)

---

## Executive Summary

Session #176 developed the theoretical framework for testing Synchronism via cluster velocity dispersion profiles. Three analyses were performed:

1. **176a**: Initial calculation (found bug in cosmological parameters)
2. **176b**: Corrected Jeans equation analysis
3. **176c**: Dynamical vs lensing mass discrepancy prediction

**Key Findings**:
- Cluster interior enhancement is SMALL (~3% at R_200, ~26% at 5 R_200)
- Environment-dependent M_dyn/M_lens ratio provides discriminating test
- New insight: Transition density ρ_t may be SCALE-DEPENDENT

---

## Background: Why Cluster Dynamics?

Session #175 concluded the Real Data Application arc with the finding that peculiar velocities are the wrong observable for testing the coherence function (MRH mismatch: kpc vs Mpc).

Session #173 showed a promising result: cluster velocity dispersion ratio (outer/inner) = 3.28, in the Synchronism direction.

This session develops the theory to interpret that observation.

---

## Theoretical Framework

### Coherence Function

The Synchronism coherence function:

```
C(ρ) = Ω_m + (1 - Ω_m) × (ρ/ρ_t)^(1/φ) / [1 + (ρ/ρ_t)^(1/φ)]
```

where:
- φ = golden ratio (1.618...)
- Ω_m = 0.3
- ρ_t = transition density

Effective gravity: G_eff = G/C(ρ)

### Key Behavior

| ρ/ρ_cosmic | C(ρ) | G_eff/G |
|------------|------|---------|
| 0.01 | 0.34 | 2.95 |
| 0.1 | 0.44 | 2.29 |
| 1.0 | 0.65 | 1.54 |
| 10 | 0.86 | 1.16 |
| 100 | 0.96 | 1.04 |
| 1000 | 0.99 | 1.01 |

---

## Session #176b: Corrected Cluster Profile

### NFW Cluster Parameters

For a Coma-like cluster (M_200 = 10^15 M_sun, c = 4):
- R_200 = 2.06 Mpc = 2063 kpc
- r_s = 0.52 Mpc = 516 kpc

### Velocity Dispersion from Jeans Equation

Modified Jeans equation with G_eff:

```
σ²(r) = (1/ρ(r)) ∫_r^∞ ρ(s) × G_eff(s) × M(<s)/s² ds
```

### Results

| r/R_200 | ρ/ρ_cosmic | σ_ΛCDM (km/s) | σ_Sync (km/s) | Ratio |
|---------|------------|---------------|---------------|-------|
| 0.01 | 406,123 | 619 | 619 | 1.000 |
| 0.1 | 22,411 | 961 | 963 | 1.002 |
| 0.5 | 976 | 931 | 941 | 1.010 |
| 1.0 | 176 | 836 | 860 | 1.028 |
| 2.0 | 27 | 717 | 773 | 1.078 |
| 5.0 | 2 | 553 | 698 | 1.262 |

### Key Prediction

Outer/Inner dispersion ratio:
- **ΛCDM**: σ(2 R_200) / σ(0.1 R_200) = **0.746**
- **Synchronism**: σ(2 R_200) / σ(0.1 R_200) = **0.803**
- **Difference**: 7.6%

**Critical Insight**: Enhancement is only significant at r > R_200 where density drops toward cosmic mean.

---

## Session #176c: Mass Discrepancy Test

### Theoretical Prediction

The key discriminating test:

```
M_dyn/M_lens = G_eff/G = 1/C(ρ)
```

This directly measures the coherence function!

### Predictions by Environment

| Environment | ρ/ρ_cosmic | M_dyn/M_lens |
|-------------|------------|--------------|
| Cluster core | 10,000 | 1.002 |
| Cluster R_200 | 200 | 1.026 |
| Cluster outskirts | 10 | 1.157 |
| Filament | 3 | 1.308 |
| Wall | 1.5 | 1.442 |
| Field | 1.0 | 1.538 |
| Void edge | 0.5 | 1.736 |
| Void interior | 0.2 | 2.045 |
| Deep void | 0.05 | 2.532 |

### Key Testable Prediction

**Void galaxies should have 2.0× higher M_dyn/M_stellar compared to cluster galaxies of same stellar mass.**

In ΛCDM: This would require dark matter halos to be denser in voids (opposite to simulations).

In Synchronism: Natural consequence of coherence function.

---

## New Insight: Scale-Dependent Transition Density

### The Problem

Galaxy rotation curves show strong "dark matter" effect at ρ >> ρ_cosmic.
Cluster dynamics shows weak effect even at ρ ~ ρ_cosmic.

How can the same coherence function explain both?

### The Solution

The transition density ρ_t may be **scale-dependent**:

```
C(ρ, MRH) = Ω_m + (1 - Ω_m) × (ρ/ρ_t(MRH))^(1/φ) / [1 + (ρ/ρ_t(MRH))^(1/φ)]
```

where ρ_t scales with the characteristic density at each MRH.

### Calibration

| Scale | MRH | ρ_t |
|-------|-----|-----|
| Galaxy | ~10-50 kpc | ~10^7-10^8 M_sun/Mpc³ |
| Cluster | ~1-10 Mpc | ~10^10 M_sun/Mpc³ (cosmic mean) |

**Ratio**: ρ_t(cluster) / ρ_t(galaxy) ~ 100-1000

This explains:
1. Strong "dark matter" in galaxy rotation curves
2. Weaker effect in cluster dynamics
3. Same coherence function form at different scales

**Consistent with RESEARCH_PHILOSOPHY.md**: "MRH must match complexity."

---

## Interpretation of Session #173 Result

Session #173 found dispersion ratio = 3.28 using:
- Inner: < 5 Mpc
- Outer: 20-60 Mpc

This compares:
- Cluster-bound dynamics (inner)
- **Cosmic field** dynamics (outer)

At 20-60 Mpc, density is ~1× cosmic mean, so G_eff/G → 1.54.

The high ratio reflects the **field-to-cluster** transition, not intra-cluster variation.

For proper cluster test, need σ(r) profile **within** R_200.

---

## Recommended Tests

### 1. Void vs Cluster Galaxy Comparison

**Data**: SDSS void catalog + spectroscopic velocity dispersions

**Method**:
- Select galaxies with similar stellar mass
- Compare M_dyn/M_stellar between environments
- Void galaxies should show 2× higher ratio

### 2. Cluster σ(r) Profiles

**Data**: Planck PSZ2 clusters with spectroscopic members

**Method**:
- Measure σ(r) from 0.1 to 3 R_200
- Compare to Jeans equation predictions
- Look for 8% excess at 2 R_200

### 3. Dynamical vs Lensing Mass

**Data**: Weak lensing + dynamics for same clusters

**Method**:
- Measure M_dyn from velocity dispersion
- Measure M_lens from weak lensing
- Compare M_dyn/M_lens as function of radius
- Should increase outward in Synchronism

---

## Files Created

| File | Description |
|------|-------------|
| `session176_cluster_dispersion_theory.py` | Initial (buggy) calculation |
| `session176b_cluster_dispersion_corrected.py` | Corrected Jeans analysis |
| `session176c_mass_discrepancy.py` | M_dyn/M_lens predictions |
| `session176_theory.png` | Visualization (buggy) |
| `session176b_corrected.png` | Corrected profiles |
| `session176c_mass_discrepancy.png` | Environment predictions |

---

## Key Takeaways

1. **Cluster enhancement is small**: ~3-8% at accessible radii
2. **M_dyn/M_lens is the key observable**: Directly measures G_eff/G
3. **Scale-dependent ρ_t explains multi-scale behavior**
4. **Void galaxies are ideal test**: 2× expected signal
5. **Session #173 measured field-cluster transition**, not intra-cluster profile

---

## Next Steps

1. Search for void galaxy velocity dispersion data
2. Compare M_dyn/M_stellar by environment using SDSS
3. Develop scale-dependent ρ_t formalism
4. Test against dwarf spheroidal galaxies

---

*Session #176 completed: December 24, 2025*
