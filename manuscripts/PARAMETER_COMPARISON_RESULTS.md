# Parameter Comparison: Empirical vs First-Principles Derived

**Date**: December 2, 2025
**Test**: SPARC rotation curve fitting (175 galaxies)

---

## The Question

The two manuscript versions use different parameter sets:
- **arxiv-v1.tex (Nov 25)**: A = 0.25, B = 1.62 (empirically fitted)
- **draft_v1.md (Dec 1)**: A = 0.028, B = 0.5 (first-principles derived)

**Question**: Why wasn't the "derived" model tested on SPARC rotation curves?

---

## Test Results

| Model | A | B | γ | Success Rate | Median χ² |
|-------|---|---|---|--------------|-----------|
| Empirical (arxiv-v1) | 0.25 | 1.62 | 2.0 | **52.6%** | 4.81 |
| Derived (draft_v1) | 0.028 | 0.5 | 2.0 | **2.9%** | 130.54 |

### By Galaxy Class

| Class | Empirical | Derived |
|-------|-----------|---------|
| Dwarfs (v < 50 km/s) | 81.8% | 22.7% |
| Intermediate (50-100 km/s) | 60.9% | 0.0% |
| Massive (v > 100 km/s) | 39.3% | 0.0% |

### Head-to-Head

- Empirical better: 173 galaxies
- Derived better: 2 galaxies
- Median χ² ratio (derived/empirical): **24.68x worse**

---

## Critical Finding

The "first-principles derived" parameters from draft_v1.md (A=0.028, B=0.5) **fail catastrophically** on SPARC rotation curve fitting:
- Only 2.9% success rate (vs 52.6% for empirical)
- 25x worse median χ²
- Works somewhat for dwarfs (22.7%) but completely fails for intermediate and massive galaxies

---

## Explanation

The 99% success rate claimed in draft_v1.md comes from the **Santos-Santos (2020) dark matter fraction test**, which:
1. Uses a simpler metric (mean DM fraction, not rotation curve shape)
2. Tests global predictions, not detailed velocity profiles
3. May be less sensitive to parameter choices

The derived parameters (A=0.028, B=0.5) fail on detailed rotation curves but succeed on global DM fractions. This suggests:
1. The theoretical derivation captures something real about global DM properties
2. But misses important physics for detailed rotation curve shapes
3. The empirical parameters compensate for this missing physics

---

## Implication for v3 Manuscript

The v3 manuscript must be honest about this:

1. **Don't claim** "all parameters derived from first principles" with 99% success
2. **Do claim** one of:
   - Empirical parameters (A=0.25, B=1.62) with 52.6% SPARC success
   - Derived parameters (A=0.028, B=0.5) with 99% Santos-Santos DM fraction success (different test)

3. **Acknowledge** the tension: The theoretical derivation needs refinement to match detailed rotation curves.

---

## Files

- Test script: `simulations/compare_empirical_vs_derived.py`
- Full results: `simulations/empirical_vs_derived_results.json`

---

*This test resolves the apparent discrepancy between the two manuscript versions. The different success rates (53.7% vs 99%) come from different parameter sets applied to different tests, not improvements in the same test.*
