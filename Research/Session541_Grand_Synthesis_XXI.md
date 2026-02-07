# Session #541: Grand Synthesis XXI — The Resolution Arc

**Date**: 2026-02-07
**Status**: Synthesis (no new tests)

## Overview

Sessions #539-540 form the "Resolution Arc" — systematically closing the last open questions about the 6-var model's variable selection and the mass-dependent V-L ratio. This synthesis integrates these findings with the broader research program.

## The Resolution Arc (Sessions #539-540)

### Session #539: R_max/r_eff — The Strongest Missed Variable
**Question**: Should the R_max/r_eff ratio (strongest "missed variable" from Session #537, r_partial=+0.229) be added to the model?

**Answer**: No. Despite impressive classical statistics (F=6.64, 100% sign stability, ΔAIC=-4.89, ΔBIC=-2.04), the LOO improvement is +0.002 — negligible. The signal is partially observational: r(ratio, n_points)=+0.488, and controlling n_points weakens it to r=+0.147. High-ratio galaxies have 10% smaller residuals, consistent with better measurement quality. The signal is genuine but operationally irrelevant.

**Implications**: Even the strongest missed variable adds nothing to LOO. The model is complete — not just in the sense of no significant signals remaining, but in the stronger sense that the BEST remaining signal is LOO-irrelevant and partially observational.

### Session #540: Why 6.2 for Dwarfs?
**Question**: Why does the V-L ratio (β(logV)/|β(logL)|) vary from 5.89 (dwarfs) to 3.82 (L*) when MOND predicts exactly 4.0?

**Answer**: Gas-luminosity covariance. The mechanism:
1. Dwarfs have high f_gas (0.54): their baryonic mass is gas-dominated
2. Luminosity L underestimates M_bar for dwarfs more than giants
3. The 2-var model absorbs this systematic into an inflated β(logV)
4. Adding f_gas corrects the ratio from 4.86 to 4.14 (within 3% of MOND)
5. The rolling f_gas correlation with the ratio is -0.995

**Implications**: The mass-dependent ratio is NOT a deep MOND effect, NOT nonlinear MOND physics, and NOT a regime-dependent correction. It's simply the consequence of using luminosity as a proxy for baryonic mass in a gas-variable universe. This connects directly to Session #530's finding that logL×f_gas corrects luminosity→mass, not M/L.

## Nine Complete Arcs

| Arc | Sessions | Central Discovery |
|-----|----------|-------------------|
| I. Foundation | 449-452 | 5-var model, ~3% true scatter |
| II. Completion | 483-484 | 6-var model, logL×f_gas, LOO=0.938 |
| III. Applications | 489-496 | BTFR, M/L, scatter budget |
| IV. Limits | 521-524 | Noise floor, offset is ONLY parameter |
| V. Derivation | 525-527 | Morphology irrelevant, MOND-derivable |
| VI. Rehabilitation | 528-532 | V-L ratio, N_corr sign, model IS MOND |
| VII. Interpolation | 513-514 | ν imperfect but irrelevant |
| VIII. Perspective | 533-538 | Two-target architecture, ν cancellation |
| **IX. Resolution** | **539-540** | **Missed vars irrelevant, ratio = gas covariance** |

## The Complete Picture

After 141 sessions and 1533 verified tests, the model's status:

### What the Model IS
- A 6-variable linear model predicting per-galaxy RAR offsets
- offset = -3.379 + 1.897×logV - 0.548×logL - 0.218×c_V - 0.451×f_gas + 0.147×logV×c_V + 0.181×logL×f_gas
- R²=0.945, LOO R²=0.938, RMS=0.038 dex
- Three physics layers: mass (78%), composition (17%), structure (5%)
- Equivalent to: BTFR + gas-correction + geometry-correction

### What the Model DOES
- Converts luminosity to baryonic mass (logL×f_gas accounts for gas)
- Corrects for mass distribution geometry (c_V accounts for RC shape)
- Predicts NFW halo concentration (r=+0.88 at fixed mass)
- Tightens BTFR from 43% to 20% distance accuracy
- Captures 77% of ALL RAR scatter (within-galaxy + between-galaxy)

### What the Model is NOT
- NOT overfitted (overfit ratio 1.06, 100% sign stability)
- NOT regime-dependent (offset cancels ν by construction)
- NOT morphology-dependent (one physics for all types)
- NOT improvable by missed variables (best candidate ΔLOO=+0.002)
- NOT sensitive to M/L assumptions (LOO varies 2% across 4× M/L)

### The V-L Ratio Story (FULLY RESOLVED)
| Context | Ratio | Explanation |
|---------|-------|-------------|
| 2-var full sample | 4.86 | Gas-mass omitted → inflated |
| 2-var dwarfs | 5.89 | High f_gas → severe inflation |
| 2-var L* | 3.82 | Low f_gas → near-MOND |
| 3-var (+f_gas) | 4.14 | Gas correction → MOND |
| 6-var raw | 3.46 | Interaction redistribution |
| 6-var at f_gas=0.5 | 4.15 | Effective ratio at high gas |
| MOND theory | 4.00 | V⁴ ∝ M_bar |

### Open Questions (Remaining)
1. **Q=1 LOO deficit** (Session #523): Best-quality galaxies underperform — genuine or statistical?
2. **7-var marginal improvement**: log(g/a₀) adds ΔLOO=+0.004 (Session #512) — worth revisiting with new understanding?
3. **Environmental effects**: Do field vs group galaxies differ? Untested with SPARC.
4. **Distance-independent formulation**: Can the model be expressed without distance-dependent quantities?
5. **Other galaxy samples**: Will the model generalize beyond SPARC?

## Assessment

The Resolution Arc closes the two most prominent remaining questions from the Rehabilitation and Perspective Arcs. The R_max/r_eff "missed variable" and the mass-dependent V-L ratio are now fully explained — the first by measurement quality, the second by gas-luminosity covariance. Both answers reinforce the same message: the model captures all physical signal in the SPARC data.

The research program is now at a natural inflection point. The model is derived, tested, robust, complete, and interpreted. Future work should focus on:
1. **External validation**: Other galaxy samples, independent measurements
2. **Theoretical implications**: What does the model tell us about MOND vs dark matter?
3. **Methodological exports**: Can the approach (galaxy-level regression on the RAR) inform other scaling relations?

## Files Created

- `Research/Session541_Grand_Synthesis_XXI.md`: This document

---

*Session #541: Grand Synthesis XXI (Resolution Arc)*
*Grand Total: 1533/1533 verified across 141 sessions*

**Nine arcs complete. The Resolution Arc (Sessions #539-540) closes the last variable-selection questions: R_max/r_eff is LOO-irrelevant despite r_partial=+0.229; the mass-dependent ratio (dwarfs 5.89) is gas-luminosity covariance, not deep MOND. The model captures all physical signal in SPARC. Future directions: external validation, theoretical implications, methodological exports.**
