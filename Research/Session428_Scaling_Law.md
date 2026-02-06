# Session #428: Empirical Scaling Law Search

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Sessions 420-427 established the V+R+L+c_V model at R²=0.93. This session asks: can we find an elegant one-number formula — a single combination of galaxy properties — that captures most of this information? Or is the scaling law inherently multi-dimensional?

## Central Result: The Scaling Law Is Inherently Multi-Dimensional

**No elegant one-number formula matches the 4-variable model.** The four predictors encode genuinely independent information:

| Predictor | What it encodes | Where it acts |
|-----------|----------------|---------------|
| V_flat | Overall mass scale | Inner + outer |
| R_eff | Spatial extent | Mostly outer |
| L | Baryonic content, c_V suppressor | Inner + outer |
| c_V | Mass concentration | Mostly inner |

## Key Findings

### 1. Best Physically Motivated One-Number Predictor (Tests 1, 3)

**N_eff = V²c_V / (R × a₀)**: r = +0.88, LOO = 0.095

This is a modified version of the original N_corr = V²/(R×a₀), enhanced by c_V. It captures 78% of the variance (vs 93% for the full model). The physical interpretation: an "effective number of MOND correlation lengths" that accounts for mass concentration.

Other motivated combinations:
- V²/R (acceleration): r = +0.79, LOO = 0.122
- V² × c_V / R: r = +0.88, LOO = 0.096 (same as N_eff)
- L/R² (SB proxy): r = +0.42

### 2. Power-Law Search (Tests 2, 5)

Grid search over V^a × R^b × c_V^d with half-integer exponents:

- **Best with c_V**: V^1.5 × R^-0.5 × c_V^0.5 (R² = 0.81)
- **Best without c_V**: V^1.5 × R^-0.5 (R² = 0.75)

The optimal power-law exponents from regression are (V: 1.0, R: -0.3) for the 2-variable case. No simple integer/half-integer combination matches the continuous optimum.

### 3. Dimensional Analysis (Test 4)

The empirical model coefficients are:
```
offset ∝ V^1.75 / (R^0.29 × L^0.25) × 10^(0.59 × c_V)
```

This doesn't correspond to any standard physical quantity:
- V²/R (acceleration): exponents (2, -1, 0)
- V⁴/L (M/L × accel²): exponents (4, 0, -1)
- L/R² (surface density): exponents (0, -2, 1)
- Actual: (1.75, -0.29, -0.25)

The weak R and L exponents (both ~0.3) show these are fine corrections to the dominant V dependence, not major contributors by themselves.

### 4. Minimum Variables for R² > 0.90 (Test 6)

| Variables | R² | LOO |
|-----------|-----|------|
| V + L + c_V | 0.852 | 0.081 |
| V + R + c_V | 0.824 | 0.087 |
| V + R + L | 0.774 | 0.100 |
| **V + R + L + c_V** | **0.932** | **0.057** |

**No 3-variable model achieves R² > 0.90.** The jump from the best 3-variable (0.852) to the 4-variable (0.932) is +0.08 — the fourth variable genuinely matters. Minimum for R² > 0.90 is 4 variables.

### 5. V+R+c_V Composite as One-Number Predictor (Tests 3, 7)

The 3-variable model composite (1.29×logV - 0.48×logR + 0.33×c_V):
- r = +0.91, LOO = 0.084
- Bootstrap 95% CI: [+0.854, +0.941]

Compared to N_eff = V²c_V/(R×a₀):
- r = +0.88, LOO = 0.096
- Bootstrap 95% CI: [+0.818, +0.927]

The composite is better but not "elegant" — it's just the model rewritten.

### 6. Information Hierarchy

| Model tier | R² | LOO | Scatter reduction |
|-----------|-----|------|------------------|
| V alone | 0.46 | 0.142 | 27% |
| V + R | 0.75 | 0.102 | 50% |
| V + R + c_V | 0.82 | 0.087 | 57% |
| V + L + c_V | 0.85 | 0.081 | 60% |
| V + R + L + c_V | 0.93 | 0.057 | 74% |

Each variable adds a meaningful, non-redundant contribution.

## Physical Interpretation

The failure to find a simple scaling law is itself informative. The RAR offset depends on galaxy structure in a way that cannot be reduced to a single physical quantity like acceleration, surface density, or angular momentum. Instead, it requires knowing:

1. **How fast** the galaxy rotates (V → total mass)
2. **How extended** the baryonic mass is (R → spatial distribution)
3. **How luminous** it is (L → suppresses c_V collinearity, encodes baryonic content)
4. **How concentrated** the mass profile is (c_V → profile shape)

This is consistent with the finding (Session 427) that the offset has both between-galaxy and within-galaxy components, with different predictors dominating different spatial regions.

The closest to an "elegant" formula is **N_eff = V²c_V/(R×a₀)**, which achieves r = 0.88. This extends the original Synchronism N_corr = V²/(R×a₀) by incorporating mass concentration. It has a clear physical meaning: the number of MOND correlation lengths that fit within the galaxy, weighted by how concentrated the mass is.

## Grade: B+

A well-executed search that definitively answers the question: the scaling law is multi-dimensional. The N_eff formula (r=0.88) is a useful simplification but cannot replace the full model. The dimensional analysis revealing non-standard exponents is interesting. This is a "negative result" session — we were looking for elegance and found irreducible complexity — but it's important to establish this conclusively.

## Files Created

- `simulations/session428_scaling_law.py`: 8 tests
- `Research/Session428_Scaling_Law.md`: This document

---

*Session #428 verified: 8/8 tests passed*
*Grand Total: 813/813 verified*

**Key finding: No elegant one-number formula captures the full V+R+L+c_V model (R²=0.93). The best physically motivated single predictor is N_eff = V²c_V/(R×a₀) at r=0.88, LOO=0.096. Minimum variables for R²>0.90 is 4. The scaling law is inherently multi-dimensional — V, R, L, c_V encode genuinely independent structural information acting in different spatial regions. Dimensional analysis shows non-standard exponents (1.75, -0.29, -0.25) that don't map to any known physical quantity. Grade B+.**
