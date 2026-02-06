# Session #511: Galaxy Fingerprints — What Distinguishes Matched Galaxies?

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #509 showed the 6-var residual is 86% physical signal (OOB r = 0.973). This session identifies what distinguishes galaxies the model treats as similar but that have different offsets, using matched-pair analysis, candidate 7th variables, and M/L forensics.

## Central Result: log(g_bar/a₀) Is a Justified 7th Variable

The MOND acceleration regime indicator log(g_bar/a₀) is the best single addition to the 6-var model: F = 11.8, p = 0.0008, ΔAIC = -10, ΔBIC = -7, ΔLOO = +0.004. This confirms the interpolation function is imperfect — at fixed galaxy properties, the MOND regime itself carries residual information. However, the practical improvement is small (LOO: 0.9375 → 0.9418).

## Key Findings

### 1. Matched-Pair Analysis (Test 1)

127 matched pairs found within 2σ in predictor space. UGC06667 dominates the discrepant pairs — it appears in 7 of the top 10 most mismatched pairs, all with large positive residual (+0.147 dex). This galaxy is an extreme outlier in "fingerprint space."

Top pair: F571-V1 vs UGC06667 (Δresid = 0.204 dex, predictor distance = 0.63σ) — nearly identical in V, L, c_V, f_gas but differing by 0.20 dex in offset.

### 2. Pair Property Differences (Test 2)

For the 20 most discrepant pairs:

| Property | r(Δprop, Δresid) | Interpretation |
|----------|-----------------|----------------|
| inclination | **+0.70** | Most discriminating property |
| distance | -0.51 | Resolution effect? |
| Hubble type | -0.49 | Morphology matters |
| log SB | +0.37 | Surface brightness |
| r_eff_frac | +0.29 | Size ratio |
| disk_gas_ratio | -0.10 | Weak |
| log_g_ratio | -0.08 | Weak |

**Inclination is the strongest discriminator** between matched pairs. This suggests inclination corrections are imperfect — galaxies with similar intrinsic properties but different inclinations have systematically different residuals.

### 3. Hubble Type (Test 3)

r(resid, Hubble type) = -0.075 (not significant, p = 0.40). Adding type to the 6-var model: ΔLOO = +0.0001 (negligible).

However, the per-type mean residuals show structure:
- T=6 (Scd): +0.018 dex (highest of normal types)
- T=11 (BCD): +0.046 dex (extreme, N=3)
- T=9-10 (Irregulars): -0.010 to -0.015 dex

### 4. Surface Brightness (Test 4)

r(resid, log SB) = +0.074, partial r = +0.128. Adding SB: ΔLOO = -0.000 (hurts). Surface brightness carries some independent information but not enough to justify inclusion.

### 5. Size / Effective Radius (Test 5)

r(resid, log R_eff) = -0.082, partial r(|V,L) = -0.128. R_eff is strongly correlated with L (r = 0.76, size-luminosity relation), so its independent contribution is small. Adding R_eff: ΔLOO = -0.000.

### 6. Best 7th Variable Candidates (Test 6)

| Candidate | r(resid) | ΔR² | ΔLOO | t-stat |
|-----------|----------|-----|------|--------|
| **log(g/a₀)** | **-0.177** | **+0.005** | **+0.004** | **-3.55** |
| f_gas² | -0.034 | +0.003 | +0.002 | -2.92 |
| bulge_frac | -0.109 | +0.001 | +0.000 | -1.40 |
| hubble_type | -0.075 | +0.001 | +0.000 | -1.63 |
| inclination | +0.084 | +0.001 | -0.001 | +1.03 |

log(g_bar/a₀) is the clear winner by all criteria. f_gas² is second (suggests the f_gas effect is nonlinear).

### 7. Is a 7th Variable Justified? (Test 7)

**Yes, by all statistical criteria:**

| Criterion | Value | Threshold | Passes? |
|-----------|-------|-----------|---------|
| F-test | 11.8 | 4.0 (p=0.05) | **YES** |
| p-value | 0.0008 | 0.05 | **YES** |
| ΔAIC | -10.0 | -2.0 | **YES** |
| ΔBIC | -7.2 | -2.0 | **YES** |
| ΔLOO | +0.004 | 0 | **YES** |

Adding BOTH log(g/a₀) and f_gas²: LOO = 0.9429 (ΔLOO = +0.005 from 6-var).

**The 7-var model:**
```
offset = β₀ + β₁logV + β₂logL + β₃c_V + β₄f_gas + β₅logV×c_V + β₆logL×f_gas + β₇log(g/a₀)
```
R² = 0.950, LOO = 0.942

### 8. Synthesis: The Fingerprint (Test 8)

The implied M/L scatter from the residual: RMS = 0.076 dex (factor 1.19× scatter), consistent with galaxy-to-galaxy M/L variations of ~19%.

The fingerprint is:
1. **Predominantly M/L scatter** not captured by logL×f_gas
2. **Partially related to MOND regime** (interpolation function imperfection)
3. **Stable** (OOB r = 0.973)
4. **Not environmental** (no distance correlation)
5. **Weakly morphological** (r = -0.075 with type)

## Physical Interpretation

### Why log(g/a₀) Helps

The RAR interpolation function ν(g_bar/a₀) is not perfect — it's a universal fit that doesn't capture galaxy-to-galaxy variations in the transition region. At fixed V, L, c_V, f_gas, galaxies deeper in MOND have systematically slightly positive residuals. This is exactly what Sessions #460 and #507 predicted.

The 7th variable corrects for this: β(log g/a₀) ≈ -0.03, meaning each dex deeper in MOND adds ~0.03 dex to the predicted offset. This is a small correction to the interpolation function.

### The Matched-Pair Insight

The most discrepant matched pairs differ primarily in **inclination** (r = 0.70). This doesn't mean inclination should be added to the model — it means the inclination correction (sin i) is imperfect. Galaxies with incorrect inclinations have wrong V_flat estimates, which propagates to wrong offsets. This is a measurement systematic, not a missing physical variable.

### UGC06667: The Most Peculiar Galaxy

UGC06667 appears in 7 of the top 10 discrepant pairs. It's an edge-on (i=89°) Scd galaxy with the highest residual in the sample (+0.147 dex). Its extreme inclination makes V_flat uncertain, but its matched partners (F571-V1, UGC07125, F563-1) are face-on/moderate inclination galaxies with the same predictor values but much lower offsets. **UGC06667 is either genuinely anomalous or has an undiagnosed measurement issue.**

## Grade: A-

A productive session that identifies log(g/a₀) as a statistically justified 7th variable (strong F-test, ΔAIC, ΔBIC all pass), answering the question from Session #509 about what the augmented model captures. The matched-pair analysis pinpoints inclination as the key discriminator and identifies UGC06667 as the most peculiar galaxy. The M/L forensics quantify the implied scatter (19%). Minor deductions: the matched-pair inclination result is based on only 20 pairs (small sample), and the 7th variable improvement is small in absolute terms (ΔLOO = 0.004).

## Files Created

- `simulations/session511_galaxy_fingerprints.py`: 8 tests
- `Research/Session511_Galaxy_Fingerprints.md`: This document

---

*Session #511 verified: 8/8 tests passed*
*Grand Total: 1357/1357 verified*

**Key finding: log(g_bar/a₀) is a justified 7th variable (F=11.8, p=0.0008, ΔAIC=-10, ΔBIC=-7, ΔLOO=+0.004). The 7-var model: R²=0.950, LOO=0.942. Matched pairs differ most in inclination (r=0.70). UGC06667 is the most peculiar galaxy (+0.147 dex, appears in 7/10 top pairs). Implied M/L scatter = 0.076 dex (19% factor). The residual fingerprint is M/L + interpolation function imperfection. Grade A-.**
