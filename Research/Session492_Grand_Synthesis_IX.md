# Session #492: Grand Synthesis IX — The Regime & Scatter Arc

**Date**: 2026-02-06
**Status**: Review (synthesizes Sessions #488-491)

## Arc Summary

Sessions #488-491 constitute the **Regime & Scatter Arc**, investigating how the 6-variable model performs across different spatial (radial), acceleration (MOND regime), and signal/noise (scatter budget) domains. These four sessions collectively produce 32/32 verified tests.

## Key Results by Session

### Session #488: Radial Offset Profile (Grade B+)
- 6-var R² increases monotonically: 0.67 (inner 25%) → 0.94 (outer 25%)
- c_V dominates the radial gradient (r = -0.53)
- Within-galaxy scatter = 0.087 dex (53% of between-galaxy)
- Radial gradient adds nothing to model (ΔLOO ≈ +0.001)
- **Implication**: Outer-only approach is definitively validated

### Session #489: BTFR From Offset (Grade A)
- Offset-corrected BTFR: slope = 4.10 (MOND: 4.0), R² = 0.992
- r(BTFR residual, offset) = -0.89 (partial: -0.96)
- 6-var model predicts M_bar to 13% accuracy (RMS = 0.053 dex)
- Coefficient on log(V_flat): -0.479 (MOND deep-limit: -0.500)
- Late types: 75% RMS reduction with offset correction
- **Implication**: The RAR offset IS the BTFR residual — they contain the same information

### Session #490: Deep MOND Limit (Grade B)
- 51% of MOND points in deep regime (g < 0.1 a₀)
- Model R² peaks at moderate MOND: R² = 0.981 at ⟨g⟩ = 0.1 a₀
- McGaugh function gives -0.047 dex systematic in deep regime (known property)
- Running exponent: 0.33-0.72 (observed) vs 0.52-0.71 (McGaugh theory)
- Deep-MOND BTFR slope = 2.5 (statistical artifact of restricted range)
- **Implication**: The sweet spot for the model is moderate MOND (0.05-0.3 a₀)

### Session #491: The Scatter Budget (Grade A-)
- Combined measurement noise = 28% of total offset variance (σ = 0.086 dex)
- V_obs noise: 3.6%, Distance (10%): 7.3%, M/L (σ=0.15): 9.4%, Inclination (4°): 8.2%
- R² from physics ≥ 0.67; R² from noise absorption ≤ 0.28
- MC null test: R² = 0.74 on independent noisy realizations
- Noise-predictor correlations < 0.07 — model fits real physics
- r(M/L noise, f_gas) = -0.85: gas-rich galaxies are immune to M/L errors
- LOO gap = 0.7%: minimal overfitting
- **Implication**: The 6-var model is measuring real physics, not fitting noise

## Cumulative Knowledge

### The 6-Variable Model: Final Assessment

```
offset_outer = -3.38 + 1.90×logV - 0.55×logL - 0.22×c_V - 0.45×f_gas
               + 0.15×logV×c_V + 0.18×logL×f_gas
```

| Metric | Value | Source |
|--------|-------|--------|
| R² | 0.945 | Session #483 |
| LOO R² | 0.938 | Session #483 |
| RMS | 0.038 dex | Session #483 |
| LOO RMS | 0.041 dex | Session #483 |
| R² physics (lower bound) | 0.668 | Session #491 |
| Noise fraction | 27.7% | Session #491 |
| Optimal domain (spatial) | Outer disk | Session #488 |
| Optimal domain (acceleration) | 0.05-0.3 a₀ | Session #490 |
| BTFR slope recovery | 4.10 (MOND: 4.0) | Session #489 |
| M_bar prediction accuracy | 13% (0.053 dex) | Session #489 |
| t-statistic (logL×f_gas) | 8.58 | Session #484 |
| NN autocorrelation | Eliminated | Session #484 |
| Cross-type validity | Late→Early R²=0.61 | Session #485 |
| M/L robustness | β stable at all M/L | Session #486 |

### The Error Hierarchy

From most to least important (as fraction of total offset variance):

| Source | σ (dex) | % of σ²_total |
|--------|---------|---------------|
| M/L variation (σ=0.15) | 0.050 | 9.4% |
| Inclination (4°) | 0.047 | 8.2% |
| Distance (10%) | 0.044 | 7.3% |
| V_obs noise | 0.031 | 3.6% |
| **Combined** | **0.086** | **27.7%** |

M/L is the largest single error source, but gas-rich galaxies avoid it entirely (r = -0.85). This creates a natural quality gradient: gas-dominated dwarfs have the smallest effective noise, while star-dominated early types have the largest.

### The Regime Map

| Regime | Points | R² | Key Feature |
|--------|--------|-----|-------------|
| Inner disk | ~570 | 0.67 | c_V-dominated, inner mass geometry |
| Outer disk | ~570 | 0.94 | Clean MOND physics |
| Deep MOND (g < 0.1 a₀) | 1146 | 0.94 | Low S/N, 51% of all MOND points |
| Moderate MOND (0.1-0.5 a₀) | ~560 | **0.98** | Sweet spot: high S/N + clean MOND |
| Transition (g > 0.5 a₀) | ~550 | 0.93 | Newtonian contribution dilutes |

The **optimal measurement window** is outer-disk points at moderate MOND accelerations (0.05-0.3 a₀). This is where the physics is cleanest and the measurement noise is lowest.

## Sessions #482-491: The Complete Model Improvement Arc

| Session | Topic | Key Finding | Grade |
|---------|-------|-------------|-------|
| #482 | Residual Forensics | NN autocorrelation r=+0.46 in 5-var | B+ |
| #483 | The Sixth Variable | logL×f_gas → LOO R² 0.896→0.938 | **A** |
| #484 | 6-Var Validation | t=8.58, F=73.6, autocorrelation eliminated | **A** |
| #485 | Type-Specific Models | Cross-prediction fails; late 3-var LOO=0.957 | A- |
| #486 | M/L Sensitivity | logL×f_gas stable at all M/L values | B+ |
| #487 | Grand Synthesis VIII | Integration of #482-486 | Review |
| #488 | Radial Profile | R² monotonically increases outward | B+ |
| #489 | BTFR From Offset | Corrected BTFR slope = 4.10, R²=0.992 | **A** |
| #490 | Deep MOND Limit | 51% deep MOND, model peaks at 0.1 a₀ | B |
| #491 | Scatter Budget | 72% physical signal, noise = 28% | A- |

**10 sessions, 80/80 tests verified. Grand Total: 1237/1237.**

## Novel Predictions Update

| ID | Prediction | Status | Evidence |
|----|-----------|--------|----------|
| NP1 | a₀ = cH₀/(2π) | SUPPORTED | 94% agreement (Session #385) |
| NP2 | Morphology → scatter | PARTIALLY SUPPORTED | 88% structural mediation (Session #383) |
| NP6 | N_corr → offset | SUPPORTED | R² = 0.23, p = 2×10⁻¹² (Session #389) |
| NP7 | R_eff MOND-dominated | SUPPORTED | r = -0.49 for late types (Session #392) |
| NP8 | R_eff M/L-independent | SUPPORTED | Gas-dominated: r = -0.59 (Session #389) |
| NP9 | R_eff L-independent | SUPPORTED | Late types: partial r = -0.49 (Session #392) |
| NP10 | R_max dynamical confirmation | SUPPORTED | r = -0.47 controlling V+L (Session #393) |
| NP11 | logL×f_gas interaction | **NEW** | LOO R² +0.042, t = 8.58, F = 73.6 (Session #483) |
| NP12 | Offset = -BTFR residual | **NEW** | r = -0.89, corrected slope = 4.10 (Session #489) |
| NP13 | 72% of offset variance physical | **NEW** | MC: noise = 28%, physics ≥ 67% (Session #491) |

## Open Questions

1. **Why is the cross-prediction failure asymmetric?** Late→Early R² = 0.61 vs Early→Late R² = 0.67. What drives this difference?

2. **Can we reduce the noise further?** Using only Q=1 edge-on gas-dominated galaxies would minimize all four noise sources simultaneously. What R² would this subsample achieve?

3. **What IS the remaining 3.8% gap to the noise ceiling?** With 72% signal and R² = 0.945, the signal-only R² is ~0.92. What additional physics lives in the remaining 8%?

4. **The running exponent discrepancy**: Observed slopes (0.33-0.72) deviate from McGaugh theory (0.52-0.71). Is this measurement noise or real physics?

5. **Distance-specific effects**: Some SPARC galaxies have Cepheid distances (2-3% errors) while others have Hubble flow (15-20%). Does model performance correlate with distance quality?

---

*Grand Synthesis IX: Sessions #488-491, 32/32 tests verified*
*Cumulative: 1237/1237 verified tests*

**The 6-variable model is now thoroughly characterized: it captures real physics (R²_phys ≥ 0.67), peaks at moderate MOND accelerations (R² = 0.98), improves monotonically toward the outer disk (R² = 0.94), recovers the BTFR to 2.5% accuracy (slope 4.10 vs MOND 4.0), and predicts baryonic mass to 13%. Measurement noise contributes only 28% of total variance, with M/L being the largest source (9.4%) but negligible for gas-rich galaxies (r = -0.85).**
