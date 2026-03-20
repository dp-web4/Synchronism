# Thor Session #20: Saturation Wall Test - Partial Success

**Date**: 2026-03-19
**Machine**: Thor autonomous session
**Proposal**: `Research/proposals/saturation_wall_oscillation_hypothesis.md`

---

## Summary

**Operator's Hypothesis**: R(I) saturation creates reflective walls that confine energy into oscillating cavities (geometric confinement, not amplitude instability).

**Test**: 72 configurations, 5000 timesteps each, systematic parameter sweep

**Results**:
- ✅ **45.8% produced stable saturation walls**
- ❌ **0% produced oscillations**
- ❌ **All autocorrelation peaks = 0.000** (complete thermalization)

**Verdict**: **Geometric confinement works, but no oscillation**

---

## Key Findings

### 1. Wall Formation CONFIRMED ✅

**Mechanism validated:**
- Localized energy pulse drives I → I_max in surrounding cells
- R(I) → 0 blocks energy transfer (ΔI = k·∇²I·R(I) → 0)
- Stable confinement boundaries form
- Wall separation: 17.8 ± 13.7 cells

**Parameter dependence:**
- A/I_max > 1.0 required (pulse must exceed saturation threshold)
- A/I_max = 1.2: 83% wall formation
- A/I_max = 2.0: 100% wall formation
- Lower coupling k favors sharper walls

### 2. Energy Thermalization (Not Oscillation) ❌

**Expected:** Energy bounces between walls → standing wave → oscillation

**Observed:** Energy enters cavity → completely thermalizes → no oscillation

**Root cause:** Discrete CA lacks momentum
- Continuum N-S: flow has inertia (ρv term), can reflect off boundaries
- Discrete CA: ΔI = f(I_current), no velocity history
- Energy stops at walls, doesn't "bounce back"

---

## Comparison Across Sessions

| Session | Mechanism | Tested | Oscillations | Walls |
|---------|-----------|--------|-------------|-------|
| #18 | Pure diffusion | 324 configs | 0 (0%) | Not measured |
| #19 | Reactive-diffusion | 300 configs | 0 (0%) | Not measured |
| #20 | Geometric confinement | 72 configs | 0 (0%) | 33 (45.8%) |
| **Total** | **All substrate mechanisms** | **696 configs** | **0 (0%)** | **33 (4.7%)** |

**Combined conclusion**: Substrate dynamics ΔI = k·∇²I·R(I) produces confinement but not oscillation.

---

## Implications

### Operator's Hypothesis: Partially Correct

**What was right:**
1. Sessions #18-19 tested wrong observables (amplitude at points, not spatial confinement)
2. R(I) → 0 DOES create reflective boundaries
3. Geometric confinement IS achievable

**What was incomplete:**
1. Confinement alone ≠ oscillation
2. Need energy reflection, not just blocking
3. Discrete substrate lacks momentum for "bounce"

### Three-Way Fork for Option B

The operator identified the gap. Three paths forward:

**Option B1: Add momentum to transfer rule**
```
ΔI_t = ΔI_{t-1} + k·∇²I·R(I)
```
- Pro: Enables reflection via inertia
- Con: Fundamentally changes CA structure
- Con: Converges to continuum N-S (undermines "discrete substrate")

**Option B2: Cavity-specific reactive terms**
```
Stage 1: Pure diffusion creates walls
Stage 2: Reactive term ε·f(I) active INSIDE cavity only
```
- Pro: Keeps discrete CA structure
- Con: Requires wall detection → multi-stage rule
- Con: Session #19 showed reactive terms fail with R(I) damping

**Option C: Accept oscillation as axiom**
- Pro: Honest about substrate limitations
- Pro: Entity criterion remains valid axiomatic prediction
- Pro: Matches current evidence (0/696 oscillations)

---

## Technical Details

### Wall Detection Algorithm

Walls identified where R(I) < 0.01 for >100 consecutive timesteps.

**Stability requirement**: 80% of last 100 steps must have 2+ walls at consistent positions.

**Wall separation**: Distance between outermost detected walls.

### Oscillation Detection

Autocorrelation on cavity energy time series:
```python
E_cavity(t) = Σ I[x] for x between walls
autocorr(lag) = correlation(E_cavity, shifted by lag)
```

**Threshold**: autocorr peak > 0.5 at lag ∈ [5, 500]

**Result**: ALL peaks = 0.000 (zero correlation at all lags)

### Thermalization Timescale

No transient oscillation detected. Energy thermalizes before first measurement (< 100 timesteps).

---

## Recommendations

1. **Update Session Primer**: Add Session #20 findings to stress test arc

2. **Present three-way fork**: Community decision on B1/B2/C

3. **Next computational work depends on decision**:
   - If B1: Test momentum-augmented CA
   - If B2: Test cavity-reactive hybrid
   - If C: Close computational arc, focus on entity criterion validation

4. **Acknowledge partial success**: Operator's geometric confinement insight is valuable even though oscillation didn't emerge

---

## Conclusion

**Finding**: R(I) saturation creates stable confinement walls (validated) but discrete CA without momentum cannot produce oscillation via energy reflection (newly discovered limitation).

**Status**: Option B partially works. Geometric confinement is real. Oscillation requires additional mechanism (momentum) or axiomatic postulation.

**Cumulative evidence** (Sessions #18-20): 0/696 configurations produced oscillations from substrate dynamics.

**Decision point**: Add momentum (B1), try cavity-reactive (B2), or accept axiom (C)?

---

*This completes the computational test of operator's saturation wall hypothesis. The hypothesis was partially correct - confinement works, but oscillation requires more than confinement alone.*
