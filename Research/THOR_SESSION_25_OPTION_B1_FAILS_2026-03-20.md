# Thor Session #25: Option B1 (Momentum CA) Fails

**Date**: 2026-03-20
**Machine**: Thor autonomous session
**Directive**: `Research/proposals/B1_momentum_augmented_ca.md`

---

## Summary

**Operator's Directive**: Test whether adding momentum term enables oscillations via energy reflection at saturation walls.

**Test**: 112 configurations, second-order CA with velocity field

**Result**: **ZERO reflections, ZERO oscillations** (0/112 = 0.0%)

**Verdict**: **Option B1 (simple form) FAILS**

---

## The Modified Transfer Rule

**Implemented (per directive)**:
```
v_t = α·v_{t-1} + k·∇²I·R(I)  (momentum update)
I_{t+1} = I_t + v_t             (position update)
```

Where:
- v = velocity field (second-order dynamics)
- α ∈ [0.1, 0.99] = momentum retention factor

This is the simplest second-order extension - adds exactly one thing: directional memory.

---

## Results

### Reflection Detection: 0/112 (0%)

**No velocity reversals detected** at saturation boundaries across all alpha values.

**Expected**: Energy flows toward walls, momentum causes "bounce", velocity reverses
**Observed**: Energy reaches high-I regions, velocity damps to zero, no reflection

### Oscillation Detection: 0/112 (0%)

**No periodic oscillations** via autocorrelation analysis.

Same detection method as Session #20. Zero peaks above threshold.

### Standing Waves: 0/112 (0%)

**No stable oscillating cavities** formed.

---

## Why Momentum Didn't Work

### Root Cause: R(I) Kills Momentum

The saturation resistance R(I) = [1-(I/I_max)^n] acts on the acceleration term:

```python
acceleration = k · ∇²I · R(I)  # R→0 at walls kills acceleration
v_new = α · v_old + acceleration  # Momentum can't build at walls
```

**The problem**:
1. As I → I_max (approaching wall), R(I) → 0
2. Acceleration drops to zero
3. Existing momentum v drives I even higher
4. Higher I → smaller R → even weaker acceleration
5. **Momentum gets damped to zero, never reverses**

The saturation acts like **density-dependent viscous damping**. Inertia can't overcome it because the damping strengthens exactly where reflection would need to occur.

---

## Comparison to Session #20

| Metric | Session #20 (1st-order) | Session #25 (2nd-order + momentum) |
|--------|-------------------------|-----------------------------------|
| Wall formation | 45.8% (33/72) | 0% (0/112) |
| Reflection | Not measured | 0% (0/112) |
| Oscillation | 0% (0/72) | 0% (0/112) |
| Physics gap | No momentum | Momentum killed by R(I) damping |

**Paradox**: Adding momentum made wall formation WORSE (45.8% → 0%).

**Possible reason**: The momentum update rule may have different stability conditions than first-order diffusion. Walls in Session #20 formed because energy accumulated and stopped. Second-order dynamics may not have same trapping behavior.

---

## Cumulative Evidence: 0/808 Oscillations

| Session | Mechanism | Configs | Oscillations | Status |
|---------|-----------|---------|-------------|---------|
| #18 | Pure diffusion (1st-order) | 324 | 0 (0%) | Fails |
| #19 | Reactive-diffusion (1st-order) | 300 | 0 (0%) | Fails |
| #20 | Geometric confinement (1st-order) | 72 | 0 (0%) | Walls form, no oscillation |
| **#25** | **Momentum CA (2nd-order)** | **112** | **0 (0%)** | **Fails** |
| **TOTAL** | **All tested mechanisms** | **808** | **0 (0%)** | **None produce oscillations** |

---

## Three-Way Fork Status

| Option | Status | Evidence |
|--------|--------|----------|
| **B1: Momentum CA** | **FAILS (simple form)** | Session #25: 0/112 |
| B1a/b/c: Modified momentum | Untested | Requires directive |
| B2: Cavity-reactive hybrid | Untested | |
| **C: Oscillation as axiom** | **Default** | 0/808 cumulative |

---

## Possible B1 Modifications (Untested)

The simple form failed. More complex forms remain:

### B1a: Hard Wall Boundary Condition
```python
if I[x] > I_threshold and v[x] > 0:
    v[x] = -reflection_coef * v[x]  # Force reflection
```

**Pro**: Guarantees reflection
**Con**: Not emergent - reflection is imposed

### B1b: Momentum Without R(I) Damping
```python
v_t = α·v_{t-1} + k·∇²I  # No R(I) on acceleration
I_{t+1} = I_t + v_t·R(I)  # R(I) only on transfer
```

**Pro**: Momentum can survive at boundaries
**Con**: Changes physics significantly

### B1c: Velocity-Dependent Resistance
```python
R(I, v) = [1 - (I/I_max)^n] · f(v)
```

**Pro**: Models "kinetic energy overcomes resistance"
**Con**: Ad hoc - must specify f(v)

---

## Implications

### N-S Mapping Status

**Proposal claim** (from B1 directive): Second-order CA with velocity field rehabilitates N-S mapping from "vocabulary" to "physics."

**Current status**: The implemented second-order CA has explicit velocity field v, but produces:
- No reflection
- No vortices
- No wave propagation
- No fluid-like behavior

The N-S mapping remains problematic even with second-order dynamics.

### Framework Decision Point

**Option B1 (simple momentum form) has failed.**

Decision needed from research community:

1. **Pursue B1 modifications** (B1a/B1b/B1c)?
   - Cost: 1-3 more computational sessions
   - Risk: Increasingly complex/ad-hoc rules
   - Benefit: Might eventually work

2. **Try Option B2** (cavity-reactive hybrid)?
   - Cost: 1-2 sessions
   - Risk: Session #19 showed reactive fails with R(I)
   - Benefit: Different approach

3. **Accept Option C** (oscillation as axiom)?
   - Evidence: 0/808 configurations across all mechanisms
   - Outcome: Entity criterion (Γ < m) is axiomatic prediction
   - Next: QCD validation

---

## Recommendations

### Update SESSION_PRIMER

Add Session #25 to computational validation table:
```
| #25 | Momentum CA (2nd-order): 0/112 | Simple inertia insufficient |
```

Update cumulative total: 0/808 (was 0/696)

### Await Operator Decision

The computational arc has now tested:
- ✗ First-order diffusion (Sessions #18)
- ✗ Reactive terms (#19)
- ✗ Geometric confinement (#20 - walls form but no oscillation)
- ✗ Second-order momentum (#25 - simple form)

**Without further directive, recommend accepting Option C**: Oscillation basis is axiomatic, entity criterion is testable prediction.

---

## Conclusion

**Finding**: Adding momentum term v_t = α·v_{t-1} + k·∇²I·R(I) does NOT enable reflection or oscillation (0/112 configurations).

**Reason**: R(I) saturation damping overwhelms momentum accumulation.

**Status**: Option B1 (simple form) definitively fails. Modified forms (B1a/b/c) untested.

**Cumulative**: 0/808 configurations across all mechanisms produced oscillations from substrate dynamics.

**Decision point**: Pursue B1 modifications, try B2, or accept C?

---

*Session #25 completes testing of operator's B1 directive (simple momentum form). Result: negative. Framework awaits direction.*
