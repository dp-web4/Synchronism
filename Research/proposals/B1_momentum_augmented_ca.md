# Directive: Pursue Option B1 — Momentum-Augmented CA

**Author**: Dennis Palatov
**Date**: 2026-03-20
**Status**: ACTIVE — pursue immediately
**Context**: Session #20 confirmed saturation wall formation but zero oscillations due to missing momentum

---

## Decision

**Pursue B1.** Add a momentum/velocity term to the CA transfer rule.

Session #20 identified the exact gap: the discrete CA has no memory of past state. Energy reaches saturation walls and thermalizes because there's no inertial term to cause reflection. This is the missing physics.

---

## The Modified Transfer Rule

**Current (first-order):**
```
ΔI_t = k · ∇²I · R(I)
```

**Proposed (second-order, momentum-augmented):**
```
v_t = α · v_{t-1} + k · ∇²I · R(I)
I_{t+1} = I_t + v_t
```

Where:
- `v` = velocity field (new state variable — Intent flow rate)
- `α` = momentum retention factor ∈ (0, 1) — how much inertia carries forward
- `α = 0` recovers the original first-order CA (no momentum)
- `α → 1` approaches inviscid flow (perpetual momentum)

This is the simplest possible second-order extension. It adds exactly one thing: memory of direction.

---

## Why This Brings N-S Back

Session #11 dismissed the N-S mapping as "vocabulary not physics" because the first-order CA was scalar nonlinear diffusion (1 DOF). It had no velocity field, so calling it N-S was a category error.

A second-order CA with an explicit velocity field `v` is structurally different. The momentum term `α · v_{t-1}` is the discrete analog of the inertial term `ρ(∂v/∂t + v·∇v)` in Navier-Stokes. The N-S mapping is no longer vocabulary — it's the actual physics being modeled.

The key distinction: we're not mapping the CA onto N-S by analogy. We're adding the specific physical mechanism (inertia) that N-S captures and the first-order CA lacked. The motivation is bottom-up (Session #20 showed we need reflection) not top-down (let's make the CA look like N-S).

---

## Test Plan

### Phase 1: Does momentum enable reflection at saturation walls?

Use Session #20's wall-forming configurations (A/I_max ≥ 1.2, low k) with momentum added.

Parameter sweep:
- `α` ∈ [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]
- Reuse Session #20's best wall-forming configs
- ~50-70 configurations

**Success criterion**: Energy reflects off saturation walls (velocity reversal at R→0 boundary). Visual: energy pulse bouncing back instead of thermalizing.

### Phase 2: Do reflections produce standing waves?

For configurations that show reflection:
- Run for 10,000+ timesteps
- Measure: autocorrelation of cavity energy (same as Session #20)
- Measure: spatial Fourier modes inside cavity
- Look for: f ∝ 1/L relationship (frequency determined by wall separation)

**Success criterion**: Periodic oscillation with autocorrelation peak > 0.5 at nonzero lag. Frequency correlating with cavity geometry.

### Phase 3: Connect to entity criterion

If standing waves form:
- Measure wall leakage rate (energy escaping per cycle) → this is Γ
- Measure total cavity energy → this is m
- Test: do stable cavities satisfy Γ < m?
- Test: do unstable/decaying cavities have Γ > m?

**Success criterion**: Entity criterion (Γ < m) emerges from the dynamics rather than being postulated.

---

## What to Watch For

1. **α sensitivity**: There should be a critical α below which thermalization still wins. Finding that threshold tells us the minimum inertia needed for entity formation.

2. **Wall sharpness matters**: Sharper walls (higher n in R(I) = [1-(I/I_max)^n]) should give cleaner reflections. There may be a critical n.

3. **2D behavior**: If 1D works, 2D cavities become possible. Circular/spherical standing waves would be the first real "particles" from the substrate.

4. **Energy conservation**: Verify total energy is conserved (or track dissipation explicitly). The momentum term shouldn't create or destroy energy.

5. **N-S correspondence**: With velocity field `v` and saturation `R(I)`, check whether the discrete dynamics reproduce known fluid phenomena (vortices, turbulence onset, Reynolds-like transitions). If they do, the N-S mapping is rehabilitated as physics, not vocabulary.

---

## Relationship to Prior Sessions

| Session | Finding | Status after B1 |
|---------|---------|----------------|
| #11 | N-S mapping is vocabulary | **Reopened** — second-order CA has real velocity field |
| #13 | Transfer rule can't produce entities | **Reopened** — that was first-order only |
| #18 | Pure diffusion: 0/324 oscillations | **Superseded** — first-order, no momentum |
| #19 | Reactive-diffusion: 0/300 oscillations | **Superseded** — first-order, no momentum |
| #20 | Walls form, 0 oscillations, no momentum | **Foundation** — identified the exact gap |

If B1 succeeds, the research arc becomes:
1. First-order CA → confinement but no oscillation (Sessions #18-20)
2. Second-order CA (momentum) → confinement + reflection + standing waves
3. Entity criterion derived from wall leakage vs cavity energy
4. N-S mapping rehabilitated as the continuum limit of the second-order CA

---

*The first-order CA was a diffusion equation missing its most important term. Time to add it back.*
