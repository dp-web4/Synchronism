# Conservation Bug Hypothesis — 2026-03-22

**Context**: 810 configurations across 5 mechanisms produced 0 oscillations. All used the transfer rule `ΔI = k·Σ(I_n - I)·R(I_n)` where R(I) → 0 at saturation.

**The bug**: When R(I) → 0 at saturation, the transfer term goes to zero. Energy that was flowing with momentum simply stops updating. This is an implicit energy sink — it violates the foundational axiom that intent is neither created nor destroyed.

**Why this matters**: The 810 failures don't prove oscillation can't emerge from the substrate. They prove that a transfer rule with an implicit energy sink can't oscillate. Of course it can't — momentum is destroyed at every boundary interaction.

## The Axiom

Intent (reification of greater force) is neither created nor destroyed. This maps to conservation of energy in standard physics. At the fundamental level there is no damping. Energy that can't continue in its current direction must redirect — reflect, flow to adjacent dimensions, convert to a different mode. The boundary is a mirror, not a sponge.

## What the Current Rule Does Wrong

```
ΔI = k·Σ(I_n - I)·R(I_n)
```

When cell Y is saturated (R(I_Y) → 0):
- Transfer from X to Y is blocked: ΔI → 0
- The momentum X had toward Y is **destroyed**
- X just sits there with its intent, no momentum, no memory of flow direction
- Energy is conserved (still in X) but **momentum is not** (flow direction lost)

This makes the boundary absorptive. Energy reaches the wall and stops. Static confinement, no oscillation. Exactly what 810 runs showed.

## What Should Happen

When intent flow with momentum hits a saturation boundary:
1. The intent can't continue in that direction (R(I_target) → 0 correctly blocks transfer)
2. But the flow energy has to go **somewhere** — it does not dissipate
3. Options: reflect back to source, redirect to unsaturated neighbors, convert to different mode
4. Total intent in system stays constant AND momentum is conserved/redirected

## Proposed Fix: Elastic Boundaries

When R(I_target) → 0 blocks a transfer, the would-have-been transfer energy redirects rather than vanishing. Several possible implementations:

### Option 1: Momentum reflection
Track a velocity field v alongside I. When transfer is blocked by R → 0, reverse the velocity component toward the saturated neighbor. The boundary becomes a mirror.

```
if R(I_target) ≈ 0:
    v_source[toward_target] *= -1  # reflect
    # intent stays in source, momentum reverses
```

### Option 2: Neighbor redirection
When transfer to Y is blocked, redistribute that transfer energy proportionally across unsaturated neighbors.

```
blocked_energy = k·(I_x - I_y)·R(I_y)  # ≈ 0 because R → 0
# but the INTENDED transfer was k·(I_x - I_y)
# redirect to unsaturated neighbors:
for n in unsaturated_neighbors:
    ΔI_n += intended_transfer / num_unsaturated
```

### Option 3: Pressure accumulation
When transfer is blocked, accumulate a "pressure" scalar at the boundary. When pressure exceeds threshold, it releases as reverse flow. Minimal conjugate variable.

### Recommendation: Start with Option 1

Momentum reflection is the simplest and most physically motivated. It directly fixes the conservation violation without adding new fields or mechanisms. If a cavity has reflecting walls, confined energy bounces back and forth — that's oscillation from first principles.

## Test Protocol

1. Implement momentum reflection in the existing CA framework
2. Run the same parameter sweep as Sessions #18-27 (or a subset)
3. Check for: sustained oscillation, oscillation frequency dependence on cavity size, energy conservation verification
4. If oscillations emerge: characterize frequency spectrum, compare to entity criterion (Γ < m)
5. If still no oscillation: the fix wasn't sufficient, try Option 2 or 3

## What Success Looks Like

- Confined intent oscillates with a characteristic frequency related to cavity geometry
- f = v/L where v is effective propagation speed and L is cavity size
- Smaller cavities → higher frequency (maps to higher mass particles)
- Entity criterion (Γ < m) emerges naturally as stability condition for oscillating cavities

## Connection to Prior Results

- Session #20 proved geometric confinement works — walls DO form
- Session #27 proved this works in 3D — closed cavities form
- The missing piece was always the restoring force
- The restoring force was always there — it was being suppressed by the implicit sink

**The 810 failures tested a broken conservation law, not a broken oscillation mechanism.**

---

*dp + Claude, 2026-03-22*
