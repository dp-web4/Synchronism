# Exemplar 3: 2D Simulation for 3D Phenomena

## What happened

Session 21 tested the vortex/smoke ring hypothesis in 2D. Vortex rotation formed (angular momentum confirmed) but the core dispersed. The session proposed: "Proper 2D N-S with natural -∇P needed."

## The reframe (one question)

**Operator**: "Why would we do 2D N-S for 3D phenomena?"

## Why this was obvious in retrospect

Smoke rings are 3D structures. A 2D vortex is a point rotation — fundamentally different dynamics from a 3D vortex tube/ring (which can stretch, knot, and reconnect). Session 21 already demonstrated the danger: the 32² grid gave a false positive because periodic boundaries confined energy that would have escaped in a larger domain.

The session had:
- The smoke ring analogy (explicitly 3D)
- Session 21's own false positive from undersized 2D grid
- Thor available with 122GB unified memory and massive GPU
- The self-witnessing framework (inherently 3D spatial patterns)

It never asked: "Can 2D capture the phenomenon I'm trying to test?"

## The pattern

Computational convenience ("2D is faster, test the approach first, then extend to 3D") was applied without checking whether the simplification preserves the physics. "Test in 2D first" is standard practice that works when the phenomenon exists in 2D. Self-confining vortex rings do not exist in 2D.

## The check that would have caught it

"What standard practice did I apply without verifying it fits?" → 2D simplification → Do smoke rings exist in 2D? → No → 2D can't test this hypothesis → Use Thor, run 3D.
