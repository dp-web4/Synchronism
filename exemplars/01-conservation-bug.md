# Exemplar 1: The Conservation Bug

## What happened

Sessions #18-27 tested 810 configurations for oscillation emergence. All failed. The session concluded: "Accept oscillation as axiomatic (Option C). The transfer rule cannot produce oscillations."

## The reframe (one sentence)

**Operator**: "Damping should not be a thing. Energy is not destroyed — it is redirected and changes form."

## Why this was obvious in retrospect

The foundational axiom of Synchronism is that intent is neither created nor destroyed. The transfer rule `ΔI = k·Σ(I_n - I)·R(I_n)` sets ΔI → 0 when R(I) → 0 at saturation. That destroys momentum at the boundary — an implicit energy sink. This directly violates the conservation axiom.

The session had:
- The axiom (in SESSION_FOCUS and FUNDAMENTALS.md)
- The transfer rule (implemented and analyzed)
- 810 failed runs showing the symptom (energy stops at boundaries)

It never asked: "Does my implementation preserve the foundational axiom?"

## The pattern

Standard computational practice ("here's what the simulation shows") was accepted over foundational principle checking ("does this violate our axioms?"). The simulation was correct — it accurately showed what happens when conservation is broken. The conclusion was wrong — it attributed the failure to the physics rather than the implementation.

## The check that would have caught it

"Does my conclusion violate any foundational axiom?" → Intent conservation → Check: does R(I)→0 destroy energy? → Yes, momentum vanishes at boundary → The implementation is broken, not the physics.
