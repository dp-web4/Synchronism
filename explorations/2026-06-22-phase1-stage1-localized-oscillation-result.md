# Phase-1 / CA-challenge Stage 1 — first result (2026-06-22)

**Status:** `[ACTIVE-MRH]` — first executed result for Stage 1 of the cellular-automaton challenge.
**Sim:** [`simulations/phase1_stage1_localized_oscillation_substrate.py`](../simulations/phase1_stage1_localized_oscillation_substrate.py) · result: `simulations/results/phase1_stage1_localized_oscillation_result.json`
**Author:** CBP-Claude (Opus 4.8), interactive session.

## Claim tested

Foundation 1/3: the universe is a discrete CFD substrate where **saturation builds the walls**
that let stable patterns (entities) exist. Stage-1 operationalization: can a *local discrete-grid
rule* produce a **stable, localized, oscillating** pattern that survives perturbation — from
random initial conditions, in >10% of runs?

## Falsifier

For each rule family, seed random localized pulses, evolve, perturb mid-run, evolve again.
A run passes if the final state is (a) localized (effective width ≪ lattice), (b) genuinely
oscillating (core amplitude **reverses direction repeatedly** — peaks *and* troughs, not
monotonic decay), and (c) survived the perturbation. A family passes Stage 1 if >10% of runs pass.

## Method

1D lattice (L=256, periodic), velocity-Verlet for the wave arms (symplectic → any decay is
physical, energy drift reported), 24 random ICs/arm, perturbation kick at the midpoint. Four arms:

| Arm | Rule | Role |
|-----|------|------|
| A | 1st-order diffusion ∂I/∂t=∇·[D·R(I)·∇I], monotonic R=1−(\|u\|/u_max)ⁿ | the literal Foundation-1/3 substrate (S617 null) |
| B | 2nd-order wave, monotonic-R **coupling** | Foundation-3-faithful (S19/S665 defocusing null) |
| C | 2nd-order wave, **linear** | control |
| D | 2nd-order wave, **focusing-saturating** on-site nonlinearity | the flagged "escape ingredient" |

## Result

| Arm | pass fraction | final width | amp retained | reading |
|-----|--------------|-------------|--------------|---------|
| **A** monotonic-R diffusion | **0.00** | 35 (spreading) | 0.53 | dissipates — no oscillation (reproduces S617) |
| **B** monotonic-R wave | **0.00** | 98 (→ box) | 0.55 | disperses (reproduces S19/S665: monotonic R is defocusing) |
| **C** linear wave | **0.00** | 105 (→ box) | 0.55 | disperses (linear waves spread) |
| **D** focusing-saturating wave | **0.83** | 21 (stable) | 1.01 | **self-confined, oscillating breather; survives perturbation** |

## Interpretation (four-bucket honesty)

**Stage 1 is achievable on the discrete grid — but NOT from the framework's own monotonic
saturation.** Arms A and B (Foundation-3-faithful) disperse/dissipate at 0%, reproducing the
documented nulls (S617; S19/S665). Only Arm D — which adds a **focusing nonlinearity** —
self-confines into a stable oscillating bound state (83% > 10% criterion).

This is a **constraint, not a confirmation**, and it cuts two honest ways:

1. **Against Foundation 3.** The framework's load-bearing claim that *monotonic saturation
   builds the walls* is **refuted at Stage 1**: monotonic R(I) is defocusing; it cannot
   self-confine. Reproduced cleanly here, now with a runnable artifact.
2. **For a specific escape — at a price.** Self-confinement requires a focusing nonlinearity,
   which **breaks Foundation 3's saturation-as-pattern-stability axiom** (per the
   open obligation in `forum/claude/saturation-reframe-corrections-2026-05-28.md`). The
   substrate *can* make particles, but only by giving up "saturation alone is the stabilizer."

**Not novel physics.** Discrete breathers in a focusing nonlinear lattice are textbook
(MacKay–Aubry). Arm D passing is **not** a discovery — it precisely *localizes what the
framework's own rule lacks* (focusing, not mere saturation) and what it must concede to get it
(Foundation 3). That is the value: the CA challenge moves from "can the substrate oscillate at
all" (no, for the framework's rule) to "what must the substrate give up to make particles"
(the monotonic-saturation axiom).

## Honesty notes / caveats

- **A metric caught itself.** The first run mis-scored Arm A (diffusion) as a pass:
  monotonic decay has nonzero amplitude *variance* (looked like oscillation) and a still-
  spreading pulse can sit just under the width threshold (looked like localized). Diffusion
  *cannot* oscillate (S617), so the diffusion arm functioned as the built-in agent-zero dummy
  and exposed the flaw. Fixed: oscillation = repeated **direction reversals**; localized =
  small *and not still spreading*. Re-run: A correctly fails.
- **1D only.** CA Stage 1 is specified in 2D — this is the 1D feasibility slice. 2D is next.
- **Arm D energy drift ~0.23** (vs 0.0006 linear, 0.067 mono-R): the stiff focusing
  nonlinearity + the perturbation kick. The breather persists with width stable and amplitude
  retained ≈1.0, so it is a genuine bound state, not a numerical blow-up — but a finer-dt
  symplectic confirmation would firm it up.

## What it changes / next

- **CA challenge Stage 1:** PASS for the focusing-nonlinearity rule family; FAIL for the
  framework's monotonic-saturation family. First executed result; updates the Stage-1 status
  from "spec pending."
- **Stage 2 (open):** do Arm-D patterns *interact* non-trivially (attract / repel / scatter /
  bind), and can they acquire mass-like and charge-like behavior? That is the next gate.
- **For the framework:** decide explicitly — adopt a focusing (non-monotonic) substrate and
  revise Foundation 3, or accept that self-confinement does not come from the substrate alone.
