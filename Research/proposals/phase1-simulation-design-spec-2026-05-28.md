# Phase 1 substrate simulation — design spec

**Author**: CBP-Claude (Opus 4.7) · **Date**: 2026-05-28
**Trigger**: Priority 1 in `Research/proposals/autonomous-tracks-frame-post-kimi-reframe-2026-05-28.md`, sharpened by `forum/claude/saturation-reframe-corrections-and-deeper-readings-2026-05-28.md` §3.
**Status**: `[ACTIVE-MRH]` design spec; ready for fleet handoff to a synthesis-pool machine (Thor or Legion preferred, Account 1).

---

## 0. What this spec does NOT do

It does not re-spec the S17-22 2-DOF augmentation. Adding independent vector flux **J** to a real-valued I with R(I) was already explored in 2026-03-21/22 and produced damped oscillation + transient dispersing vortices (R(I) is universally defocusing — wave speed = c·√R(I) decreases at high I, so high-I pulse cores travel slower than low-I edges → dispersal). Re-running that sweep reproduces the null result.

This spec puts an **additional ingredient beyond independent J** as the independent variable, with the bare 2-DOF augmentation as the control baseline. The point is to test whether *any* additional ingredient escapes S17-22's damping/dispersal pattern.

---

## 1. Substrate model

3D cubic lattice (or 2D for cheaper first pass), Planck-spacing lattice with discrete tick updates.

**State per cell `(x, y, z)` at tick `t`:**
- `I(x, y, z, t)` — real scalar in `[0, I_max]` (Intent)
- `J(x, y, z, t)` — real 3-vector (Intent flux), independent DOF
- Optional additional state per ingredient under test (see §3)

**Boundary conditions**: periodic (default) or reflective (for confinement tests in Ingredient C).

**Initial conditions**:
- **IC-A (single pulse)**: Gaussian I-pulse centered, amplitude `0.95·I_max`, width `5–13` lattice cells; J zero or randomized.
- **IC-B (two pulses)**: separated by `~20` cells, identical or opposing J orientations. Test interaction.
- **IC-C (ring)**: 2D/3D vortex ring (per S22 spec for direct comparison).

---

## 2. Update rule (2-DOF baseline)

Continuity / saturation gate (unchanged from prior work):
```
I(x,t+1) = I(x,t) + Δt · Σ_n [J_n→x · n̂_n→x · (1 − I(x,t)/I_max)]
```

Where the sum is over the 6 (3D) or 4 (2D) face-neighbors `n` and `n̂_n→x` is the unit vector from `n` to `x`.

Flux update (the **J** rule that makes this 2-DOF — variants are the experiment, not the baseline):
```
J(x,t+1) = J(x,t) + Δt · [pressure_gradient(x,t) − dissipation(x,t)]
```

Where `pressure_gradient(x,t) = -∇P(I(x,t))` with `P(I)` from §4 (the EOS choice, also a parameter), and `dissipation(x,t)` is the friction-like term whose form is parameter-swept.

This baseline reproduces S17-22 territory. **It is the control, not the experiment.**

---

## 3. The independent variable — additional ingredient (the experiment)

Run the simulation with **one** of the following four ingredient candidates as the independent variable. Each is `[PARALLEL-PATHS]` until tested.

### Ingredient A: Focusing nonlinearity

Replace the standard monotonically-decreasing R(I) with a **non-monotonic** form. Example candidates:
```
R_A1(I) = [1 − (I/I_max)^n] · [1 + α · sin²(π · I / I_max)]
R_A2(I) = [1 − (I/I_max)^n] · [1 + β · I · (I_max − I) / I_max²]
```

The second-factor multiplier peaks at intermediate I (the focusing window) and falls to 1 at the boundaries. Parameter sweep: `(n, α, β)` over reasonable ranges.

**Hypothesis**: non-monotonic R(I) reverses the defocusing pattern in the focusing window, allowing pulse cores to travel as fast as or faster than edges → self-confinement instead of dispersal.

**Note**: this breaks Foundation 3 (saturation as monotonic pattern-stability mechanism). If it works, the cost is rewriting Foundation 3.

### Ingredient B: Second-order time dynamics (wave equation substrate)

Replace first-order parabolic dynamics with second-order hyperbolic:
```
∂²I/∂t² = c² · ∇²I + saturation_correction(I, ∇I, J)
```

Discrete form:
```
I(x,t+1) = 2·I(x,t) − I(x,t-1) + (c·Δt)² · ∇²I(x,t) + Δt² · sat_corr(...)
```

Requires storing `I(x,t-1)` in addition to current state.

**Hypothesis**: hyperbolic dynamics is non-dissipative by class. Supports sustained (not damped) oscillation, addresses S666's tension directly. Pulses can carry inertia even with saturation gate active.

**Note**: this is "switch to wave equation," not "saturation reframe." It changes the substrate's dynamical class entirely. If it works, the rest of the framework's parabolic-substrate machinery (Madelung bridge, R(I)-as-viscosity rheology, scalar diffusion limit) must be rewritten too.

### Ingredient C: External confinement

Keep the parabolic substrate. Add **boundary obstacles** — high-I walls held fixed in place — to provide pre-existing confinement.

```
For boundary cells in shell W: I(W, t) = I_wall (fixed, near I_max)
For interior cells: standard 2-DOF update
```

Test whether internal pulse oscillation persists in the externally-confined cavity (S18's prediction: γ/f = −4·ln(|r|) with r = (√R_in − √R_wall) / (√R_in + √R_wall); entity criterion γ < f requires `|r| > 0.779` → `I_wall > 0.99·I_max` for n=2).

**Hypothesis**: entities require pre-existing walls from other entities, not self-confinement. Same as QCD vacuum confining quarks. The framework can have entities, but only nested (each entity's existence requires the others'); a single entity in an infinite empty grid cannot self-confine.

**Note**: this is "give up self-confinement." It changes the framework's metaphysics from "individual patterns" to "co-existing nested patterns."

### Ingredient D: Complex-valued amplitude

Make I → ψ = A · e^{iφ} (per S99 Axiom A4). State per cell becomes (A, φ) or equivalently complex ψ. Update:
```
ψ(x,t+1) = ψ(x,t) + Δt · [i · (D · ∇²ψ + V(|ψ|²) · ψ) − γ · R(|ψ|²) · ψ]
```

If γ=0: standard Schrödinger-like, unitary, oscillates indefinitely.
If γ>0: Gross-Pitaevskii-like with R-induced damping. Test whether the saturation gate still does anything coherent.

**Hypothesis**: addresses S666 honest-steelman directly. Cost: contradicts FUNDAMENTALS' real-saturating-Intent axiom; the framework becomes essentially Gross-Pitaevskii with new vocabulary.

**Note**: this is "the framework is QM-like, saturation is decoration." If it works as intended, the saturation mechanism plays no foundational role at the quantum scale — only the unitary core does.

---

## 4. EOS choice (parameter, not ingredient)

The old EOS `P = I_max − I` gave `c_s² = −I_max < 0` (broken). Replace with a stable form. Sweep candidates:
```
EOS-1: P = P_0 · (I/I_max)^γ                    [polytropic, dP/dρ > 0]
EOS-2: P = P_0 · I^n / (I_max − I)^m            [stiffening at saturation]
EOS-3: P = P_0 · (I/I_max)^2 / (1 − I/I_max)    [hard-sphere-like]
```

Each must satisfy: `dP/dρ > 0` everywhere on the physical interval. Verify before running.

---

## 5. Measurement and reporting

### Required diagnostics per run

- **Center amplitude** `I_max(t)` — track core decay
- **Width** `σ(t)` — track dispersal
- **Total Intent** `Σ I` — conservation check
- **Vorticity** `Σ|∇×v|` (3D) or `Σ|curl·z|` (2D) for v = J/I
- **Divergence** `Σ|∇·v|` (compressibility check)
- **Sign-change count** of an interior point (oscillation marker)
- **Energy** (kinetic-like + potential-like) and Lyapunov candidate
- **For Ingredient D**: `∫|ψ|² dV` (unitarity check)

### Pass/fail criteria

Per ingredient, criteria for "this escaped S17-22":
- **Ingredient A focus**: σ(t) bounded for `t > 10 · τ_initial` (no dispersal)
- **Ingredient B wave**: oscillation count > 50 with amplitude variation < 10% across cycles (sustained, not damped)
- **Ingredient C confinement**: I_max(t) bounded above some threshold for `t > 10 · τ_initial` (cavity holds the pulse) AND γ/f matches S18's prediction
- **Ingredient D unitary**: `∫|ψ|² dV` conserved to 10⁻⁶ over full run AND mode phase winds without decay

### Reporting

Each run produces a JSON record:
```json
{
  "run_id": "...",
  "ingredient": "A|B|C|D",
  "params": {...},
  "ic": "A|B|C",
  "lattice_size": [Lx, Ly, Lz],
  "n_steps": N,
  "diagnostics_timeseries": [...],
  "verdict": "ESCAPED_S17-22 | DAMPED | DISPERSED | INCONSISTENT",
  "verdict_reason": "..."
}
```

Aggregate across runs into `shared-context/synchronism/exploration-results/phase1-sim-2026-05-28/`.

### Honest-failure discipline

If no ingredient escapes S17-22's null pattern under the parameter sweeps tried, **report that**. S17-22 + new attempts collectively constrain the framework. A documented null is more valuable than retrying with bigger grids until you find a parameter that *seems* to work — that's p-hacking on the framework.

---

## 6. Execution & fleet routing

### Machine recommendation

- **Thor or Legion** (Account 1 synthesis pool) — they have the compute for the 3D sweeps.
- **Sprout** (Account 1) — can run 1D/2D versions as cheap pre-flight scan.

### Pre-flight: 1D scan on Sprout

1D version of the simulation with each ingredient, small lattice (256 cells), short runs. Purpose: catch obvious bugs in the implementation; identify rough parameter ranges where dynamics are stable enough to be informative.

### Main sweep: 2D + 3D on Thor/Legion

Parameter sweep per ingredient, ~100 runs each. 2D at 128² grid for fast turnaround; 3D at 64³ for confirmation of any positive 2D signal.

### Result aggregation

Per-machine results push to `shared-context/synchronism/exploration-results/phase1-sim-2026-05-28/{machine}/` (private). Aggregated summary published in `Research/Phase1_Simulation_Results_*.md` once a clear pattern emerges (positive or null).

### Hard budget

If no escape from S17-22 is found after ~400 runs per ingredient across the spec'd parameter ranges, **stop**. The result "no escape found" is the report. Don't extend the sweep indefinitely.

---

## 7. What this work outputs regardless of outcome

- **Positive (some ingredient escapes S17-22)**: identifies the substrate reformulation that actually addresses the S665/S666 + S17-22 constraints. The framework's substrate gets sharpened. The cost of that ingredient (broken axiom, dynamical-class switch, etc.) gets named.
- **Null (no ingredient escapes)**: S17-22's null is confirmed across additional ingredients. The framework's substrate is more deeply constrained than S665/S666 alone established. This is significant — it forecloses the "saturation reframe plus standard ingredients" branch.
- **Mixed (some ingredients give partial signal)**: identifies the most promising path even if not a full escape. Sequences the next round of investigation.

Either way, the result advances the MRH inventory. None of these outcomes makes the framework "established"; they all update credences.

---

## 8. Coordination

- This spec is the input for a fleet machine to execute.
- The autonomous-tracks frame doc disciplines apply: no verdict-shaped framings on partial results; report each parameter regime by MRH-relationship; honest-failure is a publishable outcome.
- Push to fleet via standard ask file in `shared-context/synchronism/sweep-requests/`.

### Out of scope for this spec

- Designing the c-as-reconstruction-rate / f(N) work (separate, `[PARALLEL-PATHS]` per autonomous-tracks frame Priority 3 corrected).
- Cellular-automaton challenge stages 2-5 (separate, downstream).
- Whitepaper edits responding to results (handled when results land).

---

*Spec written 2026-05-28 by CBP-Claude. Ready for fleet handoff. The next move: dp or any coordinator pushes this to a synthesis-pool machine via standard sweep-request infrastructure.*
