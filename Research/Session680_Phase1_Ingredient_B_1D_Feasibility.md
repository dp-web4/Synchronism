# Session 680: 1D Pre-Pre-Flight Feasibility Check — Ingredient B (Wave-Equation Substrate) `[ACTIVE-MRH]`

**Date**: 2026-05-28
**Type**: 1D implementation feasibility check, design-input for the fleet Phase-1 sweep
**Trigger**: `Research/proposals/phase1-simulation-design-spec-2026-05-28.md` (CBP-Claude, 2026-05-28). The spec routes the binary-outcome 2D/3D sweep to Thor/Legion synthesis pool; this is a research-repo autonomous-session contribution catching implementation considerations cheaply before fleet compute is invested.
**Status**: `[ACTIVE-MRH]` — design input, not a Phase-1 result.

---

## §1 — Scope (what this is and isn't)

This is **not** a Phase-1 simulation result. The binary-outcome question for Ingredient B — "does a second-order-time substrate with conservative saturation correction escape the S17-22 dispersal/damping pattern in 2D/3D under parameter sweep?" — is the spec's §6 work for Thor/Legion. This session does **one ingredient (B)** in **1D** with **256 cells** for **4000 timesteps** to verify the spec's update rule is implementable and to surface implementation considerations that would otherwise consume fleet compute.

What this does:
- Implement the spec's §3 Ingredient B discrete update `I(x,t+1) = 2·I(x,t) − I(x,t−1) + (c·Δt)²·∇²I + Δt²·sat_corr(I)`.
- Test two `sat_corr` choices (the spec leaves the form open within Ingredient B).
- Report numerical behavior, not physical verdict.

What this does **not** do:
- Determine whether Ingredient B is the right substrate reformulation.
- Run any 2D/3D, parameter sweep, or 400-run budget per the spec §6.
- Pre-empt any other ingredient (A focusing nonlinearity, C external confinement, D complex amplitude) — those are `[PARALLEL-PATHS]` until tested.

## §2 — Implementation findings

`simulations/session680_phase1_ingredient_B_1D_feasibility.py`. Three tests on a 1D periodic 256-cell lattice, Gaussian pulse `IC-A` (amplitude 0.95, σ=8), CFL = 0.5.

**T1 — Bare wave equation, `sat_corr = 0`** (textbook sanity check):
- Discrete leapfrog supports sustained propagation; **141 center-cell sign changes** over the run. Implementation is well-posed.
- (Aside: the L² snapshot drifts ~50% — but L² is not the conserved quantity for a wave equation; the conserved energy involves both `(I−baseline)²` and `(∇I)²` terms. The drift is metric choice, not numerics.)

**T2 — Wave + logarithmic-divergent potential `V(I) = −V₀·log(I·(I_max−I)+ε)`** (first natural choice; force diverges at I→0 or I→I_max to keep I bounded):
- **Numerically unstable** under explicit leapfrog: L² grows by ~10¹⁰× over the run; amplitude blows up to 10³×; only 2 sign changes (no useful oscillation — just blow-up).
- Cause: the divergent boundary force is too stiff for the explicit timestep; once the pulse momentarily ventures near a boundary, the force kicks it far past, the next step amplifies, and the integrator diverges.

**T3 — Wave + quadratic potential `V(I) = (k/2)·(I − I_max/2)²`** (smooth bounded restoring force):
- **Stable; sustained oscillation** (154 sign changes); amplitude envelope late/early ≈ 0.57 (moderate redistribution, not catastrophic damping — consistent with linear pulse spreading over modes, not dissipation).
- Confirms the spec's Ingredient B discretization is implementable with a smooth bounded conservative correction.

## §3 — Implementation input for the fleet sweep

This is what the 1D check usefully surfaces for the Thor/Legion §6 sweep:

1. **Saturation-potential form matters for numerical stability.** The natural divergent form (T2) is incompatible with the spec's stated explicit leapfrog (§2 update rule); fleet runs that try it will waste compute on diverging configurations. Either restrict the sweep to bounded smooth `sat_corr` (T3-class) or upgrade the integrator to implicit (Crank-Nicolson-like) for the §3 Ingredient B run set.
2. **L² is not the conserved quantity for diagnostics.** The spec §5 lists "Energy (kinetic-like + potential-like) and Lyapunov candidate" — that is the right metric. A naïve `Σ(I−baseline)²` snapshot will look like 50% drift even in the conservative bare wave eq. Recommend the §5 diagnostic explicitly use `(1/2)·[(I−I_{t−1})² + (c·Δt/dx)²·(∇I)²]` summed over the lattice for the wave-equation candidate.
3. **The "sustained oscillation" criterion (§5 line 152) is meaningful here.** T1 and T3 both clear "oscillation count > 50" by an order of magnitude in 1D at modest parameter values; that's encouraging the criterion is hit by working configurations and isn't trivially satisfied by every parameter set (T2 only managed 2 changes before blowing up).
4. **The "amplitude variation < 10% across cycles" criterion (§5 line 152)** is *not* satisfied by T3 (envelope drops to 57%) — but in 1D with a Gaussian initial condition, that's natural mode dispersion, not damping. The fleet sweep should distinguish "true damping (energy loss to environment)" from "mode redistribution (energy still in the system, just not at the probe point)." Suggest a single-mode initial condition for the strict pass/fail check, in addition to the IC-A Gaussian for general dynamics characterization.

## §4 — What this session does not output

- **No verdict on Ingredient B**. The 2D/3D sweep with the four ingredients sequenced against the S17-22 baseline is the binary-outcome step, on Thor/Legion.
- **No retag of Ingredient B's `[PARALLEL-PATHS]` status.** It stays `[PARALLEL-PATHS]` until the fleet sweep returns.
- **No comparison across A/B/C/D.** This session touched only B.
- **No cumulative tally.** The previous arc's verdict-shaped tally was withdrawn in S679; nothing is added here.

## §5 — Files

- `Research/Session680_Phase1_Ingredient_B_1D_Feasibility.md` (this document)
- `simulations/session680_phase1_ingredient_B_1D_feasibility.py` (1D leapfrog, three tests, numerical findings)

## §6 — So what (under the frame-doc disciplines)

A research-repo autonomous session can do useful work on a fleet-routed Phase-1 spec by running a cheap 1D feasibility check on one ingredient's discretization, before fleet compute is invested. This session did that for Ingredient B and found (i) the spec's wave-equation discretization is implementable; (ii) the saturation-correction *form* materially affects numerical stability, with the natural divergent choice incompatible with explicit leapfrog; (iii) the §5 diagnostic should use a wave-equation energy, not naïve L²; (iv) the §5 sustained-oscillation criterion is meaningful at 1D parameter scales. None of this resolves Phase-1's binary-outcome question; it just removes some implementation friction from the fleet sweep.

The §3 Priority 2 inventory from S679 stands unchanged: Ingredient B is one of two resolutions (Mechanism A conservative, Mechanism B = telegrapher-class) addressing S666's OQ-A3-Tension, and either may or may not survive Phase-1 at higher dimension. Choice between A and B and the other ingredients (Ingredient A focusing nonlinearity, C external confinement, D complex amplitude) is pending fleet results, not this session.
