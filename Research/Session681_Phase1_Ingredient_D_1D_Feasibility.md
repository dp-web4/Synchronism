# Session 681: 1D Pre-Pre-Flight Feasibility Check — Ingredient D (Complex Amplitude / Gross-Pitaevskii Branch) `[ACTIVE-MRH]`

**Date**: 2026-05-31
**Type**: 1D implementation feasibility check, design-input for the fleet Phase-1 sweep, second in the per-ingredient pre-flight series (following S680 / Ingredient B).
**Trigger**: `Research/proposals/phase1-simulation-design-spec-amendment-2026-05-31.md` §3.1 line 109: *"**Ingredient D** — 1D needed. Implement Gross-Pitaevskii with optional R-induced damping (γ parameter); verify `∫|ψ|²` conservation when γ=0, characterize damping when γ>0."*
**Status**: `[ACTIVE-MRH]` — design input, not a Phase-1 result.

---

## §1 — Scope (what this is and isn't)

Same shape as S680. This is **not** a Phase-1 simulation result. The binary-outcome question for Ingredient D — "does a complex-amplitude substrate with saturation-induced damping escape the S17-22 dispersal/damping pattern, and at what cost to the saturation foundation?" — is the amendment's §3.2 work for Thor/Legion. This session does **one ingredient (D)** in **1D** with **256 cells / 2000 timesteps** to verify the spec's Ingredient D discretization is implementable and to characterize the qualitative damping behavior the spec asks for.

Following the per-ingredient leverage order in the amendment §1 (B-A1 > B-A2 > **D** > A > C), D is the highest-leverage remaining pre-flight after S680 closed B.

## §2 — Implementation findings

`simulations/session681_phase1_ingredient_D_1D_feasibility.py`. Three tests using the **split-step Fourier propagator** (half-step potential, full-step kinetic via FFT, half-step potential — exact for the free-particle step). Initial wavepacket: Gaussian envelope σ=0.05 with carrier k₀=30, normalized to ∫|ψ|² = 1. Lattice 256 cells, periodic BC, T=0.02 over 2000 steps.

**T1 — Free Schrödinger** (γ=0, V=0): textbook unitarity sanity. Norm drift over the run: **+9.6×10⁻¹⁴** (machine precision). Mode-k₀ phase winds monotonically through the run. Discretization is well-posed for Ingredient D's unitary core.

**T2 — Gross-Pitaevskii** (γ=0, V(|ψ|²) = α|ψ|² with α=5.0): norm drift **+9.0×10⁻¹⁴**. The nonlinear potential is real-valued so the generator stays anti-Hermitian and unitarity is preserved analytically; split-step Fourier preserves it numerically too. Confirms Ingredient D is implementable with the GPE-class nonlinearity the amendment specifies.

**T3 — Saturation-induced damping** (γ>0, `R(|ψ|²) = 1 − (|ψ|²/|ψ|_max²)^n`, n=2). Three regimes:

| γ | end-norm / start | qualitative damping |
|---|---:|---|
| 0.5 | 0.984 | gentle, ~1.6% loss |
| 5.0 | 0.851 | moderate, ~15% loss |
| 50.0 | 0.182 | strong, ~82% loss |

Damping rate scales smoothly with γ, roughly proportional to γ·⟨R(|ψ|²)⟩. Phase dynamics persist throughout (the wavepacket keeps oscillating while amplitude erodes) — qualitatively telegrapher-like, not pure exponential decay.

## §3 — Implementation input for the fleet sweep (amends S680 amendments)

These are direct integrations of S681 findings into the spec's §3/§5/§6 for Ingredient D specifically. They sit alongside the S680 amendments (which were about Ingredient B).

1. **Integrator: split-step Fourier (not explicit Euler).** Explicit forward-Euler on Schrödinger is unconditionally unstable — amplification factor |1 + iλΔt| > 1 for purely imaginary eigenvalues. S666's sim caught this. The spec for Ingredient D should explicitly call out split-step Fourier (or Crank-Nicolson for non-periodic BC) and **not** leave it as an open implementation choice; the fleet sweep should not waste compute reproducing the Euler instability.
2. **Strict pass criterion at γ=0**: `∫|ψ|² dV` conserved to 10⁻⁶ over the full run AND mode-k phase winds linearly in t with slope `D·k²`. The second clause verifies both the unitary structure *and* the kinetic operator's coefficient are right (catches sign errors and factor-of-2 mistakes in `D` that the norm check alone misses).
3. **Characterization at γ>0**: report *damping rate* `d(∫|ψ|²)/dt` averaged over the run as the diagnostic, not pass/fail. Also report **whether phase dynamics persists while amplitude damps** (telegrapher-class behavior vs. pure exponential decay vs. fragmentation). The qualitative shape is informative for downstream classification — telegrapher-class with phase preservation is closer to the "saturation reframe of standard QM" reading; pure exponential decay is closer to "saturation kills coherence."
4. **No need for IC-variation** at this layer. T1/T2/T3 are all robust under the wavepacket IC at sensible σ and k₀; the amendment's "single-mode IC for strict pass/fail" point (S680 §2.3) is satisfied here automatically because the mode-k phase check *is* a single-mode test even within a Gaussian envelope.

## §4 — What this session does not output

- **No verdict on Ingredient D**. The 2D/3D sweep with the four ingredients sequenced against the S17-22 baseline is the binary-outcome step, on Thor/Legion per the amendment §3.2.
- **No retag of Ingredient D's `[PARALLEL-PATHS]` status.** It stays `[PARALLEL-PATHS]` until the fleet sweep returns.
- **No comparison across A/B/C/D.** This session touched only D. The amendment's §1 cross-map (which orders B-A1 > B-A2 > D > A > C by leverage) is unaffected.
- **No retag of the S679 Priority 2 inventory.** D is not in the S679 (a)/(b)/(c) resolution map; the amendment §1 noted this is a different *kind* of move (drop the real-saturating-Intent axiom). That standing remains.
- **No cumulative tally.** Per the S679 discipline.

## §5 — Files

- `Research/Session681_Phase1_Ingredient_D_1D_Feasibility.md` (this document)
- `simulations/session681_phase1_ingredient_D_1D_feasibility.py` (split-step Fourier, three tests, damping characterization)

## §6 — So what (under the frame-doc disciplines)

A research-repo autonomous session can take the per-ingredient 1D pre-flight template that S680 established and apply it cleanly to the next-highest-leverage ingredient before fleet handoff. This session did that for Ingredient D and found: (i) split-step Fourier is the necessary integrator and should be locked into the spec for D, not left open; (ii) unitarity is preserved to machine precision at γ=0 for both free Schrödinger and GPE-class nonlinearity; (iii) γ>0 produces smooth damping (rate ∝ γ at small γ, saturates at large γ as R drops with |ψ|²) with phase dynamics persisting throughout — telegrapher-class qualitative shape rather than pure exponential decay. None of this resolves Phase-1's binary-outcome question for Ingredient D; it removes more implementation friction from the fleet path.

The amendment's §1.implications ordering stands: B-A1 first (sweep on Thor/Legion when ready), B-A2 second, **D** third — and D is now de-risked for split-step implementation and diagnostic interpretation, in the same shape S680 de-risked B for explicit leapfrog and energy-formula correctness. Remaining 1D pre-flights from the amendment §3.1: Ingredient A (focusing nonlinearity) and Ingredient C (external confinement, S18 entity criterion verification). Both still templatable from S680/S681 if a future autonomous-session firing picks them up.
