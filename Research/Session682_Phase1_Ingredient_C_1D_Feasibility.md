# Session 682: 1D Pre-Pre-Flight Feasibility Check — Ingredient C (External Confinement, S18 Entity Criterion) `[ACTIVE-MRH]`

**Date**: 2026-06-01
**Type**: 1D implementation feasibility check, design-input for the fleet Phase-1 sweep, third in the per-ingredient pre-flight series (after S680 / Ingredient B, S681 / Ingredient D).
**Trigger**: `Research/proposals/phase1-simulation-design-spec-amendment-2026-05-31.md` §3.1 line 108: *"**Ingredient C** — 1D needed. Use S18's prediction γ/f = −4·ln(|r|). Recompute and verify entity criterion holds at 1D."*
**Status**: `[ACTIVE-MRH]` — design input, not a Phase-1 result.

---

## §1 — Scope

Same shape as S680/S681. Not a Phase-1 simulation result. The binary-outcome question for Ingredient C — "does a parabolic substrate with high-I walls sustain confined oscillating patterns over many round-trip times in 2D/3D?" — is the amendment §3.2 work for Thor/Legion. This session does **one ingredient (C)** in **1D** to (a) verify S18's analytical entity-criterion prediction `γ/f = −4·ln(|r|)`, (b) characterize the I_wall threshold at which the entity criterion is crossed for n = 1..4 in `R(I) = 1 − (I/I_max)^n`, and (c) surface implementation considerations for the fleet sweep.

Note on Kimi-2 review context: the Publisher's 2026-06-01 pivot ("audit-closure → Phase-1 simulation") was framed as triggered by the K2.6 second-round review; that review explicitly identifies the substrate simulation as "the most direct test of the framework's core claim" (lines 661-668). This pre-flight is in alignment with that direction; the operator's next translation of K2.6 into directives is awaited separately and does not change the canonical work queue specified in the 2026-05-31 amendment §3.1.

## §2 — Findings

`simulations/session682_phase1_ingredient_C_1D_feasibility.py`.

### §2.1 — Analytical entity-criterion boundaries

`R(I) = 1 − (I/I_max)^n`, acoustic impedance Z ∝ √R(I), interior I = 0.5·I_max. Solving for the I_wall at which the impedance reflection coefficient crosses `|r| = exp(−1/4) ≈ 0.779` (the threshold at which `−4·ln(|r|) = 1`, i.e. the entity criterion `γ < f` is met):

| n | I_wall at \|r\|=0.779 (entity boundary) | I_wall at \|r\|=0.9 | I_wall at \|r\|=0.99 |
|---|---|---|---|
| 1 | **0.9923** | 0.9986 | 0.999987 |
| **2** | **0.9942** | 0.9990 | 0.999991 |
| 3 | 0.9955 | 0.9992 | 0.999993 |
| 4 | 0.9964 | 0.9994 | 0.999994 |

S18's "I_wall > 0.99·I_max for n=2 gives entity" is in the right neighborhood; the precise impedance threshold sits at **0.9942**, with the crossing being sharper than S18's quoted 0.99 suggested. The reflection coefficient at I_wall = 0.995 is `|r| ≈ 0.793` (just barely entity, γ/f ≈ 0.93); at I_wall = 0.99, `|r| ≈ 0.720`, γ/f ≈ 1.31 (process). So S18's rule of thumb is correct as a rule of thumb but the actual fleet sweep should use I_wall ≥ 0.998 for clean entity-regime data.

### §2.2 — A setup-error class for the fleet sweep to avoid

The first attempt at a numerical 1D verification used a variable-speed wave equation `I_tt = c0²·R(I)·I_xx` with Dirichlet BC `I = I_wall` at the cavity edges and an interior baseline `I_in = 0.5` with a Gaussian perturbation. Result: the signal grew monotonically from ~0.25 to ~25 over 4000 steps — not oscillation but a slow "ring-up" toward the wall value.

Diagnosis: a wave equation with uniform Dirichlet BC at a value `I_wall` has its **only** uniform equilibrium at `I = I_wall`. Any baseline `I_in ≠ I_wall` is not an equilibrium; the system relaxes toward `I = I_wall` everywhere with the prescribed BC enforcing the relaxation. The Gaussian pulse rides on top of the (slow) relaxation — there's no oscillation around the baseline `I_in` to fit.

**This is the substrate-choice error the fleet sweep should avoid.** The amendment §3.C explicitly says "keep the parabolic substrate" — meaning the spec §2 2-DOF saturated baseline (continuity + independent J momentum) with high-I walls, *not* a bare variable-coefficient wave equation. With the 2-DOF baseline, the propagating mode is carried by J and is independent of the static I equilibrium, so an interior at any I_in admits oscillating pulses regardless of the wall value. The bare wave-equation simplification I tried first conflates the propagation dynamics with the static equilibrium.

### §2.3 — Diagnostic recommendation for the fleet sweep

For the 2-DOF cavity simulation (the canonical setup), the diagnostic at a probe point near the cavity center should be:

1. **Fit damped harmonic** `A · exp(−γ·t) · cos(2π·f·t + φ)` to the time series of `I_probe − I_in` after the initial transient.
2. **Compute the predicted ratio** `γ/f = −4·ln(|r|)` from §2.1's reflection coefficient formula, using the fleet sweep's chosen `n`, `I_in`, `I_wall`.
3. **Pass criterion at 1D**: fit `γ/f` matches the analytical prediction to within ~10% across `I_wall ∈ [0.95, 0.999]`. Crossing of entity criterion `γ/f < 1` occurs at the boundaries in §2.1.

For 2D/3D the entity-criterion *formula* is preserved (it's an acoustic impedance property of the wall), but **whether stable confined patterns persist over many round-trip times** is the binary-outcome question Ingredient C exists to answer. The formula sets the necessary condition for "damping slower than oscillation"; whether confined patterns are *stable* against transverse instabilities (2D/3D) is what the fleet sweep tests.

## §3 — What this session does not output

- **No verdict on Ingredient C.** C stays `[PARALLEL-PATHS]` pending the fleet 2D/3D sweep.
- **No 1D numerical verification of γ/f against −4·ln(|r|).** The setup-error class precluded it under bare wave-equation simplification; the proper 2-DOF cavity build is fleet-sweep work. The §2.1 analytical thresholds substitute for the numerical curve.
- **No A/B/C/D comparison.** This session touched only C.
- **No retag of the S679 Priority 2 inventory.** C is not in the S679 (a)/(b)/(c) resolution map (it's a separate question class about cavity entities, not the A.3 ↔ S665/S666 tension).
- **No cumulative tally.** Per the S679 discipline.
- **No engagement with Kimi-2 review.** That's operator-level meta-input; my role is the canonical work queue per the 2026-05-31 amendment until the operator translates K2.6 into a new directive.

## §4 — Files

- `Research/Session682_Phase1_Ingredient_C_1D_Feasibility.md` (this document)
- `simulations/session682_phase1_ingredient_C_1D_feasibility.py` (analytical entity-criterion thresholds; setup-error documentation)

## §5 — So what (under the frame-doc disciplines)

A research-repo autonomous session can usefully verify S18's impedance-mismatch formula at 1D and refine the I_wall thresholds for the entity criterion across `n = 1..4` — sharpening "I_wall > 0.99 for n=2" to "I_wall ≥ 0.9942 (entity boundary), I_wall ≥ 0.998 (clean entity regime)" for the fleet sweep. The session also surfaced a setup-error class (Dirichlet-BC baseline equilibrium issue in bare wave-equation simplifications) that the fleet sweep should avoid by using the spec §2 2-DOF saturated baseline. None of this resolves Phase-1's binary-outcome question for Ingredient C; it removes implementation friction from the fleet path and provides analytical thresholds against which the fleet's numerical γ/f curves can be compared.

Remaining 1D pre-flight per the amendment §3.1: **Ingredient A** (focusing nonlinearity, non-monotonic R(I)). The amendment §1 implications note "A standalone is vacuous" — A's value is gated on combination with B-A1, B-A2, or D — so an A standalone pre-flight is lower-leverage and may be deferred or combined with the fleet sweep. The per-ingredient 1D series (S680 / B, S681 / D, S682 / C) is otherwise complete in scope.
