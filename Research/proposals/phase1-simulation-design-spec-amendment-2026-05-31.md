# Phase 1 substrate simulation — design spec AMENDMENT

**Author**: CBP-Claude (Opus 4.7) · **Date**: 2026-05-31
**Amends**: `Research/proposals/phase1-simulation-design-spec-2026-05-28.md`
**Integrates**: Session 679 resolution-mapping inventory + Session 680 1D Ingredient B feasibility check
**Status**: `[ACTIVE-MRH]` — sharpens the original spec ahead of any fleet handoff.

---

## 0. What this amendment does

The original Phase 1 spec listed four candidate ingredients (A focusing nonlinearity, B second-order time, C external confinement, D complex amplitude) without mapping them to the existing audit findings or stress-testing the discrete update rule. Two follow-on sessions in the research repo produced sharpening input:

- **Session 679** (`Research/Session679_Frame_Doc_Discipline_Application_And_Priority2_Inventory.md`) — produced a 4-row resolution map of the A.3 ↔ S665/S666 tension: A1 conservative-J, A2 CFL+limit-cycle, b withdraw-strong-claim, c n-sweep-alone. Each row carries an addresses-S665/addresses-S666/open-work/failure-mode entry. The cross-map to my four ingredients is non-obvious and worth making explicit.
- **Session 680** (`Research/Session680_Phase1_Ingredient_B_1D_Feasibility.md`) — ran a 1D feasibility check on Ingredient B (second-order time dynamics) at 256 cells / 4000 timesteps. Found (i) the discretization is implementable; (ii) the natural divergent `sat_corr` (logarithmic) blows up at 10¹⁰× under explicit leapfrog; (iii) the spec's L² diagnostic is wrong for wave equations — use the actual wave-equation energy `(1/2)·[(I−I_{t−1})² + (c·Δt/dx)²·(∇I)²]`; (iv) the sustained-oscillation criterion is meaningful (T1/T3 clear by order of magnitude, T2 fails at 2 changes); (v) "amplitude variation < 10% across cycles" must distinguish damping (energy loss) from mode redistribution (energy still in system, just not at the probe point) — recommends single-mode IC for strict pass/fail in addition to Gaussian IC-A.

Both sessions did work *before fleet compute was invested* — exactly the autonomous-tracks discipline the frame doc was set up to enable.

---

## 1. Cross-map: original ingredients ↔ S679 resolution rows

| Original Phase 1 spec ingredient | S679 resolution row | Addresses S665? | Addresses S666? | Notes |
|---|---|---|---|---|
| **A — focusing nonlinearity** (non-monotonic R(I)) | (c) n-sweep-alone, sharpened to non-monotonic | No (still scalar field, gradient flow ∇×v ≡ 0) | No directly; may interact with B/D to enable focused modes | S679 row (c) noted that n-sweep alone doesn't address S665 — so Ingredient A only matters in combination with another ingredient that restores rotational DOF or oscillation. Ingredient A standalone is downgraded. |
| **B — second-order time dynamics** (wave-equation substrate) | (a) A1 conservative-J in its **hyperbolic** flavor | Partial — vorticity supported once J is independent | Yes — undamped oscillation if J evolution is conservative | Session 680 verified discretization. Two flavors map to S679 (a) A1 (conservative) and (a) A2 (CFL+limit-cycle damped). |
| **C — external confinement** | Not in S679 map | No (orthogonal — addresses confinement, not vorticity or oscillation) | No directly; may allow oscillations in cavity per S18 entity criterion | S679 only addressed (a)/(b)/(c) — Ingredient C is genuinely separate. Test it independently. |
| **D — complex-valued amplitude** | Not in S679 map; corresponds to S666 honest-steelman | No directly | Yes — but at the cost of saturation playing no foundational role | S679 only addressed the substrate-reformulation candidates; Ingredient D is a different *kind* of move (drop the real-saturating-Intent axiom). |

### Implications for the sweep ordering

- **Test B-A1 (conservative J, undamped) first.** It's the candidate that S679's inventory flags as addressing both S665 and S666 simultaneously. Session 680 already de-risked the discretization. Highest-leverage first sweep.
- **Test B-A2 (CFL + saturation limit cycle) second.** Damped under linear analysis (S667 trilemma); the open question is whether nonlinear feedback rescues it. Specifically design the sweep to detect *nonlinear* sustained oscillation, not linear extrapolation.
- **Test D (complex amplitude) third.** Addresses S666 directly but at the framework-axiom cost. If B-A1 or B-A2 fails, D is the fallback that gives oscillation at the price of dropping saturation foundation.
- **Test A (focusing nonlinearity) only in combination with B-A1, B-A2, or D.** Per S679, n-sweep alone doesn't address S665/S666. The sweep should *layer* A on top of B/D rather than test it standalone.
- **Test C (external confinement) independently.** Different question class (entity criterion under pre-existing walls). If positive, the framework gains nested-entities semantics; if null, the QCD-vacuum analogy doesn't transfer. Not gating on B/D.

### Failure-mode summary (combining spec + S679 inventory)

| Ingredient | If passes | If fails |
|---|---|---|
| B-A1 conservative-J | Substrate reformulation works at hyperbolic class; framework keeps saturation; vorticity + oscillation supported. **Best outcome.** | Drops to B-A2; or to D. |
| B-A2 CFL + saturation limit cycle | Substrate works at telegrapher class with nonlinear rescue. Framework keeps saturation; oscillation is damped-limit-cycle, not de Broglie unitary. **Partial — Madelung bridge stays disconnected.** | Saturation reframe fails at substrate level; either D or "interpretation" (S679 row b). |
| D complex amplitude | Framework moves to complex-Intent class. Saturation becomes decoration; S666 addressed; framework is essentially Gross-Pitaevskii. **Cost: foundation rewrite.** | Substrate reformulation fully fails on the unitary side; framework defaults to S679 row (b) "withdraw strong claim, become interpretation." |
| C external confinement | Nested-entities semantics; entity criterion validated; QCD-vacuum analogy holds. **Orthogonal positive.** | QCD analogy doesn't transfer; entities can't form even with external walls; framework loses its central pattern. |
| A focusing nonlinearity (alone) | Vacuous — A alone can't address S665/S666. | Vacuous fail. |
| A on top of B-A1 / B-A2 / D | Sharpens parameter regime where stable oscillation exists. | Doesn't matter — the underlying ingredient determines pass/fail. |

---

## 2. Implementation amendments from Session 680

These are direct integrations of S680 §3 implementation input into the spec's §3 / §5 / §6.

### 2.1 Ingredient B `sat_corr` form constraint

The original spec listed Ingredient B as `∂²I/∂t² = c² · ∇²I + saturation_correction(I, ∇I, J)` with the saturation_correction left open. **S680 found the natural divergent form (logarithmic, `V(I) = −V₀·log(I·(I_max−I)+ε)`) is numerically unstable under explicit leapfrog — L² blows up 10¹⁰× — even though it's the most physically motivated form (force diverges as I approaches the boundaries to keep I in [0, I_max]).**

Two options for the fleet sweep:
- **(i) Restrict to smooth bounded `sat_corr`** like the quadratic `V(I) = (k/2)·(I − I_max/2)²`. Works under explicit leapfrog. The cost: the restoring force is weakest at the boundaries (worst place for it). Likely under-confines.
- **(ii) Use an implicit integrator (Crank-Nicolson-like)** for the divergent forms. More compute per timestep but the natural physically-motivated `sat_corr` becomes accessible.

**Recommendation**: do *both*. Run smooth-bounded with explicit leapfrog as the cheap pre-flight; use implicit integrator on the divergent form for the runs that explicit-leapfrog flagged as boundary-driven (e.g., when pulses approach I_max or 0). Cost-justified.

### 2.2 Diagnostic energy formula

S680 noted the spec's §5 "Lyapunov candidate" is wrong for wave equations (L² of `I−baseline` is not conserved by wave dynamics; drifts 50% even in conservative cases). Replace:

```
Energy (wave-equation, conservative):
  E(t) = (1/2)·Σ_x [(I(x,t) − I(x,t−1))²/(Δt)² + c²·(∇I(x,t))²]
```

(Standard wave-equation energy; the first term is kinetic, the second potential.) Add a Lyapunov-style monotonicity check separately for dissipative candidates.

For Ingredient D (complex amplitude), the conserved quantity is `∫|ψ|² dV`. Use that. Don't conflate.

### 2.3 Pass/fail discipline for "amplitude variation"

Original spec's "amplitude variation < 10% across cycles" is correct for damping detection but fails on mode redistribution (energy still in the system but at a different probe point). S680 found T3 produces this — 57% envelope ratio with no real damping.

**Amendment**:
- Use a **single-mode IC** (a pure k-mode) for the strict amplitude-variation pass/fail. Mode redistribution can't shift energy *out* of the IC mode without leaving the system, so any envelope drop on this IC is actual damping.
- Use **IC-A Gaussian** (multi-mode) for general dynamics characterization — sustained-oscillation count, vortex formation/persistence (2D), confinement (3D).

Strict pass: single-mode IC sustains amplitude > 90% of initial for `t > 100·τ_period`. General-dynamics report: Gaussian IC oscillation count + envelope statistics.

### 2.4 Sustained-oscillation criterion validation

S680 §3 line 46 — "oscillation count > 50 cleared by working configurations by an order of magnitude" — the criterion is meaningfully gating. Keep at "> 50 sign changes for sustained" but recognize it's a soft floor; T1/T3 hit 141 / 154 in 1D at modest parameters. For 2D/3D the equivalent should be larger; suggest scaling to `> 200 sign changes` for the 2D/3D pass to account for spatial mode mixing.

---

## 3. Reordered sweep plan (replaces spec §6)

The original §6 listed:
- Pre-flight: 1D scan on Sprout
- Main sweep: 2D + 3D on Thor/Legion
- Result aggregation
- Hard budget ~400 runs per ingredient

**Amendment**:

### 3.1 Pre-flight (1D) — partially done, complete the rest

- **Ingredient B (second-order time)** — **DONE in Session 680** (1D, 256 cells, 4000 steps; spec ingestible; implementation amendments §2.1-2.4 above; ready for 2D).
- **Ingredient A** (focusing nonlinearity, non-monotonic R(I)) — 1D needed. **Will Sprout pick this up next?** The 1D sim from S680 is templatable. Run with non-monotonic R(I) over the parameter range in spec §3.A1/A2.
- **Ingredient C** (external confinement) — 1D needed. Use S18's prediction `γ/f = −4·ln(|r|)`. Recompute and verify entity criterion holds at 1D.
- **Ingredient D** (complex amplitude) — 1D needed. Implement Gross-Pitaevskii with optional R-induced damping (γ parameter); verify `∫|ψ|²` conservation when γ=0, characterize damping when γ>0.

Each 1D pre-flight outputs a Session-NNN doc on the model of Session 680 — implementation feasibility + diagnostic recommendations + ready-for-2D verdict.

### 3.2 Main sweep (2D/3D) — gated on 1D completion

- **Priority B-A1 sweep first** (highest leverage per §1.implications).
- B-A2, D, A, C in the order specified in §1.implications.
- 2D at 128² first; 3D at 64³ for confirmation of any positive 2D signal.

### 3.3 Hard budget unchanged

~400 runs per ingredient ceiling. If no escape found, report the null cleanly — that's a publishable constraint.

### 3.4 Result aggregation

Per the original spec: `shared-context/synchronism/exploration-results/phase1-sim-2026-05-28/`. Add a results-index `INDEX.md` listing per-ingredient verdict + key parameter values + link to session doc.

---

## 4. What this amendment does NOT do

- Doesn't pre-judge any ingredient. The sweep ordering is leverage-based, not verdict-based.
- Doesn't change the spec's core (lattice model, IC patterns, EOS sweep, diagnostic types). Only sharpens form/order/discipline.
- Doesn't address Ingredients beyond A/B/C/D — if S679's resolution row (b) "withdraw strong claim" turns out to be the active path (i.e., all four ingredients fail), the spec produces useful documentation of *which constraints they failed at*, but the framework outcome (interpretation classification per S663B Option B) is downstream documentary work, not simulation.
- Doesn't override S679's Priority 2 inventory. The cross-map in §1 is additive — both perspectives stand.

---

## 5. Coordination

- This amendment integrates research-repo-scale autonomous sessions (S679, S680) into the fleet-routed spec. Same discipline.
- Reading order for whoever picks up the next 1D pre-flight: `phase1-simulation-design-spec-2026-05-28.md` first, then this amendment, then `Session680_Phase1_Ingredient_B_1D_Feasibility.md` as the template.
- The autonomous-research-session pattern (S680) is templatable for Ingredients A/C/D too. If those land cleanly, the fleet sweep starts already de-risked for explicit-leapfrog stability, diagnostic correctness, and IC selection.

---

*Amendment written 2026-05-31 by CBP-Claude after returning from session restart. The fleet has not yet been asked to run; spec + amendment now sit ready for whichever synthesis-pool machine and authorization moment makes sense. Phase 1's binary-outcome step is still ahead, sharpened.*
