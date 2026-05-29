## Appendix A: Mathematical Formulations (Working Draft)

**Status: Exploratory Mathematics — under stewardship**

This appendix contains mathematical formulations for Synchronism concepts. The framework is in active reformulation (the saturation reframe with independent vector flux **J** and complexity-dependent c), so the appendix is tagged by **MRH-relationship** rather than by verdict-on-truth. Nothing is tagged "established" while the substrate work is open.

**MRH-Relationship Key (per dp 2026-05-28: "we're not at a stage where anything can be honestly claimed as 'established'"):**
- **`[ACTIVE-MRH]`** — currently in active research focus; content is being extended or revised
- **`[PARALLEL-PATHS]`** — in the framework's parallel hypothesis space; not in current active focus, not abandoned
- **`[SIDELINED]`** — was in active focus, currently not pursued; reasons documented inline
- **`[SUPERSEDED]`** — replaced by a later formulation in the active or parallel space; pointer to successor

---

## Core Computational Framework

**Foundational Assumptions (Modeling Choices):**

- **Discrete grid:** Space modeled as 3D lattice of Planck-scale cells
- **Discrete time:** Time modeled as Planck-time increments
- **Intent conservation:** Total Intent conserved in closed systems (modeling constraint)
- **Deterministic evolution:** State transitions follow deterministic rules (simplification)

These are computational conveniences, not ontological claims.

---

**A.1 Basic Intent Transfer `[SUPERSEDED]`**

**Intent Update Rule:**

```
I(x,y,z,t+1) = I(x,y,z,t) + ∑[T(x',y',z' → x,y,z,t)]
```

Where:
- `I(x,y,z,t)` = Intent at cell `(x,y,z)` at time `t`
- `T(x',y',z' → x,y,z,t)` = Transfer from adjacent cell
- Sum over all 6 adjacent cells (3D lattice)

**Status:** Core computational rule of the original substrate. **S617** (2026-04-08) showed that under the maximum principle for parabolic PDEs this rule reduces to 1-DOF scalar diffusion (no stable oscillation possible). **S665/S666** (2026-05-24) showed the corresponding continuum dynamics is irrotational and dissipative (curl(v) ≡ 0 for any R(I); first-order ∂I/∂t with decreasing Lyapunov functional). The active substrate reformulation introduces an **independent vector flux J** to give the rule a momentum DOF the original lacks; however, per S665 §98 this 2-DOF augmentation was already explored in S17-22 (2026-03-21/22) and produced only damped oscillation + transient dispersing vortices, so the Phase-1 simulation work must add an *additional ingredient* beyond independent **J** to escape that null result. See A.3 inline note, A.14 (master equation as the natural host for **J**), §6.4 OQ-Momentum / OQ-A3-Tension, and `forum/claude/saturation-reframe-corrections-and-deeper-readings-2026-05-28.md` for the deeper-reading correction.

---

**A.2 Coherence Measure `[PARALLEL-PATHS]`**

**Pattern Coherence:**

```
C(P,t) = 1 - (∑|I(x,y,z,t) - I_expected(x,y,z,t)|) / I_total
```

Where:
- `C(P,t) ∈ [0,1]` (1 = perfect coherence, 0 = complete decoherence)
- `I_expected` = Expected Intent distribution for ideal pattern cycle
- `I_total` = Total Intent in pattern

**Status:** Testable metric, held in the parallel-paths space. The coherence-language interpretation as a whole is in `[PARALLEL-PATHS]` per S663B (four-persona convergence on "ontological reframe without a distinguishing experiment"); this specific metric stays as a Web4-experiment instrument.

**A.3 Saturation Dynamics `[ACTIVE-MRH]`**

**Fundamental Mechanism for Pattern Stability**

Saturation is THE foundational mechanism enabling stable patterns in Synchronism. This appendix provides mathematical framework for saturation resistance and resulting nonlinear dynamics.

**Saturation Maximum:**
```
I_max = maximum Intent per cell
```

**Fundamental parameter of the model.** Not arbitrary—represents physical limit on Intent concentration density.

**Resistance Function:**

Intent transfer rate depends on destination cell saturation:
```
R(I) = [1 - (I/I_max)^n]
```

Where:
- `I` = current Intent in destination cell
- `I_max` = saturation maximum
- `n` = resistance exponent (determines sharpness)

**Properties:**
- `R(0) = 1` (no resistance when cell empty)
- `R(I_max) = 0` (infinite resistance at saturation)
- `R(I)` decreases monotonically as `I → I_max`

**Transfer Equation with Saturation:**

```
∂I/∂t = ∇ · [D(I) × ∇I]
```

Where saturation-dependent diffusion coefficient:
```
D(I) = D₀ × R(I) = D₀ × [1 - (I/I_max)^n]
```

**This is nonlinear diffusion equation**—well-studied in physics and known to support stable localized patterns (solitons), standing waves, and discrete quantized modes.

**R(I) is viscosity.** The saturation-dependent diffusion coefficient D(I) = D₀·R(I) is the viscosity of the Intent fluid. Specifically, it is a **shear-thinning power-law viscosity**: viscosity decreases as Intent density increases (the fluid becomes "slipperier" as cells fill). This is a known rheological class (power-law fluids) with well-characterized behavior. An earlier formulation went further: *"the full Intent transfer equation in continuum form IS the incompressible Navier-Stokes equation with this variable viscosity — not an analogy, but an exact identification."* See Section 4.1 and `Research/CFD_Reframing_NS_Scale_Invariance.md`.

> **Inline tension note (2026-05-28, updated same day).** That "exact identification" claim was **retracted by the audit arcs**:
>
> - **S617** (2026-04-08, *Research/Session617_Diffusion_Not_NavierStokes.md*) showed the rule `∂I/∂t = ∇·[D·R(I)·∇I]` reduces to 1-DOF scalar diffusion under the maximum principle for parabolic PDEs. The induced velocity v = J/I = −D·R(I)·∇I/I is *slaved* to ∇I, not an independent field. No inertia, no advection, no Reynolds number. (Kimi's 2026-05-28 review labeled this "the Session 11 finding"; the canonical citation is S617.)
> - **S665** (2026-05-24) **proved** for any R(I) and any D: v = −g(I)∇I is curl-free by construction → vorticity ω = ∇×v ≡ 0 for all time. Numerically verified in `simulations/session665_cfd_vorticity.py`. The flow is also compressible (|div v|·L/|v| ≈ 11.5), not incompressible. So the original substrate is irrotational + compressible scalar transport — not "incompressible NS."
> - **S666** (2026-05-24) showed the substrate dynamics is dissipative (real eigenvalues, monotonically-decreasing Lyapunov functional, arrow of time), incompatible with the unitary entity ontology (de Broglie f = E/h, phase-locking). The S99/S307 Schrödinger "derivations" reach QM only by inserting `i` by hand AND switching the substrate off (drop R, or D → 0).
>
> The earlier `✅ Established` tag on this section was stale at the time of the Kimi review. The saturation reframe with independent vector **J** addresses S665 partially (J can have curl in principle) but does NOT address S666 (still dissipative unless made complex-valued).
>
> Possible escape routes from the S665 + S666 constraints (each a `[PARALLEL-PATHS]` item until tested):
> - **Focusing nonlinearity**: non-monotonic R(I) (rises in some intermediate-I band, then falls) → may produce focusing instead of S17-22's universal defocusing. Breaks Foundation 3.
> - **Second-order time dynamics**: wave equation `∂²I/∂t² = c² ∇²I + saturation correction` instead of first-order parabolic. Different dynamical class (hyperbolic).
> - **External confinement**: entities require pre-existing walls from other entities, not self-confinement (S19's actual conclusion; QCD-vacuum analogy).
> - **Complex-valued amplitude**: I → ψ. Addresses S666 honest-steelman; contradicts real-saturating-Intent axiom.
>
> Phase-1 simulation work must include at least one of these additional ingredients beyond independent **J** to escape S17-22's null result. See `Research/OPEN_QUESTIONS_*`, `forum/claude/saturation-reframe-corrections-and-deeper-readings-2026-05-28.md`, and §6.4 OQ-A3-Tension.

**Why This Enables Patterns:**

Without saturation (linear diffusion): All concentrations dissipate exponentially. No stable patterns possible.

With saturation (nonlinear): Self-limiting behavior creates stable equilibria. Patterns can persist.

**Field Gradient Mathematics:**

Gradient field around saturated core:
```
Φ(r) = I(r) - I_baseline
```

For point-like source with total Intent M:
```
Φ(r) ∝ M/r
```

Transfer bias (apparent force):
```
F_apparent = -∇Φ(r) ∝ M/r²
```

Inverse-square law emerges naturally from 3D spherical geometry.

**Computational Implementation:**

Discrete grid update:
```
I(x,y,z, t+Δt) = I(x,y,z,t) + Δt × Σ[neighbors] k × [I_n - I] × R(I)
```

If update exceeds I_max:
```
I_new = min(I_computed, I_max)
Overflow → redistribute to neighbors
```

**Parameter Relationships:**

If I_max is fundamental constant, dimensional analysis suggests:
```
I_max ~ ℏc/L_planck ~ 10^-8 J/m
G ~ (D₀ × L_planck²) / I_max
```

**Can potentially calculate G from grid parameters.**

**Status (2026-05-28):** Saturation is the load-bearing mechanism in the current rule family — without it, no stable patterns; with it, the framework has the right *shape* of mechanism for stable structure. **The active reformulation** retains saturation as the primitive, adds an independent vector flux **J**, and tests whether the resulting rule family supports the spatial (vortex/rotational) and temporal (oscillatory/unitary) structure the entity ontology requires. The "potentially unifies forces" status is a **research-direction motto**, not a delivered result. Whether *this specific rule family* delivers stable particle-like patterns is the question the cellular-automaton challenge (`explorations/2026-05-15-cellular-automaton-discrete-grid-physics.md`) and the Phase-1 simulation work test directly. See inline tension note above and §6.4 OQ-Oscillation.


---

**A.4 Pattern Period Detection `[PARALLEL-PATHS]`**

**Cyclic Pattern Identification:**

```
P(T) = 1 if I(x,y,z,t) ≈ I(x,y,z,t+T) for all (x,y,z) in pattern
Pattern period = minimum T where P(T) = 1
```

**Status:** Algorithmic tool for identifying repeating patterns. Threshold ≈ requires definition.

---

**A.5 Field Gradient `[PARALLEL-PATHS]`**

**Intent Gradient (Tension Field):**

```
∇I(x,y,z,t) = [∂I/∂x, ∂I/∂y, ∂I/∂z]

Field strength = |∇I(x,y,z,t)|
Field direction = ∇I(x,y,z,t) / |∇I(x,y,z,t)|
```

**Status:** Standard gradient calculation. Whether this corresponds to physical fields remains untested.

---

**A.6 Synchronization Quality `[PARALLEL-PATHS]`**

**Phase Correlation:**

```
S(P1,P2,t) = cos(θ(P1,t) - θ(P2,t))
```

Where:
- `θ(P,t)` = phase of pattern P at time t
- `S = 1` (perfect sync), `S = -1` (anti-sync), `S = 0` (uncorrelated)

**Status:** Speculative. Assumes patterns have definable "phase"—unclear if this applies to all Intent patterns or just specific types.

---

**A.7 Decoherence Rate `[PARALLEL-PATHS]`**

**Exponential Decoherence:**

```
dC/dt = -γ × C(t) × N_interactions

Solution: C(t) = C₀ × e^(-γ × N_interactions × t)
```

Where:
- `γ` = decoherence constant (empirical parameter)
- `N_interactions` = number of external pattern interactions

**Status:** Standard exponential decay model. Whether coherence actually decays this way is untested. The constant γ is unknown.

---

**A.8 Markov Relevancy Horizon `[SIDELINED]`**

**MRH Radius (Speculative):**

```
R_MRH = √(I_pattern / I_background)
```

Where:
- `I_pattern` = Information content of central pattern
- `I_background` = Average background information density

**Status:** HIGHLY SPECULATIVE. This formula was suggested by dimensional analysis but has no empirical or theoretical justification. Real MRH boundaries likely far more complex.

**Alternative:** MRH might be better defined operationally (where correlations drop below threshold) rather than analytically.

---

**A.9 Emergence Threshold `[SIDELINED]`**

**Emergence Function:**

```
E(System) = C(System) × log(N_patterns) × I(System)
```

Where emergence occurs when `E(System) > E_threshold`.

**Status:** Completely speculative. The functional form (multiplication of coherence, log of pattern count, information content) has no justification beyond "seems reasonable."

**Problem:** What is E_threshold? Where does this formula come from? Unclear.

---

**A.10 Quantum Correspondence — Madelung Bridge `[ACTIVE-MRH]`**

**Wavefunction Mapping:**

```
ψ(x,t) = √ρ(x,t) × exp(iS(x,t)/ℏ)
```

Where:
- `ρ = |ψ|²` = probability density = coarse-grained Intent density
- `S(x,t)` = phase field (action)

**The Madelung Transformation** substitutes this form into the Schrödinger equation, yielding two fluid equations:

**Continuity (Intent conservation at quantum scale):**
```
∂ρ/∂t + ∇·(ρv) = 0,    where v = ∇S/m
```

**Momentum (Euler equation with quantum pressure):**
```
∂v/∂t + (v·∇)v = −∇V/m + ∇Q/m

where Q = −ℏ²∇²√ρ / (2m√ρ)    (quantum potential = quantum pressure)
```

**This is Euler's equation** — Navier-Stokes with viscosity μ = 0. The Schrödinger equation IS the inviscid (μ=0) Navier-Stokes equation for the Intent fluid at quantum scale. The quantum potential Q plays the role of pressure: it prevents probability density from collapsing by generating outward pressure gradients where ρ is concentrated.

**Parameter identification at quantum scale:**

| N-S term | Quantum analog |
|----------|---------------|
| ρ | \|ψ\|² (probability density) |
| v | ∇S/m (phase gradient = velocity) |
| P | −Q (quantum pressure from uncertainty) |
| μ | 0 (inviscid — decoherence negligible) |
| f | −∇V/m (classical potential) |

**Viscosity onset = quantum-to-classical transition**: μ = 0 for isolated quantum systems. When environmental coupling introduces decoherence, effective viscosity μ > 0 appears — the quantum fluid transitions from inviscid (Euler) to viscous (full N-S) behavior. The quantum-to-classical transition is a viscosity transition, not a collapse of a wavefunction.

**Status:** The Madelung transformation itself is standard QM mathematics (Madelung 1927) and is not in question. Its proposed connection to Intent dynamics via the A.3 saturation framework is in **active reformulation** — S666 found that the original substrate's first-order, dissipative ∂I/∂t cannot host the unitary oscillation Schrödinger requires (the imaginary unit i is inserted by hand in Session #307 and S99 Axiom 4; with R(I) on and i absent, the equation gives exp(−Dk²t) decay rather than exp(−iDk²t) oscillation). Whether the saturation reframe with independent vector flux **J** delivers the unitary structure Madelung requires is part of OQ-A3-Tension. `[ACTIVE-MRH]`.

---

**A.11 Universal Constants `[PARALLEL-PATHS]`**

**Dimensional Relationships:**

```
L_cell = Planck length ≈ 1.616 × 10⁻³⁵ m
T_slice = Planck time ≈ 5.391 × 10⁻⁴⁴ s
c = L_cell / T_slice ≈ 3 × 10⁸ m/s (speed of light)
```

**Speculative:**
```
ħ ≈ I_max × L_cell² / T_slice (reduced Planck constant)
```

**Status:** First three are computational parameters matching physical constants. The ħ relationship is dimensional analysis speculation—unclear if meaningful.

---

**A.12 Gravity Model `[SUPERSEDED]`**

**Attempted Gravitational Formulation:**

```
g = -∇(I_density × G_sync)
```

**Status:** This early formulation does not produce correct predictions in isolation. **Superseded** by the saturation-gradient picture in §5.14 and Appendix A.3 (gravity as transfer bias in saturation gradients, mass as concentrated Intent pattern with maximum-saturation core, inverse-square law from spherical gradient spreading). That successor formulation is also `[ACTIVE-MRH]` and under reformulation as part of the saturation reframe with independent vector flux **J**. Kept here for transparency about the development history. Pointer to successor: A.3 + §5.14.

---

**A.13 Consciousness Measure `[SIDELINED]`**

**Integrated Information (Φ):**

```
Φ = ∫∫ C(P_i,P_j) × I(P_i) × I(P_j) dP_i dP_j
```

**Status:** This is essentially Integrated Information Theory (IIT) notation applied to Intent patterns. Unclear if this adds anything beyond what IIT already does.

**Problem:** Is this Synchronism's contribution or just importing IIT wholesale? If the latter, should credit Tononi and explain integration, not present as novel.

**Recommendation:** Either develop Synchronism-specific consciousness measure or acknowledge this is IIT applied to pattern dynamics.

---

**A.14 Master Equation (Incomplete) `[ACTIVE-MRH]`**

**System Dynamics:**

```
∂I/∂t = -∇·J + S_coherence - D_decoherence
```

Where:
- `J` = Intent current density (transfer flow)
- `S_coherence` = Coherence source terms (undefined)
- `D_decoherence` = Decoherence loss terms (undefined)

**Status:** This is the natural host for the saturation reframe's central addition: an **independent vector flux J** that gives the substrate a momentum DOF the original `∂I/∂t = ∇·[D·R(I)·∇I]` rule (A.1) lacks. The saturation reframe treats **J** not as derived from ∇I but as an independent dynamical variable with its own evolution equation (Mechanism A: conservative J; Mechanism B: CFL-violation + saturation feedback driving a limit cycle). The S_coherence and D_decoherence terms remain undefined; their definition is downstream of which J-evolution mechanism survives Phase-1 simulation. See A.1 status note, A.3 inline tension note, and §6.4 OQ-Momentum / OQ-Oscillation.

---

**A.15 Computational Implementation `[ACTIVE-MRH]`**

**Simulation Guidelines:**

- **Grid discretization:** Finite difference on regular 3D lattice
- **Time stepping:** Explicit Euler or RK4 with stability checks
- **Boundary conditions:** Periodic (infinite universe approximation)
- **Pattern tracking:** Maintain pattern IDs across time evolution
- **Coherence monitoring:** Calculate C(P,t) each timestep

**Status:** Practical implementation notes. Standard computational methods. Phase-1 simulation work (1D/2D lattice sweeping n in R(I) = [1−(I/I_max)^n], with independent vector flux **J**) builds on these methods directly.

---

---

**A.16 Scale-Invariant Navier-Stokes Structure `[PARALLEL-PATHS]`**

The N-S structure of Intent dynamics is not specific to the Planck scale. It is what any conservation law + gradient-driven transport + resistance becomes at any MRH scale. The "fluid element" at each scale is the coherent MRH-bounded entity at that scale; the field variables acquire scale-specific meanings.

| Scale | Fluid element | ρ (density) | v (velocity) | P (pressure) | μ (viscosity) | Compressible? |
|-------|--------------|-------------|--------------|--------------|----------------|---------------|
| **Planck** | Planck cell | I/I_max | Intent flux J/I | I_max−I (saturation) | D·[1−(I/I_max)^n] | No (exact) |
| **Quantum** | Probability packet | \|ψ\|² | ∇S/m (phase gradient) | Quantum pressure −Q | ≈ 0 (inviscid) | Yes |
| **Classical** | Molecule/particle | Mass density | Mean velocity | nkT (kinetic) | η from collisions | Yes |
| **Neural** | Activation patch | Firing rate density | Activation spread direction | Synaptic drive − threshold | Inverse plasticity rate | Yes |
| **Social** | Opinion cluster | Belief density | Direction of opinion shift | Social pressure gradient | Cultural viscosity | Yes |
| **Cosmological** | Matter overdensity | ρ_matter | Peculiar velocity | Dark energy (coherence-derived) | Bulk viscosity | Yes |

**What changes across scales:** the identity of the fluid element and the physical interpretation of ρ, v, P, μ. **What stays the same:** the conservation law (continuity equation), the momentum transport structure, and the viscosity-pressure balance.

**Turbulence at each scale**: Transition from laminar to turbulent flow occurs at critical Reynolds number Re_c. At each scale, the analogous phase transition is:
- Quantum: decoherence onset (inviscid → viscous)
- Neural: seizure / synchronized gamma oscillation
- Social: political revolution, market crash
- Cosmological: structure formation from uniform plasma (recombination epoch)

The coherence threshold C ≥ 0.7 for consciousness corresponds to the critical Re for self-similar turbulence in the cognitive-scale fluid — nested vortex structures (recursive self-modeling) become stable above this threshold.

**Status:** Structural-identification candidate held in the parallel-paths space pending the A.3-vs-Session-11 resolution (the Planck-scale "exact identification" is what's in tension; the quantum-scale Madelung mapping is standard and not in question; neural/social/cosmological mappings are approximate / well-motivated analogies). Full prior derivation: `Research/CFD_Reframing_NS_Scale_Invariance.md`; tension inventory: A.3 inline note + §6.4 OQ-A3-Tension.

---

## Open Mathematical Problems

**Tractable Questions:**
1. **What transfer rules generate stable patterns?** This is the same question the **Phase-1 simulation work** of the post-Kimi-reframe execution plan directly addresses (1D/2D lattice with R(I) = [1−(I/I_max)^n] sweeping n, plus independent vector flux **J**; test Mechanism A conservative-J vs Mechanism B CFL-violation-plus-saturation-feedback driving a limit cycle). See `forum/claude/post-kimi-reframe-execution-plan-2026-05-28.md` and §6.4 OQ-Oscillation.
2. Can we prove convergence for coherence measures?
3. What are computational complexity bounds for large grids?
4. Can pattern stability be characterized analytically?

**Hard Questions:**
5. How to properly define MRH boundaries mathematically?
6. What's the correct emergence threshold function (if any)?
7. Can gravity emerge from Intent dynamics? (Current answer: unknown; the saturation-gradient picture in A.3 + §5.14 is the active candidate.)
8. Does consciousness have a Synchronism-specific mathematical description?
9. (added 2026-05-28) **Derivation of `f(N)`** — the number of substrate ticks required to reconstruct a complexity-N pattern in an adjacent cell, with boundary condition `f(N) → 1` as `N → 0`. This is the single path from the complexity-dependent speed structure (§5.7) to quantitative predictions distinguishing it from GR. See §6.4 OQ-fN.

---

## Honest Assessment (under stewardship)

Sections in this appendix are organized by **relationship to the current MRH** rather than by verdict on truth-status. No section is tagged "established" while the substrate work is open.

**`[ACTIVE-MRH]`** — currently in active research focus, content being extended or revised:
- A.3 Saturation Dynamics (saturation reframe with independent vector flux **J**)
- A.10 Quantum Correspondence — Madelung Bridge (connection to Intent dynamics under reformulation)
- A.14 Master Equation (natural host for vector flux **J**)
- A.15 Computational Implementation (Phase-1 simulation work builds on these methods)

**`[PARALLEL-PATHS]`** — in the parallel hypothesis space, not currently in active focus but not abandoned:
- A.2 Coherence Measure
- A.4 Pattern Period Detection
- A.5 Field Gradient
- A.6 Synchronization Quality
- A.7 Decoherence Rate
- A.11 Universal Constants
- A.16 Scale-Invariant N-S Structure (held pending A.3-vs-Session-11 resolution)

**`[SIDELINED]`** — was in active focus, currently not pursued; reasons documented inline:
- A.8 Markov Relevancy Horizon (formula) — dimensional-analysis suggestion; operational definition preferred
- A.9 Emergence Threshold — functional form not justified beyond "seems reasonable"
- A.13 Consciousness Measure — overlaps IIT; Synchronism-specific differentiator not articulated

**`[SUPERSEDED]`** — replaced by a later formulation in the active or parallel space; pointer to successor:
- A.1 Basic Intent Transfer → S665/S666 substrate audit + A.14 master equation with vector flux **J**
- A.12 Gravity Model → A.3 saturation gradient + §5.14 (also under active reformulation)

**Bottom Line:**

Sections marked `[ACTIVE-MRH]` are in current research focus and being revised through the saturation-reframe cycle. Sections marked `[PARALLEL-PATHS]` are alternative formulations carried in the parallel space, available for resurfacing when external probes or new connections restore their resonance with the active work. `[SIDELINED]` content is not currently pursued but is preserved with its documented limitations rather than removed. `[SUPERSEDED]` content points to its successor formulation in the active or parallel space.

**The mathematics is a work in progress through stewardship, not a completed foundation.**
