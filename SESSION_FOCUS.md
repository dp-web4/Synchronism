# Synchronism Session Focus

*This file contains current research state, open questions, and session priorities. Updated by both the operator and autonomous sessions.*

*Last updated: 2026-03-22*

---

## ⚠️ PRIORITY: Conservation Bug Hypothesis

**READ**: `Research/CONSERVATION_BUG_HYPOTHESIS_2026-03-22.md`

The 810 failed configurations (Sessions #18-27) may all be artifacts of a broken conservation law. The transfer rule `ΔI = k·Σ(I_n - I)·R(I_n)` implicitly destroys momentum at saturation boundaries — when R(I) → 0, flow energy vanishes instead of redirecting. This violates the foundational axiom that intent is neither created nor destroyed.

**The fix**: Make boundaries elastic, not absorptive. When transfer is blocked by R → 0, the momentum must reflect/redirect rather than vanish. Simplest: track a velocity field alongside I, reverse velocity component when hitting saturation.

**Session 17 result**: Tested. R(I) soft walls produce **DAMPED oscillation** (73 sign changes, amplitude decays 0.3→0.001). Partial reflection occurs but smooth R(I) transition absorbs energy each bounce. Hard walls give sustained oscillation; R(I) walls give decaying oscillation. The conservation bug hypothesis is PARTIALLY SUPPORTED — adding momentum (2 DOF) does produce oscillation, but damped, not sustained.

**Connection found**: Damping rate γ vs oscillation frequency f maps to entity criterion: γ < f → Γ < m → entity. If γ/f derivable from wall geometry, entity criterion follows from substrate dynamics.

**Session 18 result**: Derived. γ/f = -4·ln(|r|) where r = (√R_in - √R_wall)/(√R_in + √R_wall). Entity criterion γ < f requires |r| > 0.779 → walls need I > 0.99·I_max (for n=2). Verified against session 17 (I_wall=0.95 gives γ/f=2.62 → process, consistent with observed damping). The entity criterion IS derivable from 2-DOF cavity impedance mismatch.

**Session 19 result**: Tested nonlinear self-consistent regime. Localized high-I pulse oscillates (697 sign changes) but DISPERSES from width 13 → 81. NOT self-confined. R(I) produces DEFOCUSING nonlinearity: wave speed = c·√R(I) DECREASES at high I, so pulse edges (low I) travel faster than center (high I) → dispersal. Soliton-like self-confinement requires FOCUSING nonlinearity (opposite of R(I)). Self-consistent entities cannot form from R(I) saturation alone.

**Implication**: Entities require EXTERNAL confinement (pre-existing walls from other entities or boundary conditions), not SELF-confinement. This is analogous to how atomic nuclei confine quarks — the confinement comes from the QCD vacuum, not from the quarks themselves. The entity criterion γ/f = -4·ln(|r|) remains valid for externally confined cavities.

**Session 20 (analytical)**: Coupled oscillations also fail. Between two high-I pulses, R decreases → wave speed decreases → BARRIER, not trap. Pulses repel rather than confine. ANY monotonically decreasing R(I) is defocusing. The framework needs a focusing mechanism not currently specified. Options: (1) modified R that increases with I at some scale, (2) tensor/vector Intent with vortex confinement, (3) multi-scale R depending on coarse-grained I at different MRH. None currently in the framework.

**STATUS**: The 2-DOF dynamics (I+v) produces damped oscillation in externally-confined cavities. The entity criterion is derivable as γ/f = -4·ln(|r|). But self-confinement fails due to defocusing R(I). The framework can describe entities that EXIST in cavities but cannot explain how cavities FORM. This is the current frontier.

---

## Current Research State

### CFD Reframing (2026-03-08)

The Planck grid IS the substrate of Navier-Stokes, not merely analogous to it.

- **R(I) = [1-(I/I_max)^n] is viscosity** — shear-thinning power-law fluid. Intent transfer equation in continuum = incompressible N-S exactly.
- **Madelung bridge**: ψ = √ρ·exp(iS/ℏ) into Schrödinger → Euler equations (N-S with μ=0). Quantum potential Q = pressure.
- **Scale-invariant N-S**: Same ρ/v/P/μ structure at every MRH scale.
- **Consciousness threshold = critical Reynolds number**: C ∝ 1/μ_eff.

Full paper: `Research/CFD_Reframing_NS_Scale_Invariance.md`

---

## Open Questions

### From CFD Reframing

1. ~~**Sequential vs parallel grid**~~ **RESOLVED: parallel.** Arrow of time comes from irreversibility of update rule F, not scan direction.

2. **Oscillation basis of existence** — recurrence over tick sequences defines entities. Period = de Broglie frequency = E/h. Connection to C(ρ): coherence C may be stability measure of recurring pattern. **Currently blocked on conservation bug — see priority above.**

3. **RG as formal MRH coarse-graining**: MRH abstraction implicitly does renormalization group analysis. Making this explicit connects to well-developed math.

4. **Viscosity profile μ(scale)**: Continuous function from Planck through quantum to classical to neural to cosmic. Zeros = coherent quantum regimes.

5. **MRH → emergent spacetime**: If MRH boundary separates DOF you track from those you don't, it plays the role of a causal horizon. Same structure as holographic screen.

### Computational Validation Results (Sessions #18-27)

**Five mechanisms tested, 0/810 oscillations** — but see conservation bug hypothesis above.

1. **Pure diffusion** (#18): 324 configs, 0 oscillations
2. **Reactive-diffusion** (#19): 300 configs, 0 oscillations
3. **Geometric confinement** (#20): 72 configs, walls form (45.8%), 0 oscillations
4. **Momentum-augmented** (#25): 112 configs, R(I) damping overwhelms momentum
5. **3D synchronous** (#27): 2 configs (32³, 64³), cavities form, same failure mode as 1D

**Critical finding**: Dimensionality and update model NOT the bottleneck. Confinement works. Missing piece is elastic boundaries (conservation enforcement).

### Older Open Questions

- **OQ006**: Measurement framework integration (#250 + #291). See `Research/OPEN_QUESTION_Measurement_Framework_Integration.md`
- **OQ005**: Hot superconductor — AUDITED: standard condensed matter in η notation. 0 unique predictions.
- **OQ007**: Fractal coherence bridge — CLOSED negative. C(ρ) is classification tool, not explanatory theory.

---

## Validation State

| Claim | Status |
|-------|--------|
| γ = 2/√N_corr unification | ✅ Validated |
| Coupling-coherence sigmoid (900 runs) | ✅ Validated |
| Madelung bridge | ✅ Standard QM math (1927) — Intent identification is the claim |
| Entity criterion (Γ < m) | ⚠️ NOW DERIVABLE — γ/f = -4·ln(|r|), entity when |r| > 0.779. Follows from 2-DOF cavity impedance mismatch. Walls need I > 0.99·I_max. Caveats: requires 2-DOF (not original 1-DOF rule), linear approximation, 1D. |
| CFD/N-S reframing | Speculative — mathematically consistent, not yet formalized |
| Consciousness thresholds as Re | ❌ Untestable as stated — Re_max values differ by 440× |
| Oscillation from substrate | ⚠️ PARTIALLY RESOLVED — 2-DOF (I+v) produces DAMPED oscillation in R(I) cavity (73 sign changes). Not sustained. Damping from smooth R(I) wall absorption. Entity criterion may follow from damping rate < oscillation frequency. |
| R(I) viscosity correction | ❌ Unobservable — correction ~10⁻⁸⁰ at neutron star densities |
| Dark matter = high viscosity | ❌ Sign error vs Bullet Cluster + internal contradiction |
| Lorentz invariance from parallel update | ❌ Logical gap — no discrete 3D lattice has SO(3) |
| N-S mapping: 1 DOF vs 2 DOF | ❌ Structural problem — v derived from I, not independent |

---

## Stress Test Results (Sessions #1-7 of stress test arc)

| Issue | Status |
|-------|--------|
| Strategic ambiguity | ⚠️ 0 commitments under pressure across 7 sessions |
| "Is" vs "models" protection | ⚠️ Named — divergence regime is Planck-scale, unobservable |
| Consciousness ontology fork | ⚠️ Emergence vs panpsychism — incompatible, both invoked |
| Incompressibility error | ❌ Math error — global conservation ≠ ∇·v = 0 |
| C(ρ) superseded by C(Re) | ⚠️ Resolved — C(ρ) was proxy, C(Re) is 4-parameter |
| Oscillation basis forced fork | ⚠️ Must choose: trivial (periodic orbits) or retrocausal (formation-time prediction) |

---

## Novel Predictions Awaiting Testing

| Prediction | What it says | What would test it | Status |
|-----------|-------------|-------------------|--------|
| **Entity criterion (Γ < m)** | Decay width > mass → process, not entity | QCD exotica with Γ/m > 1 | Untested — consistent with f0(500) |
| **Grid geometry → LIV** | Cubic Planck grid → Lorentz violation at ξ₂ ~ 1 | GRB polarimetry (AMEGO/CTA/GECAM-C) | Untested — requires cubic grid commitment |
| **Formation-time bound** | Constitutive recurrence: electron t < 8×10⁻²¹ s | Ultrafast spectroscopy ~10⁴× below current | Untested — requires retrocausal commitment |

---

## Whitepaper Status

| File | Recent Change |
|------|---------------|
| `universe_grid.md` | R(I) as shear-thinning viscosity; N-S parameter table |
| `crt_analogy.md` | Simultaneity as temporal MRH construction |
| `appendix_c_consciousness.md` | Thresholds as critical Reynolds numbers |
| `mathematical_framework.md` | A.10 Madelung bridge; A.16 scale-invariant N-S |
| `executive_summary.md` | CFD reframing entry |

---

## Explainer Site

https://synchronism-site.vercel.app — auto-deploys from `synchronism-site` repo

---

## Research Posture

### Reliable, not deterministic

LLM outputs navigate a probability landscape shaped by context. Results are reliable but never mechanical. Sessions are responsive navigation, not retrieval.

### Unconfirmed ≠ unconfirmable

Novel predictions exist. They are unconfirmed because this lab can't run the experiments, not because the predictions are wrong. Discovery requires someone with the means to look, to look.

### Interactive selection

We probe, observe, adjust, reinforce. We don't create or delete — we interactively select from what's latent. This applies to the research itself: the investigation is the value, not just the results.
