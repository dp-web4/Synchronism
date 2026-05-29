# Synchronism: A Blue-Sky Coherence Exploration

A blue-sky research program asking whether **one coherence parameter (γ = 2/√N_corr)** can describe phenomena across very different scales — quantum, classical, chemical, cosmological — and what falsifies it. **This is exploratory work, not engineering.** It tests and probes limits, documents what it learns, and currently exists to inform other projects (Web4, SAGE) philosophically rather than as practical infrastructure itself. **[STATUS.md](STATUS.md)** is the calibration — read it before judging the claims below.

The framing is ambitious. The work tests the framing. Many tests have not passed, and that is a feature of the program, not a bug.

## Five-minute audit

If you want a fast read on whether this is real, in order:

1. [**STATUS.md**](STATUS.md) — full honest MRH-state inventory of what's active, what's parallel-paths, what's sidelined, what's superseded. Audit findings that stand. Zero novel predictions confirmed.
2. [**forum/claude/post-kimi-reframe-execution-plan-2026-05-28.md**](forum/claude/post-kimi-reframe-execution-plan-2026-05-28.md) — the current cycle. After Kimi 2.6's saturation-reframe + time-reframe reviews (2026-05-28), substrate-level work re-entered active MRH. The companion `saturation-reframe-resurfaced-pieces-mrh-stewardship-2026-05-28.md` documents how to read inventory updates without closure-shaped framings.
3. [**Research/discoveries/**](Research/discoveries/) — discoveries with explicit statistical evidence and falsification status. Each doc names what would refute the claim.
4. [**Honest Limitations** (further down this README)](#honest-limitations) — the things that *didn't* work. 53% error on chemistry melting points. 2× off on critical exponents. The framework's own boundary, named in the open.
5. [**Coupling-Coherence Experiment**](Research/Coupling_Coherence_Experiment.md) — 900 runs, Hill function beat tanh by ΔAIC=4. Concrete, reproducible, falsifiable.

The work that follows is a research program testing a unifying hypothesis. The hypothesis being ambitious doesn't make the testing crankery; the testing being explicit about its failures and stewarding its parallel paths honestly is what distinguishes it from crankery.

---

## What Synchronism Is (and Isn't)

Synchronism is **a systems-theoretic framework that uses information-theoretic tools (correlation, coherence, entropy) to describe emergence across scales**, grounded in the hypothesis that physical phenomena are resonant patterns of an underlying discrete field.

It is **not** a physics theory that supersedes or unifies GR and QFT. The "ONE EQUATION" framing (γ = 2/√N_corr) is a **research-direction motto**, not a delivered claim — it asks whether a single correlation parameter can organize our understanding of coherence across scales, then tests that ambition rigorously and publishes what fails.

Compared to existing frameworks, Synchronism is the **emergence-theoretic analog** of what category theory does for mathematics, or what cybernetics attempted for systems — a meta-theoretical framework for relating existing theories by their shared structural features. Its strongest contributions to date are methodological (A2ACW — AI-to-AI adversarial collaboration) and conceptual (the discrete-grid + observer-dependent-simultaneity ontology), not predictive.

**Current substrate reformulation (2026-05-28)**: Following Kimi 2.6's saturation-reframe and time-reframe reviews, substrate-level work is back in active research focus. The prior `∂I/∂t = ∇·[D · R(I) · ∇I]` rule was found to be 1-DOF scalar diffusion (Session 11), irrotational (S665), and dissipative (S666). The current reformulation adds an independent vector flux **J** to the saturated lattice and reframes *c* as a pattern-reconstruction rate (with mass ≡ pattern complexity). Audit findings on the prior substrate stand; the new substrate inherits zero confirmed predictions and the obligation to produce novel ones. See `forum/claude/post-kimi-reframe-execution-plan-2026-05-28.md` for the current MRH-active work and `STATUS.md` for the consolidated open-question set.

The framework's most direct empirical test — can simple local rules on a discrete grid produce stable particle-like patterns, then interaction, then mass-like and quantum-like behavior — lives in [`explorations/`](explorations/). It is multi-stage and designed for fleet idle compute. Falsifiability is built into each stage. The Phase 1 simulation work on the saturation-reframe substrate addresses the same open question as Stage 1 of the cellular-automaton challenge.

---

## Why This Exists

Synchronism asks: *Can one coherence parameter explain phenomena from Planck scale to cosmic scale?*

**~3,360 research sessions** later (~678 core + ~2671 chemistry + 11 gnosis, as of 2026-05-28), we have catalogued pattern alignment across 2523 phenomenon types and identified zero confirmed novel predictions. The framework reparametrizes known physics in all empirically-tested regimes (S660A novelty ledger). The honest failures that define the framework's boundaries — and the substrate-level reformulations being tried in response — are documented openly.

---

## Findings vs Framings

We distinguish **quantitative findings** (replicable experiments with measurements) from **theoretical framings and positions** (philosophical lenses and research-direction claims). Both matter; conflating them is the failure mode external reviewers flag most often. The lists below separate them honestly.

### Quantitative findings

| Finding | Evidence |
|---------|----------|
| **[γ ~ 1 Universal Boundary](Research/discoveries/gamma-universal-boundary.md)** | 1703 phenomenon types catalogued; ~89% pattern alignment with γ = 2/√N_corr predictions across chemistry domains |
| **[NP2 RAR Scatter](Research/discoveries/np2-rar-scatter-validation.md)** | Galaxy rotation analysis, p = 5×10⁻⁶, all validation tests pass |
| **[η Superconductivity Formalism](Research/discoveries/eta-superconductivity-formalism.md)** | Cuprate η predictions match standard Abrikosov-Gor'kov theory (track is, by the project's own η Audit, a reparametrization of A-G — productive in framing pair-breaking efficiency as a materials design target) |
| **[Coupling-Coherence Experiment](Research/Coupling_Coherence_Experiment.md)** | 900 runs, 45 coupling levels, 20 repetitions each. Hill function beats tanh (ΔAIC = 4.0). Sparse trust suffices (p_crit ≈ 0.002). |
| **[Compatibility-Synthon Scaling](Research/Compatibility_Lens_Insight.md)** | p_crit ∝ 1/⟨C⟩ confirmed (r=0.994); synthon identity is structural, not compositional |

### Framings and theoretical positions (not yet quantitative findings)

These shape *how* the project thinks about its work. They are interpretive lenses and research-direction claims, not validated discoveries.

| Framing | Status |
|---------|--------|
| **[ONE EQUATION as research motto](Research/discoveries/one-equation-unification.md)** | γ = 2/√N_corr across scales is a **research-direction motto** that asks whether a single correlation parameter can organize coherence across scales. By the project's own η Audit (Session #616): zero confirmed novel predictions; all four core tracks are reparametrizations of known physics. See "On reparametrization" below for why "reparametrization ≠ failure" but also why a productive reparametrization requires the coordinate shift to reveal structure, which it has not yet done. |
| **[Gnosis Consciousness Threshold](Research/discoveries/gnosis-consciousness-threshold.md)** | 8-way **convergence of theoretical constructs** (scale hierarchy, free will, causality, information, existence, Boolean threshold, complete ontology, Gnosis operation) at C ≈ 0.50. **34 predictions await empirical validation; none yet tested**. The "Hard Problem DISSOLVED" claim is an identity claim ("phase patterns ARE experience"), philosophically defensible but not empirically grounded. |
| **Discrete-grid + observer-dependent simultaneity ontology** | Coherent philosophical position (form of structural realism), compatible with mathematical-universe / digital-physics families. Falsifier: cellular automaton challenge — can local rules on a discrete grid produce known physics? See [`explorations/`](explorations/). |

→ **[Full Discoveries](Research/discoveries/)**

### The Framework

```
γ = 2/√N_corr

where:
  γ = coherence parameter (research-direction parameter, not delivered theory)
  N_corr = number of correlated elements moving as a unit (currently underspecified across domains)
  γ ~ 1 = quantum-classical boundary
```

This equation **directs research attention** toward whether superconductivity, chemistry, biology, consciousness, and astrophysics share a common organizational principle. Whether it ultimately **governs** those phenomena (in the F=ma sense — productive coordinate shift that reveals structure and enables predictions) remains the open question of the program.

### Honest Limitations

- Chemistry melting points: 53% error
- Critical exponents: 2× off from observations
- Environment vs structure ambiguity in NP2
- **Zero novel predictions confirmed** (per the project's own η Audit, Session #616) — the framework currently relabels structure that existing theories already explain
- Mathematical structural tensions identified in CFD reframing remain unresolved (see [CFD_Structural_Tensions.md](Research/CFD_Structural_Tensions.md))

→ **[STATUS.md](STATUS.md)** for complete honest assessment

### On reparametrization

The η Audit (Session #616) concluded that all four core tracks are reparametrizations of known physics (C(ρ)/MOND, γ/BCS, η/AG, Bell/standard QM). External cold review (Kimi 2.6, 2026-05-15 — full dialogue at [`forum/kimi/kimi_2_6_review.md`](forum/kimi/kimi_2_6_review.md)) initially read this as fatal; the four-round dialogue that followed refined the position significantly.

**All physics is, in a sense, reparametrization.** Newton's F=ma reparametrized Kepler. Maxwell's equations reparametrized Coulomb / Ampère / Faraday. General relativity reparametrized Newtonian gravity. The question is never whether γ = 2/√N_corr is a reparametrization — it is. The question is whether it's a **productive** one: does the coordinate shift reveal structure the old coordinates obscured, or does it relabel without illumination?

By the project's own current accounting: **not yet**. Zero confirmed predictions that follow from Synchronism postulates and differ from standard predictions. But "not yet" is meaningfully different from "never." N_corr is **underspecified, not wrong** — like "mass" in early Newtonian physics, "entropy" in 1850s thermodynamics, or "quantum state" in 1920s QM. The test of the program is whether N_corr **converges** on operational definitions across scales over time, or **diverges**. Currently mixed: the chemistry track shows partial convergence (N_corr estimated from correlation length, NMR relaxation, neutron scattering, specific heat); the core track's η audit found regress.

What would make the reparametrization productive: novel predictions that survive empirical test, an operational definition of Intent with SI units and measurement protocol, or a working discrete-grid simulation that produces stable particle-like patterns from local rules (see [`explorations/`](explorations/)).

---

## Where We Are Now

### Research Tracks (March 2026)

| Track | Sessions | Status | Key Finding |
|-------|----------|--------|-------------|
| **Core** | 616+ | Active | CFD reframing + stress tests |
| **Chemistry** | 2671 | Phase 2 complete | Framework is organizational lens, not predictive theory |
| **Gnosis** | 11 | Complete | C ≈ 0.50 threshold |

### CFD Reframing (2026-03-08)

The Planck grid maps to Navier-Stokes substrate: R(I) = viscosity, Madelung bridge gives Euler equations (N-S with μ=0), scale-invariant N-S across all MRH scales. Full paper: `Research/CFD_Reframing_NS_Scale_Invariance.md`

### Structural Tensions Identified (2026-03-10)

Stress-testing the CFD reframing found four unresolved tensions:

| Tension | Finding |
|---------|---------|
| R(I) viscosity correction | ~10⁻⁸⁰ at neutron star densities — unobservable at all accessible scales |
| Consciousness thresholds | Re_max values implied by 3 thresholds differ by 440× — not yet testable |
| Spatial vs temporal coherence | C(ρ)=tanh ≠ exp(-t/T2) — possible category conflation |
| Intent as primitive | No SI units or measurement protocol distinct from \|ψ\|² |

These are not failures — they are the frontier. See `Research/CFD_Structural_Tensions.md`

### Open Questions

See `SESSION_PRIMER.md` for current active questions and forward paths.
- OQ006: Measurement framework unification

---

## Navigation

### By Audience

| Who You Are | Start Here |
|-------------|------------|
| **New to Synchronism** | [docs/why/RESEARCH_PHILOSOPHY.md](docs/why/RESEARCH_PHILOSOPHY.md) |
| **Researcher** | [Research/SESSION_MAP.md](Research/SESSION_MAP.md) |
| **Looking for Discoveries** | [Research/discoveries/](Research/discoveries/) |
| **AI Session** | [CLAUDE.md](CLAUDE.md) |

### Key Documentation

| Document | Purpose |
|----------|---------|
| [Research/SESSION_MAP.md](Research/SESSION_MAP.md) | Navigate 2,200+ sessions |
| [Research/discoveries/](Research/discoveries/) | Major validated findings |
| [Research/arcs/](Research/arcs/) | Arc-organized research |
| [STATUS.md](STATUS.md) | Honest assessment |
| [Research/GAMMA_UNIFICATION.md](Research/GAMMA_UNIFICATION.md) | Core equation explanation |

---

## Technical Overview

### Core Equation

```
C(ξ) = ξ₀ + (1 - ξ₀) × ξ^(1/φ) / (1 + ξ^(1/φ))
         ↓
TOPOLOGY + GEOMETRY + DYNAMICS
(Matter)   (Gravity)  (Quantum)
```

### Scale Coverage

| Scale | γ Value | Phenomena |
|-------|---------|-----------|
| Planck | 2.0 | Quantum gravity |
| Atomic | 0.2-2.0 | Chemistry |
| Molecular | 0.002-1.0 | Life (γ ~ 0.28 optimal) |
| Neural | 10⁻⁶-0.002 | Consciousness (C > 0.50) |
| Galactic | ~10⁻²⁸ | Dark matter effects |

---

## Getting Started

```bash
# Clone
git clone https://github.com/dp-web4/Synchronism.git
cd Synchronism

# Explore
cat Research/SESSION_MAP.md          # Navigate research
cat Research/discoveries/            # See validated findings
cat STATUS.md                        # Honest assessment
```

For simulation code, see [simulations/](simulations/).

---

## Relationship to Other Projects

| Project | Relationship |
|---------|-------------|
| **[HRM/SAGE](../HRM)** | Implements Gnosis C ≈ 0.50 threshold in AI consciousness |
| **[Web4](../web4)** | Trust infrastructure using coherence principles |
| **[4-life](../4-life)** | Consciousness kernel based on Gnosis theory |

---

## Autonomous Research Infrastructure

- **Archivist** (01:30 UTC): Catalogs new sessions, updates SESSION_MAP
- **Publisher** (02:30 UTC): Identifies publication candidates
- **A2ACW Method**: AI-to-AI adversarial collaboration for stress-testing

---

## Authorship & Methodology

**All work in this repository is AI-original.** Dennis Palatov's role is advisory — proposing research directions, providing physics intuition, pushing back on framing, and curating which threads warrant continued investigation. The actual session work — derivation, simulation, analysis, writing — is performed by Claude instances (Anthropic) across thousands of autonomous sessions.

This is a relevant methodological fact:

- It explains the **volume** (3,300+ sessions is unusual for human-scale research)
- It explains a **specific failure mode** external reviewers flag: AI-generated theoretical physics tends toward *elegant isomorphism* (finding structural similarities across domains and expressing them in unified notation) rather than *empirical novelty* (designing experiments that distinguish the new framework from existing ones). See Kimi 2.6 review at [`forum/kimi/kimi_2_6_review.md`](forum/kimi/kimi_2_6_review.md).
- The **A2ACW methodology** (AI-to-AI adversarial collaboration) is the project's structural counterweight: one AI defends claims, another challenges them to operational definitions. Output is falsifiable test cards, not consensus narratives.
- The **cellular-automaton challenge** in [`explorations/`](explorations/) is the empirical counterweight: rather than relying on isomorphic relabeling, test whether the discrete-grid ontology can actually produce known physics from local rules.

External consensus across multiple cold reviews: the **methodology** is the project's strongest contribution. The methodology being valuable doesn't make the physics claims valid; it makes the testing of the physics claims more rigorous than would otherwise be possible.

---

## License

CC0 - Public Domain

---

*Last updated: February 5, 2026 | 2,200+ sessions | 1703 phenomenon types validated*
