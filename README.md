# Synchronism: Unified Coherence Framework

A research program asking whether **one coherence parameter (γ = 2/√N_corr)** can describe phenomena across very different scales — quantum, classical, chemical, cosmological — and what falsifies it. Open-ended; calibrated; honest about where the framework breaks down. **[STATUS.md](STATUS.md)** is the calibration — read it before judging the claims below.

The framing is ambitious. The work tests the framing. Not all tests have passed.

## Five-minute audit

If you want a fast read on whether this is real, in order:

1. [**STATUS.md**](STATUS.md) — full honest assessment of what's tested, what's pending, what failed.
2. [**Research/discoveries/**](Research/discoveries/) — discoveries with explicit statistical evidence and falsification status. Each doc names what would refute the claim.
3. [**Honest Limitations** (further down this README)](#honest-limitations) — the things that *didn't* work. 53% error on chemistry melting points. 2× off on critical exponents. The framework's own boundary, named in the open.
4. [**Research Tracks**](#research-tracks-march-2026) — current work with session counts and verdicts. "Phase 2 complete: framework is organizational lens, not predictive theory" is itself a finding.
5. [**Coupling-Coherence Experiment**](Research/Coupling_Coherence_Experiment.md) — 900 runs, Hill function beat tanh by ΔAIC=4. Concrete, reproducible, falsifiable.

The work that follows is a research program testing a unifying hypothesis. The hypothesis being ambitious doesn't make the testing crankery; the testing being explicit about its failures is what distinguishes it from crankery.

---

## Why This Exists

Synchronism asks: *Can one coherence parameter explain phenomena from Planck scale to cosmic scale?*

**2,200+ research sessions** later, we have validated predictions across 1703 phenomenon types — and honest failures that define the framework's boundaries.

---

## What We've Discovered

### Major Validated Findings

| Discovery | Evidence | Status |
|-----------|----------|--------|
| **[γ ~ 1 Universal Boundary](Research/discoveries/gamma-universal-boundary.md)** | 1703 phenomenon types, 89% validated | Validated |
| **[NP2 RAR Scatter](Research/discoveries/np2-rar-scatter-validation.md)** | p = 5×10⁻⁶, all tests pass | Strongly Supported |
| **[Gnosis Consciousness Threshold](Research/discoveries/gnosis-consciousness-threshold.md)** | 8-way convergence at C ≈ 0.50 | Theory Complete |
| **[η Superconductivity Formalism](Research/discoveries/eta-superconductivity-formalism.md)** | Cuprate predictions match | Validated |
| **[ONE EQUATION Unification](Research/discoveries/one-equation-unification.md)** | Spans 80 orders of magnitude | Theoretical |

→ **[Full Discoveries](Research/discoveries/)**

### The Framework

```
γ = 2/√N_corr

where:
  γ = coherence parameter
  N_corr = correlated particles moving as unit
  γ ~ 1 = quantum-classical boundary
```

This single equation governs superconductivity, chemistry, biology, consciousness, and astrophysics.

### Honest Limitations

- Chemistry melting points: 53% error
- Critical exponents: 2× off from observations
- Environment vs structure ambiguity in NP2

→ **[STATUS.md](STATUS.md)** for complete honest assessment

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

## License

CC0 - Public Domain

---

*Last updated: February 5, 2026 | 2,200+ sessions | 1703 phenomenon types validated*
