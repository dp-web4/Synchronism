# Session #384: Structural Analysis Arc Synthesis

**Date**: 2026-02-06
**Status**: Synthesis document (no tests)

## Arc Overview: Sessions #381-383

This arc addressed the fundamental ambiguity left by the Gas Fraction Control Arc (#376-378): is the NP2 signal (morphology → RAR scatter) driven by **environment** (Synchronism) or **intrinsic structure**?

## Arc Narrative

### Act 1 - Session #381: The Challenge
**Question**: Environment or structure?
**Method**: Composite isolation proxy from SPARC observables
**Result**: Structure favored 5-1. Within-type isolation shows zero signal.
**Meaning**: Internal SPARC proxies cannot separate the effects.

### Act 2 - Session #382: The Mechanism
**Question**: How does type → scatter work mechanically?
**Method**: Roughness decomposition, mediation analysis
**Result**: Roughness mediates 88% of the type → scatter link.
**Meaning**: Late types have messier kinematics → more scatter. This is structural.

### Act 3 - Session #383: The Residual
**Question**: What about the remaining 12%?
**Method**: Systematic offset analysis across acceleration regimes
**Result**: The offset is MOND-dominated, mass-dependent, and γ-consistent.
**Meaning**: A small but genuine signal may reflect modified gravitational coupling.

## Three-Component Scatter Model

```
╔═══════════════════════════════════════════════════════════════╗
║  RAR SCATTER = ROUGHNESS + SYSTEMATIC OFFSET + UNEXPLAINED   ║
╠═══════════════════════════════════════════════════════════════╣
║                                                               ║
║  Component 1: ROUGHNESS (R² = 0.51)                          ║
║  ├── Origin: Structural (kinematic irregularities)            ║
║  ├── Mechanism: Non-circular motions, warps, bars             ║
║  ├── Type dependence: 1.58x (late/early)                     ║
║  └── Relevance to γ: NONE (structural artifact)              ║
║                                                               ║
║  Component 2: SYSTEMATIC OFFSET (R² = 0.11)                  ║
║  ├── Origin: CONTESTED (M/L vs γ)                            ║
║  ├── Direction: Late types BELOW RAR (rotate slower)          ║
║  ├── Regime: 58% stronger in MOND regime                     ║
║  ├── Mass scaling: r(Vflat, offset | T) = +0.30              ║
║  └── Relevance to γ: POSSIBLE (maximum γ window)             ║
║                                                               ║
║  Component 3: UNEXPLAINED (R² = 0.38)                        ║
║  ├── Origin: Galaxy-to-galaxy peculiarities                   ║
║  └── Relevance to γ: Unknown                                 ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```

## Updated NP2 Assessment

| Aspect | Status | Evidence |
|---|---|---|
| Scatter difference exists | CONFIRMED (p = 5×10⁻⁶) | Sessions #376-378 |
| Difference survives confounds | CONFIRMED | Sessions #377-378 |
| Driven by environment | NOT SUPPORTED (internal data) | Session #381 |
| Mediated by roughness | CONFIRMED (88%) | Session #382 |
| Residual offset γ-consistent | POSSIBLE | Session #383 |
| Residual offset M/L-consistent | ALSO POSSIBLE | Session #383 |

**Final NP2 status: PARTIALLY SUPPORTED**
- The scatter difference is real and robust
- But it's primarily structural (roughness), not environmental
- An 11% systematic offset component is consistent with γ but not uniquely explained by it

## What We Learned About Synchronism

### Narrowing the γ Signal

Before this arc: NP2 = type → scatter (5.3% R², p = 5×10⁻⁶)
After this arc: γ window = systematic offset (11% R² of total scatter, or ~2% of total variance)

The γ signal, if it exists, is much smaller than originally appeared. It lives entirely in systematic RAR offsets, not in kinematic roughness.

### Revised γ Mechanism

The data tells us what γ DOESN'T do:
- γ does NOT increase kinematic roughness
- γ does NOT create radius-dependent scatter patterns
- γ does NOT correlate with internal isolation proxies

What γ MIGHT do:
- Shift galaxies systematically off the standard RAR
- Affect the MOND regime more than the Newtonian regime
- Scale with galaxy mass (via N_corr)
- Produce galaxy-to-galaxy offsets that correlate with morphology

### The Critical Experiment

The definitive test remains: **explicit environment data**. Specifically:
1. Cross-match SPARC with Chae et al. (2020) external field estimates
2. Test: does external field correlate with systematic offset, controlling type and mass?
3. If yes → γ (environment effect confirmed)
4. If no → M/L or baryonic model (stellar population effect)

## Arc Statistics

- Sessions: 3 (381-383)
- Tests verified: 24/24
- Galaxies analyzed: 171
- Statistical methods: Partial correlation, mediation analysis (Sobel test), path analysis, Monte Carlo, multi-variate regression, Mann-Whitney, acceleration-regime stratification
- Novel discoveries: 3 (roughness mediation, MOND-dominated offset, three-component model)
- Honest limitations: 6 (per session)
- Key negative result: Internal proxies cannot separate environment from structure

## Recommended Next Arcs

| Priority | Arc | Description | Rationale |
|---|---|---|---|
| **HIGH** | g† First Principles | Derive a₀ from γ = 2/√N_corr | Core theoretical prediction |
| HIGH | External Field Test | Cross-match Chae et al. (2020) | Definitive environment test |
| MEDIUM | Theoretical Recalibration | Reformulate γ → RAR with offset mechanism | Theory needs updating |
| MEDIUM | Literature Comparison | Compare three-component model with published work | Context needed |
| LOW | Extended Sample | THINGS/LITTLE THINGS rotation curves | Larger sample |

---

**★ STRUCTURAL ANALYSIS ARC COMPLETE ★**
*Sessions: 3/3 (381-383)*
*Tests verified: 24/24*
*Grand Total: 511/511*

**Arc discovery: The NP2 signal decomposes into roughness (51%, structural, NOT γ) and systematic offset (11%, MOND-dominated, γ-consistent but M/L-competitive). The γ theory survives but is constrained to a narrow window of systematic RAR offsets. External environment data (Chae et al. 2020) is the critical next experiment.**
