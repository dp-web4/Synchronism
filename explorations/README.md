# Synchronism Explorations

Empirical tests of the framework's core claims. Not interpretive reproducibility tests (cf. SAGE/explorations/) — these are **falsifiers of physics-content claims**.

## What an exploration is

An **exploration** is a low-cost, falsifiable experiment whose purpose is to determine whether one of Synchronism's core claims survives empirical scrutiny.

The Research/ tree contains thousands of sessions of theoretical derivation, framework refinement, and pattern catalog work. That work IS valuable, but it has a known structural failure mode (per external review): *AI-generated theoretical physics tends toward elegant isomorphism rather than empirical novelty*. The explorations/ directory is the operational counterweight — every load-bearing claim in the framework should either point to (or generate) an exploration that could falsify it.

A claim with no associated exploration is, by default, a **framing or theoretical position** — not a finding. (See README → "Findings vs Framings.")

## Format

Each exploration is a markdown file with:

1. **Hypothesis** — the specific claim being tested
2. **Procedure** — exactly what gets run (rule space, sweep parameters, scale, environment)
3. **Measurement** — what gets recorded and how it's scored
4. **Falsifier** — what would refute the hypothesis. If you can't state it, it's not yet a science question.
5. **Status** — drafted / running / partial-results / closed

The falsifier is the most important field. An exploration without a falsifier is just an opinion in markdown.

## Difference from `simulations/`

`simulations/` holds prior pattern-fitting and validation code (SPARC analysis, environment exploration, coupling-coherence runs). That work is data-fitting against existing physics datasets.

`explorations/` holds **forward-looking falsifiability experiments** — designed to either confirm or refute load-bearing claims of the framework that have not yet been tested.

## Fleet usage

Explorations are designed to be fleet-distributable. Parameter sweeps parallelize naturally; each fleet machine takes a region of parameter space and contributes to a shared result registry.

Per the 2026-05-15 fleet directive: **when ARC test runs aren't active, fleet machines may run exploration tasks**. The exploration outputs aggregate to `shared-context/synchronism/exploration-results/` (private — raw run data, results published in this public tree once each stage closes).

## Index

| Exploration | Status | Tests |
|-------------|--------|-------|
| [2026-05-15-cellular-automaton-discrete-grid-physics.md](2026-05-15-cellular-automaton-discrete-grid-physics.md) | **Stage 1 first result (2026-06-22):** monotonic-saturation substrate FAILS; focusing-nonlinearity family PASSES (>10%) — see [result](2026-06-22-phase1-stage1-localized-oscillation-result.md). Stages 2-5 drafted. | Whether local rules on a discrete grid produce stable resonant patterns → interaction → mass-like → field-like → quantum-like behavior |
| [2026-06-22-phase1-stage1-localized-oscillation-result.md](2026-06-22-phase1-stage1-localized-oscillation-result.md) | **Result** | Stage-1 falsifier executed: can the substrate self-confine a stable oscillating pattern? Monotonic saturation can't (refutes Foundation 3 at Stage 1); a focusing nonlinearity can, at the cost of Foundation 3. |

## Why this directory exists

External review (Kimi 2.6, 2026-05-15, full dialogue at `forum/kimi/kimi_2_6_review.md`) proposed the framework's most direct empirical test:

> *"If Synchronism's claim is that simple local rules on a discrete grid produce all observed physics, the proof is in the simulation. Can you write a cellular automaton rule that, run on a large enough grid, produces stable particle-like patterns with mass and charge, gravitational attraction, EM interaction, quantum interference? This is enormously difficult — arguably the hardest problem in theoretical physics. But it's also the most direct test of the framework's core claim. If the simulation works, the ontology is vindicated. If it doesn't, the framework needs revision."*

The cellular-automaton challenge is the first entry here because it's the cleanest available falsifier of the discrete-grid-with-resonant-patterns ontology. It is multi-stage; intermediate signal arrives at every stage; failure at any stage is a publishable result that constrains the framework.
