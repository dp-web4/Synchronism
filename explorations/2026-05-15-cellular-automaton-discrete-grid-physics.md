---
date: 2026-05-15
status: drafted (5 stages); Stage 1 simulator spec pending
proposed by: Kimi 2.6 (in `forum/kimi/kimi_2_6_review.md`)
fleet: idle-compute task (runs when ARC isn't running)
related:
  - Research/CFD_Reframing_NS_Scale_Invariance.md
  - Research/CFD_Structural_Tensions.md
  - Research/Coupling_Coherence_Experiment.md
  - Research/Compatibility_Lens_Insight.md
---

# Cellular Automaton — Discrete Grid Physics Challenge (CADGPC)

## The core hypothesis under test

**Synchronism's discrete-grid ontology** claims that physical reality is built from **simple local rules on a discrete grid**, with what we observe as "particles," "fields," and "interactions" being **stable resonant patterns** of an underlying Intent field. If this ontology is correct, then the structures of physics should be **derivable from a simulation** — not just *describable* by Synchronism's vocabulary.

Per external review (Kimi 2.6, 2026-05-15): *"This is enormously difficult — arguably the hardest problem in theoretical physics. But it's also the most direct test of the framework's core claim. If the simulation works, the ontology is vindicated. If it doesn't, the framework needs revision."*

This exploration tests that claim in stages. Each stage has its own hypothesis and falsifier; failure at any stage constrains the framework; success at any stage raises confidence and feeds the next.

## The saturation hypothesis (load-bearing)

A central premise of Synchronism (consistent with the Coupling-Coherence experiment's Hill-beats-tanh finding by ΔAIC=4) is that **emergent collective behavior shows saturation** — beyond some local critical mass / coupling / density, the response saturates rather than scaling linearly. Saturation is what creates **discrete emergent levels** (particles ARE patterns at saturated coherence; molecules ARE patterns at the next saturation level; etc.). Linear or unbounded response wouldn't produce discrete "things."

**This means the cellular automaton rule family being swept MUST include saturation in its update mechanism.** A rule with linear (or otherwise unbounded) response to neighbor influence is, by Synchronism's own logic, in the wrong rule family to produce stable particles.

But Synchronism has not committed to a *specific* saturation mechanism. Several are candidates, each with different consequences for emergent structure:

| Saturation mechanism | Functional form | Source / precedent | What it should favor |
|---------------------|-----------------|-------------------|---------------------|
| **Hill function** | f(x) = x^n / (K^n + x^n) | Cooperative binding kinetics; won the Coupling-Coherence comparison (ΔAIC=4 vs tanh) | Sharp phase transitions, cooperative emergence, possibly cleaner particle-like quantization |
| **Hyperbolic tangent** | f(x) = tanh(x/x₀) | Logarithmic saturation; widely used in NN activation | Smooth transitions, less cooperative, possibly fuzzier patterns |
| **Logistic** | f(x) = 1 / (1 + exp(-(x-x₀)/k)) | Sigmoidal; standard biological activation | Similar to Hill in shape, simpler analytically |
| **Hard threshold** | f(x) = 1 if x>θ else 0 | Conway's Life precedent; brittle discreteness | Sharpest cutoff; may produce Life-like discrete patterns |
| **Power-with-cutoff** | f(x) = min(x^α, C) | Fractal scaling with bounded response | Scale-free until saturation, then capped |
| **Cooperative-then-decay** | x^n / (K^n + x^n) × exp(-x/L) | Cooperative emergence with eventual quench | Spatially localized patterns by construction |

**Comparing saturation mechanisms is part of the experiment, not a pre-experiment choice.** Stage 1's parameter sweep is over BOTH the rule parameters AND the saturation mechanism. Different mechanisms may produce qualitatively different emergent behaviors; that itself is data. If only one family produces stable patterns, that constrains the framework's ontology.

## Why the simple-scalar-diffusion rule won't work

The CFD reframing's current Intent dynamics rule (∂I/∂t = ∇·[D·R(I)·∇I]) is, per `CFD_Structural_Tensions.md`, a 1-DOF scalar diffusion that cannot produce oscillating entities (maximum principle for parabolic PDEs). The cellular automaton experiment uses the **discrete-grid CFL-allowed-oscillation regime** that `CFD_Structural_Tensions.md` identifies as the framework's actual oscillation mechanism — meaning we test rules with at least one of:

- **Multi-component cell state** (e.g., scalar field + velocity field, or complex-valued field with phase + amplitude)
- **History-dependent updates** (e.g., 2-step rules: state[t+1] depends on state[t] AND state[t-1])
- **CFL-violating timescales** in the discrete update (so the continuum limit doesn't apply, and oscillation is genuine grid behavior not parabolic relaxation)

Any rule in the sweep that doesn't include one of these mechanisms is in the "guaranteed-not-to-oscillate" subspace and is included only as a control (negative result).

---

## Stage 1 — Stable Resonant Patterns (the simplest, most decisive test)

### Hypothesis

There exists a region of the (rule family × saturation mechanism × cell state space × neighborhood) parameter space such that random initial conditions on a sufficiently large discrete grid evolve to include **stable, localized, oscillating patterns** that persist >> T_lifetime and survive small perturbations.

### Procedure

**Grid**: 2D (initially — 3D is Stage 1.5 once 2D produces results). Side length L ∈ {64, 128, 256, 512}. Periodic boundary conditions.

**Cell state**: Multi-component as above. Initial sweep over:
- *Complex scalar* (ψ = re^iφ; amplitude + phase) — minimum mass-like representation
- *Real 2-vector* (analog of (ρ, J) for fluid density + flux) — preserves CFD reframing
- *Discrete N-state* (N ∈ {3, 5, 9}) — closer to Game-of-Life family, sharp quantization

**Update rule**: parameterized as `s[t+1] = f_saturation( linear_combination(neighbors, s[t], s[t-1]) )` where:
- The linear combination's coefficients are searched over a finite grid
- `f_saturation` is one of the six candidates above (Hill, tanh, logistic, hard threshold, power-with-cutoff, coop-decay)
- Neighborhood is Moore (8 cells in 2D) for first sweep; expand to von Neumann + Moore-extended in later sweeps

**Run length**: T = 10^5 update steps initially. If patterns survive, extend to 10^7.

**Initial conditions**: 100 random initializations per (rule × saturation × state-space) point. Patterns that emerge in ≥10 of 100 initializations from a given parameter point are considered "robust to initial conditions" at that point.

### Measurement (per parameter point, per initialization)

1. **Pattern detection**: After warm-up (10^3 steps), identify localized features by clustering algorithm (connected components above amplitude threshold)
2. **Persistence**: Track each pattern through subsequent timesteps. Pattern is "stable" if its centroid stays within ±N cells of mean position and amplitude variance < ε for T_observe steps.
3. **Oscillation**: Pattern is "resonant" if its internal state oscillates (Fourier peak above noise floor at some frequency)
4. **Perturbation robustness**: Inject small random noise into a stable pattern's vicinity. Measure recovery time and final state.
5. **Distinct pattern types**: How many qualitatively different stable patterns appear across the 100 random initializations? Are they discrete types or a continuum?

### Falsifier (Stage 1)

The hypothesis is **refuted** if, across the swept parameter space (6 saturation mechanisms × ~10^4 rule-coefficient points × 3 state-space variants × 4 grid sizes × 100 initializations each = O(10^7-10^8) total runs), **no parameter point produces stable resonant patterns** that survive perturbation in >10% of initializations.

The hypothesis is **supported** (but not yet vindicated) if at least one rule family + saturation mechanism reliably produces stable resonant patterns, with at least 3 distinct pattern types observed across initializations.

### What Stage 1 outcomes mean for the framework

| Outcome | Interpretation | Next step |
|---------|----------------|-----------|
| No saturation mechanism produces stable patterns | Discrete-grid-with-resonant-patterns ontology is in trouble FOR THIS RULE FAMILY. Either expand the rule family (3-step rules, longer-range neighborhoods, multi-scalar fields) or accept that the ontology needs revision. | Constraint on the framework. Publishable as the first concrete falsifier-result. |
| Only one saturation mechanism produces patterns | The framework's saturation hypothesis is real BUT specific (e.g. "must be Hill-like") rather than general | Refine Synchronism's saturation claim to be specific |
| Multiple mechanisms produce patterns, with different morphology | Different saturation classes are different "physics regimes" | Map which mechanism corresponds to which observed physics |
| Patterns appear but require very fine-tuned parameters | Patterns are not generic; ontology is fragile | Concern flag — real physics is not fine-tuned in this way |
| Patterns appear robustly across many parameters | Strong support for the discrete-grid ontology, move to Stage 2 | Most exciting outcome; Stage 2 unlocks |

### Compute budget

For a 256×256 grid running 10^5 update steps with simple per-cell updates: ~1-10 seconds per run on a modern CPU, fast on GPU. The full sweep (~10^7 runs) is large but embarrassingly parallel. Fleet partitioning:

- Each machine takes a slice of the (saturation × state-space) cross product
- Reports back per-parameter-point summary statistics (no need to ship full runs back)
- Aggregation registry in `shared-context/synchronism/exploration-results/cadgpc-stage1/`

Approximate fleet-aggregate runtime: weeks to months at low duty cycle (run only when ARC isn't), faster if dedicated bursts allowed.

---

## Stage 2 — Pattern-Pattern Interaction

Triggered only if Stage 1 succeeds.

### Hypothesis

Two stable patterns from Stage 1, placed in proximity on the grid, **interact non-trivially** — their trajectories influence each other in ways that resemble physical interactions (attraction, repulsion, scattering, merging, annihilation).

### Procedure

For each distinct stable pattern type identified in Stage 1, run a pairwise interaction sweep:
- Vary initial separation (1 to L/4)
- Vary relative phase/orientation
- Run long enough to observe interaction outcome

### Measurement

- Effective interaction force as function of separation
- Scattering angles (if patterns move)
- Whether interactions conserve a recognizable quantity (pattern count, total amplitude, phase sum)

### Falsifier

If interactions are trivial — patterns pass through, annihilate immediately, or ignore each other — no physics is emerging from the rule.

---

## Stage 3 — Mass-Like (Inertial) Behavior

Triggered only if Stage 2 succeeds.

### Hypothesis

Stable patterns exhibit **inertia**: response to external forcing produces acceleration proportional to forcing, with the proportionality constant ("inertial mass") being a property of the pattern.

### Procedure

Apply controlled external forcing (gradient in the field, time-varying source term) to a stable pattern. Measure response.

### Falsifier

If patterns respond instantly to forcing (no inertia) or don't respond at all, no F=ma-like behavior is emerging.

---

## Stage 4 — Field-Like (Long-Range) Behavior

Triggered only if Stage 3 succeeds.

### Hypothesis

Patterns produce or respond to **long-range influence** through the grid — analog of fields. Different "charge" classes of patterns interact via different field types.

### Procedure

Identify pattern internal-structure variants that could function as different "charges." Measure interaction strength as a function of distance for like and unlike charges.

### Falsifier

If all interactions are short-range (decay faster than any power law), no field-like phenomena are emerging.

---

## Stage 5 — Quantum-Like Interference

Triggered only if Stage 4 succeeds.

### Hypothesis

Single-pattern propagation through a "double-slit" analog produces an **interference pattern** at a detector array.

### Procedure

Construct a barrier in the grid with two openings; launch a stable pattern toward the barrier; record arrival distribution at detector.

### Falsifier

If arrival distribution shows no interference (purely classical), the framework's quantum claims fall apart at this scale.

---

## Output structure

Each fleet run writes to a partitioned registry. Public results (committed to this repo) summarize per-stage findings; raw data lives in `shared-context/synchronism/exploration-results/cadgpc/` (private).

```
explorations/results/cadgpc/
├── stage1/
│   ├── RESULTS.md          # Conclusion + headline findings
│   ├── per-mechanism/      # Per saturation mechanism analysis
│   └── pattern-catalog/    # Distinct stable patterns identified
├── stage2/  (only if Stage 1 succeeded)
...
```

## Open methodological questions

These need to be settled before Stage 1 simulator code is written:

1. **Rule family bounds**: How rich a rule family is "in scope"? (Number of coefficients, neighborhood size, state-space dimension.) Too narrow and we miss the right rule; too wide and the sweep is intractable.
2. **What counts as "stable"**: Single quantitative threshold or multi-criterion? (Centroid drift + amplitude variance + Fourier-peak-above-noise + perturbation-recovery — and how do we combine?)
3. **What counts as "different pattern types"**: Topological invariants? Statistical features of the local field? Need an algorithm, not a vibe.
4. **Reporting back to the registry**: What summary statistics per parameter point are sufficient to reconstruct the conclusion without shipping all runs? (Cells × timesteps × runs gets huge fast.)

These are the questions a Stage 1 simulator design doc (next deliverable) needs to answer.

## Status

- **2026-05-15**: Drafted. Stage 1 simulator design doc is the next deliverable. Fleet directive `private-context/directives/fleet-dir-cadgpc-idle-compute-2026-05-15.md` queues fleet for execution once the simulator is ready.
- After Stage 1 closes (success or failure): RESULTS.md written, findings published, next stage triggered if applicable.

## Why this matters

The framework's strongest external critique to date converges on one ask: **stop relabeling existing physics; produce new physics from the discrete-grid ontology, or accept that the ontology is conceptual rather than physical**. CADGPC is the experiment that addresses this directly. Either:

- It works → the discrete-grid ontology earns the right to be considered physics, not just description
- It fails → the framework needs revision, and we have a concrete constraint on what the revision must address

Either outcome is more valuable than the framework's current state (zero novel predictions, ambiguous physics-vs-philosophy positioning). The exploration is the operationalization of "interactive selection" applied to the framework itself — we'll see what the discrete grid actually does, not what we'd like it to do.
