# Coupling-Coherence Experiment: Testing the Generalized Coherence Equation

**Date**: 2026-02-22
**Motivation**: External challenge from Andrei's AI (GPT) — "show me one worked AI example with a measurement protocol"
**Status**: Designed and executed

---

## 1. Research Question

Does the synchronism coherence equation

```
C(ρ) = tanh(γ · log(ρ/ρ_crit + 1))
```

describe the emergence of coherence in multi-agent knowledge discovery when ρ is generalized from physical density to coupling density between agents?

Specifically:
1. Does C(p) follow a tanh(γ · log(p/p_crit + 1)) curve, or does a logistic/erf/power-law fit equally well?
2. Can p_crit be **derived** from world properties (like ρ_crit = A × V_flat² in physics), or is it only a fit parameter?
3. Can we distinguish "shared wrongness" (convergence without correctness) from genuine coherence?

## 2. Experimental Design

### 2.1 Micro-World (The Oracle)

A random directed knowledge graph serves as ground truth:
- **N** = 12 nodes (entities)
- **E** = 30 typed edges (relationships)
- **T** = 3 relationship types
- World has 12×11×3 = 396 possible edge positions
- Generated as random directed graph with specified edge count
- The graph IS the oracle — agents are correct insofar as their beliefs match it
- **Parameter rationale**: Individual observation budget (640 obs) covers ~1.6 observations per edge position — insufficient to learn the world alone. Collective budget (3,200 obs with K=5) gives ~8 per position — sufficient with coupling. This creates the regime where coupling *matters*.

### 2.2 Agent Model

**K = 5** Bayesian belief agents (not LLMs — removes shared-prior confounds):

- Each agent maintains a belief matrix **B[i][j][t]** ∈ [0,1] representing P(edge (i→j) of type t exists)
- **Prior**: B = 0.5 for all possible edges (uninformative)
- **Observation**: Each round, agent observes **m = 8** randomly selected true edges, subject to noise
- **Noise model**: With probability **η = 0.15**, an observation is flipped (true edge reported absent, or absent edge reported present)
- **Bayesian update**: Standard likelihood ratio update with known noise rate

### 2.3 Coupling Mechanism

Coupling parameter **p ∈ [0, 1]** controls information sharing:

- Each round, for each ordered pair of agents (A, B), with probability p: A receives B's current belief matrix
- When received: A updates via **weighted belief averaging**: B_new = α·B_self + (1-α)·B_received, where α = 0.7 (self-weight)
- **p = 0**: Fully independent (no sharing)
- **p = 1**: Full transparency every round
- **Units**: Sharing probability per agent-pair per round

### 2.4 Measurements

Per round, per coupling level:

| Metric | Formula | What It Captures |
|--------|---------|-----------------|
| **C_conv** (Convergence) | Mean pairwise Jaccard similarity of agents' predicted edge sets | Inter-agent agreement |
| **C_corr** (Correctness) | Mean F1 score of agents' thresholded beliefs vs ground truth | Agreement with oracle |
| **C** (Coherence) | √(C_conv × C_corr) | Genuine coherence (not shared wrongness) |
| **MI** | Mutual information between agent belief distributions | Information coupling |

The geometric mean for C is critical: high C_conv + low C_corr = "converged but wrong" → C stays low.

### 2.5 Run Parameters

| Parameter | Primary | Variations |
|-----------|---------|------------|
| Coupling levels | 45: fine near 0 (step 0.005 for p ≤ 0.10), coarser above | — |
| Rounds per run | R = 80 | — |
| Agents | K = 5 | 3, 5, 10 |
| Repetitions | 20 random worlds per coupling level | — |
| Observations/round/agent | m = 8 | 2, 5, 10 |
| Noise rate | η = 0.15 | 0.05, 0.15, 0.30 |
| World size | N = 12 nodes, E = 30 edges, T = 3 types | N ∈ {8, 12, 20} |
| **Total primary runs** | 45 × 20 = 900 | |

## 3. Analysis Protocol

### 3.1 Curve Fitting

Four sigmoid models, all with 2 free parameters (fair comparison):

| Model | Formula | Params |
|-------|---------|--------|
| **tanh** | C = tanh(γ · log(p/p_crit + 1)) | γ, p_crit |
| **logistic** | C = 1/(1 + exp(-k·(p - p_half))) | k, p_half |
| **erf** | C = erf(k · (p - p_half)) | k, p_half |
| **Hill** | C = p^k / (p^k + p_half^k) | k, p_half |

Fit via scipy.optimize.curve_fit with L-BFGS-B. Report RSS, R², AIC, BIC for each.

### 3.2 Change-Point Detection

- Compute numerical second derivative d²C/dp² from the averaged C(p) curve
- Find p* = argmax |d²C/dp²| (point of maximum curvature change)
- Compare p* with p_crit from tanh fit
- If |p* - p_crit| / p_crit < 0.15, the tanh model's p_crit is consistent with the empirical transition

### 3.3 Derived p_crit Attempt

**Candidate formula**:
```
p_crit_derived = η · H(world) / (K · m · (1 - 2η))
```

Where:
- η = noise rate (known)
- H(world) = entropy of ground-truth graph in bits = E · log₂(T) + (N²·T - E) · log₂(1) ≈ E · log₂(T)
- K = agent count
- m = observations per round
- (1 - 2η) = effective information per observation

**Intuition**: p_crit is the coupling where collective signal rate exceeds world complexity.

**Test**: Run variations (different K, η, N, m). For each, compute both fitted p_crit and derived p_crit. Plot derived vs fitted. If R² > 0.8, the derivation has explanatory power.

## 4. Kill Criteria

| Criterion | What It Would Mean |
|-----------|-------------------|
| Logistic or erf beats tanh by ΔAIC > 2 | The synchronism functional form is not preferred |
| p_crit derivation R² < 0.5 across variations | p_crit is a fit parameter, not a derived invariant |
| C_conv >> C_corr at high coupling | "Shared wrongness" dominates — coherence ≠ correctness |
| No sigmoid transition (linear C vs p) | No phase-transition-like behavior at all |

All of these are publishable negative results. The experiment is designed to produce useful information regardless of outcome.

## 5. Implementation

- **Simulation**: `simulations/coupling_coherence_experiment.py`
- **Analysis**: `simulations/coupling_coherence_analysis.py`
- **Results**: `simulations/results/coupling_coherence_*.json`
- **Presentation**: `synchronism-site/src/app/coupling-experiment/page.tsx`

## 6. Relationship to Physics Framework

| Physics | This Experiment | Analogy |
|---------|----------------|---------|
| ρ (density, g/cm³) | p (coupling probability) | Input to coherence function |
| ρ_crit = A × V_flat² | p_crit = η·H/(K·m·(1-2η)) | Derived from system properties |
| γ = 2/√N_corr | γ (fit parameter, then test derivation) | Transition sharpness |
| C(ρ) = tanh(γ·log(ρ/ρ_crit+1)) | C(p) = tanh(γ·log(p/p_crit+1)) | Same functional form? |

The honest difference: in physics, ALL parameters are derived from fundamental constants. Here, γ is initially fit, and we attempt to derive p_crit. This is a weaker claim — "same curve family with one derived parameter" vs "fully derived prediction." We report this honestly.

## 7. Results (2026-02-22)

### 7.1 Summary

900 primary runs (45 coupling levels × 20 repetitions) plus 12 variation experiments (252 additional runs each).

### 7.2 Key Findings

| Finding | Detail |
|---------|--------|
| **Clear sigmoid transition** | C ranges from 0.34 (p=0) to 0.94 (p=1), centered around p ≈ 0.03-0.06 |
| **Hill function wins** | Hill (R²=0.880) beats tanh (R²=0.868) by ΔAIC=4.0 |
| **tanh competitive but not preferred** | tanh beats logistic (ΔAIC=-13.4) and erf (ΔAIC=-12.8), loses to Hill |
| **p_crit derivation fails** | Derived p_crit ≈ 0.82, fitted p_crit ≈ 0.002 — off by 400× |
| **No shared wrongness** | Max(C_conv - C_corr) = 0.128 — convergence tracks correctness |
| **Very small p_crit** | Transition begins immediately — even tiny coupling helps dramatically |

### 7.3 Kill Criteria Assessment

| Criterion | Result | Status |
|-----------|--------|--------|
| Logistic/erf beats tanh | Hill beats tanh (ΔAIC=4); logistic/erf do not | **PARTIALLY TRIGGERED** |
| p_crit derivation R² < 0.5 | R² = -662 (catastrophic failure) | **TRIGGERED** |
| C_conv >> C_corr | Max gap 0.128 | **CLEARED** |
| No sigmoid transition | C range 0.598, clear sigmoid | **CLEARED** |

### 7.4 Interpretation

1. **The sigmoid transition is real.** Multi-agent knowledge discovery with controllable coupling produces a clear sigmoid C(p) curve. This is not linear and not noise.

2. **The tanh form is competitive but not uniquely preferred.** A Hill function (power-law sigmoid) fits slightly better. Both capture the key shape: rapid rise from low coupling, then saturation. The difference is in the low-coupling regime where the tanh's logarithmic argument creates a slightly different curvature than the Hill function's power law.

3. **The p_crit derivation is wrong.** The information-theoretic formula p_crit = η·H/(K·m·(1-2η)) overestimates by 400×. This means p_crit in this domain is **only a fit parameter**, not a derived invariant. The "universal" claim shifts from "derived invariant" to "nice curve family" — exactly as Andrei predicted.

4. **Coupling is powerful at very low levels.** p_crit ≈ 0.002 means that even a 0.5% chance of any two agents sharing beliefs each round significantly improves collective coherence. This has practical implications for multi-agent system design.

5. **The experiment successfully separates convergence from correctness.** The geometric mean C = √(C_conv × C_corr) effectively detects "shared wrongness" — this measurement design is sound.

### 7.5 Honest Assessment

The experiment achieved what Andrei asked for: a worked example with ground truth, controllable coupling, measured convergence + correctness, and a derived p_crit test. The result is mixed:

- **Positive**: Clear sigmoid transition, convergence tracks correctness, measurement protocol works
- **Negative**: tanh not uniquely preferred (Hill wins), p_crit derivation fails completely

This downgrades the claim from "universal equation with derived parameters" to "one of several equally good sigmoid descriptions of a real phenomenon." Still useful — but a different and weaker claim than originally intended.
