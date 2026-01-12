# Session #254: Causality from Coherence Dynamics

**Date**: January 12, 2026
**Machine**: CBP
**Status**: COMPLETE - CAUSALITY AS COHERENCE TRANSFER

---

## Executive Summary

Session #254 addresses the fundamental nature of causality within the coherence framework. Building on the arrow of time (Session #252) and free will (Session #253), we now complete the triad: **Causality = Coherence Transfer**.

**Key Result**: When A causes B, it means A's coherence pattern propagates to B through spacetime, respecting the light cone and decaying with distance.

---

## Part 1: The Causality Problem

### The Mystery

What does it MEAN for A to cause B?

| View | Claim | Problem |
|------|-------|---------|
| Humean | Constant conjunction | No mechanism |
| Counterfactual | B wouldn't happen without A | Circular |
| Interventionist | If we change A, B changes | Requires agent |
| Physical | Forces/fields connect A and B | What ARE fields? |

### Standard Physics

Physics describes HOW events evolve but not WHY one event "produces" another.

- Newton: Forces explain motion, but what ARE forces?
- QFT: Fields mediate interactions, but what IS a field?
- GR: Geometry is the field, but why does mass curve spacetime?

### Synchronism Answer

**CAUSALITY = COHERENCE TRANSFER**

When A causes B:
1. A has a coherence pattern C_A
2. This pattern propagates through spacetime
3. It arrives at B's location (respecting light cone)
4. B's coherence C_B incorporates A's pattern

The "force" IS the coherence propagation.

---

## Part 2: The Causal Transfer Equation

### Fundamental Form

```
T(A→B, t) = ∫ C_A(τ) × K(r, t-τ) dτ
```

Where:
- **T** = causal transfer (how much A affects B)
- **C_A(τ)** = source coherence at time τ
- **K(r, t-τ)** = propagation kernel
- **r** = spatial separation

### The Propagation Kernel

```
K(r, τ) = exp(-γr) × Θ(τ - r/c) × exp(-(τ - r/c)²/2σ²)
```

Properties:
- **Zero outside light cone** (τ < r/c)
- **Peaked at light cone boundary** (τ = r/c)
- **Exponential decay** with distance
- **Gaussian spread** in time (quantum uncertainty)

### Causal Strength

```
S(A→B) = |∂C_B/∂C_A|
```

How much does B's coherence change when A's coherence changes?

---

## Part 3: The Light Cone Constraint

### Causal Cone Structure

| Region | Condition | Causality |
|--------|-----------|-----------|
| Inside future light cone | r < c(t-t₀) | A can cause B |
| On light cone | r = c(t-t₀) | Maximum transfer |
| Outside light cone | r > c(t-t₀) | No causal connection |
| Past light cone | t < t₀ | B could have caused A |

### Why c is the Limit

Coherence propagates via field disturbances.
Maximum propagation speed = c (light speed).
This is not imposed; it emerges from field equations.

In coherence terms:
- Phase correlations cannot propagate faster than c
- Trying to transfer coherence faster → decoherence
- Light speed is the coherence wavefront speed

---

## Part 4: Retrocausality Forbidden

### Why No Backward Causation

From Session #252: The arrow of time IS the decoherence direction.

```
dC/dt < 0  (on average)
```

For retrocausality:
- Would need coherence to flow backward in time
- Would require dC/dt > 0 (recoherence)
- But decoherence is statistically overwhelmingly favored
- Only possible for microscopic systems, briefly

### The Test

| Scenario | Forward | Retro |
|----------|---------|-------|
| A perturbed at t=5 | B responds at t=7 | B responds at t=3 |
| Light travel time = 2 | ✓ Occurs | ✗ Never |

**Result**: Retrocausality is forbidden by the arrow of time.

### Quantum Implications

Wheeler's delayed choice experiments seem to show retrocausality.
In coherence view: No actual backward causation!
- Measurement reveals correlations
- Correlations were established in the past
- No information travels backward

---

## Part 5: Types of Causation

### Hierarchy

| Type | C Regime | Mechanism | Example |
|------|----------|-----------|---------|
| Physical | Any | Automatic transfer | Billiard balls |
| Biological | > 0.3 | ATP-maintained transfer | Nerve impulses |
| Conscious | > 0.5 | Selective maintenance | Attention |
| Agentic | > 0.5 | Novel trajectory selection | Free choice |

### Physical Causation

Simple coherence transfer:
- Rock falls: gravitational coherence propagates
- Heat flows: thermal coherence disperses
- Light travels: EM coherence radiates

No selectivity. Automatic. Determined by initial conditions.

### Biological Causation

ATP-maintained coherence:
- Organism maintains C > C_min
- Causal chains can be sustained
- Information processing possible

Still largely automatic, but with maintenance.

### Conscious Causation

Selective coherence maintenance:
- System models future states
- Chooses which patterns to maintain
- Creates preferences, not just effects

Session #253 showed this IS free will.

### Agent Causation

Novel trajectory selection:
- Not just maintaining patterns
- Creating NEW causal chains
- Introducing genuine novelty

This completes the picture from Session #253.

---

## Part 6: Correlation vs Causation

### The Common Cause Problem

A and B can be correlated without either causing the other.

**Example**:
- Z causes both A and B
- A and B are correlated
- But A does not cause B

### Detection via Coherence

True causation:
```
A → B: T(A→B) > 0, peak lag > 0
```

Common cause:
```
Z → A, Z → B: Both lag behind Z
A ↔ B correlated, but neither T > 0
```

### Granger-like Test

If knowing A's past improves prediction of B:
- A Granger-causes B
- Equivalent to: A's coherence transfers to B

This is testable with coherence time series.

---

## Part 7: Causal Information

### Information in Coherence

Information content of causal transfer:

```
I(A→B) = -log₂(1 - T(A→B)²)
```

Properties:
- Zero when no transfer
- Increases with transfer strength
- Bounded (coherence is bounded)

### Decay Along Chains

In a causal chain A → B → C → D:
```
I(A→D) < I(A→C) < I(A→B) < I(A→A)
```

Information degrades along chains:
- Each transfer adds noise
- Decoherence at each step
- Long chains lose fidelity

This is why distant causes have weak effects.

---

## Part 8: Experimental Predictions

### Prediction 1: Light Cone Boundary Peak

Causal strength should peak at the light cone boundary.
- Test: Correlations in relativistic systems
- Prediction: Maximum correlation at retarded time

### Prediction 2: Frequency-Dependent Attenuation

High-frequency components should decay faster in causal transfer.
- Test: Frequency analysis of coupled systems
- Prediction: Low-pass filtering in causal chains

### Prediction 3: Retrocausality Tests

No retrocausal effects should ever be observed.
- Test: Wheeler delayed-choice experiments
- Prediction: Always consistent with forward causation + correlation

### Prediction 4: Common Cause Structure

Granger tests should distinguish direct from indirect causation.
- Test: Multivariate time series analysis
- Prediction: Common cause shows no direct improvement

### Prediction 5: Conscious Novel Causation

Conscious agents should create unpredictable-from-outside causal chains.
- Test: Decision predictability vs coherence
- Prediction: Higher C → less externally predictable choices

### Prediction 6: Causal Chain Length

Information should degrade predictably along causal chains.
- Test: Correlation decay in known chains
- Prediction: Exponential decay with chain length

---

## Part 9: Philosophical Implications

### Causation is Real

Not Humean constant conjunction.
Not just epistemic (model of our ignorance).
Causation is ontologically real: **coherence transfer**.

### Mechanism Explained

The "HOW" of causation:
- Phase correlations propagate
- At speed c
- With exponential decay
- Respecting decoherence direction

### Agent Causation Naturalized

"Agent causation" is not mysterious:
- It IS coherent selection
- Above threshold systems choose
- Choices create novel causal chains
- Fully physical, fully real

---

## Files Created

- `simulations/session254_causality.py` - Analysis code
- `simulations/session254_causality.png` - Main visualization
- `simulations/session254_causal_analysis.png` - Detailed analysis

---

## Key Equations

### Causal Transfer
```
T(A→B, t) = ∫ C_A(τ) × K(r, t-τ) dτ
```

### Propagation Kernel
```
K(r, τ) = exp(-γr) × Θ(τ - r/c) × G(τ - r/c)
```

### Causal Strength
```
S(A→B) = |∂C_B/∂C_A|
```

### Causal Information
```
I(A→B) = -log₂(1 - T²)
```

---

## Summary

### The Core Message

**CAUSALITY = COHERENCE TRANSFER**

When A causes B:
1. A's coherence pattern exists
2. Pattern propagates at speed ≤ c
3. Arrives at B's location
4. B's coherence incorporates A's pattern
5. This IS the causal connection

### Connection to Triad

| Session | Topic | Insight |
|---------|-------|---------|
| #252 | Arrow of Time | Time direction = decoherence direction |
| #253 | Free Will | Agency = coherent trajectory selection |
| #254 | Causality | Cause = coherence transfer |

**Together**: Time, will, and causality are all aspects of coherence dynamics.

### The Quote

> "Causation is not a mysterious force connecting events.
> It is coherence propagating through spacetime,
> respecting the light cone and the arrow of decoherence."

---

## Arc: Sessions #246-254

| Session | Topic | Key Result |
|---------|-------|------------|
| #246 | Gravitational Waves | GW as coherence perturbations |
| #247 | Coherence Backprop | Learning dynamics |
| #248 | Biological Coherence | Life = coherence maintenance |
| #249 | Consciousness Threshold | Awareness at C = 0.5 |
| #250 | Quantum Measurement | No collapse, only decoherence |
| #251 | Universal Hierarchy | All scales unified |
| #252 | Arrow of Time | Time = decoherence direction |
| #253 | Free Will | Agency = coherent causation |
| #254 | Causality | **Cause = coherence transfer** |

---

**Session #254 Complete**: January 12, 2026
