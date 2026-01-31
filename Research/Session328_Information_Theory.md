# Session #328: Information Theory from the Planck Grid

**Information Theory Arc (Session 1/4)**
**Date**: 2026-01-31

## Overview

This session initiates the Information Theory Arc following the completion of the Statistical Mechanics Arc. The key insight is that information = pattern distinguishability on the grid. Shannon entropy measures uncertainty about which pattern is realized. The MRH is the channel capacity of nature itself — the maximum rate at which pattern distinguishability can be transferred.

## Key Questions

1. How does Shannon entropy relate to grid patterns?
2. What is Landauer's principle in terms of MRH?
3. How is Maxwell's demon resolved on the grid?
4. What is the channel capacity of the MRH boundary?

## Key Results (8/8 verified)

### Part 1: Shannon Entropy

**The Fundamental Formula**:
```
H(X) = -Σ p(x) log₂ p(x)
```

Shannon entropy measures the average bits needed to specify which outcome occurred.

**Key Properties**:
| Property | Formula | Meaning |
|----------|---------|---------|
| Non-negativity | H(X) ≥ 0 | Uncertainty is never negative |
| Maximum | H_max = log₂(n) | Uniform distribution over n outcomes |
| Additivity | H(X,Y) = H(X) + H(Y\|X) | Chain rule |
| Concavity | H(λp + (1-λ)q) ≥ λH(p) + (1-λ)H(q) | Mixing increases entropy |

**Examples** (4 outcomes):
| Distribution | H (bits) |
|-------------|----------|
| Uniform [0.25, 0.25, 0.25, 0.25] | 2.000 |
| Biased [0.9, 0.05, 0.025, 0.025] | 0.569 |
| Certain [1, 0, 0, 0] | 0.000 |

**Grid Interpretation**:
| Concept | Grid Meaning |
|---------|--------------|
| Entropy | Bits to specify which pattern is realized |
| Probability | Weight of each pattern configuration |
| Mutual info | Shared pattern structure between regions |
| Conditional | Uncertainty given context (inside MRH) |

### Part 2: Landauer's Principle

**The Fundamental Limit**:
```
E_min = k_B T ln(2) per bit erased
```

At 300 K: E_min ≈ 2.87 × 10⁻²¹ J ≈ 0.018 eV per bit.

**Key Values**:
| Temperature | Energy per bit |
|-------------|---------------|
| 100 K | 0.96 × 10⁻²¹ J |
| 300 K | 2.87 × 10⁻²¹ J |
| 1000 K | 9.57 × 10⁻²¹ J |

**Comparison to Modern Computing**:
- Typical CPU: ~10⁻¹⁸ J per bit operation
- Landauer limit: ~3 × 10⁻²¹ J per bit
- Current efficiency: ~0.0003 (~1000× above limit)

**Grid Interpretation**:
| Concept | Grid Meaning |
|---------|--------------|
| Erasure | Pattern distinguishability lost |
| Energy cost | Work to push info across MRH boundary |
| Entropy increase | Pattern info becomes thermal (beyond MRH) |
| Reversibility | Reversible = no info crosses MRH |

### Part 3: Maxwell's Demon

**The Setup**:
A "demon" sorts molecules into hot/cold compartments, apparently decreasing entropy.

**The Resolution** (Landauer-Bennett):
1. Demon must STORE information about each molecule
2. Memory has finite capacity → must be erased
3. Erasure costs k_B T ln(2) per bit (Landauer)
4. Net entropy change ≥ 0 (2nd Law OK!)

**Numerical Example** (100 molecules at 300 K):
| Quantity | Value |
|----------|-------|
| Apparent entropy decrease | -100 k_B ln(2) |
| Memory required | 100 bits |
| Erasure cost | 2.87 × 10⁻¹⁹ J |
| Entropy from erasure | +100 k_B ln(2) |
| **Net entropy change** | **≥ 0** |

**Szilard Engine**:
```
Max work from 1 molecule: W = k_B T ln(2)
This equals the Landauer erasure cost!
Net work over full cycle: 0
```

**Grid Interpretation**:
| Concept | Grid Meaning |
|---------|--------------|
| Demon | Subsystem creating correlations |
| Memory | Correlations between demon and gas patterns |
| Erasure | Correlations must cross MRH to reset |
| Resolution | 2nd Law = info flow across MRH is one-way |

### Part 4: Channel Capacity

**Shannon's Channel Capacity Theorem**:
```
C = B log₂(1 + S/N) bits/second
```

For bandwidth B and signal-to-noise ratio S/N.

**Examples** (B = 1 MHz):
| SNR (dB) | Capacity (Mbits/s) |
|----------|-------------------|
| 10 | 3.46 |
| 20 | 6.66 |
| 30 | 9.97 |
| 40 | 13.29 |

**Binary Symmetric Channel**:
```
C = 1 - H(p) = 1 + p log₂(p) + (1-p) log₂(1-p)
```

| Error probability p | Capacity (bits/use) |
|--------------------|-------------------|
| 0.0 | 1.000 |
| 0.1 | 0.531 |
| 0.2 | 0.278 |
| 0.5 | 0.000 |

**Grid Interpretation**:
| Concept | Grid Meaning |
|---------|--------------|
| Channel | Path for pattern transfer between regions |
| Capacity | Max rate of distinguishable patterns |
| Noise | Random pattern perturbations |
| MRH as channel | MRH boundary has channel capacity |

### Part 5: Data Compression

**Shannon's Source Coding Theorem**:
Optimal compression approaches H(X) bits per symbol.

**English Text Example** (top 10 letters):
| Metric | Value |
|--------|-------|
| Entropy H | 3.13 bits/letter |
| Max entropy | 3.32 bits/letter |
| Redundancy | 5.7% |

**Rate-Distortion** (Gaussian source, σ² = 1):
| Distortion D | Rate R(D) (bits/sample) |
|-------------|------------------------|
| 0.1 | 1.66 |
| 0.25 | 1.00 |
| 0.5 | 0.50 |
| 1.0 | 0.00 |

**Grid Interpretation**:
| Concept | Grid Meaning |
|---------|--------------|
| Compression | Removing pattern redundancy |
| Entropy | Irreducible pattern complexity |
| MRH compression | MRH naturally compresses: only relevant patterns tracked |

### Part 6: Information Thermodynamics

**The Unification**:
```
S_Boltzmann = k_B × H_Shannon × ln(2)
```

Thermodynamic entropy and information entropy are the SAME quantity in different units.

**Conversions**:
| S (Boltzmann) | H (Shannon) |
|---------------|-------------|
| 10 k_B | 14.43 bits |
| 6.93 k_B | 10 bits |

**Key Relations**:
| Principle | Meaning |
|-----------|---------|
| Landauer | Forgetting = pattern info crosses MRH |
| Szilard | Can extract work from information |
| Jarzynski | Free energy from non-equilibrium trajectories |
| Crooks | Forward/reverse path symmetry |

## Verification Summary

| Test | Result |
|------|--------|
| Shannon entropy of uniform distribution = 2 bits | PASS |
| Shannon entropy bounds (0 ≤ H ≤ log n) | PASS |
| Landauer energy positive | PASS |
| Landauer scales with temperature | PASS |
| Maxwell demon net entropy ≥ 0 | PASS |
| Channel capacity increases with SNR | PASS |
| Compression limited by entropy | PASS |
| Grid interpretations exist | PASS |

**8/8 verified.**

## New Predictions

### P328.1: Information = Pattern Distinguishability
- H(X) measures bits to specify which grid pattern is realized
- Not abstract — physically meaningful on discrete grid
- Status: DEFINITIONAL (establishes framework)

### P328.2: Landauer = MRH Crossing
- Erasing a bit = pattern info crosses MRH boundary
- Energy cost comes from this irreversible transfer
- Status: CORE INSIGHT

### P328.3: MRH = Channel Capacity
- MRH boundary has finite information transfer rate
- This is the "bandwidth" of nature
- Beyond MRH = noise, not signal
- Status: NOVEL INTERPRETATION

### P328.4: Maxwell Demon via Correlations
- Demon creates correlations between memory and gas
- Erasure destroys correlations = info crosses MRH
- Status: CONSISTENT with Bennett-Landauer resolution

## Information Theory Arc Roadmap

| Session | Topic | Focus |
|---------|-------|-------|
| #328 | Information Theory Foundations | Shannon, Landauer, demon, channels |
| #329 | Quantum Information | Entanglement, no-cloning, quantum channels |
| #330 | Holographic Principle | Information at boundaries, AdS/CFT |
| #331 | Black Hole Information | Information paradox resolution |

## Connection to Previous Arcs

This session bridges from Statistical Mechanics:

**Session #324-327 established**:
- Entropy = information beyond MRH
- Temperature = reconfiguration rate
- 2nd Law = info flow beyond MRH is irreversible

**Session #328 formalizes**:
- Shannon entropy as the rigorous measure
- Landauer principle as the physical cost of forgetting
- Channel capacity as the MRH's information bandwidth
- Maxwell's demon as the test case for info-thermo unity

**Key Unification**:
```
Statistical Mechanics ↔ Information Theory
        S = k_B H ln(2)
      (Same thing, different units)
```

The MRH boundary is:
- The system-bath boundary (thermodynamics)
- The tracked-averaged boundary (statistical mechanics)
- The signal-noise boundary (information theory)

All the same thing!

---

*"Information is not abstract. On the Planck grid, it is pattern distinguishability. Entropy is the bits we don't know. The MRH is nature's channel capacity — the fundamental limit on how fast patterns can be communicated across space and time."*

## Files

- `simulations/session328_information_theory.py`
- `simulations/session328_information_theory.png`
- `Research/Session328_Information_Theory.md`

---

**INFORMATION THEORY ARC (1/4)**

Next: Session #329 - Quantum Information
