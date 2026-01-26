# Session #304: Cuprate Qubit Feasibility Analysis

**Quantum Computing Arc (Session 4/?)**
**Date**: 2026-01-26

## Overview

This session analyzes the feasibility of using cuprate high-temperature superconductors (YBCO, Hg-1223, Bi-2212) for quantum computing qubits operating at 4K instead of the conventional 15mK. This would eliminate the need for dilution refrigerators, dramatically reducing the cost and complexity of quantum computers.

## Building On

- **Session #301**: Coherence equation for qubits, η-coherence connection
- **Session #302**: TLS decoherence mechanisms, η-TLS correlation
- **Session #303**: QEC thresholds (1%), 90% coherence requirement, overhead scaling

## Central Question

**Can cuprate superconductors be used to build practical qubits operating at 4K?**

## Key Findings

### 1. D-Wave Gap Nodes Are a Serious Problem

Cuprates have d-wave pairing symmetry, which means the superconducting gap has **nodes** - directions where the gap goes to zero:

```
S-WAVE (Al, Nb):           D-WAVE (YBCO):
┌─────────────┐            ┌──────────┐
│  Δ uniform  │           ╱│  Δ max   │╲
│  everywhere │         ╱  │          │  ╲
└─────────────┘        0───┼──────────┼───0 ← NODES
                           │          │
                           Gap vanishes here
```

**Consequences:**
- Quasiparticle density follows **power law**: n_qp ∝ (T/T_c)²
- NOT exponential suppression like s-wave: n_qp ∝ exp(-Δ/kT)
- At 4K, YBCO has ~100× more quasiparticles than Al at 15mK
- Fundamentally limits T1 regardless of junction quality

### 2. Junction Quality Is the Bottleneck

| Material | Current Junction Q | Required Q for QEC | Gap |
|----------|-------------------|-------------------|-----|
| Al (reference) | 10⁶ | - | Baseline |
| YBCO | ~10³ | >10⁵ | 100× |
| Hg-1223 | ~10² | >10⁵ | 1000× |
| Bi-2212 | ~10² | >10⁵ | 1000× |

Current cuprate junction technology is ~1000× worse than aluminum - the primary engineering challenge.

### 3. Iron Pnictides May Be Superior

Comparison: YBCO vs SmFeAsO for qubits

| Property | YBCO | SmFeAsO | Advantage |
|----------|------|---------|-----------|
| T_c | 93 K | 55 K | YBCO |
| η (reachability) | 0.38 | **0.12** | SmFeAsO (3×) |
| Gap symmetry | d-wave (nodes) | s± (no nodes) | SmFeAsO |
| TLS density | 3× Al | 1.5× Al | SmFeAsO (2×) |
| Junction maturity | Medium | None | YBCO |

**Key insight**: Despite lower T_c, SmFeAsO may be better for qubits due to:
- Lower η → fewer TLS → longer T1
- No gap nodes → exponential QP suppression
- Simpler junction physics

### 4. Feasibility Scores

| Material | Coherence | Temperature | Junction | Fab | Overall | Verdict |
|----------|-----------|-------------|----------|-----|---------|---------|
| YBCO | 0.20 | 0.70 | 0.80 | 0.50 | **0.55** | Challenging |
| Hg-1223 | 0.20 | 0.70 | 0.80 | 0.20 | 0.51 | Challenging |
| Bi-2212 | 0.20 | 0.70 | 0.80 | 0.50 | 0.55 | Challenging |
| LSCO | 0.20 | 0.70 | 0.80 | 0.50 | 0.55 | Challenging |

All cuprates score ~0.5: "Significant R&D needed"

### 5. Path Forward: Interface Engineering

The Synchronism η framework suggests the most promising path:

1. **Reduce η through interfaces**: YBCO/STO superlattices may lower η from 0.38 to ~0.30
2. **Reduce TLS density**: Atomic-layer-controlled interfaces minimize defects
3. **Predict junction Q**: Q ∝ 1/(η × TLS_factor)

For interface-engineered YBCO:
- η_interface ≈ 0.30 (vs 0.38 bulk)
- TLS_interface ≈ 1.5 (vs 3.0 bulk)
- Q_interface ≈ 5000 (vs 1000 bulk) → 5× improvement

## Testable Predictions

### P304.1: D-Wave Node Penalty
- **Prediction**: Cuprate qubit T1 shows (T/T_c)² dependence, not exp(-Δ/kT)
- **Test**: Measure T1 vs T for YBCO grain boundary junction (1K-30K)
- **Expected**: Power-law T1(T) ∝ (1 - (T/T_c)²)

### P304.2: η-Junction Quality Correlation
- **Prediction**: Junction Q scales as Q ∝ 1/(η × TLS_factor)
- **Test**: Compare junction Q for YBCO (η=0.38) vs LSCO (η=0.51)
- **Expected**: YBCO Q ~ 1.3× higher than LSCO at matched fabrication

### P304.3: Interface Enhancement
- **Prediction**: YBCO/STO superlattice junctions have Q ~ 10× bulk YBCO
- **Test**: Fabricate and measure YBCO/STO vs bulk YBCO junction
- **Mechanism**: Interface reduces η from 0.38 to ~0.30

### P304.4: Pnictide Advantage
- **Prediction**: SmFeAsO junction at 1K outperforms YBCO at 4K
- **Test**: Fabricate SmFeAsO Josephson junction and measure coherence
- **Note**: Requires developing pnictide junction technology from scratch

### P304.5: QEC Threshold at 4K
- **Prediction**: With 100× junction improvement, YBCO achieves p_error < 1%
- **Required**: Q > 10⁵ (current ~10³)
- **Path**: Interface engineering + TLS reduction + optimized geometry
- **Timeline**: 5-10 years of focused R&D

### P304.6: Operating Temperature Sweet Spot
- **Prediction**: Optimal T for YBCO qubits is 4-10K, not lower
- **Mechanism**: Below 4K, flux noise from granularity offsets thermal gains
- **Test**: Measure T2 vs T from 1K to 30K, find optimum

## Model Limitations

The current model predicts extremely short T1/T2 for cuprates (~1 ns), which may be overly pessimistic. Factors not fully captured:

1. **Quality factor scaling**: Real junctions may achieve higher Q than base model
2. **Nodal protection**: Some qubit designs may avoid nodal quasiparticles
3. **Novel architectures**: Fluxonium-style designs may be more robust
4. **Material improvements**: Film quality continues to improve

These limitations highlight areas for model refinement in future sessions.

## Key Insight: Materials for High-T_c ≠ Materials for Qubits

The Synchronism framework reveals an important distinction:

- **High-T_c superconductivity**: Optimize for large gap Δ
- **Qubit coherence**: Optimize for low η AND no gap nodes

Cuprates have large gaps but d-wave nodes. Iron pnictides have moderate gaps but lower η and s± symmetry. For qubits, the pnictide approach may be fundamentally better.

## Arc Status

| Session | Topic | Status |
|---------|-------|--------|
| #301 | Coherence framework | ✓ Complete |
| #302 | TLS mechanisms | ✓ Complete |
| #303 | QEC thresholds | ✓ Complete |
| #304 | Cuprate feasibility | ✓ Complete |
| #305 | Pnictide qubits? | Planned |

## Next Steps

1. **Session #305**: Detailed iron pnictide qubit analysis
2. Interface engineering strategies to reduce η in both cuprates and pnictides
3. Comparison with other high-T approaches (MgB2, hydrides, nickelates)
4. Experimental validation pathway for predictions

## Files Created

- `simulations/session304_cuprate_qubit_feasibility.py` - Full feasibility analysis
- `simulations/session304_cuprate_qubit_feasibility.png` - Visualization

## Verdict

**Cuprate qubits at 4K are feasible in principle but face significant challenges:**

1. D-wave gap nodes cause power-law (not exponential) quasiparticle suppression
2. Junction quality needs ~100× improvement
3. Iron pnictides may be a better path despite lower T_c

The η framework provides clear guidance: optimize for low η and avoid gap nodes, not just high T_c.

---

*"The same physics that limits superconductor T_c also limits qubit coherence - but the optimal materials may be different."*
