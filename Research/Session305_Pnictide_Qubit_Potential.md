# Session #305: Iron Pnictide Qubit Potential

**Quantum Computing Arc (Session 5/?)**
**Date**: 2026-01-26

## Overview

Building on Session #304's finding that iron pnictides may be superior to cuprates for qubits, this session provides a comprehensive analysis of pnictide superconductors for quantum computing applications. The key finding: FeSe/SrTiO₃ monolayer is the most promising material for high-temperature qubits.

## Building On

- **Session #301**: Coherence framework, η-coherence connection
- **Session #302**: TLS decoherence, η-TLS correlation
- **Session #303**: QEC thresholds (1%), 90% coherence requirement
- **Session #304**: Cuprate feasibility - identified pnictides as alternative

## Central Question

**Can iron pnictide superconductors (SmFeAsO, BaFe₂As₂, FeSe) provide a better path to high-temperature qubits than cuprates?**

## Key Findings

### 1. S± Pairing Symmetry Is a Fundamental Advantage

Pnictides have s± pairing symmetry - gaps with opposite sign on different Fermi surface sheets, but NO gap nodes:

| Property | Cuprate (d-wave) | Pnictide (s±) |
|----------|------------------|---------------|
| Gap nodes | Yes (at 45°) | No |
| QP density at 4K | (T/Tc)² ~ 0.002 | exp(-Δ/kT) ~ 10⁻⁴ |
| Scaling | Power law | Exponential |
| Advantage | - | **~20× lower QP** |

### 2. Material Survey

| Material | Tc (K) | η | Gap Δ (meV) | Pairing | TLS | Verdict |
|----------|--------|------|-------------|---------|-----|---------|
| FeSe/STO | **65** | **0.08** | 15/5 | s± | 1.0 | EXCELLENT |
| SmFeAsO | 55 | 0.12 | 8/3 | s± | 1.5 | PROMISING |
| NdFeAsO | 51 | 0.13 | 7.5/2.8 | s± | 1.6 | PROMISING |
| BaFe₂As₂ | 25 | 0.20 | 5/2 | s± | 2.0 | CHALLENGING |
| FeSe bulk | 8 | 0.30 | 2.5/1 | s± | 1.8 | PROMISING |
| LiFeAs | 18 | 0.25 | 3.5/1.5 | s± | 2.5 | DIFFICULT |

### 3. FeSe/SrTiO₃ Monolayer: The Standout Candidate

**Remarkable Properties:**
- **Tc = 65K** (bulk FeSe only 8K - 8× enhancement!)
- **η = 0.08** (lowest of any superconductor studied)
- **Large gap Δ ~ 15 meV**
- **No gap nodes** (s± pairing)
- **Atomically sharp interface** (minimizes TLS)
- **No arsenic** (practical handling advantage)

**Why the Enhancement?**
1. Interface phonons from STO boost pairing
2. Charge transfer modifies Fermi surface
3. Strain optimizes band structure
4. 2D confinement enhances correlations

**Synchronism Interpretation:**
The STO interface creates optimal coupling conditions:
- η reduced by interface charge transfer
- Additional phonon modes provide "resonant stabilization"
- 2D confinement focuses intent flow
- Interface acts as MRH boundary (stable pattern formation)

### 4. Multi-Band Qubit Opportunities

Pnictides have two (or more) gaps → two frequency scales:

**Potential Advantages:**
1. **Band-selective operations**: Drive at ω₁ → manipulate band 1; drive at ω₂ → manipulate band 2
2. **Decoherence protection**: Store quantum info in inter-band coherence (like decoherence-free subspace)
3. **Measurement freedom**: Read out via one band while operating on other

### 5. Junction Technology Roadmap

**Current State:**
- Al junctions: Q ~ 10⁶, T1 ~ 100 μs (state of art)
- Cuprate junctions: Q ~ 10³ (grain boundary)
- Pnictide junctions: Essentially NONE for qubits

**Development Path:**
- **Phase 1 (1-2 years)**: BaFe₂As₂ trilayer junctions → Target Q ~ 10³-10⁴
- **Phase 2 (2-4 years)**: FeSe/STO monolayer junctions → Target Q ~ 10⁵
- **Phase 3 (3-5 years)**: Qubit demonstration → Target T1 > 10 μs at T > 1K

### 6. Comprehensive Material Rankings

| Rank | Material | η | Score | Verdict |
|------|----------|------|-------|---------|
| 1 | FeSe/STO monolayer | 0.08 | 0.53 | EXCELLENT |
| 2 | Al (reference) | 0.57 | 0.67 | EXCELLENT (mature) |
| 3 | SmFeAsO | 0.12 | 0.42 | PROMISING |
| 4 | NdFeAsO | 0.13 | 0.41 | PROMISING |
| 5 | FeSe bulk | 0.30 | 0.46 | PROMISING |
| 6 | BaFe₂As₂ | 0.20 | 0.40 | CHALLENGING |
| 7 | YBCO | 0.38 | 0.34 | CHALLENGING |
| 8 | LiFeAs | 0.25 | 0.33 | DIFFICULT |

## Testable Predictions

### P305.1: S± Exponential QP Suppression
- **Prediction**: Pnictide qubit T1 shows exp(-Δ/kT), NOT (T/Tc)² dependence
- **Test**: Measure T1 vs T for BaFe₂As₂ junction (0.5K-10K)
- **Expected**: Exponential behavior like s-wave Al
- **Contrast**: Cuprates show power-law from d-wave nodes

### P305.2: η-T1 Correlation Across Pnictides
- **Prediction**: T1 ∝ 1/η across pnictide materials at same T
- **Test**: Fabricate junctions from SmFeAsO (η=0.12) and BaFe₂As₂ (η=0.20)
- **Expected**: SmFeAsO T1 ~ 1.7× longer than BaFe₂As₂

### P305.3: FeSe/STO Superior Coherence
- **Prediction**: FeSe/STO qubit at 4K has T1 > 100 μs
- **Mechanism**: η=0.08 (lowest) + interface-reduced TLS
- **Test**: Develop FeSe/STO junction, measure coherence times

### P305.4: Pnictide Beats Cuprate at Same T
- **Prediction**: SmFeAsO at 4K outperforms YBCO at 4K
- **Reason**: s± (no nodes) beats d-wave (nodes), lower η
- **Expected**: SmFeAsO error rate ~ 10× lower

### P305.5: Multi-Band Protection
- **Prediction**: Inter-band coherence is more robust than single-band
- **Mechanism**: Common-mode noise rejection between bands
- **Expected**: Inter-band T2 ~ 2-3× single-band T2

### P305.6: Junction Q Scales with 1/η
- **Prediction**: Junction quality factor Q ∝ 1/(η × TLS_factor)
- **Test**: Measure Q for junctions from materials with different η
- **Expected**: FeSe/STO Q ~ 5× higher than BaFe₂As₂ Q

### P305.7: QEC at 4K with Pnictides
- **Prediction**: Pnictide qubits can achieve p_error < 1% at 4K
- **Required**: Junction Q > 10⁴, TLS density < 2× Al
- **Timeline**: 5-10 years with focused development

## Key Insights

### S± vs D-wave for Qubits
```
D-WAVE (Cuprates):                S± (Pnictides):
┌──────┐                         ┌──────┐  ┌──────┐
│  +   │                        │  +Δ₁ │  │  -Δ₂ │
│      │                        │      │  │      │
0──────0 ← NODES                └──────┘  └──────┘
│  -   │                        Band 1    Band 2
└──────┘
Gap vanishes at nodes           Full gaps on BOTH bands
QP at ALL temperatures          QP exponentially suppressed
```

### Strategic Recommendation

**Focus R&D on FeSe/STO monolayer junction development:**

If successful, could enable:
- Quantum computing at 4K (liquid helium, not dilution fridge)
- Error rates < 0.1% (10× below QEC threshold)
- Dramatic reduction in QC infrastructure cost

## Model Limitations

1. **T1/T2 predictions depend on assumed junction quality** (Q = 10⁴)
2. **Multi-band physics not experimentally validated for qubits**
3. **Interface effects may have additional loss mechanisms**
4. **No actual pnictide qubit exists for validation**

These limitations highlight the need for experimental development.

## Arc Status

| Session | Topic | Status |
|---------|-------|--------|
| #301 | Coherence framework | ✓ Complete |
| #302 | TLS mechanisms | ✓ Complete |
| #303 | QEC thresholds | ✓ Complete |
| #304 | Cuprate feasibility | ✓ Complete |
| #305 | Pnictide potential | ✓ Complete |
| #306 | FeSe/STO deep dive? | Planned |

## Next Steps

1. **Session #306**: Detailed FeSe/STO junction physics
2. Interface engineering principles for η reduction
3. Comparison with other emerging materials (nickelates, hydrides)
4. Experimental validation pathway definition

## Files Created

- `simulations/session305_pnictide_qubit_potential.py` - Full analysis
- `simulations/session305_pnictide_qubit_potential.png` - Visualization

## Conclusion

**Iron pnictides offer a promising path to high-temperature qubits:**

1. **S± pairing eliminates gap nodes** → exponential QP suppression
2. **Lower η** → fewer TLS → longer coherence
3. **FeSe/STO is the standout** → η=0.08, Tc=65K, no arsenic
4. **Junction technology is the bottleneck** → 5-10 year R&D needed

The η framework provides clear material guidance: optimize for low η AND nodeless gaps. FeSe/STO uniquely combines both advantages with practical benefits.

---

*"The same physics that governs interface-enhanced superconductivity may enable the next generation of quantum computers."*
