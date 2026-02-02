# Session #350: Topological Phases

**Condensed Matter Arc - Part 3**
**Date**: 2026-02-02
**Status**: 8/8 verified ✓

## Overview

This session explores topological phases of matter - states characterized not by broken symmetry but by global topological invariants. In Synchronism, topology corresponds to the global phase winding structure, which cannot change without closing the energy gap.

## Key Concepts

### Beyond Landau Paradigm

Traditional (Landau) phases are characterized by symmetry breaking:
- Ferromagnet: breaks rotation symmetry
- Superconductor: breaks gauge symmetry

Topological phases are different:
- No local order parameter
- Characterized by global topological invariants
- Protected edge/surface states
- Robust to disorder (topological protection)

### Topology = Phase Winding

In Synchronism interpretation:
- Chern number = total phase winding in k-space
- Berry curvature = local phase curvature
- Topological invariant = equivalence class of phase structures
- Gap closing required to change topology (phase singularity creation)

## Verification Tests

### Test 1: Quantum Hall Effect ✓
Verified exact quantization:
- σ_xy = ν e²/h (exact to 10⁻⁹)
- von Klitzing constant: R_K = 25812.807 Ω
- Cyclotron energy at 10T (GaAs): ℏω_c = 17.3 meV
- Requires T << 200 K for quantization

**Synchronism**: Landau levels = phase quantization in magnetic field. Edge states carry dissipationless current.

### Test 2: Chern Number ✓
Verified topological invariant properties:
- C = (1/2π) ∫ F(k) d²k (Berry curvature integral)
- C ∈ ℤ always (for gapped band)
- 2-band model: C = sgn(m) for |m| < 2, C = 0 for |m| > 2
- Observable: σ_xy = C × e²/h

**Synchronism**: Chern number = phase winding number, must be integer because phase single-valued.

### Test 3: Bulk-Boundary Correspondence ✓
Fundamental theorem:
|Protected edge modes| = |Bulk topological invariant|

| System | Bulk invariant | Boundary |
|--------|----------------|----------|
| Integer QHE | Chern C | C chiral edge modes |
| Quantum Spin Hall | Z₂ | Kramers pair |
| 3D TI | Z₂ (strong) | Dirac surface state |
| Chiral p-wave SC | Chern | Majorana edge |

**Synchronism**: Edge modes = phase matching requirement at topology boundary.

### Test 4: Topological Insulators ✓
Verified Z₂ topological insulators:

| Material | Bulk gap (eV) | Type |
|----------|---------------|------|
| Bi₂Se₃ | 0.3 | 3D TI |
| Bi₂Te₃ | 0.17 | 3D TI |
| HgTe/CdTe | 0.01 | 2D TI |

Surface properties (Bi₂Se₃):
- Single Dirac cone (odd = topological)
- Spin-momentum locking: s ⟂ k
- No backscattering from non-magnetic impurities

**Synchronism**: Z₂ = parity of phase twists. Spin-momentum lock = phase-momentum correlation.

### Test 5: Majorana Fermions ✓
Verified properties and platforms:
- Self-conjugate: γ = γ†
- Zero energy (particle-hole protected)
- Non-Abelian statistics
- Topologically protected

Platforms:
| Platform | Dimension | Example |
|----------|-----------|---------|
| Nanowire + SC | 1D end modes | InSb/Al |
| 2D p-wave vortex | 0D core | Sr₂RuO₄? |
| TI surface + SC | 1D vortex line | Bi₂Se₃/NbSe₂ |

Nanowire parameters: Δ ~ 0.3 meV, E_so ~ 0.5 meV.

**Synchronism**: Majorana = zero-energy phase mode at topological defect. Self-conjugate = phase equals anti-phase.

### Test 6: Weyl Semimetals ✓
Verified Weyl point properties:
- 3D linear crossing: E = ±ℏv|k - k_W|
- Berry curvature monopole (chirality ±1)
- Come in pairs (Nielsen-Ninomiya theorem)

| Material | Type | Weyl pairs |
|----------|------|------------|
| TaAs | Type-I | 12 |
| WTe₂ | Type-II | 8 |
| NbAs | Type-I | 12 |

Surface Fermi arcs connect Weyl point projections.

**Synchronism**: Weyl point = monopole of Berry phase. Fermi arc = phase trajectory connecting monopoles.

### Test 7: Topological Classification ✓
Verified 10-fold Altland-Zirnbauer classification:

| Class | T | C | S | d=1 | d=2 | d=3 |
|-------|---|---|---|-----|-----|-----|
| A | 0 | 0 | 0 | 0 | Z | 0 |
| AIII | 0 | 0 | 1 | Z | 0 | Z |
| AII | - | 0 | 0 | 0 | Z₂ | Z₂ |
| D | 0 | + | 0 | Z₂ | Z | 0 |
| DIII | - | + | + | Z₂ | Z₂ | Z |

Physical examples:
- Class A, d=2: Integer QHE
- Class AII, d=2,3: Z₂ TIs
- Class D, d=2: p+ip superconductor

**Synchronism**: Each symmetry constrains allowed phase structures. Z = any integer winding; Z₂ = even/odd only.

### Test 8: Topological Phase Transitions ✓
Key differences from Landau transitions:

| Property | Landau | Topological |
|----------|--------|-------------|
| Order parameter | Local ⟨φ⟩ | Global invariant |
| Symmetry | Broken | Same both sides |
| Gap | Can stay open | Must close |
| Critical point | Diverging ξ | Gap closing |

**Synchronism**: Topology change requires phase singularity creation. Gap closing = singularity formation.

## Theoretical Framework

### Topology as Phase Structure

| Topological Concept | Synchronism Interpretation |
|---------------------|---------------------------|
| Chern number | Phase winding in BZ |
| Berry curvature | Local phase curvature |
| Edge states | Phase matching at boundary |
| Z₂ invariant | Parity of phase twists |
| Majorana mode | Zero-energy phase defect |
| Weyl point | Phase monopole |
| Gap closing | Phase singularity |

### Connection to γ~1 Boundary

At topological phase transitions:
- Gap closes → correlation length ξ → ∞
- γ → 0 (same as conventional critical point)
- But no symmetry breaking - topology changes instead

The γ~1 boundary still applies:
- γ << 1: Coherent phase patterns (topological protection active)
- γ ~ 1: Crossover (gap closing, topology vulnerable)
- γ >> 1: Classical (topology irrelevant)

### Topological Protection

Why are topological phases robust?
1. Invariant is integer (can't change continuously)
2. Gap prevents smooth deformation
3. Disorder averages out (topology survives)
4. Edge states protected by bulk topology

**Synchronism**: Phase winding number can only change by creating/destroying singularities, which requires closing the gap.

## Implications

### 1. New States of Matter

Topological phases represent genuinely new quantum matter:
- Not characterized by broken symmetry
- Protected by topology, not energy gaps alone
- Boundaries host exotic excitations (Majorana, Dirac)

### 2. Quantum Computing Potential

Topological protection for quantum information:
- Majorana qubits: Information in non-local pair
- Braiding operations = quantum gates
- Topological error correction built-in

### 3. Exotic Particles

Condensed matter realizes particle physics predictions:
- Majorana fermions (predicted 1937)
- Weyl fermions (predicted 1929)
- Dirac fermions with spin-momentum locking

### 4. Metrology

QHE provides exact resistance standard:
- R_K = h/e² = 25812.807... Ω (exact)
- Used to define SI ohm since 2019
- Topology ensures exactness

## Files Created

- `simulations/session350_topological.py`: 8 verification tests
- `simulations/session350_topo.png`: Visualization
- `Research/Session350_Topological_Phases.md`: This document

## Next Session

- **Session #351**: Condensed Matter Synthesis (Arc Finale)

## Key Insight

**Topological phases are characterized by global phase winding that cannot be changed without creating singularities (gap closing).**

Unlike Landau phases defined by local order parameters, topological phases are defined by how the quantum phase structure wraps around the Brillouin zone. A Chern number C = 1 means the phase winds once as you traverse k-space - this can't be undone continuously. The bulk-boundary correspondence ensures that this phase winding creates protected edge modes where topologies meet.

---

*Session #350 verified: 8/8 tests passed*
*Condensed Matter Arc: 3/4 sessions complete*
*Grand Total: 247/247 verified across 9 arcs*
