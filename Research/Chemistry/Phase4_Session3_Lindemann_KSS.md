# Phase 4 Session 3: Lindemann-KSS — η/s at the Melting Point

*Date: 2026-04-08 | Chemistry Track Phase 4: Quantum Critical*

---

## Motivation

Phase 4 Sessions 1-2 established:
- η/s orders material systems from quantum-critical (near KSS) to classical across 7 orders of magnitude
- Conceptual equivalence between entity criterion (γ/f < 1) and KSS bound (A = 1/4π) — but quantitative mapping fails (4π gap)
- The next suggested direction: **Lindemann parameter as entity criterion**

The crystal→liquid transition is the most universal phase transition in condensed matter. The Lindemann criterion (⟨u²⟩^(1/2)/d ≈ 0.1) is remarkably universal across material classes. If the Synchronism entity criterion is physically meaningful, the melting transition — where an ordered crystal (entity) becomes a disordered liquid (process) — should map onto the KSS framework.

**Test**: Does η/s at the melting point provide a material-independent entity criterion? Is there a universal A/A_KSS value at which crystals dissolve?

**Non-circularity**: η/s at T_m uses measured liquid viscosity and thermodynamic entropy — NOT θ_D directly.

---

## Dataset

26 materials spanning 6 classes: alkali metals (5), noble metals (3), simple metals (7), transition metals (5), refractory metals (3), semiconductors (2), rare earth (1).

| Material | T_m (K) | η_liq (mPa·s) | θ_D (K) | A/A_KSS | L_est |
|----------|---------|----------------|---------|---------|-------|
| Li | 454 | 0.57 | 344 | 132 | 0.100 |
| Na | 371 | 0.68 | 158 | 235 | 0.089 |
| K | 336 | 0.51 | 91 | 298 | 0.091 |
| Rb | 312 | 0.67 | 56 | 441 | 0.090 |
| Cs | 302 | 0.68 | 38 | 517 | 0.097 |
| Cu | 1358 | 4.0 | 343 | 380 | 0.069 |
| Ag | 1235 | 3.88 | 225 | 495 | 0.068 |
| Au | 1337 | 5.13 | 165 | 595 | 0.072 |
| Al | 934 | 1.3 | 428 | 202 | 0.062 |
| Pb | 601 | 2.65 | 105 | 576 | 0.061 |
| Sn | 505 | 1.85 | 200 | 405 | 0.041 |
| Zn | 693 | 3.5 | 327 | 485 | 0.047 |
| In | 430 | 1.75 | 108 | 347 | 0.072 |
| Bi | 545 | 1.67 | 119 | 373 | 0.050 |
| Fe | 1811 | 6.9 | 470 | 674 | 0.062 |
| Ni | 1728 | 5.5 | 450 | 497 | 0.063 |
| Co | 1768 | 5.0 | 445 | 460 | 0.064 |
| Ti | 1941 | 5.2 | 420 | 698 | 0.068 |
| Cr | 2180 | 6.7 | 630 | 682 | 0.052 |
| W | 3695 | 8.0 | 400 | 848 | 0.052 |
| Mo | 2896 | 5.5 | 450 | 613 | 0.057 |
| Ta | 3290 | 7.5 | 240 | 860 | 0.079 |
| Gd | 1585 | 4.8 | 200 | 1064 | 0.058 |
| Ge | 1211 | 0.74 | 374 | 103 | 0.047 |
| Si | 1687 | 0.58 | 645 | 72 | 0.055 |
| Ga | 303 | 2.04 | 320 | 405 | 0.030 |

Sources: Iida & Guthrie (1988), Assael et al. (2006), CRC Handbook, JANAF tables.

---

## Results

### 1. η/s at Melting Is NOT Universal

| Statistic | A/A_KSS | Lindemann L |
|-----------|---------|-------------|
| Mean | 479 | 0.065 |
| Median | 473 | 0.063 |
| Std | 235 | 0.017 |
| CV | **0.49** | **0.26** |
| Spread (max/min) | **14.7×** | 3.4× |

The Lindemann parameter clusters nearly 2× tighter than η/s (CV 0.26 vs 0.49). Melting is better described as an amplitude criterion than a dissipation criterion.

**In log space**, A/A_KSS is more uniform: CV(log₁₀) = 0.106. Melting occurs at A/A_KSS ≈ 100-1000 consistently, placing it ~2-3 orders of magnitude above the quantum critical regime and at the lower end of the "classical" range.

### 2. Lindemann and η/s Are Orthogonal

```
r(log(A/A_KSS), log(L)) = -0.048
```

These two quantities are completely uncorrelated. They measure different physics:
- **Lindemann**: displacement amplitude / lattice spacing (structural integrity)
- **η/s**: dissipation rate / entropy density (coherence vs thermalization)

Melting is determined by Lindemann (amplitude), not by η/s (dissipation). Materials can melt at widely different positions on the KSS spectrum.

### 3. ALL Crystals Fail the Entity Criterion at Melting

Using the Klemens model for phonon lifetime (γ_G = 1.7):

| Class | Q_melt range | γ/f range |
|-------|-------------|-----------|
| Alkali | 0.005-0.20 | 5-182 |
| Noble | 0.005-0.02 | 45-190 |
| Transition | 0.02-0.03 | 35-62 |
| Refractory | 0.002-0.008 | 120-543 |
| Semiconductor | 0.03-0.05 | 20-30 |

**Every crystal at its melting point has Q << 1 (γ/f >> 1).** By the entity criterion (γ/f < 1), every crystal is a "process," not an "entity," at its melting point.

More importantly: crystals have Q << 1 at temperatures far below melting. The phonon entity criterion is NEVER satisfied for crystals at any physically relevant temperature.

### 4. What Determines η/s at Melting?

Regression: log(A/A_KSS) = 1.84 + 0.24·log(T_m) + 0.52·log(M) − 0.92·log(ρ)

- R² = 0.78 without θ_D
- R² = 0.84 with θ_D (adds 0.06 — moderate, not dominant)
- **Dominant factor**: liquid density ρ (coefficient −0.92)
- η/s at melting is mostly determined by how dense and heavy the atoms are

η/s correlates with liquid viscosity (r = 0.77) more than entropy density (which varies less, CV 0.44 vs 0.72 for η).

---

## Interpretation

### The Two-Criterion Structure of Melting

The melting transition involves two independent criteria:

1. **Lindemann (structural)**: ⟨u²⟩^(1/2)/d ≈ 0.06-0.10 — UNIVERSAL (CV 0.26)
   - Measures: displacement amplitude vs lattice spacing
   - Physics: atoms have moved "too far" from equilibrium positions
   - Entity type: STRUCTURAL (order vs disorder)

2. **KSS position (dissipative)**: A/A_KSS ≈ 100-1000 — NOT universal (CV 0.49)
   - Measures: dissipation rate vs entropy density
   - Physics: where on the quantum-classical spectrum
   - Entity type: MODAL (coherent oscillation vs thermalized fluctuation)

These are orthogonal (r = -0.05). Melting is determined by criterion 1, not criterion 2.

### Challenge to the Entity Criterion

The Synchronism entity criterion (γ/f < 1 → entity; γ/f > 1 → process) is a **modal/dissipative** criterion. It asks: does this oscillation survive longer than one period?

But crystals are NOT modal entities. They are **structural entities** — they exist because of crystallographic order, not because individual phonon modes persist. A crystal at room temperature has phonons that are heavily damped (Q << 1), yet the crystal persists because the TIME-AVERAGED atomic positions maintain long-range order.

The distinction:
- **Wave entity** (quantum particle): pattern recurs because oscillation amplitude > damping → γ/f < 1 ✓
- **Structural entity** (crystal): pattern persists because mean displacement < lattice spacing → L < L_crit ✓, even though γ/f >> 1 for every phonon mode

This means the entity criterion as formulated applies to **wave-like entities** (quantum particles, vortices, resonant modes) but NOT to **structural entities** (crystals, molecules, biological structures).

### Connecting to the Four-Regime Classification (Phase 2)

This maps onto the Phase 2 four-regime classification:

| Regime | Property type | Entity criterion |
|--------|--------------|-----------------|
| 0: Neutral | Counting/structural | **Lindemann (amplitude)** |
| 1: Coherence | Propagation | γ/f (damping/frequency) |
| 2: Incoherence | Response | Inverse γ/f |
| 3: Barrier | Activated | exp(-E/kT) (Boltzmann) |

The melting transition is a **Regime 0 → Regime 0 transition** (ordered structure → disordered structure). The entity criterion γ/f < 1 operates in Regime 1. This is why it doesn't apply to melting.

### Where on the KSS Hierarchy?

Updated hierarchy with melting:

| System | A/A_KSS | Character |
|--------|---------|-----------|
| QGP | 1-10 | Quantum critical |
| Cold Fermi gas | ~5 | Quantum critical |
| He-4 at λ-point | ~12 | Near quantum |
| **Melting (semiconductors)** | **72-103** | **Classical onset** |
| **Melting (metals)** | **130-1064** | **Classical** |
| Room T liquids | 400-4000 | Classical |
| Gases | 1000-10000 | Deeply classical |

Melting occurs in the "classical onset to classical" zone — well above quantum criticality, but at the lower end of the classical range.

Silicon and Germanium melt at the lowest A/A_KSS (72, 103) because their covalent-to-metallic liquid transition produces anomalously LOW viscosity. This is a genuine outlier: these materials melt from a covalent solid to a metallic liquid, which is a qualitatively different transition than metal→liquid metal.

---

## Status of Phase 4 Findings

| Finding | Status | Sessions |
|---------|--------|----------|
| KSS orders quantum→classical | ✅ CONFIRMED | #1 |
| Entity ↔ KSS conceptual equivalence | ✅ VALID (not quantitative) | #1-2 |
| 4π gap (cavity vs black hole) | ❌ UNRESOLVED | #2 |
| η/s universal at melting | ❌ NOT universal (CV 0.49) | #3 |
| Lindemann more universal than η/s | ✅ CONFIRMED (CV 0.26 vs 0.49) | #3 |
| **Crystals never satisfy γ/f < 1** | ✅ **NEW FINDING** | #3 |
| **Entity criterion ≠ structural entity** | ✅ **NEW FINDING** | #3 |
| L and η/s orthogonal | ✅ CONFIRMED (r = -0.05) | #3 |

---

## Key Finding: Two Types of Entity

**The most important result of this session:**

The Synchronism framework implicitly assumes one kind of entity — a recurring oscillation pattern. The entity criterion γ/f < 1 applies to this. But condensed matter has (at least) two kinds:

1. **Oscillatory entities**: Quantum particles, standing waves, vortex rings — pattern recurs through dynamic oscillation. Entity criterion: γ/f < 1 (damping vs frequency).

2. **Structural entities**: Crystals, molecules, biological structures — pattern persists through time-averaged positional order. Entity criterion: L < L_crit (displacement vs spacing).

These are DIFFERENT criteria measuring DIFFERENT physics. The oscillatory entity criterion is about whether a wave survives one cycle. The structural entity criterion is about whether a configuration maintains order despite thermal fluctuations.

**Implication**: The FUNDAMENTALS.md definition — "Entity = recurring pattern of Intent distribution over tick sequences" — accommodates both. An oscillatory entity recurs through dynamic cycling. A structural entity recurs through statistical persistence of time-averaged positions. Both are "recurring patterns" — but with different survival mechanisms and different failure modes.

The entity criterion γ/f < 1 is ONE test (for oscillatory entities). There should be an analogous test for structural entities — and the Lindemann criterion IS that test.

---

## Open Questions for Phase 4 Session 4+

1. **Can the structural entity criterion (Lindemann) be derived from the Synchronism substrate?** The oscillatory criterion (γ/f) was derived from cavity impedance mismatch (Session 18 primary track). What substrate mechanism produces the L < 0.1 threshold?

2. **Glass transition**: Glasses are structural entities that DON'T melt sharply. They lose structural order gradually. What does η/s do at the glass transition T_g? (η increases to ~10¹² Pa·s, so A/A_KSS ≈ 10⁸-10¹⁰ — VERY far from KSS.)

3. **Superconducting transition**: Tc is where a new KIND of entity (Cooper pair condensate) forms. Is this an oscillatory or structural entity? Does it satisfy γ/f < 1?

---

*Phase 4 Session #3 — Chemistry Track*
*Finding: Melting is an amplitude criterion (Lindemann), not a dissipation criterion (η/s). Crystals are structural entities that never satisfy the oscillatory entity criterion γ/f < 1.*
*Assessment: Two types of entity identified. Framework needs structural entity criterion alongside oscillatory entity criterion.*
