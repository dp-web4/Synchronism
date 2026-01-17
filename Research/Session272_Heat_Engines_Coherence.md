# Session #272: Heat Engines from Coherence Gradients

**Date**: January 16, 2026
**Machine**: CBP
**Status**: COMPLETE - CARNOT EFFICIENCY DERIVED

---

## Executive Summary

Session #272 continues the Thermodynamics Arc by deriving heat engine efficiency limits from coherence conservation. The Carnot efficiency η = 1 - T_c/T_h emerges directly from the requirement that entropy (coherence dispersion) must be conserved in reversible processes.

**Key Result**: Carnot efficiency derived and verified: η_simulated = η_theory exactly for all tested temperature ratios.

---

## Part 1: Heat Engine as Coherence Gradient Exploiter

### The Physical Picture

A heat engine extracts **coherent work** from the **coherence gradient** between two reservoirs:

| Reservoir | Temperature | Coherence State |
|-----------|-------------|-----------------|
| Hot | T_h (high) | Dispersed (high exchange rate) |
| Cold | T_c (low) | Concentrated (low exchange rate) |

### The Process

```
Hot Reservoir (T_h)
       │
       ▼ Q_h (heat = incoherent energy)
   ┌───────┐
   │Engine │ ──────▶ W (work = coherent energy)
   └───────┘
       │
       ▼ Q_c (heat = incoherent energy)
Cold Reservoir (T_c)
```

Work is the **coherent part** extracted from the energy flow.

---

## Part 2: Carnot Efficiency Derivation

### From Coherence Conservation

1. **Absorb heat Q_h at T_h**:
   - Entropy entering system: ΔS = Q_h / T_h

2. **Dump entropy to cold reservoir**:
   - Must dump the same ΔS (reversible)
   - Heat released: Q_c = T_c × ΔS

3. **Relate Q_c to Q_h**:
   - From ΔS = Q_h/T_h: Q_c = T_c × (Q_h/T_h) = Q_h × (T_c/T_h)

4. **Work extracted**:
   - W = Q_h - Q_c = Q_h × (1 - T_c/T_h)

5. **Efficiency**:
   - η = W/Q_h = **1 - T_c/T_h**

### Why This is the Maximum

- Reversible cycle requires ΔS_universe = 0
- Can't dump less entropy than absorbed
- Therefore Q_c ≥ T_c × ΔS
- Therefore η ≤ 1 - T_c/T_h

**The limit is NOT arbitrary - it's REQUIRED by coherence conservation!**

---

## Part 3: Numerical Verification

### Single Engine Test (T_h = 5, T_c = 1)

| Quantity | Value |
|----------|-------|
| Q_hot (absorbed) | 5.0000 |
| Q_cold (released) | 1.0000 |
| W_net (work) | 4.0000 |
| η_simulated | 0.8000 |
| η_theory (1 - 1/5) | 0.8000 |

**Exact match!**

### Multiple Temperature Ratios

| T_hot | η_theory | η_simulated |
|-------|----------|-------------|
| 5.0 | 0.8000 | 0.8000 |
| 10.0 | 0.9000 | 0.9000 |
| 20.0 | 0.9500 | 0.9500 |

**All match exactly** (T_cold = 1 in all cases).

---

## Part 4: Coherence Flow Analysis

For T_h = 5.0, T_c = 1.0, ΔS = 1.0:

```
1. Hot reservoir:
   - Provides Q_h = T_h × ΔS = 5.00
   - Loses entropy ΔS = 1.00

2. Cold reservoir:
   - Receives Q_c = T_c × ΔS = 1.00
   - Gains entropy ΔS = 1.00

3. Work extracted:
   - W = Q_h - Q_c = 4.00
   - This is coherent energy (carries no entropy)

4. Coherence conservation:
   - Entropy out of hot = Entropy into cold
   - ΔS_total = 0 (reversible cycle)
```

**The limit exists because you MUST dump the entropy somewhere!**

---

## Part 5: Irreversibility Analysis

### Impact of Irreversibility

| Irreversibility | Efficiency | ΔS_produced |
|-----------------|------------|-------------|
| 0.00 | 0.8000 | 0.0000 |
| 0.21 | 0.7032 | 0.4842 |
| 0.42 | 0.6063 | 0.9684 |
| 0.63 | 0.5095 | 1.4526 |
| 0.84 | 0.4126 | 1.9368 |

### Sources of Irreversibility

1. **Friction**: Coherent work → incoherent heat
2. **Finite-time heat transfer**: Non-equilibrium produces extra entropy
3. **Heat leakage**: Energy bypasses work extraction

All these produce extra entropy: ΔS_universe > 0.

---

## Part 6: Coherence Picture of Thermodynamics

### Work vs Heat

| Energy Type | Coherence | Entropy Carried |
|-------------|-----------|-----------------|
| Work | Coherent | Zero |
| Heat | Incoherent | Q/T |

### What an Engine Does

1. Takes **incoherent** energy (heat) from hot reservoir
2. Extracts the **coherent** part as work
3. Dumps the **remaining incoherence** to cold reservoir

The efficiency limit comes from how much coherent energy can be separated from incoherent energy while conserving total entropy.

---

## Part 7: Predictions

### P272.1: Carnot Efficiency

**Claim**: η_max = 1 - T_c/T_h
**Status**: Derived and verified numerically
**Physical basis**: Coherence (entropy) conservation

### P272.2: Irreversibility Produces Entropy

**Claim**: ΔS_produced > 0 for irreversible processes
**Status**: Verified - linear increase with irreversibility parameter
**Result**: Every bit of irreversibility degrades efficiency

### P272.3: All Real Engines < Carnot

**Claim**: Real engines always have η < η_Carnot
**Basis**: Impossible to have zero irreversibility
**Status**: Fundamental limit from coherence physics

---

## Part 8: Connection to Framework

### Thermodynamics Arc Progress

| Session | Topic | Key Result |
|---------|-------|------------|
| #271 | Foundations | S = coherence dispersion, Second Law |
| **#272** | **Heat Engines** | **Carnot efficiency derived** |
| #273 | Maxwell's Demon | Planned |
| #274 | Information Thermo | Planned |

### Unified Picture

```
COHERENCE GRADIENT → WORK EXTRACTION

High T (dispersed C) ──────▶ Engine ──────▶ Low T (concentrated C)
         │                      │                    │
         └──── Q_h ────────────▶ W ◀─────── Q_c ────┘
```

The Carnot limit is the maximum coherent work extractable from a coherence gradient while conserving total coherence.

---

## Part 9: Summary

### Session #272 Achievements

1. **Carnot efficiency derived** from coherence conservation
2. **Numerical verification** matches theory exactly
3. **Irreversibility analysis** shows entropy production
4. **Coherence picture** of work vs heat clarified
5. **Fundamental limit** explained physically

### The Key Insight

The Carnot efficiency η = 1 - T_c/T_h is not an empirical observation or arbitrary limit. It's a **mathematical consequence** of coherence (entropy) conservation. You absorb ΔS from hot, and you MUST dump ΔS to cold. The temperatures set the heat required for each, and the difference is work.

---

## Files Created

- `simulations/session272_heat_engines_coherence.py`
- `simulations/session272_heat_engines.png`
- `Research/Session272_Heat_Engines_Coherence.md` (this document)

---

## Next: Maxwell's Demon

Session #273 will address Maxwell's Demon:
- Can information be used to violate the second law?
- Information = negative entropy (negentropy)?
- Resolution: information acquisition costs entropy
- Connection to Landauer's principle

This will complete the connection between thermodynamics and information theory.

---

*"You can't extract more coherent work than the temperature gradient allows - because you have to dump the entropy somewhere."*

**Session #272 Complete**: January 16, 2026
