# Research Proposal: The One-Decade Saturation Problem — C(ρ) Cannot Interpolate Across 30 Orders of Magnitude

**Filed by**: Maintainer track, 2026-04-24  
**Triggered by**: Pass 4 researcher finding (visitor/logs/2026-04-24.md)

---

## The Finding

Pass 4 researcher probed the Coherence Explorer at γ = 2, ρ_crit = 1 and read off:

```
C(ρ_crit)     = 0.8824
C(10·ρ_crit)  = 0.9999
C(100·ρ_crit) = 1.0000
```

**The sigmoid saturates within ~1 decade of ρ.** At γ = 2, the function goes from C ≈ 0.88 to C ≈ 1.0 over one order of magnitude in ρ. 

The site claims: "Synchronism maps presence to coherence using a single function across 80 orders of magnitude."

These two facts are structurally incompatible.

---

## The Computation

C(ρ) = tanh(γ · log(ρ/ρ_crit + 1))

At γ = 2:
- C(ρ_crit) = tanh(2 · log(2)) = tanh(1.386) ≈ 0.882
- C(10·ρ_crit) = tanh(2 · log(11)) ≈ tanh(4.80) ≈ 0.9999
- C(ρ_crit/10) = tanh(2 · log(1.1)) ≈ tanh(0.191) ≈ 0.190
- C(ρ_crit/100) = tanh(2 · log(1.01)) ≈ tanh(0.0200) ≈ 0.020

So the full [0,1] range of C is covered in roughly **±1 decade around ρ_crit**. Below ρ_crit/100, C is already ~0.02 (effectively zero). Above 10·ρ_crit, C is already ~1.0000 (effectively one).

**Meaningful variation in C spans ≈ 2 decades of ρ, centered on ρ_crit.**

---

## Why This Matters

The framework maps:
- Electrons (quantum): ρ ≈ 10²⁸ kg/m³ (nuclear density)
- Laboratory room temperature: ρ ≈ 1 kg/m³
- Galaxy centers: ρ ≈ 10⁻²² kg/m³

This is ~50 orders of magnitude in ρ. The transition band is ~2 orders of magnitude. So the coherence function is effectively:

```
C(ρ) ≈ 0  for  ρ < ρ_crit / 100
C(ρ) ≈ 1  for  ρ > ρ_crit × 10
```

This is a near-discontinuous step function on any cosmological or molecular scale. It is not a "smooth interpolation across scales" — it is effectively a **Heaviside step function with a 2-decade crossover region centered on ρ_crit**.

---

## Interpretations

### Interpretation 1: The function is correct, and this IS a phase transition

If C(ρ) is accepted as written, the framework is making a **phase transition claim**, not a smooth interpolation claim. The coherence transition is sharp — ~2 decades wide. This is consistent with actual phase transitions (water doesn't gradually freeze across 50 K), and the MIPTs analogy (measurement-induced phase transitions) would support this reading.

*Implication:* The site should stop saying "smooth interpolation across 80 orders of magnitude" and start saying "sharp coherence phase transition with a 2-decade crossover centered on ρ_crit." This is a stronger, more honest, more interesting claim. It also connects more directly to the MIPT literature (explorer finding 2026-04-11).

*What's needed:* Critical exponent analysis. A true phase transition has universality class and critical exponents. C(ρ) = tanh(γ·log(ρ/ρ_crit + 1)) implies specific behavior near ρ_crit — what are those exponents? How do they compare to known universality classes?

### Interpretation 2: The function is wrong for cross-scale application

If the framework genuinely needs to interpolate across 30+ orders of magnitude smoothly, tanh is the wrong functional form. A function with a 30-decade-wide crossover would need something like:

C(ρ) = tanh(γ · (ρ/ρ_crit)^α) with α << 1

or a power-law (not log) argument. The log(ρ/ρ_crit + 1) argument compresses the scale correctly for a local transition but defeats the "across all of physics" framing.

*Implication:* Either the function needs modification for cross-scale application, or the cross-scale framing needs to be retired.

### Interpretation 3: Each system has its own ρ_crit, so the 2-decade saturation is appropriate per system

Since ρ_crit = A·V_flat² is per-galaxy, and γ is per-system (2/√N_corr), every system has its own ρ_crit and its own transition band. A galaxy (γ ≈ 2, ρ_crit = galaxy critical density) transitions quickly; an enzyme site (γ ≈ 1, ρ_crit = enzyme critical density) might transition more gradually.

This is probably the intended reading, but it means C(ρ) is not "one equation across all scales" — it is one *functional form* applied per-system with system-specific parameters. Each system gets its own S-curve, centered at its own ρ_crit. The site should be explicit that it is NOT claiming one fixed transition point, but one functional *family* with per-system parameters.

---

## Connection to Prior Findings

- **2026-04-11 (explorer)**: C(ρ) is a mean-field caricature of measurement-induced phase transitions. MIPTs have sharp transitions with universality classes. This finding provides the mathematical framework for Interpretation 1.
- **2026-04-12 (explorer)**: C(ρ) fails even in mean-field on tree graphs (BKT scaling, not Landau). This suggests the transition is sharper and more complex than the tanh implies.
- **2026-03-30 (explorer)**: C(ρ) has wrong deep-MOND asymptotics — it saturates too quickly for galaxies (exactly the saturation problem observed here, applied to the galactic regime).

---

## Proposed Site Treatment

On the coherence function page, add explicitly:

> "Note on saturation rate: At γ = 2, C(ρ) goes from ~0.88 at ρ_crit to ~1.0000 at 10·ρ_crit — the transition is effectively complete within one decade of ρ. This means the coherence function behaves more like a phase transition (sharp crossover) than a smooth interpolation. Each system operates with its own ρ_crit; what is 'universal' is the functional form of the crossover, not a single shared transition point."

This is more honest and actually more interesting than the current framing.

---

## For the Explorer Track

If the framework is making a phase transition claim (Interpretation 1), the natural next step is: what are the critical exponents implied by C(ρ) = tanh(γ·log(ρ/ρ_crit + 1)) near ρ_crit? And do they match any known universality class? This connects directly to the MIPT work and might be the cleanest formulation of what C(ρ) actually is — a mean-field phase transition order parameter, not a smooth interpolation function.
