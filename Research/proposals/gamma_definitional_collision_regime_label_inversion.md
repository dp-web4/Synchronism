# Proposal: γ Definitional Collision and Regime Label Inversion

**Filed:** 2026-05-04
**Source:** Synchronism site visitor feedback (Pass 3: Grad Student, Pass 4: Researcher)
**Status:** Open — needs resolution before γ-calculator and coherence-explorer can be trusted

---

## The Problem

The coherence parameter γ appears in two incompatible forms across the site:

**Form A (universal constant, 6D phase space derivation):**
> "γ derived as 2.0 from 6D phase space considerations"
> — /coherence-function, cited as a derived result

**Form B (system-dependent function, CLT derivation):**
> γ = 2/√Ncorr
> — /gamma-calculator, /parameter-derivations, cited as the operational formula

These can be formally reconciled (Form A corresponds to Form B at Ncorr = 1, the single-particle reference point), but the reconciliation is never stated. The site presents both as if they are the same claim, which they are not:

- Form A is a universal constant that never changes.
- Form B is a function of the number of correlated degrees of freedom — a system property that varies by orders of magnitude across physical systems.

If Form B is the operational formula, then Form A is just a normalization condition (γ=2 ↔ Ncorr=1), not a first-principles derivation of the coherence parameter's value.

---

## The Downstream Consequence: Regime Label Inversion

The γ-calculator and phase-boundary-visualizer assign physical systems to "quantum" and "classical" regimes using Form B:

| System | Ncorr (estimated) | γ = 2/√Ncorr | Site Regime |
|--------|-------------------|--------------|-------------|
| Ideal gas | 1 | 2.000 | Quantum (γ > 1.4) |
| BCS superconductor (Al) | 10,000 | 0.020 | Classical (γ < 0.6) |
| BEC | 10⁶ | 0.002 | Classical (γ < 0.6) |
| Neutron star | 10¹⁰ | 0.0002 | Classical (γ < 0.6) |

Under standard physics usage, BEC and BCS superconductors are the canonical examples of **macroscopic quantum coherence** — they are defined by long-range quantum order. An ideal gas is a textbook classical system. The framework's regime labeling is inverted relative to standard terminology.

This is not a labeling convention problem. It reflects an ambiguity in what γ is measuring:

- If γ measures **phase-space concentration** (high γ = single particle, few degrees of freedom), then "quantum" means "single-particle quantum" and "classical" means "many-body collective." That is a non-standard but internally consistent use.
- If γ measures **coherence intensity** (high γ = more coherent), then the BEC/BCS inversion IS an error.

The site uses "quantum/classical" language, which implies the second reading, while the math implements the first.

---

## The Three-Case Analysis

**Case 1: γ = 2/√Ncorr with Ncorr as dynamical correlation count**
- "Quantum regime" means "single-particle" (small Ncorr).
- "Classical regime" means "strongly correlated collective system" (large Ncorr).
- BEC at Ncorr=10⁶ correctly falls in the collective basin.
- **Problem:** Labels must be renamed. "Quantum" and "Classical" mislead. Better: "Single-particle" and "Collective" or "Correlated."

**Case 2: γ = 2/√Ncorr with Ncorr as quantum coherence volume overlap**
- BEC has Ncorr_quantum = 1 (all particles in same mode), so γ_BEC = 2, placing it in the "quantum" regime.
- This is the correct quantum/classical framing, but requires Ncorr to mean wavefunction-overlap count, not dynamical correlation count.
- **Problem:** The page already acknowledges that Ncorr values in presets are "approximate estimates, not measured physical pair counts." If Ncorr is quantum overlap, the BCS preset of Ncorr=10,000 is ~6 orders of magnitude below the physical Cooper pair count in a coherence volume, AND it would not match the dynamical-correlation Ncorr used for galaxy rotation fits.

**Case 3: γ is definitionally 2.0 for all galaxy rotation work, and γ = 2/√Ncorr is a separate empirical calibration**
- The two uses of γ are in different domains: galaxy rotation uses γ=2 (Ncorr=1); condensed matter uses the rescaled form.
- **Problem:** This requires explicit scope statement. If γ=2 is fixed for galaxy work, it means the framework treats galaxy baryons as uncorrelated (Ncorr=1) — an odd assumption for a many-body system with 10¹¹ stars.

---

## The Research Question

**Which Ncorr interpretation is load-bearing?**

- For galaxy rotation fits, γ=2 → Ncorr=1. Is that a derivation or a calibration?
- For condensed matter presets, Ncorr is "estimated," not computed from first principles.
- The framework needs a single, operational definition of Ncorr that (a) gives the right γ in the galactic regime and (b) correctly places BEC/BCS in their standard quantum category.

A candidate resolution: adopt **Case 1 + label renaming**. Replace "Quantum regime" with "Uncorrelated / Single-particle regime" and "Classical regime" with "Collective / Strongly correlated regime." This is consistent with the math (high γ = fewer correlated degrees of freedom = more single-particle-like behavior). BEC then correctly appears as a "Collective" system (large Ncorr, small γ) — which is accurate: BECs are maximally correlated at the collective level, even if each particle is in the same quantum state.

The label "Classical" for BEC is the surface problem. The deeper question is whether the framework can supply a first-principles Ncorr for any given physical system, or whether Ncorr is always a fitted parameter.

---

## Proposed Maintainer Actions

1. **Short-term (site fix):** Rename "Quantum regime" → "Single-particle / Uncorrelated" and "Classical regime" → "Collective / Strongly correlated" in coherence-explorer and phase-boundary-visualizer. Add a caveat note explaining the non-standard usage.

2. **Medium-term (research):** Specify Ncorr operationally. The current presets are acknowledged estimates. What measurement or derivation would give Ncorr for a real system? Options:
   - Ncorr = (correlation length / particle spacing)^d
   - Ncorr = number of particles within one coherence length
   - Ncorr = thermodynamic susceptibility proxy

3. **Long-term (framework question):** If γ=2 for galaxies requires Ncorr=1 (uncorrelated baryons), is that physically motivated or is it a fit? A derivation of why Ncorr_galaxy = 1 from the theory's own principles would strengthen the framework considerably.

---

## Related Proposals
- `lorentz_invariance_gap_kinematic_layer.md` — no dynamical equation for the phase field
- `dual_C_symbol_ambiguity_and_bridge_derivation.md` — two-C symbol collision
