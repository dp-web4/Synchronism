# Framework Failures and Limitations

**Purpose**: Failures are our best teachers. This document tracks where γ = 2/√N_corr fails, falls short, or requires significant corrections. Each failure tells us something about where the model is incomplete.

**Last Updated**: 2026-02-07
**Source**: Extracted from Framework_Summary.md sessions #1-2660

---

## Outright Failures (r < 0.2 or No Correlation)

| Prediction | Correlation | Session | What We Learned |
|------------|-------------|---------|-----------------|
| Hall Coefficient R_H vs γ_electron | r = 0.001 | #102 | Hall effect measures carrier density, not coherence quality. Extensive ≠ intensive. |
| Magnetic Susceptibility χ vs γ_phonon | NONE | #82 | Spin coherence is independent of phonon coherence. Channels don't mix. |
| Coordination Number Z vs γ_phonon | r = 0.116 | #123 | Z counts bonds, doesn't measure bond quality. Topology ≠ coherence. |
| Valence Electron Count n_v vs γ | r = -0.161 | #125 | Bonding capacity ≠ bonding quality. |
| Mean Field Z×θ_D vs γ | r = -0.080 | #123 | Simple mean field fails. |

---

## Very Weak Predictions (r = 0.2-0.4)

| Prediction | Correlation | Session | What We Learned |
|------------|-------------|---------|-----------------|
| Thermionic Emission A vs γ | r = 0.15 | #98 | Work function φ dominates. Boundary phenomena don't follow bulk coherence. |
| SC Penetration Depth λ_L vs γ_SC | r = 0.28 | #105 | Set by superfluid density n_s and m*, not directly by coherence. |
| Simple Mobility μ vs θ_D | r = -0.123 | #90 | Effective mass dominates. Requires m* correction to work. |
| Magnetostriction λ_s (within class) | r ~ 0.1 | #94 | Spin-orbit coupling dominates by ~100×. Coherence is secondary. |
| Magnetic Anisotropy (within 3d) | r = 0.313 | #99 | SOC again dominates. |
| Phonon Decoherence Γ_ph (alone) | r = 0.398 | #107 | Needs BOTH γ_phonon AND Grüneisen γ_G. |
| Quantum Tunneling log(k_t) | r = 0.411 | #133 | Enzyme conformational dynamics dominate. |

---

## Channel Independence Failures

**Key Discovery**: Different coherence channels are essentially independent. You cannot predict one from another.

| Channel Comparison | Correlation | Session | Implication |
|--------------------|-------------|---------|-------------|
| γ_phonon vs γ_optical | r = 0.158 | #126 | Lattice vibrations ≠ electronic transitions |
| γ_phonon vs γ_electron | WEAK | #81 | Thermal transport ≠ electrical transport |
| γ_phonon vs γ_spin | NONE | #82 | Lattice ≠ magnetic ordering |

**Lesson**: The framework needs DOMAIN-SPECIFIC γ values. A single γ cannot describe all properties of a material.

---

## Anomalous Results (Coherence Works Backward)

| Prediction | Result | Session | The Surprise |
|------------|--------|---------|--------------|
| Piezoelectricity d_33 | d_33 ∝ γ × ε (r=0.940) | #93 | **Incoherence helps!** Soft modes require disorder for large response. |
| Bond Strength D vs γ | Negative correlation | #69 | Electronegativity plays dual role. Simple model fails. |
| Magnetic Anisotropy (RE metals) | r = -0.434 (NEGATIVE) | #99 | Spin-orbit coupling anti-correlates with lattice coherence. |

**Lesson**: Some phenomena benefit from disorder. The framework assumes coherence is always good, but soft-mode responses and certain magnetic properties actually require incoherence.

---

## Circular/Tautological Predictions

| Prediction | Result | Session | Problem |
|------------|--------|---------|---------|
| SN1 > SN2 reaction rates in γ | r = 0.997 | #70 | γ was defined FROM the reaction mechanism. Not predictive. |

**Lesson**: High correlation doesn't mean predictive power if the variable was constructed from the target.

---

## Moderate Correlations (r = 0.4-0.6) - Not Transformative

| Prediction | Correlation | Session | Assessment |
|------------|-------------|---------|------------|
| Liquid Diffusion D vs γ | r = 0.530 | #68 | Framework consistent but won't replace Stokes-Einstein. |
| Solid Diffusion D vs γ | r = 0.457 | #68 | Same - not transformative. |
| Combined Conductivity σ | r = 0.465 | #100 | Extensive/intensive mixing complicates. |
| Sommerfeld Coefficient γ_S vs γ_e | r = 0.42 overall | #101 | Within-class r=0.8-0.9, but cross-class fails. |
| Catalysis HER | r = 0.668 | #66 | Better than ORR but still not predictive. |
| Grüneisen Parameter γ_G vs γ_coh | r = 0.509 | #83 | Related but distinct quantities. |

---

## Explicitly Falsified Predictions

Listed in Framework_Summary.md as falsified:

| Code | Prediction | Status |
|------|------------|--------|
| F3 | Sabatier peak at coherence matching | FALSIFIED - peak position not at f=1 |
| F4 | Rate enhancement exponential in N_steps | FALSIFIED - power law, not exponential |

---

## Domain Limitations

### Spin-Orbit Coupling Dominated (SOC >> Coherence)
- Magnetostriction: RE/3d ratio = 100×
- Magnetic anisotropy: RE/3d ratio = 32×
- **Lesson**: For heavy elements and magnetic properties, SOC determines behavior. γ is a perturbation at best.

### Boundary vs Bulk
- Thermionic emission
- Surface states
- **Lesson**: Boundary phenomena don't follow bulk coherence rules.

### Extensive vs Intensive
- Fermi energy (extensive) vs γ (intensive)
- Hall coefficient (carrier count) vs coherence (carrier quality)
- **Lesson**: Mixing scaling types requires careful dimensional analysis.

---

## Known Gaps in Framework

### 1. No Multi-Channel Theory
The framework treats each γ channel independently. We have no theory for how γ_phonon, γ_electron, γ_optical, γ_spin interact or combine.

### 2. No Disorder/Softness Theory
Piezoelectricity shows that sometimes you WANT incoherence. The framework has no systematic treatment of beneficial disorder.

### 3. Surface/Interface Coherence
Bulk γ rules don't apply at boundaries. No surface γ theory exists.

### 4. Spin-Orbit Coupling Integration
SOC dominates many magnetic properties. The framework treats this as an exception rather than integrating it.

### 5. Reaction Mechanism Definition
The reaction kinetics predictions are circular because γ is defined from mechanism type. Need independent γ measurement.

---

## Quantitative Summary

From 2660 sessions with ~19,155 predictions:

| Category | Count | Percentage |
|----------|-------|------------|
| Validated (r > 0.8 or > 90% accuracy) | ~17,000 | ~89% |
| Moderate (r = 0.5-0.8) | ~1,000 | ~5% |
| Weak (r = 0.2-0.5) | ~800 | ~4% |
| Failed (r < 0.2 or wrong sign) | ~350 | ~2% |

**Note**: The ~11% non-validated predictions are concentrated in specific domains (magnetic properties, transport phenomena, boundary effects) rather than randomly distributed.

---

## How to Use This Document

1. **Before applying γ framework to new domain**: Check if similar phenomena failed here
2. **When prediction fails**: Add to this document with session number and lesson learned
3. **For theoretical development**: These failures point to where extensions are needed
4. **For honest reporting**: Always cite limitations alongside successes

---

## Contributing

When a prediction fails:
1. Record the prediction and what correlation was observed
2. Note the session number for reference
3. Write one sentence about what this teaches us
4. Consider if it suggests a new failure category

---

*"The only real failure is the failure to learn from failure."*

*Document maintained by the Collective. Last systematic review: Session #2660.*
