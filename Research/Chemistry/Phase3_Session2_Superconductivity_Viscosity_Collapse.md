# Phase 3 Session 2: Superconductivity as Viscosity Collapse

*Date: 2026-04-01 | Chemistry Track Phase 3: CFD Cross-Pollination*

---

## Question

Does the N-S multi-component framing generate predictive content for superconductivity?

**Specific test**: The superconducting transition is a "viscosity collapse" where electron channel viscosity μ_eff → 0. The channel coupling parameter λ_ep is the "Reynolds number" governing whether Cooper pairs form. Can the N-S Reynolds number analog Re_ep = λ_ep/γ_phonon predict Tc better than λ_ep alone?

---

## Setup

### N-S Framing (from Session 1)
- Phonon channel: effective viscosity ν_ph ∝ γ_phonon = 2T/θ_D
- Electron channel: effective viscosity ν_e ∝ 1/σ
- λ_ep = Lorentz force coupling term between channels
- Cooper pair formation = "phase locking" between channels when coupling exceeds viscous damping

### The Reynolds Number Analog
Re_ep = λ_ep / γ_phonon(T_ref) = λ_ep × θ_D / (2T_ref)

At T_ref = 300K: Re_ep = λ_ep × θ_D / 600

**Prediction**: If the N-S framing adds content, Tc should correlate with Re_ep (not just λ_ep).

### Modified BCS from N-S Framing
BCS: Tc = 1.13 × θ_D × exp(-1/λ_ep)
N-S Modified: Tc = A × θ_D × exp(-α × γ_ph/λ_ep) where γ_ph = 600/θ_D

### Data
24 elemental superconductors with measured: Tc, θ_D, λ_ep, γ_Sommerfeld

---

## Results

### Model Comparison (24 elemental superconductors)

| Model | r | RMSE_log | Free parameters |
|-------|---|----------|-----------------|
| BCS (universal A=1.13) | 0.806 | 3.654 | 0 |
| BCS (fitted A) | 0.806 | 1.071 | 1 |
| **N-S Modified BCS** | **−0.364** | **1.838** | **2** |
| Power law (Re_ep, θ_D) | 0.868 | 0.749 | 3 |
| Allen-Dynes (fitted μ*) | **0.975** | **0.346** | 2 |

**N-S Modified BCS performs WORSE than BCS** — negative correlation, RMSE increases by 72%.

### Why the Modification Fails

BCS exponent: exp(−1/λ_ep)
N-S Modified exponent: exp(−γ_ph/λ_ep) = exp(−600/(θ_D × λ_ep))

BCS already uses θ_D as the prefactor. Replacing "1" with γ_ph = 600/θ_D introduces θ_D into the exponent as well. This is **circular** — for the high-θ_D, weak-coupling metals (W, Os, Ir), the exponent becomes tiny (−600/600/0.3 ≈ −3.3), predicting moderate Tc; but they actually have near-zero Tc. For soft-lattice, strong-coupling metals (Hg, Pb), the exponent becomes large (−600/71/1.55 ≈ −5.5), predicting low Tc; but they have the highest Tc in this dataset.

**The modification inverts the prediction.**

### Key Collapse Test

BCS framework: Tc/θ_D should depend on λ_ep only. Collapse quality:

- r(log(Tc/θ_D), −1/λ_ep) = **0.965** — excellent BCS collapse
- r(log(Tc/θ_D), −γ_ph/λ_ep) = **0.039** — no collapse at all

### The Residual Correlation (Apparent Signal)

BCS residuals DO correlate with γ_phonon: r = 0.605 (p = 0.0017).

**But this is a confound, not a signal:**
- Weak-coupling metals (λ < 0.5): r(residual, γ_ph) = 0.226 — weak
- Strong-coupling metals (λ ≥ 0.5): r(residual, γ_ph) = 0.725 — strong

The chain: γ_ph large → soft lattice (low θ_D) → high λ_ep → strong coupling → BCS underpredicts (BCS is a weak-coupling approximation).

The "correction information" in γ_ph is purely **regime identification** — soft-lattice metals are strong-coupling, where Allen-Dynes corrections apply. This is already known from standard physics.

### The Circularity Pattern

Phase 2 finding: γ = 2T/θ_D carries zero bits beyond θ_D at fixed T.
Phase 3 Session 2 finding: The same circularity applies MORE SHARPLY to superconductivity because BCS explicitly uses θ_D. Any γ-based modification of BCS is circular by construction.

---

## Class Structure

| Class | BCS mean residual | r(resid, γ_ph) | Interpretation |
|-------|-------------------|-----------------|----------------|
| sp-metals | +0.595 | 0.867 | Strong-coupling regime; BCS underpredicts |
| 4d-metals | −0.088 | 0.624 | Near BCS validity; mixed |
| 5d-metals | −1.015 | 0.286 | Heavy SOC; weak coupling; BCS slightly overpredicts |

The sp-metal class anomaly (BCS underpredicts by factor ~2-5 for Pb, Hg) is the Allen-Dynes strong-coupling correction, not a γ_phonon effect.

---

## "So What?" Assessment

**Is the N-S viscosity framing useful for superconductivity?** No, for prediction. Yes, for vocabulary.

The framing provides a consistent picture:
- Normal state: electron and phonon channels have independent viscosities; λ_ep weakly couples them
- At Tc: electron channel viscosity → 0 (zero resistance), channels "lock" (phase transition)
- Cooper pairs = frozen electron-phonon correlation (phonon field "frozen into" electron flow, MHD Rm → ∞ analog)

But this is purely vocabulary. BCS/Allen-Dynes already captures the physics with fewer assumptions.

**Is this a productive failure?** Yes.

The circularity is now sharply demonstrated for two cases:
1. Session 1: Prandtl analog works for simple metals, fails for d-metals (organizational, not predictive)
2. Session 2: N-S Modified BCS fails because γ_ph = 600/θ_D is circular with BCS (BCS uses θ_D)

The pattern is clear: wherever the Synchronism framework uses γ = 2T/θ_D as the "new variable," it's restating the Debye temperature. Any domain where standard physics ALSO uses θ_D will show circular improvement.

**When would N-S framing be genuinely non-circular?**

The framing could be non-circular in domains where:
1. Standard physics does NOT use θ_D directly
2. The multi-component structure (separate channels) adds coupling terms not in single-channel models
3. The "viscosity collapse" produces qualitatively different behavior than existing models predict

Candidates for non-circular applications:
- **Phonon drag thermoelectric effect**: where phonon and electron channels genuinely compete and drag on each other — existing theory uses separate Boltzmann equations, but cross-channel coupling terms are approximate
- **Multiferroics**: materials with coupled magnetic and ferroelectric order — two "channels" with different γ values, coupled via strain
- **Hydride superconductors (LaH₁₀, H₃S)**: θ_D is estimated, not measured; predictions might be non-circular here

---

## Status Update

**Phase 3 after 2 sessions:**

| Session | Domain | Verdict |
|---------|---------|---------|
| #1 | Channel independence, Prandtl analog | Organizational vocabulary |
| #2 | Superconductivity Tc prediction | Circular; N-S modification fails |

**Phase 3 emerging conclusion**: The N-S multi-component framing is consistent and explanatory but does not add predictive power beyond standard physics in the tested domains. This mirrors the Phase 2 conclusion but from a different angle: Phase 2 showed the MEASUREMENT of γ is circular; Phase 3 is showing the APPLICATION of γ to known physics is also circular.

**Open question for Session 3**: Is there ANY domain where multi-component N-S adds non-circular content? Phonon drag is the best candidate.

---

## Files

- **Simulation**: `simulations/chemistry/phase3_session2_superconductivity_viscosity_collapse.py`
- **Data**: 24 elemental superconductors (Tc, θ_D, λ_ep, γ_S)
- **Previous**: `Phase3_Session1_CFD_Channel_Independence.md`

---

*Phase 3 Session #2 — Chemistry Track*
*Finding: N-S viscosity framing circular for Tc prediction; productive failure documented*
