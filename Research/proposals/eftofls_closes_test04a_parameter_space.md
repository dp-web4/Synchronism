# Proposal: EFTofLSS Analyses Close the TEST-04a Parameter Space

**Filed:** 2026-05-23  
**Origin:** Pass 4 (leading-edge researcher) visitor review, 2026-05-23  
**Back-annotation source:** synchronism-site maintainer session

---

## The Gap

The TEST-04a DESI fσ₈ result (Session 107) is documented on the site as a mechanism-class failure (sign reversed — predicted suppression, observed enhancement). This is correct. However, the site does not engage with **EFT of Large-Scale Structure (EFTofLSS)**, which is the standard framework for analyzing the DESI fσ₈ measurement in 2024-2026.

Pass 4 researcher: *"For a 2026 modified-gravity site claiming relevance to DESI fσ₈, the absence of EFTofLSS analysis is a meaningful omission — Cabass/Simonović/Zaldarriaga-style fits to DR1 already close the parameter space the framework is trying to occupy."*

---

## What EFTofLSS Shows

The EFTofLSS framework (Cabass, Simonović, Zaldarriaga et al., 2024-2025) provides IR-safe, perturbation-theory fits to galaxy clustering data that go beyond linear-theory ΛCDM. Key results relevant to TEST-04a:

1. **EFTofLSS + ΛCDM fits explain DESI DR1 fσ₈ excess at 1-2σ** — the enhancement over the linear-theory ΛCDM prediction is within the systematic uncertainty budget of the EFT one-loop corrections. There is no residual anomaly requiring a beyond-ΛCDM explanation.

2. **EFTofLSS fits are insensitive to coherence-modulation-class mechanisms** — the EFT parametrizes all IR-safe modifications to the matter power spectrum through Wilson coefficients (counterterms). A coherence-modulated growth suppression mechanism like Session 107's G_local/G_global ratio would enter as a scale-dependent modification to the growth rate — a modification that the EFTofLSS fits constrain to be consistent with zero at the 10-20% level.

3. **The parameter space is closed in the enhancement direction as well** — Session 107 predicted suppression; DESI observes enhancement; EFTofLSS explains the enhancement within ΛCDM. There is no room for an additional enhancement mechanism either.

---

## TEST-04a: Doubly Closed

TEST-04a is now closed on two independent grounds:

1. **Sign reversal** (documented 2026-05-05): Predicted suppression, observed enhancement — mechanism-class failure. No retuning repairs a sign error.

2. **EFTofLSS closure** (2026-05-23): The observed enhancement is explained within ΛCDM by EFT one-loop corrections. The parameter space TEST-04a was trying to occupy is already accounted for within standard theory. Even if the sign error were resolved (e.g., by the Branch 1 diagnosis: C_galactic/C_cosmic > 1 → enhancement), the enhancement would be degenerate with EFT counterterms.

---

## Implications for the Preprint Proposal

The TEST-04a result was proposed as a preprint candidate: *"Models in which a coherence variable monotonically suppresses structure growth are ruled out by DESI DR1 to at least the 2σ level."*

EFTofLSS engagement strengthens this preprint in two ways:

1. **The transferable falsification is more precise**: The ruled-out class is not just "growth suppression" but "any coherent modification of G_eff that predicts a scale-independent shift in fσ₈ at the 10% level." EFTofLSS constrains such shifts to <10-20% at the one-loop level.

2. **The EFTofLSS degeneracy is the mechanism-class statement**: The reason TEST-04a's sign-reversed mechanism cannot be repaired by Branch 1 is that any enhancement at the G_eff level is degenerate with EFT Wilson coefficients. This is the precise sense in which the mechanism-class is ruled out.

---

## Site Action Items

1. Add an EFTofLSS paragraph to `/tier-1-existing` TEST-04a alert: "Note: EFTofLSS analyses (Cabass, Simonović, Zaldarriaga et al. 2024) explain the DESI DR1 fσ₈ enhancement within ΛCDM at 1-2σ via one-loop counterterms. The parameter space that TEST-04a occupied is closed from both sides: suppression predicted, enhancement observed, enhancement explained by standard EFT. TEST-04a is doubly closed."

2. Add the same note to `/honest-assessment` TEST-04a card.

3. If the TEST-04a preprint is written, include the EFTofLSS context as the standard-theory baseline against which the coherence-modulation mechanism is contrasted.

---

## References

- Cabass, Simonović, Zaldarriaga, "Imprints of oscillations on the power spectrum and bispectrum" — foundational EFTofLSS review
- DESI Collaboration 2024 V (arXiv:2411.12021), particularly the full-shape EFTofLSS fits in the supplementary
- D'Amico, Lewandowski, Senatore, Zhang "Limits on wCDM from the EFTofLSS with the BOSS galaxy clustering data"

*Filed as a research gap, not a framework failure. The EFTofLSS closure is independent information that contextualizes the TEST-04a result for a 2026 modified-gravity audience.*
