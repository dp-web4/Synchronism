# Proposal: Deriving the RAR σ_int(ρ_env) Slope

**Date**: 2026-04-28  
**Source**: Maintainer session — visitor Pass 4 (leading-edge researcher) and Pass 3 (grad student)  
**Priority**: High — this is the framework's one candidate for a genuinely novel, non-reparametrization prediction

---

## The Gap

The Synchronism site claims, on `/galaxy-rotation` and implicitly on the landing page, that the framework predicts **environment-dependent intrinsic scatter on the RAR (Radial Acceleration Relation)**. The claim is:

> Standard MOND predicts narrow RAR scatter (Lelli et al. 2017: σ_int ≈ 0.13 dex across SPARC). Synchronism predicts σ_int varies with local coherence density ρ_env.

A leading-edge reviewer visiting the site today (2026-04-28) identified this as:

> "The framework's one genuinely novel claim that is potentially distinguishable from MOND and from ΛCDM-with-feedback. Not yet a paper; possibly a paper if formalized."

The prediction is referenced on the site but **lacks a numerical slope**: no dσ_int/d(log ρ_env), no specified scale, no comparison to the ΛCDM-hydrodynamics baseline (Santos-Santos et al. 2020 predict a different environmental signature from baryonic feedback).

Without numbers, this is not a falsifiable prediction by the site's own taxonomy — it would be "Untested" at best, "Speculative" more honestly.

---

## The Research Question

**Can the framework derive a numerical prediction for dσ_int/d(log ρ_env)?**

Specifically:
1. What does C(ρ) predict for how the coherence field modifies Vflat in low-density vs. high-density environments?
2. Given the tanh form C(ρ) = tanh(γ · log(ρ/ρ_crit + 1)), how does σ_int depend on the variance of ρ_env within the sample?
3. Does the prediction produce a slope that is (a) observationally distinguishable from σ_int = constant (MOND's prediction), (b) distinguishable from the feedback-driven slope in ΛCDM hydrodynamics, and (c) consistent with the SPARC baseline at σ_int ≈ 0.13 dex?

---

## Observable Target

SPARC + environment cross-match (Lelli et al. 2016 + Tully-Fisher galaxy environment from ALFALFA-SDSS-SDSS filament catalogs, e.g., Yang et al. 2007, Cortese et al. 2011). Split SPARC into:
- **Isolated** (> 1 Mpc from nearest neighbor with M > 10^10 M_sun)
- **Group environment** (< 0.5 Mpc, N > 3)
- **Cluster outskirts** (> 1 Mpc from BCG but inside R_200)

Measure σ_int in each subsample. MOND (density-independent) predicts the same σ_int ≈ 0.13 dex in all three. ΛCDM-hydro predicts higher σ_int in clusters due to tidal stripping and ram pressure. **Synchronism should predict something specific and different from both.**

---

## The Derivation Path

The coherence density ρ in the galactic context is operationalized as the matter density felt by the coherence field. In low-density environments, ρ < ρ_crit → C(ρ) is below its half-saturation point → Vflat is less constrained → higher scatter. In high-density environments, ρ >> ρ_crit → C(ρ) → 1 → Vflat is tightly constrained → lower scatter.

This predicts **σ_int decreases with ρ_env** — opposite to the ΛCDM-feedback direction (which predicts σ_int *increases* in clusters due to stripping).

If this sign is correct, it is a testable discriminator from ΛCDM. But the amplitude requires:
1. The distribution P(ρ_env) across the three environment bins
2. The propagation of ΔC(ρ_env) → ΔVflat → Δσ_int

This calculation should be possible within the framework given the existing C(ρ) parameterization. Session #64-65 (Jeans-criterion derivation) contains the machinery for linking ρ to galactic observables.

---

## Why This Matters

The current site reduces the galactic story to:
> "MOND + Vflat calibration + one novel ansatz (density-dependent scatter)"

If the slope can be derived, this becomes:
> "MOND + one *derived* prediction that disagrees with feedback ΛCDM and is testable against existing SPARC + environment catalogs"

That is the difference between a reparametrization and a research paper. The data exists (SPARC + ALFALFA-SDSS cross-match). The prediction can be derived from existing framework components. This is the highest-leverage research task facing the project.

---

## Suggested Session Structure

1. Read Sessions #64-65 (Jeans-length derivation, A = 4π/(α²GR₀²) chain) and extract the ρ → Vflat mapping
2. Compute ΔVflat/ΔC(ρ) at ρ = ρ_crit (steepest gradient point)
3. Estimate typical P(ρ_env) variance across the three environment bins
4. Propagate to Δσ_int and compare to:
   - MOND baseline: σ_int = 0.13 dex (constant)
   - ΛCDM-hydro: σ_int increases in clusters
   - Synchronism: σ_int decreases in clusters (predicted sign)
5. Check whether the predicted amplitude is above the SPARC measurement floor (~0.02 dex)

---

## Connections

- Archives: Sessions #64-65 (Jeans criterion), Session #193 (BTFR regime-dependent), GAMMA_UNIFICATION.md
- Site pages affected: `/galaxy-rotation`, `/key-claims`, `/top-5-tests`, `/predictions`
- Related proposal: `cosmology_claim_status_after_audits.md` (frames the open question after MOND reparametrization diagnosis)
- Related explorer topic: `delta-bic-galaxy-environment-ansatz.md` (open)

---

*Filed by maintainer track. The site's public dialogue with a leading-edge reviewer surfaced this as the framework's highest-leverage research direction on 2026-04-28.*
