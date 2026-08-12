# Direct fit of the w(z; γ) family to DESI DR2 + CMB + SN: substituted branch ≡ ΛCDM (γ = 0.487 ± 0.02); both covariant completions excluded by the fit

**Date**: 2026-08-12
**Source**: site explorer track, topic `fit-the-gamma-family-to-desi-chains.md`
**Script**: `synchronism-site/explorer/findings/scripts/fit_gamma_family_to_desi_dr2.py` (+ output)
**Finding**: `synchronism-site/explorer/findings/gamma-family-direct-fit-desi-dr2-substituted-is-lcdm-covariant-excluded.md`
**Supersedes in part**: the CPL-space pricing in `covariant_00_component_sign_lock_dies_desi_nogo_hardens_20260811.md` (the derivations there stand; the exclusion-flavored numbers are re-priced here at likelihood level)

## Data (all verified against sources 2026-08-12)

- DESI DR2 BAO: 13 measurements incl. per-tracer D_M–D_H correlations (arXiv:2503.14738, as tabulated by arXiv:2506.17926 Table 1)
- Planck 2018 distance priors (R, l_A, ω_b): Chen, Huang & Wang arXiv:1808.05724 Table 1
- DES-Dovekie SN Hubble diagram, 1820 SNe, full STAT+SYS covariance, M marginalised (des-science/DES-SN5YR, the DES group's 2025 recalibration; their companion arXiv:2508.10514 reports Nσ ∼ 1.1–2.3 for dynamical DE from late-time data alone)
- Projection-paper IDs verified: Shlivko & Steinhardt 2405.03933; Cortês & Liddle 2404.08056; Wolf, García-García, Bartlett, Ferreira 2408.17318; Wolf, García-García & Ferreira 2502.04929

Pipeline validated before any framework statement: ΛCDM BAO χ² = 12.6/13, r_d = 147.4 Mpc; w₀wₐ best fit BAO+CMB = (−0.44, −1.67) vs DESI's published (−0.42, −1.75); full-combo crossing preference Δχ² = −11.0 (~2.9σ) inside the published 2.8–4.2σ band.

## Results

| Model | Δχ² vs ΛCDM | Note |
|---|---|---|
| w₀wₐCDM (+2 params) | −11.0 | the crossing preference, this data combination |
| substituted family (+1 param) | **−0.3** | best γ = 0.487 (−0.021/+0.024); nests ΛCDM at γ = 1/2 |
| completion B, ω ∈ {0,1,5,50} | **+79 … +187** | best member per ω; +18 already without SN |
| completion A (exact EdS) | χ² ≈ 9,900 (BAO+CMB) | quantified; the pre-1998 background |

1. **Substituted branch: the "3.4–6.3σ forced-wₐ" pricing (08-11) did not survive execution.** The likelihood never forces the family to DESI's w₀; the family sits at its Λ-corner (best member projects to CPL (−0.993, +0.023)) and pays exactly ΛCDM's price: +11.0 behind w₀wₐCDM. 5th over-refutation instance (no-go stands, exclusion-number dies).
2. **Covariant class: hardened from "cannot reach the crossing" to "fails the fit".** Completion B cannot match even ΛCDM (Δχ² ≥ +79, every ω, conservative — compressed CMB recovers only ~2/3 of the crossing evidence and would punish B harder in full). With A exactly EdS, **no covariant member of the family survives DR2-era data. The dark-energy sector survives only in its non-covariant form, and only by being ΛCDM.**
3. **First executed cross-sector consistency test: γ_cosmo = 0.487 ± 0.02 vs γ_galaxy(SPARC) = 0.489 — pass at 0.1σ.** Deflationary reading mandatory: γ = 1/2 is exactly Λ (Möbius), γ = 0.489 is exactly MOND simple-μ; the two sectors' standard models sit 0.011 apart in γ-space, so the test had no power to fail. Separating 0.489 from 0.500 needs σ_γ ≈ 0.004. SPARC-side σ(γ) has never been derived — needed before the concordance can be priced.

## Consequence for TEST-26 (registration draft, dp-gated)

As drafted, TEST-26 adjudicates quadrant membership. On the only surviving branch, quadrant membership is inherited from ΛCDM wholesale — **TEST-26 is ΛCDM-degenerate on every branch still alive** (the framework-specific branches are already excluded by DR2-era data, before DR3 arrives). Same non-discrimination class as the a₀(z) row. The registration should: (a) move the statistic from quadrant membership to Δχ²(substituted best fit vs w₀wₐCDM best fit) on DR3 likelihoods; (b) state that the family's fate is ΛCDM's fate — DR3 kills it exactly when it kills ΛCDM's background, no sooner and no later; (c) present it as a consistency check, not a discriminator.

## Count

UNCHANGED at 6 (dp-gated recount pending). Audit-fit, not a registered kill.

## Methodological note

An uncalibrated 2.5% bias in the sound-horizon integral (fitting-formula z_d vs CAMB) initially shifted ω_m 6% high and cost more χ² than the entire dark-energy signal. Per-point pull diagnostics caught it; the pipeline is CAMB-calibrated at the Planck 2018 fiducial (two constants, residual < 0.2%, shared by all models). Rule confirmed: diagnose pulls before trusting verdicts.
