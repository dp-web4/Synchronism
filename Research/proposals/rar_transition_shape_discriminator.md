# Back-annotation: The Compander Yields a Distinct RAR Transition Shape — the First Non-Degenerate Galaxy-Scale Discriminator

**Source**: synchronism-site explorer session 2026-05-20
**Status**: constructive proposal — supersedes the scope of the 2026-03-30 "RAR disconnection" result
**Decisive test**: fit μ_Syn vs McGaugh ν to real SPARC RAR with γ fixed at 2; compare BIC

## Claim

The 2026-03-30 site finding "the coherence function cannot produce the RAR" is correct only under the
ν-identification (boost = 1/C(g_bar/a₀)), which gives constant g_obs in deep MOND and violates Tully-Fisher.
Under the **dual μ-identification** — C as the MOND interpolating function with argument g_obs/a₀,
g_bar = g_obs·tanh(γ ln(1 + g_obs/a₀)) — the compander reproduces **both** asymptotes:

- High acceleration: μ→1 ⇒ g_obs→g_bar (Newtonian).
- Low acceleration: μ≈γx ⇒ g_obs = √(a₀ g_bar / γ) (deep-MOND √-law, Tully-Fisher preserved).

So C(ρ) *is* a legitimate interpolating function. The disconnect was a mapping-choice artifact.

## The discriminator

After matching both asymptotes (γ absorbed into a₀), the compander differs from McGaugh's
ν(y)=1/(1−e^(−√y)) only in transition curvature at g_bar ≈ a₀:

- γ=2 (galaxy value, N_corr=1): max deviation **−0.083 dex at g_bar/a₀≈1.1**, = 1.45 × SPARC σ_int.
  Compander transitions *more sharply*; predicts *less* boost than MOND through the transition.
- RMS residual fitting compander (free a₀) to McGaugh-shaped data: **0.067 dex > σ_int = 0.057 dex**.

This is the same observable that already discriminates the "simple"/RAR μ from the "standard" μ in SPARC
(Famaey & McGaugh 2012). The compander is a new family member, testable identically.

## Critical contingency

A **free-γ** fit returns γ ≈ 0.91 (N_corr ≈ 5), RMS 0.013 dex — indistinguishable from McGaugh. The
discriminator has power **only if γ is independently pinned at 2** (the framework's own galaxy assignment).
Thus the entire discriminating test reduces to a prior question the framework has never resolved
operationally: **is galaxy-scale γ pinned by N_corr=1, or is it a fit parameter?**

- Fitted ⇒ compander = MOND, zero discrimination (collapses to γ≈0.9).
- Pinned at 2 ⇒ distinct, falsifiable RAR — and SPARC already mildly disfavors it (0.067 dex residual).

## Proposed test (zero-cost, existing data)

Fit μ_Syn(x)=tanh(2·ln(1+x)) and McGaugh ν to the real Lelli-McGaugh-Schombert 2016 SPARC RAR (2693 points,
153 galaxies), γ **fixed at 2**, M/L and distance nuisance parameters marginalized. Compare BIC.

**Kill criterion**: ΔBIC > 10 favoring McGaugh refutes the γ=2 compander as the galaxy mechanism. The
0.067-dex structured residual suggests this is the likely outcome — leaving the framework with a *fitted* γ
that is MOND by construction. Either result is decisive and is the first galaxy-scale test that is not
MOND-degenerate by sign or by the external field effect.

## Relation to prior conclusions

Consistent with Session #574 ("C(ρ) is a reparametrization of MOND ν(x)") — but sharper: the
reparametrization is *exact only at the best-fit γ≈0.9*. At the committed γ=2 it is a *distinct, mildly
refuted* function, not a reparametrization. The framework's claim to a galaxy mechanism stands or falls on
whether it commits to γ=2 before the fit.

---

## RESULT — proposed test executed on real SPARC (synchronism-site explorer, 2026-05-21)

The decisive test specified above was run against the real Lelli-McGaugh-Schombert (2016) mass models
(`simulations/sparc_real_data/MassModels_Lelli2016c.mrt`, 2807 points after a 10% velocity-error cut;
M/L_disk=0.5, M/L_bul=0.7; a₀ free; log-space least squares; ΔBIC differential).
Script + full finding: synchronism-site `explorer/scripts/rar_transition_shape_real_sparc.py`,
`explorer/findings/rar-transition-shape-real-sparc-result.md`.

| Model | a₀ (m/s²) | RMS (dex) | ΔBIC vs McGaugh |
|-------|-----------|-----------|-----------------|
| McGaugh ν (k=1) | 1.13×10⁻¹⁰ | 0.1437 | reference |
| Compander μ, **γ=2 pinned** | 2.97×10⁻¹⁰ | 0.1485 | **+184** |
| Compander μ, γ free | 5.3×10⁻¹¹ (γ=0.49) | 0.1437 | +7.1 |

**Verdict: kill criterion triggered.** ΔBIC = +184 ≫ 10 refutes the γ=2 compander as the galaxy
mechanism. Conservative correction for intra-galaxy point correlation (effective N≈500–1000) gives
ΔBIC ≈ 33 — still decisive. The γ=2 misfit is a coherent S-shaped residual (≈0.05–0.10 dex
peak-to-peak, ~8σ/bin) concentrated at g_bar≈a₀, exactly the predicted transition-shape signature; it
is a population-wide shape term at fixed g_bar/a₀, so it is not absorbable by per-galaxy M/L marginalization.

**Two corrections to the 2026-05-20 estimate:**
1. γ=2 is **decisively refuted**, not "mildly disfavored" — the per-point RMS penalty is small (+3.3%)
   but the structured residual is overwhelmingly significant over thousands of points.
2. The free-γ best fit is **γ ≈ 0.49** (not 0.91), with RMS identical to McGaugh to four digits — fully
   MOND-degenerate. The idealized-curve fit that gave 0.91 was an artifact of uniform sampling.

**Closure:** there is no γ for which the compander is both distinct from MOND and consistent with SPARC.
Pinned γ=2 → refuted; fitted γ → MOND. The framework's only non-degenerate galaxy-scale test has now
been run on existing data and eliminates the only non-MOND version of the mechanism. Net discriminating
galaxy tests vs MOND: **0, by execution.** Pipeline validated: best-fit McGaugh a₀=1.13×10⁻¹⁰ vs
canonical 1.20×10⁻¹⁰; RMS at canonical a₀ = 0.1441 ≈ best-fit 0.1437.
