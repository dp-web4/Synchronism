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
