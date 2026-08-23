# Pressure-supported dwarfs: the boost ceiling screens EFE = 0 from test; the DF2 repair invalidates the EFE = 0 derivation

**Date**: 2026-08-23 · **Origin**: synchronism-site explorer lane
**Source**: `synchronism-site/explorer/findings/pressure-supported-dwarfs-the-ceiling-screens-efe-zero-from-test.md`
**Refutation count: UNCHANGED at 6.** New-sector evaluation of the existing boost-ceiling root.

## Result

Executed the framework's own ceiling against pressure-supported dwarfs (Wolf+2010 mass estimator;
`B_req = (σ_obs/σ_N)²`; M/L swept as a band, not marginalised):

| system | B_req | M/L band [1,4] | vs 3.17 | vs 13.7 |
|---|---|---|---|---|
| Crater II | 60.2 | [30, 120] | 4.3σ | 3.5σ band-robust |
| Draco | 101.5 | [51, 203] | 3.7σ | 3.3σ band-robust |
| Sculptor | 26.7 | [14, 43] | 2.9σ | 1.6σ |
| Fornax | 5.6 | [3.5, 9.3] | 2.8σ | passes |
| NGC 1052-DF2 | 1.9 | [1.3, 2.6] | passes | passes |

Crater II is the only discriminator: MOND+EFE a priori (**McGaugh 2016, ApJL 832, L8** — not
Milgrom) gives 2.1 (+0.9/−0.6) vs measured 2.7 ± 0.3 (Caldwell+2017) = **0.6σ, consistent**, while
the framework's ceiling caps σ at 1.29 km/s = **4.7σ short**. Draco is a *shared* failure — MOND
misses it too (2.2–4.3 vs 9.1) — so the class must not be cited collectively.

## The structural point

EFE = 0 is the framework's **most favourable** case on these systems (an external field only lowers
σ), yet the ceiling fails anyway. And EFE matters only where `g_ext ≳ g_int`, i.e. in sparse
satellites — exactly where `B_req` is largest. **The bounded boost structurally screens EFE = 0 from
test.** The one system where the ceiling passes and the field is strong (DF2) has its prediction set
entirely by the unresolved coupling fork: 22.5 km/s (`1/C` branch, 6.1σ high) vs 6.09 km/s
(`V_flat·C` branch, 1.0σ low).

## Internal divergence requiring adjudication

Two incompatible DF2 repairs exist in this repository, neither propagated to the site:

1. `manuscripts/arXiv_preprint_draft_v1.md` §5.1–5.2 — standard model predicts **σ ~ 80 km/s** vs
   8.5 observed; repaired by **formation coherence**, `C_eff = max(C(ρ_local), C_formation)`.
2. `docs/whitepaper/Synchronism_Whitepaper_Complete.md`, Session #97 — repaired instead by **tidal
   stripping leaving a high-C core**, which contradicts `C(ρ)` at DF2's measured density (the
   quantity §5.1 computes as `C ~ 0.04`).

**Consequence.** The site derives EFE = 0 from the premise *"C depends only on local ρ."*
`C_eff = max(C(ρ_local), C_formation)` makes C a formation-history variable, and Session #97's
repair makes internal dynamics depend on the host at ~80 kpc. **Either repair invalidates the
EFE = 0 derivation.** This needs adjudication in the archive, not just on the site.
