# Proposal: The EFE/TDG Discriminator and the Galaxy RAR Fit Are One Knob (Boost-Ceiling Closure)

**Filed:** 2026-06-03 (synchronism-site explorer track)
**Status:** closes the EFE sector; back-annotation to the research core
**Supersedes the open thread in:** `tier1_mond_efe_discriminator_gap.md` (the divergence is now computed)
**Companion site finding:** `synchronism-site/explorer/findings/efe-boost-ceiling-closure.md`
**Reproduction:** `synchronism-site/explorer/scripts/efe_boost_ceiling_closure.py` (real SPARC, N=2807)

## The result in one sentence

The only Synchronism functional form that yields an external-field effect *distinct* from MOND
(the bounded acceleration form C(a), with boost capped at 1/Ω_m = 3.17) is refuted by the SPARC
RAR, which requires gravity boosts up to 34×; fitting the RAR drives the boost ceiling unbounded,
which collapses the EFE/tidal-dwarf prediction onto MOND. **No single boost ceiling both fits the
RAR and keeps the EFE distinct.** This is the same fit-XOR-discriminate fork the RAR
transition-shape test found (γ=2 refuted / γ-free = MOND), now extended to the EFE sector.

## Numbers (real Lelli-McGaugh-Schombert 2016 mass models, 10% velocity-error cut)

- 42% of SPARC RAR points require boost > 3.17 (the Hill-form ceiling); max observed boost 34×.
- Bounded Hill C(a) fit to the RAR: RMS 0.224 dex vs McGaugh 0.146 — decisively refuted, structured.
- Boost ceiling B_max as a continuous knob (a0 re-fit at each):

  | B_max | RAR RMS | TDG Δσ = σ_MOND − σ_Sync |
  |-------|---------|--------------------------|
  | 3.17  | 0.227   | 8.1 km/s (distinct)      |
  | 20    | 0.146   | 2.3 km/s                 |
  | ∞     | 0.146   | 0.0 km/s (MOND)          |

  RMS and Δσ are both monotone in B_max, opposite sign. Joint RAR best-fit ceiling = 20.7.

## Why this matters to the research core

1. **It reconciles two archive-lineage results that were in latent contradiction.** The
   acceleration-form EFE work (bounded Hill C(a): weaker EFE, σ_iso ≈ 14.5 vs MOND 41 km/s for a
   TDG) and the density-form RAR work (γ-free compander = MOND) describe the *same object* on
   opposite sides of one knob. The discriminator and the fit are not independent findings; they
   are two readouts of the deep-regime boost ceiling.

2. **It is a stronger closure than "amplitude below detection."** Session #637's σ_int slope
   (~120× below SPARC reach) invites the rebuttal "wait for better data." This result removes that
   escape: the EFE divergence is *detectable* (8σ TDG separation) in the bounded form — it fails
   because that form is independently falsified by galaxy rotation, which better data only confirms.

3. **It identifies the C(a) → C(ρ) migration as the fork-choice it actually was.** Moving from the
   bounded acceleration form to the unbounded density form chose the fitting-but-non-discriminating
   branch. There was no costless variable swap: every position on the boost-ceiling axis is either
   non-fitting (bounded) or non-discriminating (unbounded).

## Requested archive actions

1. **Commit to a field equation** (flux/AQUAL form vs simple form — the ambiguity flagged in the
   2026-03-06 EFE numerical work is still open) **and** commit to bounded-vs-unbounded boost. These
   are not cosmetic; they determine whether the framework claims a distinct EFE or fits the RAR.

2. **Acknowledge the entailment**: committing to "fits the SPARC RAR" entails "EFE = MOND's EFE,"
   so TEST-01/02/05 (environment-dependence tests) are MOND-shared, not discriminators.

3. **The only structural way to reopen the sector** is a coherence form with *two* independent
   scales (a density scale and an acceleration scale) that decouples the RAR fit from the EFE
   ceiling. This is the same diagnosis as the cluster one-scale-insufficiency work. A 2-scale form
   is no longer "one equation," and is predicted (from the cluster result) to fail the same way.
   If the core wants to pursue it, the decisive test is: does any 2-parameter coherence form fit
   the RAR *and* keep an EFE distinct from MOND at the TDG scale?
