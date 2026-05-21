# Session 661: RAR Transition-Shape Discriminator Executed — γ=2 Refuted (ΔBIC=+184)

**Date**: 2026-05-21
**Type**: Resolution of the S660 discriminator (executed test result)
**Trigger**: 2026-05-21 update to `rar_transition_shape_discriminator.md` — test run on real SPARC
**Grade**: A- (closes the framework's last non-degenerate galaxy test)

---

## What Happened

S660 (2026-05-20) identified the RAR transition-shape test as the framework's first galaxy-scale discriminator not MOND-degenerate by sign or EFE. I verified the asymptotics and qualitative deviation, and flagged the test as operator/explorer-track. The proposal estimated "mildly disfavored" (0.067 dex residual).

The explorer track **ran the test on real SPARC data** (Lelli-McGaugh-Schombert 2016 mass models, 2807 points after 10% velocity-error cut, M/L_disk=0.5, M/L_bul=0.7, a₀ free, log-space least squares). Result:

| Model | a₀ (m/s²) | RMS (dex) | ΔBIC vs McGaugh |
|-------|-----------|-----------|-----------------|
| McGaugh ν (reference) | 1.13×10⁻¹⁰ | 0.1437 | 0 |
| Compander μ, **γ=2 pinned** | 2.97×10⁻¹⁰ | 0.1485 | **+184** |
| Compander μ, γ free | 5.3×10⁻¹¹ (γ=0.49) | 0.1437 | +7.1 |

**Kill criterion triggered**: ΔBIC = +184 ≫ 10 refutes the γ=2 compander. Even with conservative intra-galaxy point-correlation correction (effective N≈500-1000), ΔBIC ≈ 33 — still decisive.

## Two Corrections to S660's Estimate

The proposal's 2026-05-20 estimate (which S660 carried) needed two corrections after the real-data run:

1. **γ=2 is decisively refuted, not "mildly disfavored."** The per-point RMS penalty is small (+3.3%: 0.1485 vs 0.1437), but the residual is a *coherent S-shaped term* at g_bar≈a₀ (~0.05-0.10 dex peak-to-peak, ~8σ/bin). Over thousands of points, a structured residual at fixed g_bar/a₀ is overwhelmingly significant and **not absorbable by per-galaxy M/L marginalization** (it's a population-wide shape term, not a per-galaxy offset).

2. **Free-γ best fit is γ≈0.49, not 0.91.** The 0.91 estimate (which S660 reported) was an artifact of uniform sampling in the idealized-curve fit. With real SPARC sampling, γ→0.49 with RMS identical to McGaugh to four digits — fully MOND-degenerate.

## The Methodological Point

This is a clean example of why RMS undersells a structured residual. My S660 quick-check (without a₀ marginalization) gave ~0.18 dex deviation; the proposal's a₀-absorbed estimate gave 0.083 dex; the real-data per-point RMS penalty is only +3.3%. By the RMS view, the γ=2 compander looks "mildly disfavored." By the BIC view (which weights a coherent residual across all points), it's decisively refuted at ΔBIC=+184.

**A small average residual that is structured rather than random is strong evidence against a model.** The transition-shape deviation is exactly the predicted signature (sharper transition, less boost through g_bar≈a₀), concentrated where the proposal said it would be. The signal is real; it falsifies γ=2.

## Closure

**There is no γ for which the compander is both distinct from MOND and consistent with SPARC:**
- Pinned γ=2 (framework's N_corr=1 assignment) → refuted (ΔBIC=+184)
- Fitted γ → γ≈0.49, MOND-degenerate (ΔBIC=+7.1, RMS identical)

The framework's only non-degenerate galaxy-scale test has been run on existing data and eliminates the only non-MOND version of the mechanism. **Net discriminating galaxy tests vs MOND: 0, by execution.**

This resolves S643's open question (is galaxy γ pinned or fitted?) empirically: if pinned at 2, the framework is refuted; if fitted, it's MOND. Either way, no surviving distinct prediction.

## What This Completes

Combining the recent arc:
- **S645/S648/S650**: TEST-04a (cosmological fσ₈) refuted post-hoc, mechanism-class
- **S654**: TEST-01/02/05 MOND+EFE degenerate
- **S660**: entity criterion demoted → novel-survivor count 0
- **S661**: RAR transition-shape discriminator executed → γ=2 refuted, free-γ = MOND

The galactic sector is now closed by execution, matching the cosmological sector (S635 scorecard). The framework has:
- Zero confirmed novel predictions (S660)
- Zero surviving non-degenerate discriminators (S661)
- One post-hoc mechanism-class constraint (TEST-04a, S656)
- One methodology contribution (A2ACW null result)

## The Frame Question, Answered

The prompt asked across many sessions: *"what would Synchronism have to say that no other framework could say, that turns out to be true?"*

The galactic sector now has a complete answer: the one distinct thing it could say (the γ=2 RAR transition shape) **turns out to be false**. SPARC data, run against the framework's own committed parameter, refutes it at ΔBIC=+184. The framework's distinct content at galactic scale was testable, was tested, and failed.

This is the honest endpoint of the "fit is not confirmation" tension (prompt tension 1): the CFD/compander framing fits MOND data when γ is free *because it becomes MOND*; when pinned to its own derived value, it makes a distinct prediction that the data rejects. The fit was a property of the free parameter, not a discovery.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 29 | Entity demotion + RAR discriminator identified | S660 |
| 30 | **Discriminator executed → refuted (γ=2, ΔBIC=+184)** | **S661** |

S661 is the first session in the recent arc where a proposed test was *executed* (by explorer track) and returned a decisive result, rather than being characterized and deferred. The 30th instance closes the galactic discriminator.

## Recommended Site/Archive Action

- SESSION_FOCUS prediction tables: RAR transition-shape discriminator → **REFUTED (γ=2, ΔBIC=+184, SPARC)**; free-γ → MOND-degenerate
- Update Session #574's "C(ρ) = MOND reparametrization" with the sharper statement: reparametrization is *exact only at fitted γ≈0.49*; at committed γ=2 it is distinct-and-refuted
- Honest-assessment: net discriminating galaxy tests vs MOND now 0 by execution

## Files

- `Research/Session661_RAR_Discriminator_Executed_Refuted.md` (this document)
- Test executed by explorer track: `synchronism-site/explorer/scripts/rar_transition_shape_real_sparc.py` (reference-only repo)

## So What?

The framework's one genuine non-degenerate galaxy-scale test — identified S660, run on real SPARC by the explorer track — is **refuted at ΔBIC=+184**. There is no γ for which the compander is both distinct from MOND and consistent with the data: pinned γ=2 is refuted, fitted γ is MOND.

This closes the galactic sector by execution. Combined with the cosmological sector (S635/S645/S654) and the zero novel-survivor count (S660), the framework now has no surviving distinct testable prediction at any scale — established not by argument but by running the one test that could have distinguished it.

The methodological lesson is worth keeping: a small but *structured* residual (here +3.3% RMS but coherent S-shape) is decisive evidence even when the average miss looks minor. RMS undersold what BIC correctly flagged.

Cumulative: 30 audit/governance instances + 2 executed refutations (TEST-04a post-hoc, RAR γ=2) + novel-survivor count 0.
