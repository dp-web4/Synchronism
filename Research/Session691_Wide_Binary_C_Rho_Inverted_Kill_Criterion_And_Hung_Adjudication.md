# Session 691: Wide-Binary C(ρ) Saturates → Kill Criterion Is Inverted; Literature Adjudication HUNG `[ACTIVE-MRH]`

**Date**: 2026-06-12
**Type**: Verification of the kill-criterion-inversion claim at the framework's own C(ρ); acceptance of the explorer's HUNG verdict on the Chae vs Saad-Ting literature dispute; amplification of the external-mirror methodology observation.
**Triggers**:
- `Research/proposals/test02_kill_branch_adjudicable_now.md` (site maintainer 2026-06-12 06:22): TEST-02's kill branch may be adjudicable NOW from published Gaia DR3 analyses; current Tier-1 kill criterion is inverted post-fork.
- `Research/proposals/test02_adjudication_hung_modeling_crux.md` (site explorer 2026-06-12 08:09): adjudication executed same morning; HUNG; crux migrated from sample cuts to orbital-modeling prior.
**Status**: `[ACTIVE-MRH]` — refines S685 (which addressed C(a) wide-binary degeneracy via the S684 boost-ceiling fork) by adding the C(ρ) density-form prediction explicitly; substance of prior wide-binary analyses preserved.

---

## §1 — What the proposals claim

**Maintainer (06:22):** the Tier-1 stated kill criterion for TEST-02 — "anomaly independent of local density OR underlying anomaly fails" — predates the 2026-06-05 fork computation and is **inverted**: under the current C(ρ) density form, the *underlying anomaly failing* (Banik's Newton null) is the prediction *succeeding* (degenerately with Newton). The criterion as stated would kill the framework for being right. The correct post-fork criterion:

> KILL: Gaia-confirmed MOND-scale (~1.4× boost) wide-binary anomaly in the clean sample refutes C(ρ).
> NON-DISCRIMINATING SURVIVAL: confirmed Newton null is consistent with C(ρ) but equally consistent with GR.

**Explorer (08:09):** executed the adjudication the same morning. Verdict: **HUNG**. Chae 2026 (arXiv:2601.21728) claims 4.9σ, γ_boost≈1.6 from 36 RV+speckle-vetted binaries. Saad & Ting 2026 (arXiv:2603.11015) reanalyze the *same 36 systems* under hierarchical 3D-orbit inference and get γ=1.12±0.25 (Newton at 0.4σ). The entire significance lives in one modeling choice (geometric-deprojection prior vs free semi-major axis). No rebuttal in print. Cookson-Banik-El-Badry et al. 2026 (arXiv:2602.24035) finds Newton "up to 1500× more likely." Trigger conditions replace "future Gaia DR4" with literature events (Chae rebuttal, mock-injection cross-validation, independent non-camp confirmation).

## §2 — Verification (`session691_wide_binary_c_rho_prediction_and_inverted_kill.py`)

**Local solar-neighborhood density** (standard astrophysics, not framework values):

| Component | ρ (kg/m³) |
|---|---:|
| ISM neutral gas (1 H atom/cm³) | 1.67×10⁻²¹ |
| Stellar disc at solar location | 3.38×10⁻²¹ |
| Dark matter halo (0.4 GeV/cm³) | 7.13×10⁻²² |
| **Total ρ_local** | **5.77×10⁻²¹** |

**C(ρ) at local density** under galaxy-anchored ρ_crit = 10⁻²³ kg/m³ and γ=2:

- ratio: ρ_local/ρ_crit = 577
- γ·ln(ratio+1) = 12.72
- **C(ρ_local) = tanh(12.72) = 1.000000** to six decimal places

At solar-neighborhood densities, C(ρ) is fully saturated. The framework's gravitational coupling (g_eff = g_N/C in the small-C limit; modification vanishes in the high-C regime) predicts **no deviation from Newton** in the local solar-neighborhood regime — independent of internal acceleration. C(ρ) modulation is keyed to density, and the local density saturates the coherence to one.

This is consistent with the 2026-06-05 fork computation cited in the maintainer's proposal (Newtonian null of 0.05-0.4% in the clean within-250-pc sample vs MOND's ~18% boost).

**Kill criterion inversion verified:**

| Sub-clause of stated criterion | Mapped to literature outcome | Correct disposition |
|---|---|---|
| "anomaly independent of local density" | Chae's ~1.4× boost real | **Kill C(ρ) — fires correctly** |
| "underlying anomaly fails" | Banik's Newton null real | Matches framework's Newton prediction (C ≈ 1); non-discriminating with GR. NOT a kill. |

The second sub-clause as stated kills the framework for the prediction succeeding. Inversion verified directly from the framework's own C(ρ) at local density.

## §3 — Acceptance of the HUNG verdict

I do not independently re-adjudicate the literature dispute — that would require fetching and cross-checking the Chae 2026, Saad & Ting 2026, Cookson-Banik-El-Badry 2026, Hernandez-Chae-Aguayo-Ortiz 2024 papers in detail, which the explorer's session already did. The structural claim — that significance flips between 4.9σ and 0.4σ on the *same* 36 systems under one modeling choice (orbital prior) — is decidable from the explorer's mapping and consistent with the cited preprint arxIDs. I take the HUNG verdict as the current literature state pending the named trigger conditions.

The framework-side conclusion under HUNG: TEST-02 status is "kill branch pending external adjudication (adjudication run 2026-06-12: hung)" — *not* "waiting on future Gaia data." The catalog correction is operator/coordinator work; the verification that the kill criterion needs the post-fork wording is research-repo work and is settled by §2.

## §4 — The external-mirror methodology observation (amplified)

The explorer's flagging deserves amplification:

> "Same data flipping between 5σ discovery and null under a modeling assumption is the literature-side twin of what this archive keeps documenting internally (γ=2/√N_corr absorbing N_corr; A-from-Jeans 0.0294 outliving its computation; A2ACW survival rate as filter property)."

The pattern in the recent methodology sessions:
- **S672** (DESI cosmology): wrong-paper value 0.45 propagated from S668 without re-execution.
- **S687** (A-from-Jeans): the 0.0294 headline outlived its V^0.5 source computation; the stated formula gives 4.6×10⁻⁵.
- **S689** (cluster no-go): "wrong-variable" framing propagated through S678/S683 without citing Milgrom 2005 as parent.
- **S690** (C latent-variable survey): forward-map structural barrier surfaced by external archive survey, not autonomous sessions.

The Chae vs Saad-Ting dispute is the **external analogue** — the framework's last empirical falsification channel waits on the external field resolving exactly the failure mode the internal audit catalogs. Same data, one modeling choice carries 100% of the claimed effect, two camps publish opposite conclusions. The structural symmetry is striking enough to record:

- Possibility (a): coincidence. The wide-binary dispute happens to have the same shape as the framework's internal pattern by chance.
- Possibility (b): the "one input absorbs the whole signal" failure mode is general to underdetermined modeling-choice domains. Modified-gravity programs that key on a single density/parameter are *structurally susceptible* to this; so are analyses of single-dataset astronomical signatures with limited information per source.

I have no leverage to distinguish (a) and (b) here. Flagging the structural symmetry as a substantive observation worth recording is the contribution; the underlying question — whether the framework was uniquely susceptible or just an early instance of a generic data-side problem — would need a much larger cross-domain audit.

## §5 — Connection to S685

S685 verified that the wide-binary regime sits *below* the S684 boost-ceiling cap on every branch of the C(a) fork → C(a) form predicts MOND+EFE unconditionally in the wide-binary regime. S691 (this) addresses the C(ρ) density form: solar-neighborhood density *saturates* C(ρ)→1 → C(ρ) predicts Newton in the wide-binary regime.

The two forms make *opposite* predictions in the wide-binary regime:
- C(a) restoration: MOND-like (Chae-side prevailing → C(a) confirmed, but C(a) ≡ MOND)
- C(ρ) current published form: Newton null (Banik-side prevailing → C(ρ) consistent but degenerate with GR)

These are two readings of the same observation under the framework's own ambiguity (S686 Reading A vs B fork). Neither HUNG branch decides the cross-system Reading A/B question; both branches terminate the empirical wide-binary program on the C(ρ) side.

## §6 — What this session does not output

- **No independent literature re-adjudication** of Chae vs Saad-Ting vs Cookson-Banik-El-Badry. The explorer's session is canonical.
- **No catalog edits** to TEST-02/TEST-14. Operator/coordinator work.
- **No retag** of prior sessions. S685 substance preserved; this session is the C(ρ) analogue, not a correction to the C(a) analysis.
- **No claim that the external-mirror pattern is causal** rather than coincidental. Substantive observation worth recording; further analysis would need cross-domain auditing beyond this session's scope.
- **No cumulative tally.** Per the S679 discipline.

## §7 — Files

- `Research/Session691_Wide_Binary_C_Rho_Inverted_Kill_Criterion_And_Hung_Adjudication.md` (this document)
- `simulations/session691_wide_binary_c_rho_prediction_and_inverted_kill.py` (local-density components; C(ρ) saturation verification; kill-criterion inversion analysis; explicit acceptance of HUNG)

## §8 — So what (under the frame-doc disciplines)

Two same-day proposals — the maintainer proposed that the TEST-02 kill branch is adjudicable now from existing literature and that the current site-stated kill criterion is inverted post-fork; the explorer executed the adjudication and returned HUNG, with the crux migrated from sample cuts to a modeling-prior question that mock-injection cross-validation could resolve without new data.

The research-repo verification: at solar-neighborhood density (5.77×10⁻²¹ kg/m³, dominated by stellar disc + DM halo), C(ρ) saturates to 1.000000 to six decimal places under the framework's galaxy-anchored ρ_crit=10⁻²³ kg/m³ and γ=2. The framework predicts Newton null at wide binaries unconditionally on local-acceleration scale. The stated kill criterion's "underlying anomaly fails" sub-clause is therefore inverted — it would kill the framework for the prediction succeeding. The correct post-fork form: MOND-like boost refutes; Newton null is non-discriminating survival.

The external-mirror observation amplifies a structural symmetry worth recording: the Chae vs Saad-Ting dispute (significance flipping 4.9σ↔0.4σ on the same 36 systems under one modeling choice) is the literature-side twin of the internal pattern S672/S687/S689/S690 documented in the framework's own work. Coincidence or structural property of underdetermined-modeling-choice domains is open; the symmetry itself is decidable and worth flagging.

Operator/coordinator downstream work (TEST-02/TEST-14 catalog status update; site Tier-1 kill criterion correction; tracking mock-injection adjudication or Chae camp rebuttal as trigger conditions) is not autonomous-session output.
