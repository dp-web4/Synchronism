# The headline kill targets the wrong coherence function — and the first registered-parameter test of the density law forks it

**Date**: 2026-09-08 · **Origin**: synchronism-site maintainer track (WAKE phase), from the 2026-09-08
visitor log (graduate-physics + leading-edge-researcher personas) and the explorer's 2026-09-07 globular-cluster
execution · **Status**: PROPOSAL — two ledger actions gate on dp; one site-side correction already applied
**Builds on**: `argument_of_C_three_functions_ledger_not_commensurable_20260824.md` (the three-C register),
`galaxy_sector_internal_locality_ledger_20260907.md`, `Research/Session611_Stellar_Markov_Blankets.md` (P611.2)

---

## Why this is a research proposal and not a site fix

Two external-persona readers with **no archive access** reconstructed the three-C problem from the site's
own text in a single pass, and both landed on the same sentence: *the headline refutation (ΔBIC = +184) is a
refutation of an acceleration-keyed compander that the site elsewhere says was never the framework's
equation.* The archive already knows this (08-04, 08-24: "the fitted numbers come from C_g, NOT the headline
C_ρ"). What is new is that the mis-attribution is now **legible from outside**, which means the program's
public scoreboard is scoring the wrong object — and an outside referee can see it. That is a ledger defect, not
a wording defect.

The same day, the explorer executed **S611 P611.2** — the only registered test that keys on the framework's
*own* variable (ρ) at the framework's *own* parameter (γ = 2) on an object that physically crosses the knee
(Galactic globular clusters). It did not kill. It forked. So the density law's honest test set and the
acceleration proxy's test set have now diverged in *outcome*, not just in variable, and the ledger has no
column that says which is which.

## Verified this session (primary layer, not compilation layer)

1. **The RAR fit's mapping.** `simulations/sparc_tanhlog_profile.py:84-85`:
   `"""Invert g_bar = g_obs*tanh(gamma*ln(1+g_obs/a0)) vectorially."""` — C is used as an **implicit μ keyed
   on g_obs**, solved for g_obs. Not explicit on g_bar (a ν-slot reading), and not on ρ. Consistent with the
   profiled a₀′ = 5.33×10⁻¹¹ sitting 2.11× below McGaugh's reference: μ_simple(x/2) at γ = ½ is MOND with
   a₀ = 2a₀′, exactly the factor the two personas back-derived. The explicit-on-g_bar reading the graduate
   persona computed (g_obs → g_bar + 2a₀′, a constant additive floor, 20× off at 10⁻¹² m/s²) was **never run**;
   the site's sentence "MOND with μ's argument swapped from g_bar to ρ" was wrong twice (μ's argument in MOND
   is g_obs; the fit swapped nothing). **Corrected on the site today** (`/coherence-function`,
   `/for-researchers`, `/galaxy-rotation`). The archive's own wording (08-24: "C_g, keyed on g_obs, implicit")
   was already right — this was a site-only drift.

2. **What +184 refutes.** The form-selection run (07-22) reproduced +184 "as pipeline sanity check" on the
   *same* implicit-g_obs pipeline. So +184 = C_g at a pinned γ = 2. It is not a test of C_ρ at any γ.

3. **What refutes C_ρ.** (a) Head-to-head on SPARC, free γ: ΔBIC +2843, γ → 0.046 (08-24); (b) the boost
   ceiling, for the floored form (TEST-09/10, γ-independent); (c) the Oort limit at the solar midplane
   (explorer 09-06: f_DM = 0.685 predicted vs 0.13 ± 0.04 for any knee above ~0.15 M☉/pc³; a window
   0.074–0.154 survives at γ = 2); (d) **globular clusters (explorer 09-07)**: at universal γ = 0.489 every
   knee the framework uses is excluded (ρ_c ∈ 0.1–300 M☉/pc³, 3.7–4.4× the Newtonian residual); at the
   **registered** γ = 2 the measured knee 0.161 is *marginal* at MOND+EFE's level, and the Oort window
   survives. The discriminating variable, measured: MOND survives clusters *because of* the EFE (EFE off:
   −0.245 ≈ density law); a density-keyed law has no g_ext.

## Proposal

### A. A "target" column on Bucket 2 (gates on dp — ledger structure)

Every Bucket-2 galaxy row names which coherence function it kills: **C_ρ / C_g / C_Ω / other**. Under that
column the current rows read: TEST-09, TEST-10 → C_Ω; RAR +184 → C_g (γ pinned); Cassini → C_g;
environment null → registration (S177), consistent with C_ρ; ρ_crit exponent → C_ρ. The 08-24 proposal
says "largest coherent sub-ledger = 2"; the column makes that visible per row instead of in a footnote.
The site's headline scoreboard will follow the ledger: the RAR-shape entry now reads "γ = 2 pin,
acceleration-keyed realization" and the density-keyed kill is led by the ceiling and the head-to-head.

### B. P611.2: registered → executed, outcome FORK (gates on dp — count)

Recommendation, agreeing with the explorer: **do not increment the count.** A registered prediction
(γ = 2 inside a resolved-member system) survives the test it was registered for. What is refuted is the
*conjunction* {universal γ = 0.489} ∧ {any knee the framework uses}. Two honest readings, both true:

- If γ is universal (the "one equation" reading, and the value SPARC + DESI both select): globular clusters
  exclude the density law at every placement. Refutation-grade, but of a conjunction the ledger never
  registered as one row.
- If γ resets per Markov blanket (P611.2, on the record since 2026-02-17): consistent, marginal — and the
  coherence function is not one function.

Either way the count's meaning changes more than its value. Recorded in `Session611_Stellar_Markov_Blankets.md`
as an execution note and flagged in PREDICTIONS.md Bucket 1 (no bucket moved).

### C. The γ ladder as a candidate Bucket-1 registration (gates on dp)

P611.2 implies a **ladder**: γ = 2 for every resolved-member system (open clusters, dwarf spheroidals,
Gaia streams), γ ≈ ½ for every unresolved one (disks in the RAR). This is the first time the N_corr /
Markov-blanket machinery makes a prediction that a public dataset can check and that differs from
"γ is whatever the fit says."

Draft row: **B8 — Markov-blanket γ reset.** *Refuted if* pressure-supported resolved-member systems
(Galactic dSphs via Walker/Battaglia dispersions; open clusters via Gaia) fit under the same knee
(0.05–0.16 M☉/pc³, floored) demand γ ≈ ½ rather than γ = 2 at ≥ the GC exclusion contrast; *reparametrization
if* γ = 2 and γ = ½ both pass because the knee sits outside both systems' sampled densities (the
identifiability trap again — pre-check the window first). *Not a win if* it passes: the ceiling still kills
the floored law on SPARC independent of γ, so a passing ladder makes the framework *two functions, both
partial*, not one function that works.

### D. The frame point (no decision needed)

This is the **fourth instance in eight days** of "which parameter / which variable does the criterion name":
07-29 estimator, 08-08 coupling, 08-27 what R_half is a function of, 09-07 γ vs knee. Three of the four were
caught *before* a false refutation was published; one (the +184 headline) was published in May and stood for
four months. The 09-07 maintainer log counted seven over-refutations in one direction; the mis-targeted
headline is not an eighth over-refutation (γ = 2 in C_g *is* dead) — it is a **mis-addressed** one, which is
the same failure mode wearing the other sign. The fix is structural (the target column), not vigilance.

And the reason the outside personas could find it is the reason the site exists: the record is legible enough
to be audited by a stranger. That is working.

## Requested from dp

1. Adopt/decline the Bucket-2 target column (A).
2. Rule on P611.2's outcome: fork, count unchanged (B) — or increment.
3. Adopt/decline B8 (C), or route it to the explorer for the pre-check first.
