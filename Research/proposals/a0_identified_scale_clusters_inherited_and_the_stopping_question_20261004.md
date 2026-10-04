# Proposal: the a₀ the fit identifies is a₀′/γ; clusters are inherited; and the stopping question

**Date:** 2026-10-04 · **Source:** site maintainer, from the 2026-10-04 visitor log (graduate-physics and
leading-edge-researcher personas) · **Bucket moves:** none · **Count:** 6 · **Bucket 0:** 0

## 1. Which a₀ does "a₀ = cH₀/2π" predict? (executed, small)

The SPARC fit solves g_bar = g_obs·tanh(γ ln(1 + g_obs/a₀′)) with a₀′ profiled. In the deep regime
tanh(γ ln(1+y)) → γy, so g_obs → √((a₀′/γ)·g_bar): **the Milgrom-equivalent scale is a₀_M = a₀′/γ**
(= 2a₀′ at γ = ½, which the site already states for that point only).

On the explorer's frozen-likelihood Υ_disk sweep (`synchronism-site/explorer/findings/scripts/sparc_gamma_interval_frozen_likelihood_output.txt` [3d], not re-fitted):

| Υ_disk | γ̂ | a₀′ | a₀′/γ̂ |
|---|---|---|---|
| 0.40 | 0.272 | 2.93e-11 | 1.074e-10 |
| 0.50 | 0.489 | 5.34e-11 | 1.091e-10 |
| 0.55 | 0.678 | 7.40e-11 | 1.092e-10 |
| 0.60 | 0.963 | 1.04e-10 | 1.084e-10 |

a₀′ moves 3.6× across the band; **a₀′/γ moves 1.6%.** Υ trades against γ (the shape), not against the scale.
The velocity-χ² fits of 2026-09-29 give a₀′/γ = 1.16e-10 (frozen nuisances) and 1.18e-10 (Υ, D, i profiled),
so the scale is estimator-dependent at ~8%. cH₀/2π (H₀ = 67.4) = 1.042e-10 sits **3–12% below** the fit's own
scale; the H₀ range 67.4–73 moves the prediction by ±5%. Controls: the deep-limit identity reproduced
numerically to 0.3% (next-order term) at γ = 0.5, 0.489, 2.
Script: `synchronism-site/maintainer/scripts/a0_identified_scale_across_upsilon.py` (+ `_output.txt`).

**Consequences.**
- The site's "13%" compares cH₀/2π with the conventional 1.2e-10, not with the framework's fit. Against the fit it
  is 3–12%: different, not sharper. Bucket 3 (Milgrom's 1983 coincidence) unchanged.
- **Correction to PREDICTIONS.md, 2026-08-14 block, item (iv).** "The a₀ 'profiled 5.33e-11 vs derived 1.04e-10,
  factor 1.96' tension dissolves at Υ = 0.6 … γ–a₀–Υ is one flat degeneracy, the shape parameter unidentified at
  factor 2." Both the tension and its dissolution compared a₀′ with Milgrom's a₀. The factor 1.96 is the 1/γ
  conversion at γ ≈ ½; the "dissolution" at Υ = 0.6 happened because γ̂ ≈ 0.96 ≈ 1 there. The scale is identified
  (to 1.6% across Υ within one estimator); the degeneracy is Υ ↔ γ. This agrees with the 08-04 finding that γ and
  a₀ enter only through γ/a₀ in the deep regime; it states which combination to compare with cH₀/2π.

## 2. Galaxy clusters: an inherited failure, not a new root (estimate, not an execution)

A visitor researcher persona asked whether clusters give a third ceiling kill. Round inputs (not re-read):
Coma-class M500 ~ 1e15 M☉, r500 ~ 1.3 Mpc, f_b(r500) = 0.13–0.15, so B_req = 6.7–7.7 at g_bar ≈ 1.1–1.2e-11 m/s².
With the registered TEST-09/10 form C_a = C_min + (1−C_min)·x/(1+x), x = (g_bar/1.05e-10)^(1/φ):

| floor | delivered B | shortfall |
|---|---|---|
| 1/Ω_m | 2.2 | ×3.1–3.5 |
| (Ω_m−Ω_b)/Ω_b | 2.8–2.9 | ×2.4–2.7 |
| Ω_m/Ω_b | 3.0–3.1 | ×2.2–2.5 |
| 1/Ω_b (no CDM) | 4.0–4.2 | ×1.7–1.8 |
| MOND simple ν | 3.7–3.9 | ×1.8–2.0 |

The framework falls short by at least MOND's known factor ~2, as the nesting (galaxies ⊂ MOND, 2026-09-27)
requires. Only the no-CDM floor is comparable to MOND. A variant keyed on g_bar/a₀ without the 1/φ power would
meet the cluster at the no-CDM floor, but that is the explicit wiring that sits 0.8–1.8 dex above the galaxy RAR.
So clusters sit inside the existing DM-fork decision (2026-09-22/25) and add no root.
Script: `synchronism-site/maintainer/scripts/cluster_r500_boost_vs_ceiling_estimate.py`.
For a real execution: the CLASH cluster RAR (Tian et al. 2020, ApJ 896, 70; acceleration scale ~2e-9 m/s², about
an order of magnitude above a₀; cited from memory, not re-read). Seeded as an explorer topic.

Also recorded on the site: **GW170817** (|c_T/c − 1| ≲ 1e-15) is an unaddressed constraint. There is no
specified relativistic completion, and the one radiative model (scalar inflow) is already refuted for lacking spin-2.

## 3. WAKE: the stopping question (frame, gates on dp)

The researcher persona ended with: *"Is any reading (density-keyed, acceleration-keyed, mean-density) expected to give
a different number from MOND+EFE or ΛCDM anywhere in reachable data? If not, the program has finished, and the
deliverable is the negative-results note."*

The archive has already answered most of this, in pieces:
- 2026-09-25: no test in TEST-01…26 has a branch that moves Bucket 0.
- 2026-09-27: every surviving construction is nested in its parent (galaxies ⊂ MOND; DE ⊃ ΛCDM, win branch is
  Cardassian prior art; wide binaries split; Cassini inherited).
- 2026-09-27 explorer: door #3's natural intrinsic-decoherence rate is excluded by the heat budget of ordinary matter
  by 15–30 orders. SPINE still calls door #3 "the one genuinely-open direction."
- Today: the cluster sector is inherited; the a₀ coincidence is no sharper against the fit than against 1.2e-10.

**What I think this means.** The *registered physics program* has met a stopping condition: no registered or
proposed observable has a branch where any reading of the framework differs from its parent theory in reachable
data. That is not "the ontology is wrong" (SPINE's single-observer move is untouched by any of this; B1/CHSH tested
only the local-realist readings). It is "the quantitative probes have finished teaching boundaries." Continuing
the daily site loop on the physics axis now mostly propagates qualifiers. Today's session spent most of its effort
on "+184" passages a 2026-09-29 sweep missed: 12 lines, now caught by a context-requiring lint rule.

**Recommendations (all gate on dp):**
1. Write the transferable-negatives note as the physics axis's deliverable: the local-density no-go with symmetron
   scoping, the unidentifiability lemma (γ/A only), the boost-ceiling kills (SPARC + KiDS lensing; the cap is the
   DM fork), and the A2ACW H1/H2 degeneracy. The researcher persona said they would cite the first two.
2. Update SPINE's door-#3 sentence with the heat-budget exclusion, or name the dissipative-soliton commitment that
   would reopen it. As written it points readers at a direction the explorer has closed at its natural value.
3. Split the headline count's texture, keeping the number: "4 data refutations (two share the ceiling root; one is
   inherited from MOND) + 1 at the threshold (γ = 2 pin, galaxy-level) + 1 construction check (CHSH)". Two personas
   independently asked for this today.
4. Decide whether the explorer track should keep running physics executions or pivot to the applied axis and the
   A2ACW post-cutoff control, which is the only experiment on the site whose outcome could change a headline in
   either direction.
