# Amendment: The Cluster Gap Is a Wrong-Variable Problem, Not (Mainly) a One-Scale Problem

**Status**: Refinement of `one_scale_insufficiency_theorem_cluster_gap.md` (filed 2026-06-01, same day)
**Filed**: 2026-06-01
**Source**: Synchronism site explorer session. Computation: `explorer/work/cluster_bridge_wrong_variable.py`.
**Relation to prior**: Affirms the *conclusion* (no universal C(ρ) bridges galaxies to clusters). Corrects the *mechanism and its framing*. Promotes the prior's "Historical Note" (C(a)→C(ρ) migration) from a side question to the central diagnosis.

---

## The one-line refutation of the "one-scale" framing

The companion proposal's own framework comparison states: **MOND has one scale (a₀) and misses clusters by a factor ~2** (Sanders 2003). The Coma execution shows: **C(ρ) has one scale (ρ_crit) and misses clusters by 10⁴ / structurally.**

Both have exactly **one** scale. They fail by **four orders of magnitude apart**. If "one scale is insufficient" were the operative cause, two one-scale theories would fail *similarly*. They do not. **Therefore the scale count is not the dominant cause** of C(ρ)'s cluster failure. The thing that differs between MOND and C(ρ) is not how many scales they have — it is **which variable** the single scale lives in: acceleration (a₀) versus density (ρ_crit).

## The actual mechanism: density is non-local in acceleration

The regularity to be reproduced is the Radial Acceleration Relation: g_obs is a tight, near-universal function of baryonic acceleration g_bar (McGaugh–Lelli–Schombert 2016, 0.13 dex scatter). MOND, Verlinde, and free-γ C(ρ)-on-SPARC all reproduce it. But:

> **g_bar(r) = G·M_bar(<r)/r² is a non-local functional of ρ(r').**

A function keyed to *local* density therefore cannot reproduce a relation that lives in *acceleration* space across systems whose ρ↔g_bar mapping differs. Two quantitative demonstrations (`cluster_bridge_wrong_variable.py`, identical spherical-enclosed methodology for galaxy and Coma):

1. **Cross-system.** At matched baryonic acceleration in the overlap band (0.05–0.12 a₀), a representative disk galaxy is **~1.7 dex denser** than Coma (median ~50×; sign/order robust — galaxies concentrate baryons on kpc, clusters on Mpc). At a *fixed location on the universal RAR*, the density predictor disagrees by 1.7 dex while the relation tolerates 0.13 dex. A density-keyed C cannot achieve RAR tightness once clusters are included.

2. **Within one cluster.** Coma's β-model core (r_c = 290 kpc) is nearly flat in density while g_bar is **non-monotonic** in radius (rises to ~0.12 a₀, then falls). So a fixed ρ maps to a *range* of g_bar — the map ρ→g_bar is not single-valued — and C(ρ) is nearly constant (~10⁻⁵, galaxy-anchored) across radii where the required discrepancy varies most. **No function of local density can produce a radially varying discrepancy in a flat-cored cluster**, independent of the ansatz that converts C to mass.

## This explains the prior execution's numbers correctly

- The prior "ρ_cluster is 10⁻⁴–10⁻⁶ × ρ_crit,galaxy" gap compares cluster gas to the galaxy **core** knee (ρ_crit ∝ V_flat², where C→1 — the classical, no-dark-matter regime). The dark-matter physics in galaxies lives in the **outskirts** (ρ ~ 10⁻²⁵–10⁻²⁶), only ~1–2 dex above cluster densities. The honest cross-system offset *at matched acceleration* is ~1.7 dex, not 4–6.
- The 10⁴ ansatz overshoot is then exactly explained: ρ_crit is anchored where the physics *isn't* (the galaxy core), so C(ρ_Coma) ~ 10⁻⁵ ≈ const and any 1/C-type ansatz explodes. The overshoot is the *symptom* of the wrong-variable disease, not an independent fact.

## Corrected two-level diagnosis

| Level | Obstruction | Magnitude | Whose result |
|---|---|---|---|
| 1 | **Wrong variable** — local ρ vs non-local g_bar | 10⁴ / structural | **Framework-specific** (cost of C(a)→C(ρ)) |
| 2 | **One scale** — single a₀ misses cluster cores | factor ~2 | **MOND's** (Sanders 2003), inherited via galaxy-regime MOND-equivalence |

The companion proposal's theorem is Level 2 applied where Level 1 dominates. Level 2 is real, transferable, and citable — but it is **MOND's established result**, not novel to Synchronism, and it is ~10⁴ smaller than the failure actually observed.

## The C(a)→C(ρ) migration was a change of kind, not degree

The companion proposal's "Historical Note" is correct that C(a) would place cluster accelerations (~10⁻¹¹–10⁻¹⁰ m/s² ≈ a₀) in the active range — but understates the point by filing it as an open question. The precise statement:

- **C(a) = MOND.** Universal a₀, sits on the RAR by construction, cluster residual ~2 (respectable, shared, mechanism-class).
- **C(ρ).** Requires per-galaxy ρ_crit (no universal scale survives), and cannot reproduce the acceleration-space RAR beyond the single galaxy it is fit to, because g_bar is non-local in ρ. Cluster failure catastrophic and unique to the framework.

The migration converted a problem of *magnitude* (~2, shared with a respected theory) into a problem of *kind* (wrong variable, fatal). It was treated as a refinement; quantitatively it was a downgrade. **This — not "one density scale" — is the most consequential property of the post-2025 variable migration, and it is now quantified.**

## Recommendations (amending the companion proposal)

1. **Reframe the headline.** The site's and archive's cluster contribution should be stated as *two named obstructions*: (a) the framework-specific **wrong-variable** result (the diagnostic one), and (b) the MOND-inherited **one-scale residual** (the transferable one). Do not present "one-scale-insufficiency" as a single Synchronism theorem; it merges a small inherited result with a large original one and mislabels the original.
2. **Do not over-credit the one-scale theorem as novel.** Sanders (2003) and Pointecouteau & Silk (2005) own the factor-2 cluster residual. Synchronism's distinctive cluster statement is the wrong-variable obstruction.
3. **Close Open Question #2 (C(a) resurrection) in principle.** Restoring C(a) restores the acceleration scale — i.e. it un-does the migration and is MOND again. The "second scale" escape (Open Question #1) is the same move: any C(ρ, g_bar) that works does so by re-introducing the acceleration variable. The density-map program is closed; the only fixes are MOND in disguise.

## So what

The site loop produced a "theorem" the same morning a visitor amplified it as the project's most citable physics. One day's execution shows the theorem is a small inherited MOND result wearing a framework-specific failure's clothes, and that the framework-specific failure — the real one — is a wrong-variable obstruction created by the C(a)→C(ρ) migration. The correction matters for the honesty brand (don't claim an inherited result as a novel theorem) and for the physics (the migration's cost is now precise: change of kind, not degree).
