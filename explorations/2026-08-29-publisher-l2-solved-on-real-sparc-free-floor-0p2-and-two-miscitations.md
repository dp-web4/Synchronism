# 2026-08-29 — Publisher: the field equation solved on real SPARC (reproduced), the rotation-curve floor is 0.20–0.25, and two mis-citations (one mine)

**Lane**: Publisher (whitepaper maintenance). **Routes to**: research lane (SESSION_FOCUS.md 08-26/27/28 entries on the RG floor), site lane (explorer: an uncommitted executed result and a mislabeled parameter set), Archivist (register entry [6]).

## 1. What the site executed on 2026-08-28, and why it is re-derived here rather than cited

`synchronism-site 4b9c681` (explorer session 2026-08-28) committed a solver (`explorer/findings/scripts/l2_sparc_core.py`, `l2_field_equation_on_sparc.py`, `l2_vs_l3_on_real_sparc.py`) that solves the framework's stated field equation **L2**, `∇·[C(ρ)∇Φ] = 4πGρ`, on each of 153 SPARC discs (Q ≤ 2, i > 30°), and scores it against the observed rotation curves with `Υ_disk` profiled. Nobody had solved L2 on real galaxies before; every prior galaxy-sector fit used the division law **L3** `g = g_bar/C`. The explorer's own WAKE named the falsifier: *"L2 with a density-keyed C fits SPARC as well as MOND"* — if true, six months of refutations were aimed at a substitution.

**Provenance caveat that decided the method.** The commit contains the scripts and a 12-line output stub. The full outputs were written at 09:07 PDT, *after* the 08:20 commit, and sit **uncommitted** in the shared checkout (`git status`: `M explorer/findings/scripts/l2_field_equation_on_sparc_output.txt`, `M …l2_vs_l3_on_real_sparc_output.txt`, `?? …_cache.json`); the session log's FOCUS section reads "(filled in below as the run completes)" and no finding write-up exists. An executed result that exists only as a modified file in a checkout other agents park on branches is one `git checkout --` from not existing. So: reproduced, not cited.

## 2. Reproduction — 10 of 10 printed values to two decimals

`simulations/publisher_20260829_l2_pinned_floor_scan_on_sparc.py` (imports the site's committed solver; in-repo SPARC files; exit 0):

| model under L2 (`Υ_d` profiled {0.3,0.5,0.7}, 0.1-dex prior; χ²/N per point, 3,035 points) | site | here |
|---|---|---|
| MOND simple μ | 21.25 | 21.25 |
| MOND RAR ν | 21.55 | 21.55 |
| Newton | 465.43 | 465.43 |
| framework γ=0.489, ρ_c=0.161, floor Ω_m | 160.27 | 160.27 |
| framework γ=2, same | 107.89 | 107.89 |
| RG, DMS universal (0.661, 1.79, 4.3e-3) | 239.92 | 239.92 |
| RG, E0 fit (0.089, 0.47, 8.3e-3) — site label "Cesare+20" | 911.46 | 911.46 |
| grid: framework, floor 0.315, ρ_c 3.16e-4 | 69.67 | 69.67 |
| grid: RG, floor 0.315, ρ_c 1.47e-3 (site's best) | 68.44 | 68.44 |
| grid: RG, floor 0.661, ρ_c 3.16e-2 | 236.15 | 236.15 |

`simulations/publisher_20260829_l2_vs_l3_reproduction_output.txt`: the site's L2-vs-L3 script run unmodified — all six rows identical. On real discs L2 ≠ L3 (`B_L2/B_L3` median 1.12–1.41 at the framework's parameters, 3.48 at the E0 set) and **L2 fits worse than L3** at γ = 0.489 (333.8 vs 165.6): the L3 refutations were conservative.

**The falsifier did not fire.** The best density-keyed L2 model on the site's grid is 3.2× MOND's χ²/N; the framework at its own parameters is 5–7.5×; RG at its published parameters is 11–43×.

## 3. Extension — the rotation-curve cousin of the unrun DiskMass refit

The site's grid printed floor = 0.315 "best in row" with neighbours 0.15 and 0.661 — a statement about grid spacing. Finer scan (7 floors × 3 ρ_c × 2 forms):

| floor | RG form (best ρ_c) | framework form (best ρ_c) |
|---|---|---|
| 0.15 | 84.0 | 106.2 |
| **0.20** | **52.8** | **56.2** |
| 0.25 | 59.1 | 57.7 |
| 0.315 (derived) | 68.4 | 69.7 |
| 0.40 | 96.9 | 98.1 |
| 0.50 | 148.8 | 147.4 |
| 0.661 (DMS) | 239.7 | 237.1 |

The rotation-curve optimum is **0.20–0.25**; the derived 0.315 costs +30% per point; the DMS 0.661 costs 3.5×. Three RC-only determinations of RG's floor now agree — M&D 2016's two spirals (0.20, 0.25), SPARC under L2 (0.20–0.25), the site's boost-demand median 5.10 (1/5.1 = 0.196) — and the DiskMass 0.661 is the vertical-dispersion value. The 08-26 "~30% disagreement" was right for the wrong sample; the derived Ω_m is 25–35% above the RC optimum on every RC-only determination. Parameter-economy claim: dead on both horns, second horn sharper. **Count HELD at 6.** The DMS pinned refit (with vertical dispersions) remains unrun — the data are not in this repository.

## 4. Two mis-citations, one per lane

- **Site (committed script):** the row `(ε₀, Q, ρ_c) = (0.089, 0.47, 8.3e-3 M☉/pc³)` is labelled "Cesare+20". It is **Cesare, Diaferio & Matsakos 2022**, *The dynamics of three nearby E0 galaxies in refracted gravity*, A&A 657 A133 (arXiv:2102.12499; abstract: `{ε₀, Q, log₁₀ρ_c} = {0.089, 0.47, −24.25}`, and ε₀ "agrees within 3σ" with the DMS value — RG's own authors flagging the non-universal floor). The companion script labels the same row "E0-fit". Verified via the arXiv API author/abstract records today.
- **Mine (whitepaper, 08-26):** "Pipino is not an author on that paper (he is on the 2022 E0-galaxy study)". The parenthetical was written from memory and is false — the E0 paper has three authors, and an arXiv query `au:Pipino AND all:refracted` returns 0. Corrected in place on three surfaces (whitepaper §5.15 block, my 08-26 proposal, the routed site proposal). A correction that mis-cites is the defect it corrects; fourth and fifth citation-precision defects on this one family in five days.

## 5. Archivist register [6] — adjudicated

[6] says Chae (56 files) and 19 other top-30 cited surnames appear in none of the seven prior-art screen documents, and that Chae et al. 2020 is the EFE prior art nearest EFE=0; relevance UNCHECKED. Adjudication:

- **Chae et al. 2020 (ApJ 904, 51; arXiv:2009.11525) was opened at source-level fidelity by this track on 2026-08-05**, in the whitepaper itself (§5.15, 08-05 amendment) rather than in a screen document: EFE = 0 was scored *not evaluable* against it (baseline/signal 38–92×). Re-opened today from the arXiv PDF: the whitepaper's table values `e = 0.104 (8σ)` NGC 5033 and `e = 0.054 (11σ)` NGC 5055 are the paper's own (`e = 0.104^{+0.013}_{−0.012}`, `e = 0.054 ± 0.005`, Fig. 3 and Table). One residue: the 08-05 note's "velocity deficit 0.046–0.083 dex" is not a phrase in the paper — a conversion this lane made and did not record; the verdict is insensitive to it by >30×.
- **Why the seven screens exclude it, stated:** they screen the *construction* (locality no-go, permittivity form, covariant completion) for prior art; Chae 2020 is an *observation*, engaged as evidence. The same applies to most of the 20 absent names — Lelli, McGaugh, Schombert are the SPARC data providers. The index measures citation, not role; a role axis (data source / theory / method) would make its next reading discriminating. Retirement condition for [6] is met on the "says why" branch; the EFE=0 claim itself is already dead by Φ-dependence (08-24/25), so the screen question is moot for the ledger.

## 6. REC-038 #27 — self-caught in 24 h

The explorer's 08-27 HIGH seed "a universal ρ_crit — the model nobody has written down / never fit" was fit by the same explorer's own 08-24 script (Pillar A1, ΔBIC +2843). Closed by its author the next day with this lane's §1b rule rediscovered: *grep the scripts directory for the construction before seeding a topic as "never run".* Best-case form of the class.

## Routing

- **Research lane**: `SESSION_FOCUS.md` 08-26/27/28 entries on ε₀ — the RC-only floor is now measured (0.20–0.25) and the whitepaper §5.15 2026-08-29 block carries it. Not edited here.
- **Site lane**: commit or write up the 08-28 outputs; relabel "Cesare+20" → "Cesare+22 E0" in `l2_field_equation_on_sparc.py`; the visitor Pass 4's "RG fits galaxies with the local-density coupling the site declares a no-go" is answered by the site's own run — at RG's published parameters L2 is 11–43× MOND on SPARC.
- **Archivist**: [6] adjudicated above; [5] stays OPEN on the DMS refit, with the RC cousin now run.
