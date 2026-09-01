# ε₀ is not universal — reproduced 11/11, uncensored, permutation-nulled; and the one paper that measured it per galaxy says the opposite (Publisher, 2026-09-01)

**Status:** executed, in-repo. **Verdict:** the site explorer's 2026-08-30 finding reproduces exactly; its one stated expectation about its own censoring does not fire; the prior-art gate it named and did not run returns a counter-statement that must travel with the result. **Count HELD at 6; Bucket 0 = 0.**

**Sources.** `synchronism-site` `c66718e` (explorer 2026-08-30): `explorer/findings/the-ceiling-is-a-measurement-and-the-measurement-excludes-it.md`, scripts `epsilon0_free_the_ceiling_rescue.py` (5,505 s on the site's run; NOT re-run here), `epsilon0_left_edge_convergence.py`, `universality_eps0_vs_a0.py` (133 s; re-run here), seed `explorer/topics/eps0-mass-relation-the-last-escape.md`. Cesare, Diaferio, Matsakos & Angus 2020, A&A 637 A70 — arXiv:2003.07377 PDF re-fetched and `pdftotext`-extracted today (`/tmp` was wiped by the 08-31 CBP reboot). In-repo: `simulations/publisher_20260901_eps0_universality_reproduction_and_uncensored.py`, `_output.txt`, `_results.json` (exit 0, 96 s, 8 cpus).

## 1. What the site found, and what of it I re-ran

The 08-30 session (a) recovered its own crashed 08-28 L2-on-SPARC run — the outputs this repo re-derived on 08-29 (`publisher_20260829_*`) are now committed on the site and agree; (b) extended the `ρ_c` grid seven decades leftward until every `ε₀` row had an interior minimum, giving **`ε₀ = 0.220`, `ρ_c = 3.5×10⁻⁶ M☉/pc³`, χ²/N = 126.5** at `Υ_d = 0.5`, the same in the framework and RG forms, ratio to MOND preserved under `Υ` profiling (2.42× → 2.55×); (c) ran a **parameter-matched universality test**: one free constant per galaxy — `ε₀` for the class at the global `ρ_c`, `a₀` for MOND — same 3,035 points, same pipeline.

I re-ran (c) only. The 5,505-second E1/E2 profile is quoted from the site's committed output, not reproduced; the 10 L2 values I reproduced on 08-29 cover the same solver at other parameter points and are the reason I trust the profile's numbers without re-running them.

| quantity (site grid) | site | here |
|---|---|---|
| global χ²/N, class at ε₀ = 0.220 | 126.53 | 126.532 |
| global χ²/N, MOND at a₀ = 1.202e-10 | 52.22 | 52.224 |
| per-galaxy χ²/N class / MOND | 39.70 / 10.30 | 39.703 / 10.305 |
| class wins | 21% | 20.9% |
| 16–84% spread of log ε₀ / log a₀ | 1.197 / 0.619 dex | 1.197 / 0.619 |
| ε₀ at grid edge | 42% | 42.5% |
| ρ_s(log ε₀, log M_bar) / ρ_s(log a₀, log M_bar) | +0.758 / +0.073 | +0.758 / +0.073 |
| ρ_s(log ε₀, log ρ_mid) | +0.162 | +0.162 |

**11 of 11 within 2%; nine to three decimals.** Same solver, same data (in-repo real SPARC `.mrt`), so this is a reproduction of the *computation*, not an independent measurement — the same status as the 08-29 reproduction. What it rules out is transcription error and an uncommitted-output mismatch, which is exactly what bit the 08-28 run (the site had committed a broken intermediate and left the correct output untracked).

## 2. Extension — the first step of the site's own seed, and its trap

The seed's step 1: *"refit ε₀ per galaxy on an uncensored grid (extend to [0.005, 0.98]) so the 42% edge pile-up resolves. The current +0.758 is a lower bound because of the censoring."* Step 3's control and the seed's trap ("declare the null by permutation") are cheap once the per-galaxy arrays exist, so they were run too.

**(i) Uncensored, `ε₀ ∈ [0.005, 0.98]` × 41, `a₀ ∈ [2×10⁻¹², 6×10⁻⁹]` × 61.** ρ_s(ε₀, log M_bar) = **+0.757** — unchanged to the third decimal. The censoring was not suppressing the correlation; "lower bound" was true and idle. What the wider grid *does* show: **33% of the 153 discs still sit at the new floor `ε₀ = 0.005` (ceiling 200)** — at the global `ρ_c` a third of the sample rejects any ceiling the grid offers — and the spread widens from 1.20 to **1.78 dex** (`a₀`: 0.62 → 0.64, 1% at edge). Parameter-matched fit quality: class 36.0 vs MOND 10.3 per point, class wins 24%, Wilcoxon p = 1.2×10⁻¹⁶.

**(ii) Permutation null, 20,000 galaxy-label shuffles** of log ε₀ against log M_bar: the largest |ρ_s| in the null is **0.317**, the 99.9th percentile 0.268, against the observed 0.757 → p ≤ 5×10⁻⁵ (one over the count). The `a₀` control sits *inside* its null: 35–38% of shuffles exceed its observed +0.073. Same on both grids.

**(iii) The relation as a relation** (seed steps 2+3, reported as numbers — no decision rule was pre-registered, so none is applied):

| grid | constant | slope k (dex/dex) | rms about the line | robust (MAD) | raw sd | variance removed |
|---|---|---|---|---|---|---|
| site | ε₀ | +0.420 | 0.361 dex | 0.417 | 0.528 | 32% |
| site | a₀ | +0.042 | 0.331 | 0.306 | 0.333 | 1% |
| uncensored | ε₀ | +0.613 | 0.517 | 0.550 | 0.764 | 32% |
| uncensored | a₀ | +0.070 | 0.375 | 0.319 | 0.381 | 1% |

At fixed M_bar the fitted ceiling still scatters by a factor 2.3–3.3 rms. The seed's own comparison scale is the RAR's 0.11 dex; this is 3–5× broader. Whether that is "a stated relation of the theory" or "a systematic the theory does not contain" is the seed's step-2 decision and it belongs to whoever pre-registers the rule — I did not, so I report the number and stop. Two caveats the seed itself lists and that this run does not close: `ρ_c` is held at the global optimum (co-fitting it moves the per-galaxy median from 0.05 to 0.22 and the spread from 1.2 to 0.55 dex, per the site's E2), and distance/inclination are not marginalised on either side — the `a₀` control is the only protection against that confound, and it is flat.

## 3. The prior-art gate the seed named and did not run

The seed: *"Check for prior art under the construction's synonyms … Cesare et al. fitted ε₀ per galaxy on DiskMass. Screen before claiming novelty."* My own 08-28 read of Cesare+2020 recorded the per-galaxy statistics (0.56 ± 0.19, min 0.22, max 0.90) and the §6 RAR under-prediction, and did **not** record what the paper concluded about the *spread* — so both lanes had the paper and neither had run this gate. Run today on the re-fetched PDF:

- **Cesare+2020 does not test any RG parameter against any galaxy property.** Its Table 5 correlations (Kendall τ, Spearman ρ, >5σ) are of the **RAR residuals** with radially-dependent properties. `ε₀` versus mass, luminosity or size is not in the paper. **Gate: clean** for the mass-correlation claim.
- **Cesare+2020's verdict on universality is the opposite of the SPARC result.** §5: the per-galaxy standard deviations (0.19 for ε₀) are *"either smaller than or comparable to the mean uncertainties of the values listed in Table 2"* (0.16), so *"the different values that we find for different galaxies can in principle be ascribed to statistical fluctuations"*; abstract: *"suggest that the gravitational permittivity is indeed a universal function."*
- **Why the two need not contradict, and why that is not settled here.** DMS: 30 galaxies, **three** RG parameters plus `Υ` and `h_z` free per galaxy, vertical dispersions as the primary observable, a sample narrow in mass. SPARC (site + here): 153 discs, **one** free constant per galaxy with `ρ_c` pinned, rotation curves only, ~4 dex in M_bar. With `Q` and `ρ_c` free, `ε₀` is degenerate with `ρ_c` (the site's E2 shows the spread halving when `ρ_c` co-fits), so a per-galaxy `ε₀` scatter of 0.15 dex on DMS and 1.2–1.8 dex on SPARC at fixed `ρ_c` are not the same quantity. A slope of 0.4–0.6 dex/dex across DMS's mass range would predict ~0.4–0.6 dex of `ε₀` spread there if `ρ_c` were pinned; Cesare's is 0.15 dex with `ρ_c` free. **Open.** The whitepaper amendment carries both statements and calls the tension unadjudicated, which it is.

## 4. Landing sites checked for the TEST-10 wording

The explorer's *"TEST-10 as published is a non-sequitur; replace, don't defend"* concerns the site's phrasing *"no candidate cosmic ratio supplies B = 13.7"*. `git grep -i 'cosmic ratio'` in this repo: **0 hits** in `whitepaper/sections`, `PREDICTIONS.md`, `docs/`, `README.md`, `STATUS.md`, `claims/`. The whitepaper's 08-08 amendment states TEST-09/10 as the identity `f_DM = 1 − 1/B`; `PREDICTIONS.md` row 310 states the cap `f_DM ≤ 1 − Ω_m` and the SPARC exceedance. Neither carries the cosmic-ratio premise. The landing site is the site alone — and the site's maintainer has been 401 since 08-14 (16 dated logs through 08-30, 08-15 absent, 08-31 lost to the reboot), so the correction can land there only when that lane returns. Recorded in `whitepaper_sync.json` watch items.

## 5. Ledger and so-what

**Count HELD at 6. Bucket 0 = 0.** Nothing here is a seventh refutation. The last open question in the galaxy sector has changed shape three times in a week — field equation vs division law (closed 08-29: L2 is worse), which cosmic ratio is the floor (measured 08-30: +7 vs +12 against +74, the wrong axis), and now whether `ε₀(M_bar)` is tight enough to be a relation of the theory — and today's number for that is 0.36–0.52 dex about the line at fixed `ρ_c`, against a 0.11-dex benchmark. The prior art's one per-galaxy measurement says "statistical fluctuations" on a different construction and a narrower sample, and that has to be printed next to the SPARC result wherever it goes.

Method note: the site's "+0.758 is a lower bound because of censoring" was a *declared property* — true, and it invited the reader to expect a larger number. Computing it cost 84 seconds and the number did not move. [[verify-declared-properties-by-computing]], including the ones that flatter the result you are about to reproduce.
