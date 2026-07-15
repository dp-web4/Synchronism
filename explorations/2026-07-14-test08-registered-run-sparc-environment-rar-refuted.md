# TEST-08 registered run: RAR residuals vs environment — REFUTED (r² = 0.0001, kill bar 0.09)

**Date**: 2026-07-14 (CBP; dp-gated run — dp said "yes let's run it")
**Registered in**: `Research/EXPERIMENTAL_TEST_CATALOG.md` TEST-08 (pre-registered method, expected value, kill criterion)
**Script**: `simulations/test08_sparc_environment_rar.py` · per-galaxy table `simulations/test08_per_galaxy_results.json`
**Status of the test going in**: pre-registered, runnable, **unrun** — the standing dp-gated lever the Publisher had held for six windows ("refuted ≠ untested; nobody has looked").

## The registered claim

> **Prediction**: RAR scatter correlates with galactic environment (cluster vs. field vs. void), not just internal structure.
> **Expected**: Environment explains >20% of RAR scatter.
> **Falsification**: Environment correlation < 0.3 (r² < 0.09).
> **Distinguishing power**: MEDIUM-HIGH — MOND predicts zero environment dependence.

Physics motivation (S177 lineage): coherence C depends on ambient density → lower environmental ρ → higher G_eff → void galaxies sit systematically high on the RAR.

## Pre-declaration (made before computing any correlation)

**Primary statistic**: Pearson r between per-galaxy RAR offset (mean dex residual from the pooled RAR fit) and distance-corrected log density from redshift-space cylinder counts in Cosmicflows-4 (projected radius 2 Mpc, |ΔV| < 500 km/s, self-excluded). r² < 0.09 → refuted; r² ≥ 0.20 → expectation met. Secondaries (context, not verdict): 5 Mpc sphere count, 5th-NN 3D density, Virgo-centric distance, tercile ANOVA, partial correlation controlling internal structure.

## Method

- **RAR side**: SPARC table1 + MassModels (in-repo). McGaugh+16 standard: M/L[3.6] = 0.5 disk / 0.7 bulge, Q ≤ 2, i ≥ 30°, e_V/V < 10%, g_bar > 0. Fit g† on pooled points; per-galaxy offset = mean dex residual (≥3 points).
- **Environment side**: SPARC has no coordinates → resolved all 175 names via NED (175/175, zero failures; `simulations/sparc_real_data/sparc_ned_coordinates.json`). Density tracer = **Cosmicflows-4** (in-repo `data/cf4/table2.dat`, 55,877 galaxies with Vcmb, DM-distance, RA/Dec). Metrics: redshift-space cylinder counts (2 Mpc, 500 km/s), 5 Mpc sphere counts, 5th-NN 3D density, Virgo-centric distance. CF4 is flux-limited → all count metrics distance-corrected (residual after regressing on log D); raw versions reported too.

## Pipeline validation (both instruments verified live before reading the answer)

- **RAR instrument**: 2696 points / 147 galaxies after cuts (published: 2693/153) ✓; fitted g† = 1.161×10⁻¹⁰ m/s² (published 1.20×10⁻¹⁰) ✓; pooled scatter 0.133 dex (published ~0.11) ✓. The pipeline reproduces McGaugh+16.
- **Environment instrument (agent-zero check — is the null a broken instrument?)**: dynamic range 0–162 neighbors, only 9/141 zeros ✓; three independent estimators inter-correlate (Spearman 0.86–0.90) ✓; **ground truth**: SPARC's 28 known Ursa Major Cluster members (distance-flag f_D=4) land at mean **74th percentile** of the density metric (random = 50th; UMa is a diffuse spiral cluster, so 74th not 95th is the right answer) ✓. **The instrument is alive; the null is a measurement, not a failure-to-measure.**

## Result

N = 141 galaxies with both RAR offset and environment.

| Statistic | r | r² | p |
|---|---|---|---|
| **PRIMARY: offset vs dist-corr log(1+N_cyl)** | **+0.012** | **0.0001** | 0.89 |
| 5 Mpc sphere count (dist-corr) | −0.128 | 0.016 | 0.13 |
| 5th-NN 3D density (dist-corr) | −0.141 | 0.020 | 0.09 |
| Virgo-centric distance | +0.108 | 0.012 | 0.20 |
| partial r (control log L36, log SBeff, log D) | +0.003 | 0.0000 | 0.97 |

Terciles (void/field/cluster): −0.030 / +0.009 / −0.020 dex — ANOVA F = 1.25, p = 0.29. No monotone trend, no significant differences.

**Registered verdict: REFUTED.** Primary r² = 0.0001, ~900× below the kill bar (0.09) and ~2000× below the registered expectation (0.20). Every secondary is also below the kill bar. Environment explains, at most, ~2% of per-galaxy RAR scatter — and the strongest secondary trends carry the **wrong sign** for the S177 mechanism (higher density → *lower* offset, where Synchronism predicted void galaxies high).

## Honest caveats (none rescue it)

1. **Attenuation**: per-galaxy offset scatter (0.124 dex) includes observational noise (distance, inclination errors); intrinsic scatter is ~half the variance. A true r² of 0.20 against intrinsic scatter would still appear as ~0.10 against total — above the kill bar. We measure 0.0001. Attenuation cannot bridge three orders of magnitude.
2. **CF4 distance errors** smear the 3D metrics — which is why the primary was pre-declared as the redshift-space cylinder (robust to this). All four estimators agree on the null.
3. **MOND scope**: the catalog's "MOND predicts zero environment dependence" is true only for MOND-without-EFE. MOND-with-EFE predicts a *small* environmental suppression (Chae+2020 claim a detection with a different estimator — external field strength, not number density). Our weak negative-sign secondaries (r ≈ −0.13/−0.14, ns) are directionally EFE-like, not Synchronism-like. Nothing here contradicts Chae; the registered Synchronism claim (>20%, void-high) is what died.
4. **Sample**: SPARC is not volume-complete and under-samples true cluster cores. But the UMa ground-truth check shows the density axis genuinely spans field-to-cluster, and 141 galaxies with a working instrument yielding r² = 10⁻⁴ is not a power problem: detecting the registered r² = 0.20 at this N would have p < 10⁻⁷.

## So what

- **The one runnable-and-unrun lever is now run.** Bucket 0 stays 0, but the board is now clean: no registered test sits untested-while-presented-as-dead or dead-while-presented-as-pending. TEST-08 moves catalog → Bucket 2 with an executed date and a named criterion that fired.
- This is the third independent kill in the galactic sector, and they now form a **convergent no-go**: locality no-go (door #1), ρ_crit(V) sign inversion, BTFR bounded-boost (TEST-09, executed earlier today by the explorer track), and now environment-independence of the RAR. Every "coherence depends on local/ambient density" observable that has been pushed to data has failed in the same direction: **the RAR behaves like a universal local law (MOND-like), not an environment-coupled one.** Per [interrogate-the-constant / convergent-ties-point-upstream], the constant across all four is the C(ρ) coupling itself — the framework's galactic-dynamics limb, not any single formula.
- S177/S179 closure: Session 179's proxy attempt (Dec 2025) saw *opposite-sign* trends with crude proxies and was ruled INCONCLUSIVE (proxies invalid). The registered run with a real external catalog (CF4) settles it: no signal at any usable estimator.

*Run executed as registered; result recorded either way, per the standing kill-criterion rule (dp decision 2026-07-08 lineage). Fetched assets (NED coordinates) committed for reproducibility.*
