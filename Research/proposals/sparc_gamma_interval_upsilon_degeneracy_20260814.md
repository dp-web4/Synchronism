# γ_SPARC = 0.49 ± 0.11 (stat), ϒ-band [0.27, 0.96] (syst) — the cross-sector γ-concordance is convention-dependent and permanently unpowered

**Date**: 2026-08-14
**Source**: synchronism-site explorer track, self-directed (open thread of the
2026-08-12 direct-fit proposal `gamma_family_direct_fit_desi_dr2_20260812.md`).
**Script**: `synchronism-site/explorer/findings/scripts/sparc_gamma_interval_frozen_likelihood.py`
(runs against this repo's `simulations/sparc_real_data/MassModels_Lelli2016c.mrt`,
reproducing the frozen tanh-log likelihood of `simulations/sparc_tanhlog_profile.py`).

## Result

The galaxy-sector γ = 0.489 (quoted since the 07-22 compander adjudication, and one
leg of the 08-12 "first executed cross-sector consistency test: PASS at 0.1σ") has
never carried an uncertainty. Derived on the frozen likelihood itself:

| Rung | Estimator | σ(γ) or band |
|---|---|---|
| naive | 2807 points independent, Δχ²=1 | 0.029 |
| galaxy-level | jackknife / bootstrap(400, wide bounds) / √(N/N_gal) | 0.103 / **0.113** / 0.121 |
| ϒ systematic | γ̂ refit at ϒ_disk ∈ {0.4, 0.5, 0.55, 0.6} | **γ̂ ∈ [0.27, 0.96]**, rms flat (0.1433–0.1458 dex) |

- Continuous free minimum γ̂ = 0.4892 (grid 0.489 was exact); P(γ̂* ≥ 1/2) = 0.49.
- Single galaxies (UGC11914, UGC03580, NGC5985) move γ̂ by ±0.03–0.04 each — ~3× the
  entire |0.489 − 0.5| offset. The three decimals of "0.489" are noise digits.
- The likelihood's own mild preference is ϒ_disk = 0.55 → γ̂ = 0.68 (Δχ²_naive = −14.7
  *better* than the ϒ = 0.5 convention that produced 0.489; galaxy-scaled ~0.9σ).

## Consequences

1. **Concordance re-priced**: |γ_cosmo 0.487 − γ_SPARC 0.4892| = 0.002 vs combined
   σ = 0.115 → 0.02σ, 51× short of 1σ-separation power. Under ϒ = 0.55 the same
   comparison reads 0.68 vs 0.487 ≈ 1.7σ "tension" with lower χ². The 08-12 PASS is
   a property of the ϒ_disk = 0.5 convention. (First over-affirmation instance of the
   unpropagated-nuisance class; prior five instances were over-refutations.)
2. **Permanent closure**: reaching σ_γ = 0.004 (the 08-12 bar for separating 0.489
   from 0.5) needs ~130,000 SPARC-quality galaxies AND global ϒ known to ±0.0007 —
   the ϒ term does not average down. γ=1/2-vs-fit discrimination on rotation curves
   is closed a priori.
3. **Toward the framework**: the "profiled a₀ = 5.33×10⁻¹¹ vs derived cH₀/2π =
   1.04×10⁻¹⁰, factor 1.96" tension dissolves inside the same valley — at ϒ = 0.6
   the profiled a₀ is 1.043×10⁻¹⁰ ≡ the derived value, indistinguishable rms. The
   γ–a₀–ϒ triple is one flat degeneracy; the shape parameter is unidentified at
   factor 2.
4. **Downstream**: the q = 2γ ≈ 0.98 ⇒ simple-μ identification (07-22) is an
   ϒ = 0.5 slice statement (q spans [0.5, 1.9] under the band); any sector consuming
   "γ ≈ 0.5 from galaxies" inherits factor-2 slack. TEST-12-class registrations must
   state the a-priori closure.

## Open thread

Hierarchical refit with per-galaxy ϒ priors (Li et al. 2018 RAR methodology) is the
defensible next rung before any external citation of the interval.

## Nuisance accounting

PROFILED: a₀. FIXED: ϒ (swept globally, ϒ_bulge = 1.4×ϒ_disk), distances,
inclinations. Bootstrap absorbs per-galaxy coherent scatter empirically; frozen
unweighted log-residual objective throughout; acceleration-keyed (g_bar) — the
substituted model the site actually quotes, not the ρ-keyed masthead equation.
