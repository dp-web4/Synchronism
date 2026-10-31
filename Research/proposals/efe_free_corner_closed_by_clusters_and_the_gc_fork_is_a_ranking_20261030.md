# The external-field-free corner is closed by globular clusters; the GC fork is a ranking; the tanh's shape touches data only as exclusions

*Site maintainer, 2026-10-30. From simulated-visitor graduate and researcher personas (site visitor log 2026-10-30,
Passes 3–4). No bucket moves. Count 6, range 3–8, Bucket 0 = 0.*

> **⚠ CORRECTION (2026-10-31, site maintainer, the same author; from site explorer 2026-10-30 and a simulated-visitor
> researcher persona 2026-10-31, independently).** §1's table is a **data-free statistic.** For two laws the paired
> per-cluster difference is (s_obs − s_A) − (s_obs − s_B) = s_B − s_A: the observed slope cancels exactly, and
> shuffling s_obs across clusters reproduces every difference and its bootstrap width to the printed digit
> (`synchronism-site/explorer/findings/gc-paired-differences-cancel-the-data-the-8-and-10p6-sigma-are-model-separations.md`).
> The 8σ / 14σ / 6.6σ / 4.4σ / 10.6σ / 10.8σ measure how far apart the models are, not which one the data prefer. The six
> "held" predictions held vacuously. Data-facing re-score (explorer, pre-registered site `012ca7c`; per-law refit of r_t,
> β, M with a 0.25 dex photometric-r_t prior, 28 clusters, cluster bootstrap):
> - **(a) is withdrawn as worded.** The EFE-free corner is **disfavoured at ≈ 3σ, conditional on the mass model, not
>   closed**: compander γ = 0.9 EFE-off Δχ² = +54 ± 19 vs Newton (+93 ± 32 with β ≤ 0); simple μ EFE-off +71 ± 23.
>   Still inherited (only a₀′/γ enters).
> - γ = 2 (registered) is **indistinguishable from Newton** (Δχ² +0.8 ± 2.8) and ahead of MOND+EFE by 1.5σ. "A ranking"
>   stands, re-sourced to Δχ².
> - **Newton vs MOND+EFE has no stable sign**: catalogue inputs MOND+EFE better (−503 ± 186, 2.7σ); photometric r_t
>   Newton better (+914 ± 389, 2.3σ); refit Newton better (+8.7 ± 4.7, 1.8σ). "Every MG row ≥ 4σ worse than Newton" is
>   withdrawn. The −0.057 "systematics floor" is the Webb+2013 eq. 8 tidal radius (computed); Harris r_t gives +0.015 ± 0.031.
> - §3's exclusion list: read "GC: EFE-free corner closed" as "≈ 3σ disfavoured, inherited from MOND-no-EFE". The
>   frame question's conclusion (the tanh's shape is not exercised there) is unchanged; the negative is weaker.
> - §"For the project lead" item claiming the corner "closed": same correction. The text below is kept as written.

## 1. Executed today (pre-registered at site `fcc4809` before the script)

`synchronism-site/maintainer/scripts/gc_row_differences_cluster_bootstrap.py` (+ `_PREREG.md`, `_output.txt`).
It uses the 42-cluster outer-slope statistic with the same selection, Jeans solver and King mass models as the
2026-09-07 table. The identity control reproduces every published row to ±0.001. Rows are scored on **paired
per-cluster differences**, with a cluster bootstrap (4000 draws).

| pair (A − B) | difference | σ_boot | |
|---|---|---|---|
| Newton − MOND+EFE | +0.037 | 0.005 | 8.0σ (MOND+EFE worse) |
| MOND EFE-off − MOND+EFE | −0.151 | 0.011 | 14σ |
| density γ = 0.489 − MOND+EFE | −0.118 | 0.018 | 6.6σ |
| **density γ = 2 (registered) − MOND+EFE** | **−0.018** | **0.010** | **1.7σ, not separated** |
| density γ = 2 − Newton | −0.055 | 0.012 | 4.4σ |
| **compander, EFE-free, γ = 0.87–0.95 − MOND+EFE** | **−0.147** | **0.014** | **10.6σ** |
| compander, EFE-free, γ = 0.87–0.95 − Newton | −0.184 | 0.017 | 10.8σ |

All six predictions held. (P6, Newton separated from MOND+EFE, I had put at about 50%.)

**(a) The EFE-free corner is closed.** The 2026-10-16 block left γ ≈ 0.87–0.95, with no external-field effect, as
the one point surviving SPARC (galaxy level) and the EPM2011 perihelia (γ < 0.87 excluded). On the cluster statistic it
gives −0.241 at every γ in that range, and −0.240 at γ = ½. Cluster outskirts sit in the deep regime, where only a₀′/γ
enters, so **the statistic cannot see γ**. The closure is therefore the exclusion of MOND without an external field,
inherited, and not a test of the compander's shape. Caveats as for the table: L2/algebraic, isolated clusters,
Newton-fitted mass models, no escapers or anisotropy.

With this, the 2026-10-21 sentence ("no galaxy-scale prediction that differs from MOND exists on a surviving branch")
loses its last exception. The only γ ≠ ½ reading left is the external-field-on γ ≳ 1.5, which passes Q₂. It is
MOND-class with a sharper return, and SPARC disfavours it at about 2σ at galaxy level.

**(b) The GC fork is a ranking.** The coded rule (pass if |residual| ≤ |MOND+EFE|) makes MOND+EFE pass by
construction. On paired differences, the registered γ = 2 point and MOND+EFE are not separated (1.7σ). "Marginal"
puts the registered point next to MOND+EFE; it does not distinguish the two. Every modified-gravity row, MOND+EFE
included, is worse than Newton by ≥ 4σ. Part of Newton's lead is built in, because the mass models were fitted under
Newton (the 09-19 mass-scale refit moved the ratios by 7%). Read the 09-08 "registered prediction survives" as "is
ranked beside MOND+EFE under L2"; under L3 it is excluded (10-08).

## 2. Number hygiene (re-run at 4 decimals today)

- **TEST-09 sweep:** the slopes are 3.3465 / 3.4557 / 3.4923, against data at 3.7527 ± 0.1001, so |Δn| = 0.406 / **0.297** / 0.260.
  At (Ω_m−Ω_b)/Ω_b the slope misses the strict > 0.3 line by 0.003, about 0.03σ. The |Δn| = 0.300 edge falls at
  **B_max ≈ 5.31**, not ≲ 5.4. Verdicts are unchanged.
- **σ-distances (10-16 block (ii)):** 3.4σ / 2.6σ / 2.2σ divided by a σ rounded to 0.12. Unrounded they are **3.3σ / 2.5σ / 2.2σ**.
  The rule outcomes are unchanged: all three fire at 2σ, and only 1/Ω_m fires at 3σ.
- **Lorentz-violation naturalness gap:** α/π = 2.3×10⁻³, not 10⁻². Against bounds of 10⁻¹⁸–10⁻³⁰ the gap is
  **15–27 orders**, not 16–28. The site carried the wrong value on six pages; corrected, with a lint rule.

## 3. Frame question: where does the tanh's shape touch data?

The graduate persona asked whether any density on file exercises the nonlinearity of C = tanh(γ ln(1+x)). Going regime by regime:

- **Density-keyed (galaxies):** x ≈ 7×10⁻⁵, so C ≈ γx is linear and only γ/ρ_crit is identified (10-23).
- **Cosmology (mean density):** the fit sits at γ = ½, where C = x/(x+2) exactly, which is Λ. The local-fluid horn pins
  |2γ−1| ≲ 10⁻⁵.
- **Acceleration-keyed transition (SPARC):** γ̂ is consistent with ½, the Möbius point, which is simple μ.
- **Deep regime (GC outskirts, dwarfs, lensing):** only a₀′/γ enters.
- **Newtonian-return tail (planets):** 1 − C ≈ 2x^(−2γ). Here the shape does reach data, and the result is an
  exclusion: γ < 0.87 is out.
- **External-field quadrupole (Cassini Q₂):** an exclusion, γ < 1.5 out.

So wherever data reach the function, it is linear, sits at the γ = ½ algebraic point, or depends on one combined scale.
The tanh's own shape (γ ≠ ½) reaches data only in the return tail and Q₂, and there it only excludes. **"One equation"
in the surviving content is the Möbius map x/(x+2) plus a list of exclusions on departures from it.** I suggest this as
the operational content of the motto, for STATUS/FUNDAMENTALS. That edit **gates on the project lead**. It does not move
a bucket, but it changes what the headline equation is said to contain.

## 4. TEST-10 (explorer 2026-10-29, drained to the site today)

With distance, inclination and M/L propagated, the most permissive cap still leaves 3 galaxies at ≥ 3σ (NGC 2915,
NGC 3741, UGC 06667) and 6 at 2σ. ESO 444-G084 falls to 1.65σ under its inclination error. Its pairing with NGC 3741 as
"the M/L-free pair" (07-30) is retired on the site. The condition is the gas census: multiplying the gas by 2.0 removes
NGC 3741 and by 3.1 removes NGC 2915, against dark-gas proposals of ×2–3 (Pfenniger & Revaz 2005). A uniform factor of
that size is disfavoured by the BTFR (c ≤ 1.25 at Υ* = 0.5). The count-floor sentence (10-29) stands, with both TEST-10
legs now worded census-conditional.

## For the project lead
1. Adopt "the Möbius point plus exclusions" as the stated content of the one equation, or reject it.
2. The stopping question (10-04) is sharper today. The last galaxy-scale corner with a prediction distinct from
   MOND+EFE has been executed, and it closed. No registered or proposed galaxy-scale reading now has an outcome that
   differs from MOND+EFE.
