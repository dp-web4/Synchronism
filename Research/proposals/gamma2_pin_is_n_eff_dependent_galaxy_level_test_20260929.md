# The γ = 2 pin is N_eff-dependent: at the galaxy level the "+184" kill is ΔBIC ≈ 11 and 1.6–2.2σ (site maintainer, 2026-09-29)

**Status:** executed, pre-registered (site commit `976e3f9` before the script). Count unchanged (6); Bucket 0 = 0. Whether root 2 of 5 ("γ = 2 pin") keeps its place in the headline gates on dp.

## What was asked
A visitor graduate-physics persona (site log 2026-09-28) noticed that the site scores the γ = 2 kill on the SPARC RAR under two effective-N conventions: /galaxy-rotation divides the +184 by ~5.6× ("ΔBIC ≈ 33, still decisive"), /core-idea divides the +2843 density-vs-acceleration number by ~20× (one point per galaxy). Applied to the +184, the 20× convention gives ≈ 9, under the site's own ΔBIC > 10 rule. PREDICTIONS.md 2026-09-24 (iii) says "the +184 and +2843 verdicts survive any of the three conventions." For the +184 that is arithmetically false at 20×. Nobody had estimated N_eff from the data.

## What was run
`synchronism-site/maintainer/scripts/gamma2_pin_galaxy_level.py` (+ `_PREREG.md`, `_output.txt`, `_diag.py`, `_diag_output.txt`). Frozen 2026-07-22 pipeline (`simulations/sparc_tanhlog_profile.py`: same row cut, Υ_disk = 0.5, Υ_bulge = 0.7, log-space SSR, a₀′ profiled), with the galaxy name kept per row. 2,807 points, **166 galaxies**.

Identity controls reproduced before any new number was read: γ = 2 vs McGaugh ν **ΔBIC = +184.0**; free γ̂ = **0.489**, ΔBIC vs McGaugh **+7.1**.

## Results against the registered rules
| | Registered prediction | Result | Verdict |
|---|---|---|---|
| P1 sign test | γ = 2 worse in > 65 % of galaxies, p < 10⁻³ | worse in **94/166 = 56.6 %**, sign p = 0.051, Wilcoxon p = 0.068 | **REFUTED** |
| P2 galaxy-level 10-fold CV | free γ beats γ = 2 by > 3σ (galaxy-level SE) | held-out ΔlnL = +76.9 ± 48.9 ⇒ **1.6σ** | **REFUTED** |
| P3 ΔBIC with N → N_gal | 5–15 | **+10.9** | HELD |
| P4 galaxy-block bootstrap of the full-N ΔBIC | 95 % lower bound > 50 | median 179, **95 % [23, 351]**; sd 82 vs ~19 expected under point independence ⇒ **N_eff ≈ 150** | NEITHER (lower bound 23) |
| P5 signal in galaxies with median g_bar/a₀ < 1 | > 70 % of the paired excess | 51 % | REFUTED |

**Registered verdict rule:** the pin stays "refuted, convention-free" iff P1 and P2 both hold. **Neither does. The γ = 2 pin's refutation is N_eff-dependent.**

Diagnostic (not pre-registered): the net paired excess (+3.95 in log10² units) is +9.52 from galaxies where γ = 2 is worse minus −5.58 where free γ is worse; **seven galaxies carry 98 % of the net** (UGC 11914 21 %, NGC 2841 16 %, UGC 02953 15 %, NGC 5985 14 %, IC 2574 11 %, DDO 161 11 %, UGC 03205 10 %), three of them high-acceleration spirals (median g_bar/a₀ = 4.0, 2.6, 0.9). The ΔBIC ≤ 10 retained γ interval is (0.425, 0.60) at N = 2807, (0.35, 0.85) at N/5.6, and **(0.3, 1.7) at N_gal = 166**.

## What this does and does not change
1. **The +184 is a point-independence number.** The empirical replication unit is the galaxy (bootstrap spread ⇒ N_eff ≈ 150, matching the 08-14 σ(γ) result that the fit is galaxy-limited). At that N the same comparison is ΔBIC ≈ 11 (Wilks p ≈ 10⁻³, ~3σ if the rescaling is taken literally), 2.2σ cluster-robust in-sample, 1.6σ out-of-sample. "Decisive" (the site's word) and "~8σ per bin" are not supported. "Disfavoured at the threshold" is.
2. **The 2026-09-24 (iii) sentence "the +184 and +2843 verdicts survive any of the three conventions" is withdrawn for the +184.** The +2843 (density- vs acceleration-keyed) survives every convention (2843/17 = 167).
3. **TEST-25's "robust empty intersection under every recorded BIC convention" is convention-dependent in the same way:** every recorded convention was point-independent. At N_gal the SPARC-retained interval reaches γ = 1.7, and the 2026-09-17 block records that the compander passes Cassini post-hoc at γ ≳ 1.5–2. So the in-house squeeze is not robust to N_eff. **TEST-25's defensible content is unchanged**, because since 2026-09-17 it rests on the published, marginalized 8.7σ (Desmond, Hees & Famaey 2024) for the RAR-preferred functions, which MOND shares. What moves is the *in-house* instrument's claimed independence from conventions.
4. **This is not a rescue of γ = 2.** Nothing here makes γ = 2 preferred; it is worse in 57 % of galaxies and at the ΔBIC threshold. The Desmond et al. analysis marginalizes per-galaxy Υ, distance and inclination, which this pipeline (like the frozen artifact) holds fixed; that is the stronger instrument and it should carry the SPARC side of the squeeze. I could not re-read it (no network); its SPARC-side exclusion of fast-return functions is the number to cite.
5. **Method lesson (fourth instance of the same family):** a ratio, a σ or a ΔBIC quoted on 2,807 points from ~170 curves inherits an N_eff assumption. The site's own 08-20 galaxy-CV method existed and was not applied to the headline kill. Rule for the site: any SPARC ΔBIC gets its galaxy-block bootstrap interval beside it.

## Recommendations (gate on dp)
- Ledger Bucket 2 row "C(ρ) ⇒ MOND (γ free)": append the galaxy-level numbers; relabel the γ = 2 pin "disfavoured at the ΔBIC threshold (N_eff ≈ 150), not decisive".
- Root count: the 09-17 finding that the over-refutation share of corrections rose applies here. Recommend the headline sentence read "2 framework-specific mechanism roots (boost ceiling; γ = 2 pin, the latter at the threshold once galaxies are the unit)". Whether it stays a root is dp's call; I am not recounting.
- Register the per-galaxy-nuisance refit (Υ, D, i free per galaxy, under each law, covariance reported) as the next execution on this row. It is the visitor's request and the published method.
