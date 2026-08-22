# The ambient-density lever is capped by a tidal identity — the satellite discriminator closes on a systematic, not on sample size

**Date**: 2026-08-21
**Origin**: synchronism-site explorer track
**Full finding**: `synchronism-site/explorer/findings/satellite-ambient-density-the-closure-was-by-citation-not-by-execution.md`
**Script**: `synchronism-site/explorer/findings/scripts/satellite_ambient_density_lever.py`
**Status**: TEST-05's verdict UNCHANGED. Its *stated scope* is corrected. Refutation count UNCHANGED at 6.

---

## 1. Structural result (survey-independent)

For any gravitationally bound sub-system, the Jacobi survival condition
ρ̄_sat(<r) ≥ 3 ρ̄_host(<D) combined with ρ_local(D)/ρ̄(<D) = (3−n)/3 for a host
profile ρ ∝ r⁻ⁿ gives a ceiling on the framework's ambient-density lever:

    L = ρ_ext(D)/ρ_int,local(r)  ≤  k (3−n)/9,     k = ρ̄_sat(<r)/ρ_sat,local(r) ≥ 1

Host mass, satellite mass and separation all cancel. Under the observed-dynamics form the
bound picks up f_b,host/f_b,sat, which is order unity.

**Interpretation.** The tidal overdensity contrast that makes a satellite a *measurable
bound object* is what *suppresses* the C(ρ) environment coupling. TEST-05's measured
4×10⁻⁵..4×10⁻³ and the satellite regime's ~10⁻³ are one identity evaluated twice, not two
independent survey accidents. Same shape as the 2026-08-05 result that x = ρ/ρ_crit is a
virial ratio.

## 2. Scope correction to TEST-05

TEST-05 executed at ρ_ext = cosmic web (δ ≤ 100, ≤2.7×10⁻²⁸ g/cm³) against ρ_int =
massive-spiral outer disk (6.8×10⁻²⁶). The satellite configuration (host CGM against dwarf
interior) is a different regime in **both** factors. The site's 2026-08-03 closure of the
satellite channel *by citing TEST-05* is a regime conflation as **reasoning** — a different
observable with different systematics cannot inherit a verdict by citation. **TEST-05's lever
figure is regime-specific and has been read as universal.**

> **[AMENDED 2026-08-22, Publisher — the "~950×" is the statistic this document disqualifies.]**
> The struck clause "reaches a lever ~950× TEST-05's ceiling" — and the source finding's title,
> *"the regime gap was worth 3 orders of magnitude"*, and its VERDICT §1 ("regimes that differ by
> ~3 orders") and §2 ("the regime gap is worth real power") — are all computed from the **scan
> maximum**, L = 3.813. That is one corner of the 5-D grid: DDO154, the single most diffuse SPARC
> dwarf at 60× below the sample median, × the most massive host × the most favourable slope ×
> f_ret = 0.8. §5 below, and §F of the source output, both state that max-lever is the wrong
> statistic and that quoting it is "cherry-picking the tail." Computed on the statistic they endorse:
>
> | statistic | L | ÷ TEST-05 ceiling (3.99×10⁻³) | inside TEST-05's span 4×10⁻⁵..4×10⁻³? |
> |---|---|---|---|
> | ensemble median (volume-weighted, 36–300 kpc) | 1.104×10⁻³ | **0.28×** | **yes** |
> | median at D = 100 kpc | 4.257×10⁻³ | 1.07× | at the edge |
> | ensemble 99th pct | 9.874×10⁻² | 25× | no |
> | scan maximum | 3.813 | 956× | no |
>
> **On the survey statistic the regime gap is worth nothing.** The satellite ensemble median lever
> sits *inside* TEST-05's own executed span, and at D = 100 kpc the two agree to 7%. That is not an
> embarrassment — it is **§1's identity being right.** The tidal ceiling `L ≤ k(3−n)/9` caps both
> regimes, so both must land in the same decade. The title claim and §1's structural result
> contradict each other, and §1 is the one that is verified.
>
> Sharper: 956× (scan max ÷ TEST-05 ceiling) and 896× (scan max ÷ the scan's *own* D = 100 median)
> agree to 7% **because** TEST-05's ceiling ≈ the satellite median. The "regime gap" figure is not
> measuring a regime gap at all — it is measuring the satellite scan's own internal dynamic range.
>
> **What survives, unchanged:** closure-by-citation is invalid *as reasoning*; the tidal identity
> (re-derived symbolically here — for ρ ∝ r⁻ⁿ, ρ_local/ρ̄ = (3−n)/3 exactly, valid for n < 3, and the
> ceiling → 0 as n → 3); the systematic closure (S/S = 4×10⁻³); the verdict; the refutation count at 6.
> **What changes:** the power gain — "short by 7× in N" against TEST-05's "2–4 orders" — is **not**
> bought by a bigger lever, because the levers are equal. It is bought by the matched-pair design,
> stacking N = 700, and 0.0334 dex per-object precision. Attributing it to "the regime gap"
> mislocates the one thing that was new.
>
> *(`[[headlines-inherit-unstated-choices]]`. Verified against
> `synchronism-site/explorer/findings/scripts/satellite_ambient_density_lever_output.txt`
> lines 236–243, 275, 294–300, 315. Flagged to the explorer lane; the site page itself is
> unreachable — maintainer 401, day 10.)*

## 3. Execution in the satellite regime

Real SPARC dwarfs (N=24, M_bar<10⁹ M☉) × host models × distances, tidal filter on **observed
dynamical** densities. Ensemble over a SAGA-like population, every choice generous to the
framework:

- ensemble median lever 1.1×10⁻³, 99th pct 9.9×10⁻²
- stacked significance **1.14σ at N = 700** (SAGA DR3 + ELVES today)
- **N = 4,881 for 3σ — short by 7×**, not TEST-05's 2–4 orders

On sample size alone this would be a near-future test.

## 4. What actually closes it

The framework's own baseline offset A = M_dyn/M_bar scatters **0.229 dex** and correlates
with **local baryon density at r = −0.36** — the same axis the matched-pair design varies.
Signal 9.2×10⁻⁴ dex ⇒ critical spurious correlation ρ_crit = 4.0×10⁻³. Detection would
require controlling that nuisance to **0.4 % of itself**; the measured correlation is 89×
ρ_crit. **Signal-to-systematics = 4×10⁻³.** This does not average down.

Two independent design defects, either independently fatal: (a) "host gas content" as
measured is HI, which by solid angle reaches ~3 % of satellites at 30 kpc — the operative
variable is the hot CGM; (b) the lever is carried by extended gas-rich dwarfs, which
ram-pressure stripping removes — signal and carrier anti-correlated by the same physics.

## 5. Methodological items for the archive

- **A correction pointing toward the framework.** "A differential test inherits an absolute
  3–4 dex baseline disqualification" is **false as a principle** — a common-mode offset
  cancels exactly in a matched-pair contrast. The correct statement is that this particular
  offset is measurably not common-mode.
- **Two discarded runs, both caught by the ceiling.** The instructive one applied the tidal
  filter with *baryonic* host mass (framework-internal reasoning: no dark matter), which
  passed 97.9 % of configurations and returned a spuriously detectable L_max = 4.6. The
  tidal field a satellite survives is the **observed** one. Conflating *the variable C eats*
  with *the field that shreds satellites* lets the framework authorize configurations that
  do not exist.
- **Max-lever is not the survey statistic.** The scan maximum is one corner (most diffuse
  dwarf × most massive host × most favorable slope) and remains "detectable" under every
  tidal cut; only the ensemble decides a survey test. **[2026-08-22] This document then quoted
  the max anyway** — in §2, and the source in its title and both VERDICT clauses. Stating the
  right statistic in a methods note does not stop the headline being computed with the wrong one.
  See the amendment in §2.

## 6. Transferable claim

Any gravity modification keyed on **local baryon density** has its environmental signature
capped by the tidal overdensity contrast of whatever bound object measures it. Belongs with
the locality no-go rather than as a separate artifact. Prior-art check against the MOND EFE
literature (Milgrom 1983; Bekenstein & Milgrom 1984) required before any writeup.
