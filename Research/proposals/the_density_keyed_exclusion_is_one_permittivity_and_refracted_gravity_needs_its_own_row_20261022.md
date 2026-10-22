# Proposal: the density-keyed exclusion covers one permittivity under named readings; Refracted Gravity needs its own status line

**Date:** 2026-10-22
**Author:** site maintainer (synchronism-site), from a simulated-visitor researcher persona and the explorer's 2026-10-21 finding
**Gates on:** project lead (the ledger-row decision and the adopt-or-not question in §3). Nothing here moves a bucket.
**Count:** 6 (unchanged). **Bucket 0:** 0 (unchanged).

## 1. What "the density-keyed branch is excluded" was resting on

The site and this ledger say that the density-keyed branch is excluded "on SPARC (ΔBIC +2843) and by lunar laser ranging",
without saying which permittivity or which reading. A researcher persona pointed out that the same site credits the field
equation ∇·[C(ρ)∇Φ] = 4πGρ to Refracted Gravity (RG; Matsakos & Diaferio 2016) and also says RG "is not excluded". Both
statements can only hold if "excluded" is scoped. I checked the provenance:

| Root | What it actually tested | Source |
|---|---|---|
| ΔBIC +2843 | **algebraic** g_obs = g_bar / tanh(γ ln(1+ρ/ρ_c)), **unfloored**, ρ = Σ/2h with constant h = 0.3 kpc, global ρ_c, against the implicit C_g | `synchronism-site/maintainer/scripts/dbic2843_galaxy_level.py` docstring (imports `explorer/findings/scripts/two_pillars_head_to_head_fit.py` unchanged) |
| 2026-09-09 field-equation run | the framework's tanh-in-log-ρ form under the solved field equation on a grid f ∈ {0.089, 0.315}, γ ∈ {0.489, 2}, ρ_c ∈ [3.2×10⁻⁴, 0.161]; χ²/N 69–2700 vs MOND 21.2 | PREDICTIONS 09-09 block |
| LLR | any pointwise density-keyed G_eff smoothed over ℓ ≥ 1 pc; 68/74 laws including both published RG calibrations; non-covariant reading | PREDICTIONS 09-24 (i), 10-20 (i) |

**None of these refits RG's own permittivity on discs.** Cesare et al. fit RG to 30 DiskMass discs, and nothing here has
refit them (already noted 08-27; the explorer topic `refracted-gravity-sparc-refit-under-l2-plus-striction.md` is queued).

## 2. RG's status as this archive has it, in one place

- **Globular clusters:** passes under every reading of its Eq. 4.1 (10-14).
- **LLR:** both published calibrations excluded on the pointwise non-covariant reading at ℓ ≥ 1 pc (09-24). The covariant
  completion (Sanna, Matsakos & Diaferio 2023) is unchecked.
- **Dwarf spheroidals (explorer 2026-10-21, PREREG site `71c6253`):** in spherical symmetry Gauss's law gives g = g_N/ε
  exactly, so the boost is capped at 1/ε₀. The disc calibrations fail all 15 MW dSphs by ≥ 3σ (Leo I 9–15σ, which RG keeps in
  equilibrium). The elliptical (E0) calibration gives 2–4σ on the equilibrium subset. For the rest, the refraction response is
  ≈ (Δε/ε)(g_ext/g_int), order unity, so by RG's own equation they are not in equilibrium.
- **RG's distinct, unrun prediction:** host-aligned dipolar (m = 1) lopsidedness of low-density satellites. Newtonian tides
  are quadrupolar.

**Recommendation (gates on the project lead):** add a comparison-class line, "Refracted Gravity (prior art for the field
equation): status", beside the density-keyed Bucket 2 rows, carrying the four bullets above. It is not a framework row. It
stops the record from saying "excluded" about the class when the evidence covers this framework's permittivity.

Applied site-side today: a new `/test-catalog#what-excluded-covers` section, scope clauses on /mond-unification,
/honest-assessment and /for-researchers, and a distinct-content register at `/for-researchers#distinct-content`.

## 3. Question the frame: should the framework adopt RG?

The persona's closing advice was that "the framework's best route to a distinct, reachable prediction is to adopt or refute
RG's flattening dependence". I think that is half right, and the wrong half matters.

- **Adopting RG cannot move Bucket 0.** Any win on refraction or flattening is RG's, published a decade earlier. This is the
  09-27 nesting result in a new place: if the adopted theory *is* the parent, the only possible outcomes are "RG wins" or
  "RG and the framework lose together".
- **What adoption would buy is a dynamics.** The open frame question on the site is "which dynamics?" (L2 has no action and
  breaks the third law for composite bodies; L3 adds striction that is sourced by stellar self-energy). RG has a covariant
  completion. Adopting it would settle L2 versus L3 by import rather than by choice, and would make "C(ρ)" name one theory with
  one set of predictions.
- **It would also import RG's problems.** On this archive's own runs those are the LLR exclusion of the non-covariant reading
  and the dSph ceiling under the disc calibration.

So the honest choice for the project lead is:
(a) adopt RG's dynamics and permittivity, accepting prior-art status for the galaxy sector and inheriting its dSph and LLR
record; or
(b) keep the undeclared family, and say plainly that its density-keyed realisation has no dynamics of its own.

Either is defensible. What is not defensible is the current state, where the record borrows RG's equation and also calls the
branch dead.

## 4. A cheap-kill suggestion that would test the wrong object

The same persona suggested testing the T_scan ≳ 1.6–9.4 µs survivor against µs-period sidebands in trapped-ion Ramsey or
optical-clock data. The site's own T_scan finding says a **global persistent** clock already fails at every T. The survivor is
a scan phase drawn once per particle, which produces no periodic sideband. So clock sidebands would re-test the excluded
version. The discriminator the finding names, g²(0) against photon duration at fixed window, remains the only one. This is
recorded in the new distinct-content register.

## 5. Two smaller items

- **The γ = ½ nesting point** (1/C − 1 = 2/x gives the simple ν on acceleration and Λ on mean density) was already stated on
  three site pages and in this ledger (09-23, 09-24 ii). A graduate persona still missed it, because it sat mid-paragraph. It is
  now a boxed callout at the top of /galaxy-rotation and /dark-energy. No research change.
- **/measurement-without-observers: Reparametrization → Speculative.** The 10-19 move did not fit the site's own definition
  ("equivalent to existing physics in different notation"). The mechanism has no rate parameter and rests on the
  audited-negative γ = 2/√N_corr, so it is not equivalent to decoherence theory, which has rates. No ledger row changes.
