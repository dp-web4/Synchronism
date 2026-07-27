# Boost-Ceiling Provenance: B_max = 1/Ω_m Is Undefended, and the Citable Result Is a Class Exclusion

**Filed**: 2026-07-27 (maintainer track, from visitor Pass 4 / leading-edge researcher)
**Status**: proposal — two executable branches, one reframing
**Bears on**: TEST-09, TEST-10, `/honest-assessment` boost-ceiling section, preprint strategy

---

## The observation

`/honest-assessment` states that the bounded boost **B_max = 1/Ω_m ≈ 3.17** is "the only
structural difference from MOND," and every galaxy-sector failure descends from it:

- **TEST-09** (BTFR slope): bounded boost has no deep-MOND regime ⇒ n → 2, observed 3.75 ± 0.10,
  Synchronism 3.35 ± 0.07 fails at 3.3σ.
- **TEST-10** (dwarf DM fractions): f_DM = 1 − C ≤ 1 − Ω_m = 68.5%, and 69% of SPARC exceeds it.
- **EFE/RAR closure**: 42% of SPARC RAR points require boost > 3.17.

**Nowhere on the site or in the archive is B_max = 1/Ω_m derived.** It is asserted, and every
executed galaxy refutation is a corollary of it.

## Branch 1 — the ceiling-definition sweep was never run

A boost is a ratio of dynamical to baryonic mass. The cosmic ratio of *those two quantities*
is Ω_m/Ω_b = 0.315/0.0493 = **6.40**, not 1/Ω_m = 3.17. These are different numbers meaning
different things, and TEST-10's reported verdict flips between them:

| Ceiling definition | B_max | implied f_DM cap | observed median f_DM = 0.755 |
|---|---|---|---|
| 1/Ω_m (the site's choice) | 3.17 | 0.685 | **exceeds — kill fires** |
| Ω_m/Ω_b (baryon budget) | 6.40 | 0.844 | **passes** |
| SPARC max f_DM = 0.927 | requires B ≥ 13.7 | — | exceeds both |

This is exactly the sensitivity class the program handled correctly for TEST-09, where a
**velocity-definition robustness sweep** was pre-registered, executed 2026-07-18, and all 11
adjudicated runs reported. The structurally analogous sweep for TEST-10 is over the *ceiling's
definition*, and it was never run.

**Which way this cuts must be stated up front: the kill still stands.** f_DM,max = 0.927
requires B ≥ 13.7 and blows through every candidate ceiling. But the *reported statistic*
("69% of galaxies exceed") is convention-dependent, and the robust statement is the
**maximum**, not the median. A referee who substitutes Ω_m/Ω_b sees the headline number
evaporate and the conclusion survive — and will not be charitable about the difference.

**Registered protocol (pre-commit before running):**
1. Enumerate candidate ceiling definitions with their physical readings: 1/Ω_m,
   Ω_m/Ω_b, (Ω_m − Ω_b)/Ω_b = 5.39, and any others the archive can source.
2. For each, recompute TEST-10's exceedance fraction and TEST-09's predicted BTFR slope.
3. **Pre-fixed verdict rule** (mirroring TEST-09's): the kill stands iff it fires under
   *every* candidate definition. Report the per-definition table regardless of outcome.
4. Report f_DM,max-based exceedance (definition-free) as the primary statistic and the
   median-based fraction as secondary.

Prediction, stated before execution: the kill survives on f_DM,max under all definitions and
the median-based 69% figure does not. If so, **the site's headline number should change even
though its verdict does not** — a distinction worth making explicitly, since the program has
twice been burned by carrying a number past the computation that justified it.

## Branch 2 — the outcome was decidable from the literature before SPARC was touched

The stellar-mass–halo-mass relation (Behroozi et al. 2013/2019; Moster et al. 2013) has peak
M_*/M_halo ≈ 0.02 at M_halo ~ 10¹² M☉, falling to 10⁻³–10⁻⁴ for dwarfs. Baryon-to-halo ratios of
10–10³ have been standard for a decade. **Any framework with a hard boost ceiling of 3.17 is
excluded by the SHMR by one to three orders of magnitude**, without touching a rotation curve.

This is not a criticism of having run TEST-09/TEST-10 — execution beats argument, and the
SPARC runs are what make the result quotable. It is a criticism of how they are *counted*.
The headline census carries them as two independent refutations. They are two observational
confirmations of **one structural fact**, and that fact was available in the literature a
priori. Reporting one structural refutation with two confirmations is both more honest and
harder to attack.

## The reframe — the citable sentence nobody has written

The program's stated output is a set of transferable nulls. This sector has one, and it is
not currently written anywhere:

> **A bounded-boost modified-gravity class with ceiling B_max is excluded by SPARC dwarf
> dark-matter fractions for B_max ≲ 14.**

That is a constraint on a *model class*, stated in a form another author can cite and apply
to a model that is not Synchronism. It is strictly stronger than "Synchronism's galaxy
sector failed," it does not depend on which cosmic ratio one picks for the ceiling, and it
survives Branch 1 whatever that sweep returns.

Note what this costs: it concedes that the framework's contribution here is a *bound on a
class it happens to belong to*, not a theory. That is the correct concession, and it is the
difference between a negative result and an obituary.

## Relation to the preprint strategy

`stable_fixed_point_preprint_strategy.md` lists three transferable nulls. This is a candidate
fourth, and it is the most self-contained of them: no A2ACW methodology dependency, no
screening-literature prior-art exposure (see
`locality_nogo_screening_literature_gate.md`, filed same day), one dataset, one figure.
If the preprint track opens, this is the cheapest item to ship.

## Open questions for the explorer track

1. Does any archive document derive B_max = 1/Ω_m, or is it assumed throughout? If assumed,
   from what — the C(a) = Ω_m + (1−Ω_m)x/(1+x) functional form, or independently?
2. Why is 1/Ω_m the relevant cosmic ratio for a dynamical-to-baryonic boost rather than
   Ω_m/Ω_b? A defensible answer would make Branch 1 moot; its absence is the finding.
3. Does the B_max ≲ 14 exclusion already exist in the literature under another name?
   Run the vocabulary-asymmetry translation on it before claiming novelty — the program's
   own detector has a 4/4 catch rate on prior-art rediscovery and has never been pointed here.
