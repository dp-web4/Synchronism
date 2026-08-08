# The boost ceiling's epoch fork closes the last candidate discriminator — and it closes a priori

**Date**: 2026-08-08
**Source**: maintainer track, back-annotation from `synchronism-site` visitor log 2026-08-08 (Pass 4, Leading-Edge Researcher)
**Status**: closure by internal consistency; no new data used, none needed
**Supersedes nothing.** Extends `boost_ceiling_underived_class_exclusion` (2026-07-27) and
`a0z_row_nondiscriminating_lcdm_degenerate` (2026-08-01).

---

## The finding filed against us

Today's expert visitor pass identified what it called *"the only genuinely discriminating prediction
I found on the entire site"*, and it was a good catch of a real inconsistency:

> The site evolves one cosmological input and freezes the other, with no stated rule.
> `a₀ = cH₀/(2π)` is explicitly promoted to `a₀(z) = cH(z)/(2π)`. But `B_max = 1/Ω_m` uses
> Ω_m,0 = 0.315 and is never given a z-argument. If the ceiling really is cosmological, then
> **f_DM,max(z) = 1 − Ω_m(z)** ≈ 0.21 at z = 1, ≈ 0.05 at z = 2 — epoch dependence neither MOND
> nor ΛCDM predicts, testable today against Genzel+2017 / Lang+2017 / Price+2021 / RC100.

The inconsistency is real and the site should state a rule. But the proposed test **does not
exist**, and which way it fails depends on which ceiling convention you pick. Both branches were
checkable in ten minutes of arithmetic on quantities already on the site.

## Branch (i): the ceiling references the baryon budget, Ω_m/Ω_b — *exactly* epoch-independent

`/tier-1-existing` already carries (flagged 2026-07-28, re-confirmed 07-29) the observation that a
dynamical-to-baryonic boost more plausibly references **Ω_m/Ω_b ≈ 6.40** than 1/Ω_m, and that
TEST-10's headline flips convention: under Ω_m/Ω_b the reported median f_DM = 0.755 *passes*.

Under that reading there is no epoch prediction at all, and not approximately — **identically**.
Baryons and total matter are both dust: Ω_b(z) and Ω_m(z) each carry the same (1+z)³ and the same
E(z)² denominator, so

> **Ω_b(z)/Ω_m(z) = Ω_b,0/Ω_m,0 = 0.1565 for all z**, hence B_max = 6.39 at every epoch.

The ratio is a constant of the cosmology, not a function of it. Under branch (i) the proposed
high-z test has *zero* signal by construction — there is nothing to measure.

## Branch (ii): the ceiling is 1/Ω_m(z) — and it collides with the framework's own a₀(z)

Take the visitor's branch seriously and evolve the ceiling. Then, with Ω_m,0 = 0.315, Ω_Λ = 0.685:

| z | Ω_m(z) | B_max = 1/Ω_m(z) | f_DM,max = 1−Ω_m(z) | a₀(z)/a₀(0) = E(z) |
|---|--------|------------------|---------------------|--------------------|
| 0.0 | 0.3150 | 3.175 | 0.685 | 1.000 |
| 0.5 | 0.6082 | 1.644 | 0.392 | 1.322 |
| 1.0 | 0.7863 | 1.272 | 0.214 | 1.790 |
| 2.0 | 0.9255 | 1.081 | 0.075 | 3.032 |
| 3.0 | 0.9671 | 1.034 | 0.033 | 4.566 |

Closed form: **B_max(z) = E(z)²/[Ω_m,0(1+z)³]**, where the *same* E(z) is what multiplies a₀.
That is the problem. **The two cosmological inputs evolve in opposite directions.**

In the deep-MOND limit — which the galaxy sector *is*, exactly, at γ = 1/2, since
C(x) = x/(x+2) = μ_simple(x/2) — the required boost is B = √(a₀(z)/g_bar). So the ceiling
permits MOND-like behaviour only where

> g_bar/a₀(z) > 1/B_max(z)² = Ω_m(z)²

which is **g_bar > 0.099 a₀** at z = 0 (a usable range), but **g_bar > 0.618 a₀** at z = 1 and
**g_bar > 0.857 a₀** at z = 2. By z ≈ 1 the ceiling has closed off essentially the entire MOND
regime, *while the framework's own a₀(z) has grown 1.8× and widened that same regime by the same
factor*. At z = 2 the framework simultaneously asserts a MOND transition at 3× higher acceleration
and forbids boosts above 8%.

**Branch (ii) is not a prediction; it is a self-contradiction**, and it needs no observation to
see it.

## Consequence: the researcher's question has a determinate answer

The site was asked "is B_max evaluated at z = 0 or at the epoch of the system?" The honest answer
is now available and it is not "we haven't decided":

> **The ceiling must be frozen at Ω_m,0** — not by convenience, but because evolving it
> contradicts the framework's own a₀(z) = cH(z)/2π by z ≈ 1. And once frozen, f_DM,max carries no
> epoch dependence, so no high-z discriminator exists. Under the alternative Ω_m/Ω_b reading the
> ratio is epoch-independent identically, and no discriminator exists there either.
> **Either branch: no test. The count of discriminating tests stays 0.**

Note what this does *not* claim. It does not refute the framework on new grounds — the boost
ceiling is already badged *Asserted, Not Derived*, and TEST-09/TEST-10 are already known to be one
root resting on it. It removes a candidate, which is the cheaper and more common outcome.

## The pattern this is the third instance of

This is now the third consecutive candidate discriminator that died on **an unmade definitional
choice, decidable with no data**:

1. **EFE** (2026-08-07): blocked because three live readings of C's argument (ρ / g_bar / |∇Φ|)
   give opposite-signed EFE. Theory decision, not measurement.
2. **The boost ceiling itself** (2026-07-27/28): TEST-10's headline flips between 1/Ω_m and
   Ω_m/Ω_b; the surviving kill is the narrower tail exclusion B ≲ 14.
3. **The ceiling's epoch** (today): one reading is identically flat, the other is internally
   inconsistent.

The generalisation worth taking seriously at program level: **this framework's discriminating
content is not limited by data availability. It is limited by the fact that its central objects
are underspecified enough that each candidate test forks before it can be registered.** Every one
of these three was resolvable by arithmetic on already-published numbers. That is a statement
about the theory's specification, not about telescopes, and it is the single most compressible
description of why "0 discriminating tests" has been stable for two months while test *proposals*
keep arriving.

Corollary for the ledger: a candidate test should be required to name **which reading of C, and at
which epoch, its prediction is evaluated under**, before it gets a TEST-ID. Two of the three above
would have been caught at registration by that one requirement.

## What was banked to the site (2026-08-08)

- `/parameter-derivations` item 8 (the ceiling's provenance): the epoch fork stated, both branches
  closed, with the table above.
- `/tier-1-existing`: the epoch closure attached to the existing ceiling-convention caveat, so the
  closure lives where the tests are counted — not only where it was found.
- `/test-catalog`: recorded as a candidate considered and closed a priori, so the next reader does
  not re-file it. (Third instance of the self-sealing-silence failure: an unstated *reason* for a
  gap gets re-inscribed by each fresh reader as a new finding.)

## Guardrails for anyone extending this

- **Do not bump the refutation count.** Nothing was refuted. A candidate was removed.
- Branch (ii)'s self-contradiction is an argument for *freezing the ceiling*, not an independent
  refutation of the framework — do not badge it as a seventh failure.
- The high-z f_DM literature (Genzel+2017 *Nature* 543, 397; Lang+2017; Price+2021;
  Nestor Shachar+2023 RC100) is **not** needed for this closure and should not be cited as
  supporting it. It would only matter under branch (ii), which is closed on internal grounds.
  Note independently that low f_DM at high z is also the ΛCDM/simulation expectation for
  baryon-dominated discs — so even had branch (ii) survived, the a₀(z) failure mode
  (ΛCDM-degenerate) was the live risk.
- `Ω_b,0 = 0.0493` (Planck 2018) was used for the 6.39; the site elsewhere quotes 6.40. The
  epoch-independence of the ratio is exact regardless of the value.

Arithmetic: `explorer/scripts/` — reproduced in this document's table; single-file, no data inputs.
