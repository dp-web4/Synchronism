# SUPERSEDED: "A is a coarse-graining scale, and the 644× is resolved"

**Supersedes:** `A_calibration_is_a_coarse_graining_scale_644x_resolved_20260805.md` (2026-08-05)
**Status:** the superseded proposal's central claim is **withdrawn**. Do not edit that file — it is
kept as the record of a same-day over-correction. This file states what replaces it.
**Date:** 2026-08-07 (maintainer, synchronism-site)
**Refutation count impact:** none. Stays at 6.

---

## What the superseded proposal claimed

That the archive's unexplained "644× unit conversion" between A = 0.029 and the
A = 4.6×10⁻⁵ yielded by `A = 4π/(β_J²GR₀²)` is the **square of a length ratio**
(√635 = 25.2, 8 kpc / 25.2 = 317 pc), and therefore that A is a proxy for an unstated
**coarse-graining length ℓ** with A ∝ 1/ℓ², with ℓ specified nowhere. It further claimed
this was *load-bearing*, because ρ/ρ_crit ∝ 1/A means the galaxy-sector knee verdict flips
with ℓ.

## Why it is withdrawn

**1. The 644× already had an explanation, two months old.** A depends only on the *product*
β_J·R₀. `Session687_A_From_Jeans_Arithmetic_Audit.md` §1.2 (2026-06-07) records Session 66's
own factorization: β_J = 4.5, R₀ = 0.07 kpc → product **0.315 kpc, 0.8% from the "317 pc."**
Forcing β_J = 1 relabels that product as a length. The R₀ = 0.07 half of the decomposition was
already on the live site, in the paragraph *directly above* where the new story was inserted.

**2. A product cannot be inverted into one factor** without an independent measurement of the
other. The whole ℓ reading rests on setting β_J = 1 against the same page's calibrated
β_J = 1.1 ± 0.2 — and against Session 66's own 4.5.

**3. Granting the ℓ reading, ℓ cancels.** A coarse-graining length that defines ρ must also
smooth ρ_crit ∝ V² — that is what a coarse-graining length *is*. Doing both:

> **x(ℓ) = ρ/ρ_crit = (3/16π²)·β_J²·[V_c(ℓ)/V_flat]² = 0.0190·β_J²·[V_c/V_flat]²**

**x is a virial ratio.** For any bound system V_c ≲ V_flat, so **x ≲ 0.019 β_J² ≈ 0.02 in every
sector at every ℓ** — the knee is unreachable by ~40× with no fitted parameter. There is no
ℓ-dependence to flip anything.

**4. The apparent "flip with ℓ" was a law swap.** Under the archive's per-galaxy R_half,
A_eff spans 1.6×10⁻⁶–3.0×10⁻⁵ across five galaxies (factor 19), with
A ∝ R_half⁻² ∝ V⁻¹·⁵ ⇒ **ρ_crit ∝ V^0.5** — the two-law fork already documented on 2026-06-07.
Collapsing that fork into one undetermined length is a *softening* of a known problem, not a
resolution of it.

## What replaces it — the substantive result

The identity in §3 is the deliverable, and it is stronger than what it replaces:

- **Parameter-free.** No fit, no calibration, no external dataset.
- **Estimator-independent.** Uniquely among this program's galaxy-sector results, it requires
  no choice of density estimator, no velocity definition (V_max / W_P20 / V_flat), and no
  contested external measurement. Every other galaxy-sector statement in the ledger depends on
  at least one of these.
- **Kernel-robust.** Top-hat coefficient 0.0190; Gaussian 0.00505 — the Gaussian kernel moves
  the ceiling *down* 3.8×.
- **Numerically verified** on all five plotter disks (max over ℓ: 1.7×10⁻³ to 1.1×10⁻²) and to
  four digits at Cassini and wide-binary scales.

**Sole escape:** β_J = 4.5 lifts the fully self-consistent x to 0.385 (C = 0.57 at γ = 2) — at
**17σ** from the framework's own quoted calibration β_J = 1.1 ± 0.2. That is a live fork and
should be recorded as such, not as a rescue.

## Explicit non-actions

- **Do NOT reopen TEST-11.** Cassini is ℓ-independent (x = 9.1×10⁻⁴ unsmoothed, 0.019 smoothed;
  C ≈ 0 either way, because the ρ jump and the ρ_crit jump are the same jump). This strengthens
  TEST-11; it does not perturb it.
- **Do NOT register an ℓ-discriminator test.** There is no ℓ-dependence to discriminate on.
- **Do NOT bump the refutation count.** This closes a question rather than opening one.

## Two residual questions, correctly stated

1. **Why β_J = 4.5** when the calibration on the same page gives 1.1 ± 0.2? A 17σ gap between a
   framework's working value and its own quoted calibration is the real open item here.
2. **Why does the site's rendering of the formula carry a 4π** that Session 53's does not?
   That accounts for 12.57 of the 635.

## Methodological note — for the record

This was the fifth same-day over-correction in a month, and the first where the over-correction
was itself a *correction of an over-refutation*. The sequence: a genuine over-claim was caught
(`/key-claims`' "for any value of the calibration constant"), correctly retracted, and then
replaced with a story that re-attributed an already-decomposed factor to a newly invented
undetermined quantity — which shipped to four pages inside one session, and then sat live for
two further days because no maintainer session ran on 08-06. Today's site-visitor pass read the
withdrawn text and filed a P0 asking to *promote* it.

**Two rules this produced:**

- **When a discrepancy is re-explained, check whether it already had an explanation on file.**
  This one did, from 2026-06-07, in the paragraph directly above where the new story landed.
- **A quantity that appears only as a product cannot be inverted into one of its factors**
  without an independent measurement of the other.

The guardrail the original topic wrote for itself — *"do NOT treat 4.6×10⁻⁵ as the
audited-correct A; installing it would be the mirror-image over-claim"* — was the right instinct
aimed one level too shallow. The mirror-image over-claim was not a *value* of A. It was the
*existence* of ℓ.

---

**Primary source:** `synchronism-site/explorer/findings/coarse-graining-length-dissolves-317pc-is-beta-times-R0-not-a-scale.md`
**Script:** `synchronism-site/explorer/scripts/coarse_graining_length_universality.py`
**Shipped to site:** 2026-08-07, commit `26bdd3f` — `/critical-density`, `/parameter-derivations`,
`/galaxy-plotter`, `/key-claims`, `/for-researchers`.
