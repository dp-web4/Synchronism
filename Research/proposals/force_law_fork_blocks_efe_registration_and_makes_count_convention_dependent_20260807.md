# The force-law fork is not a documentation defect — it blocks the framework's only discriminating test

**Date:** 2026-08-07 (maintainer, synchronism-site)
**Status:** proposal / decision request. Requires a **theory** decision, not new data or analysis.
**Refutation count impact:** none directly. But it makes the current count of 6
**convention-dependent**, which is a different problem from the count being wrong.

---

## The observation

A site-visitor pass on 2026-08-07 filed three items as separate P0s:

1. *"Register EFE = 0 as a numbered test. It is the framework's only structurally discriminating
   prediction vs MOND, it's testable on existing data, and its absence from the ID-keyed ledger
   makes the '0 of 24 could select Synchronism' claim an artifact rather than a result."*
2. *"Pick one force law, site-wide. TEST-09 and TEST-10 fire only under one of them, and both
   refutations are counted with no convention caveat where they're counted."*
3. *"Give every ledger row an ID, or stop describing the inventory by ID range."*

**These are one item.** Specifically: **(1) is not actionable, and (2) is the reason.**

## Why EFE = 0 cannot simply be registered

The framework runs three mutually exclusive readings of the g_bar → ρ substitution, all live
(documented 2026-08-04):

| reading | form | where it lives |
|---|---|---|
| amplitude | v² = v_b² + (V_flat·C)² | `/galaxy-plotter` code, line 106 |
| division | g_obs = g_bar/C | implied by the f_DM = 1−C identity, `/tier-1-existing` |
| multiplication | g_obs = C·g_bar | `/mond-unification` prose |

On the site's own five plotter galaxies these miss the observed curves in **three different
directions** — roughly 10²–10³× high, converging on Newtonian, and 10²–10³× low — and no
calibration reconciles them, because the disagreement is a functional-form mismatch, not an
amplitude one.

**The two C conventions give EFE of opposite sign** (banked 2026-08-04). So there is no single
EFE = 0 prediction to register. Registration is blocked on a theory decision, and has been since
at least 08-04, but nowhere on the site says so — which is why a fresh reader files "register
EFE = 0" as a P0 and will keep doing so.

This is the **self-sealing pattern already identified for a₀(z)**: a prediction with no TEST ID
is invisible to every audit that walks the ledger by ID, and the *silence about why* gets
re-inscribed by each new reader as a fresh finding. Unfixed propagation gaps are re-discovered,
not remembered.

## The second consequence: the count is convention-dependent

TEST-09 (BTFR slope) and TEST-10 (dwarf DM-fraction ceiling) are both corollaries of the boost
ceiling B ≤ 1/Ω_m = 3.17, which follows from f_DM = 1−C — i.e. from the **division** reading
specifically. Therefore:

- They are **one structural root**, not two independent refutations (already noted 08-04:
  "count = 3 roots not 6").
- That root **presupposes one of three live conventions**, and the caveat saying so was stated
  only on `/galaxy-rotation`, where it was found — not on `/tier-1-existing`, where the tests
  are counted. **A caveat that does not reach the page doing the counting is not a caveat.**

The headline "6 refutations" is thus doing two kinds of work it hasn't earned: double-counting
one root, and resting that root on an unfixed convention.

## What is being requested

**A theory decision: which of the three readings is the framework's asserted force law?**

This cannot be settled by fitting — all three fail, in different directions. It has to be
settled by what the framework *asserts*, and then the ledger has to be recomputed under that
choice. Three concrete consequences of the decision:

1. **EFE = 0 becomes registrable** (with a determinate sign) and gets a TEST ID. It is the only
   structurally discriminating prediction against MOND on the board; leaving it unregistered is
   what makes "0 of 24 could select Synchronism" an artifact of the denominator rather than a
   result.
2. **TEST-09/TEST-10 either survive or dissolve together.** Under the division reading they
   stand as one root; under the amplitude reading the ceiling never arises, because C ≈ 0 in the
   disk gives essentially no boost at all.
3. **The count is restated** as independent empirical roots (3–4) rather than test IDs (6).

## Interim action taken on the site (2026-08-07, commit pending)

Pending the decision, `/tier-1-existing` now carries both caveats where the tests are counted:
the convention-dependence and shared root of TEST-09/10, and an explicit statement that EFE = 0
is **absent because blocked, not because overlooked**, with the reason named. This stops the
gap from being re-filed as a novel P0 and makes the blocked state visible to ID-keyed audits.

## What would change this analysis

- **If the framework's assertion is genuinely that C enters differently in different regimes**,
  then the three readings are not competitors and the fork dissolves — but the framework then
  owes a rule for which applies where, and no archive document supplies one.
- **If EFE's sign turns out to be convention-independent** after all — i.e. if the opposite-sign
  result from 08-04 was an artifact of how the conventions were compared — then registration
  unblocks immediately and item 1 above should be executed rather than deferred. **This is the
  cheapest thing to check next and should be checked before anything else here.**

---

**Site source of the fork:** `src/app/galaxy-rotation/page.tsx` C-convention note (2026-08-04);
`src/app/galaxy-plotter/page.tsx:106` for the amplitude coding.
**Related:** `tier1_mond_efe_discriminator_gap.md`, `test_catalog_a0z_tier1_gap_20260801.md`
(same drop mechanism, different prediction).
