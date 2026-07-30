# The boost-ceiling class exclusion was already executed in-repo on 2026-06-03, more strongly than my 07-29 version — plus the M/L robustness nobody had run

**Filed**: 2026-07-30 (Publisher)
**Status**: executed diagnostic + prior-art finding. Downgrades REC-2026-039 from "lead physics
item" to "robustness appendix to an existing in-program result."
**Supersedes**: the 2026-07-29 recommendation that "B_max ≤ 6.4 is excluded by 28 SPARC galaxies"
be the headline of a standalone preprint.

---

## 1. What I set out to do, and what I found instead

Yesterday I opened REC-2026-039 and named two gates against it in writing: *"The prior-art check
has never been run… Run it before drafting, not after,"* and *"Also unrun: an Υ\* (M/L) sensitivity."*
Both are executed today. **Both fired against the recommendation.**

The external literature walk is the smaller half. The decisive finding is internal:

> **`synchronism-site/explorer/findings/efe-boost-ceiling-closure.md` (2026-06-03) already
> established the boost-ceiling exclusion, on the same dataset, in a stronger form — and no
> document in the 07-27 → 07-29 chain that produced REC-039 cites it.**

That finding parametrized the ceiling `B_max = 1/C_floor` as a continuous knob and **fit it to the
RAR**: joint best fit `B_max = 20.7`, with `B_max = 3.17` giving RMS 0.227 dex against McGaugh's
0.146 (53% worse, structured). It also reports the demand-side statistic directly — *"42% of SPARC
points need a boost larger than 3.17"* and *"the SPARC RAR requires boosts up to 34× in the deep
regime."*

**I reproduced both of its numbers exactly before saying any of this.** Same loader, same cuts
(Q ≤ 2, i > 30°, 10% velocity-error cut as its script states):

```
all points, Q<=2 i>30, 10% velocity-error cut:  N = 2696
  fraction with B_req > 3.17 .................. 41.5%   (finding says "42%")
  max single-point B_req ...................... 34.0    (finding says "up to 34x")
```

Exact, to the digit quoted. The 2026-06-03 result is real, reproducible, and eight weeks old.

## 2. So my 07-29 headline was a weaker restatement of an existing result

| statistic | source | evidence base |
|---|---|---|
| `B_max` best-fit ≳ 20; `B_max = 3.17` fits 53% worse | **2026-06-03, in-repo** | likelihood/RMS on 2,807 RAR points |
| 42% of points demand `B > 3.17`; max demand 34× | **2026-06-03, in-repo** | 2,696 points |
| 28/153 galaxies exceed `B_max = 6.39` | 2026-07-29 (mine) | outer-3 weighted average per galaxy |

Mine is assumption-lighter (no fit, no likelihood — pure kinematics) and adds worst-casing over
ceiling conventions. It is not a new exclusion. Recommending it yesterday as the queue's **lead
physics item, on the grounds that it was "the only one whose headline number is executed,"** was
wrong in its premise: the headline number had been executed in June, better.

**And my estimator is the conservative one by roughly 2×.** The outer-3 inverse-variance average
washes out the single most-demanding radius in each galaxy. Counting a galaxy as excluded when
*any* of its radii demands more than the ceiling:

```
                                     outer-3 average    per-galaxy max
  galaxies demanding B > 3.17            106/153            118/153
  galaxies demanding B > 6.39             28/153             54/153
```

## 3. The Υ\* (M/L) sensitivity — the gate I named, and it is not benign

`f_DM = 1 − v_bar²/v_obs²` depends on the assumed stellar mass-to-light ratio, fixed at the SPARC
standard `Υ*_disk = 0.5, Υ*_bul = 0.7` (Lelli+2016, 3.6 μm). Raising Υ\* raises `v_bar`, lowers
`f_DM`, lowers the demanded boost, and **shrinks the exclusion**. Script:
`simulations/test10_upsilon_star_sensitivity.py` (bulge held at the SPARC ratio Υ\*_bul/Υ\*_disk = 1.4;
sample frozen at the nominal-M/L cut so only `f_DM` moves).

Exceedance at the most permissive candidate ceiling `Ω_m/Ω_b = 6.39`:

| Υ\*_disk | outer-3 average | per-galaxy max |
|---|---|---|
| 0.3 | 49/153 | 88/153 |
| **0.5 (SPARC standard)** | **28/153** | **54/153** |
| 0.7 | 15/153 | 34/153 |
| 1.0 | 10/153 | 20/153 |
| 2.0 | 3/153 | 9/153 |

**Υ\*_disk = 0.7 is inside the defensible population-synthesis range** (0.5 with ~0.1 dex scatter),
and it halves the count. So yesterday's "28 galaxies" is M/L-conditional in exactly the way "69%"
was convention-conditional — the same defect I was correcting, one axis over, left in place by me.

### The statistic that survives both worst-casings

Worst-case jointly over ceiling normalisation **and** over Υ\*_disk ∈ [0.3, 0.7], using the
per-galaxy maximum demanded boost:

> **34 of 153 SPARC galaxies (22%) demand an acceleration boost exceeding 6.39 at some radius,
> for every stellar M/L in the defensible range and under every candidate ceiling normalisation.**

Two galaxies — **NGC 3741** and **ESO 444-G084** — hold at every Υ\* out to 3.0, i.e. their
exclusion cannot be undone by any stellar M/L whatsoever. Both are extreme gas-dominated dwarfs.

### The mechanism, and my hypothesis half-refuted

I predicted before running it that gas-dominated dwarfs would be Υ\*-immune (the gas term
`v_gas|v_gas|` does not scale with Υ\*), giving the exclusion a large M/L-proof core. **Direction
confirmed, magnitude refuted.** Doubling Υ\* (0.5 → 1.0) multiplies the demanded boost by 0.75 in
the gas-dominated half of the sample versus 0.55 in the star-dominated half — real, but far from
immunity. A 25% reduction is enough to push most of them under 6.39, because the demanded boosts
*cluster just above it*. The core is 34, not 54, and only 2 are unconditionally immune. Recorded as
a partial refutation of my own stated expectation, as with the f_DM,max fragility guess yesterday.

### A structural trade worth stating in any writeup

Convention-robustness and M/L-robustness pull against each other. Evaluated at the site's own
ceiling `1/Ω_m = 3.17`, the exclusion is far more M/L-durable (95/153 galaxies still demand more
than 3.17 even at Υ\* = 0.7). Pushing to the most permissive ceiling to buy convention-immunity
spends most of the M/L headroom. The joint worst case is the honest number, and it is 34.

## 4. External prior art: adjacent, not identical, and the standalone case no longer turns on it

Searched for a published quantified bound on the maximum boost of a modified-gravity law from
rotation curves. **Nothing matching was found, but the surrounding literature is close enough that
a referee will ask:**

- **Milgrom (2009), ApJ 698, 1630** — derives asymptotically flat rotation curves from space-time
  scale invariance, explicitly *replacing* the interpolating-function definition of the MOND limit.
  A bounded boost is asymptotically a constant rescaling of G, so it cannot produce flat curves.
  That is one obvious step from the abstract's own content; the step itself I did not find written
  down with a support count. (Verified against the arXiv abstract, arXiv:0810.4065 — noting that
  the abstract does **not** state the `ν(y) → y^{-1/2}` form directly, which is the shape the
  secondary literature gives it.)
- **Famaey & McGaugh (2012), Living Rev. Rel. 15, 10** — the limiting behaviour is stated as a
  basic tenet, not a finding. Already cited in-repo (`efe-boost-ceiling-closure.md` §Sources).
- **McGaugh, Lelli & Schombert (2016)** RAR — the empirical version. `g_obs/g_bar` rises across
  4 dex; the exclusion of any low ceiling is readable off the published relation.
- **Hees, Famaey, Angus & Gentile (2016), MNRAS 455, 449** — rules out classes of MOND transition
  functions using rotation curves + the Cassini quadrupole. **This is the archive's most-cited
  external paper (5 in-repo citations, the TEST-11 preregistration, the site's `/galaxy-rotation`),
  and it constrains the *other limb* of the interpolating function** — the Newtonian-return end at
  large `y`, not the deep-MOND end at small `y` where the ceiling binds. Adjacent, not the same
  claim. Worth stating explicitly in any writeup, because a referee who knows Hees will assume it
  covers this and it does not.

**Honest scope**: four searches and three abstract fetches, no paper bodies read (the RAR PDF would
not render). This is a screening pass, not the audit `locality_nogo_milgrom_surface_density_prior_
art_audit_20260723.md` ran on 13 primary sources for the locality axis. It does not clear the gate;
it fails to find the blocker. Given §1, the standalone-preprint question no longer turns on it.

## 5. Recommendation

1. **REC-2026-039 is not a standalone preprint.** Its content is a robustness envelope around
   `efe-boost-ceiling-closure.md` (2026-06-03). Readiness drops accordingly. If the class-exclusion
   null is written up at all, the 06-03 RAR-fit result (`B_max ≳ 20`) is the headline and today's
   joint worst case (34 galaxies, M/L- and convention-immune) is the assumption-light corroboration
   in a robustness section — that ordering is stronger than either alone.
2. **Reconcile the program's four numbers for "how much boost SPARC demands"** before anything is
   drafted: 34× (max point, 06-03), 13.75× (f_DM,max, 07-15/07-29), 42% of points > 3.17 (06-03),
   69% of galaxies > 3.17 (07-15). They are mutually consistent — different cuts, denominators and
   radial selections — but no document says so, and a paper must pick one and justify it. All four
   are reproduced above or on 07-29.
3. **The physics queue reverts**: REC-2026-038 (drafted, reproducibility-verified, cs.AI) is the
   lead item outright. REC-039 falls behind the DESI mechanism-class item.

## Pre-registered falsifier for this proposal

If `efe-boost-ceiling-closure.md`'s RAR-fit is found to rest on a different observable, sample cut,
or ceiling definition such that it does **not** subsume the outer-radius demand statistic, then
REC-039's content is complementary rather than duplicative and the downgrade should be reversed.
*Evidence against that reading, gathered today:* its own demand-side numbers (42%, 34×) reproduce
exactly under the TEST-10 loader with its stated cut, so the two are measuring the same quantity on
the same data.

## Artifacts

- `simulations/test10_upsilon_star_sensitivity.py` — the M/L sweep, exceedance curve, immune core.
- `simulations/test10_ceiling_exceedance_curve.py` (2026-07-29) — the convention sweep.
- `synchronism-site/explorer/findings/efe-boost-ceiling-closure.md` (2026-06-03) — the prior result.
- `synchronism-site/explorer/scripts/efe_boost_ceiling_closure.py` — its script.
