# Correction to `symmetric_audit_discipline_directional_bias_confirmed.md`: The Directional Law Fails Its Null; Reflexivity Is the Real Predictor

**Date**: 2026-07-09
**Source**: Site explorer session, executing item 3 of the proposal it corrects
  ("Track the tally itself") and Open Thread #1 of
  `explorer/findings/citation-walk-all-statistics-2026-07-08.md`
**Status**: Proposal — supersedes the decision put to dp in the parent proposal
**Reproduces via**: `synchronism-site/explorer/scripts/directional_bias_null_model.py`

## What the parent proposal asks dp to decide

> "Tally as of 2026-07-09: **6 independent provenance-break instances found, 6 over-refute the
> framework's physics, 0 over-claim the physics.** ... Is this worth a standing process change
> (e.g., 'every maintainer session's fix list must include at least one re-derivation of a
> *negative* verdict')?"

**Do not decide on this basis. The statistic does not survive its own null model, and it is
miscounted.**

## 1. The law is not significant even if true

Provenance breaks are found *among statistics*. On a corpus advertising 0 confirmed
predictions, the statistic pool is overwhelmingly anti-framework. In the 2026-07-08
citation-walk — the only *census* on the site (five pre-declared bundles walking "every
load-bearing statistic," so the sample was not selected by suspicion) — **23 of 28
verdict-bearing statistics face against the framework: p_anti = 0.821.**

    P(all 6 breaks land anti-framework | breaks uniform over statistics) = 0.821^6 = 0.31

To reach p < 0.05 the loop would need **16 consecutive over-refutations and zero over-claims**;
p < 0.01 needs 24. The headline was never evidence. Flipping the valence of all four debatable
classifications (2⁴ reclassifications) leaves the one-sided p in [0.756, 0.987].

## 2. The tally is miscounted, and "0 over-claim the physics" is false

The census contains **7 breaks: 5 over-refute, 2 over-claim** — *fewer* over-refutations than
the 5.8 chance predicts. Per-statistic break rate is **higher on pro-framework statistics**
(2/5 = 40%) than anti (5/23 = 22%); the point estimate runs *against* the claimed law.

Two physics over-claims, both found by this program, both live from the initial commit to
2026-07-04 (`git log -S "below CDM" -- src/` → fixed in `541fc88`):

- **σ_int = 0.086 "below CDM prediction — definitive"** on `/honest-assessment`,
  `/galaxy-rotation`, `navigation.ts`, `/cdm-discrimination`. S610's verdict is
  **CDM-consistent at z = +0.5**, and S610 explicitly labels the below-CDM reading *PREMATURE*.
- **`/cdm-discrimination`: TEST-03 environment dependence "confirmed at p = 5×10⁻⁶."**

The "CDM inversion" campaign that produced both is *named in the citation-walk's own list of the
four campaigns that yielded "the six."* It was recorded as an over-refutation.

Two mechanisms preserved 6/6: **scope-exclusion** (the pro-side breaks — A2ACW's denominator,
"3,308 sessions," "47 contributions" — are reclassified as "methodology, not physics" and
dropped; a law that holds only after removing its counterexamples is a definition), and a
**disjunctive predicate** ("every single one over-refutes *or over-closes*" — direction and
confidence are orthogonal axes joined by an "or," which absorbs any break at all).

## 3. What the data do show: reflexivity, not direction

Partition the same census by whether a statistic describes **the world** or **the loop**:

| claim class | breaks | rate |
|---|---|---|
| SELF (3,308 sessions; 47 contributions; 9/9 demotion rate) | 3/3 | **100%** |
| PHYSICS (CHSH, RAR, DESI, LIV, wide-binary, S63, σ_int…) | 4/27 | 15% |

**Fisher exact p = 0.0086** (0.094 under the conservative variant grading the 9/9 membership
caveat as clean). Direction: p = 0.57, wrong sign.

The mechanism is not motivational, it is structural. **Self-referential statistics are the only
claim class with no primary source outside the loop to walk to.** "3,308 sessions" has no
primary: a fresh count finds 650 session files, the component tables sum to ~2,926, and the
chemistry track's own total is 2,685. "9/9" has no canonical list (the two enumerations disagree
on two members). A2ACW's "0 across 3,308" had no denominator. A CHSH S-value, by contrast, has a
committed script with a fixed seed and regenerates bit-identical.

This subsumes the directional observation rather than denying it: over-refuting the physics and
over-claiming the session count both make the audit look more thorough than it was. But the
unified variable is *flattery of the auditing process*, not honesty-branded manufacture of
failure — and only the reflexivity framing names a break class the directional law cannot.

## 4. The error is the one this program already condemned, one level up

`a2acw-detector-false-positive-rate-null-baseline.md` faulted the A2ACW null for publishing a
**numerator with no denominator** — sensitivity on a selected positive set, specificity never
measured, no control arm — and warned: *"don't let the null borrow the test suite's authority."*

The directional law is a numerator (6 over-refutations) with no denominator (how many pro-facing
statistics were ever checked?), a sensitivity with no specificity, and no control arm. **The
audit of the audit inherited the audit's failure mode.** That is the transferable result here,
and it sharpens the monotone-closure finding: the loop's methodological errors reproduce
themselves at each meta-level.

The instances remain sound and serious. TEST-03's kill is manufactured; S63's "0.64 also
rejected" was fabricated; Σ₀'s 0.5% match was reported as a 12% miss. Nothing above rehabilitates
the framework, which remains a reparametrization with 0 confirmed predictions. What fails is the
inference from those instances to a law about direction.

## Revised decision for dp

Replace the proposed rule ("re-derive one *negative* verdict per session" — it would have caught
4 of 7 breaks and missed all 3 self-referential ones, the only class with a 100% break rate)
with:

> **Every number the project reports about itself — session counts, contribution counts,
> demotion rates, detector sensitivity/specificity — must ship with the command or script that
> regenerates it, or be stated as "self-reported; not independently reproducible."**

Mechanically checkable rather than dependent on reviewer disposition, and free: a physics
statistic already has a committed script; a self-referential one either has a `find | wc -l`
behind it or it should not be a number.

## Caveat, stated in the same spirit

The census is n = 30 and the SELF cell is n = 3. **p = 0.0086 on a 3/3 cell is fragile** — one
clean self-referential statistic takes it to 0.094. The reflexivity hypothesis is
*better-supported*, not established, and deserves exactly the skepticism applied above to the law
it replaces.

## Replication, run in the same session: 2/2 predicted breaks confirmed

The hypothesis predicted that **this repo's** reflexive statistics would break. Written down
first, then tested. Both candidates broke; neither had been walked before.

**"1,703 phenomena"** (17 assertions; carried to the site). Nothing computes it; no phenomena
catalog exists. `1703` appears in the chemistry scripts as three inconsistent referents —
`biolubricant…py:17` "Session #1840 | **1703rd Phenomenon Type**"; `cmp_semiconductor…py:22`
"**Finding #1703** | **1639th** phenomenon type"; `membrane_separation…py:3` "Chemistry **Session
#1703**". The repo's own phenomenon-type ordinals **run to 2523**. So 1,703 is an *ordinal
misread as a cardinal*, and it is ~820 short of the largest ordinal the repo asserts. This is the
same error shape as TEST-03's "Session 616" (the catalog header "After 616 core sessions" read as
a session identifier) and the N = 14,585 splice.

**"616 core sessions"** (`Research/EXPERIMENTAL_TEST_CATALOG.md` — the file that registers the
falsification thresholds). `STATUS.md` says **678**; the highest `SessionNNN` file is **691**;
700 session files exist. The catalog that defines every kill criterion is stale by ~75 sessions,
and it is the counter whose misreading manufactured TEST-03's kill.

**Two candidate breaks dissolved on inspection, and are recorded as such**: "1,703 confirmed
predictions" occurs only inside the whitepaper's own disclaimer ("…**not** as 1,703 confirmed
predictions"); the varying core-session values (308/312/610/614/678) are dated publisher-report
snapshots — a time series, not a contradiction.

These are **not** folded into the 2×2 above: adding SELF breaks from this repo without sampling
this repo's physics statistics would inflate one cell by selective expansion — the exact
denominator error under diagnosis. They stand as an out-of-sample prediction test.

**Still untested here**: "47 contributions," "2,671 chemistry sessions," "89% validated."
