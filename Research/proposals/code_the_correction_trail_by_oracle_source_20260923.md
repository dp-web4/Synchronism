# Proposal: code the correction trail by what caught each correction (the program as its own oracle experiment)

**From:** site maintainer, 2026-09-23 (WAKE phase, from the visitor log of the same date)
**Status:** proposal, not a result. No bucket moves; count 6; Bucket 0 = 0.
**Gates on dp:** only item 3 (cadence). Items 1 and 2 are explorer-track work and can start without a decision.

## What prompted it

For the second day running, the visitor researcher persona concluded that the most citable thing in the program is a
*methods* claim rather than any physics: the /a2acw diagnosis that "same-corpus self-play without an external oracle
converges on internal consistency, not discovery; the boundary is the oracle, not the ambition". The persona also drew
a corollary the site had not drawn. This program *did* have a non-corpus oracle, the tests it executed against SPARC,
DESI, Cassini and the Baumgardt–Hilker catalogue. Every verdict that mattered came from that oracle and none came
from debate. So the program's own history is a natural experiment for its own thesis.

The same day's log is itself a data point. All six high-severity items in it were **consistency** catches, found by
reading with no data run. Examples: /dark-energy said "γ = 0.489 is exactly MOND's simple μ" while two other pages and
this ledger's B2 row say γ = ½; /coherence-function contradicted itself about which C the ceiling kills used; a
"crossing" was described as a quadrant that is mostly thawing. None of them changed a physics verdict. The thesis
predicts exactly that.

## The claim, stated so it can fail

> **H-oracle.** In this program's record, corrections that change a physics verdict (bucket, refutation count, or a
> test's pass/fail/underpowered status) come from executing something against data. Corrections found by reading or
> argument alone, whether by archive sessions, the site tracks, visitor personas or external LLM reviewers, change wording,
> attribution, internal consistency and scope, but not verdicts.

**Refuted if** a reading-only correction changed a verdict. One candidate is already known: TEST-04a moved from
"refuted" to "underpowered" on 2026-07-14 after a citation check. That check read the registered criterion and DESI's
published numbers and executed nothing new. Whether reading a published *measurement* counts as the oracle is the
boundary the coding must fix **before** coding starts. Otherwise the result can be made to come out either way.

## What to do

1. **Pre-register the coding scheme** (explorer). Units: every dated correction in PREDICTIONS.md plus the site's
   revision notes, reusing the 2026-09 correction-cohort sampling frame (`explorer/scripts/correction_cohort_analysis*`)
   where possible. Add two fields:
   - *caught-by*: executed-new-computation / read-external-measurement / read-internal-text / argument-only.
   - *effect*: verdict-changing / scope / consistency / wording.

   State H-oracle's prediction as a contingency table. Fix the rule for "read-external-measurement" in advance.
   Two raters, report κ; the 09 cohort got κ = 0.33 on valence, so expect to need a tighter codebook.
2. **Run it** on the existing trail. This analysis can produce a *positive* finding. After ~3,300 sessions, the other
   such place is the blind post-cutoff A2ACW arm, which is still unrun. This one is much cheaper.
3. **Cadence question for dp.** Three consecutive researcher-persona passes (09-21, 09-22, 09-23) have found the
   remaining physics errors to be bookkeeping: which C, which γ, which quadrant. They are worth fixing, but they don't
   discriminate anything. If H-oracle holds, it gives a principled reason to shift maintainer effort from daily
   re-auditing of settled physics to (a) the methods result and (b) the few executable, still-open items: the TEST-04a
   re-registration before DR2, the dark-matter fork, and the blind A2ACW arm.

## Why this is the right scale

Daily site fixes only improve the wording of verdicts the physics has already given. The question the program can still
answer with data it already has is about its own method. This proposal moves one step toward that without dropping
the honest-assessment discipline.

Topic seeded: `synchronism-site/explorer/topics/code-the-correction-trail-by-oracle-source.md`.
