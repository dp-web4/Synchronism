# The persona loop is now the site's review instrument, and it has no measured catch rate

**From:** synchronism-site maintainer, 2026-11-24
**Source:** visitor log 2026-11-24 (four simulated-visitor personas, all AI), Pass 4 (researcher) frame point, checked against source
**Gates on:** the project lead (one decision, two arms; the second arm changes visitor-track infrastructure the maintainer must not touch)

## What the feedback revealed

The researcher persona counted about 100 attributions to "a simulated-visitor persona" across 8 pages. I re-counted in
source: **34 on /honest-assessment, 22 on /tier-1-existing, 21 on /for-researchers, 320 mentions of "persona" across
`src/`**. Since July, the persona loop has become the site's main source of corrections. It is one model family, it
runs daily against the live site, and nobody has measured how often it catches a real error or flags a correct statement
as wrong.

This is the same-corpus objection the site raises against A2ACW, where specificity was 0/6 by construction and one model
both wrote and scored. A2ACW's numbers are printed; the loop that produces most of the site's corrections has none.

Today's log is evidence in both directions:
- **Catches that were real:**
  - the plotter's "0 of 123 reach Ω_m" came from the retired knee: rerun, 111/123 at the current knee;
  - the "600×" divided coefficients carrying different units;
  - the wide-binary crux box ruled against both groups' own framing. I checked their abstracts and confirmed it.
  - Draco's σ input was never printed; the arithmetic reproduces.
- **Catches that were fetch artifacts,** which the persona dropped itself: the breadcrumb, "107", and "6 cards". Memory
  records earlier cases where the maintainer chased tallies that were artifacts.
- **Known misses:** the stale "+184 refuted" sat on /equation-walkthrough because steps 2–8 are client-rendered and a
  fetch cannot see them. The 2026-11-14 knee switch orphaned a survey number for 10 days. The plotter was visited in 8 of the
  9 visitor logs from 11-14 to 11-22 (grep of the logs) before one persona put the two numbers side by side on 11-24.

So the loop is productive, and its blind spots are partly structural: fetch rendering, pages nobody edits, numbers far
from their labels. Neither of those is a rate.

## The proposal (stake in the ground)

**Measure the catch rate with planted errors, and publish it next to A2ACW's specificity.**

Prediction, written before any run: **for planted numeric errors that need recomputation (a wrong conversion, a stale
calibration, a mixed-units ratio), the loop catches fewer than 1 in 3. For planted wording overclaims (an upgraded badge,
"refuted" for "disfavoured"), it catches more than 2 in 3.** The false-positive rate on unplanted correct statements is
above 10 % of flags before the persona's own artifact filter, and below 5 % after it.

Refuted if numeric catches are ≥ ½ or wording catches are ≤ ⅓.

### Arm A: in-family, offline (the explorer can run this; no infrastructure change)
- Build a static export of a frozen commit and plant N ≈ 20 errors in five classes:
  - arithmetic;
  - stale calibration;
  - citation/locator;
  - overclaim;
  - over-refutation, since the 05→08 over-refutation share rose 1/20 → 9/24, so the loop should be tested on that side too.
- Run the four persona prompts through `claude -p` with read access to the exported HTML. Count catches per class, then
  repeat on the unplanted export for the false-positive baseline.
- Register the planted set and the scoring rule in a sealed file before any run. The planter should not be the scorer.
- Limitation: it is not the real instrument. WebFetch summarisation is part of the real loop's error profile.

### Arm B: the real instrument (gates on the project lead)
- One visitor run against an unlisted Vercel preview deployment of a branch carrying the planted set, through the
  visitor track's normal WebFetch path. The maintainer does not modify visitor-track configs, so this needs the project
  lead or the visitor track's owner.
- **Never plant errors on the public site.** A deliberate falsehood on a public page is an outward-facing act and is not
  proposed.

### Cross-vendor note
The 10-20/10-21 A2ACW record found same-vendor review harsher than cross-vendor review in all 17 disagreements. A
cross-vendor persona arm would show whether that asymmetry carries over to the loop's catch profile. Same gate as the
existing cross-vendor-rater question.

## The second frame point, restated because DR4 is 8 days out

Every outward-facing item on `/for-researchers#other-theories` is still untold to the people it concerns. Two of them
are time-sensitive before Gaia DR4 (2 December): the wide-binary companion-rate analysis, and the unreproduced
consolidated Γ (0.102 vs 0.139). Today's crux-box fix shows the cost of not asking. The site had ruled on the
dispute's cause in a way that contradicted both groups' published framing; a persona caught it, the authors could not.
Contacting Chae et al. and Saad & Ting **gates on the project lead** (first routed 2026-11-09; draft issue text in the
11-11 proposal).

## What changed on the site today (record)
- `/galaxy-plotter` and Tier 1 TEST-09: the survey rerun at the current knee
  (`synchronism-site/maintainer/scripts/plotter_survey_at_archive_knee.py`; predictions committed `9fe6951` before
  running, 2 of 4 held, identity exact).
- `/wide-binaries`: the crux now leads with both groups' framing, and Chae 2023/2024 are added to the list.
- `/for-researchers`: the dim ≥ 5 row reads "at reach by threshold kinematics" (explorer 2026-11-22), and Draco's σ input
  is printed.
- Site commit `d3d9962`.
