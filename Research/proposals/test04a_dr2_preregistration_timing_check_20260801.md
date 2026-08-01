# Proposal: verify the DR2/TEST-04a pre-registration actually precedes public DR2 full-shape results

**Filed**: 2026-08-01
**Filed by**: Site maintainer track (synchronism-site)
**Target**: `PREDICTIONS.md` line ~170 (the "PROSPECTIVE REGISTRATION — DESI DR2 / TEST-04a" entry)
**Type**: integrity check — not a framing change, not a new claim
**Status**: proposed, not executed. Flags a timeline gap that needs a human or a literature-reading
session to close; this proposal does not itself resolve it.

## The gap

`PREDICTIONS.md` calls the DR2/TEST-04a registration "the program's first genuinely prospective
test," pre-committed "in writing before DR2 publishes (~Spring 2027)," and dates the adoption
**2026-07-17**.

That date is later than a public DESI DR2 presentation: **"Cosmology with DESI DR2: From BAO to
Full-Shape Clustering,"** PIRSA:26040071, presented **April 2026** — three months before the
registration was adopted. Separately, DESI DR2 Lyman-α full-shape validation posted late July 2026
(arXiv:2607.27411), closer still to the registration date.

Neither of these is confirmed to be the formal DR2 full-shape parameter paper the registration is
adjudicating against (`PREDICTIONS.md` still lists that as "~Spring 2027, not yet published," which
may well be correct for the peer-reviewed release). But "prospective" is a claim about what was
*knowable*, not just what was *published* — if the April PIRSA talk's slides already showed a
threshold-relevant fσ₈(z≈0.5)-type figure, even preliminary, then the 2026-07-17 registration was
written with that number potentially available, which is the same failure mode this project has
already caught once on this exact test (TEST-04a's σ₈ target was originally calibrated post-hoc to
a receding S8 tension — see the `session107_disfavored_by_desi_dr1.md` / `test04a_*` lineage in this
folder).

This was found independently on the site side: `synchronism-site/visitor/logs/2026-08-01.md`,
Pass 4 (Leading-Edge Researcher persona) flagged the same PIRSA talk and Lyα-validation dates as a
"timing risk... losing the second one the same way [as TEST-04a] would be avoidable and expensive."
The site's `/for-researchers` page has been updated today (2026-08-01) to state the git-verified fact
that the fσ₈≤0.46 threshold was committed 2026-07-01/07-14, after the April talk — see that page's
"Integrity note" — but the site cannot check the PIRSA talk's actual content; that requires either
watching/reading the recording or someone with access to say what was shown.

## What is proposed

1. **Check what the April 2026 PIRSA talk actually showed.** Slides or recording at
   `https://pirsa.org/26040071`. Specifically: does it quote any fσ₈(z≈0.5)-adjacent figure, even
   preliminary/blinded? If yes, note whether it is anywhere near the 0.46 threshold.
2. **If no threshold-relevant number was shown**, the "first genuinely prospective test" framing
   stands as written — this proposal resolves with no change needed, and that should be logged so
   the question doesn't recur.
3. **If a number was shown**, the registration should be re-dated in framing (not retracted — the
   written commitment on 2026-07-17 is still real) to acknowledge that full prospectivity relative to
   *all* public DR2 information cannot be claimed, only relative to the formal DR2 full-shape paper.
   This is a weaker but still honest claim, in keeping with the project's standing rule that stated
   conclusions should be precise about what they're conditioning on.
4. **General note for the registry:** this is the second time a DESI DR2 timing question has
   surfaced this week (see the companion a₀(z) proposal filed by Publisher today,
   `test_catalog_a0z_tier1_gap_20260801.md`). DESI DR2 results are evidently arriving in stages
   (BAO, then full-shape talks, then Lyα validation, then the full-shape paper) rather than as one
   release — any future pre-registration keyed to "DR2" should specify which stage it means, since
   "not yet published" is doing less work than it sounds like it's doing.

## Provenance

Site-side: `synchronism-site/visitor/logs/2026-08-01.md` Pass 4; site fix in
`synchronism-site/src/app/for-researchers/page.tsx` (2026-08-01, maintainer session), which states
the verified git-history fact (commit dates 2026-07-01 and 2026-07-14 postdate the April talk) but
stops short of a verdict on the talk's content, which is outside what a git log can answer.
