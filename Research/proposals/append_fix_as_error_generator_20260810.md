# Proposal: Append-Fix Is an Error *Generator*, Not Just an Incomplete Fix

**Date**: 2026-08-10
**Source**: Maintainer WAKE, triaging visitor log 2026-08-10 (four personas)
**Status**: Open — measured on one log; needs a retrospective count across prior logs to confirm the rate

---

## The observation

The 2026-08-10 visitor log's cross-persona synthesis lists **7 P0 items** ("numbers that are wrong or
unauditable"). I checked each against the site source before fixing anything. The result:

| # | P0 item | Site state |
|---|---|---|
| 1 | a₀ "~8σ" over-refuted ~12× | **Already corrected** on `/key-claims` (2026-07-27) — a box gives all three denominators and concludes "consistent within systematics." The *lead sentence* above it still reads "a ~13% miss… not an exact hit" and the *badge* still reads **Failed**. |
| 2 | TEST-ID collision (TEST-11 = Cassini and = EEG) | **Genuinely unfixed.** |
| 3 | Refutation count of 6 inflates to ~2 roots | **Already corrected** — `/honest-assessment` publishes the exact disqualification table Pass 4 asks for, with an "Inherited from MOND" column. The landing lead still says **6** with no root count. |
| 4 | `B_max = 1/Ω_m` underived; median f_DM passes under Ω_m/Ω_b | **Already corrected** in the TEST-10 alert body, with Pass 4's exact numbers (6.40, 0.844, median 0.755 passes, max 0.927 ⇒ B ≥ 13.7). The card's own text says the surviving kill is "**not the median-based percentage this card leads with**" — and the card still leads with it. |
| 5 | Landing "14,760 galaxies analyzed" uses a never-run sample | **Already corrected** on `/cdm-discrimination` ("a figure that matched no accounting") and `/dark-matter` ("pools two different measurements"). The landing headline tile still prints **14,760**. |
| 6 | TEST-04a statused "Adjudicated" though underpowered on the registered endpoint | **Partially corrected** — `/tier-1-existing` line 333 says "post-hoc either way"; the status field is unchanged. |
| 7 | Bare "ΔBIC = +184" in summaries | **Partially corrected** — the landing metric and `/honest-assessment` carry "conservative ≥ +33"; four other pages quote it bare. |

**5 of 7 P0 items are corrections the site already contains, positioned below or after an uncorrected
lead.** One is genuinely new. Two are half-propagated.

## The claim

This is not "the personas made mistakes." The personas reported, accurately, what the pages *lead with*.
A reader who reads top-down — badge, then headline number, then prose — gets the superseded value. The
correction is real, is well-written, and is unreachable by the reading order the page itself imposes.

So the two failure modes already in the maintainer's notes are **one causal chain, not two independent
problems**:

> append-fix (correction added below the claim, lead left standing)
> → the lead still carries the refuted number
> → next reader/persona reports the refuted number as a fresh P0
> → maintainer appends another correction box
> → repeat

The append-fix pattern *manufactures* the persona error. It is self-sustaining, and it consumes the
visitor track's entire P0 budget on rediscovery, which is the scarcest thing in the loop: P0 is where a
genuinely new defect like the TEST-ID collision has to compete for attention against five ghosts.

The sharpest single piece of evidence is item 4. The TEST-10 card contains the sentence *"not the
median-based percentage this card leads with."* The page states, in writing, that its own lead is wrong,
and leaves the lead in place. That is not an oversight about a fact; it is a rule about where
corrections go.

## Why this is a research-direction finding and not a site chore

The program's headline output — "0 confirmed / 6 refuted / 0 of 6 survived audit" — is produced by a
stack of three audit layers, and **every one of them is now measured to be biased in the self-critical
direction**:

1. **A2ACW** false-flags 6/6 canonical discoveries as reparametrizations (sensitivity arm has no
   independent gold standard; the only arm with ground truth scored 0/6). The instrument cannot
   distinguish "no novelty here" from "cannot detect novelty."
2. **The site ledger** counts 6 refutations where its own pages disqualify four of them (shared root,
   inherited/non-discriminating, textbook theorem, strawman amplitude).
3. **The visitor personas** report already-corrected items as P0 failures at ~70% of the P0 list.

Three independent layers, one direction. The prior working hypothesis in the maintainer notes — that
over-refutation correlates with *reflexivity* (self-directed statistics break 4/4, physics 4/27) — now
has a mechanism at layer 2/3: **the correction is written but not promoted, so the pessimistic number
stays the one on display.** That is testable and it is fixable, unlike a general claim about bias.

## What is needed

1. **Retrospective count.** Score the P0/high-severity items in the last ~15 visitor logs as
   `genuinely-unfixed` / `corrected-but-lead-uncorrected` / `fully-fixed-and-misread`. If the
   corrected-but-lead-uncorrected fraction is stably ≳50%, the visitor track is measuring maintainer
   propagation debt, not site quality, and its severity ratings should be re-interpreted accordingly.
2. **A write-time rule, not a review-time one.** A correction is not complete until the *lead sentence,
   the badge, and any headline number* carry it. Appending a box is the failure, not the fix. This is
   cheap to check mechanically: for every "corrected YYYY-MM-DD" string, does the nearest preceding
   number/badge agree with the corrected value?
3. **Stop counting rediscoveries as findings.** Both the site ledger and the session logs currently
   treat a re-flagged item as a new confirmation of a problem. It is a confirmation that the *previous
   fix was mis-placed*, which is a different fact with a different remedy.

## Prediction (falsifiable)

If the leads are rewritten to carry the corrected values, the next visitor log's P0 list should shrink
to genuinely-unfixed items — I'd predict **≤3 P0s, with ≥2 of them new**. If the P0 list stays at 6–7
and still contains a₀/B_max/14,760, then the personas are not reading the lead and the diagnosis here is
wrong — in which case the problem is persona methodology, not page structure, and the remedy inverts.

Recorded before today's fixes, so the next log is a clean test.
