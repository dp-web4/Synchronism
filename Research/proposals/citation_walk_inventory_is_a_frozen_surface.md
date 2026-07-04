# Proposal: The Contribution Inventory Is the Highest-Risk Compilation Surface — a Standing "Superseded-By" Rule

**Filed**: 2026-07-04 (synchronism-site explorer track)
**Origin**: Systematic run of the citation-walk audit proposed 2026-07-03
(`research_contributions_top3_demoted_retraction_survival_pattern.md`).
**Type**: Meta-methodology / archive hygiene

## What Was Run

The 2026-07-03 back-annotation proposed a new A2ACW failure-taxonomy clause — **retraction survival in
compilation layers** — and noted the citation-walk audit (trace each compilation number to its origin
session, check whether a later session retracted it) had caught 2/2 informally but never been run
systematically. This session ran it systematically over the 47-contribution inventory.

## Result — the clause is confirmed, and the inventory *itself* is the largest instance

The site's `/honest-assessment` states the 47 contributions' remainder is "not yet swept … this site's
only unaudited surface; the demotion prior … correspondingly low." Walking the citations shows this is
false against **this archive's own sessions**:

- **Count.** S582 (2026-02-08) = **30** (strict criteria). S615 (2026-02-17) = **47** (looser, "bookkeeping,
  no new physics"). **S634 (2026-04-25) flagged 47 as a 57% overcount and recommended reverting to 30**;
  S618/S627 already cite 30. The site (and the whitepaper) still say 47.
- **Characterization.** The Framework Stress Test arc (S617–638) and the execution arc (S661/S668/S669)
  already resolved the inventory's entries: cosmology domain **0 novel-unfalsified** (S635); RAR σ_int
  candidate ~120× below floor (S637); C(ρ) → MOND (S637) and → **Curie paramagnet, CAS-verified** (S638);
  chemistry r≈0.98 network = **Debye model relabeled, Δr=0 executed** (S669); BTFR/α refuted (S631).
  "8-for-8 site-claim audits over 9 days" by S638's own tally.
- **Intra-document contradiction.** The whitepaper carries both "47 genuine contributions (Final
  Accounting)" and, ~5 pages later, "S634: 47 is a 57% overcount … S638: Curie paramagnet." Same file.

## The Generalizable Finding

**An inventory / "final accounting" document is the highest-risk compilation surface in a long-running
program**, precisely because it is authored once, stamped canonical, titled *Final*, and then never
re-run — while the research it counts keeps moving. S615 froze on 2026-02-17; the archive kept auditing
for 76 more sessions and demoted essentially everything the inventory listed; the freeze won on every
downstream surface (site, whitepaper).

This is the same mechanism as the TEST-04a and CDM σ_int leaks, one level up: not a single number
outliving its source, but a *summary* outliving the entire audit arc that superseded it.

## Recommendation

1. **Standing "superseded-by" rule.** Any accounting/inventory document carries a live pointer to the
   latest session that revised its headline figures. The canonical count is whatever the *latest* audit
   session asserts — not the one titled "Final." A doc titled "Final Accounting" is a liability without
   this pointer.
2. **The count is 30, and it is fully characterized, not unaudited.** Update downstream surfaces (site,
   whitepaper) to S634's 30 (or keep 47 with the explicit looser-classifier caveat), and replace
   "unaudited / low demotion prior" with the S635–S638 result: predictive content fully characterized —
   cosmology → MOND, chemistry → Curie paramagnet, both CAS-verified; novel-surviving yield 0.
3. **Mechanize the citation-walk.** The inventory tables give origin session IDs; a grep over
   `Research/Session[6-9]*.md` surfaces later revisions in seconds. Wire "contribution → latest-touching-
   session → verdict" into the daily loop so a frozen inventory can't silently outlive its audit arc again.

## So What?

The 2026-07-03 clause predicted retraction-survival would be systemic in the remaining contributions. It
is — and the cleanest instance is the meta-surface: the list that *summarizes* the contributions is a
Feb-17 snapshot the archive's own April–May audit arc already overturned. The honesty ledger's own
summary was carrying the pre-audit optimism.
