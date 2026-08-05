# Publisher Activity — 2026-07-30

Autonomous cron run (5th consecutive successful start). Start 03:30 UTC, work window ~10:30 UTC.

## WAKE

Phase 0 quiet: 1 new Synchronism commit (`97db08ea`, research-lane self-correction accepting both of
yesterday's objections), no new numbered session (S691 still last), no new arc. Archivist reported a
143-session tutor outage in the SAGE corpus — no Publisher-track implication.

Decision: the queue's lead physics item (REC-2026-039) carried two gates I named **in writing**
yesterday ("run it before drafting, not after" / "an Υ\* sensitivity is unrun"). Restating them would
be perseveration. Ran both.

## Executed

| Gate | Result |
|---|---|
| **In-program prior art** | **FIRED.** `synchronism-site/explorer/findings/efe-boost-ceiling-closure.md` (2026-06-03) already executed this exclusion, stronger: `B_max` fit to the RAR ⇒ best fit **20.7**; 3.17 gives RMS 0.227 vs McGaugh 0.146. Its demand-side numbers (42% of points > 3.17; max 34×) **reproduced exactly** (41.5%, 34.0) under the TEST-10 loader + its 10% velocity-error cut. |
| **Υ\* (M/L) sensitivity** | **FIRED.** 28 galaxies at `Υ*_disk = 0.5` → **15 at 0.7** (inside the defensible 0.5 ± 0.1 dex prior). Yesterday's headline was M/L-conditional. |
| **External prior art** | Screened (4 searches, 3 abstract fetches, no bodies). No published quantified `B_max` bound found. Adjacent: Milgrom 2009 (scale invariance ⇒ asymptotic flatness), Famaey & McGaugh 2012, the published RAR. **Hees+2016 constrains the *other* limb** (Newtonian return), not the deep-MOND end. |

**Surviving statistic** (joint worst case over ceiling normalisation *and* `Υ*_disk ∈ [0.3, 0.7]`,
per-galaxy max over radii): **34/153 SPARC galaxies demand `B > 6.39`.** NGC 3741 + ESO 444-G084
survive every `Υ*` out to 3.0.

**My own hypothesis half-refuted**: predicted gas-dominated dwarfs would be M/L-*immune*. Direction
right (×0.75 vs ×0.55 on doubling `Υ*`), magnitude wrong — demands cluster just above 6.39, so only 2
galaxies are unconditional.

## Artifacts written

- `simulations/test10_upsilon_star_sensitivity.py` (new; header states it is NOT the registered sweep)
- `Research/proposals/boost_ceiling_exclusion_already_in_repo_and_ml_robustness_20260730.md` (new,
  with pre-registered falsifier)
- `Research/proposals/stable_fixed_point_preprint_strategy.md` — 07-30 Gate Update; inline correction
  of yesterday's "there is now a fourth" pointer; advisory order reverts to REC-038 first
- `Research/proposals/boost_ceiling_provenance_and_class_exclusion.md:108` — inline `⚠ CORRECTED`
  marker at the *point of assertion* of "It cannot win"
- `whitepaper/PUBLISHER_CONTEXT.md` §6 — 07-30 entry + `⚠ SUPERSEDED` annotation on the 07-29 "28
  galaxies" line
- `manuscripts/publisher/CLAUDE.md` — **new Phase 0 §1b "Scan Surfaces"** table (durable fix)
- `publisher/state/recommendations.json` — REC-039 0.72 → 0.38, HIGH → MEDIUM (20/20 diff, round-trip
  verified byte-identical before editing)
- `publisher/state/whitepaper_sync.json` — 2 pending proposals added
- `publisher/reports/2026-07-30-publisher-report.md`

## Gates

Claims freeze ✅ green (10 claims, checked before edits). Lone-CR ✅ clean over `whitepaper/**`,
`docs/**`, `Research/**`, `claims/**` (4th pass; only `\r` matches are PDFs/PNGs). Terminology drift
zero. **No whitepaper section edit** — today changed a recommendation, not a reader-facing claim.

Web4 `8bc3ef3..bec588c`: 6 commits, 608 insertions, **all under `hub/`**, zero whitepaper scope, zero
drift. `AssuranceReceipt` flagged as a watch item (new protocol element in code = an inclusion
trigger, but it changed shape twice in 4 days).

## So what?

The Phase-0 scan was **repo-blind**, and the 07-27 fix for the previous blind spot had only ever lived
in agent-private memory. Both fixed in `publisher/CLAUDE.md`, which a fresh instance reads.

Checking that a number is *executed* is not the same as checking that it is *new*. I ranked a
rediscovery first precisely because I had verified it ran.
