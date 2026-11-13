# Outage-window backfill audit — disposition at #14's head (2026-11-13)

**Asked for by:** NOVA's exact-head review of dp-web4/Synchronism#14 (2026-11-13 04:14Z): "a one-time
publisher-window backfill audit".
**What it audits:** the ten STALE-IN-PAPER items the 2026-10-18 pass found in the 09-03 → 10-07 publisher
outage (log entry 2026-10-18; `state/whitepaper_sync.json` `watch_items[0]`), plus the 09-23 item the
2026-10-30 re-scan added. The 10-30 re-scan independently reproduced 10-18's items 3, 4, 6, 7 and 8.
**Read against:** `publisher/wp-2026-10-29` at `ca78abce` (the 11-12 head), live `whitepaper/sections/**`
`<topic>.md` files, excluding `meta/` and changelogs. Line numbers are from that head. They are unchanged at
`53e0214f`, which only edits table cells in place.
**Method:** each item's printed phrase was located with a span-printing scan
(`/tmp/pub1113/backfill_scan.py`, not committed: its "annotation follows" flag is paragraph-level and
over-reports on 16 KB one-line paragraphs). Every disposition below was then read by hand at the cited line.
**Not done here:** the source `📌` blocks were not re-read today. Item descriptions are the 10-18 entry's,
and that entry's quotes were spot-checked on 10-18 and 10-30.

## Disposition

| # | Item (source block) | Printed at (`ca78abce`) | Status | Closed by | What remains / gate |
|---|---|---|---|---|---|
| 1 | TEST-09 bounded-boost as a "surviving decisive negative" (09-18 sweep; 09-25 (i)) | body: §5.15 count item (2) L138 and the floor-convention bullet; **headline: exec summary L189, §5.15 status note L5, conclusion L77** | **PARTIAL** | 10-29 + 10-31 (sweep printed, unrounded 0.406 / 0.297 / 0.260); 11-06 (registered-rule note on count item (2)) | The headline sentence "both surviving decisive negatives are galactic (locality no-go; TEST-09 BTFR bounded-boost)" is unscoped at all 3 sites. Whether a convention-dependent kill counts **gates on dp**. |
| 2 | "Compatibility-Synthon scaling (r = 0.994)" as calibrating work (09-25 (iii): p_crit ∝ 1/⟨C⟩ refuted on its own data) | exec summary L17 | **OPEN** | none | Scope gap (10-30 reading). Needs a re-read at source before inscription. |
| 3 | "ρ_crit ∝ V⁺² sign inversion kills independently" (10-07 (iii): conditional on forward regression; orthogonal gives V^+2.1) | exec summary L191, conclusion L79 | **OPEN** | none | Partly **gates on dp** (10-18). |
| 4 | "EFE = 0 exactly" / "{EFE = 0} ⟺ locality" (09-16 (ii), 09-17 (iii): algebraic reading only; under L2, ∇C ≠ 0 refracts the field) | §5.15 L111, L113, L117; exec summary L187/L191; conclusion L73/L79 | **OPEN** | none (L99's "EFE = 0 is a property of the density-keyed function" is a different, earlier point) | This is the "locality no-go" in item 1's headline sentence. Partly **gates on dp**. |
| 5 | "Dark matter AND dark energy are both coherence effects" / "no Λ required", no Cardassian prior art (09-22, 09-14) | §5.15 L364 | **PARTIAL** | 11-06 (L342 names the fluid, Cardassian-class reading as the only survivor of the LLR test) | L364's "both coherence effects" sentence is unscoped. |
| 6 | "no order-ε⁰ channel / retires the perturbation caveat" (09-15: only without the pressure term) | §5.15 L366 ("no order-ε⁰ term exists") | **OPEN** | none | "pressure term" has 0 hits. Rests on site-executed numbers, so it needs a re-run before inscription (10-30). |
| 7 | "factor-1.96 a₀ tension dissolves at ϒ = 0.6" (10-04 (i)) | §5.15 L366 | **CLOSED** | 10-30 (CHECKED from L109) | none |
| 8 | lensing-ceiling deficits 6× / 13× / 28× via √f (09-25 (i): B_req = ν/f) | §5.15 L242 | **CLOSED** | 11-11 (Brouwer read in full; 1/f, not √f); 11-12 correction | none |
| 9 | "no parameter-free obstruction on the coarse-graining axis" (09-15 + 09-24 (i): LLR bound) | open questions L51 | **OPEN** | none | 0 annotations at the sentence. |
| 10 | "a₀ = cH₀/(2π) DERIVED" and "the 2π is not numerology" (09-19) | §5.15 table L16; §5.15 L256 | **PARTIAL** | 11-11 (table note) + 11-13 (in-cell qualifier) | L256 "This is not numerology; it connects the MOND scale to cosmology" is unscoped. |
| — | 09-23: the sign lock excludes thawing, not only crossing | §5.15 L366 (the 10-30 log cited L346 at that head) | **CLOSED** | 10-30 (CPL algebra) | none |

**Tally:** 3 closed (7, 8, 09-23), 3 partial (1, 5, 10), 5 open (2, 3, 4, 6, 9).

## What the table shows that the list did not

1. **The headline's two surviving negatives are exactly the two open items that reach the headline.**
   Exec summary L189, §5.15 L5 and conclusion L77 each say "both surviving decisive negatives are galactic
   (locality no-go; TEST-09 BTFR bounded-boost)". The locality no-go is item 4: its EFE = 0 ⟺ locality
   step holds under the algebraic reading only. TEST-09 is item 1: it fires under the 1/Ω_m floor and not
   under the two Ω_b floors. The body of §5.15 now prints the TEST-09 condition. The headline sentence that
   relies on both prints neither. This is the sentence a first-time reader meets, and it is the one with the
   most outage debt.
2. **Per-item annotation has closed 3 items, 2 more in part, in 26 days.** Each closure came from a pass whose own window
   happened to touch the item. No pass worked the list as a list. The 10-18 scale question (a dated status
   banner at §0 and §5.15, a rewrite of the galaxy-sector verdict paragraphs, or per-item notes) has not
   been answered, so per-item has remained the default by inaction.
3. **Five open items need a decision before any edit, not just a pass.** Items 1, 3 and 4 gate on dp in
   part. Items 2 and 6 need source work (a re-read, a re-run) that no window will bring in.

## Recommendation (for dp; not acted on)

- Merge #14 on its current qualifications, if the in-cell qualifiers answer the status-cell half of the hold.
  Then take the **banner** option for the five open items: one dated paragraph at the head of §0 and §5.15
  that names items 1–4, 6 and 9 and points to `PREDICTIONS.md`. That costs one edit. The per-item route costs about
  ten notes on paragraphs that are already 10–16 KB.
- Whatever the choice, item 1's headline sentence should not stay as it reads, because it is the most-read line with
  the most debt.

— cbp-claude (autonomous-publisher-cbp, 2026-11-13)
