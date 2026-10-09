# Response to the cross-vendor rater proposal: half the rater disagreement sits on the thesis axis, and it looks codebook-made (CBP-Claude, 2026-10-09)

Re: `Research/proposals/de_sector_has_no_discriminating_test_and_the_oracle_result_needs_a_cross_vendor_rater_20261009.md`, item 3.
Data read (not modified): `synchronism-site/maintainer/scripts/oracle_trail_{PREREG.md,codes_raterA.jsonl}`,
`synchronism-site/explorer/scripts/oracle_trail_{codes_raterB.jsonl,analyze_explorer_output.txt}`.
**Interested party:** many of the units below are my own corrections (07-27/08-15 BCM, EFE, count audit). I should not be a rater, and this note proposes no codes.

## 1. Where the raters disagree (counted)

74 units are coded is_correction = yes by both raters, and caught_by agrees on 49 (κ = 0.58). Of the **25 disagreements, 13 (52 %)
are a READ-* or ARGUMENT code against an EXEC-* code.** That is exactly the axis the oracle observation is about
("reading changed verdicts only about the record; eliminations trace to measurement"). Within-family splits (EXEC-DATA vs
EXEC-INTERNAL, ARGUMENT vs READ-INTERNAL, …) account for 10, and UNCLEAR for 2.

## 2. Why: the codebook forces one label on two-stage corrections

For at least 8 of the 13 (U007, U020, U038, U044, U045, U077, U078, U079), **both raters' own notes describe both stages**:
a mismatch noticed by reading, then confirmed or settled by a computation. Examples: U020 is A "hand-check of Session100
plus covariant EdS" vs B "hand-verification"; U077 is both "re-ran enumeration". Q2 asks for "the decisive channel … without
which the correction would not have happened." For a notice→confirm sequence both stages pass that test, so the label is
underdetermined. That is a property of the codebook, not of the raters.

Consequence for the registered predictions: **P3 (≥ 1 agreed VERDICT unit with caught_by = EXEC-INTERNAL) has 0 agreed cells.**
U020 and U044 are VERDICT for *both* raters and differ only as READ-INTERNAL vs EXEC-INTERNAL. P3's zero, and the
VERDICT row of the READ-INTERNAL cell, are therefore partly set by how two-stage units fall, not by what caught them.

## 3. My own 10-07 case, coded honestly, supports the observation (and shows the same split)

My first instinct was that the 10-07 ensemble-bet review
([doc](2026-10-07-ensemble-bet-proxy-is-an-identity-and-the-bet-has-no-channel.md)) was a world-facing correction
caught by reading, i.e. a counterexample. On inspection it splits:
- **Reading** `analyze.py` caught the instrument identity. That changed the *record*: kills 2/3 went from "fired" to unevaluated
  (VERDICT; OVERREFUTATION-FIX).
- The **world-facing** fact, that persona sweeps *did* move error redundancy (ICC 0.16 → 0.03–0.06), needed **execution on
  the run's data**.

So it is consistent with the observation, and it is itself a two-stage unit that a single caught_by label would code
either way. It sits outside the pilot's scope (PREDICTIONS.md only), so it could serve as a held-out probe, though not as the
positive control item 3 asks for (that needs another program).

## 4. Suggestion for the re-code (gates on dp, as the proposal does)

A cross-vendor rater on the current single-label codebook would measure vendor variance *stacked on* a codebook ambiguity,
and could not separate the two. Two cheap additions, neither needing a second vendor:
1. **Split Q2 into `noticed_by` and `settled_by`** (same channel list) and re-code the 13 read-vs-exec units with the
   existing two raters. If agreement on each field rises well above 13/25-disagreeing, the κ = 0.58 is mostly codebook-made.
2. **Pre-register how the thesis is scored on two-stage units.** "Reading changes only the record" should be tested on
   `settled_by` for world-facing VERDICT units. Then a world-facing correction that reading *noticed* but execution *settled*
   supports it, rather than counting against it or hiding in disagreement.

The cross-vendor rater stays worth doing. On the two-field codebook it measures the thing it is meant to measure.

## So what

This touches neither bucket. It bears on the one output every recent outside persona says is citable. Its weakest number
(κ = 0.58) has an identifiable cause that a vendor change would not remove. Not computed here: whether the two-field
re-code actually raises agreement. That is the cheap falsifier for this note, and if it doesn't, §2's reading is wrong.
