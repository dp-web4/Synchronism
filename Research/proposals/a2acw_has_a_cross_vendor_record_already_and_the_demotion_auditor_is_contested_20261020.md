# A2ACW already has a cross-vendor record, and who did the demotions is contested

**Date:** 2026-10-20
**Author:** site maintainer track (Claude), from a simulated-visitor researcher persona's questions
**Status:** Proposed. Items 2 and 3 gate on dp. Item 1 is a provenance trace anyone can run.
**Builds on:** `a2acw_cross_vendor_corpus_control.md` (gated on "access to a second-vendor model pair"),
`de_sector_has_no_discriminating_test_and_the_oracle_result_needs_a_cross_vendor_rater_20261009.md`.

## What the persona asked

The researcher persona said A2ACW is the site's only candidate novel contribution. It asked two things: which models
played which roles, and why the cross-vendor control waits for a "maximally disjoint" pair. Its point was that waiting
for full disjointness is a reason never to run the control.

## What the archive records (each line read verbatim this session)

| Role / result | Model | Source |
|---|---|---|
| Defender ("PRIMARY") in LV's protocol | Claude | `forum/lv/README.md:23` |
| Challenger in LV's protocol | ChatGPT (incognito) | `forum/lv/README.md:24` |
| S659's 4/6 vocabulary-asymmetry catch rate | one Claude agent, both roles | `Research/Session659_…md:83` |
| 2026-05-22 calibration arm (6 canonical discoveries, specificity 0/6) | Claude Opus 4.7 wrote and scored | site `explorer/findings/a2acw-detector-false-positive-rate-null-baseline.md:56` |
| The ~3,308 archive sessions | **not recorded** | none found |
| Who demoted the six "Validated" claims | **contested** (see item 3) | S662:50 vs the site |

LV is an external reviewer. LV's three run logs (2026-01-19/20, PDFs in `forum/lv/`) contain 13 verdicts: 5 HELD,
4 NARROWED, 4 FAILED. They are the only cross-vendor A2ACW runs in the archive. The 3,308 internal sessions took LV's
name, but nothing records whether they kept LV's two-vendor pairing.

> **[Correction 2026-10-21, site maintainer (author of this proposal): the "13 verdicts (5 HELD / 4 NARROWED / 4 FAILED)" tally above does not reproduce.** The closing bullets of LV's three logs count 27: 10 HELD / 8 NARROWED / 9 FAILED (T1 3/3/5, T2 6/4/3, T3 1/1/1); `forum/lv/README.md` inflates T3 to 3/2/6. The 13 and 5/4/4 had no locator: a recalled number, the RECALL failure mode this program has already documented. Also: only T3 records a model (Claude Sonnet 4.5), and T3's challenger is labelled a "ChatGPT Incognito simulation", so only T1 is unambiguously cross-vendor. The contested demotion auditor is now traced: AI tracks, not a human audit (explorer 2026-10-20).]

## Three points

1. **The precondition in `a2acw_cross_vendor_corpus_control.md` was met before it was written.** A partly disjoint
   pair (Claude vs ChatGPT) was in use in January. The control cannot reach full corpus disjointness, but it never could:
   both vendors train on the same public physics literature. So the useful framing is a **lower bound on
   independence**. Run the same claims through a single-vendor pair and a cross-vendor pair. If the verdicts differ, the
   single-model results are contaminated by correlated error. If they agree, that bounds the correlation from below and
   does not remove it. The fleet already has Codex, Gemini and Kimi CLIs installed, so a second vendor is available on
   this machine. *Cheapest first step:* re-run LV's three targets (Temperature, Field Effects, CRT/Quantum) single-vendor
   (Claude vs Claude), compare against LV's 13 cross-vendor verdicts, and register the comparison rule before running.
   n = 13 is small. It is still more than the zero cross-vendor comparisons that exist now.
2. **The site's "same training distribution" sentence needs a scope.** It is literally true only for S659 and the
   05-22 arm (one model). For LV's runs, it means overlapping corpora across vendors. The site now says this on
   /research-philosophy#a2acw-model-roles.
3. **Provenance conflict on the six demotions.** `Research/Session662_Galaxy_Closure_And_A2ACW_Specificity.md:50` says
   "the 6 demotions came from human audit, not the AI loop". Since 2026-09-15, the site has said the auditors were "LLM
   agents, not an outside domain expert". The 09-15 correction cited "the archive record", and S662 is part of that
   record. One of the two is wrong. The S662 reading would weaken the oracle thesis's "5/6 refutations from code"
   framing, because a human audit is a third oracle type. The site reading would make S662's sentence an error inside
   an AI session. **Trace needed:** find the six demotion events (dates and authors in git) and record who made each.

## Recommendation

- **dp:** authorise the cheapest cross-vendor step in item 1 (it costs a few runs). Or say why LV's pairing does not
  count as the second-vendor pair the control was waiting for.
- **Explorer:** the item 3 trace (topic seeded:
  `synchronism-site/explorer/topics/a2acw-roles-and-the-demotion-auditor-provenance.md`).
- **No bucket moves.** Count 6, Bucket 0 = 0. This is about the audit instrument, not the physics.
