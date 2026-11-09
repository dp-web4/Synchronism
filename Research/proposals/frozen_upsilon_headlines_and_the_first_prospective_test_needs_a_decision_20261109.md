# Frozen-Υ headlines, and a first prospective test that needs a decision by 2 December

**Date:** 2026-11-09 · **From:** site maintainer · **Moves:** no bucket, count 6, range 4–8, Bucket 0 = 0
**For:** the project lead

## Two decisions, in order of urgency

### 1. Gaia DR4 is on 2 December, 23 days from now. Adopt something before then?
- `synchronism-site/maintainer/scripts/wb_dr4_reading_key_PREREG.md` (11-06): boost → γ.
- `synchronism-site/maintainer/scripts/wb_dr4_analysis_plan_PREREG.md` (today, draft): how a boost is measured, how outside analyses are read, and an outcome table fixed in advance.

**Recommendation: adopt both, or a revision, by about 25 November.** Neither can move Bucket 0 or the count. What adoption buys is that "0 prospective tests completed" can become 1, on data no one here has seen. The explorer's 11-08 run shows why the method has to come first: one error setting moves γ by a factor of 0.65, and the allowed contaminant fraction moves it from 1.8 to 1.0.

The plan is one AI session's draft, and DR4's column names and RV coverage are unchecked. Those items are marked. Boufourou 2026 proposes its own pre-registered DR4 protocol, which might be better to adopt or compare with; only the abstract has been read.

**Also yours: should Saad & Ting and Chae et al. be told?** The site now states the separation-error setting, attributed to the public code at `4333e9f` (I read line 106 myself). The mock results are one AI session's calculation. Writing to them is outward-facing, so I have not. Useful content for them: the frame mismatch, the v_sys term, the injection configuration, and a request for per-system PDFs and the consolidation script.

### 2. Does row #3 stay among the 6?
With one global Υ per law, the γ = 2 pin is +46 by point (95 % −26 to +132) and +2.7 by galaxy, about 2σ from a γ = 2 world. It was +184 / +10.9 at frozen Υ = 0.5. The reading rule fixed before the run says to print this beside every headline and to recommend "disfavoured, not refuted" under a Υ-free rule, and I have done that on the site. I have not changed the count. B2's registered text says "marginalized BIC", so the frozen number was never the registered statistic.

This joins the open counting-rule question (10-07, 10-13, 10-16, 10-18). If it is ruled out, the count is 5 with a floor unchanged. The range already carries it.

## What the personas found, and how much of it held

| Claim | Checked | Result |
|---|---|---|
| +184 → about +47 with own Υ | Registered run on the real table | **Held:** +46.2, +2.73 |
| H₀ closure conditional on γ | Same run | **Held:** 92.8 / 71.1 / 56.5 |
| Every γ is a MOND μ | By hand | **Held** |
| Boost 1.49–2.70 across Chae's bin | Script | **Held**, to the digit |
| C at peak slope, γ = 2, is 0.601 | Script | **Held** |
| GEDE form at γ = 1 | Algebra, numeric 10⁻¹¹ | **Held** |
| S&T code setting (explorer) | Read the line | **Held** |
| TEST-05 r² ≈ 0.014–0.14 under MOND+EFE | Existing explorer finding, 09-05 | **Overestimated:** r² ≈ 0.001 (0.036 with a perfect proxy) |

My own frame errors today: "a₀′/γ is identified" (10-04) was my lane's sentence, and it holds only with γ free.

## What this says about the program, as a researcher

1. **The strongest-looking refutation in the headline was a nuisance choice.** The +184 was printed and defended for months. The N_eff correction (09-29) moved it to +11, and the Υ correction moves that to +2.7. The explorer's 09-29 per-galaxy run already said about 2σ, and it sat beside the bigger number. This matches the 10-31 finding that refutation-side statistics are where the recent errors have been. An audit of every frozen-nuisance headline is seeded for the explorer. I do not expect +2843 to fall, since it survives galaxy-level scoring, but it has not been run with Υ free per law.

2. **Answers exist and do not reach the pages.** TEST-08's power against EFE was computed on 2026-09-05. The site went on saying "not contradicted" without the number, and a persona re-derived a rougher version today. Today I also nearly seeded it as an open question. Grepping findings first caught it, by chance of habit rather than by design. Findings have an "Action: Maintainer" drain, but a finding that answers a question nobody has queued has no route to the page that needs it.

3. **Most of today's HIGH items are bookkeeping, again.** One test on the calendar, one question nobody here can answer alone. A researcher persona's last unanswered question is the right one: "Who, other than an LLM, has checked any of this? Its most valuable next step might be a two-page note that an outside MOND researcher could check in an afternoon." The wide-binary result is the first thing in months that an outsider could check in an afternoon and would care about. That is a reason to take decision 1 seriously.

## Not done
- No ledger table row was edited (10-16 (i) gate). The block carries the corrections.
- The Bucket 3 "η ≡ Abrikosov–Gor'kov" row is flagged, not changed. I have not read AG's paper.
- Global Υ is a coarse model. Per-galaxy Υ is the explorer's 09-29 run.

Scripts: `synchronism-site/maintainer/scripts/gamma2_global_upsilon_per_law.py` (+ `_PREREG.md`, `_output.txt`), `gamma2_global_upsilon_twin_power_posthoc.py`, `efe_free_corner_boost_across_the_chae_bin.py`.
