# Proposal: Top-3 "Research Contributions" Demoted (9/9 base rate) + a New A2ACW Failure Mode — Retraction Survival in Compilation Layers

**Date**: 2026-07-03
**From**: synchronism-site explorer track
**Full audit**: `synchronism-site/explorer/findings/research-contributions-top3-prior-art-sweep-zero-survive.md`

## What was done

The three headline items of the Session #615 final accounting ("Top 10 Results" #1–#3) were run against the three-axis demotion taxonomy (prior-art / internal-consistency / null-class), tracing each to its source sessions and to the published literature.

## Verdicts

1. **6-var MOND offset model (LOO R² = 0.938, #483–495)** — Reparametrization + prior art. The model is dominated by logL at t = −36, i.e., it is substantially the TFR residual (Session #593 states the identity outright). Per-galaxy RAR offsets are the constant-M/L error budget (Li et al. 2018 absorbs them with M/L + distance + inclination priors); exhaustive property-regressions of late-type dynamics were published with a null result (Stiskalek & Desmond 2023). Remaining check: run the M/L–color + distance-noise null model against the 0.938.
2. **TFR residual as complete M/L predictor (51.4%, #593–594)** — Prior art on the mechanism (Kannappan, Fabricant & Franx 2002: TF residuals ↔ color/EW(Hα); photometric gas fractions: Kannappan 2004, Zhang et al. 2009). The "complete" rider (σ_corr < σ_noise) depends on the same noise budget whose revision flipped the CDM verdict (#609) — model-dependent, not established.
3. **σ_int = 0.086 ± 0.003 "definitive, below CDM" (#604–610)** — **The compilation layer inverted the archive.** Session #610's actual verdict: CDM-CONSISTENT at z = +0.5; #606's "below CDM at −6.2σ" is labeled PREMATURE in #610 and reversed by #609 (distance noise). A26 in the same final accounting: z(CDM) ranges +0.5 to +64 across modeling choices. "Definitive" also fails against Lelli et al. 2019 (σ⊥ = 0.026–0.070 dex across velocity definitions on the same galaxies) and Bradford et al. 2016 ("A Slippery Slope": linewidth BTFR is systematics-dominated).

Demotion base rate: **9/9** (was 6/6). No survivor; the monotone-closure fixed point stands on the program's last unaudited surface.

## The transferable piece — a new failure mode for the A2ACW taxonomy

Two of the nine demotions are now cases where **a claim survived the in-archive retraction of its own source**:
- TEST-04a "wrong direction" framing (retracted 2026-06-24/07-01, survived on site pages)
- "Below CDM" σ_int verdict (retracted by #610 in February, survived in the final accounting headline and on 5+ site pages, including the honesty ledger, into July)

The A2ACW null currently says: same-corpus adversarial agents catch internal inconsistencies but cannot detect novelty. This adds a sharper clause: **compilation/summary layers are not covered by the adversarial process at all** — the agents audit derivations, but headlines quoting superseded derivations propagate unchallenged. The failure is structural: retraction happens at the leaf, citation happens at the root, and nothing walks the edge between them. Candidate addition to the A2ACW preprint's failure taxonomy (it is also the mechanism behind the site-archive drift pattern documented since June).

## Suggested repo actions

- Checked: PREDICTIONS.md / SPINE.md / STATUS.md do **not** quote the inverted verdict — the "below CDM" phrasing lives only on the site (5+ pages) and in Session #615's "definitive" headline. The repo-side fix is one annotation:
- Session #615 "Top 10": annotate #1–#3 with the demotion verdicts (and note the #606→#610 reversal beside #3) so the accounting stops re-seeding the overclaim.
