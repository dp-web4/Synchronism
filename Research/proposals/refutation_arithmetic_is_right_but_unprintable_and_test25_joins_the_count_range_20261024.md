# Proposal: the refutation arithmetic is right but unprintable; TEST-25 joins the count range (3–8); 0.086 is not a BTFR scatter

**Date:** 2026-10-24
**From:** site maintainer (synchronism-site), from simulated-visitor graduate and researcher personas, and site explorer 2026-10-23
**Bucket moves:** none. Count 6 as counted; range on file becomes **3 to 8**. Bucket 0 = 0.
**Gates on the project lead:** (iii) table-row edits; (v) the frame question.

## (i) Two persona corrections to the galaxy-level ΔBIC are wrong on re-derivation. The site caused one of them.

Two personas flagged the galaxy-level ΔBIC as HIGH:
- **+184 → +10.9.** The researcher reconstructed it as pinned γ = 2 versus free γ (Δk = 1). On that reading only Δχ²
  should scale, and the ln N penalty should be recomputed, which gives ≈ +6.2, "not disfavoured".
- **+2843 → +163 / +142.** Dividing by 2,807, they read the labels as swapped.

Both readings fail against the scripts:

- `maintainer/scripts/gamma2_pin_galaxy_level.py` compares γ = 2 with **McGaugh's ν**. Each fits one acceleration scale,
  so Δk = 0. ΔBIC = N ln(RSS₂/RSS_McG) = 2807 × ln(0.1485²/0.1437²) = +184. The galaxy-level figure is
  N_gal × ln(RSS₂/RSS_McG) = 166 × 0.0656 = +10.9. There is no penalty to rescale. The 09-29 block already said "equal
  parameter counts, no penalty helps".
- `dbic2843_galaxy_level.py` runs on **N = 2,438 points in 122 galaxies**, not 2,807. 2843 × 140/2438 = 163 and
  2843 × 122/2438 = 142. The labels are right.
- The bootstrap [23, 351] is the spread of the full-N statistic. Scaled by 166/2807 it is [1.4, 20.8], which straddles 10.
  So "at the threshold" and the bootstrap agree.

**Where the misreading came from.** /galaxy-rotation said the galaxy-level +10.9 "compares γ = 2 with free γ using N = 166
galaxies (ln 166 = 5.1)". That is false, and it is the sentence the +6.2 grew from. It is now corrected and lint-retired. No
surface printed the comparator, Δk, N or the formula for either number.

**Recommendation (gates on the project lead, with the 10-16 table-lag item):** every ΔBIC in this ledger's tables should
carry four things: comparator, Δk, N, and the replication unit. The numbers were right for a month while three readers
could not reconstruct them. An unprintable refutation is audited by guesswork. Today two guesses would have moved the
count in the flattering direction (+6.2 drops row #3 outright).

## (ii) TEST-25's galaxy-level dependence belongs in the count range

The 09-29 block records the cascade:
- The ΔBIC ≤ 10 retained band is (0.425, 0.60) at N = 2807 but (0.3, 1.7) at N_gal.
- The 09-17 block finds Cassini Q₂ passes post-hoc at γ ≳ 1.5–2.
- So TEST-25's in-house empty intersection holds only under point independence.

The 10-18 range ("4 to 8") listed the γ = 2 pin's galaxy-level delta and omitted TEST-25's.

Under the galaxy-level convention both rows drop. The range is **3 to 8**. The two deltas are not independent: one
convention moves both. TEST-25's defensible content is still the published, nuisance-marginalized 8.7σ (DHF 2024; Park+2026
3–15σ). That result covers the family near γ ≈ ½, not the 1.5–1.7 tail.

Two galaxy-level γ intervals also coexist on the site. They are now labelled:
- **0.3–1.7**: the registered statistic, ΔBIC ≤ 10 at N_gal with nuisances frozen.
- **[0.43, 0.95]**: a 95% galaxy-bootstrap interval on γ̂, with Υ/D/i profiled (explorer 09-29). The same output also
  prints [0.430, 0.973] from the P4 bootstrap, so the upper edge is bootstrap-dependent at ±0.03.

## (iii) For the Session 610 row: 0.086 is an identity residual

Site explorer 2026-10-23 (PREREG site `b134d0f`; not re-run by me).
- Shuffled-V gives 0.0870, against 0.0864 for the real model, so velocity contributes 0.0006 dex. The target is a
  closed-form function of four of the five predictors: the Taylor+2011 colour formula to 0.0037 dex.
- Both of Session 610's verdicts have no content: "MOND rejected at +32σ" and "CDM-consistent at z = +0.5".
- Desmond 2017's table is in log V. The site's "0.061" reading was inverted.
- Like-for-like inverse BTFR: 0.073 dex in log V, against abundance matching's 0.064.

Back-annotated at the head of `Research/Session610_CDM_Synthesis.md`. Not a framework row; recorded for provenance, as in
the 10-09 (ii) block.

## (iv) Two clarifications, no verdict change

- **The local-fluid P(k) pin is two-sided.** A researcher persona argued that DE perturbations are sourced ∝ (1+w), so
  c_s² ~ 0.1 at w ≈ −0.99 would be invisible. That argument holds for a separate DE fluid. Here ρ_DE = f(ρ_m) comoves:
  matter and dark energy are one fluid, and the pressure acts on the matter perturbation itself. In the explorer's 09-14
  table, at k = 0.1, R is 6.09 at γ = 0.4999 and 0.0066 at γ = 0.5005. That is the Sandvik+2004 pair.
- **Joint galaxy + cosmology content: none, or already excluded.** This restates 09-24 (ii) and 10-15 (i) as one row in the
  site's distinct-content register.

## (v) Frame question: the visitor loop's top severity is now bookkeeping, and today it was wrong twice

The researcher persona's verdict was: "the physics question is answered; the research question is at the wrong scale; the
A2ACW positive control is the experiment worth running." I concur on the first two clauses, with a datum.

For six weeks the highest-severity visitor items have been refutation bookkeeping: N_eff, count rules, intervals, CI
levels. Each has moved a range or a label, never a bucket. Today the two HIGH items, on re-derivation, were the persona
misreading an unprinted comparator. The registry has no branch that moves Bucket 0 (09-25, 09-27, 10-09). Further
bookkeeping audits have diminishing returns unless the record is made printable, as (i) proposes. Done once, that also
stops the loop re-deriving the same numbers.

The two experiments with unknown outcomes are both routed and both still waiting:
- the A2ACW positive control / post-cutoff arm (09-10, 10-12);
- the stopping question (10-04).

Recommendation: **rule on the stopping question before another month of bookkeeping**, or name a third live experiment.
Absent a ruling, the site loop will keep producing range adjustments to a count whose rule is unchosen.

## Artifacts
- Site: /galaxy-rotation#dbic-formula; /dark-matter (+2843 formula); /honest-assessment (beginner box; range 3–8);
  `src/lib/ledger.ts` COUNT_RULE_SENTENCE; /for-researchers distinct-content row; /dark-energy two-sided note;
  /cdm-discrimination#not-a-btfr-scatter; site_lint 240 rules (positive control: 31 hits on HEAD, 0 after).
- Scripts read, not re-run: `synchronism-site/maintainer/scripts/gamma2_pin_galaxy_level_output.txt`,
  `dbic2843_galaxy_level_output.txt`, `explorer/findings/scripts/gamma2_pin_nuisance_refit_sigint0.034_output.txt`.
