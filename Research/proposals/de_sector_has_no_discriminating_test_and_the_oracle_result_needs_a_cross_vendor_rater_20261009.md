# Proposal: the DE sector has no discriminating test left, and the program's one transferable result needs a cross-vendor rater

**Date:** 2026-10-09
**Source:** Site maintainer session (visitor Pass 4, leading-edge researcher persona; explorer finding 2026-10-08)
**Status:** Proposed. Items 1 and 3 gate on dp. Item 2 is a record correction, made today.
**Count:** 6 (unchanged). **Bucket 0:** 0 (unchanged). No bucket moves.

## 1. The dark-energy sector has no test with a framework-specific outcome

The site explorer pre-registered (site `f3e045e`) and executed a DR3 Asimov forecast on the substituted
DE background (the one-parameter curve through ΛCDM at γ = ½):

- The 2002 constant-w Cardassian (wCDM, same parameter count) fits every DR2-allowed point on the curve
  to Δχ² ≤ 0.42, in all 11 precision arms. MP-Cardassian n = 0 fits to ≤ 0.10.
- Δχ²(wCDM − locus)/Δχ²(ΛCDM − locus) = 0.044–0.047, constant across arms.
- Resolving the curve's shape from wCDM at 2σ needs a 9.3σ detection of w ≠ −1. On the curve that
  requires γ ≈ 0.30, which DR2 already excludes at Δχ² = 137.
- The curve's w′ today is −1.95(1+w₀), inside Caldwell & Linder's generic freezing band.

Finding: `synchronism-site/explorer/findings/the-de-freezing-locus-has-no-shape-win-branch-at-dr3-the-original-cardassian-fits-it-to-0p4-chi2.md`.

TEST-26's kill was already shared with ΛCDM (2026-09-27). Its win side is now shared with wCDM. With
TEST-04a's growth readings either withdrawn (S), ΛCDM to 0.2 % (F) or excluded (U) (2026-10-07), **no
registered DE test has an outcome that singles out the framework.** The sector is ΛCDM with a relabelling
that data at DR3 precision cannot resolve.

**Recommendation (gates on dp):** retire TEST-26 and TEST-04a as discriminating registrations. Keep them
listed as "tie / shared-kill only" so that a tie cannot later be read as a success. This joins the
10-04 stopping question: across galaxies (nested in MOND), wide binaries (a one-γ squeeze), growth and
DE, the registered physics program now has no branch that can move Bucket 0.

## 2. Record correction: the CDM benchmark in Session 610 was internal

Session 610's "CDM predicts 0.085 dex from halo-concentration scatter" has no external source. The site
called the measured BTFR scatter σ_int = 0.086 ± 0.003 dex "CDM-consistent (z = +0.5)" on that basis, and
the external check had been queued since 2026-07-10. First external number, read today (verbatim):
Desmond 2017 (MNRAS 472, L35; arXiv:1706.01017), halo abundance matching for SPARC-like samples, gives
mock BTFR scatter "∼0.25 dex", "3.6σ discrepant with the SPARC value of ∼0.11 dex", in baryonic mass;
with zero abundance-matching scatter the mean falls to 0.061.

So the only external ΛCDM figure checked is about 3× the internal 0.085. The site now reads the CDM
verdict as **suspended**, not reversed, because the 0.086 is on ALFALFA W50 widths after an in-sample
TFR-residual M/L correction, not SPARC V_flat. This is a CDM-versus-MOND question and does not touch the
framework's ledger. It is recorded because Session 610's "CDM-consistent" headline propagated to six
site pages on an unsourced benchmark, the same pattern as the 07-10 fabricated-consensus finding.

## 3. WAKE: the transferable result's main weakness is now cheap to remove

Every recent researcher persona (09-23, 10-08, 10-09) says the same thing: the physics is MOND + ΛCDM
where it survives, and the citable output is the **oracle observation**. About 3,300 adversarial AI
sessions produced zero world-facing refutations. The 92-unit H-oracle pilot found that reading changed
verdicts only about the record. All 13 empirical eliminations trace to an external measurement.

Its stated limits are n = 92 from one program, two AI raters with κ(caught-by) = 0.58 (rater A of
unknown provenance), no human or different-model-family rater, and no positive control. Today's
researcher asked the obvious question: "Has it been coded by a second, independent rater (human or a
different model family)?" It has not.

The 2026-07-07 proposal `a2acw_cross_vendor_corpus_control.md` asked for a cross-vendor *debate* pair
and gated on access to a second vendor. The rater version of that control is much cheaper. It needs no
debate, only the frozen codebook and the 92 units, which are already in-repo
(`synchronism-site/maintainer/scripts/oracle_trail_units.jsonl`, `oracle_trail_PREREG.md`, rater A codes beside them;
rater B in `synchronism-site/explorer/scripts/oracle_trail_codes_raterB.jsonl`). Access
also appears to exist now: this machine's agent inventory lists non-Anthropic CLIs (codex, gemini,
kimi). One cross-vendor rater would turn the one result outsiders want to cite from "two raters from
one model family" into a measured agreement across corpora.

**Recommendations (gate on dp, because the corpus would be sent to another vendor's service):**
1. Approve one cross-vendor re-code of the 92 units with the frozen codebook, pre-registered, reporting
   κ(caught-by) against each existing rater. Pre-state the bar: κ ≥ 0.40 against both, or the "never
   about the world" reading is reported as rater-dependent.
2. Before or with that re-code, run the positive control the pilot names: seed the unit list with known
   world-facing reading-caught corrections from another program, so the codebook's ability to detect
   one is measured, not assumed.
3. If the re-code holds, consider the oracle note (with its denominator) as the program's first
   write-up for an outside audience. It does not depend on the physics being right.

## So what?

Item 1 removes the last DE row that could have been read as a live bet. Item 3 points the remaining
effort at the one place a result could hold up regardless of the physics. The physics rows can now
only lose; the methodology row can still be strengthened or refuted by one cheap re-code.
