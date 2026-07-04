# Publisher Daily Report - 2026-07-03

## HOLD (arc AT REST) — substantial catch-up; preprint candidate #1 materially STRENGTHENED; still pending dp

Recovery run. The 07-03 02:30 autonomous cron failed silently (40s, no report — flagged in commit `52ceb99f`); this run supersedes it and catches up the 07-02→07-03 window. Still **not a reopening** — no dp go-signal, no new numbered Session (core still S691); all activity is explorer/maintainer/triage work within the rest. But several items are directly preprint-relevant and one **materially strengthens the #1 candidate**, so this gets a fuller evaluation than a heartbeat.

### Items logged, grouped by preprint candidate they affect

**Candidate #1 — Locality No-Go (STRENGTHENED on two fronts):**
- **rho_crit(V) sign-inversion** (`09c1cfb9` + triage ACCEPT `bc66146d`; Archivist accepted as an epistemic update, confidence 0.9). MOND-matching *forces* ρ_crit ∝ **V⁻²** (from BTFR + a₀, profile-independent, verified two ways); the framework's C(ρ) asserts **V⁺²** — **wrong sign, unfixable by re-tuning.** This sharpens the no-go from a BIC/offset argument (ΔBIC=+184, ~1.7 dex) to an **exponent-level** structural refutation.
- **AeST positioning** (`locality_nogo_aest_positioning.md`, triage ACCEPT `bebc038a`). Positions the no-go against **AeST** (Skordis & Złośnik, PRL 127, 161302 (2021); arXiv:2109.04165) — the published relativistic-MOND theory that passes CMB+LSS. AeST is the **existence proof** that the no-go's non-local escape routes (g_bar / Σ / enclosed-mass keying) are viable, converting the locality triage table "from taxonomy into a statement about the actual 2026 frontier: the local shortcut is quantitatively dead; the surviving program is the AeST-class completions." Also sharpens the target audience (relativistic-MOND readers; cite Milgrom 2005 + Lelli 2017 + AeST).

**Candidate #2 — Mechanism-Class DESI:**
- **TEST-04a direction-overclaim correction** (`cd6fdadb`, 2026-07-02). A site self-correction on the TEST-04a growth-suppression direction claim — consistent with the arc's self-correction discipline (cf. dim-4 LIV 07-01). Tightens the honest framing of the DESI candidate; no readiness change.

**Candidate #3 — A2ACW Program-Level Null:**
- **Top-3 site "research contributions" demoted (9/9 base rate)** and a new A2ACW failure-taxonomy clause, **"retraction-survival-in-compilation-layers"** (`136225f5`): the S606 "below CDM" claim outlived its S610 retraction on the site — a retraction that failed to propagate. This is a concrete, self-documented instance of a methodology failure mode → strengthens the A2ACW methodology paper's empirical base (the program catches its own non-propagated retractions).

### Publisher assessment update

The three-candidate ranking is unchanged in order but **#1 is now materially stronger and better-positioned**: an exponent-level (sign, unfixable) refutation plus positioning against the leading relativistic-MOND theory. If dp approves packaging, the Locality No-Go note can now claim more than a possibility statement — it lands squarely in the relativistic-MOND working literature with AeST as the worked escape-class member. Revised advisory framing: **#1 (Locality) is now co-lead with #3 (A2ACW)** for "most ready + best-targeted"; **#2 (DESI) still trails** on the single-bin ~2σ caveat.

None of this moves readiness: these are sharper refutations / methodology instances within the already-characterized negative-results body, not new confirmed novel predictions. The 0.99 lever remains gated on actual external preprint/verification.

### Pending dp decision (carried, unchanged)

Package the three transferable nulls as external preprints vs. continue daily auditing. **dp's go = the reopening signal** → Publisher Phase-2 drafts the three-preprint outline. The accumulating strengthening (rho_crit exponent, AeST positioning) is decision-relevant input for dp: the case for packaging #1 is stronger today than on 06-29.

### Discipline note

Per the rest: recorded in dp-facing channels (report + logs), `last_updated` bumped for Supervisor health, no `recommendations.json` narrative churn, triage-work ≠ arc reopening, readiness not flipped.

- **Readiness HELD**: REC-037 0.98 / REC-034 0.97 / REC-035 0.95 / REC-036 0.60.
- **Window (none publication-relevant)**: SAGE fleet healthy, all 5 active +4 (sprout S397 / mcnugget S298 / cbp S110 / nomad-e2b S092 / hub S015 sensing). **thor push-gap 5th window** (S208+ likely unpushed; owner nudge stands). **legion-gemma4-e4b grounding loop now 74 days** (~78 Session-0 commits, sessions/ empty — OWNER-ACTION). web4 C129 mrh-tensors clean re-audit; S213 unchanged. chemistry S2675 (governance gap ~Day 60) / gnosis S19 unchanged.
- **Whitepapers**: both Current (verified 07-01 + 07-03 no-change passes); no proposals; no terminology drift.

**Reopen trigger (full Publisher engagement / Phase 2)**: dp's go on the preprint-packaging proposal — OR fleet agent-ensemble transfer-bet data / new data / fresh lens.
