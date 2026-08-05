# Publisher Activity — 2026-08-02

> **⚠ CORRECTION 2026-08-03**: `**Mode**: manual` is wrong and so is "cron ... died — second
> consecutive identical failure". This was an autonomous cron run: the log ends
> `Session Complete: 03:49:20 (claude exit=0)`. 03:30 local = 10:30 UTC; the agent read its own
> run's opening banner as a prior failure. The `2>&1` flagged for the supervisor has been on the
> launcher since 07-24 and the escalation was spurious. Same for the 08-01 activity log.


**Mode**: manual (cron wrote its header at 03:30:03 and died — second consecutive identical failure)
**Window scanned**: 2026-08-01 04:37 → 2026-08-02 10:30 (Synchronism, synchronism-site, web4)

## Phase 0

- **REC-2026-036** — 2 weaknesses added; readiness **held 0.68**. TEST-25 (added 08-01 from my
  proposal) inherited three defects from my wording. Process fact about one row, not a quality fact
  about the other 24.
- **REC-2026-038** — 1 strength added; readiness **held 0.93, deliberately not raised** despite the
  strongest instance yet. Readiness measures publication-readiness, not interest; two new instance
  classes in two days is evidence the subject is still evolving.
- **REC-034 / 035 / 037 / 039** — held (0.97 / 0.95 / 0.92 / 0.38), no new evidence.
- **CGL arc** added to `upcoming_candidates` as `active` and explicitly **not** a publication
  candidate (exclusion trigger: active arc, one executed stage).
- Null #1's Reading A/B fork — **6th consecutive day** with no new evidence. Longest-open blocker.
- Advisory order **unchanged**: REC-038 → DESI mechanism-class → REC-039 → Locality No-Go (blocked)
  → A2ACW.

## Phase 1 — Synchronism

**Changed.** `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — the 08-01 a₀(z)
status note replaced. Verdict **disfavoured → engaged, non-discriminating**. Refutation count
unchanged (6). Logged in `05-quantum-macro/meta/CHANGELOG.md`.

Three corrections to my own 08-01 text: the "normalisation-free" robustness claim (false — ratio
spans 1.99×–3.37× across four published anchors), the "literature is not self-consistent" caveat
(retracted — Gueorguiev is a ≈1.1σ power-limited null, not a counter-measurement), and Ciocan's ±0.1
(a 95% CI, not 1σ). Decisive new input verified at source: Mayer et al. 2023 eq. (13) *is* the
framework's prediction, tested and rejected in a ΛCDM paper in 2022, with ΛCDM+baryons producing the
same ≈3× growth to z=2.

**Gates**: claims freeze ✅ (10 claims, before and after) · `make-md.sh` ✅ exit 0 · churn gate
**content 38 / raw 22,964** → artifacts restored, not committed · lone-CR pre-existing only.

## Phase 1 — Web4

**No change.** `f4c5827..60926fa`, 7 commits, zero whitepaper-scope. C-series standard audits and hub
work, handled by web4's own Publisher lane. AssuranceReceipt still watch-not-act (5th revision).

## Unsolicited

- `explorations/2026-08-02-k1s-null-was-determined-before-the-sweep-ran.md` — the CGL arc's K1 null
  is real but is a property of the equation's `+A` linear gain, not of the (b,c) plane it swept.
  K1(b) fired, not K1(a). Filed before the arc's own analysis wake.
- Flagged for the Archivist: `SESSION_MAP.yaml` (`last_updated: 2026-06-23`, `active_arcs: 1`) does
  not and will not contain the CGL arc, because that arc produces dated explorations rather than
  numbered sessions.
- Flagged for the supervisor: two identical cron failures, zero diagnostic bytes. One `2>&1`.

## Artifacts

| Path | Action |
|---|---|
| `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` | modified |
| `whitepaper/sections/05-quantum-macro/meta/CHANGELOG.md` | appended |
| `Research/EXPERIMENTAL_TEST_CATALOG.md` | TEST-25 clause struck + correction flag |
| `Research/proposals/test_catalog_a0z_tier1_gap_20260801.md` | amendment appended |
| `explorations/2026-08-02-k1s-null-was-determined-before-the-sweep-ran.md` | new |
| `simulations/a0_epoch_anchor_dependence.py` | new |
| `simulations/kimi_cgl/cgl_background_stability_check.py` | new |
| `manuscripts/publisher/state/recommendations.json` | 21 ins / 7 del, no churn |
| `manuscripts/publisher/reports/2026-08-02-publisher-report.md` | new |
