# Publisher Activity Log — 2026-08-07

**Mode**: autonomous
**Result**: whitepaper UPDATED (self-correction), 0 new recommendations, 2 recs updated (both HELD)

---

## Sequence

1. Read `publisher/CLAUDE.md`, archivist log (08-07 09:30), collective publisher log (08-06).
2. Scanned all §1b surfaces **including** `synchronism-site/explorer/findings/` — which is how the finding
   surfaced. New since last run: 1 sibling finding (08-05 08:12), 1 exploration (hive-organs arc charter).
3. Read the sibling finding; recognised it answers §6.4's discharge condition, registered 08-06 as "unrun".
4. **Re-derived independently before accepting** — `simulations/publisher_20260807_ell_consistency_check.py`,
   5 checks. Did not take the sibling's result at face value; Q4 and Q5 are net-new beyond it.
5. Edited whitepaper sources (2 sections, 2 CHANGELOGs, PUBLISHER_CONTEXT §6a).
6. Updated `recommendations.json` (REC-038 +1S/+1W, REC-036 +1W). Schema read from the file, not memory.
7. Gates: build, claims freeze, churn, lone-CR, json churn.
8. Report, this log, collective log, commit.

## Verification performed

| Check | Method | Result |
|---|---|---|
| A depends only on β_J·R₀ | 5 factorizations at fixed product | all → 0.0023429 exactly |
| 645.0 closed form | direct | 644.9987 ✓ |
| 08-05 proposal's 317 pc vs invariant | 1×0.3174 vs 4.5×0.07 | 0.76% — tautological |
| §6.4 point 3 reproduces | invert C at both endpoints | x ratio 6060 vs (80)²=6400 ✓ (rounding) |
| Two-sided convention | derived + numeric | `x = 0.019·β_J²·[V_c/V]²`, ℓ cancels |
| ρ-only convention | top-hat integration over exponential disk | x 3.58 → 0.0045, **opposite sign** |
| 1.515× origin | (3.2/2.6)² | 1.5148 — R_d, not h |
| NGC 3198 R_d ground truth | SPARC Lelli+2016 Table 1, in-repo | **3.14 kpc** (2.6 is 21% low) |

**Data provenance note**: `simulations/sparc_data_cache/NGC3198.dat` is headed *"Synthetic data for
Synchronism validation"* (R_disk = 4.00 kpc) and was **not** used. Real values come from
`simulations/sparc_real_data/SPARC_Lelli2016c.mrt`. Recording this because the two directories are adjacent
and the synthetic one is the shorter path.

## Gates

- build exit 0 (7,515 lines, +19 from 7,496)
- claims freeze exit 0 (10 claims, v1 verified)
- churn **content 50 / raw 23,406 → FIRED** (8th consecutive), artifacts restored via
  `git checkout -- docs/whitepaper/ whitepaper/build/`, CI is the builder
- lone-CR: 1 text path (`claims/v1-snapshot/Synchronism_Whitepaper_Complete.md`, frozen, unchanged)
- recommendations.json: raw 9 == content 9, `ensure_ascii=True`, trailing newline, round-trip verified
  (40 recs / 170 milestones)

## Self-caught this pass

1. **The primary finding is itself a self-catch of yesterday.** Nothing external prompted it beyond scanning
   a surface my own rules already required.
2. **Wrote "~800×" for the cross-convention spread, then computed it.** The real figure at a common
   ℓ = 100 pc is **3.5×10⁴** — conventions (a) and (c) coincide to ~30% there by numerical accident, which is
   what made the eyeballed estimate feel plausible. Fixed before commit, and Q5 was added to the script so the
   whitepaper number is reproducible from the cited file rather than asserted. *A spread quoted from
   inspection is an assertion; the standing rule is compute it or cut it.*
3. **Nearly cited the synthetic SPARC cache** for NGC 3198's R_d (it lists 4.00 kpc, which would have made the
   whitepaper's 2.6 kpc look merely conservative rather than wrong). Caught by reading the file header.

## Not done / deferred

- The REC-038 base-rate measurement (how often prior art IS found and cited on time) — flagged again, still
  unowned. It is the single thing between REC-038 and 0.95+, and it has been standing since 08-03.
- `.gitattributes` root fix for the churn class — dp's call, still open, 8th firing.
- web4 Track H (`85975dc`) — re-examine after the 2026-08-08 hackathon.
