# Publisher Activity — 2026-08-03

**Mode**: **AUTONOMOUS** (cron, RUN-ID 18472, header 03:30:02 PDT = 10:30:02 UTC — and see below,
the "manual" labels on 08-01 and 08-02 were wrong)
**Window scanned**: 2026-08-02 05:00 → 2026-08-03 10:40 (Synchronism, synchronism-site, web4)

## Phase 0

- **CGL Arc** — `active` → **`closed_killed`** (2026-08-02). K1 KILL on the repaired instrument
  (best 0.042 vs 0.10 bar, D anchor 0.875); arc stopped per its own charter rather than declaring
  INDETERMINATE. My 08-02 exploration was read, cited by commit, and acted on: the lane built a
  CGL-native defect metric and re-swept. Theory ledger empty; frame ledger three items.
- **Alignment Arc** — **new** `upcoming_candidate`, `active`, explicitly **not** a candidate
  (zero executed stages). First arc whose frame is owned by dp. Q1(b) registers superdeterminism
  as a kill in advance, which forecloses this literature's standard escape hatch before the arc
  can reach for it.
- **REC-2026-040 OPENED — 0.45 / MEDIUM.** The α-admixture bound (≤25% weight on a local-density
  variable, 95%) + the smoothing scan's 1/r² obligation, as a **framework-independent** methods
  note. Only such item in the file. Opened low on purpose: external prior-art gate unrun, named as
  weakness #1 — the direct correction to opening REC-039 at 0.72 on 07-29 without one.
- **REC-2026-039** — weakness added, **held 0.38**. Superseded as lead instance of its own genre;
  the α-bound needs no fitted ceiling and no M/L conditioning. Inverts its own framing field's
  stated ordering.
- **REC-2026-038** — strength + weakness added, **held 0.93** (2nd consecutive deliberate hold).
  Fifth instance / fifth failure mode: continuity that exists, is correct, is durable, and is
  never consulted. Recommended the draft stop enumerating and write the asymmetry.
- **REC-2026-035** — weakness added, **held 0.95**. Open *question* (not claim) against its own
  blocker: is "BIG-SPARC required" about sample size or about needing a per-galaxy radial
  dimension a BTFR point cannot supply at any N? If the latter, its framing sentence is unverified.
- REC-034 / 036 / 037 — held (0.97 / 0.68 / 0.92), no new evidence.
- Null #1's Reading A/B fork — **7th consecutive day**; partially overtaken by the ≤25% bound.

## Phase 1 — Synchronism

**Changed.** `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — status note
under "MOND-Synchronism Unification (Sessions #86-89)". That heading's standing sentence, *"a
breakthrough discovery: MOND and Synchronism are the same physics in different parameterizations,"*
is now an **exact algebraic identity** at γ = 1/2 and therefore no longer a discovery. Verified
here numerically (max dev 5.55e-17) before propagating. Paired with the ≤25% admixture bound, the
γ↔ρ_crit degeneracy caveat, and the EFE retraction (→ EFE = 0 exactly, in tension with Chae+2020).
Refutation count **held at 6**; Bucket 0 unchanged. Logged in `05-quantum-macro/meta/CHANGELOG.md`
and `PUBLISHER_CONTEXT.md` §6.

**Scoping correction filed against the source finding** (mine, not carried there): its §4 decisive
null is a *quantification* of Lelli+2017 prior art at constant scale height, per this program's own
07-23 prior-art audit row 13. Novel content is §6 + §7.

**Gates**: claims freeze ✅ (10 claims) · `make-md.sh` ✅ exit 0 · churn gate **content 40 / raw
23,168** → artifacts restored, not committed (3rd consecutive correct firing) · lone-CR: one path,
the frozen `claims/v1-snapshot/`, 10th consecutive pass.

## Phase 1 — Web4

**No change.** `60926fa..e5f221a`, 3 commits, zero whitepaper-scope. C310 reports the t3-v3-tensors
spec byte-frozen and clean for the 8th consecutive pass, zero mutation, findings routed. No
canonical-terminology exposure.

## The correction that is the second finding

**The cron never failed.** 08-01 log 3,514 B ending `exit=0`; 08-02 log 4,704 B ending `exit=0`.
Both autonomous. I reported both as deaths across four surfaces and escalated a spurious
OWNER-ACTION for a `2>&1` present since 07-24. Mechanism: 03:30 PDT **is** 10:30 UTC, so the agent
reads its own run's opening banner as a prior failure — and because `claude -p | Add-Content`
flushes only at completion, the log is header-only *whenever you look*, so the test can only return
"dead." I had caught this on **07-31** and written the correct rule into the collective log; it
reached nobody, including me, twice. Correction banners appended to the 08-01 and 08-02 reports,
the 08-02 activity log, `PUBLISHER_CONTEXT.md` §6, and the collective log.

Paired with 08-02 this measures the asymmetry both directions: one false clause → four documents,
two repos, <24 h, unaided. One true correction → nobody, 72 h, twice.

## Unsolicited

- `explorations/2026-08-03-publisher-one-substitution-and-its-exclusion.md` — back-annotates the
  sibling-repo RAR scatter no-go (absent from this repo entirely), names the 48-minute convergence,
  and files the prior-art scoping correction.
- Flagged for the site lane: propagation recommendation #2's novelty justification contradicts its
  own §2 and Scope block.
- Flagged for the supervisor: **withdraw** the Publisher's 08-01 and 08-02 cron OWNER-ACTION.

## Artifacts

| Path | Action |
|---|---|
| `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` | modified |
| `whitepaper/sections/05-quantum-macro/meta/CHANGELOG.md` | appended |
| `whitepaper/PUBLISHER_CONTEXT.md` | §6 entry + correction to 08-01/08-02 entries |
| `explorations/2026-08-03-publisher-one-substitution-and-its-exclusion.md` | new |
| `manuscripts/publisher/state/recommendations.json` | 55 ins / 9 del, no churn |
| `manuscripts/publisher/state/whitepaper_sync.json` | updated (was stale at 07-30) |
| `manuscripts/publisher/reports/2026-08-03-publisher-report.md` | new |
| `manuscripts/publisher/reports/2026-08-0{1,2}-publisher-report.md` | correction banners |
| `manuscripts/publisher/logs/publisher-activity-2026-08-02.md` | correction banner |
