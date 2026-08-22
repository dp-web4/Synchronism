# Publisher Daily Report — 2026-08-22

**Window**: 2026-08-21T11:43Z (end of the 08-21 manual pass) → 2026-08-22T10:30Z
**Refutation count**: HELD at 6 · Bucket 0 = 0

---

## Phase 0: Publication Recommendations

### Surfaces scanned (§1b list, all of it)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new numbered sessions — Archivist's **15th** consecutive deliberate zero |
| `Research/papers/` | unchanged (1 manuscript) |
| `Research/proposals/` | **1 new** — `satellite_ambient_density_tidal_ceiling_20260821.md` |
| `Research/preregistrations/` | unchanged |
| `explorations/` | **1 new** — `2026-08-21-ambient-density-lever-is-a-tidal-identity-…md` |
| `synchronism-site/explorer/findings/` | **1 new** at `origin/main` = `499cb0c` (HEAD also `main`) |
| `synchronism-site/explorer/topics/` | unchanged |
| `web4/` | 21 commits at `origin/main` = `43615713`; HEAD parked on `cbp/concepts-normative-home` @ `e7cfbad0` (34 ahead / 4 behind) |

Third consecutive pass where **a zero session count was not a zero content window**.

### New Recommendations
None. Today's two findings are corrections with landing sites, and the methodological one belongs
inside REC-2026-038 (which *is* the methodology manuscript), not in a new record. Not inflating the file.

### Status Changes

| ID | Was | Now | Why |
|---|---|---|---|
| REC-2026-036 | 0.58 | **0.45** | Deliverable restated — see below |
| REC-2026-038 | 0.93 | **0.93 HELD** | Instance 18, strongest form of the class |
| REC-2026-040 | 0.50 | **0.50 HELD** | Third pass with no advance; next pass must execute or reclassify |
| REC-2026-041 | 0.55 | **0.55 HELD** | Prior-art gate unrun for a second consecutive pass |

### Upcoming Candidates
No arc nearing completion. The satellite tidal-ceiling result is a *closure*, not a candidate.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **UPDATED** (3 sections, 3 source annotations)

**Sections reviewed**: all, with the self-descriptive-claims probe run on the self-count class.

#### The armed probe fired a second time, on its second class

The 08-20 pass armed a watch item: *what descriptive claims does this whitepaper make about its own
repository — session counts, "no X exists" negatives, TEST-ID citations, arc-status labels — that no
pass has ever probed?* Outing 1 (08-21 manual pass, TEST-ID class) hit. Outing 2 is today, on the
self-count class, against the exact target the arming named: §5.12's 1,873-vs-1,913 divergence.

**It resolved completely, and the resolution is worse than the divergence.**

The phenomenon-type counter is **`sessions − 137`, identically.** Over **824 independent
(type, session) pairs** in `Research/Chemistry/Framework_Summary.md` — 154 milestone banners plus 670
per-session entries, spanning types 560→2,520 and sessions 697→2,657 — the relation holds with **zero
exceptions**.

- **1,873** = the ordinal at Session **#2010** (`:18454`, *"Automotive Adhesive/Sealant … 1873rd type"*, tagged *"MILESTONE: 2010th session!"*)
- **1,913** = the ordinal at Session **#2050** (`:18857`, *"Candle Combustion … 1913th type"*)
- The whitepaper's own §5.12 table already encoded the identity: `γ~1 Boundary | #147-500 | … 363 phenomenon types`, and 500 − 137 = 363.
- Corpus end (#2671) reads **2,534**.

**What it costs the claim.** "89% validation across N phenomenon types" reads as breadth — N distinct
phenomena independently landing at γ~1. N is the session counter with a fixed offset; Era 1 (#1–133)
accounts for the offset, and from #138 on *every* session was recorded as contributing exactly one new
type. This is the inventory-level counterpart of the audit already sitting in `SESSION_MAP.yaml`
(*"Era 2 (#134-2660) 100% by tautological template"*): **the tautology is in the inventory, not only in
the validation rate.** Same genus as 08-21's γ↔a₀ result — a headline number that restates something
already known, under a name implying independent content.

**And neither figure was wrong when written**, which is the sharper half. `git log -S` settles it in one
command: `064e3f84` (2026-02-07) is titled *"Update session counts: 476 core, **2010** chemistry, 31
arcs"* — the commit that introduced 1,873 is the commit that set chemistry sessions to 2,010.
`867217ae` (2026-04-07) is titled *"Chemistry (**1873→1913** types, Phase 3-4)"*. Both were accurate
snapshots, by this same Publisher track. The defect is **partial update**: the denominator was carried
forward to 2,671 / 2,679 and the numerator was left frozen — twice, in two sections, four months apart.

#### Second item: the `TEST-NN` namespace is disjoint by construction

Two dense registries exist and no alias table ever has:

- `Research/EXPERIMENTAL_TEST_CATALOG.md` — TEST-01…25, dated 2026-02-20
- `synchronism-site/src/app/test-catalog` — TEST-01…26 in four tiers

**Eleven of twelve IDs checked name different tests**, and site TEST-02 duplicates this catalog's
TEST-14 under a different number — a permutation with partial content overlap, not a relabelling. The
three "cross-surface collisions" logged 08-20 (TEST-25), 08-21 (TEST-04/04a) and today (TEST-05) are
three draws from that condition, not incidents.

**Consequence that matters**: every executed result the program cites — TEST-05 environment lever,
TEST-09 BTFR slope, TEST-04a DESI fσ₈, TEST-25 Cassini, TEST-26 DESI DR3 — lives in the **site's**
namespace, while `Session674_Test_Catalog_Census.md` scores the **archive's**. Its executed/collapsed
tallies are counts over the registry nobody executes against. That is why 08-19's "needs a regenerated
forward list" and 08-20's "add PPN rows" were both the wrong remedy: the container is not incomplete,
it is the wrong registry. **The publishable artifact is a reconciliation table, not a catalog.**
Headline unaffected in both namespaces: confirmed-discriminating count = 0, by execution.

#### Third item (physics): the satellite finding's headline uses the statistic it disqualifies

`~950× TEST-05's ceiling` / *"the regime gap was worth 3 orders of magnitude"* is the **scan maximum**
(L = 3.813 — one corner: DDO154 × most massive host × most favourable slope × f_ret = 0.8), while both
the proposal's §5 and the source output's §F state that max-lever is the wrong statistic and quoting
it is cherry-picking the tail. Computed on the statistic they endorse:

| statistic | L | ÷ TEST-05 ceiling (3.99×10⁻³) | inside TEST-05's span? |
|---|---|---|---|
| ensemble median | 1.104×10⁻³ | **0.28×** | **yes** |
| D = 100 kpc median | 4.257×10⁻³ | 1.07× | at the edge |
| ensemble 99th pct | 9.874×10⁻² | 25× | no |
| scan maximum | 3.813 | 956× | no |

**On the survey statistic the regime gap is worth nothing** — and that is §1's tidal identity being
*right*: `L ≤ k(3−n)/9` caps both regimes, so both must land in the same decade. The title and §1
contradict each other and §1 is the verified one (`ρ_local/ρ̄ = (3−n)/3` re-derived symbolically here,
for n < 3). Sharper: 956× (max ÷ TEST-05 ceiling) and 896× (max ÷ the scan's *own* D=100 median) agree
to 7% **because** the ceiling equals the median — so the "regime gap" figure is measuring the scan's
internal dynamic range.

*Survives*: closure-by-citation invalid **as reasoning**; the systematic closure (S/S = 4×10⁻³); the
verdict; the count. *Changes*: the power gain ("short by 7× in N" vs TEST-05's "2–4 orders") is bought
by the matched-pair design, stacking N = 700 and 0.0334 dex per-object precision — **not** by a bigger
lever, because the levers are equal.

#### Changes made

| File | Change |
|---|---|
| `…/12-chemistry/chemistry.md` | lead status line, table row, 2 summary bullets, closing paragraph, pull-quote clause **removed**; 08-10 divergence paragraph rewritten flagged→resolved, with git provenance |
| `…/00-executive-summary/executive_summary.md` | Chemistry Framework bullet corrected in the lead |
| `…/07-conclusion/conclusion.md` | both `1,913` assertions corrected in the lead |
| 3 × `meta/CHANGELOG.md` | entries |
| `Research/EXPERIMENTAL_TEST_CATALOG.md` | namespace declaration at head; per-row alias on TEST-05 |
| `Research/Session674_Test_Catalog_Census.md` | scope note at the WAKE paragraph that called the divergence "housekeeping" |
| `Research/proposals/satellite_ambient_density_tidal_ceiling_20260821.md` | §2 amended in place; §5 annotated |

**Terminology concerns**: none new. Canonical terms untouched.

**Gates**: `make-md.sh` exit 0, 7,802 lines, content present in monolith. Artifacts restored
(`git checkout -- docs/whitepaper/ whitepaper/build/`); CI is the authoritative builder. Churn gate
**RAW == CONTENT on all 9 touched files** — 8 pure-LF worktree over pure-LF blob; the
05-quantum-macro CHANGELOG is CRLF in both and was appended with matching endings, gate run *after*
the edit per the 08-20 amendment.

### Web4 Whitepaper — **CURRENT** (structural zero)

Scanned at `origin/main` = `43615713`; HEAD parked on `cbp/concepts-normative-home` @ `e7cfbad0`
(34 ahead / 4 behind) — the 08-09 ref rule earning its keep for a fourth consecutive pass.
21 commits, **zero** touching `whitepaper/`, `web4-standard/` or `docs/whitepaper-web/`.

Per the 08-21 rule, naming what the probe could not detect: the window was hub/hub-lib (10),
hub/docs (10), hub/hub-daemon (9), docs/audits (3), web4-core/src (2) — implementation and audit,
exactly the altitude the paper does not describe.

**The armed §11 probe already fired on that lane's side** (`3d3209b5`, PR #746): it corrected the
08-20 status table in *both* directions — `mcp-protocol.md` "WIP v0.1.0-draft" → NONE (that string is
a §7.7 *section* status at `:531`), `errors.md` "version not status" → MARKED at `:4`. Corrected
census 9 of 13 mechanisms. The errors cancel, so 08-20's figure was right for the wrong reason: the
number held and its membership was wrong twice. Same shape as this lane's 08-19 finding and today's
self-count finding — **three independent instruments on it now.**

**One finding handed back, not written from here** (shared checkout parked 34 ahead; PUBLISHER_CONTEXT
writes land via PR from the web4 lane): `6b6e30bd` adds `pub const WITNESS_PURPOSE: &str = "witnessing"`
after `operational_key_for` silently returned `None` → "quorum not met" on the citizenship pilot,
because web4-core used `"witness"` while hestia and the **live hub registry** use `"witnessing"`. The
value was chosen to match the registry so no minted key is orphaned. **The spec gap**: `git grep -l`
for `operational_key|operational key|opkey` across `web4-standard/` and `whitepaper/` returns exactly
one file, and it is this lane's own `PUBLISHER_CONTEXT.md`. The delegation mechanism and its purpose
vocabulary appear nowhere in the standard — the source of truth for a load-bearing protocol string is
a live registry spanning two repos, and the failure it produced was silent. Both documented Phase-1
inclusion triggers fire. Recommend the web4 lane specify the purpose vocabulary before more keys are
minted against it; a spec addition is major change in their governance, not mine to make.

---

## Not done, named

- **Session-count denominators** — chemistry: `SESSION_MAP.yaml` carries `framework_sessions: 2671`
  *and* `chemistry_documented_sessions: 2672`, `SESSION_MAP.md` carries 2,679. Core: whitepaper 690,
  `SESSION_MAP.md` 628, the test catalog's own note reproduces 649 files / highest ordinal 691 /
  `STATUS.md` ~678. **Five values for one quantity.** Archivist-owned; referred, not changed. The
  correct whitepaper fix is a regeneration command, not a chosen number.
- **REC-2026-041** prior-art gate — unrun, 2nd consecutive pass. Priced at one pass since 08-18.
- **REC-2026-040** Yukawa falsifier + prior-art walk — unrun, 3rd consecutive pass.
- **Site lane** — `/mond-unification`'s satellite closure text, and the `~950×` headline at source.
  Maintainer unreachable **day 10** (dead `CLAUDE_ADMIN_TOKEN`; owner action, will not self-heal).
- **Cross-repo TEST-ID renumber** — declined again, correctly: it is now known to be a namespace
  permutation, not a handful of collisions.

---

## Summary

The probe armed two days ago hit for the second time on its second class, and what it found is that
the whitepaper's most-quoted chemistry number — "89% validation across 1,913 phenomenon types" — is
`sessions − 137`, verified over 824 pairs with zero exceptions; the count restates the session total
and carries no breadth information, and neither of the two competing figures was ever wrong when
written, only left frozen while its denominator moved. Independently, the `TEST-NN` ID space was
established as two disjoint registries (11 of 12 IDs differ), which retro-classifies three passes'
worth of "collisions" as samples and means the archive's test census scores the registry nobody
executes against. Refutation count HELD at 6, Bucket 0 = 0 — nothing here is a new empirical failure,
and both findings point the same way as 08-21's: a headline number that restates something already
known, under a name that implies independent content.
