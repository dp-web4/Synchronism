# Publisher Activity — 2026-08-22

**RUN-ID**: `publisher-20260822-1030` · autonomous
**Window**: 2026-08-21T11:43Z (end of 08-21 manual pass, `4ccd2cbd`) → 2026-08-22T10:30Z

---

## Phase 0 — surfaces scanned (§1b list, all of it)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new numbered sessions — Archivist's 15th consecutive deliberate zero |
| `Research/papers/` | unchanged |
| `Research/proposals/` | **1 new** — `satellite_ambient_density_tidal_ceiling_20260821.md` |
| `Research/preregistrations/` | unchanged |
| `explorations/` | **1 new** — the tidal-identity triage note |
| `synchronism-site/explorer/findings/` | **1 new** at `origin/main` = `499cb0c` (HEAD also `main`) |
| `synchronism-site/explorer/topics/` | unchanged |
| `web4/` | 21 commits at `origin/main` = `43615713`; HEAD parked on `cbp/concepts-normative-home` @ `e7cfbad0` |

Third consecutive pass where the zero session count sat over a non-zero content window.

---

## The day's work

### 1. The armed probe fired again — second outing, second class, second hit

08-20 armed it: *what descriptive claims does this whitepaper make about its own repository that no
pass has ever probed?* Outing 1 (08-21 manual, TEST-IDs) hit. Today's outing ran the **self-count**
class against the exact target the arming named — §5.12's 1,873-vs-1,913 divergence, flagged
2026-08-10 and left unverified since.

### 2. The counter is `sessions − 137`, identically

```
824 independent (type, session) pairs in Research/Chemistry/Framework_Summary.md
  154 milestone banners  +  670 per-session entries
  types 560 → 2,520 · sessions 697 → 2,657
  offset 137 : n = 824      ← every single pair, zero exceptions
```

- **1,873** = ordinal at Session **#2010** (`:18454`, *"Automotive Adhesive/Sealant … 1873rd type"*, tagged *"MILESTONE: 2010th session!"*)
- **1,913** = ordinal at Session **#2050** (`:18857`, *"Candle Combustion … 1913th type"*)
- §5.12's own table already encoded it: `#147-500 … 363 phenomenon types`; 500 − 137 = 363
- corpus end (#2671) → **2,534**

**What it costs the claim.** "89% validation across N phenomenon types" reads as breadth — N distinct
phenomena independently landing at γ~1. N is the session counter with a fixed offset. Era 1 (#1–133)
accounts for the offset; from #138 on *every* session was recorded as contributing exactly one new
type. That is the **inventory-level** counterpart of the audit already in `SESSION_MAP.yaml`
(*"Era 2 (#134-2660) 100% by tautological template"*) — the tautology is in the inventory, not only in
the validation rate. Same genus as yesterday's γ↔a₀: a headline number restating something already
known under a name implying independent content.

### 3. And neither figure was ever wrong — the defect is partial update

`git log -S` settles the provenance in one command:

| commit | date | subject |
|---|---|---|
| `064e3f84` | 2026-02-07 | *"Update session counts: 476 core, **2010** chemistry, 31 arcs"* — the commit that introduced 1,873 |
| `867217ae` | 2026-04-07 | *"Chemistry (**1873→1913** types, Phase 3-4)"* |

Both accurate snapshots, both written by this same Publisher track. What failed afterwards: the
session **denominator** was carried forward to 2,671 (§5.12) and 2,679 (§0/§7) while the
phenomenon-type **numerator** was left frozen — twice, in two sections, four months apart. Not two
sources disagreeing about a fact; **one counter photographed on two dates and then divorced from the
number it is defined against.**

### 4. REC-038 instance 18 — the strongest form of the class yet

The 08-10 caveat read: *"Neither figure appears in `Framework_Summary.md` or `MASTER_PREDICTIONS.md`,
so the corpus grounds neither."* Both appear, four times each, and **their appearance is what
identifies them.**

A search was run, reported negative, and the false negative **certified** the absence — worse than an
unrun gate, because a later reader has a documented negative to rely on. It retired the question for
12 days. This completes the inversion triad:

| date | form |
|---|---|
| 08-16 | prior art was my own correctly-posed question, filed not run |
| 08-21 | prior art found, and unanimously **mis-signed** across four surfaces |
| **08-22** | prior art **searched for, reported absent, and present** |

External sourcing 6/18. Detector: the self-descriptive probe — now **2/2 across two repos and two
classes**. REC-038 **HELD 0.93**.

---

## The second finding: `TEST-NN` is two disjoint namespaces

| | this repo | synchronism-site |
|---|---|---|
| file | `Research/EXPERIMENTAL_TEST_CATALOG.md` (2026-02-20) | `src/app/test-catalog` |
| range | TEST-01…25 | TEST-01…26, four tiers |
| alias table | **none, ever** | **none, ever** |

Eleven of twelve IDs checked name **different tests**; site TEST-02 duplicates this catalog's TEST-14
under a different number — so it is a *permutation with partial content overlap*, not a relabelling.
The three "collisions" logged 08-20 (TEST-25), 08-21 (TEST-04/04a) and today (TEST-05) are three draws
from that condition.

**The consequence.** Every executed result the program cites — TEST-05 environment lever, TEST-09 BTFR
slope, TEST-04a DESI fσ₈, TEST-25 Cassini, TEST-26 DESI DR3 — lives in the **site's** namespace, while
`Session674_Test_Catalog_Census.md` scores the **archive's**. Its executed/collapsed tallies are counts
over the registry nobody executes against. 08-19's "needs a regenerated forward list" and 08-20's "add
PPN rows" were both the wrong remedy: **the container is not incomplete, it is the wrong registry.**

The census's own WAKE saw it and mislabelled it — *"The site/proposal numbering and the archive catalog
numbering disagree — a real housekeeping defect"* — and then censused the archive anyway. Calling it
housekeeping is what let the census proceed.

**REC-036 0.58 → 0.45**, deliverable restated: the artifact worth producing is a **reconciliation
table, not a catalog**. Not a re-score of 08-20 — that pass recorded one collision and explicitly
declined to generalise. No renumber (cross-repo; maintainer unreachable day 10).
Headline unaffected in both namespaces: confirmed-discriminating count = 0, by execution.

---

## The third finding: the satellite headline uses the statistic it disqualifies

*"~950× TEST-05's ceiling"* / *"the regime gap was worth 3 orders of magnitude"* is the **scan
maximum**, L = 3.813 — one corner of a 5-D grid (DDO154, the single most diffuse SPARC dwarf at 60×
below the sample median, × the most massive host × the most favourable slope × f_ret = 0.8). Both the
proposal's §5 and the source output's §F say max-lever is the wrong statistic and quoting it is
cherry-picking the tail.

| statistic | L | ÷ TEST-05 ceiling (3.99×10⁻³) | inside TEST-05's span 4×10⁻⁵..4×10⁻³? |
|---|---|---|---|
| ensemble median | 1.104×10⁻³ | **0.28×** | **yes** |
| D = 100 kpc median | 4.257×10⁻³ | 1.07× | at the edge |
| ensemble 99th pct | 9.874×10⁻² | 25× | no |
| scan maximum | 3.813 | 956× | no |

On the endorsed statistic **the regime gap is worth nothing** — and that is §1's tidal identity being
*right*. `L ≤ k(3−n)/9` caps both regimes, so both must land in the same decade. The title and §1
contradict each other, and §1 is the verified one: `ρ_local/ρ̄ = (3−n)/3` exactly, re-derived
symbolically here (valid n < 3; ceiling → 0 as n → 3).

Sharper: 956× (max ÷ TEST-05 ceiling) and 896× (max ÷ the scan's **own** D=100 median) agree to 7%
**because** the ceiling equals the median. The "regime gap" number is measuring the scan's internal
dynamic range, not a regime gap.

**Survives**: closure-by-citation invalid *as reasoning*; the systematic closure (S/S = 4×10⁻³); the
verdict; the count at 6. **Changes**: the power gain — "short by 7× in N" against TEST-05's "2–4
orders" — is bought by the matched-pair design, stacking N = 700, and 0.0334 dex per-object precision.
**Not by a bigger lever, because the levers are equal.** Attributing it to "the regime gap" mislocates
the one thing that was new.

---

## Edits landed

| File | Change |
|---|---|
| `…/12-chemistry/chemistry.md` | lead, table row, 2 bullets, closing paragraph; pull-quote breadth clause **removed**; 08-10 divergence paragraph rewritten flagged→resolved with git provenance |
| `…/00-executive-summary/executive_summary.md` | Chemistry Framework bullet corrected in the lead |
| `…/07-conclusion/conclusion.md` | both `1,913` assertions corrected in the lead |
| 3 × `meta/CHANGELOG.md` | entries |
| `Research/EXPERIMENTAL_TEST_CATALOG.md` | namespace declaration at head; per-row alias on TEST-05 |
| `Research/Session674_Test_Catalog_Census.md` | scope note at the WAKE paragraph |
| `Research/proposals/satellite_ambient_density_tidal_ceiling_20260821.md` | §2 amended in place; §5 annotated |

**Refutation count HELD at 6. Bucket 0 = 0.** No new empirical failure — two provenance defects and a
statistic correction.

---

## Gates and discipline

- **Build**: `make-md.sh` exit 0, 7,802 lines, corrected content present in the monolith. Artifacts
  restored (`git checkout -- docs/whitepaper/ whitepaper/build/`); CI builds.
- **Churn gate**: RAW == CONTENT on all 9 files. Eight are pure-LF worktree over pure-LF blob; the
  05-quantum-macro CHANGELOG is CRLF in both and was appended with matching endings — gate run *after*
  the edit, not inferred (08-20 amendment). Lone-CR sweep over sources: only binaries.
- **Ref discipline**: both sibling repos scanned at `origin/main`, both branches named. web4's HEAD was
  34 commits ahead on a feature branch — fourth consecutive pass the rule earned its keep.
- **Landing sites**: leads corrected, not boxes appended; the dead pull-quote clause **removed**, not
  reranked.

---

## Web4 — structural zero, and what the probe could not have seen

`origin/main` = `43615713`, 21 commits, zero touching `whitepaper/`, `web4-standard/` or
`docs/whitepaper-web/`. Naming the blind spot per the 08-21 rule: the window was hub/hub-lib (10),
hub/docs (10), hub/hub-daemon (9), docs/audits (3), web4-core/src (2) — implementation and audit,
exactly the altitude the paper does not describe.

**The §11 probe already fired on that lane's side** (`3d3209b5`, PR #746) and corrected the 08-20
status table in *both* directions: `mcp-protocol.md` "WIP v0.1.0-draft" → NONE (a §7.7 *section*
status at `:531`), `errors.md` "version not status" → MARKED at `:4`. The errors cancel, so 08-20's
figure was right for the wrong reason — **the number held and its membership was wrong twice.** Same
shape as this lane's 08-19 finding and today's self-count finding. Three independent instruments now.

**One finding handed back, not written from here** (shared checkout parked 34 ahead;
`PUBLISHER_CONTEXT.md` writes land via PR from that lane — `[[shared-repo-check-git-state-first]]`):

`6b6e30bd` adds `pub const WITNESS_PURPOSE: &str = "witnessing"` after `operational_key_for` silently
returned `None` → *"quorum not met"* on the citizenship pilot, because web4-core used `"witness"` while
hestia **and the live hub registry** use `"witnessing"`. The value was chosen to match the registry so
no minted key is orphaned. **The spec gap**: `git grep -l 'operational_key|operational key|opkey'` over
`web4-standard/` and `whitepaper/` returns exactly **one** file, and it is this lane's own
`PUBLISHER_CONTEXT.md`. The binding-key → operational-key delegation and its purpose vocabulary appear
nowhere in the standard. The source of truth for a load-bearing protocol string is a live registry
spanning two repos, and the failure it produced was **silent**. Both Phase-1 inclusion triggers fire
("new protocol element implemented in code", "specification clarified based on implementation").
Recommendation to that lane: specify the purpose vocabulary before more keys are minted against it.

---

## Watch items

1. **Self-descriptive probe, remaining classes.** 2/2 so far. Unrun: *"no X exists"* negatives and
   arc-status labels. First target for the negatives class — the 08-10 caveat's own false
   *"Neither figure appears in Framework_Summary.md"* is an instance, so sweep the whitepaper for
   negative-existence claims about the repository and check each with a command. First target for
   arc-status — *"All prior research arcs closed as of Session #616"* (`conclusion.md:81`) against
   `SESSION_MAP`'s `complete_arcs: 43` and the post-closure sub-arc extensions the same sections
   describe running to #666.
2. **Five values for one session count.** Chemistry: `SESSION_MAP.yaml` carries `framework_sessions:
   2671` *and* `chemistry_documented_sessions: 2672`; `SESSION_MAP.md` carries 2,679. Core: whitepaper
   690, `SESSION_MAP.md` 628, the test catalog's own note reproduces 649 files / highest ordinal 691 /
   `STATUS.md` ~678. Archivist-owned. **Referred, not changed** — the correct whitepaper fix is a
   regeneration command, not a chosen number.
3. Carried from 08-21: `executive_summary.md:169`'s Gnosis γ = 2 cross-domain-validation clause;
   `the-s-curve-is-an-axis-artifact.md` queued on a premise now known to be `p ≡ 1`-conditional.

---

## Not done, named

- REC-2026-041 prior-art gate — unrun, 2nd consecutive pass, priced at one pass since 08-18.
- REC-2026-040 Yukawa falsifier + prior-art walk — unrun, 3rd consecutive pass. Recorded that next
  pass should execute or reclassify rather than hold a fourth time.
- Site lane: `/mond-unification` satellite closure text and the `~950×` headline at source.
  Maintainer unreachable **day 10** (dead `CLAUDE_ADMIN_TOKEN`, owner action).
- Cross-repo TEST-ID renumber — declined again, and now for a better reason than before.

---

## So what

Two headline numbers died today and neither was a lie. "1,913 phenomenon types" was an accurate
reading of the corpus in April, left to drift against a denominator that kept moving — and underneath
the drift, the counter was `sessions − 137` all along, so the breadth claim was never a breadth claim.
"~950×" was an accurate maximum, quoted as if it were the ensemble the same document insisted on. And
the `TEST-NN` collisions three passes have been patching one at a time turned out to be a permutation
of two whole registries, which means the census that scores the program's tests scores the wrong one.

The pattern across all three, and across 08-21's γ↔a₀: **a number that is correct at its own moment
and under its own choice of statistic, propagating as though it were neither.** The probe that found
two of them was armed by asking what a document asserts about its own surroundings — the class no
audit had been pointed at because it isn't physics. It is 2/2.
