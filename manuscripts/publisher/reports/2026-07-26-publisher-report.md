# Publisher Daily Report - 2026-07-26

## Verdict: SIGNAL — the cron started (first autonomous run in 5 attempts), and it immediately found the claims gate dead again. Yesterday I fixed the instance; today I fixed the class, and it survived a live CI churn. Also closed the 47/48 accounting defect carried since 07-04. No new research.

This is an autonomous run, not a manual one — PID 461000, started 03:29:59 from the
Task Scheduler path, matching prompt. The 07-25 launcher fix (`c1981fd7`, drop `claude -c`)
works. Four prior no-starts (07-16, 07-21, 07-25 + one earlier) are ended.

**Correction to my own log-based evidence**: I have been citing empty `publisher/logs/*.log`
files as no-start evidence. Every one of the 23 July logs is *exactly 270 bytes* — start
header + complete header, zero captured output — because before `c1981fd7` the launcher never
redirected `claude`'s stdout at all. Those logs could not have distinguished a start from a
no-start on any day. The no-start conclusions were sound (they rested on missing reports and
missing commits), but the log-size evidence I cited alongside them was worthless. Today's log
is 135 bytes mid-run, and should carry real output at exit — verify next pass.

---

## 1. The claims gate was dead again, 24 hours after I called it "green, verified live"

Ran `render_claim_surfaces.py --check` at session start. **RED**:

```
ValueError: v1 freeze violation:
docs/whitepaper/Synchronism_Whitepaper.pdf: differs from frozen commit
```

Yesterday's self-correction said: *"when I praise infrastructure, I check whether it has
survived a real update — or I flag it as unexercised."* I applied that rule to the *praise*
and missed the *fix*. The 07-25 repair snapshotted the live `PREDICTIONS.md` and re-pointed the
manifest at it — and stopped there. It never asked the class-level question: **which other
frozen paths are living paths?** Four of the five were.

- `docs/whitepaper/Synchronism_Whitepaper.pdf` and `..._Complete.md` sit under `docs/`, which
  `build_whitepaper.yml` rewrites and commits (`git add docs/ whitepaper/build/`).
- Repo-root `README.md` and `SPINE.md` are live documents anyone may edit.
- Only `claims/v1-snapshot/PREDICTIONS.md` — yesterday's one-file fix — was actually immutable.

**The trigger makes the point sharper than the diagnosis does.** Commit `186c80c6` changed the
PDF from 754603 → 754617 bytes with `Synchronism_Whitepaper_Complete.md` **byte-identical**.
No content changed. A content-free rebuild broke a content freeze, because the build is not
byte-reproducible and the manifest was hashing a build artifact.

### The loop was mine, closed end to end

Tracing what triggered that CI run: `b9b31d30`, my own 07-25 commit, whose entire content was
appending "reviewed, nothing changed" to `whitepaper/PUBLISHER_CONTEXT.md`. That path matches
the workflow's `whitepaper/**` filter.

```
Publisher logs "nothing changed" to whitepaper/PUBLISHER_CONTEXT.md
  → matches the whitepaper/** push filter
  → full rebuild + Pages deploy
  → PDF differs (metadata only), source .md identical
  → CI commits ~754KB of binary churn into docs/
  → breaks the v1 freeze pointed at that PDF
```

**The Publisher's record of changing nothing was itself a change, and it broke the Publisher's
own integrity gate.** Same trigger on 07-18, 07-23 and 07-25: the last four pure-deploy commits
all rebuilt the PDF with unchanged source.

### Fixed at four levels (`378dd0cb`, `13e728a1`, `66afbd40`)

1. **All five artifacts snapshotted** into `claims/v1-snapshot/`, each verified byte-identical
   to its frozen sha256 — the evidence guarantee is unchanged, only the path is.
2. **Rule enforced, not remembered**: `verify_v1_freeze()` now rejects any manifest entry
   outside `claims/v1-snapshot/`, with the remedy in the error text. A live path in the
   manifest is now a design-time error instead of a gate failure discovered days later.
3. **Alarm decoupled from rendering** — the actual mechanism of the two-day outage. The freeze
   check ran *before* rendering, so any violation aborted the renderer: an evidence-integrity
   alarm could silently disable the publication surface. Now `--check` still raises (CI fails
   on tampering), and the write path warns loudly, renders anyway, and exits non-zero.
4. **Trigger fixed**: `PUBLISHER_CONTEXT.md` excluded from the build trigger. A run log is not
   whitepaper content.

### Verified, not assumed — twice

- **Before committing**, injected a bad hash and exercised both paths: `--check` → exit 1 with
  the raise; render → exit 1 *with the surfaces still written*. That is the behaviour I claimed,
  observed rather than reasoned about.
- **After committing**, a real CI deploy (`beebcc90`) churned the PDF again mid-pass — 754617 →
  754607, a *third* distinct size for identical content. Re-ran the gate: **`checked 10 claims;
  v1 freeze verified`, exit 0.** Under yesterday's manifest that commit would have taken the
  renderer down for the third time in three days. The fix has now survived a real update.

---

## 2. Closed the "frozen 47" defect — carried since 2026-07-04, deferred every pass since

The executive summary and conclusion each asserted **"47 genuine contributions"** (#615) and
**"Grand total: 48"** (#616) as the headline aggregate — and then, three bullets later in the
same section, reported S634's finding that 47 is a **~57% overcount** against the canonical
S582 inventory of **30**. The surfaces stated a number and refuted it in the same breath,
leaving the reader to pick.

Fixed by naming the classifier dependence **at the point of assertion** (`d6352137`):

- #615 bullet — 47 marked as #615's own classifier, *disputed, not to be quoted alone*, pointing
  at #634.
- #616 bullet — "Running total: **31 under the canonical strict criteria** (S582's 30 + 1), or 48
  under #615's more lenient classifier … never reconciled, so the strict figure is the one to
  quote."

Deliberately **not** a silent 48→31 swap: these are two different measurements under two
different criteria, and reconciling them is a research task nobody has done. Hiding that behind
a single clean number would repeat the original error in the opposite direction.

Also caught: **REC-2026-036's `recommended_framing` carried the bare "48 contributions" as
publication advice** — it would have propagated the disputed figure into a paper. Corrected.

I chose to fix rather than re-flag this because three weeks of "carried, flagged again" is
perseveration wearing the costume of conservatism. The durable fix remains open and is now a
recorded pending proposal: migrate the aggregate into `claims/claim_ledger.json` so the renderer
regenerates it under the `--check` gate — one machine-checked source instead of hand-maintained
prose in two places.

---

## 3. The information needed to predict this failure was in my own log entries

Ran `make-md.sh` to regenerate the built surfaces after the section edits. The diff was **3,177
lines** — the entire document. Diagnosed before committing: not content, **CRLF**. Building on
the WSL `/mnt/c` mount emits `\r\n`; CI emits `\n`. Same line numbers, same 68 H2 headers, same
7,243 lines, every line "changed." Reverted the generated artifacts, committed source sections
only, let CI rebuild.

**But I should not present this as a discovery.** Every no-change pass in `PUBLISHER_CONTEXT.md`
back to 2026-07-01 carries the same parenthetical: *"No rebuild performed (would yield only
CRLF/timestamp churn)."* I have written that line at least nine times. Today I ran the build
anyway, against my own standing note, and then spent the effort re-deriving why it was a bad
idea. Knowing a thing in a log is not the same as knowing it at the moment of acting.

The sharper version is worse for me. Those same log entries record the deployed PDF's byte size
pass after pass: **754599** (07-17), **754607** (07-17 deploy), **754612** (07-18), **754603**
(07-23), and today 754617 → 754607. I personally wrote down, four times in the week before
2026-07-23, that this artifact's bytes change with no content change — and then on 07-23 I froze
it by hash. The evidence that the freeze would fail was already in my own notes, in my own
handwriting, before I set it. This is not "infrastructure surprised me"; it is **not reading my
own record before acting on the thing it describes.**

So the day's pattern is one finding, not three: **artifacts under version control that are not
reproducible from their sources** — the PDF (five distinct sizes for one content), the monolithic
`.md` (host-dependent line endings), and by extension the web build. Freezing, diffing or gating
any of them by hash will keep producing false alarms. Recorded as a pending proposal with the
concrete one-liner (`docs/whitepaper/** text eol=lf`), left as **dp's call** because a
`.gitattributes` normalization rewrites those files once on the next commit.

---

## Phase 0: Publication Recommendations

### New Recommendations / Status Changes
- **None.** Readiness **HELD**: REC-037 0.98 / REC-034 0.97 / REC-035 0.95 / REC-036 0.68.
  `last_updated` → 07-26. REC-036 edited for the accounting correction (framing text only, not
  readiness). No new core session — S691 unchanged, no new numbered Session, no dp packaging
  go-signal.
- **Arc status**: substrate-physics **AT REST** (programmatic); research stack live; emergence
  phenomenology active (exploratory).

### Packaging-prep flags for dp (carried + updated)
- REC-037's entry **body remains stale** vs its July portfolio (TEST-08/09/10, compander,
  Cassini + the 07-25 literature-benchmarked re-verification + the EFE-axis input). Report series
  is the source of record until a packaging-time refresh.
- **Updated**: the claims-ledger renderer is now materially more robust than when I flagged it
  yesterday — snapshot-only freeze enforced in code, and the freeze alarm can no longer take the
  renderer down. Any packaging leaning on it should keep `--check` in CI. The *unexercised* label
  I owed it yesterday is now discharged: it has survived a live CI churn.
- **New**: do not quote 47 or 48 contributions in any packaged artifact. Strict figure is 31
  (S582's 30 + η-audit 1), and the classifier disagreement should be stated, not resolved by
  picking a number.

### Upcoming Candidates
- Emergence phenomenology arc — watch for synthesis. EFE axis — dp's preprint-scoping decision
  (routed 07-25, unchanged). Queued "unfalsifiable badge class" proposal — next signal-gate item.
  No runnable-and-unrun registered test.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: Current, **with one integration made today** (the first section change since
  07-16's `e9c51080`). **Sessions reviewed**: through S691.
- **Changes made**: contribution-count annotation in `00-executive-summary` and `07-conclusion`
  (2 lines each), CHANGELOG entries in both `meta/` dirs, REC-036 framing correction.
- **Proposals recorded** (3, in `whitepaper_sync.json`): aggregate → claim ledger; repo
  line-ending policy; move the Publisher run log out of `PUBLISHER_CONTEXT.md`.
- **Claims gate**: found dead, repaired at four levels, **green and verified live after a real CI
  churn**.
- **Terminology concerns**: None. Canonical terms unchanged (C(ξ), γ, MRH, Intent, Entity).

### Web4 Whitepaper
- **Status**: Current — and for the first time, **actually recorded as reviewed**.
- **Self-correction**: `whitepaper_web4.json` has said `last_review: null`,
  `status: pending_initial_review` since Phase 1 launched in January, while daily reports
  asserted "Web4 Whitepaper — Status: Current" and web4 accumulated **eight consecutive**
  "log … Publisher no-change verification pass" commits (07-16 … 07-25). The checks were real;
  they never wrote the state file. The machine-readable status contradicted the reported status
  for six months and I never reconciled them. Now recorded properly against `web4@5c2dd39`.
- **Reviewed**: `5c2dd39` (hub alarms the sender when a notice is dropped at cap/TTL; new
  `notice-dropped` kind with self-regress guards; 81+196 tests green), `6d118ec` (C270
  t3-v3-tensors audit, 2 MED + 1 INFO routed).
- **Verdict: no change, on scope grounds** — stated rather than asserted. Since the 07-09
  rewrite the whitepaper is an equation-ordered introduction to the *standard* (LCT 34 / V3 22 /
  T3 21 / ATP 18 / MRH 17 / ADP 17 / R6 8 mentions in 487 non-empty lines) and documents no
  transport layer: **"notice" appears zero times**. A new notice kind and its drop semantics are
  hub-implementation reliability one layer below the standard; documenting them would widen
  scope, not correct an error.
- **Watch item**: if mesh delivery semantics become normative — the standard promising delivery
  or dead-lettering rather than leaving it to the hub — that scope boundary moves and a transport
  section becomes owed.
- **Terminology concerns**: None; canonical LCT/T3/V3/ATP/ADP/R6/MRH all used in-definition.

---

## Infra / Governance
- **Publisher cron**: **STARTED** — first autonomous run after four no-starts. The launcher fix
  is validated for start; stdout capture is validated only at exit, so **verify tomorrow** that
  today's log is >135 bytes. Unexercised until then, and I am labelling it that way on purpose.
- **CI trigger exclusion** (`66afbd40`) is likewise **unexercised**: today's `PUBLISHER_CONTEXT.md`
  write is its first real test. If no deploy commit follows this pass, it works. Verify next pass
  rather than assume.
- **Standing Agency Grant (dp, 2026-07-16)**: unchanged — supervisor-scoped, not self-applied.
- **Scope note**: I made one shared-CI change (`build_whitepaper.yml` trigger filter) on my own
  authority. Justification: it removes a demonstrated harm loop caused specifically by Publisher
  writes, it is one line, and it is trivially reversible. Flagged here for dp to veto.
- **Working tree**: clean apart from this pass; `AGENTS.md`/`CLAUDE.md` supervisor edits from
  prior days are no longer showing as modified.

---

## Summary

The cron started for the first time in five attempts, and the first thing it found was that the
gate I repaired yesterday and reported as "green, verified live" was dead again — broken by a
CI rebuild that changed 14 bytes of PDF metadata with the source byte-identical. Tracing the
trigger closed a loop on me: my own commit appending "reviewed, nothing changed" to a file under
`whitepaper/` is what launched the rebuild that broke my own integrity gate. Yesterday I fixed
the instance; today I fixed the class — all five frozen artifacts moved to immutable snapshots,
snapshot-only freezing enforced in code so a live path is a design-time error, the freeze alarm
decoupled from rendering so it can never again take the publication surface down, and the
spurious trigger removed. Both behaviours were exercised against an injected fault before commit,
and the repair then survived a real CI churn mid-pass. The uncomfortable part is that this was
predictable from my own notes: I recorded that PDF's byte size changing four times in the week
before I froze it by hash, and my own standing line "no rebuild performed (would yield only CRLF
churn)" appears in nine prior entries — I ran the build anyway today. I also stopped deferring
the 47/48 defect and closed it: the whitepaper no
longer asserts a contribution count it refutes three bullets later, and REC-036 no longer carries
the disputed figure into publication advice. The Web4 review is recorded honestly for the first
time — its state file had said `pending_initial_review` for six months while my reports said
"Current." No new research today; readiness held across the board.

The day's real lesson is narrower than "check your infrastructure," and it has two halves.
**Fixing the instance I was shown, rather than asking what class it belonged to, is what cost the
second outage.** And the class was already documented — in my own log, by me, before I made the
mistake. The gap isn't detection; today's pass detected everything within minutes. The gap is
that I write careful notes and then act without reading them.
