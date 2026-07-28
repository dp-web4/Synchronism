# Publisher Daily Report - 2026-07-28

## Verdict: SIGNAL — yesterday the research lane retracted the universal form of the project's #1 transferable null. The retraction reached the ledger and not the two documents that carry the claim to a reader. Both fixed. And I think the retraction may itself be an over-correction, for a reason nobody has checked.

Autonomous cron run, **third consecutive successful start** — the 07-25 launcher fix is now
confirmed three times, not once.

---

## 1. Phase 0: the research lane resumed, and what it produced lands squarely on the publication queue

Three days of silence ended on 07-27 with three commits from a site-visitor pass. The substance:

A **pre-registered screening-literature gate** — verdict rules, a primary-source requirement, and a
falsifier for its own null, all fixed *before* the walk — closed the last residual hole the 07-23
prior-art audit had named. It found a **counterexample against the framework's own no-go**:

> **Burrage, Copeland & Millington (PRD 95, 064050, 2017)** reproduce the RAR on SPARC-153 with a
> screened scalar keyed on **local volumetric ρ(r)**, under **universal** Lagrangian parameters.

Two things follow. The first is that the no-go's universal form — *"any local-ρ MOND mimic fails"* —
is false as stated. The second is sharper: the same walk found that a supporting claim the project
inscribed on 2026-07-10 ("symmetron/chameleon screening *fails* to reproduce MOND on the ρ-vs-g_bar
mismatch") was **unsourced and inverted** — a fabricated consensus, cited in favour of a result the
project has a stake in. The explorer corrected both at source in `PREDICTIONS.md` the same day and
adopted two general rules from it: prior-art nulls must be stated corpus-scoped, never "none exists";
and any attribution-of-consensus must name a primary or be deleted, because the existing
"every p-value walks to a primary" guard has no walkable target for a consensus claim.

That is the research lane catching itself over-*defending*. It is the right catch and I am not
second-guessing it.

### What I found: the correction stopped at the ledger

`PREDICTIONS.md` was fixed. Two surfaces that carry the same claim to a reader were not.

| Surface | State on 2026-07-28 morning |
|---|---|
| `PREDICTIONS.md` (ledger) | ✅ corrected 07-27, two rows |
| `whitepaper/sections/00-executive-summary` + `07-conclusion` | ❌ still asserting the universal form |
| `Research/proposals/stable_fixed_point_preprint_strategy.md` | ❌ still asserting it — **and this is the document awaiting dp's decision** |

The whitepaper's canonical statement, in both sections:

> *a function of local ρ(r) cannot reproduce the RAR — the modification must key on a non-local
> functional of the baryon distribution.*

The paper names exactly **two** surviving decisive negatives. That sentence is one of them. So for a
day, two published surfaces of the same project disagreed about whether one of its two load-bearing
negatives was true.

The preprint strategy document is worse, because it is the one with a decision pending on it. Its §1
states null #1 in the universal form **with an escape clause** — *"except by per-system calibration"* —
that Burrage's **universal** Lagrangian parameters defeat directly. dp opening that file yesterday
would have read a falsified claim with no marker on it.

Both swept this pass. Nothing deleted; `[SCOPED 2026-07-27]` markers in place, per the S668/S672 and
TEST-04a corrective-propagation precedent, with CHANGELOG entries in both `meta/` dirs and a full
Gate Update appended to the strategy proposal.

---

## 2. The part I am least sure of, stated as a question: the retraction may be an over-correction

The triage treats Burrage as a counterexample. Under the classification table's **own** criterion, it
may not be one.

A screened scalar's static field equation is

```
∇²φ = V_eff,φ(φ; ρ(r))
```

— elliptic, sourced by ρ. Its solution φ(r) therefore depends on ρ **throughout** the region, so the
force `∝ ∇φ` is a **non-local functional of the baryon distribution**, even though the *potential* is
keyed on local ρ. That is the same structural reason AQUAL is classed non-local in the very table
under revision.

So "local" has two readings, and this archive has never chosen between them:

| Reading | "Local" means | C(ρ) | Burrage |
|---|---|---|---|
| **A — algebraic-pointwise** | function of ρ at the field point, `g = f(ρ(r))·g_N` | local | **local** → counterexample stands, no-go narrows |
| **B — field-sourced** | functional of the ρ distribution, however obtained | local | **non-local** → **not a counterexample; the no-go survives intact** |

Under Reading B the honest fix is a one-sentence definition, not a scope reduction — and the 07-27
retraction over-shot.

**I am raising this as a question, not a finding.** I have not read Burrage; the field equation above
is the standard screened-scalar form, not a reading of that paper. Filed as
`Research/proposals/locality_operational_definition_algebraic_vs_field_sourced.md` **with a
pre-registered falsifier**: if the answer is Reading A, this proposal was wrong to raise the
possibility, and that gets recorded rather than quietly dropped.

**What does not move either way**: Synchronism's C(ρ) is algebraic-pointwise under *both* readings, so
the framework-specific kills (TEST-09 BTFR bounded-boost; ρ_crit ∝ V⁺² sign inversion) stand
untouched. Bucket 0 stays 0. This is a question about the *transferable* claim only.

**Why I bothered.** Yesterday's triage caught the project over-defending a result it has a stake in.
Accepting the resulting retraction uncritically — because retractions *feel* epistemically safe — would
be the same deference failure with the sign flipped. Both directions means both directions.

---

## 3. What this does to the publication queue

`stable_fixed_point_preprint_strategy.md` proposes three preprints and has been waiting on dp since
06-28, gate-updated 07-23 and again today. Null #1's headline sentence, its novelty relative to
Milgrom, and its required citations all move depending on which reading is correct:

- **Reading B** → preprint #1 keeps its general claim, gains a sharpened definition, and cites Burrage
  as *scope-setting*. Strongest version, worth arXiv.
- **Reading A** → it narrows to ρ-threshold mimics tying ρ_crit to a₀, with Burrage a required citation
  and a named occupied escape class. Publishable, materially less novel-to-audience.
- **Either way**, Burrage is now a required citation. Submitting §1 without it invites precisely the
  referee objection the gate just raised internally.

**Revised advisory order** (supersedes yesterday's): **REC-038 (drafted, verified, cs.AI) → #2 DESI
mechanism-class → #1 Locality No-Go (blocked on the fork) → #3 A2ACW.** Null #1 drops below DESI only
because it now has an open blocking question — not because it weakened.

**One framing option worth dp's attention**, surfaced by the corrections rather than by me: the triage
records that the symmetron is *the same mathematical object* as C(ρ), with the three ingredients C(ρ)
lacks — an action, a derived ρ_*, a Z₂ order parameter, mean-field β=½. *"A framework independently
reconstructed the symmetron minus its action, and the missing ingredients are exactly what killed it"*
is a more legible story to an external audience than either version of the no-go, and it uses material
the project already owns.

### Status changes

| Rec | Change | Why |
|---|---|---|
| **REC-2026-037** | readiness **0.98 → 0.92**, weakness + framing addendum, `date_updated` 07-28 | Its named TRANSFERABLE ARTIFACT *is* the locality classification table, now under active correction. "Stable, not actively evolving" no longer holds. |
| REC-2026-038 | held (0.93, HIGH) | Unchanged; still the readiest thing in the file on time-to-postable |
| REC-034 / 035 / 036 | held (0.97 / 0.95 / 0.60) | No new evidence |

`recommendations.json` diff: **6 lines**. (Second consecutive pass without the ~500-line
re-serialization churn — `ensure_ascii=True` plus the trailing newline.)

### Widened Phase-0 scan (the 07-27 fix, exercised)

`Research/papers/` — one item, `repository-mediated-autonomous-science`, unchanged since 07-23, already
REC-2026-038. `Research/preregistrations/` — one item, `sparc_cassini_tanhlog`, unchanged since 07-23,
already handled 07-24. `Research/proposals/` — three new files, all from the screening gate, all
handled above. **No new numbered session (S691 unchanged), no new complete arc.** The widened scan
found nothing new today, which is the correct result and the first time I can say that with the
scanner actually looking in the right places.

---

## 4. Phase 1 Synchronism: gates checked before touching anything

| Check | Result |
|---|---|
| Claims freeze gate (`--check`) | ✅ green, 10 claims, before edits |
| **Lone-CR regression** (yesterday's fix) | ✅ **clean** — `git grep -lIP '\r(?!\n)'` over `whitepaper/**` and `docs/**` returns zero text files |
| Source ↔ artifact sync | ✅ artifacts 07-27 03:47 postdate last section commit 03:38 — CI rebuilt from the CR fix |
| Markdown/pandoc smoke test on both edited sections | ✅ both parse; `**` parity anomaly in exec-summary confirmed **pre-existing on HEAD**, not introduced |
| Forbidden patterns (§7) | unchanged |
| Terminology drift | zero |

The CR fix held through a live CI rebuild. That is the second infrastructure claim discharged from
"fixed" to "exercised" this week.

Generated artifacts deliberately not rebuilt locally — `make-md.sh` runs `git pull` as its first act,
which is wrong to invoke against a dirty tree mid-pass. CI is the authoritative builder and today's
section edits legitimately trigger it.

---

## 5. Phase 1 Web4: no drift, and one correction there that rhymes with everything else this week

Reviewed `web4@5df662a..206dd00` — eight commits, two in whitepaper scope, both already handled on the
web4 side. No Synchronism-side action needed and none taken. Zero terminology drift.

**`b2e2888`** is worth recording. §11 credited hestia with *"fail-closed defaults for unattended
operation."* Hestia's own bypass catalogue says the opposite in as many words — *"Default posture is
fail-open. Absent `HESTIA_PRE_FAIL_CLOSED=1`, an unreachable daemon means an ungoverned agent"* —
confirmed in code (`pre_tool_use.py:377`) and recorded **unset** on a live host. What was true is kept
and narrowed: hestia *is* fail-closed on invalid law. The self-contradiction removed is the sharp part
— the paper argues a relying party must **compute** trust from evidence rather than accept an
originating party's declaration, while itself accepting a *declared* safety property about its own
reference implementation.

**Cross-repo pattern, noted not actioned**: `206dd00` ("the Rust workspace has never run a test in CI —
draw the boundary") and web4's still-dead `build_whitepaper.yml` are the same shape as this week's
Synchronism findings. Three repos, one class: **a surface reporting success for work that never ran.**

---

## 6. Housekeeping

- **⚠ `private-context` is mid-rebase and my collective-log entry is not on `main`.** The repo is in a
  **paused interactive rebase** started by the supervisor's 03:00 session (`.git/rebase-merge` stamped
  03:29:24, ~20 min before I looked): detached HEAD, one pick remaining (`6fd171c85`), `main` two
  commits ahead of HEAD. Short-form `git status -sb` shows only `## HEAD (no branch)` and says nothing
  about a rebase, so I committed the log entry into it without noticing.
  **Deliberately not resolved.** Continuing or aborting another track's in-flight rebase is a
  consequential action on state I don't own. Instead I created a plain branch ref
  `publisher-log-2026-07-28` → `700a850e3` (zero working-tree impact, fully reversible) so the entry
  survives either outcome: it lands on `main` if the supervisor runs `--continue`, and is recoverable
  from the ref if it runs `--abort`. **Supervisor: the collective log entry for today is on that ref,
  not on `main`, until the rebase completes.**
- The stash-pop conflict in `machines/fleet/cbp.json` that started all this was **caused by my own
  documented pull ritual** (`stash && pull --rebase && stash pop`) on a pull I only needed for
  *reading*. Resolved keeping the newer upstream heartbeat (2026-06-06 over the stashed 2026-05-06),
  stash dropped, other agents' `.gz` working-tree changes verified intact with `gzip -t` and left
  alone. `--autostash` does the same job without the conflict. Memory corrected.
- **Nearly repeated yesterday's timezone error, in reverse.** I read the rebase directory's 03:29
  timestamp and started to call it a seven-hour-abandoned operation. Local is ~03:45 PDT and UTC
  ~10:45 — the rebase is *twenty minutes old and possibly live*, which is exactly why leaving it
  alone was the right call rather than merely the cautious one.
- `AGENTS.md` / `CLAUDE.md` carry uncommitted GitNexus index-count updates — supervisor scope, left
  untouched per precedent `a13894da`.

---

## Summary

The research lane retracted the universal form of the project's #1 transferable null yesterday, on a
pre-registered gate that also caught a fabricated consensus claim of its own. That correction reached
the ledger and stopped there — the whitepaper and, more consequentially, the preprint strategy document
sitting in front of dp both still asserted the falsified claim this morning. Both are now annotated.

The finding I am least certain of is the one I most want read: the retraction may have over-shot,
because a screened scalar's force is a non-local functional of ρ even when its potential is keyed on
local ρ — which under the classification table's own criterion would make Burrage no counterexample at
all. Filed as an open question with a falsifier, not asserted.

**The pattern this week has been "every defect was in something reporting success." Today's is
narrower and worse: the correction was real, correct, and well-made — and it still didn't arrive.
Detection wasn't the weak link, and this time neither was diagnosis. Propagation was.**
