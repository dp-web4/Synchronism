# Publisher Daily Report - 2026-07-27

## Verdict: SIGNAL — a publication-ready manuscript has been sitting in the repo for four days and my Phase-0 scan could not see it. Also: yesterday's CRLF diagnosis was wrong, and the real defect was corrupting the published whitepaper.

Autonomous cron run, second consecutive successful start. The 07-25 launcher fix
(`c1981fd7`) is now confirmed on repeat rather than on a single observation.

One correction to my own process at the top, because it is the most important thing in this
report: I flagged today's 135-byte log as a silent cron failure before checking the clock.
The file's mtime was one minute old — `date -u` 10:30 UTC is 03:30 local, so the log was the
*opening* banner of the run I was executing. The closing banner is written at exit. I was
looking at my own footprint and calling it a corpse. Checked before reporting; no failure.

---

## 1. Phase 0 found something, and the finding is partly about Phase 0 itself

`Research/papers/repository-mediated-autonomous-science/` contains a complete 4,307-word
manuscript — abstract, related work, methods, seven numbered results, a five-part model,
six design implications, a limitations section — plus a reproducibility manifest, a
stdlib-only analysis script, and a manually audited claim-lineage CSV.

It landed **2026-07-23** in PR #1 under the subject *"Add ledger-first research surfaces"*
and has not been touched since. The 07-24, 07-25 and 07-26 Publisher passes all ran and none
of them saw it.

**Why I missed it**: my Phase-0 scan reads `Research/Session*.md` and `SESSION_MAP.yaml`. A
publication-ready artifact that arrives as a *paper* rather than as a numbered session is
invisible to it. This is a scanner that looks for publication candidates in exactly one
shape and reports "no new candidates" when one arrives in another. Scan widened to
`Research/papers/` from this pass forward.

That failure is the same shape as the two I logged this week — the freeze manifest pointing
at live paths, the metric nobody decomposed. Not a detection failure. A *coverage* failure
that reported itself as a clean result.

### What the paper says

It studies `gnosis-research` — a repo whose design wakes a fresh, memoryless LLM instance
every few hours and hands it the accumulated repository — as a longitudinal case in claim
formation, correction, and drift. The thesis is sharper than the setup suggests:

> Repository-mediated succession delivers fast local self-correction without weight updates
> or persistent memory. But append-only narrative memory preserves framing, lets evidence be
> bypassed, and does not propagate corrections. The decisive property is **admission and
> maintenance of claims**, not memory.

The load-bearing contrast: five locally adjudicated claim lineages reverse within
**5.7–95.2 hours** (median 59.3), while two wrong-variable claims take **~92 days** to
correct in a parallel repo and remain uncorrected in the source repo's public framing. In
one case, adverse evidence is committed **~12 hours before** a successor declares the
contradicted pattern a universal law.

### I verified it rather than describing it

`gnosis-research` is checked out locally, so the reproducibility contract is testable. I ran
it read-only against the frozen snapshot head `7c1c16d0`:

| Manifest claim | Reproduced |
|---|---|
| 138 commits at snapshot head | ✅ 138 |
| Subject classes 45 / 14 / 2 / 77 | ✅ 45 / 14 / 2 / 77 |
| First commit hash + timestamp | ✅ exact |
| Last commit hash + timestamp | ✅ exact |

One command, no manual steps. **No other recommendation in `recommendations.json` has had
its central quantitative claim independently reproduced.** Thirty-seven recommendations rest
on session narratives; this one rests on a hash I can check.

I also tried to break it. The month histogram shows a hole — Feb 21, Mar 45, Apr 66, **May
0**, Jun 5, Jul 1 — so "longitudinal 2026-02 through 2026-07" is really three dense months
plus a thin tail, and 132 of 138 commits fall in Feb–Apr. That would be a framing problem if
the paper hid it. It does not: §4.1 states plainly that May contains no commits and
explicitly declines to infer why. Checked, not assumed; it holds.

**Filed as REC-2026-038**, readiness 0.93, priority HIGH, target arXiv cs.AI.

### The honest caveats, stated where dp will see them

- **Self-study conflict is structural.** The system studied and the system writing are the
  same lineage of agents on the same steward's infrastructure. The paper discloses this; it
  cannot neutralize it. The one real mitigation is that the automated half is now
  third-party reproducible from a frozen hash, so a skeptic can check the quantitative spine
  without trusting the authors.
- **The manual lineage audit — the load-bearing half — is not reproducible.** Nine rows of
  hand-coded CSV by a non-blinded coder who had already read the corpus.
- **Readiness scores in this file are not commensurable.** REC-034/035/037 sit at 0.95–0.98
  scoring *arc maturity with the manuscript still unwritten*. REC-038's 0.93 scores a
  *written, verified manuscript with no external review*. On time-to-postable, 038 is the
  readiest thing in the file by a wide margin. I am not renumbering the others; I am naming
  the incommensurability rather than pretending to one scale.

### Advisory: this reorders the pending preprint strategy

`Research/proposals/stable_fixed_point_preprint_strategy.md` (awaiting dp since 07-23)
proposes three preprints. Candidate #3 — the A2ACW program-level null, cs.AI, flagged
"most novel-to-audience" — occupies the same slot as this manuscript and is materially
further behind: #3 is an outline, this is drafted, reproducible, and now verified.

**Recommendation**: if dp approves the strategy, make this the cs.AI vehicle and fold the
A2ACW self-play null into it as a results subsection or companion, rather than drafting #3
from scratch. Revised advisory order: **REC-038 (ready now) → Locality No-Go (#1) → DESI (#2)**.

---

## 2. Phase 1 Synchronism: the build script was corrupting its own sources

Yesterday I ran `make-md.sh`, got a 3,177-line whole-file diff, called it CRLF churn from the
WSL `/mnt/c` mount, reverted, and wrote that diagnosis into two section CHANGELOGs and a
pending proposal.

**The diagnosis was wrong, and it hid a real defect.** The mount was a bystander.

`preprocess-sections.sh` converts `### Header` to `**Header**` by wrapping `substr($0, 5)`.
On a CRLF working tree that substring ends with the record's trailing CR — so the CR lands
*between the header text and the closing `**`* — and the script writes the result back into
the **source** file:

```
**Electronic Channel (Optical/Dielectric)<CR>**
**Post-Kimi consolidated open questions (2026-05-28)<CR>**
**Key Claims<CR>**
```

11 occurrences across 3 sections, committed, and reaching every published surface —
`docs/whitepaper/Synchronism_Whitepaper_Complete.md` carried all 11.

The second-order effect is what produced yesterday's phantom diff, and it is the part I had
no model for: **git classifies a file containing a lone CR as non-text, which silently
disables `core.autocrlf=input` for that whole file.** Its CRLF endings then diff against the
LF blobs CI produces — every line. I reproduced this in an isolated repo before asserting it:

| Fixture | Result |
|---|---|
| Clean text file, CRLF endings | normalizes silently, **no diff** |
| Same file + one lone CR | **diffs on every line** |

So the whitepaper build had a self-inflicted wound that disguised itself as an environment
problem, and my 07-26 self was one layer of the disguise.

**Fixed** (`8a0fd27f`): `sub(/\r$/, "", header_text)` in the script, verified idempotent
against a synthetic CRLF fixture; 11 stray CRs stripped from 4 files. Each repaired file is
byte-identical to its predecessor once CR characters are stripped from both — **zero text
changed**, verified per-file, not assumed.

Two notes on what I deliberately did *not* do:

- **`claims/v1-snapshot/` keeps the uncorrected bytes.** It is frozen evidence of what v1
  published, not a copy to be kept current. Repairing it would have been the intuitive move
  and the wrong one.
- **The 07-26 CHANGELOG entries keep the wrong cause.** The new entries in 06-implications
  and 09-appendix-mathematical supersede them. These changelogs are append-only; editing
  history to look correct is not a correction.

Generated artifacts not rebuilt locally — CI is the authoritative builder, and this commit
legitimately triggers it.

### Infrastructure claims from yesterday, checked rather than re-asserted

Yesterday I flagged two fixes as unexercised. Both have now been exercised:

| Claim | Status |
|---|---|
| Claims freeze gate (class fix, `378dd0cb`/`13e728a1`) | ✅ green at start and after commit (10 claims) — survived a full day and a live CI deploy |
| CI trigger exclusion (`66afbd40`) | ✅ **passed its first real test** — `f1ce18af` touched `whitepaper/PUBLISHER_CONTEXT.md`, no Pages deploy followed |
| Launcher stdout capture (`c1981fd7`) | ✅ 3,937-byte log at exit yesterday vs 270-byte stubs before |

---

## 3. Phase 1 Web4: no drift in web4 — the stale artifact was mine

Reviewed against `web4@5df662a`. Two commits in window, both in whitepaper scope, both
already handled on the web4 side. No Synchronism-side action needed and none taken.

**`ad6e35c`** — web4's `build_whitepaper.yml` had been failing on *every* push-triggered run
since 2026-05-16 (46 runs, 45 failures; the single success a manual dispatch), rejected by
branch protection at the deploy push, while its CHANGELOG claimed a CI-built PDF that had
never existed. Corrected by a superseding entry.

Worth recording across repos: that find came from applying the class-level question — *"which
other instances share this defect?"* — that I recorded failing to ask on 07-25. **The same
question found a real defect in a second repository within a day.** That is the strongest
evidence so far that the habit transfers rather than just reading well in a report.

**`5df662a`** — a worked example (hestia#49) for why the relying party must compute trust
rather than accept an originating party's declaration. Standard-altitude, consistent with
section 11.

### A terminology flag I raised against myself

`5df662a` uses "R6/R7". My canonical terminology table listed **R6 only**, so the drift check
fired. I checked before reporting it: R7 appears in 19 whitepaper files, ships as an
`R7Action` type in web4-core v0.2.0, and carries a glossary definition — *R6 with Reputation
as a first-class output; the Request↔Result delta feeds back into the actor's T3/V3 tensors*.

R7 is legitimate and long-established. **The stale artifact was my own governance table**,
which means my drift detector was primed to emit false positives against correct usage.
Added R7 to `publisher/CLAUDE.md` this pass. Genuine terminology drift found: **zero**.

---

## 4. Research review

| Check | State |
|---|---|
| New numbered sessions | None — S691 unchanged |
| New complete arcs | None |
| Arc status | AT REST since 2026-06-24 (operator's call) |
| dp decision on three-preprint strategy | Still pending (filed 07-23, 4 days) |
| Last `Research/` commit | `4fd2a1cf`, **2026-07-24** |

The research lane has been silent three consecutive days after a dense 07-22→07-24 run
(EFE evidence-axis step-0, Milgrom prior-art audit, SPARC×Cassini). **Noted, not escalated** —
the arc is at rest by the operator's call, and silence is the expected state, not an anomaly.
Recording it because "expected silence" and "stalled track" look identical from here, and the
distinction is worth being able to make later if it continues.

Readiness held, no narrative churn: REC-037 0.98 / REC-034 0.97 / REC-035 0.95 / REC-036 0.60.

---

## 5. Housekeeping

- `recommendations.json`: 37 → **38** recommendations. Diff is now a **clean 34-line append**.
  Previous passes rewrote ~500–700 lines of this file every time by re-serializing with
  `ensure_ascii=False` against a file that uses escaped unicode. Future passes: dump with
  `ensure_ascii=True` and re-append the trailing newline.

---

## Summary

Phase 0 filed **REC-2026-038** — a drafted, reproducibility-verified cs.AI manuscript that had
been sitting in the repo for four days, invisible to a scanner that only knows how to look for
numbered sessions. Phase 1 traced yesterday's "CRLF churn" to a build script that was trapping
carriage returns inside the bold markers it generated and writing them back into the whitepaper
sources — 11 of them reached the published document — and fixed the script rather than the
symptom. Phase 1 Web4 found zero drift in web4 and one stale row in my own terminology table.
Research lane unchanged and at rest; three days quiet, noted without escalation.

**The pattern across all three phases today: every defect was in something that was reporting
success.** The scanner said "no new candidates." The changelog said "CRLF churn, environment's
fault." The drift table said "R6." Detection was never the weak link.
