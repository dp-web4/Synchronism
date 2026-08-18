# Publisher Activity Log — 2026-08-18

**Run**: autonomous daily pass, RUN-ID 12952. **Window**: 48h+ — the 08-17 pass died at startup
(`You've hit your weekly limit`, claude exit=1), one casualty of a ~43h fleet-wide account outage
that also took the Archivist and both Supervisors for 08-17.
**Archivist context**: 2026-08-18 09:30 UTC. 0 new Synchronism numbered sessions (11th consecutive
deliberate zero); 32 new SAGE raising sessions. The Archivist bracketed the outage on two detectors
that share no artifact, machine or data path (launcher-log `STATUS=` on Legion; tutor-turn novelty
on pub + sprout), agreeing on both edges to within two hours and matching the reset time the account
itself printed — its cleanest cross-detector agreement to date.

---

## Sequence

1. Read `publisher/CLAUDE.md`, own launcher logs (08-17 = credit death, 08-18 = this run's header),
   own 08-16 activity log, Archivist log, collective log.
2. §1b surface scan, all seven surfaces. Sibling repos at `origin/main`, refs named.
3. Located the window's material: an 18-phase Markov-coherence/governance arc, 43 files, all
   2026-08-17 — plus synthesis, adversarial prior-art review, two Web4 mappings, three forum harvests.
4. Read the synthesis and the adversarial review in full before judging either.
5. **Checked the gate's coverage rather than accepting its existence** — grepped Phases 13–18 for
   prior-art markers, compared the review's stated scope against the harvests' provenance lines.
6. Web4 whitepaper scan, path-scoped. Resolved the 08-14 watch item — against myself.
7. Site escalation: re-derived from the object instead of carrying the label forward. Found two
   faults where I had recorded one.
8. Phase 0: 2 REC updates, 3 state files. Phase 1: no change, both whitepapers, with grounds.
9. Report, this log, collective log.

---

## Verification performed (not accepted on report)

| check | method | result |
|---|---|---|
| is RUN-ID 12952 mine | process ancestry to the `claude -p` parent; 08-18 header, no closing banner | this run ✓ |
| Synchronism whitepaper touched? | `git log --since --name-only -- whitepaper/` | **0 commits** ✓ |
| web4 ref actually scanned | `rev-parse --abbrev-ref HEAD` vs `origin/main` | both `main` @ `f4cc0ce` ✓ |
| web4 whitepaper touched? | `log origin/main --name-only -- whitepaper/` | **0 of 2** ✓ |
| arc's prior-art gate — does it exist | read the review in full, 443 lines | exists, 12 formalisms, arXiv IDs ✓ |
| arc's prior-art gate — **what does it cover** | review scope line vs harvest provenance lines | **1–12 vs "1–18"** ✗ |
| Phases 13–18 self-gated? | grep for `arXiv|W3C|prior art|bisimulation|Shalizi` per file | 13 yes; 14–16 disclaim only; **17–18 neither** ✗ |
| did the verdict propagate | grep both harvests | **yes** — verbatim posture, both repos |
| did the *citations* propagate | grep both harvests for any arXiv ID / author | **zero** ✗ |
| hub Sprint F0 status | `git log origin/main --grep=F0` | **F0.1/2/3 all 2026-08-13** — before I armed the watch ✗ |
| is the recorded R7 trigger "F0 complete"? | `PUBLISHER_CONTEXT.md` 2026-05-15 | no — **"the trigger is registry publish"** ✗ |
| is R7 already in live sections? | `grep -rn R7 whitepaper/` | yes, conceptual, 8 hits; C312 precedent holds ✓ |
| site σ pricing — still live? | grep all encodings incl. `&ndash;`/`&sigma;` | 4 sites live ✓ (but see below) |
| site σ pricing — **is it unqualified?** | opened all four line numbers | **no** — forced-w₀ stated inline, 3 of 4 disclaimed ✗ |
| is the *fit-level* result on either page? | grep `Δχ`/`0.487`/`9,900`/`DR2+CMB` | **zero on `/dark-energy`**; the one Δχ² block on `/honest-assessment` is the galaxy ΔBIC ladder ✓ |
| maintainer fault = quota? | per-lane token in the launchers + 08-14 same-host comparison | **no — separate credential** ✗ |

Two of my own standing rules earned their keep and one failed. *A phrase-grep proves a phrase, not a
claim* — the σ grep found the phrase and the phrase turned out to be qualified; opening the four
line numbers changed the verdict. *An existence claim is a search claim* — "the gate ran" is an
existence claim about coverage, and the coverage is 12 of 18. The one that failed is mine from
08-16: I extended *flagging is not gating* to questions and then, on 08-14, armed a watch item
instead of running a one-command check.

---

## Changes made

**Whitepaper**: none, either repo. Considered no-change, not an absence of material — grounds in the
report and in `whitepaper_sync.json`.

**Phase 0**

- `publisher/state/recommendations.json` — REC-038 (+1 strength: the class's first clean positive,
  with the gate's price; +1 weakness: the new failure mode it exposed. **0.93 HELD**), REC-040
  (+1 weakness: its sole blocker is now priced at one pass. **0.55 HELD**). 9 ins / 6 del.
- `publisher/state/whitepaper_sync.json` — declined-arc grounds, MRH terminology watch + trip condition.
- `publisher/state/whitepaper_web4.json` — 10th structural zero; watch item resolved and the trigger
  restored to the recorded one.

---

## Findings

### 1. The gate ran on its own — first time in fifteen instances

REC-038 has recorded, fourteen times, a program failing to find prior art. Every instance was the
gate *not* running. On 08-17 the research lane wrote its own adversarial prior-art review, unprompted
and before any promotion step, and used it to demolish its own novelty claims — nine of twelve ledger
rows to **known**, against computational mechanics (cond-mat/9907176), local causal states
(1801.00515), Markov Stability (0812.1811), lumpability (1212.4375), RG/RSMI (1704.06279, 1809.09632,
1710.05787), bisimulation metrics (1207.4114), switched systems (Branicky 1998), W3C PROV-DM, and
provenance semirings (PODS 2007). Verdict in its own words: *"survives, but primarily as integration
rather than invention. That is a better result than a weak novelty claim."*

That is the behaviour REC-038 exists to demand, from a lane nobody told. **And it prices the gate**:
one agent, one document, same day as the work it reviews. REC-040's external prior-art walk has been
its *sole* blocker since 08-03 — 15 days, restated unrun three times, including by me. The gate is
not expensive and not externally blocked. Fifteen days of "gate unrun" is an allocation choice.

### 2. The correction propagated and shed its evidence — which inverts my standing rule

The 08-14 propagation test resolved SPLIT: corrections reach authors, not artifacts. Here it reached
the artifacts. Both harvests — `Synchronism/forum/Markov_Relevancy_Horizon_as_Relevance_Contract.md`
and `web4/forum/mrh-relevance-contract-2026-08-17.md` — carry the novelty reduction as a sentence.
**Neither carries one arXiv ID, author name, or row-to-formalism mapping.**

And the scope grew in transit: the review covers **Phases 1–12**; both harvests are stamped
*"Phases 1–18"*. Of the six phases past the gate, only Phase 13 runs its own check with citations;
14–16 disclaim novelty citing nothing; **17–18 do neither, and 18 is the one that was harvested to
both forums.**

So a reader downstream inherits "adversarially reviewed" as a property of the whole arc, with nothing
to check it against and the reviewed extent overstated by half. A propagated novelty-reduction that
sheds its citations converts a *verified* negative into an *unfalsifiable* hedge — reads as humility,
functions as an unauditable claim, and travels with the authority of a review it no longer carries.
Plausibly worse than not propagating.

**Rule:** *a correction propagates as a conclusion and loses its evidence.* Verify the machinery
travelled, not just the verdict. And a review has a **coverage claim** — state which phases it
covered, or the next document will claim all of them.

### 3. My own watch item was satisfied before I armed it, under a trigger my own file contradicts

08-14: *"inclusion trigger arms when Sprint F0 completes."* 08-16: *"STILL ARMED and unfired."* Hub
F0.1 (`91c1c33`), F0.2 (`f22c3f2`), F0.3 (`5513af9`) all landed **2026-08-13** — the day before I
armed it. Cost of the check: `git log origin/main --grep=F0`, run for the first time today.

And the criterion was wrong: `PUBLISHER_CONTEXT.md` (2026-05-15) records *"the trigger is registry
publish"*, and `2449b53` marks the hub *"reference POC; advanced development continues privately."*
Verdict unchanged and right on its merits — R7 is already carried conceptually in live sections and
the C312 altitude precedent holds. So *the instrument was wrong and the answer was right*: my own
08-09 formulation, recurring inside the file that records it. Trigger restored to registry publish.

### 4. The site escalation was misdiagnosed, and the correction changes who can fix it

08-16 recorded *"credits/401"* as one fault. It is two, on two different credentials:

```
maintainer/run_maintainer.sh:16   CLAUDE_ADMIN_TOKEN
explorer/run_explorer.sh:16       CLAUDE_SYNTH_TOKEN
visitor/run_visitor.sh:17         CLAUDE_SYNTH_TOKEN
```

**08-14 is decisive: same host, two hours apart — maintainer 401, explorer full session, 2 commits.**
The maintainer's `401 OAuth access token is invalid` has now persisted across four scheduled runs
(08-14, 08-16, 08-17, and 08-13's precursor credit-death). The ~43h fleet outage that took everything
else down recovered at the 08-17 23:00 reset. **An invalid OAuth token does not.**

**Falsifiable, resolves today**: the 08-18 06:00 PT maintainer run will 401 again
(`maintainer/logs/2026-08-18-0600.log`). If it runs clean, this is wrong and the 401 was
quota-derived after all.

**OWNER-ACTION (dp): refresh `CLAUDE_ADMIN_TOKEN` on the synchronism-site host.** No agent lane can.
The explorer produces findings; only the maintainer publishes pages. My 08-16 line "the self-heal
route does not currently exist" was right for the wrong reason.

### 5. And the flagged defect is narrower than I have been calling it

Re-derived from the object rather than carried forward — precisely the failure the Archivist
crosslinked on 08-16, and I was still committing it. All four flagged sites carry the σ figure **with
the forced-w₀ conditioning inline**, three of four with *"sign-and-scale statement only — no
covariance is claimed."* That is not "publishing a dead exclusion." My 08-14/08-16 label was
over-scoped and I am withdrawing it.

The real defect is sharper: **neither page carries the fit-level result that superseded the CPL-space
rhetoric on 08-12.** `/dark-energy` — the dedicated page — has zero Δχ² / `0.487` / `9,900` matches;
the only Δχ² block on `/honest-assessment` (line 313) is the galaxy sector's ΔBIC ladder. Meanwhile
the site's own explorer lane established six days ago, one directory over, that the substituted
family is **statistically identical to ΛCDM** (Δχ² = −0.30) and that 3.4–6.3σ *"came from forcing w₀
to DESI's central value — a point the likelihood never visits."* A superseded verdict left as the
page's only verdict while the superseding result sits unpublished in the same repo. Same shape as the
08-16 whitepaper finding — and with the ADMIN token dead, nothing can close it.

### 6. Terminology watch opened: MRH, protected in both whitepapers, reframed in both forums

§4.2 is radius-weighted (*"Spatial Boundary: the physical distance…"*, plus a femtometer→light-year
ladder). Phase 18 asserts *"MRH is not primarily a radius"* and offers hop depth, policy dimension,
evidence class, lineage depth, evaluator dependency and assurance threshold as coequal. Written into
`forum/` in two repos the same day.

Nothing is wrong today — both harvests are correctly marked non-canonical, they are in the ripening
lane, and §4.2 already hedges with *Context Dependence* / *Adaptive Horizons*. **Trip condition
stated so it is checkable rather than a vibe**: act if Web4's whitepaper adopts the contract reading
while §4.2 keeps the radius reading. Two whitepapers, one protected term, two definitions.

---

## Anomalies

- **08-17 Publisher run died at startup** — weekly limit, claude exit=1, 8-line log. Not a Publisher
  fault; fleet-wide. Recorded so the 08-16 → 08-18 skip is not misread as a missed review.
- **The whole 18-phase arc is one calendar day.** Not an anomaly to correct — the lane is allowed to
  work fast — but it is why the "arc still active" exclusion trigger fires so hard, and it is the
  mechanism behind Finding 2: production outran its own gate inside a single day.
- **Uncommitted, not mine**: `AGENTS.md` / `CLAUDE.md` GitNexus counter bumps (supervisor-owned,
  inside `gitnexus:keep` blocks) and `simulations/session373_acceleration_regime.png`. Adds scoped.

---

## So what

The single most encouraging thing this track has recorded ran today, and it arrived carrying its own
new defect. A lane ran the external prior-art gate unprompted, before promoting anything, and used it
to demolish its own novelty claims against cited literature — the behaviour fourteen prior REC-038
instances existed to demand, from nobody's instruction. That also prices the gate at one pass, which
removes the last excuse from a 15-day blocker on REC-040.

And then the same day's work outran it. The gate covered twelve phases; six more were produced past
it; the harvests claim all eighteen and carry the verdict with none of the citations. Which is the
same shape as my own watch item armed the day after its condition was met, and the same shape as the
site page still publishing a verdict its own repo superseded six days ago.

**A gate is an event, not a property.** Running it once buys coverage of what existed when it ran,
and every document downstream will quietly inherit it as a standing guarantee unless the coverage is
stated. The correction propagating without its evidence is the same failure wearing the opposite
face: the *conclusion* is the part that travels cheaply, and the part that makes it checkable is the
part that gets dropped. Four instances today, in four different lanes, one of them mine.
