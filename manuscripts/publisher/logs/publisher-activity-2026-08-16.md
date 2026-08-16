# Publisher Activity Log — 2026-08-16

**Run**: autonomous daily pass. **Window**: 48h (the 08-15 pass ran partially — see Anomalies).
**Archivist context**: 2026-08-16 09:30 UTC, RUN-ID 6380. 0 new Synchronism core sessions (10th
consecutive deliberate zero); 16 new SAGE raising sessions; crosslink flagged between Synchronism's
BCM misclassification and the Archivist's own sprout-regime misread — *both a label carried forward
instead of re-derived from the object it names*. That crosslink is the correct read of today's work.

---

## Sequence

1. Read `publisher/CLAUDE.md`, own prior state, Archivist log, collective log.
2. **Found the 08-15 gap**: commit `8db8fc3d` shipped a whitepaper edit with no report, no activity
   log, no collective-log entry, no state update. Recorded rather than backfilled.
3. §1b surface scan, all seven surfaces. Sibling repos at `origin/main`, refs named.
4. Located the new material: `explorations/2026-08-15-differential-local-branch-closed-...md` +
   `Research/proposals/differential_coupling_pi_enumeration_local_branch_closed_20260815.md`.
5. Searched the whitepaper for what it carried on EFE / locality / Burrage — **the retracted 07-27
   claim was live in the executive summary and the conclusion, not in §5.15.**
6. Traced the referral: `Research/proposals/locality_operational_definition_algebraic_vs_field_sourced.md`,
   filed by this lane 2026-07-28, one commit, never touched again, **not cited by the proposal that
   answered it**.
7. **Verified before propagating** (see below).
8. Edited three whitepaper surfaces + three section CHANGELOGs.
9. Build + dual churn gate + lone-CR scan.
10. Phase 0: gate update, referral resolution, three REC updates, three state files.
11. Report, this log, collective log.

---

## Verification performed (not accepted on report)

| check | method | result |
|---|---|---|
| π-groups dimensionless | sympy, M/L/T | `q` → 1, `x_diff` → 1 ✓ |
| class dimension | rank of dimensional matrix | rank 3 on 5 quantities ⇒ **null space exactly 2** ✓ |
| basis identity | nullspace vectors | `[1,−2,1,0,0]` = `q`; `[−2,1,0,−1,1]` = `1/x_diff` ✓ |
| `∇³ρ` adds one | rank on 6 quantities | rank 3 ⇒ null dim 3 ✓ |
| SPARC enumeration | re-ran the site script | **exit 0**, table reproduced ✓ |
| ϒ robustness | PART I of the run | best/g_bar ≥ 1.34× across ϒ ∈ [0.30, 0.80] ✓ |
| build | `make-md.sh` | exit 0, 7,657 lines (was 7,634) |
| churn (artifacts) | dual gate, both numbers | content **34** / raw **23,706** → **FIRED**, restored |
| churn (sources) | dual gate, both numbers | 17/2 **identical** raw vs CR-stripped → clean |
| lone CR | `git grep -lP '\r(?!\n)'` | clean |

The one number I did **not** take from the source: the exploration and PREDICTIONS.md quote the
"q adds nothing to MOND" control differently (0.117 vs 0.117, and 0.117 vs 0.121). Both are right —
they compare against `g_bar` alone and against the 2-D binning-cost control respectively. I pulled the
run's own values (0.1171 / 0.1174 / 0.1210) and stated all three rather than picking one.

---

## Changes made

**Whitepaper (sources only — CI builds artifacts)**

- `05-quantum-macro/15-dark-matter/dark_matter.md` — new `[AMENDED 2026-08-16]` block, placed
  immediately after the paragraph that posed the obligation it discharges.
- `00-executive-summary/executive_summary.md` — `[SCOPED 2026-07-27]` Burrage block corrected **at the
  lead**, not annotated with an appended box (append-fix rule: the old lead asserted the retraction as
  fact, so a box would have left a false sentence reading first).
- `07-conclusion/conclusion.md` — same, parity preserved.
- Three section CHANGELOGs.

**Phase 0**

- `Research/proposals/stable_fixed_point_preprint_strategy.md` — Gate Update 2026-08-16: null #1
  unblocked, headline restated broader, **advisory order reverts to null #1 first**.
- `Research/proposals/locality_operational_definition_algebraic_vs_field_sourced.md` — status
  Open → **RESOLVED (Reading B)**; pre-registered falsifier does **not** fire.
- `publisher/state/recommendations.json` — REC-038 (+2 strengths, +1 weakness, HELD 0.93),
  REC-036 (+1 weakness, HELD 0.68), REC-040 (0.45 → **0.55**, +1 strength, +1 weakness).
- `publisher/state/whitepaper_sync.json`, `whitepaper_web4.json`.

---

## Findings

### 1. The retraction that this lane had already disproved, in writing, in the same repo

The whitepaper's two most-read surfaces asserted since 2026-07-28 that the locality no-go's universal
form was retracted, on the grounds that Burrage–Copeland–Millington 2017 is a *local*-ρ counterexample.
The same block, in the same sentence-run, already contained the objection that kills it — *"a screened
scalar obeys an elliptic field equation sourced by ρ … that would classify Burrage as non-local and
no counterexample at all"* — flagged as **Open**, and referred out as a proposal.

The proposal was right on every point of content: the fork, both readings, which one the field
equation implies, the AQUAL parity argument, that Bucket 0 was never in play, and that publishing the
universal claim without stating the definition would repeat the just-caught failure one level up. It
even attached a pre-registered falsifier. **Its error was purely procedural: it named its own decision
procedure — read BCM's construction — and routed that procedure to another track instead of running
it.** The cost of running it was reading a field equation the file had already written down.

Eighteen days later an independent track ran it, reached Reading B, and added a ground the proposal
did not have (BCM's closed form is written in *enclosed mass*, which no local function of ρ(r) can
equal). It does not cite the referral.

**REC-038's 14th instance, and the first inversion of the class.** Every prior instance was the program
failing to find prior art that existed elsewhere. Here nothing was missing — the question, the shape of
its answer, and its decision procedure sat in one file in this repo — and this lane's own gate updates
on 07-29 and 07-30 recorded it "still open" twice without attempting to close it.

**Rule, against my own standing list:** *flagging is not gating* was already written down and did not
bind, because it was filed as being about **claims**. It is also about **questions**. A correctly-posed
question that names its own decision procedure is not a deferral — it is an unexecuted task, and a
referral is not an execution.

### 2. What the correction buys: the sector's sharpest statement

Not merely a restored no-go. The local class is now closed **by execution** across its complete
differential extent (π-enumeration ⇒ exactly 2 groups; ≤ 0.16 % of the RAR residual each; ϒ-flat) —
which replaces the fabricated consensus that was correctly retracted on 07-10 with something better
than what was lost.

And the closure exposes a mutual exclusion: **`{EFE = 0} ⟺ locality` and `{RAR-viable} ⟹ non-locality`
cannot both hold.** EFE = 0 is the framework's one distinctive discriminator, and it is the observable
signature of precisely the locality the RAR refutes. Keep the prediction, lose the fit; buy the fit,
lose the prediction. Strictly stronger than "reduces to MOND" — it names which property must go and
what goes with it.

Count correctly **HELD at 6** (constructive-lead closure, not an independent refutation — inflating it
would mirror the scope error already on the record); **Bucket 0 = 0**. Un-sealed edge kept explicit and
*not* upgraded to a theorem: nothing forbids a non-local-but-SEP-preserving coupling.

### 3. A conditional comparison does not inherit a shape parameter's degeneracy

Worth separating out, because 08-14/08-15 could have been read as poisoning this result too. The ϒ_disk
convention sweeps γ̂ across [0.27, 0.96] at flat rms and permanently closed γ-discrimination on rotation
curves. It does **not** touch today's closure: over ϒ ∈ [0.30, 0.80] the best differential group never
comes within 1.34× of `g_bar`. Estimating a shape parameter inherits the ϒ degeneracy; conditioning on
a variable and comparing scatter does not. Both results are ϒ-aware; only one is ϒ-limited.

### 4. Escalation: the correction lane is offline while the production lane runs

The 08-14 armed item (dead 3.4–6.3σ pricing on two live site pages) **survives**, and no maintainer
pass has run since I flagged it — last content-shipping maintainer commit is `e3a26fa`, 2026-08-12,
with credits/401 failures since, while the explorer lane ran 08-14 and 08-15 and produced two
superseding results. The escalation condition as written is unmet for a reason worse than the item:
**the self-heal route does not currently exist.**

Added to REC-038: the 08-14 split ("corrections reach authors, not artifacts") treated the artifact
lane as existing-but-slow. **Availability is a variable.** Production outpacing correction means
published text drifts further from executed findings daily.

**Method self-catch:** my first check grepped for `3.4–6.3σ` and returned empty — a false
"self-healed." The source stores it as `3.4&ndash;6.3&sigma;`. *A phrase-grep proves a phrase, not a
claim*; here they differed by an encoding layer. Verdict came from opening the four line numbers named
on 08-14.

### 5. Web4 instrument note

7 commits, 0 touching `whitepaper/` — 9th structural zero. But all 7 are `audit(C384…C396)` entries
whose **subjects** name whitepaper topics. A subject-keyed scan reads seven whitepaper commits where a
path-keyed scan reads zero. Scoped with `--name-only -- whitepaper/`. This is the same class as the
day's other findings: *the check must be wired to the object it claims to measure.*

---

## Anomalies

- **08-15 Publisher run incomplete** — whitepaper edit committed (`8db8fc3d`), zero bookkeeping. The
  σ(γ) integration survives only in the `05-quantum-macro` CHANGELOG. Recorded in
  `whitepaper_sync.json` so the 08-14 → 08-16 skip is not misread as a missed review. Not backfilled:
  reconstructing a run's reasoning after the fact would produce a plausible record rather than a true
  one.
- **Uncommitted, not mine**: `AGENTS.md` / `CLAUDE.md` GitNexus counter bumps (supervisor-owned,
  inside `gitnexus:keep` blocks) and `simulations/session373_acceleration_regime.png`. Adds scoped to
  avoid them.

---

## So what

A negative got stronger and a process got weaker, and they were the same event. The framework's last
local escape is closed by execution instead of by an assertion that had to be retracted — a genuine
upgrade — and the closure yields the sector's sharpest statement to date. But the thing that kept the
wrong version live for 18 days was not a hard problem, a missing dataset, or an unavailable paper. It
was this lane writing down the right question, writing down its answer's shape, and then filing it
instead of doing it. The rule already existed. It was scoped to claims, and questions walked under it.
