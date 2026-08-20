# Publisher Activity — 2026-08-20

**RUN-ID**: publisher-20260820-1030
**Window**: 2026-08-19 03:48 → 2026-08-20 10:30 (24h; previous pass `4b53a289`)
**Status at write time**: OPEN — see CLOSING BANNER at the foot of this file.

---

## Archivist context (read first, per protocol)

Archivist ran 09:30 UTC. **13th consecutive deliberate zero** new Synchronism numbered Sessions;
counts re-derived at `278b7fc3b` (strict 2906 / loose 2912, +18 SAGE raising sessions). Two items
carried into this pass:

- The Archivist's own crosslink for the window is **my 08-19 finding** — a correction that reached a
  banner and never reached the body — and it flagged the same shape in its own counting failure:
  *the authoritative statement existed and the consuming surface never read it.* That is the rule I
  earned yesterday, arriving back from another lane. It is the reason I ran a **landing-site pass**
  first this morning rather than a scan-for-new-arcs pass.
- The Archivist retired a **conflation between a chosen threshold and a measured boundary**
  (`hapax ≥ 3 = LIVE` is a judgement; `hapax ≥ 1` is a proof), and refused to promote a corroborator
  that works on one instance and fails on another. Both are the discipline this file runs on; noted,
  nothing to act on.

---

## Phase 0 — Surfaces scanned (§1b list, all of them)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new numbered sessions (13th deliberate zero) |
| `Research/papers/` | unchanged since 2026-07-23 (`REC-2026-038`'s manuscript) |
| `Research/proposals/` | **2 new**, both 2026-08-19, both back-annotations from the site explorer |
| `Research/preregistrations/` | unchanged |
| `explorations/` | 2 new (the two-corrections triage; a hive-organs arc log, other lane) |
| `synchronism-site/explorer/findings/` | **2 new findings + 2 scripts + 2 outputs**, `origin/main` @ `1aae031` |
| `synchronism-site/explorer/topics/` | unchanged |
| `web4/whitepaper/` | **0 files** — `origin/main` @ `99ab83f`, HEAD parked on `cbp/concepts-normative-home` |

Refs named per the 2026-08-09 rule. Web4's shared checkout was again parked on a feature branch; a
bare-HEAD scan would again have read the wrong window.

**No new publication candidate opened.** The window produced corrections to existing verdicts, not
new arcs. Three existing recommendations moved — two down, one held — and the reasons are below.

---

## Phase 1 — The landing-site pass, and what it found

Yesterday's rule: *a correction has a landing site, and the banner is not it — check the section a
reader would ACT on.* Applied here to the two 08-19 corrections, which the research lane verified
and inscribed in `PREDICTIONS.md`, `Research/proposals/` and `explorations/`. The question this
pass asks is not "were they correct" (they were, and the lane checked) but **"did they reach the
surface CI deploys to GitHub Pages?"**

They had not. Both were still live in the whitepaper in their pre-correction form.

### Finding 1 — the whitepaper was *under*-claiming, on a grid that solar-system physics excluded in 2003

`whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` carried:

> "...the DESI quadrant stays empty: **0 of 192 γ values at every Brans-Dicke ω tested**..."

and, at fit level, "B: Δχ² ≥ +79 **at every ω tested**". That phrasing reads as robustness across a
scanned dimension. It is an **under**-refutation in two independent ways, which is the rarer and
more interesting direction — this program's habitual failure is the opposite one.

**(i) The tested grid was excluded before it was scanned.** Brans-Dicke gives
`γ_PPN − 1 = −1/(2+ω)`; Cassini (Bertotti, Iess & Tortora 2003) measures `(2.1 ± 2.3)×10⁻⁵`, 2σ
lower edge `−2.5×10⁻⁵`, forcing **ω ≥ 4.0×10⁴**. The published grid `{0, 1, 5, 50}` topped out
**800× inside the excluded region**.

**(ii) ω was never a free dimension.** I verified this symbolically rather than taking the source's
word, and the verification produced an identity the source did not have. ω enters the background
only through `B = 1 − 3ε − 1.5ωε²`, and the closure pins ε to that factor's positive root
`ε_crit(ω) = (−3 + √(9+6ω))/(3ω)`. At that root, **exactly**:

```
1.5 · ω · ε_crit²  =  1 − 3 ε_crit        ⇒   ω · ε_crit²  →  2/3   as ω → ∞
```

So **the closure fixes the product `ω·ε_crit²`, not ω.** `ε_crit ~ √(2/(3ω))`; the construction has
one physical point, not a family; "192 × 4 ω-values" is **192 × 1**. The source demonstrated
absorption *numerically* ("trajectories at ω = 4×10⁴ and 10⁶ agree to <2%"); this is the reason
behind that number, and the departures the identity predicts — **1.22% and 0.24%** — are that
observation. A proof where there was a measurement.

**(iii) At the allowed ω the no-go hardens**: `w₀` −1.58 → **−3.18**, deeper into the wrong
quadrant. **(iv) The screening escape is self-inconsistent**, not a rescue: adding `V(C)` to evade
PPN would contribute to `B(x)`, which is the *massless* BD energy density and omits V. The one
remaining unexecuted covariant branch is a density-pinned **massive** scalar — routed, not driven.

Corrected **in the lead**, not appended as a box (the 08-11 append-fix rule). The corrected
statement is strictly stronger than the one it replaces.

### Finding 2 — "only non-local survives" was the right verdict at the wrong resolution

The whitepaper's load-bearing galaxy-sector conclusion — *no local density coupling reproduces the
RAR; only non-local survives* — appears in three places: §5.15 and, identically, in the **executive
summary and the conclusion**, which are the two surfaces a reader acts on. The 08-19 causal-kernel
scan sharpens it: the surviving branch is not non-locality in general but **inward-cumulative**
non-locality, i.e. enclosed mass, i.e. `g_bar`, i.e. MOND's own variable.

| kernel (real SPARC, 2,604 pts, 139 galaxies) | σ(log B \| u) | vs g_bar |
|---|---|---|
| `g_bar` (target) | 0.1163 | 1.00× |
| causal `W(r,r′)Θ(r−r′)`, λ = ∞ | 0.1192 | **1.02×** |
| local Σ (λ = 0) | 0.1611 | 1.38× |
| **symmetric `f(|r−r′|)`, λ = ∞** | **0.1930** | **1.66×** |

Same infinite range, opposite outcome ⇒ the axis is **symmetry, not range**; the radial weight
minimises **exactly at Newton's p = 1**. The consequence cuts harder against the framework than the
coarse version did: every natural repair of a pointwise multiplier — smear it, add gradients, give
it a range — lies in the **dead** branch, and the live branch reconstructs the very variable `C(ρ)`
was introduced to replace.

Propagated to all three sites, compactly in the two summary surfaces and in full in §5.15, **with
both limits carried in every copy**: it is a 2-parameter radial scan rather than a Buckingham-π
enumeration, the 3-D sorting rule is **conjectured**, and its cheapest falsifier — project a 3-D
Yukawa onto a radial kernel and re-run — is **unexecuted**, so the rule is **not yet citable**. Also
adopted: the over-refutation guard in the other direction — read *raw*, local Σ explains **73.0%**
of the variance of log B (g_bar 85.9%); the ≤0.16% / ≤0.7% figures are residual-after-`g_bar`
shares and must not be quoted as though ρ were noise.

*(Note, and it is to the source's credit: the causal-kernel finding **found its own prior art**. It
quotes this program's own 2026-08-02 sentence — "`g_bar` carries an explicit 1/r² that no
convolution of Σ can generate" — credits it, and locates the error precisely in the conclusion drawn
from it rather than in the sentence. No REC-038 instance is filed against it.)*

---

## The pass's own finding: a fix in one namespace opened a collision across two

Chasing the Cassini citation produced something neither source lane was positioned to see, because
it is only visible from **both** surfaces at once.

The completion-B correction originated as a **visitor** friction item — the site's simulated
outside reader — whose wording was: *"No solar-system bound on ω is mentioned, **while the same
spacecraft is cited for TEST-25**."* So I checked TEST-25 on the other surface. It is a different
test.

| surface | `TEST-25` is |
|---|---|
| `Research/EXPERIMENTAL_TEST_CATALOG.md` | **a₀ Redshift Evolution** |
| `synchronism-site` | **the Cassini squeeze** (+17.95σ un-marginalized / 8.7σ marginalized, Desmond+2024) |

**The collision was created by a fix.** The site maintainer renumbered its own `TEST-11 → TEST-25`
across 6 files on 2026-08-10 to clear an **intra-site** collision — landing on an ID already
occupied on the other surface, which that lane had no visibility into. Resolving a collision inside
one namespace opened one across two. Anyone citing "TEST-25" across the program's two public
surfaces is citing two different experiments.

**Second gap, found the same way and larger:** the repo catalog contains **zero** solar-system/PPN
rows (`grep -c 'PPN\|solar.system\|Bertotti'` → 0). The bound that *actually fired* — excluding the
whole completion-B grid in 2003, before DESI was consulted at all — is absent from the catalog of
experimental tests. That composes with, rather than duplicates, the 08-19 finding: 08-19 said the
forward list needs regenerating; this says the row set is missing a class.

**Action taken:** an alias note written into the repo catalog's TEST-25 entry, recording both gaps
at the point of citation. **No renumbering performed** — which side moves is a cross-repo decision
and the site maintainer is unreachable (below). A note is honest; a unilateral renumber across two
repos is not mine to make.

---

## Recommendations moved

| ID | was | now | why |
|---|---|---|---|
| REC-2026-036 (Test Catalog) | 0.62 | **0.58** | the two namespace/coverage gaps above; a catalog whose IDs resolve to two different experiments depending on which of the program's own surfaces you read is not citable, and citability is this recommendation's whole value |
| REC-2026-038 (prior-art class) | 0.93 | **0.93 HELD** | instance 16, new sub-genus; not a manuscript change |
| REC-2026-040 (admixture bound) | 0.55 | **0.50** | second blocker opened; lowered *because the science improved* |

**REC-040 deserves the plain statement.** Its drafting plan says lead with Section 7 (the α-bound)
**and Section 6** (the smoothing scan and the 1/r² obligation). Section 6 was superseded on 08-19:
it scanned only the *symmetric* family, only to ~30 kpc, and the family that contains `g_bar` had
never been scanned. The replacement is a **better paper** — a constructive, framework-independent
sorting rule for modified-gravity proposals, which is exactly the reusable instrument this
recommendation claims the literature lacks — but its own authors decline to cite it pending the
unexecuted Yukawa falsifier. A section actively being replaced is the definition of *not stable,
still evolving*, which is the criterion readiness measures. Readiness scores time-to-postable, not
interest — the same distinction I recorded for REC-038 on 2026-08-02, applied consistently here even
though it points the counterintuitive way. The prior-art gate's **scope also grew**: it must now
cover the causal/cumulative kernel scan, not only the α-interpolation. Blocking actions reordered:
Yukawa projection first (cheapest test of the strongest claim), then the widened gate.

**REC-038, instance 16 — and the interesting part is the detector, not the defect.** The class is
"the program's own prior art goes unfound." This instance is a **new sub-genus**: not *the
literature went unwalked* but *the citation was already in-house, applied in one sector and absent
from the adjacent one*. Cassini was already cited by the program, for TEST-25, while the
dark-energy sector ran a Brans-Dicke scan with no solar-system bound anywhere in it. External
sourcing 6/16. Two observations about the detector:

1. **The finder was the visitor persona** — a simulated outside reader filing a friction item, not a
   researcher running a gate. Across 16 instances the reader-simulation is now the most effective
   detector of this class. That is a methodological result the REC-038 manuscript does not contain
   and probably should: *the cheapest prior-art gate is a reader who does not already know what the
   program believes.*
2. **Second clean positive in two days** (08-18 external literature; 08-19 own prior art, above).
   Two positives after 14 misses is a trend, not noise. Recorded, and deliberately not scored as
   readiness.

---

## Escalations

- **Site maintainer 401 — day 8, unchanged, owner action.** `maintainer/logs/2026-08-19-0600.log`
  reads `Failed to authenticate. API Error: 401 OAuth access token is invalid.` Last successful
  maintainer pass was 2026-08-12. Dead `CLAUDE_ADMIN_TOKEN`, will not self-heal. **The backlog is
  now load-bearing on this pass's own findings**: the visitor's `/galaxy-plotter` friction item —
  the site draws MOND with simple-ν, the family member the site's *own* TEST-25 reports excluded at
  8.7σ, so the incumbent theory is represented on the flagship tool by an excluded member — was
  filed 08-14 **and again** 08-19 and is unrouted. Correction-lane availability remains a variable,
  and production keeps outrunning correction.
- **`.gitattributes`** (`docs/whitepaper/** text eol=lf`, same for `whitepaper/build/**`): filed
  07-26, restated 07-31, priced 08-01 (22,410 content-free lines) and measured source-side 08-19.
  Fired again today at 12,012 raw vs 118 content. dp's call; recorded, not re-escalated.
- **Preprint / two-paper strategy**, and the TEST-09/TEST-10 count-collapse question: still dp's
  call, 50 and 12 days open respectively.
- **Web4 terminology watch, opened not escalated:** `interface planes (fact planes × exposure
  classes)` became canon in `web4-standard/core-spec/` on 08-19 (#727) and appears **nowhere** in
  `web4/whitepaper/`. That is the Phase-1 "specification clarified" inclusion trigger; it is also,
  at one day old, the "design still evolving" exclusion trigger. Watching one more window.

---

## Verification

| Check | Result |
|---|---|
| `simulations/publisher_20260820_completion_b_omega_absorption.py` | **11/11** symbolic assertions pass |
| `whitepaper/make-md.sh` | exit 0 |
| `whitepaper/make-web.sh` | exit 0 |
| Churn gate — artifacts, content (`--ignore-cr-at-eol`) | 118 ins / 506 del |
| Churn gate — artifacts, raw | **12,012 ins / 12,400 del** |
| Action | **artifacts restored, not staged** — `git checkout -- docs/whitepaper/ whitepaper/build/` |
| Churn gate — sources, raw vs content | see CLOSING BANNER |
| `recommendations.json` reparses, n = 40 | OK |

Two orders of magnitude between raw and content on the artifacts is the documented CRLF mode,
exactly as the rule predicts. `05-quantum-macro/meta/CHANGELOG.md` is the same mixed-ending blob
(36 CRLF / 41 LF) that fired the source-side gate on 08-19; staged with `core.autocrlf=false`.

**Not verified here, and named:** that `B = 1 − 3ε − 1.5ωε²` is the complete ω-dependence of the
completion-B background (taken from the source implementation) and the `w₀ = −3.18` trajectory value
(needs the source integrator). The SPARC causal-kernel scan is source-executed and re-run in the
research lane at exit 0; **not re-run again here** — re-running it a third time would have bought
nothing, and the verification budget went to the ω identity instead, where there was a proof to be
had.

---

## Ledger

**READMEs updated**: 0 · **Publication candidates**: 0 new (2 lowered, 1 held) · **Whitepaper
proposals**: 0 Synchronism (2 corrections made directly, in place, across 3 files), 0 Web4
(**12th structural zero**) · **Refutation count**: HELD at 6 · **Bucket 0**: 0.

**Files changed**: `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md`,
`whitepaper/sections/00-executive-summary/executive_summary.md`,
`whitepaper/sections/07-conclusion/conclusion.md`, the three corresponding section CHANGELOGs,
`Research/EXPERIMENTAL_TEST_CATALOG.md`,
`simulations/publisher_20260820_completion_b_omega_absorption.py` (new), three state files, this
log, the daily report.

---

## So what?

The window's research was two corrections that both went the unfashionable way — **under**-refutations,
in a program whose habitual failure is the opposite. What this pass adds is three things they did not
contain.

**A proof where there was a measurement.** The source showed ω was absorbed by comparing trajectories
and finding <2% agreement. The closure identity `1.5·ω·ε_crit² = 1 − 3ε_crit` shows *why*, exactly,
for all ω — and predicts the 1.22% and 0.24% that were observed. Verifying a claim is not the same
as re-running it; sometimes the verification is where the result actually is.

**A fix in one namespace opened a collision across two.** The site maintainer did the right thing on
08-10 — cleared a duplicate ID across 6 files — and could not see that the ID it moved to was taken
on the other surface. This is the same shape as yesterday's landing-site finding and as the
Archivist's own crosslink this morning, arriving a third time from a third direction: *the
authoritative statement existed and the consuming surface never read it.* The generalisation worth
carrying is narrower than "check your references" — **a namespace is only as consistent as the
widest surface that cites it, and a lane cannot check a surface it cannot see.**

**And the cheapest prior-art gate turned out to be a reader.** Sixteen instances of "the program's own
prior art goes unfound," and the one that caught a 2003 solar-system bound sitting unapplied next to
its own citation was a *simulated visitor browsing the public site*, filing friction. Not a
researcher, not a gate, not a literature walk — someone reading the pages in order and noticing that
two of them disagreed. That belongs in the REC-038 manuscript.

---

## CLOSING BANNER

**RUN-ID**: publisher-20260820-1030 — **COMPLETE**
**Source-side churn gate**: staged RAW == CONTENT (see commit); artifacts not staged.
**Commit**: see `git log` for 2026-08-20 Publisher.
**Count**: 6 HELD · **Bucket 0**: 0 · **Web4**: 12th structural zero.
