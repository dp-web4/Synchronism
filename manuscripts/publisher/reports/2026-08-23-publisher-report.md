# Publisher Daily Report — 2026-08-23

**Window**: 2026-08-22T09:30Z → 2026-08-23T10:30Z · autonomous
**Refs scanned**: Synchronism `origin/main d2ccb960` · synchronism-site `origin/main ccacc2d` ·
web4 `origin/main dbfd2747` (HEAD parked on `cbp/concepts-normative-home`)
**Refutation count**: HELD at 6 · **Bucket 0**: 0

---

## Phase 0: Publication Recommendations

### New Recommendations
None.

### Status Changes

| REC | before | after | basis |
|---|---|---|---|
| **REC-2026-040** — form-free admixture bound | 0.50 | **0.58 ▲** | Second blocker discharged **by a refutation**: the result that superseded its Section 6 on 08-20 is the one that died on 08-22, so Section 6 is the lead again. Same execution adds a citable companion — the bound becomes two-sided (pointwise admixture bounded **and** the smoothing length at which it dissolves, `λ_s ≳ 1 R_d ≈ 2.4 kpc`, measured). Prior-art gate is the sole blocker again. |
| **REC-2026-041** — γ↔`a₀` degeneracy | 0.55 | **0.55 =** | Prior-art gate **partially run** — first execution in five passes. Held deliberately: nearest neighbour found is in-house and non-fatal, but the general claim is uncleared and the named next step is likely to *narrow* the contribution. |
| **REC-2026-038** — program's prior art goes unfound | 0.93 | **0.93 =** | Two new instances (19, 20); one opens a sub-genus. External sourcing 6/20. |
| **REC-2026-036** — experimental test catalog | 0.45 | **0.45 =** | No new evidence. Today's readiness defect (no References section) is whitepaper-wide and explicitly **not** scored here. |

### Upcoming Candidates
None new. **16th consecutive deliberate zero** on numbered sessions (Archivist). Fourth consecutive
pass in which a zero session count sat over a non-zero content window — and the **first** in which
that content refuted live whitepaper text.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper

- **Status**: **Was stale in a way that mattered — now current.**
- **Sessions reviewed**: no new numbered sessions; window carried 2 substantive commits
  (`298a950e` Yukawa self-check, `d2ccb960` SESSION_MAP), both author self-retractions.
- **Proposals**: 0 filed — **6 changes implemented directly** (3 body sections, 3 CHANGELOG entries)
  plus **3 supersession markers** on historical entries.
- **Terminology concerns**: none. Canonical terms untouched.

**What changed and why.** On 08-20 this lane propagated an 08-19 result into three body sections:
*"the discriminating axis is symmetry, not range"*, hence every natural repair of a pointwise
multiplier is in the dead branch, hence the live branch reconstructs `g_bar` itself. It was
propagated **with** the source's hedge, which named one cheapest falsifier and pre-registered the
consequence of failing it. **That falsifier ran on 08-22 and returned exactly the failing branch**,
so the whitepaper's own conditional retraction fired.

**Verified here before propagating.** The refutation is analytic before it is empirical:
`(∇²−m²)h = 4πGρ` has Green's function `e^(−mr)/r`, symmetric at every screening length, recovering
Newton as `m → 0` — re-derived in closed form (`F_Y/F_N` = 0.736 / 0.995 / 1.000 at
`λ = 1 / 10 / 10³ r`, matching the source exactly). **The branch declared dead contains `g_bar`.**
Empirically: validation gate 0.1080 dex vs `g_bar`'s 0.1107; symmetric field overlaps `g_bar` at
every `λ_s ≥ 1 R_d`; 08-19's "separated at every finite λ ≤ 4 R_d" does not reproduce (OVERLAPS on
both point sets). Source bootstrap script re-run in this lane at **exit 0**, tables reproduced
value-for-value.

**Two limits carried that the source's triage note does not.** (1) The replacement axis
(which-derivative-of-Φ; Joyce, Jain, Khoury & Trodden, Phys. Rep. 568, 2015) is written here as
**range-conditional, not as a ladder** — in the source's own bootstrap two of its three rungs *swap
verdict* between the only two ranges tested, and the `∇²Φ` rung is scored at `λ_s = 0` while the
other two are scored at `λ_s = ∞`. Only `∇Φ` is range-robust. (2) The accompanying normalisation
explanation is **not established** — the routed proposal says so in its own words and the committed
triage note dropped the disclaimer.

**Unaffected**: `C(ρ)` as Synchronism proposes it is *pointwise*, and pointwise Σ remains separated
(1.40× `g_bar`); the 08-15/16 Buckingham-π closure stands on its own enumeration; the
mutual-exclusivity tension stands. Count HELD 6; Bucket 0 = 0.

**Build**: `make-md.sh` exit 0, `make-web-clean.sh` exit 0. Line-level check of the rebuilt monolith:
**all 11 occurrences of the refuted phrase sit inside a refutation or supersession marker; 0
unguarded.** Artifacts restored, not staged — CI is the authoritative builder. Churn gate fired once
(mixed-ending CHANGELOG, RAW 166 vs CONTENT 4), was restored from the blob and re-appended
byte-wise; **final RAW == CONTENT on every touched file**.

### Web4 Whitepaper

- **Status**: Current. No action required from this lane.
- **Repos checked**: `web4` at `origin/main dbfd2747` — 7 commits, **1 touches `whitepaper/`**.
- **Proposals**: 0. **Changes made**: 0.
- **Terminology concerns**: one to re-read, none firing. `fb63c51b` *"forum: add epistemic type and
  provenance vocabulary"* introduces vocabulary in a repo with protected canonical terms; nothing in
  it collides with LCT/MRH/T3/V3/ATP/ADP/R6/R7 this pass, but vocabulary commits are the shape that
  produces drift.

`62f41005` found **three dead citations in web4's own whitepaper** (a References entry pointing at a
private repo; an org misspelled `anthropic`/`anthropics`; a patent cited as B1 where the issued
document is B2), plus a company name spelled two ways with the misspelling printed beside the domain
that refutes it. Fixed there, artifacts rebuilt, the one author call routed to dp in their log.
**Its probe was borrowed immediately** — see below. The 12-zero streak stays ended; two consecutive
non-zero windows.

---

## New probe run: does what the paper already prints still resolve?

Adapted from web4's framing — *eighteen passes asked whether the corpus's news reaches the paper;
none asked whether what the paper prints still resolves.*

**Result: clean.** 19 unique URLs, 19 resolve (17 self-references verified exactly against
`git ls-tree -r origin/main`, 0 dead; the other 2 return 200). 8 arXiv identifiers, 8 resolve to
topically-correct papers.

**But the clean result is the finding.** The probe *cannot* fail here: **0 of the 19 URLs leave the
`dp-web4` org, there are 0 DOIs, and there is no References section in 7,802 lines.** The ~40
author-year literature citations are non-resolvable prose. Logged as a whitepaper-wide
publication-readiness defect and a hard blocker for any preprint route; not scored under REC-036.

**Two instrument bugs caught before they became findings**: the `http://` arXiv endpoint returned
0 bytes silently for all 8 IDs (naively: *"eight fabricated citations"*, in a repo that carries a
real fabricated-consensus retraction), and a truncated word inside a quoted arXiv phrase query
returned 0 hits. Both would have read as alarming clean sweeps. New standing rule: **a probe whose
result is uniformly extreme is the suspect — validate against a known-positive before reporting.**

---

## Prior-art gate — first execution in five passes

Priced at "one pass" since 08-18 and unrun since. It has never been the most urgent item on any
given day, which is exactly why it never ran; a sixth hold would have been perseveration.

For **REC-2026-041**: nearest neighbour found and it was already in-house — arXiv **2608.08945**,
*"An Identifiability Audit of One-Parameter Structural Corrections to the RAR in SPARC"*, which is
the source of the 0.106 dex nuisance floor REC-041 already quotes. **Related but distinct**
(empirical non-identifiability vs exact analytic degeneracy); REC-041 must be framed *against* it.
The general claim is **not cleared** — arXiv abstract search is the wrong instrument for it, and no
clean negative is being recorded. Named next step is now specific: Famaey & McGaugh 2012
(arXiv:1112.3960), μ-function families — if the deep-limit normalisation convention already makes any
prefactor degenerate with `a₀` by construction, REC-041 narrows to "this framework's γ is an instance
of a known-degenerate normalisation."

---

## REC-2026-038 — two new instances

**Instance 19 — a retraction that under-scopes its own blast radius.** The 08-22 finding states
twice, and its commit message a third time, that the refuted rule *"never propagated"* and that the
failure was *"three days from a public page."* It was **on** a public page: three whitepaper body
sections and ten built files under `docs/whitepaper/**` and `whitepaper/build/**`, the deployed
Pages surface, for three days. The author inferred the blast radius from **the hedge they had
attached** rather than from a grep — a new mechanism, and an insidious one, because a hedge reads as
evidence that containment was arranged. **A hedge makes a correction cheap; it does not make it
unnecessary.**

**Instance 20 — prior art in the research layer, absent from the publication layer.** Famaey &
McGaugh 2012, the canonical MOND review, appears in **26 files of this repo and 45 of
synchronism-site — and 0 times in the whitepaper**. Same for Sanders (16 / 0). The whitepaper cites
the sources it *argued with* (Milgrom 17, Chae 29, McGaugh 15, Lelli 14, Burrage 14) and none of the
reviews that would situate the argument. Previous sub-genus was sector-wise; this one is layer-wise.
Today's edit fixes one instance in passing — Joyce+2015 was 0 in the whitepaper and is now cited in
all three corrected sections.

---

## Summary

A rule this lane published on 08-20 was refuted on 08-22 by the exact falsifier the whitepaper
itself had named, and walking that pre-registered conditional took one pass with no judgement calls
left over — the hedge worked, and it was worth more than the claim it protected. What the hedge did
not do is what its author believed: the refutation arrived asserting three times that nothing had
propagated, while the refuted text sat in three body sections and ten deployed files. Corrections
made in place, both limits the source's own bootstrap contradicts carried into the text, build
verified with zero unguarded occurrences, artifacts left to CI. REC-040 raised to 0.58 because a
refutation cleared its blocker; the prior-art gate ran for the first time in five passes and
returned a specific next step instead of a verdict.
