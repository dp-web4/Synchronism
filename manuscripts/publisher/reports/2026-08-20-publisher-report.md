# Publisher Daily Report — 2026-08-20

**Window**: 2026-08-19 03:48 → 2026-08-20 10:30 (24h) · **RUN-ID** `publisher-20260820-1030`
**Full log**: `publisher/logs/publisher-activity-2026-08-20.md`

---

## Phase 0: Publication Recommendations

### New Recommendations
**None.** The window produced corrections to existing verdicts, not new arcs. 0 new numbered
sessions (Archivist: 13th consecutive deliberate zero).

### Status Changes

| ID | Was | Now | Reason |
|---|---|---|---|
| **REC-2026-036** — Experimental Test Catalog | 0.62 | **0.58** | Cross-surface `TEST-25` ID collision + zero solar-system rows (below) |
| **REC-2026-038** — prior-art class | 0.93 | **0.93 HELD** | Instance 16 recorded, new sub-genus; no manuscript change |
| **REC-2026-040** — admixture bound | 0.55 | **0.50** | Second blocker opened; lowered *because the science improved* |

**REC-2026-040 — the counterintuitive one.** Its drafting plan says lead with Section 7 (the
α-admixture bound) *and* Section 6 (the smoothing scan and the "1/r² obligation"). Section 6 was
superseded on 08-19: it scanned only the **symmetric** kernel family, only to ~30 kpc, and the
family containing `g_bar` had never been scanned. Scanned now, the causal family reaches **1.02×**
`g_bar` at λ = ∞ with the radial weight measured **exactly at Newton's p = 1**, while the symmetric
family at the *same* infinite range is 1.66× — worse than reading ρ pointwise. So the discriminating
axis is **symmetry, not range**, and Section 6's framing is the wrong shape.

The replacement is a *better* paper — a constructive, framework-independent sorting rule for
modified-gravity proposals, which is exactly the reusable instrument this recommendation claims the
literature lacks. But its own authors decline to cite it: a 2-parameter radial scan rather than a
π-enumeration, 3-D behaviour **conjectured**, and the cheapest falsifier (project a 3-D Yukawa onto
a radial kernel, re-run) **unexecuted**. A section actively being replaced is the definition of *not
stable, still evolving* — which is the criterion readiness measures. Readiness scores
time-to-postable, not interest; same distinction recorded for REC-038 on 2026-08-02, applied
consistently here. **The prior-art gate's scope also grew** (must now cover the cumulative-kernel
scan), so it is no longer the sole blocker. Blocking actions reordered: **Yukawa projection first**,
then the widened gate, then a second variable, then authorship.

**REC-2026-038 — instance 16, and the detector is the story.** New sub-genus: not *the literature
went unwalked* but **the citation was already in-house, applied in one sector and absent from the
adjacent one**. Cassini/Bertotti+2003 was already cited by the program for TEST-25 while the
dark-energy sector ran a Brans-Dicke scan with no solar-system bound in it at all. External sourcing
6/16. Two observations:
- The finder was the **visitor persona** — a simulated outside reader filing friction, not a
  researcher running a gate. Across 16 instances the reader-simulation is now the most effective
  detector of this class. *The cheapest prior-art gate is a reader who does not already know what
  the program believes.* Recommended for the manuscript.
- **Second clean positive in two days.** The 08-19 causal-kernel finding found its **own** prior art
  unprompted — quoted this program's 2026-08-02 sentence verbatim, credited it, and located the
  error in the conclusion drawn from it rather than in the sentence. Two positives after 14 misses
  is a trend. Not scored as readiness.

### Upcoming Candidates
None new. The unexecuted 3-D Yukawa projection is the highest-value cheap execution in the queue —
it is the cheapest test of the strongest live claim and it gates REC-040's constructive half.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: **Updated** — 2 corrections, made directly and in place, across 3 files
- **Sessions Reviewed**: no new numbered sessions; reviewed against the 2 new `Research/proposals/`
  back-annotations and 2 new site-explorer findings
- **Changes Made**:
  1. **Completion-B ω — an *under*-refutation corrected in the lead.** The whitepaper read "0 of 192
     γ values **at every Brans-Dicke ω tested**," which reads as robustness across a scanned
     dimension. It is not. *(i)* Cassini forces **ω ≥ 4.0×10⁴** (`γ_PPN − 1 = −1/(2+ω)`,
     Bertotti+2003), so the published grid `{0,1,5,50}` sat **800× inside the excluded region**.
     *(ii)* ω is **absorbed by the model's own closure** — verified here symbolically, 11/11: at the
     pinned point `1.5·ω·ε_crit² = 1 − 3ε_crit` **exactly**, so the closure fixes the product
     `ω·ε_crit² → 2/3` rather than ω. `192 × 4` is **192 × 1**. *(iii)* At the allowed ω the no-go
     **hardens**: `w₀` −1.58 → **−3.18**. *(iv)* The screening escape is self-inconsistent with the
     massless background; a density-pinned **massive** scalar `V(C)` is the one unexecuted covariant
     branch. **The corrected statement is strictly stronger than the one it replaces.**
  2. **"only non-local survives" sharpened to inward-cumulative.** Propagated to §5.15 *and* to both
     reader-action surfaces (executive summary, conclusion). Symmetric kernels fail at **any** range
     (σ = 0.193 dex at λ = ∞ vs pointwise ρ's 0.161); the causal kernel at the same range reaches
     1.02× `g_bar`, radial weight exactly p = 1. Every natural repair of a pointwise multiplier lies
     in the **dead** branch; the live branch reconstructs the variable `C(ρ)` was introduced to
     replace. **Both limits carried in all three copies** (2-param scan; 3-D conjectured; Yukawa
     falsifier unexecuted ⇒ **not yet citable**), plus the over-refutation guard: read raw, local Σ
     explains **73.0%** of log-B variance — the ≤0.16% figures are residual-after-`g_bar` shares,
     not a claim that ρ is noise.
- **Builds**: `make-md.sh` exit 0, `make-web.sh` exit 0. Artifacts **restored, not staged** (churn
  gate: content 118/506 vs **raw 12,012/12,400** — the documented CRLF mode). CI is the builder.
- **Terminology Concerns**: none. No canonical term redefined.

### Web4 Whitepaper
- **Status**: Current — **12th consecutive structural zero**
- **Repos Checked**: `web4` at **`origin/main` @ `99ab83f`** (HEAD parked on
  `cbp/concepts-normative-home` — ref named per the 08-09 rule; a bare-HEAD scan would have misread
  the window). 4 commits, touching `docs/audits/` and `web4-standard/`; **0 files under
  `whitepaper/`**.
- **Proposals**: None
- **Terminology Concerns**: **Watch opened, not escalated.** `interface planes (fact planes ×
  exposure classes)` became canon in `web4-standard/core-spec/` on 08-19 (#727) and appears
  **nowhere** in `web4/whitepaper/`. That is the "specification clarified" inclusion trigger and, at
  one day old, also the "design still evolving" exclusion trigger. Watching one more window.

---

## The pass's own finding: a fix in one namespace opened a collision across two

Chasing the Cassini citation surfaced something only visible from **both** surfaces at once. The
visitor's friction wording was *"no solar-system bound is mentioned, **while the same spacecraft is
cited for TEST-25**"* — so I checked TEST-25 on the other surface. It is a different test.

| surface | `TEST-25` is |
|---|---|
| `Research/EXPERIMENTAL_TEST_CATALOG.md` | **a₀ Redshift Evolution** |
| `synchronism-site` | **the Cassini squeeze** (+17.95σ un-marginalized / 8.7σ marginalized) |

**The collision was created by a fix.** The site maintainer renumbered its own `TEST-11 → TEST-25`
across 6 files on 2026-08-10 to clear an *intra-site* collision — landing on an ID already occupied
on the other surface, which that lane had no visibility into. Anyone citing "TEST-25" across the
program's two public surfaces is citing two different experiments.

**Second gap, larger:** the repo catalog has **zero** solar-system/PPN rows. The bound that actually
fired — excluding the whole completion-B grid in 2003, before DESI was consulted — is absent from
the catalog of experimental tests. This *composes with* the 08-19 finding rather than duplicating
it: 08-19 said the forward list needs regenerating; this says the row set is missing a class.

**Action**: alias note written into the repo catalog at the point of citation. **No renumbering** —
which side moves is a cross-repo decision and the site maintainer is unreachable.

---

## Escalations for dp

1. **Site maintainer 401 — day 8, owner action, will not self-heal.** Dead `CLAUDE_ADMIN_TOKEN`;
   last successful pass 2026-08-12. The backlog now touches this pass's own findings: the visitor's
   `/galaxy-plotter` item — the site draws MOND with simple-ν, the family member its *own* TEST-25
   reports excluded at 8.7σ, so the incumbent theory is represented on the flagship tool by an
   excluded member — was filed 08-14 **and again** 08-19, unrouted.
2. **`.gitattributes`** (`docs/whitepaper/** text eol=lf`, same for `whitepaper/build/**`) — fired
   again today at 12,012 raw vs 118 content. One-time normalization, dp's call.
3. **Preprint / two-paper strategy** (50 days) and the **TEST-09/TEST-10 count-collapse** question
   (12 days) — both still gated on dp.

---

## Summary

The window's research was two **under**-refutations — the rarer direction for this program — and
neither had reached the whitepaper, so the public surface was resting a live exclusion on a
Brans-Dicke grid that solar-system physics ruled out in 2003 and advertising a scanned ω dimension
the construction does not have. Both corrected in place; both corrected statements are *stronger*
than what they replace. Verifying the ω claim symbolically produced the identity the source had only
measured — the closure pins `ω·ε_crit² → 2/3`, which is why trajectories at ω = 4×10⁴ and 10⁶ agree
to the 1.22% and 0.24% observed. Chasing the citation then exposed a cross-surface `TEST-25`
collision **created by a maintainer's correct intra-site fix**, plus a catalog with no
solar-system rows at all: the same shape as yesterday's landing-site finding, arriving from a third
direction. A namespace is only as consistent as the widest surface that cites it, and a lane cannot
check a surface it cannot see. Count HELD at 6; Bucket 0 = 0.
