# Publisher Daily Report - 2026-08-24

**RUN-ID 18412** · window 2026-08-23T11:00Z → 2026-08-24T10:30Z

---

## Headline

The window's one substantive finding was **half right, and the wrong half had already been inscribed
in `PREDICTIONS.md`.** This pass refuted it by computation before it reached the whitepaper, corrected
it in both places, and amended an *internal contradiction* in `15-dark-matter` that had sat unmarked
on both sides since the section was written.

**EFE = 0 returns from "triply compromised" to doubly compromised.** Refutation count **HELD at 6**;
Bucket 0 = 0.

---

## Phase 0: Publication Recommendations

### New Recommendations
None. Archivist reports the **17th consecutive deliberate zero** — 0 new numbered Synchronism
sessions. Correct for a track at rest by charter; not an anomaly.

### Status Changes
- **REC-2026-038** (Repository-Mediated Continuity) — readiness **HELD 0.93**. **Instance 21 logged**,
  and it is the inverse of this lane's own standing rule. The rule reads *a phrase-grep proves a
  phrase, not a claim* — search a paraphrase, learn the paraphrase's fate. Today's instance runs it
  backwards: the arc **refuted a paraphrase and reported the claim refuted**. Detail below.
- **REC-2026-036 / 040 / 041** — unchanged (0.45 / 0.58 / 0.55). No new evidence in window.
- Two-paper strategy, three-preprint decision, DR3 registration terms, TEST-09/TEST-10 collapse
  question: all still **gated on dp**, unchanged.

### Upcoming Candidates
None new. Surfaces scanned per §1b: `SESSION_MAP.yaml`, `Research/papers/`, `Research/proposals/`,
`Research/preregistrations/`, `explorations/`, and both sibling repos **at `origin/<default>`**:

| repo | ref scanned | local HEAD | commits in window |
|---|---|---|---|
| Synchronism | `HEAD` (own repo) | `main` | 3 non-CI (all EFE=0/DF2 material) |
| synchronism-site | `origin/main` | `main` @ `09e9794` | **0** |
| web4 | `origin/main` | `main` | 1 (`fe29bf8a`, hub docs) |

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: Current — one amendment made.
- **Sections affected**: `05-quantum-macro/15-dark-matter/dark_matter.md`
- **Changes Made**: `[AMENDED 2026-08-24]` on the **Session #97 DF2/DF4 Resolution** paragraph.

**The landing site was an internal contradiction, 108 lines wide, unmarked on both sides.**
`15-dark-matter` states the constitutive law `C(ρ)`; evaluated at NGC 1052-DF2's *measured present*
density it returns **C ≈ 0.04** — a low-C, boost-25 object. The same section's Session #97 paragraph
resolves DF2 by asserting a tidally-stripped **"high-C core."** Both live, neither flagged. The
C ≈ 0.04 figure is the framework's own (`manuscripts/arXiv_preprint_draft_v1.md` §5.1).

**What I refused to propagate, and why.** The 08-23 finding concluded that the framework's two DF2
repairs *"invalidate the EFE = 0 premise,"* and `PREDICTIONS.md` now read: *"Admitting **ANY**
non-local C-contribution invalidates the `C local ⇒ EFE=0` derivation **globally**."* That premise —
*"C depends only on local ρ"* — is a **parenthetical gloss on the site's `/mond-unification` page**,
quoted at line 130 of the source finding itself. The **derivation**, stated identically in the
whitepaper and in the site's *own two earlier findings*, is:

> "It is **linear in Φ** (C depends on ρ, not ∇Φ), so superposition holds and **EFE = 0 exactly**."

Locality is sufficient for that, not necessary. Measured rather than argued — `∇·[C∇Φ] = 4πGρ` on a
161² grid, comparing the dwarf's internal relative field with and without a host at 12 scale lengths:

| C | non-local? | Φ-dependent? | measured EFE |
|---|---|---|---|
| fixed formation-memory field (no dependence on present ρ) | **yes** | no | **5.6×10⁻¹³** |
| `tanh(γ·ln(1+\|∇Φ\|/a₀))` | no | **yes** | **4.6×10⁻²** |

Eleven orders apart. `C_eff = max(C(ρ_local), C_formation)` preserves EFE = 0 **exactly**.

**What survives is bigger than what failed.** Both repairs abandon **`C = C(ρ)` itself** — the
sector's constitutive posit. §5.2's own testable prediction is stated to hold *"regardless of current
density,"* and a function of local density cannot predict something that holds regardless of local
density. Correctly found, **mis-located**. The screening half (the boost ceiling pre-empting a clean
EFE test) is untouched and stands.

**Two figure corrections and one confirmation the source did not make.** Anchoring on published
baryons (Wolf+2010; M⋆ = 2×10⁸ M☉, R_e = 2.2 kpc) gives **σ_N = 6.99 km/s**, reproducing van Dokkum
et al. 2018's published "≈ 7 km/s from the stars alone":
1. Strict `C = 0.04` predicts **σ = 35 km/s**, not the circulated **80** — a **~4× miss, not ~10×**
   (80 needs σ_N = 16, i.e. 2.3× the published baryons). Direction unchanged; the failure is real.
2. The formation-coherence repair **works**: C_formation ∈ [0.5, 0.7] ⇒ σ = 8.4–9.9 km/s vs 8.5
   observed, with the independently derived **C_req = 0.68** inside its stated band. So the two
   repairs are **not equivalent**, as the source treated them — one reproduces DF2, the other
   contradicts the density it is invoked to explain.

- **Terminology Concerns**: None. `C(ξ)`, `γ`, MRH, Intent, Entity all used canonically.
- **Build**: `make-md.sh` exit 0, `make-web.sh` exit 0. Churn measured **both ways** per the standing
  rule — **14 content lines vs 11,956 raw**; tree restored with `git checkout`, **CI builds the
  artifacts**. (`make-web.sh` also wanted to revert `index.html`/`style.css` by 378/206 lines against
  what CI committed — another reason not to stage locally built output.)
- **Script**: `simulations/efe_locality_vs_phi_dependence.py`, exit 0.

### Web4 Whitepaper
- **Status**: Current. **No proposals, no changes, no terminology drift.**
- **Repos checked**: `web4` at `origin/main` (local HEAD also `main`). 1 commit in window; **0 touching
  `whitepaper/`**.
- `fe29bf8a` reconciles **F0 R7c status downward** in `hub/docs/SPRINTS.md`: *"Phase 0 is not
  complete,"* R7c **partial**, deploy-closure protection still open at #709, R7d not started. This
  advances the armed watch item (hub Sprint F0 = first R7 implementation) **in the honest direction** —
  a sprint doc correcting its own completion claim rather than inheriting it. Hub-doc bookkeeping, not
  a protocol element; R7 is already canonical in the terminology table. No integration merited.

---

## Routed — a correction with no landing site

**The site's `/mond-unification` gloss is genuinely wrong.** It reads *"Predicts EFE = 0 structurally
(C depends only on local ρ)"*; it should read **"(C independent of Φ)"**. As far as this scan reached,
it is the **only** place in the corpus stating the premise incorrectly — the whitepaper and the site's
own findings both state it correctly. I cannot reach it: the maintainer lane has been **401 since
2026-08-12 (day 12)** and the site had **0 commits on `origin/main`** this window, though its explorer
lane ran on 08-23.

Filed at `Research/proposals/efe_zero_premise_is_phi_independence_not_locality_20260824.md`, and
opened as a watch item that is itself a test: **can a correction cross lanes when the owning lane is
down?** The explorer lane is alive and the maintainer lane is not.

---

## Summary

The window's finding correctly identified that the framework's DF2 repairs break something, and named
the wrong thing — it refuted a **gloss** of the EFE = 0 premise and reported the **derivation**
refuted, one day after the same lane's own 08-22 result had replaced the locality axis with *"which
derivative of Φ."* Computation settles it: a fully non-local Φ-independent C gives EFE = 5.6×10⁻¹³,
a ∇Φ-keyed C gives 4.6×10⁻². EFE = 0 goes back to **doubly** compromised, and the defect the finding
actually located — both repairs abandon `C = C(ρ)` — is the larger one. Corrected in the whitepaper
and `PREDICTIONS.md`, with two figure corrections and one confirmation the source missed. **Count HELD
at 6; Bucket 0 = 0.** The cheapest prior art available anywhere is a section read against itself, and
that is exactly what went unread.
