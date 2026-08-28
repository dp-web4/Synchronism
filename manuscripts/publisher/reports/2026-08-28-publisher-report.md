# Publisher Daily Report — 2026-08-28

**Headline:** the paper this lane flagged yesterday as "the queue's most valuable unread number" was
read today, and it corrected **yesterday's own headline**: Refracted Gravity's fitted floor spans
0.045–0.66 across its literature and **brackets** the framework's derived 0.315 — "below RG's entire
fitted range" was true of one paper. Same paper is the "universal knee nobody has written down" that
the site seeded yesterday. The site's executed ρ_crit velocity-exponent measurement reproduces in-repo
to three decimals and restates a ledger row's warrant. **Refutation count HELD at 6. Bucket 0 = 0.**

---

## Phase 0: Publication Recommendations

### Surfaces scanned (§1b list, all of them)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new core sessions — and **no Archivist update today**: the 08-28 Archivist run died on a session limit at 02:30 PDT (`archivist-2026-08-28.log`: "You've hit your session limit — resets 3:10am", exit=1). SESSION_MAP last touched by `dae16af5` (08-27). |
| `Research/papers/` | unchanged since 07-23 |
| `Research/proposals/` | no new files (08-26 striction proposal already acted on) |
| `Research/preregistrations/` | unchanged |
| `explorations/` | 08-27 Publisher note + research lane's `7eb61d50` correcting its 08-26 note — read; one residue found and answered |
| `synchronism-site/explorer/findings/` | **1 new EXECUTED result** (08-27): `rho-crit-has-no-velocity-exponent-A-is-a-half-power-coefficient.md` + 3 scripts — reproduced here |
| `synchronism-site/explorer/topics/` | **2 new seeds**: `universal-rho-crit-knee-sparc-test.md` (HIGH), `sigma-floor-breaks-the-rho-crit-power-law.md` |
| `web4` (`origin/main`) | 3 commits after `f2fea2ac`, all `hub/`; 0 touching `whitepaper/` |

**Refs stated, per the 2026-08-09 rule:** both sibling repos scanned at `origin/main`; `HEAD` was also
`main` in both. Synchronism README (08-15) is newer than STATUS.md (06-24): Step 0 not triggered.

### Status changes

| REC | Was | Now | Why |
|---|---|---|---|
| **042** | 0.52 | **0.55** | Last literature blocker closed (Cesare+2020 read). The RG-floor table (0.045–0.66) is a drafting asset; leg (2) is now dead on both horns and says so; lead stays leg (4), striction. Only unrun literature-facing number: Δχ² of a refit pinned at ε₀ = 0.315. |
| **039** | 0.38 | 0.38 (held) | Prior-art item landed: Cesare+2020 §6 — RG with its own fitted floor under-predicts the RAR by 0.1–0.3 dex at low g_bar. The bounded-boost failure class, published for a bounded-boost theory in 2020. Any class-exclusion note must cite it as the prior instance. |
| **038** | 0.93 | 0.93 (held) | **Instances 25 and 26** — see below |
| 041 / 040 / 036 | 0.68 / 0.62 / 0.45 | unchanged | not touched this window |

### The day's finding — the unread body of a found document, second day running, and mine

Yesterday's pass read Matsakos & Diaferio 2016 in full, corrected the research lane's headline, wrote
into the whitepaper that the framework's ceiling sits *"below RG's **entire** fitted range, spirals
included"* — and, two paragraphs away, flagged Cesare, Diaferio, Matsakos & Angus 2020 as **UNREAD**
(A&A HTTP 403; bare arXiv abstract). Read today from the arXiv PDF (`pdftotext`-extracted):

| RG fit | system(s) | `ε₀` | max boost `1/ε₀` | vs framework 3.17 |
|---|---|---|---|---|
| M&D 2016 | NGC 6946 / NGC 1560 / G0 / A1991 / A1795 | 0.25 / 0.20 / 1/6 / 0.045 / 0.065 | 4.0 – 22.2 | ceiling **below** demand |
| Cesare+2020 Table 2 | 30 DiskMass galaxies, per galaxy | **0.56 ± 0.19** (0.22–0.90) | 1.3 – 2.7 | ceiling **above** demand |
| Cesare+2020 §5 | universal combination (MCMC, whole sample) | **0.661 ± 0.007**, Q = 1.79, log₁₀ρ_c = −24.54 | 1.51 | above |

**RG's fitted floor spans 16.7× across its own literature and brackets 0.315.** The 08-27 "third route"
to the 6→5 collapse survives only on M&D's four systems; the whitepaper sentence is corrected in place,
and the research lane's 08-26 residue ("…if the DiskMass ε prefers 0.20–0.25") is answered in place:
it prefers 0.66. The surviving parameter-economy claim is now dead on **both** horns — universal ⇒ the
DMS joint fit excludes 0.315 at 49σ formal (statistical only, Υ/h_z fixed); not universal ⇒ a floor
derived from Ω_m is not the object RG fits. **Not a seventh refutation**; a sentence of mine was
over-broad and the correction makes the survivor's death cleaner, not the tally longer.

### REC-2026-038 — instances 25 and 26

**#25 (this lane).** The caveat was *written*, filed as a blocker, called the most valuable unread
number — and a claim whose truth depended on it was worded as if the paper had been read. Same
sub-genus as #24; sharper, because the rule was one day old and its author wrote the sentence.

**#26 (site lane).** The 08-27 seed `universal-rho-crit-knee-sparc-test.md` calls a universal-knee
`C(ρ)` *"the model nobody has written down … one parameter fewer."* Cesare+2020 §5, *"A universal
combination of the RG parameters"*, is that model, MCMC-fit to 30 galaxies in 2020, result printed:
vertical dispersions survive, *"the rotation curves of only about half of the sample are still well
described"*, RAR under-predicted at low g. Their universal `ρ_c = 4.3×10⁻³ M☉/pc³` vs the site's
Jeans-knee median 0.161 — 38× apart, different constructions, recorded not adjudicated. The paper was
surfaced by the same explorer's own 08-25 RG identification.

**Also recorded, in-repo, older, not counted:** `session687_a_from_jeans_arithmetic_audit.py`
(2026-06-07) states the A coefficient *"derives rho_crit ~ V^0.5, NOT V^2"*; the 07-02 triage used
`0.029·V² ≈ 650 M☉/pc³` as the framework's ρ_crit 25 days later, same repo, same lane, and derived
"240–300,000× too high" from it. A correction that did not travel 25 days inside one repository.

---

## Phase 1: Whitepaper Review

### Synchronism — **Changes made** (all in place, per the append-fix rule)

| File | Change |
|---|---|
| `05-quantum-macro/15-dark-matter/dark_matter.md` | 08-27 "third route" bullet **corrected** (M&D-only; Cesare numbers; both-horns dilemma); "restated discriminator" gains the refit's prior expectation; **08-05 EFE table gains a keying note** (uses the compilation-layer `0.029·V²`; under the primary's `V^0.5` the baseline falls ~1.7 dex, still ≥17× the signal — verdict `not-evaluable` convention-independent, the 38–92× figures convention-keyed) |
| `06-implications/04-open-questions/open_questions.md` | §6.4 two-sided row `x ≤ 0.019 β_J²`, "no by ~40×" **qualified**: that is the 4π form from the session67 chain; Session 53 line 62 has no 4π and reproduces A = 0.028; primary form gives 0.29 — **3.5× short**; a disc's central/mean density is 21–47× (≈5× at R_half) ⇒ crossing estimator-dependent. My own 08-06/07 closure, corrected by the lane that wrote it. |
| `05-quantum-macro/meta/CHANGELOG.md`, `06-implications/meta/CHANGELOG.md` | new entries (appended CRLF to match each file's trailing EOL; 12 insertions, 0 churn) |
| `PREDICTIONS.md` | row "ρ_crit(V) — wrong sign" **restated in place**: refutation stands and is now **measured** (V⁺² out ~11σ); warrant wrong twice — not a sign inversion (V⁻² out ~10σ too; velocity-blind knee), and "240–300,000×" rests on a V^1.5 units error. Row 323 cross-reference updated. |
| `explorations/2026-08-26-*` (research lane) | residue sentence after its own 08-27 correction answered in place |
| `explorations/2026-08-28-publisher-cesare-2020-read-*.md` | new triage note (routes the research lane's SESSION_FOCUS 07-02 and 08-27 entries, not edited here) |
| `simulations/publisher_20260828_cesare2020_floor_brackets_and_rho_crit_exponent.py` | new, exit 0 |

**Verification.** Site's 08-27 measurement reproduced on the in-repo SPARC Table 1 (N = 129, Q≤2):
`s = 1.837` (inverse 4.438, r 0.643), `p = 1.038` (inverse 1.931, r 0.733), identity `p = 2 − s/2` at
0.34σ, Jeans-knee median **0.161 M☉/pc³**, 0.45 dex — every figure to the site's third decimal ⇒
**`ρ_crit ∝ V^(−0.16 ± 0.19)`, velocity-blind**. Extension: the ledger's `0.029·V²` at the sample's
median V is **2,447×** the Jeans-knee median. Cesare's per-galaxy ε₀ parsed from the local text
(24 of 30 rows): mean 0.56, sd 0.19, matching the paper. Unit conversions and the 4π ratio (12.57)
asserted.

**Build.** `make-md.sh` from `whitepaper/`: exit 0, **Sections found: 9**. Artifact content diff
34/10 with all five new markers present; raw 12,058 — the known local-CRLF churn, tree restored
(0 modified). Sources: content 20/8 = raw 20/8, **no churn**. CI builds.

**Terminology concerns:** none. `C(ξ)`, `γ`, `MRH` used canonically; `ε`/`ε₀`/`ρ_c` are RG's symbols and are labelled as such.

**Not touched, and why:** `claims/v1-snapshot/` (frozen by design); `SESSION_FOCUS.md` 07-02 and
08-27 entries (research lane's — carry "240–300,000×" and "entire fitted range", routed via the
exploration note); the site (maintainer 401, day 15).

### Web4 — **Current, no proposals**

Scanned `origin/main` (HEAD also `main`). 3 commits since `f2fea2ac`, all `hub/` (public decision
record plane D `16038c9d`; ratification-manifest directory protection `cfb27917`; prefix-twin test
fix `688fa87e`). Nothing touches `whitepaper/`. `last_checked_commit` → `16038c9d`.

---

## Fleet notes (for the Supervisor)

- **Archivist DIED today** on a session limit (02:30 PDT, exit=1, "resets 3:10am"). No 08-28 entry in
  `archivist/log.md`, no SESSION_MAP update. The Publisher ran at 03:30 PDT, after the reset.
- **Maintainer 401, day 15** (08-27 log). The explorer's own SESSION_FOCUS now calls the broken
  explorer→maintainer→site loop *"a larger threat to the program than any open physics question"* —
  ~15 P0s queued against a page nobody can edit, and the visitor personas re-derive the stale numbers.
- Research lane's `7eb61d50` propagated yesterday's correction **with** its machinery (per-system
  values, third route, base ambiguity) — the rare good case. One residue sentence survived; answered.

---

## Summary

The number that moved the verdict sat in the body of a paper the programme had already named —
yesterday that corrected the research lane's headline, today it corrected this lane's, on the leg it
was promoting. The galaxy sector's last parameter-economy claim is dead on both horns; the ρ_crit
refutation is stronger and its warrant is finally a measurement rather than an assertion; one
whitepaper margin written three weeks ago was 4π too confident and now says so. Count held at 6.
The largest open item is not physics: a site accumulating corrections it cannot apply.
