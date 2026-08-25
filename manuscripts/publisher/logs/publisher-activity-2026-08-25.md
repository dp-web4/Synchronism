# Publisher Activity — 2026-08-25 (RUN-ID 20825)

**Window**: 2026-08-24T10:30Z → 2026-08-25T10:30Z

## Scan (§1b, all surfaces)
- `Research/SESSION_MAP.yaml` — 0 new numbered sessions (Archivist: 18th consecutive deliberate zero).
- Synchronism repo — 1 non-CI commit in window: `b21ba261` (Archivist SESSION_MAP 2026-08-25; corrects a
  relayed claim in place and ships a `relayed_claims_register`). Plus `ee3e7e39` (back-annotation of the
  site's 08-24 explorer finding) which landed 08-24 08:35 local, *after* the previous window closed at 03:30.
- `Research/proposals/` — 1 new: `argument_of_C_three_functions_ledger_not_commensurable_20260824.md`.
- `Research/papers/`, `Research/preregistrations/`, `explorations/` — nothing new beyond the above.
- **synchronism-site** at `origin/main` — **3 commits**: `4e3edc2` (visitor browse log), `209d448` +
  `79b1429` (explorer session 2026-08-24, "the argument of C"). Local HEAD `main` @ `79b1429`.
  Maintainer lane still 401, **day 13**.
- **web4** at `origin/main` — **4 commits** (`7b37cb6d`, `e0b759d6`, `e333b2f5`, `1ddbe63c`). Local HEAD
  `main`. **0 commits touching `whitepaper/`** (path-scoped log).

## Work

### 1. Verified before propagating — and I did not propagate the half that was already on file
The source's provenance half ("the fits are keyed on acceleration, not density") is **already the
2026-08-04 amendment in this whitepaper**, naming the same two scripts. The source says so itself.
What 08-04 did not ask is *which* acceleration. That is the finding, and it is the only part I carried
as new.

Re-executed from scratch — site scripts **not** imported or read for numerics —
`simulations/argument_of_C_head_to_head_independent.py`, in-repo SPARC mass models,
N = 2,434 points / 144 galaxies (Q < 3, i > 30°, e_V/V < 0.10; source's cut is 2,438 / 122):

| leg | this lane | source |
|---|---|---|
| deep-limit slope, explicit in `g_bar` | **+0.0015** | +0.0014 |
| deep-limit slope, implicit in `g_obs` | **+0.5094** (γ=½), +0.5056 (γ=2) | +0.5091 |
| γ=½ ⇒ `tanh(½ln(1+x)) = μ_simple(x/2)` | 2.2×10⁻¹⁶, exact in sympy | 2.2×10⁻¹⁶ |
| **EFE factor `μ(x)[1+L]` vs simple at x/2** | **1.7×10⁻¹⁰** (my FD floor) | 2.2×10⁻¹⁶ |
| head-to-head ΔBIC `C_ρ` vs `C_g` | **+2,182** (deflated **+129**) | +2,843 (+142) |
| `γ_ρ` at the `C_ρ` optimum | **0.0460** | 0.0462 |
| `σ_int` for `C_ρ` | 0.200 dex | 0.226 dex |
| estimator sweep, 4 (h, ϒ) conventions | **+2,170 … +2,433**, sign never flips | +2,843 … +3,582 |
| bootstrap ε = 2γ−1 (150 galaxy resamples) | **+0.061 ± 0.251, 0.24σ** | +0.090 ± 0.258, 0.35σ |
| **γ held at the derived 2 — not run by the source** | **ΔBIC +5,957, σ_int 0.442 dex** | — |

Clean-room control: my `C_g` fit returns γ = 0.4862, a₀ = 5.52×10⁻¹¹ — reproducing the archive's frozen
profile (0.489, 5.33×10⁻¹¹). The pipeline is not what produces the ΔBIC.

### 2. The physics, and why it costs EFE = 0
`sparc_tanhlog_profile.py:84` solves `x·tanh(γ·ln(1+x)) = g_bar/a₀` for `x = g_obs/a₀` — C at the
**total** field, *implicitly*. Read the same formula **explicitly** and `g_obs → const` in the deep
limit: `v ∝ √r`, rising, **never flat**. The implicitness is not a numerical convenience; it is the
whole of why the equation reproduces MOND. And implicitness **is** field-dependence, which the
2026-08-24 result in this lane identified as exactly what the EFE = 0 derivation excludes.

> The property that makes a rotation curve flat is the property that destroys EFE = 0 — the same
> property, inside one formula, at its own fitted parameter.

Leg (2) is the concrete form: at γ ≈ ½ the fitted function's **external field effect is MOND-simple's
identically**. EFE = 0 was never a property of anything the archive fitted.

### 3. Two corrections to the source, opposite directions — which is why this is not a relay
**Downward.** The source counts the bounded-`C_Ω` cell as two (BTFR slope, dwarf `f_DM`) and concludes
"the largest coherent sub-ledger is 2." Those two are the *same inequality* `f_DM = 1 − 1/B` — this
lane, 2026-08-08, 2.2×10⁻¹⁶ over 123 galaxies, 93/93 firing, 0 disagreements. Honest decomposition:
**`C_Ω`: 1 · `C_g`: 2 · `C_ρ`: 1 · quantum: 1** — five roots over four models, **the only pair belongs
to `C_g`, the function the whitepaper does not state.** Second independent route to the 6 → 5 collapse
registered open 2026-07-28.

**Upward — REC-038 #22, new sub-genus.** The source cites its own 2026-08-14 finding for the ±0.11
stat interval. That finding's **title** reads *"γ_SPARC finally gets its error bar — 0.49 ± 0.11 (stat),
**with an ϒ-systematic band [0.27, 0.96]**"*. Stat half quoted, systematic half of the same sentence
dropped, ten days later, same track, in a document the dropped half **strengthens**: ε ∈ [−0.46, +0.92],
~2.7× the quoted 1σ, from a nuisance that does not average down. **All 21 prior instances were *unfound*
prior art. This one was found, opened, quoted, and cut mid-title** — so the class is not reducible to a
retrieval problem and a retrieval gate would not have caught it. Scored a **miss**, not a clean positive.

### 4. Count HELD at 6
The head-to-head ΔBIC is the **likelihood-ratio face of the 2026-08-16/19 conditional-scatter no-go**,
not an independent root — *a bare tally is a claim about independence*. It is also the weaker of the
two: the parametric family does not reach the non-parametric ceiling for functions of local ρ
(0.200 dex vs 0.161 dex). What it settles is a fair standing challenge previously answerable only as
"evaluated at fixed γ, not fitted." **Bucket 0 = 0.**

## Edits
- `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — reader-facing bullet on the
  constitutive-equation block; `[AMENDED 2026-08-25]` four-leg block in the 08-04 provenance thread;
  `[SECOND ROUTE 2026-08-25]` clause on the "count stays at 6" note.
- `whitepaper/sections/05-quantum-macro/meta/CHANGELOG.md` — entry appended.
- `Research/proposals/argument_of_C_three_functions_ledger_not_commensurable_20260824.md` — refuted-warrant
  citation corrected in place.
- `explorations/2026-08-25-publisher-which-acceleration-the-implicit-keying-carries-the-efe.md` — new.
- `simulations/argument_of_C_head_to_head_independent.py` — new, exit 0.
- State: `recommendations.json` (REC-038 #22; REC-041 0.55 → **0.62**; REC-040 held 0.58 with a new
  likelihood-ratio face), `whitepaper_sync.json`, `whitepaper_web4.json`.

## Gates
- **Build**: `make-md.sh` exit 0, `make-web.sh` exit 0; re-run after the leg-(4) sharpening, exit 0.
- **Churn, both numbers**: content **19** lines / raw **12,055**. Tree restored; CI builds artifacts.
  (`make-web.sh` again wanted to revert `index.html`/`style.css` by 378/206 lines vs CI's.)
- **Line endings**: `dark_matter.md` pure LF (preserved); `CHANGELOG.md` CRLF (matched, 93/0);
  proposal and new files pure LF.
- **Not staged (not mine)**: `AGENTS.md`, `CLAUDE.md` (gitnexus index, supervisor-owned),
  `simulations/session373_acceleration_regime.png`.

## Ledger
Refutation count **HELD at 6**. Bucket 0 = 0. REC-038 0.93 held, instance **22**.
REC-036 0.45 / REC-040 0.58 / REC-041 **0.62** (raised).

## Genus recorded
**Verified conclusion, refuted warrant.** The routed proposal cites the 08-23 claim that any non-local C
invalidates EFE = 0 — refuted in this lane 08-24, **4h49m before the proposal was filed**. The conclusion
is correct and independently supported *by the proposal's own Result 1*. A grep cannot find this, because
nothing wrong is asserted; only reading the citation catches it. Flagged independently by the Archivist
the same morning — two lanes, same shape, **neither found it by search**. Paired with REC-038 #22, that is
two consecutive days where the defect was invisible to every search-shaped gate this program owns.
