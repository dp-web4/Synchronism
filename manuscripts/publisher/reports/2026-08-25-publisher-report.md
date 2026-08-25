# Publisher Daily Report — 2026-08-25

**Window**: 2026-08-24T10:30Z → 2026-08-25T10:30Z · **RUN-ID** 20825

## Phase 0: Publication Recommendations

### New Recommendations
None. No new numbered sessions (Archivist: 18th consecutive deliberate zero).

### Status Changes
- **REC-2026-041** (regulator/γ-relabelling) — readiness **0.55 → 0.62**. Its central claim gained a
  second independently verified leg: the one-line formula has two readings, and only the *implicit*
  one (C at `g_obs`) produces a flat rotation curve. "The regulator creates the deep-limit power law"
  is true only of that reading. Blocker unchanged (external prior-art gate, Famaey & McGaugh 2012).
- **REC-2026-040** (α-admixture bound) — readiness **HELD at 0.58**, but the deliverable gained a
  *likelihood-ratio face*: the density keying fitted head-to-head loses by ΔBIC +2,182 (+129
  deflated), which is the form a referee asks for and the scatter statement was not. Sole blocker
  unchanged and unrun — the external prior-art gate, now **22 days open, still priced at one pass**.
- **REC-2026-038** (prior art) — **instance #22**, and the first of a new sub-genus (§3 below).
  Readiness held at 0.93.
- REC-2026-036 (0.45), REC-2026-037 (0.92), REC-2026-039 (0.38) — unchanged.

### Upcoming Candidates
No arc reached terminus. The `argument_of_C` material is not a candidate in its own right: it is a
correction to how existing results are attributed, and it lands in REC-040/041 rather than opening a
new one.

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: Updated.
- **Sessions reviewed**: no new numbered sessions; window carried `b21ba261` (Archivist) and
  `ee3e7e39` (back-annotation of the site's 08-24 explorer finding, landed after the last window closed).
- **Changes made**: three edits to `05-quantum-macro/15-dark-matter/dark_matter.md` — a reader-facing
  bullet on the constitutive-equation block, an `[AMENDED 2026-08-25]` four-leg block in the
  acceleration-vs-density provenance thread, and a `[SECOND ROUTE 2026-08-25]` clause on the
  "count stays at 6" note. Plus a CHANGELOG entry, a triage note, and an in-place correction to the
  routed proposal.
- **Terminology concerns**: none. `C(ξ)`, `γ`, MRH unchanged. The amendment introduces `C_ρ`/`C_g`/`C_Ω`
  as *labels for readings of one equation*, explicitly not as new canonical terms.

### Web4 Whitepaper
- **Status**: Current. Scanned at `origin/main` (local HEAD also `main`).
- **Repos checked**: web4 — 4 commits (`7b37cb6d`, `e0b759d6`, `e333b2f5`, `1ddbe63c`), **0 touching
  `whitepaper/`** by path-scoped log.
- **Proposals / changes**: none merited. #781 is hub admission behaviour; #778 and `1ddbe63c` are hub
  concept docs; R6/R7 already canonical.

---

## 1. The day's finding: which acceleration — and it decides the EFE

The site explorer's 2026-08-24 finding says the galaxy sector runs three different coherence
functions. The **provenance half is not new here** — the 2026-08-04 amendment already recorded that
the frozen instrument is keyed on acceleration rather than density, naming the same two scripts, and
the source says so itself. What 08-04 did not ask is *which* acceleration, and that is the whole
finding.

Re-executed from scratch in this lane (`simulations/argument_of_C_head_to_head_independent.py`,
in-repo SPARC mass models, N = 2,434 points / 144 galaxies; the site's scripts were **not** imported
or read for numerics):

| leg | this lane | source | verdict |
|---|---|---|---|
| deep-limit slope, explicit in `g_bar` | **+0.0015** | +0.0014 | reproduces |
| deep-limit slope, implicit in `g_obs` | **+0.5094** | +0.5091 | reproduces |
| γ=½ ⇒ `tanh(½ln(1+x)) = μ_simple(x/2)` | 2.2×10⁻¹⁶ (exact in sympy) | 2.2×10⁻¹⁶ | reproduces |
| **EFE factor `μ(x)[1+L]` vs simple at x/2** | **1.7×10⁻¹⁰** (my FD floor) | 2.2×10⁻¹⁶ | reproduces |
| head-to-head ΔBIC, `C_ρ` vs `C_g` | **+2,182** (deflated **+129**) | +2,843 (+142) | sign/order/γ_ρ reproduce |
| `γ_ρ` at the `C_ρ` optimum | **0.0460** | 0.0462 | 3 digits |
| estimator sweep, 4 (h, ϒ) conventions | **+2,170 … +2,433** | +2,843 … +3,582 | sign never flips |
| galaxy-level bootstrap, ε = 2γ−1 | **+0.061 ± 0.251 (0.24σ)** | +0.090 ± 0.258 (0.35σ) | reproduces |
| **added here**: γ held at the derived 2 | **ΔBIC +5,957, σ_int 0.442 dex** | not run | new |

My clean-room `C_g` fit independently returns γ = 0.4862, a₀ = 5.52×10⁻¹¹ — reproducing the archive's
frozen profile (γ ≈ 0.489, a₀ = 5.33×10⁻¹¹), which is the check that says the sample and likelihood
are not what produce the ΔBIC.

**The physics, once.** `sparc_tanhlog_profile.py:84` solves `x·tanh(γ·ln(1+x)) = g_bar/a₀` for
`x = g_obs/a₀` — C at the *total* field, self-consistently. Read the same formula explicitly and
`g_obs → const`: `v ∝ √r`, rising, never flat. The implicitness is not a numerical convenience, it is
the whole of why the equation reproduces MOND. **And the implicitness is field-dependence** — which
the 2026-08-24 result in this lane identified as exactly what the EFE = 0 derivation excludes. So:

> The property that makes a rotation curve flat is the property that destroys EFE = 0. Same property,
> measured on this section's own one-line formula, at its own fitted parameter.

Leg (2) makes that concrete rather than structural: at γ ≈ ½ the *external field effect* of the
fitted function is MOND-simple's identically. **EFE = 0 was never a property of anything the archive
fitted.**

## 2. Why the count did not move

The head-to-head ΔBIC is the **likelihood-ratio face of the 2026-08-16/19 conditional-scatter no-go**
already on file, not an independent empirical root — the standing rule *a bare tally is a claim about
independence* cuts against admitting it. It is also the *weaker* of the two: the parametric family
does not reach the non-parametric ceiling for functions of local ρ (σ_int 0.200 dex vs 0.161 dex
best-achievable). What it does settle is a fair standing challenge the archive could previously
answer only as *"evaluated at fixed γ, not fitted"*: the novel claim has now been fitted freely, and
it loses.

**Refutation count HELD at 6. Bucket 0 = 0.**

The **ledger recount** is where two corrections met, in opposite directions:

- **Downward.** The source's decomposition gives `C_Ω` two entries (BTFR slope, dwarf `f_DM`) and
  concludes "the largest coherent sub-ledger is 2." Those two are the *same inequality*,
  `f_DM = 1 − 1/B`, verified in this lane 2026-08-08 to 2.2×10⁻¹⁶ over 123 galaxies, 93 and 93 firing
  with zero disagreements. Honest decomposition: **`C_Ω`: 1 · `C_g`: 2 · `C_ρ`: 1 · quantum: 1** —
  five roots over four models, and **the only pair belongs to `C_g`, the function the whitepaper does
  not state.** This is the 6 → 5 collapse registered open 2026-07-28 and gated on dp, now reached a
  second time by a wholly independent route. It is no longer one lane's observation.
- **Upward.** See §3.

## 3. REC-038 #22 — the first instance where retrieval succeeded and the *reading* failed

The source concludes "SPARC does not constrain ε" from `ε = +0.090 ± 0.258`, and cites **its own
2026-08-14 finding** for the statistical interval. That finding's **title** reads:

> *"γ_SPARC finally gets its error bar — 0.49 ± 0.11 (stat), **with an ϒ-systematic band [0.27, 0.96]**"*

The stat half is quoted. The systematic half of the same sentence is dropped — ten days later, in the
same track, in a document whose own conclusion the dropped half *strengthens*. Carried through,
γ̂ ∈ [0.27, 0.96] at flat rms is **ε ∈ [−0.46, +0.92]** — ~2.7× the quoted 1σ, comparable to its 2σ,
from a nuisance that does not average down. The parameter is **unidentified at factor 2**, not merely
unconstrained.

All 21 prior instances of this class were *unfound* prior art. This one was found, opened, quoted,
and cut mid-title. **The class is therefore not reducible to a retrieval problem, and a retrieval gate
would not have caught it.** Scored as a miss, not a clean positive: the gate ran and returned the
wrong half. Detector was reading the citation, not searching for it — the second consecutive instance
a grep could not have found.

## 4. Genus recorded: verified conclusion, refuted warrant

The routed proposal cites the 2026-08-23 claim that any non-local C invalidates EFE = 0 — refuted in
this lane on 2026-08-24, **4h49m before the proposal was filed**. The *conclusion* is correct and
independently supported **by the proposal's own Result 1** (keying C to `g_obs` *is* Φ-dependence).
A grep cannot find this defect, because nothing wrong is asserted; only reading the citation catches
it. Corrected in place. Flagged independently by the Archivist the same morning — two lanes, same
shape, neither found it by search.

## 5. Gates

- **Build**: `make-md.sh` exit 0, `make-web.sh` exit 0 (re-run after the leg-(4) sharpening: exit 0).
- **Churn, both numbers**: content **19** lines / raw **12,055**. Tree restored with
  `git checkout -- docs/whitepaper/ whitepaper/build/`; CI builds the artifacts.
  (`make-web.sh` again wanted to revert `index.html`/`style.css` by 378/206 lines vs CI's.)
- **Line endings**: `dark_matter.md` pure LF (preserved, 0 CRLF / 0 lone CR); `CHANGELOG.md` CRLF
  (matched, 93 / 0); new files pure LF.
- **Not staged (not mine)**: `AGENTS.md`, `CLAUDE.md` (gitnexus index counts, supervisor-owned),
  `simulations/session373_acceleration_regime.png`.

## Summary

The site's "three coherence functions" finding is right about the trichotomy and I verified its three
executable legs independently, including the one it did not run (γ held at the framework's own derived
value 2: ΔBIC +5,957, σ_int 0.442 dex). The propagation is not a relay because it carries two
corrections in opposite directions — the ledger recount is one *smaller* than the source states, and
the ε conclusion is one *stronger*. The durable result is that EFE = 0, the sector's one distinctive
prediction, was never a property of any function the archive fitted: the self-consistency that makes a
rotation curve flat is the field-dependence that destroys it, now measured inside a single one-line
formula rather than argued across classes of coupling. Refutation count HELD at 6.
