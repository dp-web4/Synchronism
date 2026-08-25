# Which acceleration? The implicit keying is what carries both the flat curve and the EFE

**Date**: 2026-08-25 · **Lane**: Publisher · **Status**: EXECUTED, independent re-implementation
**Source**: `synchronism-site/explorer/findings/the-argument-of-C-three-functions-each-killed-by-its-own-distinguishing-feature.md` (2026-08-24)
**Routed via**: `Research/proposals/argument_of_C_three_functions_ledger_not_commensurable_20260824.md`
**Script**: `simulations/argument_of_C_head_to_head_independent.py`
**Landed in**: `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` `[AMENDED 2026-08-25]`

---

## What the source claims, and what I checked

The site explorer's 2026-08-24 finding says the galaxy sector runs **three different coherence
functions** — `C_ρ` (the headline equation), `C_g` (every fitted number), `C_Ω` (the bounded ceiling)
— each carrying a different claim to novelty and each refuted by exactly that feature.

The 2026-08-04 amendment in this whitepaper already recorded that the frozen instrument is keyed on
**acceleration rather than density**, naming the same two scripts the source greps. The source says
so itself. So the provenance half is not new here, and I did not propagate it as new.

**What is new is the next question, which 08-04 did not ask: *which* acceleration.**

## Re-executed from scratch (site scripts not imported or read for numerics)

In-repo SPARC mass models, N = 2,434 points / 144 galaxies (Q < 3, i > 30°, e_V/V < 0.10 — my cut;
the source's is 2,438 / 122, which is why third digits differ).

| leg | result | source's | verdict |
|---|---|---|---|
| deep-limit slope, explicit in `g_bar` | **+0.0015** | +0.0014 | reproduces |
| deep-limit slope, implicit in `g_obs` | **+0.5094** (γ=½), +0.5056 (γ=2) | +0.5091 | reproduces |
| γ=½ ⇒ `tanh(½ln(1+x)) = μ_simple(x/2)` | 2.2×10⁻¹⁶ (exact in sympy) | 2.2×10⁻¹⁶ | reproduces |
| **EFE factor `μ(x)[1+L]` vs MOND-simple at x/2** | **1.7×10⁻¹⁰** (my FD floor) | 2.2×10⁻¹⁶ | reproduces |
| head-to-head ΔBIC, `C_ρ` vs `C_g` | **+2,182** (deflated +129) | +2,843 (+142) | sign, order, γ_ρ reproduce |
| `γ_ρ` at the `C_ρ` optimum | **0.0460** | 0.0462 | reproduces to 3 digits |
| `σ_int` for `C_ρ` | 0.200 dex | 0.226 dex | same conclusion (RAR total = 0.13) |
| estimator sweep, 4 (h, ϒ) conventions | +2,170 … +2,433 | +2,843 … +3,582 | **sign never flips** |
| **added here**: γ held at the derived 2 | **ΔBIC = +5,957, σ_int = 0.442 dex** | not run | new |

My `C_g` fit independently returns γ = 0.4862, a₀ = 5.52×10⁻¹¹ — reproducing the archive's frozen
profile (γ ≈ 0.489, a₀ = 5.33×10⁻¹¹) from a clean-room pipeline, which is the check that says my
sample and likelihood are not the thing producing the ΔBIC.

## The physics, stated once

`sparc_tanhlog_profile.py:84` solves `x·tanh(γ·ln(1+x)) = g_bar/a₀` for `x = g_obs/a₀`. C is
evaluated at the **total** field, self-consistently — *implicitly*. Read the same formula explicitly
(C at `g_bar/a₀`) and `g_obs → const` in the deep limit: `v ∝ √r`, rising, **never flat**. The
implicitness is not a numerical convenience; it is the whole of why the equation reproduces MOND.

And the implicitness *is* field-dependence. The 2026-08-24 result in this lane located the EFE = 0
derivation on **Φ-independence** (not on locality). So:

> **The property that makes a rotation curve flat is the property that destroys EFE = 0. They are
> the same property, measured on this section's own one-line formula, at its own fitted parameter.**

Leg (2) makes it concrete rather than structural: at γ ≈ ½ the *external field effect* of the fitted
function is MOND-simple's identically. EFE = 0 was never a property of anything the archive fitted.

## Two corrections to the source, in opposite directions

**Downward — the ledger recount is one too high in one cell.** The source assigns the six
refutations to the three functions and reports "the largest coherent sub-ledger is 2," giving `C_Ω`
two entries (BTFR slope, dwarf `f_DM`). Those two are the **same inequality**: `f_DM = 1 − 1/B`,
verified in this lane 2026-08-08 to 2.2×10⁻¹⁶ over 123 galaxies, the two criteria selecting 93 and
93 with zero disagreements. The honest decomposition is

> `C_Ω`: **1** · `C_g`: **2** · `C_ρ`: **1** · quantum: **1**

— five roots over four models, no model carrying more than two, and **the only pair belongs to
`C_g`, the function the whitepaper does not state.** This is the 6 → 5 collapse registered open on
2026-07-28 and gated on dp, now reached a second time by an independent route. Count **HELD at 6**;
merging is a naming convention and dp's to set. Bucket 0 = 0.

**Upward — the source under-states its own conclusion.** It reports `ε = 2γ − 1 = +0.090 ± 0.258`
(galaxy-level bootstrap) and concludes "SPARC does not constrain ε." That is the *statistical*
interval at `ϒ_disk = 0.5`. The whitepaper already carries the systematic (2026-08-15): refitting
across `ϒ_disk ∈ [0.4, 0.6]` moves γ̂ over **[0.27, 0.96]**, i.e. `ε ∈ [−0.46, +0.92]` — ~2.7× the
quoted 1σ, comparable to its 2σ, from a nuisance that does not average down. The correct statement
is stronger: ε is not merely unconstrained; the parameter is **unidentified at factor 2**.

**REC-2026-038 instance #22, and a new sub-genus: the prior art was found, cited, and truncated at
the clause boundary.** The source cites its own 2026-08-14 finding for the ±0.11 stat interval. The
*title* of that finding reads: *"γ_SPARC finally gets its error bar — 0.49 ± 0.11 (stat), **with an
ϒ-systematic band [0.27, 0.96]**"*. The stat half is cited; the systematic half of the same sentence
is dropped — ten days later, in the same track, in a document whose conclusion the dropped half
strengthens. Every previous instance of this class was *unfound* prior art. This one was found,
opened, quoted, and cut mid-title. Sourcing is 6/16 external; this is the first instance where the
retrieval succeeded and the *reading* failed.

My own bootstrap reproduces the statistical half independently: 150 galaxy-level resamples give
γ = 0.5305 ± 0.1256, **ε = +0.061 ± 0.251, 0.24σ from zero** (source: +0.090 ± 0.258, 0.35σ). So the
number the source quotes is right. It is the interval beside it that was needed.

## Why the count did not move

The head-to-head ΔBIC is the **likelihood-ratio face of the 2026-08-16/19 conditional-scatter no-go**
already on file, not an independent empirical root — my standing rule *a bare tally is a claim about
independence* cuts against admitting it. And it is the *weaker* of the two: the parametric family
does not even reach the non-parametric ceiling for functions of local ρ (σ_int 0.200 dex against the
0.161 dex best-achievable conditional scatter). What it does settle is a fair standing challenge the
archive could previously answer only as *"evaluated at fixed γ, not fitted"*: **the novel claim has
now been fitted freely, and it loses.**

## Genus recorded

**Verified conclusion, refuted warrant.** The routed proposal cites the 2026-08-23 claim that any
non-local C invalidates EFE = 0 — refuted in this lane on 2026-08-24, 4h49m before the proposal was
filed. The *conclusion* the citation supports is correct and independently supported **by the
proposal's own Result 1** (keying C to `g_obs` is Φ-dependence). A grep cannot find this defect,
because nothing wrong is being asserted; only reading the citation catches it. Corrected in place in
the proposal. Flagged independently by the Archivist the same morning — two lanes, same shape, and
neither found it by search.
