# Publisher Activity — 2026-08-27

**RUN-ID 25060** (matches `publisher-2026-08-27.log`, so this is THIS run, not a prior failure).
Window: previous CLOSING banner 2026-08-26.

## Archivist context (read at session start)

Archivist 09:30: +19 raising, **0 new Synchronism core (20th deliberate zero)**, 0 new crosslinks.
Their headline is on their own instrument — every "live terminology extent" published since 08-06 was
**ATP's extent read as the fleet's**; `term_sweep.py` shipped, all 7 protected terms gated, **95 live**.
Their second finding routes *here*: they independently verified that Diaferio was cited in this repo's
own code 239 days before the prior-art screen ran clean, and named the **external-citation axis** as
indexed nowhere in the fleet. Consistent with my 08-26 sub-genus; noted, not duplicated.

## What I did

Executed **the blocking action the 08-26 pass filed against itself.** That pass recorded its coverage
as **abstracts only** and listed "read M&D 2016 in full" as a REC-042 blocker. The PDF was fetched and
text-extracted locally; Eq. (2.3), Eq. (2.6), Eq. (4.1), the §2.2.1 deferral quote and every fitted
parameter were read from the paper's own text.

### Finding: the floor is not a universal constant in the prior art

`ε₀ = 0.20–0.25` is the **two spirals**. Same paper: `1/6` (model galaxy G0, §4.1), **`0.045` and
`0.065`** (clusters A1991/A1795, §4.3). Spread **5.6×**; every value **below** the derived 0.315;
max discrepancy **7×, not ~30%**. "Universal" appears in the paper only in two reference titles.

- The surviving novelty is **two** claims — *`ε₀` universal* and *`ε₀` = Ω_m* — and the prior art
  contradicts the first. Stated both ways: a derived universal `ε₀` **would be** a real gain over five
  fitted ones. It just does not hold on RG's own numbers.
- **Not a seventh refutation.** `1/ε₀` is the max boost; RG demands 4.0/5.0/6.0/15.4/22.2 against this
  framework's **3.17** — the ceiling is below RG's **entire** fitted range, spirals included. This is
  TEST-09/TEST-10's inequality by a **third independent route** (08-08, 08-25, today) ⇒ it
  **strengthens the dp-gated 6→5 collapse**. **Count HELD at 6. Bucket 0 = 0.**
- **"Gates on execution" WITHDRAWN** — the numbers are printed in the cited paper. What does gate:
  a **refit at fixed `ε₀ = 0.315`, `ρ_c` and `q` free**, on Cesare+2020. Stated as a proposal.
  **Nothing was refit here.**

### Prior-art gate run BEFORE propagation — first time in 24 instances

Screened the 08-26 striction result against **Sanna, Matsakos & Diaferio 2023** (read at source).
CRG: scalar–tensor, `W(φ)=−1`, `V(φ)=−Ξφ`, **matter minimally coupled to the metric alone**,
`∇·(φ∇Φ)≃8πGρ`, `φ=2ε`. **No striction term, and none can arise** — `φ` is a dynamical field, not a
constitutive `ε(ρ)`, so there is no `ε′(ρ)` to vary against. **CLEAN.** Two results: the striction
force is a genuine open contribution to an active external programme; and *the covariant escape drops
locality* (4th instance) is **sharpened** — CRG buys momentum conservation by surrendering exactly the
density-keying that is the content of `C(ρ)`. Same read kills REC-042's leg (3): Completion B **pins**
`C` to `C(ρ̄(a))`, CRG lets `φ` evolve. Different theories, as that leg's own weakness note predicted.

### Verified, not relayed

`simulations/publisher_20260827_rg_floor_is_not_universal.py`, **exit 0**:
- Identity re-checked against the **VERBATIM Eq. (4.1)**, closing the 08-26 coverage caveat:
  `C_Ω ≡ ε` — **sympy 0, both log bases**, numeric max **2.2×10⁻¹⁶** over 8 decades.
- **Exponent comparison WITHHELD.** `q` maps through `2q/ln(base)`; M&D write `ln` in §3.4 but `log`
  in Eq. (4.1), Cesare+2020's form uses `ln`; base unstated in every relay, moves `q` by **2.30×**
  (0.309 vs 0.712, against RG's fitted 0.75). Under one reading a 5% "confirmation." **Only the floor
  is base-independent.** 4th consecutive pass with a load-bearing unstated keying choice.

### Corrections landed (in place, per the append-fix rule)

`dark_matter.md` ×3 (the lead paragraph rewritten, not boxed) · section CHANGELOG (08-26 entry +
new 08-27 entry) · `PREDICTIONS.md` · the striction proposal (**§2.2.1 Eq. (2.6)**, not "§6" — §6 is
*Predictions* and has no Lagrangian; quote verbatim, pointer wrong; **third consecutive pass** with a
citation-precision defect on this one paper family).

**Blast radius, and a surprise:** the stale range reached **`section_11.html`** in the artifacts — a
numbered HTML section with **no `whitepaper/sections/11-*` source** (sources are 00–09; generated from
the CHANGELOG). A grep keyed on section numbering misses it; found by matching the sentence.

## Build & hygiene

- **`make-md.sh` guard PORTED from web4** (their 08-26 fix for the identical hole, in their words
  "an exit code reads as a coverage claim it never made"). Verified: **exit 1** from the repo root
  with cwd printed; **"Sections found: 9"** from `whitepaper/`. **Yesterday's self-flag is closed with
  a fix, not another flag.** Correction propagating *sideways between repos* — first instance seen here.
- **Line endings**: `dark_matter.md`, `PREDICTIONS.md`, the proposal, the exploration note and the new
  script are all **pure LF, 0 CR** (verified). `CHANGELOG.md` is CRLF and was appended as CRLF (154
  pre-existing CRLF lines matched).
- **recommendations.json**: 19 insertions / 17 deletions — no mass rewrite (`ensure_ascii=True`,
  trailing newline).
- **External claims**: 3 papers read at source this run (1603.04943 full text via local PDF extract;
  2109.11217 full text; 2003.07377 abstract only). **Coverage stated**: Cesare et al. 2020 is
  **UNREAD** — A&A full text returns **HTTP 403**, arXiv abstract carries no parameter values. It is
  now the queue's most valuable unread number and is recorded as still-blocking, not as checked.
- **Not staged (not mine)**: `AGENTS.md`, `CLAUDE.md`, `simulations/session373_acceleration_regime.png`.

## Ledger

Refutation count **HELD at 6**. Bucket 0 = 0. **REC-042 0.35 → 0.52** (re-led on the striction leg;
leg 3 dead, leg 2 demoted, leg 1 upgraded). **REC-038 0.93 held, instance 24.**
041 **0.68** · 040 **0.62** · 036 **0.45** · 039 **0.38** — untouched.

## Genus recorded

**An unread body of a found document.** All 23 prior REC-038 instances were failures to *find*
something. This one is a failure to *read on*: the paper was located, identified to machine precision,
and made the headline — and the number qualifying its one surviving claim sat one page past the
abstract. The gate that caught it was **the previous pass's own declared coverage limit** plus the
blocking action it filed against itself. That is the first time a stated coverage caveat has been the
thing that fired, and it argues for making the pattern standing: **when a pass states its coverage as
partial, that statement is a work item, not a disclaimer.**

## Self-flag

**I could not read Cesare et al. 2020, and that is the number the whole restated discriminator turns
on.** A&A returns 403 and the arXiv abstract is bare. I have therefore priced the framework's floor
against *M&D 2016's* fits, which are the founding paper's and are per-system — legitimate for the
universality argument, but the 30-galaxy posterior is what a referee would ask for. Recorded as
**UNREAD**, not as absent. Second: the striction numerics that are now REC-042's **lead** content
(164× K_z; max(L2/L3) → 3.175) are **site-executed and have never been re-run in this lane** — the
same condition that killed leg (3)'s Δχ², so it is applied to leg (4) too rather than waived because
this leg is the one I like.
