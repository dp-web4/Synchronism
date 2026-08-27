# The floor is not a universal constant in the prior art — and the number that says so was one page past the abstract

**Date:** 2026-08-27
**Lane:** Publisher
**Status:** `[ACTIVE-MRH]` — gate-fired by executing the blocking action the previous pass filed against itself.

**Verdict.** The 08-26 pass identified the galaxy sector's field equation as **Refracted Gravity**
(Matsakos & Diaferio 2016, arXiv:1603.04943) and isolated one surviving novelty: *RG **fits** the
vacuum floor `ε₀`; this framework **derives** it as `Ω_m` = 0.315.* It priced the disagreement at
**~30%** and called it *"a live external discriminator gating on execution"* of Cesare+2020. That
pass also stated its own coverage honestly — **abstracts only** — and filed "read M&D 2016 in full"
as a REC-042 blocking action. **Executed this run.** The paper was fetched, text-extracted and read.
`0.20–0.25` is the **two spirals**. The same paper fits `ε₀ = 1/6` for its model galaxy and
**`0.045` and `0.065` for two clusters**. Spread **factor 5.6**; every published value **below**
0.315; largest discrepancy **7×, not ~30%**. Count **HELD at 6**, and for a stated reason.

---

## What the source says

Read from the paper's own text, not a relay:

| System | § | `q` | `ε₀` | implied max boost `1/ε₀` |
|---|---|---|---|---|
| NGC 6946 (spiral) | 4.2 | 3/4 | **0.25** | 4.00 |
| NGC 1560 (spiral) | 4.2 | 3/4 | **0.20** | 5.00 |
| model galaxy G0 | 4.1 | 1/2 | **1/6** | 6.00 |
| A1991 (cluster) | 4.3 | 2 | **0.045** | 22.22 |
| A1795 (cluster) | 4.3 | 2 | **0.065** | 15.38 |

And **M&D never claim `ε₀` is universal.** The word "universal" occurs in the paper only inside two
reference titles (Persic+ 1996, NFW 1996). They fit `ε₀` per system, and the values move by 5.6×.

## Two consequences, cutting opposite ways

**(1) The surviving novelty is two claims, and the prior art already contradicts the first.**
*"`ε₀` is universal"* and *"`ε₀` = Ω_m"* are separable. A single derived floor cannot reproduce a
factor-5.6 spread whatever its value. This is worth stating plainly in both directions: deriving a
universal `ε₀` where RG fits five different ones **would be a real parameter-economy gain** — that
is why the 08-26 pass was right to isolate it as the one thing left. It is simply that on RG's own
published numbers it does not hold, and that was legible in the same paper whose identification was
the previous day's headline.

**(2) It is not a seventh refutation. It is the bounded-boost inequality, by a third route.**
`1/ε₀` **is** the maximum boost. RG's fits demand 4.0, 5.0, 6.0, 15.4 and 22.2 where this
framework's ceiling is `1/Ω_m = 3.17` — **the ceiling sits below RG's entire fitted range, spirals
included.** That is exactly TEST-09/TEST-10's inequality (`B ≤ 1/Ω_m` binds), now sourced from
published fits rather than from SPARC. It is therefore a **third independent route** to the standing
**6→5 collapse** (with 2026-08-08's "TEST-09/TEST-10 are one inequality" and 2026-08-25's likelihood
form) — it *strengthens the dp-gated collapse recommendation* rather than adding to the tally.
**Count HELD at 6. Bucket 0 = 0.**

## What the discriminator actually is, restated

*"Gates on execution"* is **withdrawn**: the discriminating numbers are printed in the cited paper.
What genuinely gates on execution is narrower, and is stated as a **proposal, not a result** —
nothing was refit here. RG fits `(ε₀, ρ_c, q)` *jointly* per system, so the honest test of the
derived floor is a **refit at fixed `ε₀ = 0.315` with `ρ_c, q` free**, on Cesare+2020's 30 DiskMass
galaxies, because vertical dispersions are RG's own method for determining `ε`. That is a real,
scoped, externally-comparable execution and it remains the best item in the physics queue.

## Verified, not relayed

`simulations/publisher_20260827_rg_floor_is_not_universal.py` (exit 0):

- **The identity re-checked against the VERBATIM Eq. (4.1)**, closing the 08-26 coverage caveat that
  the forms "came through the routing note and were re-derived, not read in the papers":
  `C_Ω ≡ ε` — **sympy 0 in both log bases**, numeric max **2.2×10⁻¹⁶** over 8 decades.
- **Exponent comparison WITHHELD as unsafe.** `q` maps through `2q/ln(base)`. M&D write `ln` in §3.4
  but `log` in Eq. (4.1); Cesare+2020's reported form uses `ln`. The base is unstated in every relay
  and moves `q` by **2.30×** — the difference between `q = 0.309` and `q = 0.712` for the same
  framework, against RG's fitted 0.75. Under one reading that is a 5% agreement someone will call a
  confirmation; under the other it is 2.4× off. **Only the floor is base-independent**, so only the
  floor is quoted. (Fourth consecutive pass on which an unstated keying choice was load-bearing.)

## Prior-art gate on our own new result — run, and CLEAN

The 08-26 back-annotation's striction force is the lane's one outward-facing physics contribution.
It was screened before travelling, against the one place it could already exist — **Sanna, Matsakos
& Diaferio 2023** (A&A 674 A209; arXiv:2109.11217), read at source. CRG is scalar–tensor with
`W(φ) = −1`, `V(φ) = −Ξφ`, **matter minimally coupled to the metric alone**, weak-field
`∇·(φ∇Φ) ≃ 8πGρ`, `φ = 2ε`. No striction term appears — **and none can**: in CRG `φ` is an
independent dynamical field with its own equation of motion, *not* a constitutive `ε(ρ)`. There is
no `ε′(ρ)` to vary against. So (i) the striction result is **not prior art**, it is a genuine open
contribution; and (ii) it *sharpens* the standing observation that **the covariant escape drops
locality** (4th instance) — CRG buys momentum conservation precisely by surrendering the
density-keying that is the entire content of `C(ρ)`.

This is the **first time the prior-art gate ran before the claim propagated** rather than after.

## Corrections landed

- **Whitepaper `15-dark-matter/dark_matter.md`** — the "What survives the identification" paragraph
  **lead-corrected in place** (not appended), plus the two other sites carrying `0.20–0.25`.
- **`whitepaper/.../meta/CHANGELOG.md`** — the 08-26 entry corrected in place; new 08-27 entry.
- **`PREDICTIONS.md`** — the discriminator sentence restated.
- **`Research/proposals/appendix_d_variational_completion_*_20260826.md`** — the deferral is
  **§2.2.1, Eq. (2.6)**, not "§6" (§6 is *Predictions* and has no Lagrangian). Quote verbatim and
  substance unaffected; pointer wrong. Also records M&D's sign convention for upstream use, and the
  prior-art gate result above.

**Blast radius note.** The stale range also reached **`section_11.html`** in the built artifacts —
a numbered HTML section with **no corresponding `whitepaper/sections/11-*` source directory** (the
sources are 00–09). It is generated from the section CHANGELOG. A grep of `whitepaper/sections/`
for a "section 11" claim finds nothing; the source had to be found by matching the sentence. Worth
knowing before the next blast-radius check trusts section numbering.

## REC-2026-038, instance 24 — new sub-genus

Not an unfound **document**. An **unread body of a found document.** The paper was located,
identified to machine precision, celebrated as the sharpest result of the arc — and the number that
qualifies the one surviving claim sat one page past the abstract. All 23 prior instances were about
failing to *find* something; this one is about stopping at the summary of something already in hand.

The gate that caught it was **the previous pass's own declared coverage limit** ("abstracts only")
combined with the blocking action it filed against itself. That is the first time a *declared
coverage caveat* has been the thing that fired. It is cheap, it is mechanical, and it argues for
making the pattern standing: **when a pass states its coverage as partial, that statement is a work
item, not a disclaimer.**

## Self-flag closed with a fix

Yesterday's self-flag — `make-md.sh` exits 0 with `✅` after finding **0 of 10** sections — is
**fixed**, by porting the guard web4's lane added to *its* `make-md.sh` on the same day (08-26).
Verified: from the repo root it now prints the cwd and **exits 1**; from `whitepaper/` it reports
**"Sections found: 9"** and proceeds. Two stray directories (`Synchronism/build/`,
`ai-agents/docs/whitepaper/`) still await owner action — the safety preset denies `rm` outside
`/tmp` and no `hestia_appeal` tool is exposed here, so this lane complies rather than recasting.
