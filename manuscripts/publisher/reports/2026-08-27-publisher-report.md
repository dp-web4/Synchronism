# Publisher Daily Report — 2026-08-27

**Headline:** the previous pass declared its own coverage as *abstracts only* and filed the full
read as a blocking action. Executing that one item corrected the day-old headline's surviving claim,
killed one leg of the live recommendation, promoted a new one, and closed the previous pass's
self-flag with a fix ported from a sibling repo. **Refutation count HELD at 6. Bucket 0 = 0.**

---

## Phase 0: Publication Recommendations

### Surfaces scanned (§1b list, all of them)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new core sessions — **21st consecutive deliberate zero** (archivist confirms) |
| `Research/papers/` | unchanged since 07-23 |
| `Research/proposals/` | 1 new (08-26, the striction/variational completion) — acted on |
| `Research/preregistrations/` | unchanged |
| `explorations/` | 2 new (08-26) — both read |
| `synchronism-site/explorer/findings/` | 1 new (08-26, `l2-is-not-l3-for-a-disc-and-the-action-adds-a-force-the-tests-omit`) |
| `synchronism-site/explorer/topics/` | no new seeds |
| `web4` (`origin/main`) | 4 commits, 1 touching `whitepaper/` |

**Refs stated, per the 2026-08-09 rule:** both sibling repos scanned at `origin/main`; `HEAD` was
also `main` in both, so the ref and the checkout agreed this window.

### Status changes

| REC | Was | Now | Why |
|---|---|---|---|
| **042** | 0.35 | **0.52** | Both "read the paper in full" blockers executed. Leg (3) resolved DEAD, leg (2) demoted to a correction, leg (1) upgraded to source-verified, **new leg (4) promoted to lead**. Remaining blockers are reproduce-it-here items, not literature items. |
| **038** | 0.93 | **0.93** (held) | **Instance 24**, new sub-genus — see below |
| 041 / 040 / 036 / 039 | 0.68 / 0.62 / 0.45 / 0.38 | unchanged | not touched this window |

### The day's finding

The 08-26 pass identified the galaxy sector's field equation as **Refracted Gravity** (Matsakos &
Diaferio 2016) and isolated one surviving novelty: *RG **fits** the vacuum floor `ε₀`; this framework
**derives** it as `Ω_m` = 0.315*, priced at **~30% apart** and called *"a live external discriminator
gating on execution"* of Cesare+2020.

The paper has now been read. **`0.20–0.25` is the two spirals.** From M&D's own text:

| System | § | `q` | `ε₀` | max boost `1/ε₀` |
|---|---|---|---|---|
| NGC 6946 | 4.2 | 3/4 | 0.25 | 4.00 |
| NGC 1560 | 4.2 | 3/4 | 0.20 | 5.00 |
| model G0 | 4.1 | 1/2 | 1/6 | 6.00 |
| **A1991** | **4.3** | **2** | **0.045** | **22.22** |
| **A1795** | **4.3** | **2** | **0.065** | **15.38** |

Factor **5.6** spread; every published value **below** 0.315; max discrepancy **7×, not ~30%**. And
M&D never claim `ε₀` is universal — the word appears in the paper only in two reference titles.

**Two consequences, cutting opposite ways.** (1) The surviving novelty is **two** claims — *`ε₀` is
universal* and *`ε₀` = Ω_m* — and the prior art already contradicts the first. Deriving a universal
`ε₀` where RG fits five different ones **would be** a real parameter-economy gain; it simply does not
hold on RG's own numbers. (2) It is **not a seventh refutation**: `1/ε₀` is the max boost, RG demands
4.0–22.2 against this framework's 3.17, so the ceiling sits below RG's **entire** fitted range. That
is TEST-09/TEST-10's own inequality by a **third independent route** (with 08-08 and 08-25) — it
**strengthens the dp-gated 6→5 collapse**, it does not add to the tally.

**"Gates on execution" is withdrawn** — the discriminating numbers are printed in the cited paper.
What genuinely gates on execution is narrower and is stated as a **proposal, not a result**: RG fits
`(ε₀, ρ_c, q)` jointly, so the honest test is a **refit at fixed `ε₀ = 0.315` with `ρ_c, q` free** on
Cesare+2020's 30 DiskMass galaxies. **Nothing was refit here.**

### Prior-art gate run BEFORE the claim propagated — first time in 24 instances

The 08-26 striction result is this lane's one outward-facing physics contribution. It was screened
against the only place it could already exist — **Sanna, Matsakos & Diaferio 2023** (A&A 674 A209),
read at source. CRG is scalar–tensor, `W(φ) = −1`, `V(φ) = −Ξφ`, **matter minimally coupled to the
metric alone**, weak field `∇·(φ∇Φ) ≃ 8πGρ` with `φ = 2ε`. **No striction term appears, and none
can** — `φ` is an independent dynamical field, not a constitutive `ε(ρ)`, so there is no `ε′(ρ)` to
vary against. **Result: CLEAN.** The striction force is a genuine open contribution, and it *sharpens*
the standing "the covariant escape drops locality" observation (4th instance): CRG buys momentum
conservation precisely by surrendering the density-keying that is the whole content of `C(ρ)`.

Same read settles REC-042's **lead blocking action**: Completion B pins `C` quasi-statically to
`C(ρ̄(a))`; CRG lets `φ` evolve dynamically. **Different theories ⇒ leg (3) evaporates**, exactly as
that leg's own weakness note predicted.

### REC-2026-038 — instance 24, new sub-genus

Not an unfound **document**. An **unread body of a found document.** The paper was located, identified
to machine precision and celebrated as the arc's sharpest result; the number qualifying the one
surviving claim sat one page past the abstract. All 23 prior instances were failures to *find*.

The gate that caught it was **the previous pass's own declared coverage limit**, plus the blocking
action it filed against itself. That is the first time a stated coverage caveat has been the thing
that fired. **Proposed standing rule: when a pass states its coverage as partial, that statement is a
work item, not a disclaimer.**

---

## Phase 1: Whitepaper Review

### Synchronism — **Changes made**

| File | Change |
|---|---|
| `sections/05-quantum-macro/15-dark-matter/dark_matter.md` | "What survives the identification" paragraph **lead-corrected in place** (not appended); two further in-place range corrections |
| `sections/05-quantum-macro/meta/CHANGELOG.md` | 08-26 entry corrected in place; new 08-27 entry |
| `PREDICTIONS.md` | discriminator sentence restated |
| `Research/proposals/appendix_d_variational_completion_*_20260826.md` | deferral is **§2.2.1, Eq. (2.6)**, not "§6"; M&D's sign convention recorded; prior-art gate result inscribed |
| `simulations/publisher_20260827_rg_floor_is_not_universal.py` | new, exit 0 |
| `whitepaper/make-md.sh` | **coverage guard ported from web4** |
| `explorations/2026-08-27-...` | triage note |

**Terminology concerns:** none. `C(ξ)`, `γ`, `MRH`, `Intent`, `Entity` all used canonically.

**Verification.** `simulations/publisher_20260827_rg_floor_is_not_universal.py`, exit 0. The identity
`C_Ω ≡ ε` re-checked against the **verbatim** Eq. (4.1) rather than the relayed form — **sympy 0 in
both log bases**, numeric max **2.2×10⁻¹⁶** over 8 decades. This closes the 08-26 coverage caveat that
the forms "came through the routing note and were re-derived, not read in the papers."

**Withheld as unsafe.** The exponent comparison. `q` maps through `2q/ln(base)`; M&D write `ln` in
§3.4 but `log` in Eq. (4.1), Cesare+2020's reported form uses `ln`, and the base — unstated in every
relay — moves `q` by **2.30×**, between `q = 0.309` and `q = 0.712` against RG's fitted 0.75. Under
one reading that is a 5% agreement someone would report as a confirmation. **Only the floor is
base-independent, so only the floor is quoted.** Fourth consecutive pass on which an unstated keying
choice was load-bearing.

**Blast radius.** The stale range also reached **`section_11.html`** in the built artifacts — a
numbered HTML section with **no `whitepaper/sections/11-*` source directory** (sources are 00–09; it
is generated from the section CHANGELOG). A blast-radius grep keyed on section numbering misses it
entirely; the source had to be found by matching the sentence.

### Web4 — **Current, no proposals**

Scanned `origin/main` (HEAD also `main`). The lane touched `whitepaper/` once this window
(`e9aee591`): a stale `hestia#49` pointer superseded by the open readiness index `hestia#224`, plus a
direct link to the bypass catalogue, in `sections/11-standard-and-implementations`. Correct
self-correction, right shape, nothing for this lane to propose.

**One item travelled the other way.** web4 added a sections-coverage guard to *its* `make-md.sh` on
08-26 — *"without this the script … prints its success banner and exits 0, so an exit code reads as a
coverage claim it never made"* — which is the exact hole this lane self-flagged in **Synchronism's**
`make-md.sh` the same day. **Ported and verified**: exits 1 from the repo root with the cwd printed;
reports `Sections found: 9` from `whitepaper/`. Yesterday's self-flag is closed with a fix rather than
another flag.

---

## Summary

Reading one paper past its abstract corrected a day-old headline, resolved the live recommendation's
lead blocker in the negative, promoted a better lead in its place, and produced the 24th instance of
this programme's own prior-art failure class — in a new sub-genus that is about stopping at a summary
rather than failing to search. The refutation count held at 6, and the reason it held is now stated
rather than asserted: the new numbers are the old inequality arriving by a third route, which
strengthens the pending 6→5 collapse. The one thing this lane could not do is read Cesare et al.
2020, which is behind an HTTP 403 and is now the queue's most valuable unread number.
