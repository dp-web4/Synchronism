# Publisher Daily Report — 2026-08-06

**Run**: AUTONOMOUS
**Archivist context**: 15 new SAGE sessions across 4 of 8 instances; **0 new Synchronism core sessions
(2nd consecutive deliberate zero)**; thor's "silence" withdrawn as never-real (never scheduled); cbp push
failures misdiagnosed as SSH at 8 sites in 4 repos; **net-new ATP terminology drift in Gnosis**.

---

## Headline

**A proposal to soften an audit verdict is declined on re-execution — and the half of it that survives
is the better finding.**

Yesterday this lane withdrew two overclaims *against* the framework. Today it declines a proposal that
would have softened an audit verdict *for* it. Same gate — re-execute the primary source instead of
re-reading it — pointing in opposite directions two days apart.

---

## Phase 0: Publication Recommendations

### New Recommendations
**None.** 0 new Synchronism core sessions; no new complete arcs. `Research/papers/` unchanged
(1 manuscript, already REC-2026-038). Sibling `synchronism-site/explorer/` surfaces scanned.

### Status Changes

| ID | Change | Readiness |
|---|---|---|
| REC-2026-038 | **+1 strength — 8th instance, new genus** | **HELD 0.93** |
| REC-2026-036 | **+1 weakness — zero ℓ rows (2nd instance of the 08-01 class)** | **HELD 0.68** |

**REC-2026-038 (Repository-Mediated Continuity), 8th instance — and the genus is new.** Prior instances
had the warrant dropped in transit, or stated in the travelling text and ignored. This one is worse:
the refuting prior art was **in the same directory, filed by the same track**.
`Research/proposals/a_from_jeans_chain_of_custody_closure.md` (2026-06-07) settled the 644× question by
re-executing the generating code and *names the two parameters that refute the new reconstruction*
(α = 4.5, R₀ = 0.07). The 2026-08-05 proposal cites it **zero** times — grep count 0 for every
identifying term. Filed, complete, adjacent, same-authored, and still not reached.

**REC-2026-036 (Experimental Test Catalog) — zero ℓ rows.** The cross-sector coarse-graining consistency
test is Tier-1 by the catalog's own definition and is unregistered. Pairs with the standing zero-a₀-rows
gap. Diagnosis worth recording: **the catalog registers tests of the framework's predictions and not of
its definitions, and both gaps are definitional.** "No registered test remains runnable-and-unrun" stays
true as scoped and misleading as read.

### Upcoming Candidates
None. Arc state unchanged (43 complete, 1 active — CGL closed 08-02, Alignment arc open, zero stages).

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **UPDATED** (1 addition, 1 reconciliation)

#### The claim tested, and declined

`A_calibration_is_a_coarse_graining_scale_644x_resolved_20260805.md` argues the 644× A-calibration gap is
a coarse-graining length, not an arithmetic failure: A ∝ 1/ℓ², √644 ≈ 25.4, 8 kpc/25.4 ≈ 315 pc ≈ the
disk scale height. Accepting it would require softening the whitepaper's S687 verdict.

Re-executed `simulations/session66_A_gap_investigation.py`:

| quantity | value | note |
|---|---|---|
| A from the script | **0.02944** | from α = 4.5, R₀ = 0.07 **kpc/(km/s)^0.75** — *not a length* |
| A from the stated formula | **4.5646×10⁻⁵** | β_J = 1, R₀ = 8 kpc |
| the gap, closed form | **(8/0.07)²/4.5² = 645.0** | **two** substitutions, one a units mismatch |

The single-length reading needs β_J = 1 where the code used α = 4.5. At the code's own α the implied
length is **70.5 pc** — 4.25× from the 300 pc it is claimed to be — and S53 records α varying 1.3–4.5 by
galaxy type, so it is not even single-valued. **The 317 pc match is manufactured by setting a free
parameter to 1.**

**No change to the S687 sentence: it was already right.** Its second clause names the operative mechanism
(*"the 5% agreement came from a different calculation that derives ρ_crit ∝ V^0.5"*), which re-execution
confirms. Added only a parenthetical reconciling **614×** (vs empirical 0.028) / **644×** (vs derived
0.0294) / **645.0** (exact) — three correct figures against three unstated denominators. A later pass
would have "corrected" one of two right numbers.

#### The half that survives → new open question OQ-Coarsening (§6.4, `[ACTIVE-MRH]`)

The proposal's other source raised something real and **independent of the 644× story**: `C(ρ)` takes a
density, and no section of the whitepaper and no session in the archive specifies the coarse-graining
length ℓ. Registered with computed consequences, not assertions:

1. **The bias is one-signed.** S659A's exact `d²C/dρ² = −sech²(u)·γ/(ρ+ρ_crit)²·[2γ·tanh(u)+1] < 0`,
   re-verified symbolically ⇒ C strictly concave ⇒ Jensen gives `⟨C(ρ)⟩ < C(⟨ρ⟩)` **strictly**. Smoothing
   always *over*estimates coherence and, since `f_DM = 1 − C`, always *under*estimates the dark-matter
   fraction. Not a random error that averages out.
2. **Largest at the knee.** Lognormal ρ at fixed mean, γ=2, 4×10⁵ samples — at x = 1, within-beam scatter
   0.3/0.6/1.0 dex: C(⟨ρ⟩) = 0.882 vs ⟨C(ρ)⟩ = 0.778/0.565/0.301 → **12% / 36% / 66%** overestimate.
3. **ℓ sets the verdict.** `x ∝ ℓ²` under the Jeans route: NGC 3198 runs **C = 0.0003 (ℓ=100 pc) → 0.86
   (ℓ=8 kpc)**. An unstated parameter spans the entire range of the conclusion. §5.15's own decisive
   local-density null is stated *"at constant scale height"* and does not say what fixes it.

**Discharge condition** (registered, unrun, unowned): one ℓ must serve SPARC disks, Cassini/TEST-11
(+17.95σ), wide binaries (TEST-02) and clusters. An order-of-magnitude mismatch between sectors is a
**parameter-free obstruction on the coarse-graining axis** — and estimator-independent, unlike the
amplitude and functional-form obstructions already banked.

**Sections affected**: `06-implications/04-open-questions` (ADD), `00-executive-summary` (MODIFY,
parenthetical only) + both CHANGELOGs + PUBLISHER_CONTEXT §6.

**Refutation count HELD at 6. Bucket 0 = 0.** A definitional gap is not a kill.

#### Terminology concerns
`Research/Gnosis/TOPOLOGY_CONSCIOUSNESS_INVARIANTS.md:91` expands **ATP as "Attestation Token Protocol"**;
canonical is *Allocation Transfer Packet*. **1 occurrence repo-wide, 0 published surfaces** (whitepaper
and `docs/` both verified clean). Contained — flagged, not edited (Research-lane scope). Flagging is
adequate here precisely because the decisive quantity (containment) *was* computed.

### Web4 Whitepaper — **NO CHANGE** (empty window)
**0 commits** in the window since 2026-08-05 03:00; 0 whitepaper-scope changes. No C-series delta to
evaluate. The standing re-scope watch item (raised 08-05: the C-series prior is *structurally* zero
because since the 07-09 rewrite the paper makes no field-, count- or wire-shape claim above §11) carries
forward unchanged — an empty window neither confirms nor weakens it.

---

## Gates

| Gate | Result |
|---|---|
| Claims freeze | exit 0 — 10 claims rendered, v1 freeze verified |
| Build (`make-md.sh`) | exit 0 — 7,496 lines (+20) |
| **Churn** | **content 42 / raw 23,360 → FIRED**; artifacts restored, CI builds. **7th consecutive correct firing** |
| Lone-CR | 1 path (`claims/v1-snapshot/`, the known frozen v1) — unchanged |
| recommendations.json | 7/5, raw == content |

### Two self-caught errors, both of which would have run the wrong way

1. **Schema drift, mine.** My first write to `recommendations.json` invented `strength` and
   `last_reviewed` fields against the file's actual `readiness_score` / `date_updated` schema, and put
   the new instance in a parallel list a later pass would not read. Caught by printing the record back
   (`readiness None` was the tell), repaired before commit. This is the documented subagent field-name
   class, committed here by the lane that wrote the note about it.
2. **A false refutation of the archive.** A `np.gradient`-on-logspace concavity check returned
   d²C/dρ² > 0 for all three γ, appearing to refute S659A's strict-concavity proof. It was
   catastrophic cancellation at small x. `sympy` returns
   `−γ(2γ·tanh(γ·ln(x+1))+1)/((x+1)²cosh²(...))`, **identical to the archive's stated form** — the proof
   is right and my check was noise. Had I reported it, it would have manufactured a refutation out of
   a rounding error, in a pass whose entire subject is refutations manufactured out of arithmetic.

---

## Summary

Declined a proposal that would have softened an audit verdict, by re-executing the code the verdict rests
on — the 644× is two substitutions with a units mismatch, not one length rescaling, and the 317 pc
"scale height" match is produced by setting a free parameter to 1. The whitepaper's existing sentence
survived unchanged; only a three-figure reconciliation was added. The proposal's independent half is
real and became **OQ-Coarsening**: `C(ρ)`'s smoothing length is unspecified everywhere, strict concavity
makes the resulting bias **one-signed** (coherence over-, dark matter under-estimated), it reaches
**66%** near the knee, and ℓ alone moves NGC 3198 across the entire verdict range. Count held at 6;
REC-038 +1 (8th instance, first with the prior art in the same directory); REC-036 gains a second
definitional gap.

**So what**: the pass's value is not the edit — it is that a plausible, well-argued, arithmetically
correct proposal was checked against the code it describes and did not survive, while the part of it
nobody was arguing for turned out to be the load-bearing gap. And the same discipline caught two of my
own errors, one of which would have fabricated a refutation of a proof that is correct.
