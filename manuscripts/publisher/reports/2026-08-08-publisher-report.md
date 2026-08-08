# Publisher Daily Report — 2026-08-08

**Headline**: the whitepaper asserts a refutation count of **6** four times and, until today, enumerated it
**nowhere** — and two of its six members are described nowhere in the whitepaper at all. When the collapse
everyone has been asserting is actually measured, it is **6 → 5**, not the "3–4" now circulating. Count
**HELD at 6**, no row merged; the naming decision remains dp's.

---

## Phase 0: Publication Recommendations

### New Recommendations
**None.** No new complete arc; Archivist reports the **4th consecutive deliberate zero** on Synchronism core
Sessions. All of today's input arrived as proposals and sibling-repo logs, not numbered sessions.

### Status Changes

| Rec | Change | Readiness |
|---|---|---|
| **REC-2026-038** (repository-mediated continuity manuscript) | **+1 strength, +1 weakness** — 10th instance of the documented class, and the **first externally-sourced double** | **HELD 0.93** (3rd day) |
| **REC-2026-039** (bounded-boost class exclusion) | **+1 strength, +1 weakness** — its class-exclusion framing is now the *measured*-correct level of description; it also inherits a counted single-dataset/single-M/L dependence | **HELD** |
| **REC-2026-036** (experimental test catalog) | **+1 weakness** — 3rd instance of the ID-keyed-audit blind spot, and the first about the catalog's own headline rather than a missing row | **HELD 0.68** |

### Upcoming Candidates
None new. §1b scan of all seven surfaces: `Research/papers/` unchanged since 07-23,
`Research/preregistrations/` unchanged since 07-30, `synchronism-site/explorer/findings/` newest is
08-07 08:13 (read, via the proposal that carries it), `explorer/topics/` nothing new, `explorations/` one
new file (08-07 kimi saturation-horizon synthesis, hive-organs arc wake 1 — arc-charter scope, not
publication scope).

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **1 integration**

**Status**: updated. **Sessions reviewed through**: no new core Sessions (4th deliberate zero); inputs were
three 08-07 proposals, one exploration, and two sibling-repo lane logs.

**Section affected**: `05-quantum-macro/15-dark-matter` (a five-paragraph amendment appended to the
2026-08-03 note), plus its CHANGELOG.

**What the amendment does.**

1. **Enumerates the six for the first time anywhere in this whitepaper**, each with the statistic it fired
   on and the data it used: (1) RAR transition shape at γ=2, SPARC, ΔBIC=+184; (2) BTFR slope (TEST-09),
   SPARC, |Δn|=0.41 vs a 0.3 bar; (3) dwarf DM-fraction ceiling (TEST-10), SPARC, 69% above
   f_DM = 1−Ω_m = 0.685; (4) RAR environment dependence (TEST-08), SPARC × Cosmicflows-4, r²=0.0001 vs a
   registered >0.20; (5) SPARC × Cassini joint squeeze (TEST-11), +17.95σ; (6) observer-relative Bell/CHSH
   substrate bet (B1), S=1.98 ≤ 2 on both arms.

2. **Records that two of them appear nowhere in the whitepaper body** — verified by grep over
   `whitepaper/sections/`: `3.17`, `68.5`, `0.927`, "boost ceiling" → **0 occurrences each**; `CHSH`,
   `no-signaling`, `Kuramoto` → **0 occurrences each** (the only "Bell" is §5.4's exposition of Bell's
   theorem and a prior-art specificity list).

3. **Measures the collapse** (`simulations/publisher_20260808_test09_test10_independence.py`, in-repo SPARC
   mass models, 123 galaxies at Q ≤ 2, i > 30°; pipeline conventions copied verbatim from
   `test09_btfr_bounded_boost_real_sparc.py` and `test10_dwarf_dm_fraction_ceiling.py` so any difference is
   a difference in question, not in pipeline):

   - **TEST-09 and TEST-10 are the same inequality**, not two corollaries of a shared root.
     `f_DM = 1 − (V_bar/V_obs)² = 1 − g_bar/g_obs = 1 − 1/B`, so `f_DM > 1−Ω_m` ⟺ `B > 1/Ω_m`.
     `max|f_DM − (1 − 1/B)| = 2.2×10⁻¹⁶`; the two criteria select **93 and 93** galaxies with **0
     disagreements**, and **88 and 88** with 0 disagreements on the error-weighted outer-3 point.
   - **TEST-09's registered kill does not fire where the ceiling is not binding.** Deleting every TEST-10
     firing galaxy: observed n = 3.39 ± 0.23 vs framework 3.38 ± 0.18, deviation **0.01 (0.0σ)** against a
     0.3 bar. Full sample: 3.79 ± 0.10 vs 3.35 ± 0.07, deviation 0.44 (3.5σ). Exceeders only: 0.33 (2.6σ).
   - **The leave-out test is underpowered, and the amendment says so.** σ on the non-exceeder deviation is
     0.29, so the 0.3 threshold sits **1.0σ** out; the mass lever arm shortens 3.45 → 3.06 dex. It could
     not have demonstrated independence had independence existed. Verdict stated as **"the burden of the
     two-row convention is unmet, not disproved."**

4. **Corrects two figures circulating in the sibling repo**, both asserted rather than computed, both
   overstating the collapse: *"honest independent-root figure: 3–4, not 6"*
   (`synchronism-site/maintainer/logs/2026-08-07.md`) and *"at most four independent roots"*
   (`…/visitor/logs/2026-08-07.md`). On the axis actually measured it is **6 → 5**, one pair.

5. **Names the sharper axis, which nobody was counting: the dataset.** Five of the six run on SPARC; three
   of the scripts behind them were opened and share a fixed Υ_disk = 0.5, Υ_bulge = 0.7
   (`test09…`, `test10…`, `simulations/sparc_tanhlog_profile.py`). TEST-09 and TEST-10 additionally share
   *every* parameter of one force law (`C(a)=Ω_m+(1−Ω_m)x/(1+x)`, Ω_m=0.315, φ=1.618034, a₀=1.05×10⁻¹⁰).
   Every non-SPARC channel enters only jointly or by reference — Cosmicflows-4 as TEST-08's predictor,
   Cassini as TEST-11's second constraint, the textbook CHSH violation as B1's comparison — and the sixth
   analyses **no external dataset at all**, so the site's *"executed on external data… astronomical,
   ephemeris, and laboratory"* is loose on its third data type. **An M/L shift moves five rows together.**

**Deliberate non-action.** Count **HELD at 6**, no row merged, **Bucket 0 = 0**. Whether TEST-09/TEST-10 are
*counted* as one refutation or two is **open question 4 of
`Research/proposals/boost_ceiling_provenance_and_class_exclusion.md`, registered 2026-07-28 and explicitly
gated on dp** — a naming convention, and dp's to set. What the edit supplies is what that question asked for
and did not have: an enumeration and a measurement.

**Terminology concerns**: none. No canonical term redefined; MRH, C(ξ), γ, Intent, Entity untouched.

**The other three 08-07 proposals, adjudicated:**

| Proposal | Whitepaper action |
|---|---|
| `A_calibration_coarse_graining_reading_SUPERSEDED_20260807.md` | **None needed — convergent.** Its §3 identity (`x = (3/16π²)β_J²[V_c/V_flat]²`, ceiling 0.019β_J²) is the same identity this lane derived independently on 08-07, and its β_J·R₀ = 0.315 kpc / "0.8% from 317 pc" matches this lane's 0.76%. §6.4's correction already landed 08-07. Two lanes, two derivations, one result. |
| `efe_sign_is_forkless_the_blocker_is_a_category_error_and_the_test_has_no_power_20260807.md` | **None — and it reinforces the 08-05 withdrawal.** EFE = 0 under all three force laws, executed at e ∈ {0, 0.05, 0.5, 5.0}; the deciding fork is C's *argument*, not how C multiplies. Its §7 says explicitly "refutation count stays at 6". Its archive correction (Session215's EFE=0 *argument* is a non-sequitur while its conclusion holds) touches an archive session, not whitepaper text. Its §5 point — EFE=0 is the SEP, which GR and ΛCDM also predict, so it is refutation-only — is worth carrying if EFE=0 is ever registered; nothing in the whitepaper currently claims it selects. |
| `force_law_fork_blocks_efe_registration_and_makes_count_convention_dependent_20260807.md` | **Count axis integrated (above); force-law axis NOT integrated.** Its EFE half is superseded by the proposal above it, same day. Its count half is the item this report acts on — but note the amendment does not adopt its "3–4", which the measurement contradicts. |

### Web4 Whitepaper — **no change, verified**

**Status**: current. **Repos checked**: `web4` (all branches, post-fetch).

- **3 commits in window**, all `docs/audits/` **only** — C332 entity-types, C334 errors, C336
  security-framework, 8th deltas. **Zero** whitepaper scope, **zero** `web4-standard/` mutation (C336
  records the 8th consecutive zero-mutation pass across that rotation).
- **PAIRED-CHANNELS fired by judgement and returns no action, on a new channel.** C336's finding is a
  retired crypto token surviving in an outward-facing artifact (`P-256` where the spec says `ECDH-P256`,
  fixed 65 days ago; plus `AES-256-GCM` vs `AES-128-GCM`, `SHA-3-256`). Swept the whitepaper for all of it —
  `P-256|ECDH|AES-(128|256)-GCM|SHA-3|W4-(BASE|FIPS|IOT)-1|cose:` over `whitepaper/sections/`: **0 hits.**
  The only matches in the repo's whitepaper tree are in `whitepaper/archive/sections-2026-07-09-pre-rewrite/`
  (frozen) and in PUBLISHER_CONTEXT prose. **Fourth independent confirmation of the structural-zero prior**,
  and the first from the crypto-suite channel specifically: the Web4 whitepaper defines no algorithm tokens,
  so a suite-level drift cannot reach a sentence in it.
- Last commit touching whitepaper sources remains **`57625ed` (2026-07-30)**.

**Terminology concerns**: none. LCT, MRH, T3, V3, ATP, ADP, R6/R7 all intact.

---

## Gates

| Gate | Result |
|---|---|
| Build (`make-md.sh`) | **PASS** — 7,529 lines, 680K |
| Claims freeze (`claims/render_claim_surfaces.py --check`) | **exit 0** — checked 10 claims; v1 freeze verified |
| Edit present in local build | `count of test IDs` → **1/1** (build monolith / docs monolith) |
| Artifact churn (both numbers, per the 08-01 rule) | content **30 lines**; raw **23,412 lines** ⇒ **churn — tree restored, CI builds.** 9th firing |
| Lone-CR sweep (`git grep -lP '\r(?!\n)'`) | hits are **binary only** (PDF/PNG); no text path |
| `recommendations.json` diff | **14 insertions / 9 deletions** — no mass rewrite |

---

## Summary

No new arcs and no new recommendations; the whole pass is one integration and one verified no-change. The
integration is a defect in my own bookkeeping: this whitepaper has asserted "refutation count unchanged at
6" four times without ever printing the list, and two of the six were never described in it — so the number
was unauditable from the document that states it. The list is now there, and the collapse everyone has been
asserting is measured rather than asserted: TEST-09 and TEST-10 fire on **the same galaxies by algebraic
identity** (93/93, 0 disagreements, residual 2.2×10⁻¹⁶) and TEST-09's registered kill has **no power** where
the ceiling is not binding (0.01 vs a 0.3 bar) — but the leave-out test's own σ is 0.29, so the honest
verdict is *burden unmet, not collapse proved*, and the collapse is **6 → 5**, not the circulating 3–4. The
more decision-relevant finding is the one nobody was counting: **five of six rows run on SPARC at a fixed
M/L, and the sixth analyses no external dataset at all** — so an M/L shift moves five "independent"
refutations together.

**So what?** The transferable rule is not about this count. It is that the shared root was **on file since
07-15, registered as open on 07-28, and still re-filed as novel twice on 08-07 by two independent non-author
readers** — one of whom diagnosed exactly this mechanism the same day, for a different row. Re-discovery
tracks whether a claim can be **walked**, not whether it is remembered. So: *a bare tally is a claim about
independence — publish the list or drop the number.* And the direction matters for REC-2026-038, whose
standing worry was that its instance count is an observer artifact: today's instance is externally sourced,
doubled, and against an item already adjudicated open, which is the first evidence that the rate is a
property of unenumerated claims rather than of one observer's attention.
