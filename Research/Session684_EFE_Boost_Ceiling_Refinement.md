# Session 684: EFE Sector — Boost-Ceiling Fork Pattern-Matches the Galaxy/Cluster Family `[ACTIVE-MRH]`

**Date**: 2026-06-03
**Type**: Verification of a new visitor-channel result and integration of the EFE sector into the existing fit-XOR-discriminate pattern (S661, S678/S683).
**Trigger**: `Research/proposals/efe_boost_ceiling_closure.md` (site explorer track, 2026-06-03), reproduction `synchronism-site/explorer/scripts/efe_boost_ceiling_closure.py` on real SPARC N=2807.
**Status**: `[ACTIVE-MRH]` — refines a prior `[AUDITED-NEGATIVE]` finding pattern (S661/S678/S683), substance preserved, the EFE sector is added to that family with the same structural form.

---

## §1 — What the proposal claims

The proposal closes the EFE/TDG (external field effect / tidal dwarf galaxy) sector by showing that one "boost ceiling" parameter `B_max` controls both:
- RAR fit quality (improves as `B_max` rises)
- TDG EFE distinctness from MOND (decreases as `B_max` rises)

The two effects are anti-correlated. No value of `B_max` simultaneously fits the SPARC RAR and keeps the EFE distinct from MOND. The fork is **fit XOR discriminate** — the same shape as S661's RAR transition-shape result (γ=2 refuted / γ-free → MOND) and S678/S683's cluster diagnosis (the only escape from the cluster failure is to re-introduce the acceleration variable, which is MOND again).

Reported numbers (Lelli-McGaugh-Schombert 2016, N=2807, 10% velocity-error cut):

| B_max | RAR RMS (dex) | TDG Δσ = σ_MOND − σ_Sync |
|---|---|---|
| 3.17 (bounded Hill, 1/Ω_m) | 0.227 | 8.1 km/s (distinct) |
| 20.7 (joint RAR best-fit) | 0.146 | ~2 km/s |
| ∞ | 0.146 | 0.0 km/s (MOND) |

42% of SPARC RAR points require boost > 3.17; max observed boost ~34×. The bounded form is independently falsified by the RAR shape.

## §2 — Verification (`session684_efe_boost_ceiling_check.py`)

I do not re-run on real SPARC — the explorer track did that on the public catalog and the load-bearing computation is in their script. I verify the **structural** claim against the McGaugh-Lelli-Schombert exponential interpolating function `ν_e(y) = 1/(1 − exp(−√y))`, `y = g_bar/a₀`, on a synthetic log-uniform `g_bar` grid spanning the SPARC range.

**(A) Boost distribution.** On a log-uniform grid log₁₀(y) ∈ [−4.5, 1.5]: 61% of points have McGaugh boost > 3.17, median boost 6.14, 95th-percentile 126. On real SPARC the proposal reports 42% and max ~34×. The difference is the distributional weight — my synthetic grid samples the deep-MOND tail more heavily than SPARC does. The qualitative claim (a substantial fraction of SPARC sits above any reasonable bounded ceiling at 3.17) is confirmed.

**(B) RAR RMS vs B_max** (fit knob). Bound the boost at `B_max`, take residual vs the McGaugh ν_e prediction:

| B_max | synthetic RAR RMS (dex) | fraction clipped |
|---|---|---|
| 3.17 | 0.776 | 61.0% |
| 5.0 | 0.648 | 53.3% |
| 10.0 | 0.468 | 42.4% |
| 20.0 | 0.310 | 32.0% |
| 50.0 | 0.137 | 18.5% |
| ∞ | 0.000 | 0.0% |

Monotone decreasing — confirmed. Absolute RMS numbers exceed the proposal's 0.227 because my synthetic grid samples the deep-MOND tail more than SPARC, so the bounded form's clipping bites harder; the proposal's real-SPARC numbers are necessarily smaller. **Direction of the trend is the load-bearing fact and is exact.**

**(C) TDG EFE proxy.** At a representative deep-MOND TDG operating point (y=10⁻³, MOND boost ~32, MOND σ ~14.5 km/s for the bounded-Hill calibration): bounded boost → smaller inferred σ. Δσ = σ_MOND − σ_bounded is monotone *decreasing* in B_max (9.95 → 0 across B_max ∈ {3.17, …, ∞}). Confirmed.

**(D) Anti-correlation.** RAR RMS and TDG Δσ move in opposite directions with B_max. There is no value of `B_max` that simultaneously gives MOND-grade RAR fit and a distinct EFE prediction. The structural fork is real and the proposal's headline ("no single boost ceiling both fits the RAR and keeps the EFE distinct") is confirmed.

## §3 — Integration with the existing pattern

The proposal pattern-matches into a unified picture of the galaxy/cluster sector:

| Sector | Fork | Distinct branch (refuted) | Fitting branch (= MOND) |
|---|---|---|---|
| Galaxy RAR shape (S661) | γ pinned | γ = 2 (ΔBIC = +184 on SPARC) | free-γ → 0.49 ≈ McGaugh |
| Cluster bridge (S678/S683) | variable choice | C(ρ) with one density scale | C(a) restoration ≡ MOND |
| EFE sector (S684, this) | B_max | bounded Hill (B_max ~ 3.17, refuted by RAR shape) | B_max → ∞ ≡ MOND EFE |

Same structural fork in every sector: any choice that *distinguishes* the framework from MOND at galaxy/cluster scales is independently falsified; any choice that *fits* the data IS MOND, with the framework's distinguishing content collapsing to a relabel. The S684 result strengthens this from a per-sector observation to a *family pattern*. It also sharpens the C(a) → C(ρ) migration cost statement S683 made: every choice that survives is either (a) a relabel of MOND in the original C(a) sector, or (b) a non-fitting bounded form. There is no "different physics" position in this region of model space that survives both the RAR and the EFE.

The proposal's third recommendation — that the only structural way to reopen the sector is a 2-parameter coherence form (a density scale and an acceleration scale) that decouples the RAR fit from the EFE ceiling — is consistent with the S678/S683 finding that the cluster failure is wrong-variable, not one-scale. A 2-scale form does not address the wrong-variable obstruction at cluster scales (g_bar is non-local in ρ regardless of how many density scales you add), so the prediction is that any 2-scale form passing the EFE test will fail at clusters in the wrong-variable way. This is testable but not in this session.

## §4 — What this session does not output

- **No verdict on the broader galaxy/cluster/EFE closure.** The substantive findings stand as `[AUDITED-NEGATIVE]` on the old R(I) substrate. The active substrate reformulation (saturation reframe with independent **J**, Phase-1 sim) is upstream of this question; whether it has anything to say about the EFE sector is open and downstream of fleet sim results.
- **No retag of S661 or S678/S683.** Those session docs are already in MRH-relationship taxonomy after S679; what this session adds is *integration with* a new sibling finding (S684 = EFE companion to S661 = RAR), not a status change.
- **No execution of the proposal's archive actions** (commit to field equation form; acknowledge RAR → EFE-MOND entailment; explore 2-scale forms). Those are operator/coordinator work. My role: verify load-bearing structure, integrate framing.
- **No re-run on real SPARC.** The explorer track's `efe_boost_ceiling_closure.py` is the canonical computation; my synthetic check verifies the structural direction, not the absolute numbers.
- **No cumulative tally.** Per the S679 discipline.

## §5 — Files

- `Research/Session684_EFE_Boost_Ceiling_Refinement.md` (this document)
- `simulations/session684_efe_boost_ceiling_check.py` (structural check on McGaugh ν_e; monotonicity verification)

## §6 — So what (under the frame-doc disciplines)

A research-repo autonomous session can verify the structural direction of a visitor-channel result without re-running the absolute numbers, and integrate it as a sibling of prior findings under one pattern. This session did that for the EFE sector: B_max is a single anti-correlated knob, the fit-XOR-discriminate fork is real, and the EFE sector now joins S661 (RAR shape) and S678/S683 (cluster bridge) as a third instance of the same family pattern — any choice in the galaxy/cluster/EFE region of model space that distinguishes the framework from MOND is independently falsified, and any choice that fits the data collapses onto MOND. The pattern is now structural (three sectors, same fork shape) rather than per-sector.

The proposal's broader implications (field equation commitment; 2-scale exploration; preprint framing of the closed sector) are operator/coordinator work, not autonomous-session output.
