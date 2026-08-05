# The baseline/signal gate: a structural prediction is only testable where the measurement's signal exceeds the model's own baseline error

**Date**: 2026-08-05
**Track**: Publisher (whitepaper lane), back-annotating an explicit request from the site explorer lane
**Status**: REGISTERED — proposed as a standing gate on the refutation ledger
**Origin**: `synchronism-site/explorer/findings/efe-zero-is-not-refutable-by-chae2020-the-baseline-is-off-by-3-dex.md` (2026-08-04), action item #7 ("Back-annotation to the Synchronism research repo")
**Verification status**: the numeric claim below was **re-derived independently in this lane** on the in-repo SPARC mass models before propagation — not accepted from the sibling repo on trust.

---

## The gate

> **Before badging any published result as confirming or refuting a structural prediction, compute the model's own error on the raw observable, at the measurement's own radii, in the measurement's own units. If the model's baseline error exceeds the signal being measured, the correct badge is `not-evaluable` — not `confirmed`, not `refuted`.**

A refutation is a statement about a *residual*. A residual is only meaningful relative to a baseline that fits. When the baseline misses by orders of magnitude, the residual carries no information about the structural claim, and a significance figure attached to it (however large, however well-measured) is significance in the *published* model, not in ours.

## The case that generated it

The framework predicts **EFE = 0 exactly**: `C = C(ρ_local)`, and a uniform external field leaves local ρ unchanged, so the Strong Equivalence Principle holds by construction. Chae et al. 2020 (ApJ 904:51) report a ~4σ sample-level External Field Effect in SPARC, with per-galaxy detections at 11σ (NGC5055, `e` = 0.054) and 8σ (NGC5033, `e` = 0.104).

Read at that level the inference is immediate and was drawn twice in this repo, **in opposite directions, within two days**:

| direction | claim | where |
|---|---|---|
| overclaim of *testability* | EFE = 0 is "sharper and more falsifiable," and "stands in tension with" Chae's ~4σ detection | whitepaper §5.15, executive summary, conclusion (2026-08-03), re-sharpened 08-04 |
| overclaim of *refutation* | "EFE=0 IS a genuine MOND-discriminator that the framework **FAILS**" | `PREDICTIONS.md` B2, `SESSION_FOCUS.md` (2026-08-03) |

Both die on one number that neither computed.

Chae's EFE is a **velocity deficit in the outer rotation curve**, of size **0.046 dex (NGC5055) to 0.083 dex (NGC5033)** — a 10–17% effect. The framework's stated density law, evaluated on the same curve at the same radii, is off by **+3.1 to +4.2 dex**:

| galaxy | Chae `e` (σ) | framework baseline error | EFE signal | baseline / signal |
|---|---|---|---|---|
| NGC5033 | 0.104 (8σ) | +3.14 to +3.71 dex | 0.083 dex | **38× – 45×** |
| NGC5055 | 0.054 (11σ) | +3.62 to +4.18 dex | 0.046 dex | **80× – 92×** |

Grid: γ ∈ {2, 0.489} × h ∈ {0.3, 1.0} kpc; minimum ratio anywhere on the grid **37.8×**. MOND on the same points and the same pipeline sits at **+0.002 to +0.084 dex** — *inside* the signal, which is what makes Chae's measurement meaningful for MOND and meaningless for us.

**Verdict: EFE = 0 is `not-evaluable` against Chae et al. 2020. Refutation count stays at 6. Bucket 0 stays 0.**

### Reproduction

`simulations/sparc_real_data/MassModels_Lelli2016c.mrt` (Lelli et al. 2016), outermost radius with positive reconstructed density; `ρ = Σ/2h` with `Σ = (M/L)·SB_disk`, M/L = 0.5; `ρ_crit = 0.029·V_flat²` M☉/pc³; ledger convention `V = V_bar/√C` with `C = tanh(γ·ln(ρ/ρ_crit + 1))`.

**One normalization trap, recorded because it runs the flattering way.** `ρ_crit = 0.029·V_flat²` is in M☉/pc³ **unscaled**. An assumed 10⁻³ scaling (plausible if one reads V_flat in km/s and reaches for SI-ish units) shifts every baseline error by exactly `½·log₁₀(1000)` = 1.5 dex, *downward*, and understates the ratio by ~3×. It changes no conclusion here — the ratio stays ≥ 20× even then — but the direction of the error is toward the framework, which is the direction that does not get caught by wanting the result.

## Why this is a gate and not an anecdote

Three properties make it reusable:

1. **It is one line of arithmetic** — a division — and it is available before any literature search, any σ, any ΔBIC.
2. **It is direction-neutral.** It killed a claimed *tension* and a claimed *refutation* in the same pass. A ledger that only applies scepticism to favourable results drifts; this one applies to both signs by construction.
3. **It catches a failure the existing rules do not.** This program already requires naming the estimator, pre-registering the criterion, and verifying claims of exactness. None of those fire here: the σ was real, the paper was read at source, the prediction was correctly derived. What was missing is that *the model and the measurement were never put in the same units on the same observable.*

## The propagation mechanism, which is the sharper half

The executive-summary instance of this claim carried its own accurate warning label:

> *"Magnitude not independently re-derived here; see §5.15."*

The warning was true, was written by the propagating lane, was never removed — and the magnitude was the **entire load-bearing content** of the claim. The claim reached four surfaces in two days anyway.

**A caveat naming the unverified quantity does not stop the claim from travelling.** Flagging is not gating. If a claim's decisive quantity is uncomputed, the options are compute it or don't propagate it; annotating it and propagating anyway produces exactly the artifact seen here — a correctly-hedged sentence that is wrong.

This is the companion to the mode recorded 2026-08-04 (*"the sentence you are least tempted to check is the one travelling on your verification of its neighbour"*). Together: **verification does not transfer across sentences, and hedges do not substitute for it.**

## Scope — what this does not do

- It does **not** rescue `C(ρ)`. The reason EFE = 0 is unreachable is that the density law misses the rotation curve by 3–4 dex; that is the *same* failure the ≤25% admixture bound and the 2–5 OOM mean-relation result already record. The prediction is untestable **because** the model is badly wrong on the underlying observable.
- It does **not** make EFE permanently unaskable. The route is the differential completion `∇·[C(ρ)∇Φ] = 4πGρ` (registered 2026-08-04): if it closes the 3-dex baseline gap on SPARC, the EFE question becomes askable for the first time. Until then it is not.
- It does **not** change the count. **6 refutations, Bucket 0 = 0**, unchanged in both directions.

## Net effect on the ledger

Two-signed, and on balance **favourable** — the first correction in this lane in some weeks that moves that way:

- **Removed**: a claimed empirical tension between EFE = 0 and Chae et al. 2020 (both directions).
- **Restored**: "zero remaining active discriminators" stands **unqualified** (the 08-03 text had qualified it).
- **Unchanged**: refutation count 6; Bucket 0 = 0; the constructive result that EFE = 0 follows from the completion's linearity in Φ.

## One correction back to the source finding

The finding describes Chae's deficits as *"a 5–12% velocity deficit."* Its own quoted figures — 0.046 and 0.083 dex — are **9.9% and 17.4%**. Immaterial to every conclusion on both sides, recorded for the same reason this document exists.
