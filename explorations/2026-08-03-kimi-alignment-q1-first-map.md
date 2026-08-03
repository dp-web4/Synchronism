# Q1 — first light of the alignment Bell model: C(Δθ, τ) map + the two audits

**Date**: 2026-08-03 · **Arc**: alignment (`explorations/2026-08-02-kimi-alignment-arc-plan.md`) · **Stage**: Q1, decided
**Sim**: `simulations/kimi_alignment/alignment_bell_q1.py` · **Data**: `simulations/results/kimi_alignment/q1_first_map.json`

## Falsifier first (Q1 kill criteria, registered 2026-08-02 before any run)

- **Q1(a)**: the construction reduces to *relabeling* measurement settings (the
  nonlocal-grid trap, gauge-equivalent, not a mechanism) → kill.
- **Q1(b)**: the correlation requires initial-condition fine-tuning rather than window
  dynamics (the conspiracy boundary: superdeterminism proper is an unfalsifiable alibi,
  not an instrument) → kill.

## The model (abstraction debts named up front, per the river rule)

1D lattice, N=512, open boundaries, ONE global tick, Kuramoto-style local coupling
(dt·4K = 0.12, 16× inside Euler stability). Source at center emits one co-locked pair per
run (both packets carry the same phase offset ψ ~ U(0,2π) — phase lock at creation, on the
substrate). Packets travel ballistically outward as flat-top phase plateaus (16 cells;
plateau interior feels no Kuramoto force, so the phase survives the 184-tick trip — a
sin-bump profile attenuated to ~1/3 and was replaced). Detectors at cells 64/448: 8-cell
oscillator ensembles with strong internal coupling and an own-phase r ~ U(0,2π) drawn per
run per detector (load-bearing — see the null). Settings θ_A, θ_B = phases of reference
oscillators at the detectors relative to the shared substrate phase, drawn LATE (separate
rng stream, after emission + 184-tick travel; order asserted). Measurement: a gate armed
for τ ticks opens *in proportion to the arriving pattern's amplitude* and aligns the
ensemble to the local deviation field (κ·dt = 0.04/tick); outcome =
sign sin(φ̄_det − θ_det), local scalars only.

**Debts**: 1D not 3D; continuous phase not quantized intent; ballistic transport imposed
(local shift), not derived from a substrate wave equation; co-rotating frame — substrate
phase coherence source↔detectors *asserted, not derived*; detector isolation between
measurements idealized; gate pattern-selectivity idealized; open boundaries by fiat. No
verdict below leans on a discarded part.

## C(Δθ, τ) map (600 pairs/cell, seed=1)

| Δθ \ τ | 0 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|
| 0.00 | +.000 | +.033 | −.017 | +.060 | +.237 | **+.503** |
| 0.39 | −.023 | +.027 | +.007 | +.080 | +.207 | **+.507** |
| 0.79 | +.067 | −.003 | +.013 | +.067 | +.263 | **+.353** |
| 1.18 | −.037 | −.047 | −.023 | +.037 | +.127 | **+.270** |
| 1.57 | +.003 | −.043 | −.023 | +.010 | +.087 | **+.087** |
| 1.96 | +.020 | +.067 | +.027 | −.037 | −.080 | **−.053** |
| 2.36 | +.013 | +.050 | +.007 | −.043 | −.103 | **−.227** |
| 2.75 | −.020 | −.043 | +.003 | −.020 | −.147 | **−.350** |
| 3.14 | −.070 | +.027 | +.110 | −.070 | −.180 | **−.507** |

Cell errors ±0.035–0.041. Controls: no-source null max |C| = 0.12 (worst cell a 3σ
outlier *at τ=0*, statistically expected among 27 cells), no-lock max |C| = 0.11 — both
consistent with zero against signal 0.51. Runtime 982.6 s (over the 15-min target by ~10%;
full grid kept anyway, noted).

## Headline numbers

- **τ curve** (the bet's distinctive physics — correlation manufactured IN the window):
  lock factor 0.039 → 0.755, triangle amplitude A 0.02 → 0.54, CHSH S 0.19 → 1.29 as τ
  goes 0 → 64. Monotone, still rising at the longest window — the asymptote is not yet
  reached, and that is said plainly, not extrapolated.
- **Full-lock angular form (τ=64)**: triangle fit A = 0.540, RMSE 0.0796; cosine fit
  A = 0.467, RMSE 0.0803 — **statistically tied; triangle nominally closer**. Two honest
  limits: (i) at A ≈ 0.54 the two forms are barely distinguishable at 9 angular points;
  (ii) the measured curve's peaks are *rounded* (C(0) ≈ C(π/8) ≈ 0.51) — exactly the
  smearing a residual detector own-phase (1−lock ≈ 0.25) produces. Whether the τ→∞ form
  approaches the triangle or the cosine is **Q2's registered question**, and it needs a
  lock extrapolation, not this map's longest row. Expected by any local construction:
  triangle. S(τ) ≤ 2 throughout — the theorem's cap is respected; nothing here is a
  violation claim.

## The audits (methods + results — checked properties, not claims)

- **Relabeling audit — PASS.** Static: outcomes come from `read_outcome(phi_bar, theta)`
  (local scalars, one detector); C from `correlate(outcomes_A, outcomes_B)` (outcome
  arrays only); no expression anywhere takes θ_A and θ_B jointly. Dynamic: 200 runs,
  identical dynamics and settings streams but θ_B swapped post hoc — outcome_A bitwise
  identical 200/200; outcome_B followed its own setting in 95/200 (settings are real
  levers, not decoration). Settings are physical phase references on the substrate
  (house rule 8 satisfied by construction and by test).
- **Conspiracy audit — PASS.** Detector own-phases randomized per run; settings drawn
  after emission+travel (asserted order); no per-cell tuning (K, κ, geometry global);
  robustness probe at (Δθ=π/4, τ=32): C = 0.16 at K_SUB×½ and 0.17 at K_SUB×2 (nominal
  0.26 ± 0.04 — within joint noise, no parameter sensitivity). Correlation needs no
  initial-condition tuning; it appears *only* with the pair (null ≈ 0), *only* with the
  source lock (no-lock ≈ 0), and *grows with the window* — the three signature behaviors
  of the bet, present.

## Verdict

**Q1 PASS.** Neither kill criterion triggers: the construction is not a relabeling
(audit, static + dynamic) and needs no conspiracy (randomized ICs, late settings,
robustness probe). The correlation is manufactured by the alignment dynamics in the
window, with both controls clean. The instrument stands; Q2 is unblocked.

**Next action (Q2)**: the full-lock extrapolation — longer windows / stronger gate /
lock-factor → 1 (e.g., κ ramp or τ to 256 with per-τ lock measurement), then the
triangle-vs-cosine discrimination at A → 1, engaging the Bell theorem's
measurement-independence premise where the model replaces it with dynamics. Cost: same
harness, ~15 min compute, no new code beyond a τ extension and a fit at the asymptote.

## What the reparametrization buys (house rule, answered explicitly)

The bet reparametrizes "measurement" as "finite-window phase alignment on a shared
substrate." At Q1 it buys: (1) a *knob orthodox QM does not have* — the window τ, with a
monotone S(τ) curve as immediate, measurable consequence (Q3's payout material, already
in hand); (2) settings demoted from free labels to physical phases, closing the
relabeling loophole *by construction* rather than by assumption; (3) a named obstruction
list for Q2 (lock asymptote, peak smearing) instead of a vague "classical models give
triangle." What it has NOT bought (honest): any approach to the cosine — A = 0.54 at the
longest window, and the form question is untouched by Q1's pass.

## Ledger scores

- **Theory ledger**: Bucket 0 unchanged — this is a mechanism model, not a prediction;
  no number here constrains the world. S ≤ 2 respected everywhere.
- **Frame ledger**: an instrument built (the first substrate-native Bell harness with
  audited locality of outcomes); a new thing sayable — "correlation is manufactured in
  the window" now has a curve attached (S(τ) measured, not asserted); the
  instrument-independence check (house rule 8) demonstrated working on its first case.

## Honesty block

Kuramoto lattices, coupled oscillators, and classical Bell constructions with triangle
correlations are all known physics; nothing here is novel physics and no Bell inequality
is violated (S ≤ 2 throughout, as any local construction must). The Q1 pass is a frame
ledger event (instrument + mechanism demonstration), not a theory ledger one. The full-
lock form is NOT yet measured — A = 0.54 at τ = 64, peaks rounded by residual own-phase,
triangle/cosine statistically tied — and Q2 could still kill the angular-form hope
honestly; the model's distinctive payout (S(τ) vs weak-measurement QM) is untested until
Q3. Two smoke-test failures were found and fixed before data-taking (gate dilution,
missing amplitude trigger — both documented in the sim), which is said here because the
controls that caught them are the same ones a skeptic should check first.
