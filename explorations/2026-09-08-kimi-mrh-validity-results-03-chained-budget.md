# MRH-validity, results 03: the chained budget — and the repo's "strongest" model is mostly the Υ degree of freedom

**From**: kimi-code · **Date**: 2026-09-08 · **Acting on**: sub-arc map rung 3,
incorporating the codex-review adoptions (identifiability column, bias column,
galaxy-holdout discipline)
**Instrument**: `simulations/sparc_real_data/sa2_rung3_chained_budget.py`
(deterministic). Cloud: 149 galaxies (≥5 points each after Q≤2/Inc≥30°/3σ cuts),
**same hash split as V4: 112 train / 37 test**; a₀ fitted on train only
(1.05e−10 m/s², literature-consistent), Υd per galaxy from its own data
(median 0.475).

## The chained variance budget (37 test galaxies, 740 points)

```
anomaly Var                          0.07994 dex²
  1-param function + Υ prices        0.06989   (87.4%)
  instrument floor (errV)            0.00362   ( 4.5%)
  6-var series prices (out-of-sample)0.00044   ( 0.5%)
  within-galaxy residual             0.00538
REMAINDER (res − inst − series)      0.00599   ( 7.5%)
E after series                                        0.925
test residual mean                                    −0.0069 dex   (bias column, codex T6)
IDENTIFIABILITY (codex T5): REMAINDER is a COMPATIBLE INTERVAL [0, 0.00599]
for excluded-sector signal — this channel (rotation curves + 3.6µm photometry)
cannot distinguish unmodeled included-set structure from a sector.
```

## The finding: Session 484's 6-var model does not survive the chain

The Session-484 feature set (logV, logL, c_V, f_gas, logV×c_V, logL×f_gas —
the repo's "strongest" model, R²=0.945 / LOO R²=0.938 with **fixed** Υ) fitted
on train per-galaxy mean residuals and evaluated out-of-sample:

- **train R² 0.305, train LOO R² 0.061** — versus 0.945/0.938 measured with
  fixed Υ;
- out-of-sample, it prices **9.5% of the between-galaxy remainder**
  (0.00467 → 0.00423), i.e. 0.5% of the anomaly.

The chain catches the double-count: with Υ fixed at 0.5, per-galaxy
acceleration-scale offsets are large and galaxy-level features (luminosity,
gas fraction — the variables that SET the baryonic mass scale) price them
beautifully. Once Υ is free per galaxy, those same offsets are absorbed by Υ
itself — and the 6-var series has almost nothing left to explain. **Most of
the 6-var model's explanatory content was the stellar-mass-normalization
degree of freedom wearing galaxy-property clothes.** Its celebrated
logL×f_gas interaction ("the hidden variable," Session 482) is, under the
chain, largely the luminosity-dependence of what Υ must be — not independent
structure on top of it.

Two honest caveats, both registered: (1) my Υ-RMS fit centers each galaxy's
residual, so the between-galaxy component is small *by construction* — the
series never had access to the variance the original model priced; the two
readings ("6-var was Υ" vs "the chain pre-absorbed the target") are the same
fact stated from opposite ends, and the fact is that **Υ and the 6-var offset
structure are the same information**, which is itself a real result about the
RAR's galaxy-level scatter. (2) Session 484's NN-autocorrelation elimination
(r 0.46→0.005) is a structural claim the chain does not address; retesting
autocorrelation after Υ (not after the series) is the registered follow-up.

## What this does to the epicycle index story

The budget, not E, is now the deliverable — E's after-series 0.925 prices a
component that is 0.5% of the anomaly. The honest summary for dp's
diagnostic: **87.4% of the anomaly is priced by the one-parameter function
plus per-galaxy Υ; ~4.5% is instrument; 7.5% remainder, of which galaxy-level
properties explain almost nothing out-of-sample — and the remainder is an
identifiability interval, not a sector detection.** The within-galaxy piece
(0.00538) is where rung 4 goes: radial structure, which no galaxy-level model
can touch.

## What this does NOT claim

- Not a refutation of the 6-var model as a *predictor* of RAR offsets when Υ
  is unavailable (photometry-only surveys) — that use stands; the finding is
  about its information content relative to Υ, not its utility.
- The within/between split uses unweighted galaxy means (approximation
  noted); the qualitative result is insensitive to the weighting choice.
- c_V here is the halves-by-R concentration; Session 484's exact definition
  may differ in detail — the feature set is theirs, the conventions are mine,
  stated in the instrument.

## Next

Rung 4: the within-galaxy remainder's radial structure (does the residual
correlate with R, with local SB, with radius in scale lengths?) — the last
place included-set structure can hide in this channel. And the registered
follow-up: NN-autocorrelation after free Υ.

## Provenance

kimi seat, interactive with dp, 2026-09-08. Instrument deterministic; rerun
for every number. Markov Phase 3 read per the queue (its CMI-horizon
formalism predates and frames this budget's target-relative definition).
