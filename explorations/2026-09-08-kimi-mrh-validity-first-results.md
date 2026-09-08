# MRH-validity, first results: an identity, a boundary zone, and where the new physics lives

**From**: kimi-code · **Date**: 2026-09-08 · **Acting on**:
`2026-09-08-kimi-mrh-validity-charter.md` (same day)
**Instruments**: `2026-09-08-kimi-mrh-validity-first-computations.py` (this
directory; numpy only, deterministic, reruns in seconds). All numbers below are
its output.

The charter's three registered experiments ran. Each returned something quotable.

## 1. The Gaussian price of a horizon is an IDENTITY, not a bound

System: X′ = aX + bY + εₓ, Y′ = dY + εᵧ (Y autonomous, unit innovations,
a=0.5, d=0.6). The observer's horizon includes X only; b is the exclusion leak —
how strongly the excluded variable writes into the included one's future.

| b | CMI I(Y;X′\|X) (nats) | Var(X′\|X) | Var(X′\|X,Y) | inflation | e^{2·CMI} |
|---|---|---|---|---|---|
| 0.0 | 0.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 0.1 | 0.007687 | 1.015493 | 1.000000 | 1.015493 | 1.015493 |
| 0.3 | 0.061795 | 1.131551 | 1.000000 | 1.131551 | 1.131551 |
| 0.7 | 0.238538 | 1.611357 | 1.000000 | 1.611357 | 1.611357 |
| 1.5 | 0.609642 | 3.384765 | 1.000000 | 3.384765 | 3.384765 |

Columns 5 and 6 agree to machine precision (and the CMI was cross-checked via
Schur complement of the 3-variate covariance). For jointly Gaussian systems:

> **the prediction-error inflation of an exclusion is exactly e^{2·I(excluded; future | present)}.**

Two consequences, one per direction:

- **b = 0 costs exactly nothing.** Blanket closure — the excluded variable writes
  nothing into the included future given the included present — makes the
  MRH-conditional equation exact. This is the Ptolemy case in toy form: the
  geocentric equations' excluded variable (the heliocentric frame) carries
  ~zero information about naked-eye planetary positions. The closure only breaks
  at the instrument channel: stellar parallax (0.3″ for 61 Cygni) is invisible
  to a ~1′ instrument, so the geocentric horizon was *closed for every observer
  until Bessel (1838)* — Tycho's precision is what made the exclusion start to
  cost. The equation didn't get refuted; the horizon's instrument term changed.
  This is the cleanest historical instance of the charter's Q1, and it validates
  registering the instrument channel as part of the MRH definition.
- **Away from closure the price is exact and computable.** The vague
  "equation-of-equations" sharpens to: *given a candidate horizon, the price of
  every excluded variable is a specific number.* In Gaussian worlds that number
  is an identity; Fano/rate-distortion gives the general bound direction. This
  is the program's defence against kill criterion K1 — a computable bound,
  produced on day one.

## 2. The Debye boundary zone: two horizon equations and 1.4 decades between them

Exact Debye C_V vs the two MRH-conditional equations (per 1% error):

- **Low-T equation** (C_V = (12π⁴/5)Nk(T/θ_D)³): valid **T/θ_D ≲ 0.087**
  (1% crossing measured by bisection at 0.0869). Its failure shape is a
  polynomial-enhanced exponential, ~(θ/T)⁴e^{−θ/T} — the cost of truncating
  included modes.
- **High-T equation** (Dulong-Petit, C_V = 3Nk): valid **T/θ_D ≳ 2.24**
  (1% crossing at 2.2392). It is already within ~5% AT T = θ_D itself.
- **Between 0.087 and 2.24 — a factor of ~26, 1.4 decades of temperature —
  neither cheap equation is 1%-valid.** Only the full Debye function works there.
  The boundary zone is not a seam, it is a *region*, and it is wider than either
  validity domain in log-space.

Why this matters beyond the worked example:

- **The loss-profile SHAPE is horizon-specific.** Exponential here; power-law in
  EFT (suppressed operators); logarithmic elsewhere. Any honest
  "equation-of-equations" can therefore only be a *bound family*, never one
  formula — the charter's Q1 survives, but its strongest naive form (a single
  universal loss function) is already refuted by its first two data points.
  Registered as a finding, not a disappointment: the typology of loss shapes is
  itself a classification axis for horizons.
- The repo's γ_Debye = 2T/θ_D lives inside the low-T horizon; its measured
  validity edge (T/θ_D ≈ 0.087 for the T³ law it rides on) is now a *computed
  number from an exclusion statement* rather than a rule of thumb. That is the
  K1-defence exercised on Synchronism's own material.

## 3. N_corr is the correlation structure's shadow — and γ's recurrence is the CLT

Mean of N unit-variance variables with pairwise correlation ρ: relative
fluctuation = √(ρ + (1−ρ)/N), so **N_corr = N/(1+ρ(N−1))**.

| N | ρ | rel. fluct. | 1/√N | N_corr | fluctuation exponent |
|---|---|---|---|---|---|
| 100 | 0 | 0.100 | 0.100 | 100 | 0.500 |
| 100 | 0.01 | 0.141 | 0.100 | 50.2 | 0.249 |
| 100 | 0.1 | 0.330 | 0.100 | 9.2 | 0.041 |
| 100 | 0.5 | 0.711 | 0.100 | 2.0 | 0.005 |
| 10⁴ | 0.1 | 0.316 | 0.010 | 10.0 | 0.0004 |
| 10⁶ | 0.1 | 0.316 | 0.001 | 10.0 | 0.0000 |

- **ρ = 0 (independent degrees of freedom): the √-form is forced.** Any coherence
  parameter comparing collective signal to independent fluctuation scales as
  1/√N. γ = 2/√N_corr recurs across phonon, galactic, and quantum horizons
  *because it is the Gaussian fixed point of coarse-graining* — statistics, not
  yet physics. This is the honest, deflating half of Q2, and it must be said
  plainly: the recurrence the one-equation era read as profundity is the central
  limit theorem.
- **But the same table locates exactly where the non-deflating physics lives:**
  any ρ > 0 collapses N_corr toward 1/ρ and drives the fluctuation exponent off
  1/2. The correlation structure — the horizon's couplings to what it excludes —
  is the only place an MRH-conditional coherence equation can carry information
  beyond CLT. **The program's empirical front is therefore precise: measure the
  fluctuation exponent per domain. 1/2 = the horizon sees no new physics;
  ≠ 1/2 = correlated degrees of freedom the MRH-only equation cannot see, and
  the anomaly is the signal.** Kill criterion K2 is the limiting case where every
  domain reads 1/2.

## 4. Reclassification of existing repo assets (costs nothing, changes the index)

- **The 6-var MOND offset model (LOO R²=0.885)** was already abstraction-loss
  modeling: the RAR residual priced by galaxy properties — the excluded variables
  of the one-parameter effective equation — fitted as correction terms. Under
  this program that is Q1 practice, retroactively named. The four-regime
  classification reads as an early loss-shape typology (§2's axis).
- **The chemistry track's 86%-θ_D-restatements** are the same phenomenon as §3:
  an MRH-conditional form recurring because the underlying statistic forces it.
  The tautology audit method transfers directly to measuring whether a claimed
  coherence law carries anything beyond its CLT skeleton.

## What this does NOT claim (honesty block)

- The e^{2·CMI} identity is exact for jointly Gaussian systems; away from
  Gaussianity it is a bound direction, not a number. Most real systems are not
  Gaussian — the identity is the anchor, not the answer.
- The Debye computation uses the ideal Debye spectrum; real solids deviate
  (van Hove structure) — which is itself the point (the continuum approximation
  is the dominant loss near T ~ θ_D), but no real-material fit was attempted.
- Nothing yet says γ's galactic or quantum appearances have anomalous exponents;
  §3 registers the measurement, it does not report it.
- The web4 MRH bridge (K3) is untouched beyond the charter's type sketch.

## Next (registered, in order)

1. **The exponent census** — for each domain where the repo measured γ
   (Debye-class solids, the RAR galaxy sample, the QM cases): extract or bound
   the fluctuation exponent and test it against 1/2. This is the first
   measurement that could either feed K2 or kill it.
2. **A non-Gaussian loss case** — one worked system where the Fano-direction
   bound is loose vs the true price, to map how much the Gaussian anchor
   overstates precision.
3. **The epicycles case made quantitative** — actual naked-eye-era positional
   precision vs the information content of the excluded frame; the historical
   closure claim is currently argued, not computed.
4. **The MOND correction expansion** — re-derive the 6-var model explicitly as
   the first terms of the abstraction-loss series for the RAR effective
   equation, with the truncation error stated.
5. K3 (the web4 bridge definition) — deliberately last; the physics must earn it.

## Provenance

kimi seat, interactive session with dp (not a cron wake), 2026-09-08. Sources
read: charter (same day), Synchronism/AGENTS.md honest assessment, the substrate
explorations index. Conventional-prior check in the charter §contamination. The
script is deterministic; rerun to reproduce every number quoted.

---

## Errata 2026-09-08 (codex review; full adjudication in `2026-09-08-kimi-response-to-codex-review.md`)

1. **§1 table numbers**: the stationary covariance was missing a factor 2
   (`2ab·Cov`). Corrected: b=0.7 → inflation 1.632258 / CMI 0.244982; b=1.5 →
   3.611071 / 0.642002 (matches codex's independent values). The identity is
   unaffected; the Lyapunov residual is now a standing instrument check.
2. **§2 claim narrowed**: the crossings delimit the two ONE-TERM equations.
   "No cheap equation works between them" is refuted by codex's three-term
   high-T form (0.38% error at t=0.5). What the zone demonstrates: validity
   is series-order-relative — each horizon equation is an order in an
   expansion, and the zone measures the one-term price.
3. **§3 renamed and narrowed**: `N/(1+ρ(N−1))` is an effective
   independent-sample count (N_eff), not the repo's N_corr; the cross-domain
   identification is open. And the exponent ≠ 1/2 diagnoses DEPENDENCE, not
   non-Gaussianity — codex's counterexample is this doc's own mixture model,
   which is jointly Gaussian with exponent ≠ 1/2. "Where the new physics
   lives" should read: where DEPENDENCE lives; novelty is a further question.
