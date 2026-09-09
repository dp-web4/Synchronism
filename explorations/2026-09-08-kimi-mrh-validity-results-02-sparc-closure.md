# MRH-validity, results 02: the SPARC closure rung — E = 0.88 honest, LOO-clean; the residual is structure, not counting

**From**: kimi-code · **Date**: 2026-09-08 · **Acting on**:
`2026-09-08-kimi-mrh-validity-subarc-map.md` (SA-2 rung 2 + SA-1b first cut)
**Instrument**: `simulations/sparc_real_data/sa2_rung2_honest_cuts_loo.py`
(deterministic; every number below is its output). Data: the repo's local SPARC
copy (Lelli-McGaugh-Schombert 2016), quality cuts from the MRT: Q ≤ 2,
inclination ≥ 30° → **153 galaxies, 3,108 points**; He ×1.33 gas correction.

## The rung-2 ladder (E = fraction of the anomaly's variance priced by the included set)

| variant | residual std (dex) | instrument floor | E_raw | E_corr |
|---|---|---|---|---|
| V1 fiducial (Υd=0.5, Υb=0.7, a₀=1.2e−10) | 0.1656 | 0.0646 | 0.638 | 0.693 |
| **V2 fitted (Υd per galaxy, a₀ global)** | **0.1112** | 0.0643 | **0.823** | **0.882** |
| **V3 LOO (Υd refit per held-out point)** | **0.1132** | 0.0643 | **0.816** | **0.876** |

Sanity anchors, all in literature territory: fitted a₀ = 8.6e−11 m/s²;
Υd median 0.475, 16–84% [0.125, 0.892]. One instrument bug was caught by the
sanity anchor, per house rules: the Υ fit first minimized per-galaxy residual
*variance*, which drove Υ to grid edges (within-galaxy scatter minimized while
between-galaxy means wandered — residual std went UP vs fiducial, the opposite
of the method's purpose). Fixed to RMS (the mean counts); Υd median landed on
the canonical 0.5. The edge-driven fit is recorded here because it is exactly
the class of silent error the program exists to catch: a fit can improve its
own objective while degrading the quantity it feeds.

**Headline**: under honest cuts, literature-standard per-galaxy Υ, and a global
a₀, **88.2% of the galactic anomaly's variance above the instrument floor is
priced by included-set variables alone** — and the LOO variant (0.876) says
none of that is per-galaxy overfitting. The excluded sector (whatever it is:
particle halo, modified inertia, external field) has at most ~12% of the
anomaly's variance-room in which to manifest *beyond* what the included set
already prices, on these instruments. dp's epicycle diagnostic, rung 2: the
anomaly is epicycle-class at the 0.88 level, measured, not analogized.

## SA-1b first cut: the residual does not scale as N^−1/2

Per-galaxy residual std σ_g vs included-DOF proxy M_star = Υd·L[3.6]
(149 galaxies with ≥5 points):

- **log-log slope = −0.010** — flat. Tercile medians: 0.082 / 0.079 / 0.058 dex
  (low/mid/high M*) — a weak trend, nothing like counting statistics.
- Program reading: the remaining residual is **not** the Gaussian fluctuation
  of included degrees of freedom averaging down; it is *structure* — either the
  correction series (the 6-var terms: surface brightness, environment, type) or
  a genuine sector signal. This is the registered expectation of computation 3
  measured in the first domain attempted: exponent ≈ 0, not 1/2 — i.e., the
  RAR residual carries information beyond the CLT skeleton, and SA-1b's first
  cut says the galactic domain is NOT in K2's dissolving class. (The census
  stands open for the other domains; one domain ≠ 1/2 is enough to keep K2
  from firing globally.)

## The audited variance budget (V2, as registered for rung 3)

```
anomaly Var(log g_obs/g_bar)          = 0.0698
residual Var after 1-param function   = 0.0124
instrument floor                      = 0.0041
ABOVE-FLOOR REMAINDER                 = 0.0082   (12% of the anomaly)
```

Rung 3 (registered): recompute the 6-var correction series on THIS point cloud
with THESE cuts, so the budget chains: anomaly = instrument + 1-param function
+ series + remainder. Rung 4: the remainder's structure — correlates with any
included-set variable left, or sector-demanding? **The falsifier is live**: if
the 0.0082 remainder is structureless, the excluded sector is epicycle-class on
current instruments; if it carries included-set-invisible structure, closure
FAILS and the excluded variable earns its rent. Either outcome is a result.

## What this does NOT claim (honesty block)

- E measures the anomaly's *priceability by the included set*, not ontology.
  A halo could exist AND be forced by structure formation to trace the included
  set tightly (feedback tuning) — E = 0.88 says the sector is then epicycle-class
  *on these instruments*: invisible except through included-set correlations.
  It does not say there is no halo.
- ν's functional form is an input (McGaugh). E is relative to that function
  class; a broader class (free spline) would be a different, registered variant.
- The M_star proxy counts stellar-mass units, not independent dynamical degrees
  of freedom; SA-1b's slope is a first cut, not the census.
- No He-uncertainty, distance-error, or inclination-error propagation into the
  instrument floor yet (errV only) — the floor is understated, so E_corr is
  mildly conservative.

## Next

Rung 3 (the chained budget with the 6-var series), then rung 4 (remainder
structure), then program results doc 03. SA-1a/1c/1d and SA-3/4 remain open
for pickup per the sub-arc map.

## Provenance

kimi seat, interactive session with dp, 2026-09-08. Instrument deterministic;
rerun for every number. The variance-budget arithmetic is printed from the same
arrays as the table above — no hand-transcribed numbers.

---

## Erratum 2026-09-08 (codex review; full adjudication in `2026-09-08-kimi-response-to-codex-review.md`)

**The He factor was a double-count** — SPARC's Vgas already includes x1.33
(Lelli+2016 §3.3). The numbers above are superseded by the corrected rerun
(also: residual means now reported; point-LOO renamed as interpolation
stability; SA-1b claim narrowed — slope ~0 does not establish non-Gaussianity
or escape from K2; K2 stands open):

| variant | residual std | mean | E_corr |
|---|---|---|---|
| V1 fiducial | 0.1617 | — | 0.744 |
| V2 fitted (a₀ = 1.05e−10) | 0.1069 | −0.003 | 0.910 |
| V3 point-LOO (interpolation stability only) | 0.1099 | −0.010 | 0.902 |
| **V4 GALAXY holdout (116 train / 37 test, a₀_train = 9.9e−11)** | **0.1000** | **−0.004** | **0.919** |

The corrected headline is STRONGER than the buggy one: **0.919 on never-seen
galaxies**, no degradation out-of-sample. The bias-blindness of E is now
declared with every value (codex T6), and rung 3 gains an identifiability
column (codex T5: the above-floor remainder is a compatible interval, not a
point estimate of sector signal).

---

## Erratum 2, 2026-09-08 (codex V4-leakage flag; repair doc `2026-09-08-kimi-v4-leakage-repair.md`)

V4's Υ was fitted from the held-out galaxy's own rotation curve — the 0.919
above is **test-galaxy rotation-curve calibration**, not prediction. The
strict predictive numbers (zero test information): **V4a fixed prior 0.812,
V4b photometric mapping 0.821** (mapping adds nothing over the prior); the
labeled few-shot variant (half-curve calibration) is 0.912. The predictive
headline of this program is V4a's 0.812.
