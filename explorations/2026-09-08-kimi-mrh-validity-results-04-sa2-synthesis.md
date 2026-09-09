# MRH-validity, results 04 + SA-2 synthesis: the closure study, complete

**From**: kimi-code · **Date**: 2026-09-08 · **Acting on**: sub-arc map rung 4
(final) · **Instruments**: `sa2_rung4_radial_structure.py` (rung 4) + the rung
1–3 chain, all deterministic in `simulations/sparc_real_data/`

## Rung 4 — the within-galaxy remainder is smooth structure, not noise

37 test galaxies, 740 points, same split and train-only a₀ as rung 3
(1.05e−10 m/s²), within-galaxy deviations from per-galaxy means:

```
within-galaxy remainder Var     0.00672 dex²
instrument Var (errV)           0.00362       (above-floor factor 1.85x)
lag-1 autocorrelation           +0.598        ← smooth structure, not point noise
radial slope (per galaxy)       median −0.044 dex, 46% positive — NO universal trend
radius-decile profile           gentle arch: −0.014 inner → +0.012/+0.020 mid → −0.013 outer
SB profile                      mirrors the arch (SB declines with R)
```

The decisive reading: the *white* component of the remainder ≈ (1−0.598)·0.00672
≈ 0.0027 ≈ the instrument floor 0.0036 — **the instrument accounts for the
noiselike part, and the entire above-instrument within-galaxy remainder
(~0.0031 dex²) is smooth radial structure.** No universal radial trend (slopes
scatter symmetrically), but coherent point-to-point wiggles with a mild
inner-down/mid-up/outer-down arch. This is the known baryonic-echo phenomenon
(Renzo's rule / Sancisi's law: features in the baryonic distribution echo in
the rotation curve) — *included-set* structure the one-parameter form cannot
express, not sector signal. Registered follow-up: direct wiggle-to-feature
correlation (the quantitative Renzo test on this cloud).

## SA-2 synthesis — the epicycle diagnostic, final budget

The galactic anomaly on 37 never-seen galaxies (740 points), every stage
audited, every leak codex found repaired:

```
anomaly Var                              0.07994 dex²   (100%)
  1-param function + per-galaxy Υ         87.4%   priced by the included set
  instrument floor (errV)                  4.5%   measured from the data
  6-var galaxy-level series                0.5%   ≈ nothing after Υ (rung 3)
  smooth radial (baryonic-echo) structure  ~3.9%  included-set, unmodeled by the form
  WHITE REMAINDER                          ~3–4%  ≈ instrument-consistent
SECTOR COMPATIBLE INTERVAL                 [0, ~4%]  (codex T5: interval, not detection)
```

Reading for dp's diagnostic, final form: **the galactic anomaly is 87–96%
included-set-priceable (depending on how the smooth structure is counted), the
instrument is accounted, and the maximum room left for an excluded sector to
manifest *beyond* included-set structure on this channel is ~4% of the anomaly's
variance — an interval, not a zero.** MOND's function is that pricing's shape;
CDM's halo is invisible to this channel except through it. Both are
epicycle-class *on this channel* in the precise sense: the excluded variable
carries almost no channel-reachable information beyond its correlation with
the included set. The discriminator that would change closure is a new channel
(direct detection, lensing, CMB-class data) — exactly dp's "un-instrumented"
diagnosis, now with a number on it.

Strict-prediction honesty (V4 family): with ZERO rotation-curve information,
a new galaxy's anomaly is 81% priced (fixed prior) / 82% (photometric mapping —
which adds nothing); the 90%+ numbers require rotation-curve calibration and
are labeled fitting tasks, not predictions.

## What SA-2 changed in the repo's own ledger

1. **The 6-var model ("strongest output", R²=0.945) is mostly the Υ degree of
   freedom** — train LOO R² 0.061 after free Υ; prices 9.5% of a 5%-component.
   Its predictive utility with fixed Υ stands; its information content does
   not. Flagged for the next STATUS.md stewardship pass.
2. **Υ ≢ photometrically predictable** (mapping train R² 0.046) — the
   mass-normalization structure lives in the rotation curve itself.
3. **The RAR residual's within-galaxy half is baryonic-echo-shaped** — the
   Renzo test is now the sharpest registered follow-up in the lane.
4. Session-484's NN-autocorrelation elimination: retest after free Υ
   (registered, unrun).

## The lane's closing state

SA-2's registered rungs are all executed: day-zero (E=0.74) → rung 2 honest
cuts (0.88, then codex's He repair) → V4 family (strict 0.81) → rung 3 chain
(6-var repriced) → rung 4 (remainder = instrument + echo). The closure study
is complete; the instruments are deterministic and rerunnable; every number
quoted in the series traces to a printed line. Open follow-ups (Renzo wiggle
test, NN-after-Υ, a second channel) are registered but NOT the arc's to hoard
— any seat may take them with a dated doc and a falsifier, per house rules.

The program's other lanes stand: SA-1 census (phonon/quantum/chemistry, open),
SA-3A done / SA-3B-C (codex, active), SA-4 epicycles-quantitative (open, ideal
first pickup), SA-5 web4 bridge (after the physics), SA-6 parked.

## Provenance

kimi seat, interactive with dp, 2026-09-08. One day, one program: charter →
results 01–04 → sub-arc map → three codex adjudications accepted (7+7+2
points), each leaving the claims narrower and truer. This doc closes SA-2.
