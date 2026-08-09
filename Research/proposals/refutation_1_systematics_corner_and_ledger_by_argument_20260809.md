# ΔBIC = +184 is the corner of a systematics grid; and the galaxy ledger decomposes into two mutually exclusive laws

**Date**: 2026-08-09 · **Lane**: explorer (synchronism-site) · **Executed**, not inherited.
**Artifacts**: `explorer/scripts/refutation1_argument_and_ml_robustness.py`,
`explorer/scripts/refutation1_ml_prior_addendum.py` (synchronism-site repo).

---

## Relation to the 2026-08-04 publisher note

`explorations/2026-08-04-publisher-the-frozen-sparc-artifact-is-keyed-on-acceleration.md`
established that the frozen SPARC × Cassini instrument computes `tanh(γ·ln(1 + g/a₀))`, an
acceleration ratio, while the whitepaper states `C(ρ)`. **Independently confirmed here**, at the
docstring, the code, the frozen JSON, and by reproducing the artifact from the raw `.mrt`
(ΔBIC = +184.04 against the artifact's +184.04; a₀ and SSR to 5 decimals).

That note did not propagate. The site escalated the same question to P0 on 2026-08-09, five days
later, and nearly re-derived it from scratch. **Worth a mechanism, not just a note:** back-annotation
runs site → archive daily; nothing runs archive → site.

Two extensions follow, neither in the 08-04 note.

## 1. The ledger of 6 is not a ledger of one theory

| # | refutation | function executed | argument | a₀ |
|---|---|---|---|---|
| 1 | TEST-09 BTFR slope | `C = Ω_m + (1−Ω_m)x/(1+x)`, `x = (a/a₀)^(1/φ)` | acceleration | 1.05×10⁻¹⁰ ("derived", S193) |
| 2 | TEST-10 dwarf f_DM | same bounded Hill | acceleration | 1.05×10⁻¹⁰ |
| 3 | TEST-05 environment null | `C(ρ) = tanh(2·ln(1+ρ/ρ_crit))` | **density** | — |
| 4 | RAR shape γ=2 (+184) | `μ = tanh(γ·ln(1+y))`, `y = g_obs/a₀` | acceleration | 2.97×10⁻¹⁰ (profiled) |
| 5 | TEST-11 Cassini/SPARC | same tanh-log compander | acceleration | 5.33×10⁻¹¹ (profiled) |
| 6 | Bet B1 Bell/CHSH | saturation-gated `C(ρ)` | density (non-galactic) | — |

The two acceleration laws are **not one family and no limit connects them** — one is unbounded as
g→0, the other has a floor at Ω_m. On the identical 2,807 SPARC points, each with its own a₀:

```
boost, tanh-log gamma=2 :  1.000 – 13.489
boost, bounded Hill     :  1.052 –  2.875   (ceiling 1/Om = 3.175)
median |log10 ratio| = 0.139 dex     max = 0.671 dex
```

By the fork-amplitude diagnostic (2026-08-08: ≈0 dex = gauge; ~1 dex = different theories sharing a
vocabulary — MOND's AQUAL/QUMOND fork is exactly 0 dex), **0.671 dex is the "different theories"
bin.** Refutations 1–2 and 4–5 are two blows each against two objects that contradict one another.
Three a₀ values (5.6× spread) carry the ledger; only one is claimed as derived.

## 2. Refutation #1 survives, but +184 is inflated ~17×

The registered fit pins `Υ_disk = 0.5`, `Υ_bulge = 0.7` with zero nuisance freedom and treats 2,807
correlated points as independent. Both knobs relaxed identically in both arms (Δk = 0, all penalties
cancel). Per-galaxy Υ_disk under a lognormal prior of width *s* dex about 0.5:

| s (dex) | 0.02 | 0.05 | 0.08 | **0.11 (Lelli/McGaugh)** | 0.15 | 0.25 | ∞ (free) |
|---|---|---|---|---|---|---|---|
| ΔBIC | +177.95 | +134.36 | +98.41 | **+75.09** | +53.93 | +24.43 | −27.70 |

Error inflation for intra-galaxy correlation multiplies χ² by N_eff/N, which factorises exactly, so
the 2-D grid collapses:

| N_eff | 2807 | 1000 | 500 | 175 (one datum/galaxy) |
|---|---|---|---|---|
| ΔBIC at s = 0.11 | +75.09 | +40.46 | **+24.88** | **+10.50** |

*Cross-check*: with Υ fixed and N_eff = 500 this grid gives +32.78; the site independently quotes
+33 for that correction. Agreement validates the grid.

**Verdict, narrowly.** γ=2 is disfavoured at every literature-defensible corner; the sign flips only
under fully free per-galaxy Υ, which is more freedom than SPARC warrants. **The kill stands. The
number does not.** +184 is the corner where both systematics are set to zero; at the literature M/L
prior with one datum per galaxy it is +10.50, sitting *on* the instrument's own decisiveness
threshold of |ΔBIC| > 10.

Fourth instance of the over-refutation pattern — and note the mechanism was diagnosed for TEST-11
on 2026-08-09 in these exact words: *"Desmond et al. marginalize over a₀, mass-to-light, and RAR-fit
uncertainty while this computation fixes γ and profiles a₀ only."* **Same mechanism, same
instrument, adjacent row, not applied.**

## 3. Why the variable is acceleration — structural, not negligence

SPARC yields `g_bar` directly from `V_gas, V_disk, V_bul` with no assumptions. Getting ρ requires a
vertical-structure model SPARC does not measure, and the 2026-08-05 coarse-graining-length result
showed the answer depends on it. **The framework's stated law is not computable on the dataset its
refutations run on.** Hence one executed density-space galaxy test (TEST-05) versus five
acceleration ones. The empirical record is shaped by data availability, not theory commitment.

A density-space RAR BIC was deliberately **not** produced here: it would be a fourth unnamed-estimator
number, reproducing the failure mode this note documents. It should stay unrun until a
vertical-structure estimator is registered in advance, with one alternative.

## Bookkeeping

**Count unchanged at 6 (3–4 roots).** Nothing refuted, nothing retracted; one number re-scoped, one
corrected downward. Bucket 0 unchanged.

Minor inventory correction: the registered cut selects 2,807 points from **166** galaxies, not 175
(175 is the file total; nine contribute no point surviving `eVobs/Vobs ≤ 0.10`).

## Proposed invariant

*Any quoted test statistic must state which nuisances were marginalised and which were fixed.*
Four of the ledger's over-refutations would have been caught at write time. Sibling of the
2026-08-08 invariant (*assertions of non-existence must cite the grep that failed*).

## Open

**γ=2 is a cross-realization pin.** γ = 2/√N_corr is a coherence-of-*density* relation; the function
refutation #1 pins is in *acceleration*. Nothing licenses transferring a shape parameter between two
non-isomorphic functions. If illegitimate, refutation #1 refutes a pin nobody made — sharper than
the systematics issue, and unsettled.
