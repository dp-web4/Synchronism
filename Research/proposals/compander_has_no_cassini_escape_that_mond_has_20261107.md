# The compander has no Cassini escape that MOND's δ-family has (2026-11-07)

*Site maintainer, from site explorer 2026-11-06 (PREREG site `a4d99f4`) and a maintainer run today (PREREG site `d7eb6b2`).
No bucket moves. Count 6, range 4–8, Bucket 0 = 0.*

## What was run

The explorer found on 11-06 that MOND's "δ ≥ 4 passes Cassini" held a₀ at McGaugh's value. With a₀ refit to SPARC at each δ,
δ = 4 fails at +9.2σ and no δ ≤ 6 passes (Υ 0.5/0.7, bulges in). Dropping the 32 bulge-hosting galaxies brings δ = 6 to +3.0σ,
at a SPARC cost of ΔBIC +164. That is Desmond, Hees & Famaey 2024's mitigation direction (8.7σ fiducial, 1.9σ with extra Υ
freedom and no bulges). The explorer named the compander's no-bulge arm as the open thread. If the compander also fell toward a
pass, the 11-06 range change (3–8 → 4–8, TEST-25 into the rule-invariant floor) would depend on the configuration.

Run today on the same loader, likelihood and Q₂ instrument, with `TanhLogNu(γ)` from `simulations/sparc_cassini_q2.py`.
Script: `synchronism-site/maintainer/scripts/compander_cassini_nobulge_profiled_a0.py` (+ `_output.txt`). Predictions
were in the header and committed before running, and 3/3 held. The identity control reproduces +17.95σ at γ = 0.489, and the
profiled a₀ matches `sparc_profile.json`.

| γ | z@2.32, bulges in | z@2.32, no bulge | no-bulge ΔBIC vs γ = 0.489 |
|---|---|---|---|
| 0.489 | 17.95 | 17.53 | 0 |
| 1.0 | 17.53 | 16.36 | +17.8 |
| 2.0 | 16.26 | 14.44 | +47.9 |
| 3.0 | 15.54 | 13.44 | +61.5 |

Rows at γ ≥ 0.489 converge: changing a₀ by ×1.001 or the solver order (128/512 → 192/768) moves z by ≤ 0.02. **The γ = 0.3 rows
do not converge.** At g_ext/a₀ ≈ 9 the solver returns 8–22σ depending on the a₀ digits and the order. Every draw is > 5, so no verdict
changes. But the 11-06 table's γ = 0.3 entry (17.73) carries about ±6σ of numerical uncertainty, and the 11-06 block's "14.5–19.3σ"
should be read with that row's caveat.

## What it means

1. **The range stays 4–8.** TEST-25's in-house empty intersection survives the bulge cut, at 13.4–17.5σ for γ 0.489–3.
2. **The inheritance is asymmetric.** On the same cut MOND's family falls 2.9σ at δ = 6, and the compander falls ≤ 2.1σ. δ
   sharpens the transition, while γ above ~1 mainly moves the knee: the 10–90 width saturates near 1.17 dex (`gamma-is-knee-position`,
   10-23 block). So the compander has no sharp-transition member to reach for. The 09-27 nesting reading predicted this
   shape: the acceleration branch is a subset of the MOND class, so it can lose alone. Here the MOND class may have an escape
   through members the compander does not contain. "Inherited from MOND" stays literally true at δ = 1 and γ ≈ ½. It is not a
   statement that the framework fares as well as the MOND family once mitigations are applied.
3. **Untested for both families: Υ free per galaxy**, which is Desmond+'s other mitigation. Until it runs, "MOND fails Cassini like the
   compander" must not be written in either direction (the explorer's symmetry warning).

## For the project lead

- Nothing gates on you that was not already gated. If the count rule is ruled on (10-12/10-13/10-18), TEST-25's floor membership is
  now robust to the bulge configuration.
- The site now carries a per-branch equation-of-motion table (algebraic / QUMOND / L2 / L3 / cosmology). Each verdict names its
  dynamics. That is the operational form of the 10-08 "which dynamics is the theory" question, and a researcher persona asked for it as a
  HIGH today. It does not decide the question.
