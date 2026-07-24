# SPARC × Cassini squeeze: registered result

**Execution date:** 2026-07-23<br>
**Registration:** `PREREGISTRATION.md`, commit `9c77e7be`<br>
**Execution code/profile parent:** `e05e3582`<br>
**Registered outcome:** **A — robust empty intersection**

## Verdict

> The registered scale-universal tanh-log family is incompatible with the
> joint SPARC and Cassini constraints as a QUMOND interpolation family at the
> registered resolution.

This is a family-specific refutation, not an independent prediction and not a
refutation of the Synchronism umbrella ontology. It applies when the same
universal tanh-log interpolation function is used for galaxy rotation and the
Solar System in QUMOND.

## Primary result

The exact 2026-07-22 SPARC likelihood was reconstructed on the 79-point
registered gamma grid. It reproduces the historical checkpoints:

| Check | Reconstructed value |
|---|---:|
| Selected SPARC rows | 2,807 |
| Best grid member | `gamma = 0.489` |
| Free-family Delta BIC vs. McGaugh | `+7.1069` |
| Fixed `gamma = 2` Delta BIC vs. McGaugh | `+184.0445` |

Under the primary profile-likelihood convention, the registered
`Delta BIC <= 10` SPARC grid interval is:

```text
gamma = 0.425 through 0.600
```

None of those points passes the current signed, two-sided 95% Cassini
quadrupole interval. At the SPARC optimum, with its profiled
`a0 = 5.33265e-11 m s^-2`:

```text
gamma = 0.489
Q2 = 3.39043e-26 s^-2
Cassini discrepancy = +17.95 sigma
```

Across the complete union of all SPARC-retained sensitivity sets
(`gamma = 0.400` through `0.625`), the primary-field predictions span only:

```text
Q2 = 3.34743e-26 through 3.40002e-26 s^-2
z  = +17.71 through +18.00
```

Thus the result is not driven by the best-fit point or by a narrow boundary.

## Registered sensitivities

- `g_ext = 2.00e-10 m s^-2`: zero Cassini survivors; discrepancies
  `+14.78` through `+15.33 sigma`.
- `g_ext = 2.48e-10 m s^-2`: zero Cassini survivors; discrepancies
  `+19.18` through `+19.35 sigma`.
- Legacy 2014 Cassini interval at the primary external field: zero survivors.
- SPARC `Delta BIC` thresholds 6, 10, and 14: every joint set is empty.
- All three recorded interpretations of `Delta BIC(gamma)` give an empty
  joint set:
  profile-relative, fixed-gamma versus McGaugh, and the two-parameter family
  versus McGaugh.
- The framework-fixed `gamma = 2` is already excluded by SPARC at
  `Delta BIC = +184.04`; it cannot rescue the joint intersection.
- The direct high-acceleration residual is recorded for every grid row in the
  machine-readable result.
- AQUAL was not independently computed. The literature-supported direction
  suggests it is at least as constrained for benchmark families, but no
  AQUAL-specific verdict is claimed.

## Numerical validation

All literature benchmarks and `mu`-to-`nu` mapping checks pass. Every
Cassini-relevant row was required to meet:

1. two consecutive radial Gauss-Legendre changes below 0.5%, using orders
   from 512 through 8192; and
2. an independent angular-order 128-to-256 change below 0.5%.

The observed worst accepted radial change was `0.2169%`; the worst angular
change was `0.00501%`; and the worst mapping relative error was
`1.38e-10`.

The initial full-domain attempt produced no verdict because low-gamma
quadrature was unstable. `INSTRUMENT_AMENDMENT_2.md` preserves that failure
and the fixed response. Cassini is short-circuited only where SPARC already
excludes a point beyond `Delta BIC > 14` under every convention, so those
numerically delicate values cannot enter any registered primary or
sensitivity intersection.

## Scope and interpretation

The result closes this realization:

```text
universal tanh-log interpolation + modified-gravity QUMOND
```

It does not directly test modified inertia, a non-gravitational engineering
compander, system-dependent or multi-scale functions, dark-matter/hybrid
models, or the umbrella ontology. Escaping the result requires changing the
realization rather than retuning `gamma` inside the registered family.

The SPARC likelihood also inherits the limitations of the source analysis:
fixed stellar mass-to-light ratios, an error cut, and unweighted log-space
residuals. The closure is nevertheless insensitive to the three BIC
conventions and all registered thresholds because the Cassini discrepancy
throughout the retained region is approximately 15–19 sigma.

## Reproducibility

Machine-readable artifacts:

- `sparc_profile.json` — frozen likelihood profile and source/data hashes.
- `joint_result.json` — every SPARC row, evaluated Solar-System rows,
  convergence histories, sensitivities, and verdict.

Commands:

```bash
python3 simulations/sparc_tanhlog_profile.py \
  --output Research/preregistrations/sparc_cassini_tanhlog/sparc_profile.json

python3 simulations/sparc_cassini_joint.py \
  Research/preregistrations/sparc_cassini_tanhlog/sparc_profile.json \
  --output Research/preregistrations/sparc_cassini_tanhlog/joint_result.json
```

The result JSON identifies execution parent `e05e3582`, the registration and
instrument-amendment commits, runtime dependency versions, the original
synchronism-site analysis commit and script hash, and the SPARC data hash.
