# Instrument amendment 1: fixed radial quadrature

**Date:** 2026-07-23<br>
**Applies after registration commit:** `9c77e7be`<br>
**Scientific thresholds changed:** No<br>
**Target-family mapping changed:** No

## Trigger

The first post-registration target-family run reached SciPy's adaptive
subdivision ceiling for `gamma = 0.489` and `0.5`. Under preregistered outcome
branch D, those values could not be used for a scientific result.

The warning arose in the one-dimensional log-radius integration surrounding
the analytically known field-cancellation location:

```text
v^2 = e_N.
```

Splitting the adaptive integral at that point did not remove the warning. The
reported values were stable, but warning suppression by raising the
subdivision ceiling would not provide a transparent convergence argument.

## Fixed change

The angular Gauss-Legendre quadrature and stabilized QUMOND integrand are
unchanged. The radial integral now uses deterministic Gauss-Legendre
quadrature in log-radius, split at `v^2 = e_N`:

```text
primary radial order = 512 per segment
convergence order     = 256 per segment
```

The reported numerical error is the absolute difference between those two
orders. The target row passes only if this difference is below 0.5% of
`|q|`. A doubled-order run at 1024 is retained as a separate sensitivity.

This change was chosen from the geometry of the published integral, not from
which result it produced.

## Revalidation

After the amendment, the three literature benchmarks remain within the
registered 2% tolerance:

| `e_tilde` | Published `|q|` | Computed `|q|` | Relative difference |
|---:|---:|---:|---:|
| 1.0 | 0.094 | 0.0947124440 | 0.758% |
| 1.5 | 0.159 | 0.1594854552 | 0.305% |
| 2.0 | 0.221 | 0.2209973126 | 0.0012% |

For the benchmarks, doubling radial order changes `|q|` by less than
`2.5e-10` relative. Doubling angular order changes it by less than `5.3e-9`
relative. Expanding the log-radius bound from 14 to 16 changes it by less
than `5.8e-13` relative.

For the slowest target member tested, `gamma = 0.489`, the combined
angular-256/radial-1024 result differs from the primary
angular-128/radial-512 result by approximately 0.035%. It therefore clears
the registered 0.5% numerical threshold.
