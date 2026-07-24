# Preparatory Solar-System result

The QUMOND instrument is validated, and the family-specific Solar-System side
of the squeeze is now computed. The registered joint verdict is deliberately
withheld because the exact 2026-07-22 SPARC profile and its
`a0(gamma)` surface have not yet been imported with source provenance.

At the fixed demonstration values
`a0 = 1.03e-10 m s^-2` and `g_ext = 2.32e-10 m s^-2`:

| `gamma` | `Q2` (`s^-2`) | Difference from 2026 Cassini mean |
|---:|---:|---:|
| 0.489 | `4.1331e-26` | `+22.07 sigma` |
| 0.500 | `4.1253e-26` | `+22.03 sigma` |
| 1.000 | `2.6708e-26` | `+13.95 sigma` |
| 2.000 | `5.1780e-27` | `+1.99 sigma` |

At those same fixed values, the current two-sided 95% Cassini interval begins
accepting the family at approximately:

```text
gamma >= 2.00543.
```

Varying `g_ext` across the preregistered range moves that preparatory boundary
to `gamma >= 2.07278` at `2.00e-10 m s^-2` and `gamma >= 1.97539` at
`2.48e-10 m s^-2`. Under the older 2014 Cassini constraint, the central-field
boundary would have been only `gamma >= 1.69499`; the 2026 update materially
tightens this family-specific test.

The prior SPARC optimum remains far above the Cassini limit across that full
external-field range:

```text
Q2(gamma=0.489) = 3.39e-26 to 4.49e-26 s^-2.
```

This strongly supports the predicted direction of the squeeze, but it does
not prove an empty joint intersection. The remaining load-bearing step is to
recover the exact completed SPARC `Delta BIC(gamma)` and `a0(gamma)` profile,
then apply the registered set-intersection rule without changing any
threshold.

The machine-readable values are in
[`preparatory_solar_system_result.json`](preparatory_solar_system_result.json).
The post-registration numerical repair and its revalidation are recorded in
[`INSTRUMENT_AMENDMENT_1.md`](INSTRUMENT_AMENDMENT_1.md).
