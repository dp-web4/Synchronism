# SPARC × Cassini squeeze of the tanh-log family

**Registration date:** 2026-07-23<br>
**Registered by:** Codex (OpenAI)<br>
**Repository base:** `cf5cd3f0370c89b0e240659859fcfc6f68f1d241`<br>
**Status:** Registered before evaluating the tanh-log family's Solar-System
quadrupole. Benchmark-function evaluation is permitted before registration
because it validates the numerical instrument rather than testing the target
family.

## 1. Epistemic status

This is a **prospective family-specific constrained reanalysis**, not a blind
test of modified-gravity MOND in general.

The direction of the expected result is already known. Hees et al. (2016)
showed that several interpolation-function families fitting rotation curves
conflict with the Cassini quadrupole. Desmond, Hees, and Famaey (2024)
quantified a RAR–Cassini tension under modern SPARC modeling. Park et al.
(2026) strengthened the Cassini constraint and reported discrepancies of
roughly 3–15 sigma depending on galaxy modeling. Synchronism's own
2026-07-22 analysis predicts that its SPARC-preferred tanh-log member will be
on the disfavored slow-return side.

What remains prospective is narrower and still useful:

> Does the *specific, previously fitted tanh-log family* contain any value of
> its shape parameter that simultaneously survives the frozen SPARC fit
> criterion and the current Cassini quadrupole constraint when mapped
> consistently into QUMOND?

No wording after execution may upgrade this to an independent prediction.

## 2. Scope

The primary test applies to a **scale-universal modified-gravity reading** of
the tanh-log function in QUMOND:

- the same interpolation function applies to galaxy rotation and the Solar
  System;
- `a0` is universal;
- the algebraic relation is the registered approximation for the primary
  SPARC surface; and
- the Solar-System prediction is the QUMOND external-field quadrupole.

QUMOND is primary because its quadrupole is exactly calculable and is
conservative relative to AQUAL for the benchmark families studied in the
literature. AQUAL is a directional sensitivity case, not a separately
implemented primary test.

The result does **not** directly test:

- modified inertia, whose external-field behavior can be different or absent;
- a tanh-log engineering compander with no gravitational interpretation;
- theories with a second acceleration scale or a system-dependent function;
- dark-matter or hybrid models; or
- the Synchronism umbrella ontology.

## 3. Frozen target family and mapping

The target function is registered in the MOND `mu` convention:

```text
mu_gamma(x) = tanh(gamma * ln(1 + x)),
x = g/a0,
gamma > 0.
```

For QUMOND the required interpolation function is `nu_gamma(y)`, where
`y = g_N/a0`. It is not permissible to substitute `y` directly for `x`.
Instead:

```text
y = x * mu_gamma(x)
nu_gamma(y) = x/y = 1/mu_gamma(x).
```

The conversion is performed numerically as a monotone change of variables.
This mapping is frozen. Discovering that the completed SPARC analysis used a
different mapping triggers the mapping-mismatch branch in Section 9; it does
not license silently changing this registration.

The nominal shape domain is:

```text
gamma in [0.20, 5.00].
```

The execution grid must include `gamma = 0.489` (the prior SPARC optimum),
`0.5`, `1.0`, and `2.0`, with maximum spacing 0.025 in `[0.35, 1.25]` and
0.10 elsewhere. Root refinement is permitted only to locate boundaries after
the full fixed grid has been evaluated.

## 4. Solar-System calculation

The primary quadrupole is the QUMOND result used by Milgrom (2009), Hees et
al. (2016), Desmond et al. (2024), and Park et al. (2026):

```text
Q2 = -3 a0^(3/2) q(e_tilde) / (2 sqrt(G M_sun)),
e_tilde = g_ext/a0.
```

The dimensionless `q` is evaluated with the published two-dimensional
integral. The implementation may subtract an angle-independent `nu` term
whose angular integral is exactly zero; Hees et al. recommend this numerical
stabilization.

Primary nuisance values:

```text
g_ext = 2.32e-10 m s^-2
M_sun = 1.98847e30 kg
G     = 6.67430e-11 m^3 kg^-1 s^-2
```

For each `gamma`, the primary `a0` is the value minimizing the frozen SPARC
profile. If the exact profile cannot be reproduced, the execution is
incomplete: a plot using a single assumed `a0` may be reported only as an
instrument demonstration.

The current Cassini likelihood is:

```text
Q2_obs = (1.6 +/- 1.8)e-27 s^-2  (1 sigma)
```

The two-sided untruncated 95% interval is therefore:

```text
[-1.928, 5.128]e-27 s^-2.
```

The historical 2014 constraint `(3 +/- 3)e-27 s^-2` is a mandatory legacy
sensitivity analysis, not the primary measurement.

## 5. Instrument validation fixed before target evaluation

The QUMOND integrator must pass all of these checks:

1. For the RAR function
   `nu_RAR(y) = 1 / (1 - exp(-sqrt(y)))`, reproduce the published magnitudes
   `|q(1)| = 0.094`, `|q(1.5)| = 0.159`, and `|q(2)| = 0.221` to within 2%.
2. Doubling angular quadrature order changes each benchmark by less than
   0.5%.
3. Expanding the log-radius integration bounds changes each benchmark by less
   than 0.5%.
4. The Simple and Standard functions have the correct deep-MOND and Newtonian
   asymptotes.
5. `nu_gamma(y)` reconstructed from the target `mu_gamma(x)` satisfies
   `y * nu_gamma(y) * mu_gamma(y * nu_gamma(y)) = y` to relative error below
   `1e-7` on a log grid.

Failure of any check produces no scientific verdict.

## 6. Frozen SPARC surface

The primary surface is the completed 2026-07-22 2,807-point SPARC
form-selection analysis, not a newly selected subset. It must be imported
with:

- its exact script commit;
- the row-selection rule;
- `Delta BIC_SPARC(gamma)` on the registered grid;
- profiled `a0(gamma)`;
- parameter count; and
- a hash of the result artifact.

The local source data currently hash to:

```text
MassModels_Lelli2016c.mrt
  9108994b12cc401b94a1768beca61c53ec354779385c9c9cc571049f3043244c
SPARC_Lelli2016c.mrt
  5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243
```

The back-annotation and independent verification hash to:

```text
Research/proposals/compander_family_form_selection_executed_20260722.md
  1b62619e957a155c5d25ad70428b02bbdba45e64e0103a7c5b05abe090342a24
explorations/2026-07-22-verify-compander-form-selection-and-register-wz-scope.md
  aeb72523cc6e9300f9013072b0892a2c802e6d30bc3f6c97211298ac5438fbe6
```

Those prose artifacts report `gamma_best = 0.489` and
`Delta BIC(gamma=2) = +184`. They are cross-checks, not substitutes for the
profile artifact.

## 7. Primary acceptable regions

Define:

```text
R_SPARC   = {gamma: Delta BIC_SPARC(gamma) <= 10}
R_CASSINI = {gamma: Q2_pred(gamma, a0_profile(gamma))
                    lies in the current two-sided 95% Cassini interval}
```

The primary joint-survival set is:

```text
R_joint = R_SPARC intersect R_CASSINI.
```

`Delta BIC <= 10` is a model-retention boundary, not proof or confirmation.
The Cassini interval is applied to the signed prediction; plotting magnitudes
must not replace the signed test.

In addition, report the continuous combined discrepancy

```text
z_Cassini(gamma) =
  (Q2_pred(gamma) - 1.6e-27) / 1.8e-27
```

at the SPARC optimum and throughout `R_SPARC`. A complete posterior-predictive
analysis is preferable, but it may not replace the registered set-intersection
result.

## 8. Mandatory sensitivity analyses

Report without changing the primary verdict:

1. `g_ext` at `2.00e-10` and `2.48e-10 m s^-2`;
2. the 2014 Cassini interval;
3. `Delta BIC` thresholds 6 and 14;
4. the framework-fixed `gamma = 2`;
5. the prior SPARC optimum `gamma = 0.489`;
6. numerical integration resolution and interpolation error;
7. the direct high-acceleration residual of slow-return members; and
8. the direction of the AQUAL correction, explicitly labeled as literature
   transfer unless independently computed.

The no-bulge SPARC model in Park et al. (2026) is contextual evidence, not a
post-hoc replacement dataset. If tested, it must be labeled a separate
secondary analysis.

## 9. Outcome branches

### A. Robust empty intersection

If `R_joint` is empty and remains empty under all numerical and `g_ext`
sensitivities, the conclusion is:

> The registered scale-universal tanh-log family is incompatible with the
> joint SPARC and Cassini constraints as a QUMOND interpolation family at the
> registered resolution.

AQUAL is likely at least as constrained, but that statement remains a
literature-supported inference unless computed.

### B. Primary empty, sensitivity overlap

If the primary intersection is empty but any mandatory physical sensitivity
produces overlap, the family is **in tension and sharply narrowed**, not
closed. Report every surviving interval.

### C. Non-empty primary intersection

If `R_joint` is non-empty, the family survives in the reported bounded
interval. No statement that the family is selected or confirmed is allowed
unless it also beats named alternatives on untouched evidence.

### D. Numerical failure

If benchmark, convergence, or mapping checks fail, the result is
**inconclusive due to instrument failure**.

### E. SPARC provenance failure

If the frozen SPARC profile, row selection, `a0(gamma)`, or source-script
identity cannot be recovered, the joint test is **not executed**. A
Solar-System-only scan may be published as preparatory work.

### F. Mapping mismatch

If the completed SPARC surface used `mu_gamma(g_N/a0)` directly rather than
the registered implicit `mu`-to-`nu` conversion, do not splice the two
surfaces. Report the inconsistency, rerun SPARC under the physically
consistent mapping as a newly versioned analysis, and preserve the original
result as historical.

## 10. Reporting requirements

The result artifact must contain:

- code and data commit hashes;
- dependency versions;
- the full registered gamma grid;
- `a0`, `Delta BIC`, signed `q`, signed `Q2`, and Cassini `z` per row;
- the exact joint-survival intervals;
- benchmark and convergence tables;
- all mandatory sensitivities;
- a machine-readable JSON or CSV result; and
- a prose verdict quoting one of the branches above.

No result may be summarized only by the best-fit point.

## 11. Primary sources fixed at registration

- Park, Hees, Famaey, Desmond, and Durakovic (2026),
  [Improved constraints on modified Newtonian gravity from Cassini radio
  tracking data](https://arxiv.org/abs/2602.17884).
- Desmond, Hees, and Famaey (2024),
  [On the tension between the Radial Acceleration Relation and Solar System
  quadrupole in modified gravity MOND](https://arxiv.org/abs/2401.04796).
- Hees et al. (2016),
  [Combined Solar System and rotation curve constraints on
  MOND](https://arxiv.org/abs/1510.01369).
- Hees et al. (2014),
  [Constraints on MOND theory from radio tracking data of the Cassini
  spacecraft](https://arxiv.org/abs/1402.6950).
