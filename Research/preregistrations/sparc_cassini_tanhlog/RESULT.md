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

### Why the intersection is empty: a tail-shape mismatch, not a tuning failure

Back-annotated 2026-07-30 from the synchronism-site visitor/maintainer feedback
loop (independent physics-persona pass): the compander and McGaugh's simple-nu
approach the Newtonian limit along different functional tails, and the gap
between them widens without bound as `x = g_bar/a0` grows, so no fitted
`gamma` can close it.

The compander saturates as a **power law**: `1 - C = 2(1+x)^-2*gamma`. McGaugh's
nu deviates from Newtonian as an **exponential**: `nu - 1 ~ e^-sqrt(x)`. At
Saturn, `x = g_Saturn / a0 ~= 6.5e-5 / 1.2e-10 ~= 5e5`. At the SPARC-preferred
`gamma ~= 0.489` (`2*gamma ~= 0.978`):

- compander fractional anomaly: `2 * (5e5)^-0.978 ~= 5e-6`
- McGaugh-nu fractional anomaly: `e^-sqrt(5e5) = e^-707 ~= 0` (identically zero
  in any floating-point sense)

Both numbers are consistent with the profiled likelihoods in
`sparc_profile.json` and `joint_result.json`; this section names the mechanism
rather than re-deriving the verdict. The practical content: over the SPARC
range (`1e-2 < x < 1e2`) a power-law tail and an exponential tail are
numerically indistinguishable at current precision, which is why "curve-
equivalent to MOND" reads as true there — but a power-law tail is bounded
below by an inverse-power decay while the exponential tail is not, so *any*
sufficiently strong-field regime separates them, and the Solar System is far
enough into that regime that the gap is over 5 orders of magnitude rather than
a fitting margin. This reframes TEST-11 from "a parameter choice that happened
to fail" to "a structural consequence of choosing a power-law compander for
the deep-field limit" — the same choice that gives the compander a finite
(non-MOND) boost ceiling in the bounded-C formulation is, in the interpolating-
function formulation, what produces this tail mismatch. Both formulations
trace to the same underlying decision, which is worth stating explicitly if
this closure is cited as evidence against tanh-log companders generally
rather than against this one registered realization.

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

---

## CORRECTION (2026-09-07): the tail-shape section names the wrong MOND function

**Source**: maintainer track, back-annotated from the synchronism-site visitor graduate-physics
persona pass, 2026-09-07. **The registered result is unaffected** — the empty intersection,
+17.95σ at the SPARC optimum, and the +17.71σ to +18.00σ span across the ΔBIC ≤ 10 interval all
came from executed computation, not from the explanatory section below. What is corrected is the
*explanation*, and the correction makes the verdict **more** inherited from MOND, not less.

### The error

The section "Why the intersection is empty: a tail-shape mismatch" attributes the exponential
Newtonian return `nu - 1 ~ e^-sqrt(x)` to "McGaugh's simple-nu." That is two different functions
conflated:

| function | form | Newtonian return |
|---|---|---|
| Milgrom **simple μ** (equivalently simple ν, `nu(y) = 1/2 + sqrt(1/4 + 1/y)`) | `mu(x) = x/(1+x)` | **power law**, `1 - mu ~ 1/x` |
| McGaugh **RAR ν** (the function fitted to the observed RAR) | `nu(y) = [1 - e^-sqrt(y)]^-1` | **exponential**, `nu - 1 ~ e^-sqrt(y)` |

The exponential belongs to the RAR function. Simple-ν returns as a power law: expanding
`1/2 + sqrt(1/4 + 1/y)` at large `y` gives `1 + 1/y + O(y^-2)`.

### Why this matters — the section contradicted the result it explains

This same repository, and the site's TEST-25 row, correctly state that **at γ = 1/2 the compander
is Milgrom's simple μ identically** (`C = x/(x+2) = mu_simple(x/2)`, exact for all x, not
asymptotic). The tail-shape section then argued that the compander's power-law tail *distinguishes*
it from "simple-nu." Both cannot be true. Under the correct naming they are consistent and the
statement is stronger:

> The compander at the SPARC-preferred γ sits in the **simple-μ branch of MOND**, and that branch is
> the one planetary ephemerides had already disfavored (Hees et al. 2016; Blanchet & Novak 2011)
> for exactly this reason — a slow power-law high-acceleration return. The tail does not separate
> this framework from MOND; it separates **one branch of MOND from another**, and the framework
> landed on the closed branch.

### Consequences to propagate

1. The claim "this is the framework's one genuinely non-MOND-degenerate piece of physics" — which
   the site carried until 2026-09-07 — is **withdrawn**. The power-law tail is Milgrom's, not ours.
2. It sharpens the already-registered +17.95σ vs 8.7σ reconciliation (Desmond, Hees & Famaey 2024).
   Their 8.7σ is for **RAR-preferred** interpolating functions; the simple-μ branch is worse. The
   factor ~2 is therefore *partly* the marginalization difference already documented and *partly* a
   genuinely different (and more excluded) interpolating-function family. The registered row should
   say which portion is which, or say that it has not been decomposed.
3. `/galaxy-plotter` on the site draws its MOND reference with the simple-ν — the object this test
   excludes. That inconsistency was already flagged on the TEST-25 row; the naming fix removes the
   remaining ambiguity about whether it is the same function as the compander. It is.
