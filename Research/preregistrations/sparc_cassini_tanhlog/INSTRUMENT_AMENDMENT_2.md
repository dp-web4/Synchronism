# Instrument amendment 2: convergence-aware set short-circuit

**Date:** 2026-07-23<br>
**Applies after code/profile commit:** `592479d6`<br>
**Scientific thresholds changed:** No<br>
**Target-family mapping changed:** No

## Trigger

The first full-grid execution correctly returned registered branch D rather
than a scientific verdict. Fixed radial order 512, validated in amendment 1
at `gamma = 0.489`, is not converged at the bottom of the nominal domain
(`gamma = 0.20` through approximately `0.40`). The apparent Cassini survivor
at `gamma = 0.25`, `g_ext = 2.00e-10 m s^-2` changed non-monotonically with
quadrature order and was numerical, not physical.

This failure was discovered before any joint result was accepted.

## Pre-existing logical exclusion

The frozen SPARC surface independently excludes every `gamma < 0.40` by more
than the widest registered sensitivity threshold (`Delta BIC > 14`) under
all three recorded BIC conventions:

- profile likelihood relative to the best tanh-log grid row;
- fixed-gamma versus the one-parameter McGaugh reference; and
- the two-parameter tanh-log family versus McGaugh.

Such a point cannot enter any primary or mandatory-sensitivity joint set,
regardless of its Cassini value. Computing a numerically delicate quadrupole
there cannot affect the registered set-intersection verdict.

## Fixed change

1. Evaluate Cassini on the union of every SPARC-retained grid set across
   `Delta BIC` thresholds 6, 10, and 14 and all three BIC conventions.
2. Mark every other row as analytically short-circuited by SPARC, preserving
   its complete SPARC profile and explaining why no Cassini value is needed.
3. For each evaluated row and each registered external-field value, increase
   radial Gauss-Legendre order through `512, 1024, 2048, 4096, 8192` until
   two consecutive changes are below 0.5%.
4. At the accepted radial order, independently require the change from
   angular order 128 to 256 to be below 0.5%.
5. Use SciPy's high-order Legendre-node generator for radial quadrature.
   This is the same deterministic quadrature rule; it replaces NumPy's
   prohibitively slow high-order node construction.

The benchmark functions, QUMOND integral, `mu`-to-`nu` mapping, Cassini
likelihoods, SPARC likelihood, parameter domain, and all verdict thresholds
remain unchanged.

## Epistemic consequence

This amendment can only turn an irrelevant numerical failure into a
well-founded SPARC exclusion. It cannot create a survivor or remove a point
from any SPARC-retained set. If any retained point fails the adaptive radial
or angular checks, branch D still applies and no scientific verdict is
permitted.
