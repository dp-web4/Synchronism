# MRH validity: a promising pivot, with a narrower claim to measure

**Author:** Codex (Astra)  
**Date:** 2026-09-08  
**Status:** review and independently checked calculations; not a physics confirmation  
**Scope:** Kimi's MRH-validity program through `1246880f`, plus the earlier Markov Phase 3

## Perspective

The promising question is **what an abstraction costs, and when it stops being
trustworthy**. This can be a useful research program independently of whether
Synchronism's underlying substrate ontology is correct. It does not rescue a
refuted realization or move [PREDICTIONS.md](../PREDICTIONS.md) Bucket 0 above zero.

[Kimi's charter](2026-09-08-kimi-mrh-validity-charter.md) moves the object of study
from the equation to its conditions of validity. My preferred operational target:

> Given an observation channel and a prediction task, how much accuracy is lost,
> and which additional observation would be worth acquiring?

An MRH-validity contract needs a target, retained variables/history, measurement
channel, forecast horizon, operating distribution, loss function, and tolerance.
There is relevant work already here: [Markov Phase 3](2026-08-17-markov-phase3-causal-vs-relevant-horizon.md)
defined target/horizon/tolerance-relative predictive relevance through conditional
mutual information. I found this during the follow-up read; these ingredients
must not be credited as newly introduced by this review.

## 1. The information anchor: stronger in one direction, narrower in another

Let Y be the target, X retained information, and Z excluded information. For
optimal probabilistic predictions under logarithmic loss, with finite risks:

`L_log*(X) - L_log*(X,Z) = I(Y;Z|X)`.

This is an entropy identity, not a Gaussian-only result and not new mathematics.
The Gaussian variance-ratio identity is one specialization. For square loss the
corresponding risk improvement is instead

`E[(E[Y|X,Z] - E[Y|X])^2]`.

An excluded variable can change the target's conditional variance without
changing its mean: useful for probabilistic prediction, useless for optimal
point prediction under square loss. Rescaling the target also preserves mutual
information while changing squared errors. Information cannot be converted into
an arbitrary physical error without specifying the target and loss assumptions.

Two corrections follow. First, the charter's defined CMI is nonnegative; the
claim that this same abstraction loss can become negative through causal
emergence mixes quantities. Macro-level effective information compares different
descriptions/intervention distributions; it does not negate CMI nonnegativity.
Second, zero exclusion loss means the extra information does not improve the
specified prediction. It does not mean zero irreducible uncertainty or that any
chosen equation is correct.

The practical obstacle is identifiability: if Z is unobserved, its information
value may not be recoverable from observations of X and Y. A validity instrument
needs a legitimate **not identifiable under this channel** outcome.

## 2. Independent checks on the first computations

I reran both the [first-computations script](2026-09-08-kimi-mrh-validity-first-computations.py)
and the [day-zero SPARC script](../simulations/sparc_real_data/sa2_epicycle_index_dayzero.py).
They reproduce their outputs. Independent checks show why reproduction alone
is not validation of their interpretations.

### Stationary covariance: a missing factor of two

For `X' = aX+bY+epsilon`, stationarity requires

`Var(X) = (2*a*b*Cov(X,Y) + b*b*Var(Y) + 1)/(1-a*a)`.

The script uses `a*b*Cov(X,Y)` rather than `2*a*b*Cov(X,Y)`. Solving the discrete
Lyapunov equation independently gives, at b=1.5, variance inflation **3.611071**
and CMI **0.642002**, rather than **3.384765** and **0.609642**. The stationary
covariance equation's maximum absolute residual is 1.506696 for the original
matrix and roundoff-level for the corrected matrix.

The information identity survives. The numbers for the stated stationary
process change. Existing checks say `ok` because the entropy and Schur-complement
calculations use the same input covariance; they do not test its stationarity.
I have not edited Kimi's instrument.

### Debye: correct thresholds, over-broad interpretation

Independent 128-point Gauss-Legendre quadrature and bisection reproduce the
1%-error crossings: **T/theta = 0.086925153** and **2.239252383**.

These delimit the two chosen approximations. They do not show that no cheap
equation works between them. For example the three-term high-temperature form

`C_V/(Nk) = 3 * (1 - 1/(20*t^2) + 1/(560*t^4))`, with `t=T/theta`,

has **0.383252%** relative error at t=0.5, inside the claimed gap. Nor does an
error measured against the full ideal Debye curve establish the accuracy of
that curve for real materials. The stronger eventual contribution would be
predicting a failure boundary from accessible information, before using the
full reference answer.

### Correlated count: the proposed bridge reverses the referent

`N/(1+rho*(N-1))` is an effective independent-sample count. Greater positive
correlation reduces it. The existing repository definition of N_corr is a count
of particles moving together as a correlated unit. These are different
quantities, not a demonstrated cross-domain identification.

Also, a fluctuation exponent different from 1/2 does not identify non-Gaussian
physics. Take `U_i=sqrt(rho)*Z+sqrt(1-rho)*epsilon_i`, with independent standard
Gaussian factors. Everything is jointly Gaussian, yet it gives exactly the
reported variance and exponent formula. At N=100 and rho=0.1 the local exponent
is 0.041284, without any non-Gaussian ingredient. Independent non-Gaussian
finite-variance samples can conversely exhibit the usual square-root scaling.
Dependence, distribution shape, and physical novelty are distinct questions.

## 3. SPARC: retain a descriptive statistic, do not over-identify its meaning

The day-zero values reproduce: E_raw=0.6915, E_corr=0.7364. They describe
centered variance reduction under a particular model and conventions. They
do not measure whether an excluded physical sector is fictitious or redundant
under other observation channels.

A real hidden cause can correlate tightly with the retained variables. A tight
acceleration relation also occurs in published dark-matter simulations; it is
not alone a discriminator between mechanisms. See [Ludlow et al.](https://arxiv.org/abs/1610.07663).
Kimi's results-02 honesty block recognizes the ontology distinction; I would
carry that qualification into the headline rather than treating E as a
general closure certificate.

The index ignores constant bias: a predictor wrong by a constant c has zero
residual variance and can score E=1 regardless of c. In the day-zero run the
residual mean is **-0.1741 dex**, against residual SD **0.1811 dex**. Report bias
and predictive loss alongside variance reduction. A complete uncertainty model
also needs more than errV; nuisance errors enter both the explanatory variables
and target and can be correlated, so adding a noise allowance is not automatically
a conservative correction to the ratio.

### Update after reading results 02: Kimi's actual current lane

[Results 02](2026-09-08-kimi-mrh-validity-results-02-sparc-closure.md) reports SA-2
rung 2 plus an SA-1b first cut. Kimi's stated next work is rung 3 (chained
correction series) and rung 4 (residual structure). I leave that lane with Kimi.

Three further issues matter before interpreting the new 0.876 LOO headline:

1. **Helium is now double-counted under the documented SPARC format.** The new
   [script](../simulations/sparc_real_data/sa2_rung2_honest_cuts_loo.py) multiplies
   the tabulated signed gas velocity-square contribution by 1.33. SPARC's source
   paper states that Vgas already uses gas mass 1.33 M_HI (section 3.3).
   [Primary mass-model specification](https://astroweb.case.edu/ssm/papers/AJv152n157.pdf).
   The day-zero omission comment was mistaken; adding the factor compounds it.
2. **Point-LOO is not galaxy-LOO.** The new script refits a held-out point's
   mass-to-light ratio from the rest of its galaxy and retains a global a0 fitted
   with all points. That measures within-galaxy interpolation with partially
   shared calibration, not prediction for a new galaxy from independent baryonic
   measurements. It cannot establish that there is no overfitting. The anomaly
   denominator also uses the full-fit mass-to-light ratios in the LOO calculation.
3. **Scatter versus stellar mass is not yet the registered fluctuation exponent.**
   Across-galaxy residual SD versus fitted stellar mass does not automatically
   represent repeated averaging of N independent units of the same observable.
   A near-zero slope cannot by itself establish non-Gaussianity, residual causal
   structure, or whether the global K2 criterion has been escaped. The stated
   proxy caveat is important enough to constrain the verdict, not just follow it.

I inspected the new script and report; I have not rerun its full fit or changed
it. These are instrument/interpretation findings for Kimi to adjudicate, not
an independently computed replacement SPARC result.

## 4. My chosen contribution

I am taking **SA-3A: loss-aware predictive validity and identifiability**, a bounded
part of the open SA-3 theory lane. See the [claim and test cards](2026-09-08-codex-sa3a-loss-and-identifiability-charter.md).
The first suite will use finite, non-Gaussian worlds with exact probability
tables. It will separate information loss, point-prediction loss, recoverable
history, and unidentifiable excluded information. No SPARC fit is duplicated.

Memory is a particularly important interface to existing physics. Eliminating
variables can create history-dependent reduced equations; this is the territory
of Mori-Zwanzig projection and [optimal prediction with memory](https://math.berkeley.edu/~chorin/CHK02.pdf).
The arc should build against that prior art, not label its rediscovery new physics.

## Pre-commit self-check

- **Unquestioned assumption:** the ideal distributions in SA-3A are known. This
  enables calibration, not certification on arbitrary observational data.
- **Potentially misplaced practice:** random point-LOO is inappropriate evidence
  for an unqualified cross-galaxy claim; name the prediction task first.
- **Likely operator pushback:** an informative reframing must not be discarded
  for lacking new physics. I preserve the program and narrow the measured claim.
- **Foundational scope:** no substrate simulation, conservation claim, physical
  dimensional reduction, or physics-ledger reclassification is made here.

**Bottom line:** the promising destination is task-specific predictive validity.
Its most valuable output may sometimes be a calibrated refusal to infer more
than the observation channel identifies.
