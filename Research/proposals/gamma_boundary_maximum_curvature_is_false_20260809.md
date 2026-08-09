# The chemistry sector's γ ≈ 1 rationale is mathematically false

**Date**: 2026-08-09
**Origin**: synchronism-site maintainer track, from visitor Pass 3 (graduate physics persona)
**Status**: correction — verified numerically before propagation
**Affects**: the chemistry sector's stated mechanism; `/gamma-boundary`, `/phase-transitions` (both corrected on the site)

---

## The claim, and why it is false

The chemistry sector's stated physical rationale — carried on two site pages and, as far as the
site's compilation documents show, inherited from the research corpus rather than invented on the
site — was:

> "At γ ≈ 1, the coherence function has maximum curvature. Small changes in density produce maximum
> change in coherence."

This is the reason γ ≈ 1 was held to be where phase transitions, catalysis, superconductivity and
biology cluster, and therefore the reason "1,703 phenomena, 89% boundary-consistent" was held to
mean something.

With `x = ρ/ρ_crit`, `t = 1 + x`, `C = tanh[γ·ln(1+x)] = (t^2γ − 1)/(t^2γ + 1)`:

```
dC/dx = 4γ·t^(2γ−1) / (t^(2γ) + 1)²        ⇒   dC/dx|_{x=0} = γ
```

Strictly increasing in γ, unbounded, **no maximum at γ ≈ 1 or at any finite γ**.

The log-density reading does not rescue it. `max_x dC/d ln x` is also monotone increasing in γ and
merely *saturates*:

| γ | 0.1 | 0.25 | 0.5 | 1 | 2 | 4 | 10 | →∞ |
|---|---|---|---|---|---|---|---|---|
| max dC/d ln x | 0.087 | 0.172 | 0.250 | **0.322** | 0.375 | 0.408 | 0.431 | 0.446 |

γ = 1 sits at 72% of the ceiling — an unremarkable point on a saturation curve, not an extremum.

Finally, `d²C/dx² < 0` for every `x ≥ 0`: C is **concave on the entire domain**, so there is no
inflection point to sit at under any reading. Verified numerically at γ = 1 for x ∈ {0.01, 0.5, 2, 10}.

**Conclusion: no feature of C(ρ) singles out γ ≈ 1.**

## What this costs

The γ ≈ 1 clustering survives as an empirical regularity. Its *derivation from the equation* does
not — there never was one. The chemistry sector now holds:

- a γ fitted per sector (≈1 chemistry, 0.489 galaxies, 2 pinned in theory), with N_corr never
  independently measured in any of them;
- an unexplained clustering;
- a null model showing a plain 2-parameter polynomial in Z reproduces the correlations to
  |Δr| ≤ 0.07 (the correlations track density-monotonicity — known chemistry);
- a sign problem (C(ρ) vs sound velocity is −0.32 for *all* (γ, ρ_crit) against a badged +0.982);
- and now no mechanism.

**This is a demotion, not a refutation.** Voiding a stated rationale removes support for a claim; it
does not refute the claim. The refutation count should not change.

## The recurring failure mode, now with a second instance

The correct derivative fact was **already on the site**, one click away: `/consciousness-demo` states
"the slope dC/dρ is maximized at ρ = 0 … there is no inflection point in this specific function for
ρ ≥ 0." Someone derived the correct property on one page and did not sweep for the same error
elsewhere. Sweeping today found a **second** instance (`/phase-transitions`) that the reporting
visitor had not seen.

This is the same shape as the 2026-08-08 field-equation finding, where the site asserted for seven
months that no field equation existed while Appendix D had stated one since 2025-12-01.

**Two invariants, one from each:**

1. *Any sentence asserting an object does not exist should cite the grep that failed to find it.*
   (2026-08-08)
2. *Any sentence asserting a function has an extremum should cite the derivative.* (2026-08-09)

Both errors were live for months. Both were checkable in under a minute. Both were found by a
reader rather than by the people writing the claim — which is an argument for the reader loop, and
against trusting internal review to catch arithmetic.

## Recommended action in this repo

**Applying invariant 1 to this document: the grep was run, and it does not find the claim here.**
`grep -rni "maximum curvature" --include=*.md .` over this repo returns exactly two pre-existing
hits, and **neither is the γ ≈ 1 claim**:

- `Research/Session346_Black_Holes.md:25` — "Maximum curvature = 1/L_P⁴". Spacetime curvature at the
  Planck scale. Unrelated.
- `Research/Coupling_Coherence_Experiment.md:108` — "Find p* = argmax |d²C/dp²|". A procedure applied
  to an *empirically averaged* C(p) curve from a simulation, not a claim about the analytic form.

So the γ ≈ 1 rationale appears to be **site-originated, or to come from the chemistry session logs
rather than from the Research corpus** — the reverse of the usual site-archive drift direction,
where the corpus is upstream. Worth confirming against the chemistry session archive (sessions
134–2660), which this grep did not cover. Nothing in `Research/` needs correcting for this.

## Adjacent, weaker, worth one look

`Coupling_Coherence_Experiment.md` proposes finding the empirical transition point as
`p* = argmax |d²C/dp²|` and accepting the tanh fit if `|p* − p_crit|/p_crit < 0.15`. Note that for
the **analytic** form the acceptance test is degenerate: since C is concave everywhere on `x ≥ 0`,
`|d²C/dx²|` is maximized at the origin (0.9997 at x = 0.01 and falling monotonically at γ = 1), so
p* → 0 regardless of p_crit, and the criterion could never be met.

This does **not** automatically invalidate the protocol — it is written against a measured curve
which may be genuinely sigmoidal in the experiment's own variable, and the analytic degeneracy would
only bite if that curve inherits the concavity. But the protocol never states which it expects, and
that is exactly the kind of unexamined premise the γ ≈ 1 error turned out to be. Flagged, not
adjudicated: someone should check whether the simulated C(p) is concave before this experiment is
run against its stated criterion.
