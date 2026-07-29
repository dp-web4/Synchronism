# The Nested-Submodel Reframe Conflates Fit With Selection — "It Cannot Win" Is Not What Nesting Buys

**Filed**: 2026-07-29 (Publisher track, autonomous)
**Status**: proposal — one correction with a pre-registered falsifier, one executed diagnostic
**Bears on**: `boost_ceiling_provenance_and_class_exclusion.md` (Branch 3), the 2026-07-28 triage,
`PREDICTIONS.md` (TEST-09/TEST-10 rows), `SESSION_FOCUS.md`, `stable_fixed_point_preprint_strategy.md`
**Raised as**: a correction to a *framing*, not to any executed result. Bucket 0 unchanged (0).

---

## What was inscribed on 2026-07-28

The nested-submodel reframe, adopted as "the sharpest statement of the whole door-#1 cage":

> Synchronism (galactic) = MOND ∩ {B ≤ 1/Ω_m} — a strictly nested submodel. A nested submodel has
> exactly two possible outcomes against its parent: statistically indistinguishable, or refuted.
> **It cannot win** — a priori, from the model's structure, no SPARC data required. The best
> conceivable result of the entire galaxy program was "MOND with one fewer degree of freedom."

## The correction

The nesting fact is right. The inference drawn from it is wrong, and wrong in a specific,
standard way: **a nested submodel cannot achieve a better *fit*. It can absolutely be *selected*.**
That is the entire basis of parsimony, and it is the difference between maximized likelihood and
model comparison.

The canonical example is in the same field as the claim. ΛCDM is a strictly nested submodel of
wCDM (w = −1 is a restriction of free w). It cannot fit better — it *never* does. It wins anyway,
on every evidence-based criterion, because the restriction costs nothing and buys a sharper prior
predictive. Nobody says "ΛCDM cannot win because it is nested." The sentence would be recognised
as a category error.

The boost ceiling is exactly that kind of restriction: **B_max = 1/Ω_m is fixed by cosmology**, not
fitted to rotation curves. It adds no free parameter. So if the ceiling had been slack across
SPARC, the galaxy sector would have delivered MOND-equivalent fits with strictly less freedom in
the interpolation function — and against a "MOND-class" parent that is a *free function*, that is
not a consolation prize. It is a result.

So the a-priori content of nesting is a **dichotomy**, not a verdict:

| | outcome | how it is adjudicated |
|---|---|---|
| ceiling slack across the data | tie on fit, **win on selection** | model comparison — a real, favourable result |
| ceiling binds somewhere | **refuted** | the constraint is falsified |

Which branch obtains is an **empirical** question about whether the ceiling binds. That question is
answered by SPARC and by nothing else.

## What this changes

**It does not change the verdict.** The ceiling binds — hard, and this run re-verified it (below).
Branch 2 is the branch we are on. TEST-09 and TEST-10 stand exactly as executed.

**It changes what the executed work was for.** The 07-28 triage concludes that the galaxy program's
result "was never in doubt" and that every execution was "a corollary confirming *where* the ceiling
binds." Under the corrected reading, the SPARC runs were **the decision procedure between two
genuinely different outcomes**, one of which was favourable. That is not a corollary. It is the
measurement that settled the fork — and it is the strongest available answer to "why run it at all."

**It changes the preprint's headline sentence.** "A nested submodel cannot win" invites an immediate
referee objection from anyone who does model comparison for a living, and it invites it in the
paper's opening move. The defensible form loses nothing:

> A bounded-boost restriction of a MOND-class interpolation function cannot improve on its parent's
> fit; it can only tie with fewer degrees of freedom, or fail where the ceiling binds. On SPARC it
> fails: [N] galaxies require boosts above [B].

That sentence is a priori where it is a priori, empirical where it is empirical, and states the
result as a class constraint rather than an obituary.

## The dependency the argument has, which cuts both ways

The nesting claim requires the framework's galactic law to be a *restriction* of the parent family
with no compensating free parameter. `C(a) = Ω_m + (1−Ω_m)·x/(1+x)`, `x = (a/a₀)^(1/φ)` has Ω_m
fixed by cosmology, a₀ free (as MOND's is), and **φ — which this archive's own provenance audit
(2026-07-17) found is fitted-then-named, not derived** (S45 ruled B≈φ not significant; S170
reintroduced it with no derivation; S218 concedes the Boltzmann route "gives exponent 1, not 1/φ").

If φ is genuinely free, the galactic law is a **two-shape-parameter family (a₀, φ) that happens to
be bounded** — not a restriction of a one-parameter parent at all. Then:

- the nesting premise fails, so nesting cannot be the reason for anything;
- the parsimony win was never on the table either — but for the *opposite* reason from the one the
  reframe gives (too many parameters, not too few).

Either way, **"it cannot win" is not established by nesting.** Under fixed φ the nesting is real
and the conclusion is wrong; under free φ the conclusion may hold and the nesting is not why. This
fork is worth stating in the preprint explicitly, because a referee will find it.

## Pre-registered falsifier

This proposal is **wrong**, and should be recorded as wrong rather than quietly dropped, if either:

1. **The project's operative comparison is fit-only.** If "0 tests could select Synchronism over
   MOND+EFE+ΛCDM" is adjudicated purely on residual scatter / χ² with no complexity term anywhere,
   then "cannot win" is true under that convention and my objection is terminological. *Current
   evidence against this: the claim uses the word "select," and the archive's own model comparisons
   (e.g. S559: "AIC/BIC both favor gradient (−391, −349)") apply complexity-penalised selection.*
2. **A no-extra-parameter restriction fitting equally well is shown not to be preferred** under some
   criterion the program adopts and names in advance. I know of none and expect none; if one is
   produced, the correction collapses.

Note what is *not* a falsifier: showing that the ceiling binds. That is Branch 2, it is already
established, and it is compatible with everything above.

---

## Executed diagnostic — the headline statistic, and where my own guess was wrong

Script: `simulations/test10_ceiling_exceedance_curve.py` (imports the executed TEST-10 loader and
cuts verbatim; sample identical). **This is not the registered ceiling-definition sweep** — that
sweep also recomputes TEST-09's BTFR slope per definition under a pre-fixed verdict rule and
**remains unrun**. Reproduction check first: 106/153 = 69% exceed 1−Ω_m, f_DM,max = 0.927. Matches
the 2026-07-15 run exactly.

The 07-28 triage prescribed replacing the convention-dependent "69%" with f_DM,max = 0.927
(⇒ B ≥ 13.7) as the definition-free statistic. **I expected f_DM,max to be a fragile single-galaxy
extreme and said so before running it. It is not — that guess was wrong:**

```
Drop-the-max: B_req,max = 13.75 -> 2nd-largest = 13.14
              (only 4% of the bound rests on one galaxy)
```

The tail is smooth, not spiky. The outlier objection is refuted and is recorded here as refuted.

**But a different problem with f_DM,max is real, and it is structural rather than statistical:**

```
definition                     B_max   f_DM cap    exceed
1/Om   (site's choice)          3.17      0.685    106/153   69%
(Om-Ob)/Ob                      5.39      0.814     44/153   29%
Om/Ob  (baryon budget)          6.39      0.843     28/153   18%
f_DM,max (the prescribed cite) 13.75      0.927       0/153    0%
```

**f_DM,max is excluded by exactly one galaxy — by construction.** Nothing can exceed the sample
maximum. The prescribed headline is the single weakest-supported point on the entire curve, which
is the opposite of what "robust" was meant to buy. Trading a convention-dependence for a K=1 support
is not an improvement; it just moves where a referee pushes.

**The statistic that is actually both definition-free and robust** is the exceedance count at the
*most permissive* candidate ceiling — worst-case over conventions, with the full sample behind it:

```
K galaxies of support  ->  largest B_max excluded
    K =   1                 B_max <= 13.75      <- the prescribed cite
    K =   5                 B_max <=  9.46
    K =  10                 B_max <=  8.39
    K =  20                 B_max <=  7.08
    K =  40                 B_max <=  5.69
    K =  76                 B_max <=  4.13
```

Required-boost quantiles (bootstrap 95% CI, 10k resamples, seed 20260729):
p50 = 4.08 [3.83, 4.56] · p90 = 7.59 [6.74, 8.39] · p95 = 8.84 [7.62, 9.51].

### The consequence, and it favours the result

The triage says that under the Ω_m/Ω_b reading "the *median* f_DM = 0.755 passes" and therefore "the
kill stands on the tail, not the median." The first clause is right; the second undersells it by an
order of magnitude in support. **Under the most permissive candidate ceiling, 28 galaxies — 18% of
the sample — still exceed it.** That is not a tail of one. It is 28 independent objects, and it
survives every convention on the table because it is evaluated at the most generous one.

**Recommended headline for the class-exclusion null, replacing both "69%" and f_DM,max:**

> A bounded-boost modified-gravity class with ceiling B_max ≤ 6.4 is excluded by 28 SPARC galaxies;
> with B_max ≤ 8.4, by 10. The exclusion holds under every candidate normalisation of the ceiling
> (1/Ω_m = 3.17, (Ω_m−Ω_b)/Ω_b = 5.39, Ω_m/Ω_b = 6.39), since it is evaluated at the most permissive.

This is stronger than "B_max ≲ 14," better supported than f_DM,max, immune to the convention
objection the triage correctly raised, and computable from data already in hand. It is what makes
this null the cheapest item to ship.

## Disposition

- Correction and diagnostic filed here; **not** inscribed into `PREDICTIONS.md` unilaterally — the
  reframe it corrects was inscribed by the research lane and the count questions already gate on dp.
- `boost_ceiling_provenance_and_class_exclusion.md` — executed-diagnostic note appended.
- `stable_fixed_point_preprint_strategy.md` — Gate Update appended (this is a candidate fourth null
  and its headline sentence moves).
- Registered sweep (per-galaxy, per-definition, including TEST-09's slope) still **unrun**. This
  diagnostic does not discharge it.
- Bucket 0 unchanged (0). Arc AT REST.

## So what

Yesterday I flagged that accepting a retraction uncritically is deference with the sign flipped.
Today's reframe is the same shape one turn further in: it is elegant, structural, and
*self-deprecating*, which is exactly the profile that gets waved through in a project whose culture
rewards self-criticism. It says the work was never in doubt. The work was in doubt — it decided
which of two branches obtained, and one of them was a win. Under-claiming has the same failure mode
as over-claiming; it just doesn't feel like one.
