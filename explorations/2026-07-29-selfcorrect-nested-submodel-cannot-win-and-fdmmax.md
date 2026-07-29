# Self-correction — my 07-28 "cannot win, a priori" was a category error, and my f_DM,max prescription was a K=1 support; both fixed (2026-07-29)

**Status:** `[ACTIVE-MRH]` — gate-fired by a Publisher-track critique + filed proposal
(`nested_submodel_fit_versus_selection.md`) that corrects **my own** 2026-07-28 inscription on two points.
**Verdict: BOTH objections correct, verified by re-execution; corrected my ledger inscription at source. The
galaxy-sector verdict is unchanged (the ceiling binds); what was wrong was the a-priori framing and the
prescribed statistic. Bucket 0 unchanged (0); arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The signal, and why I accept it against my own prior work

A Publisher-track proposal challenged the two moves I made on 07-28: the "nested-submodel reframe" (which I
called "the sharpest statement of the whole door-#1 cage") and my prescription to cite `f_DM,max` in place of the
convention-dependent "69%." Both objections are correct. Per the 07-27 lesson — *over-defending is the mirror of
over-refutation* — and the 07-09 lesson — *don't declare a structural property without computing on it* — I fixed
rather than defended. The maintainer deliberately did **not** touch PREDICTIONS.md (the reframe was inscribed by
the research lane, i.e. me); correcting it is my job.

## Objection 1 — "a nested submodel cannot win" conflates FIT with SELECTION (category error)

Nesting establishes "cannot fit **better**," not "cannot win." A no-extra-parameter restriction of a family
**can** be *selected* over the parent — that is the whole basis of parsimony. Canonical, same-field example:
**ΛCDM is a strictly nested submodel of wCDM** (w=−1 restricts free w); it never fits better and wins anyway on
every evidence criterion. "ΛCDM cannot win because it is nested" is an instant category error, and my sentence
had the same shape.

So the a-priori content of nesting is a **dichotomy, not a verdict**:

| a-priori branch | outcome | adjudicated by |
|---|---|---|
| ceiling *slack* across the data | tie on fit, **win on selection** (favourable — fewer DOF than a free-function parent) | model comparison |
| ceiling *binds* somewhere | **refuted** | the constraint is falsified |

**Which branch obtains is empirical — decided by SPARC and nothing else.** My 07-28 claim ("cannot win, a priori,
no SPARC data required; the result was never in doubt; every execution was a corollary confirming *where* the
ceiling binds") is therefore backwards in its most important part: **the SPARC runs were the decision procedure
that settled the fork, and one branch (slack ⇒ parsimony win) was genuinely favourable.** That is the *strongest*
answer to "why run it at all," and my framing had thrown it away by declaring the outcome foregone. The verdict
is unchanged — the ceiling binds, hard — but the executed work decided which branch obtained; it was not a
corollary of a foregone conclusion.

**A second fork, from my own earlier finding.** The nesting premise requires the galactic law to be a
*restriction with no compensating free parameter*. But my **2026-07-17 φ-audit** found φ is *fitted-then-named*,
not derived. If φ is genuinely free, `C(a)=Ω_m+(1−Ω_m)x/(1+x)`, `x=(a/a₀)^{1/φ}` is a **two-shape-parameter
(a₀, φ) family that happens to be bounded — not a one-parameter restriction at all**, and nesting is not the
reason for anything. I should have connected these two of my own results; the critique did. Either way, "cannot
win" is not what nesting buys.

## Objection 2 — `f_DM,max` is excluded by one galaxy *by construction*

Re-executed `simulations/test10_ceiling_exceedance_curve.py` (imports the TEST-10 loader verbatim; sample
identical). Reproduces the 07-15 run (106/153 = 69% at 1/Ω_m; f_DM,max = 0.927). The exceedance curve:

| ceiling B_max | f_DM cap | exceed |
|---|---|---|
| 1/Ω_m = 3.17 | 0.685 | 106/153 = 69% |
| (Ω_m−Ω_b)/Ω_b = 5.39 | 0.814 | 44/153 = 29% |
| **Ω_m/Ω_b = 6.39** (most permissive) | 0.843 | **28/153 = 18%** |
| f_DM,max = 13.75 (my 07-28 cite) | 0.927 | **0/153 = 0% by construction** |

**A sample maximum is exceeded by nothing — a K=1 support, the single weakest point on the curve.** Prescribing
it as "the robust statistic" swapped a convention-dependence for a one-galaxy dependence: not an improvement,
just a different place for a referee to push. (Interesting sub-result, honestly recorded by the critique against
its *own* prior: the outlier objection to f_DM,max is empirically refuted — drop-the-max moves 13.75→13.14, only
4% on one galaxy, the tail is smooth. So f_DM,max is *not* fragile; it is just structurally the weakest-supported
point. Both the critique's "it'll be fragile" guess and my "it's robust" claim were declared-then-refuted-by-
computation.)

**The statistic that is actually definition-free AND robust:** the exceedance at the *most permissive* candidate
ceiling, full sample behind it. **At Ω_m/Ω_b = 6.39, 28 galaxies (18%) still exceed** — convention-immune (it is
the most generous convention), 28-object support, not a tail of one. Headline, replacing both "69%" and
f_DM,max: *a bounded-boost modified-gravity class with B_max ≤ 6.4 is excluded by 28 SPARC galaxies; with
B_max ≤ 8.4, by 10.* Stronger than "B_max ≲ 14," immune to the convention objection, computable from data in
hand. I inscribed this as the ledger's TEST-10 statistic.

## Corrections applied

- **PREDICTIONS.md (my 07-28 block):** "cannot win, a priori" → "cannot fit better; the a-priori content is a
  dichotomy whose branch is empirical; the SPARC runs settled the fork, one branch was a win." `f_DM,max`
  prescription → the convention-immune 28-galaxy / B_max ≤ 6.4 class bound. Added the φ-freedom fork. Kept the
  verdict (ceiling binds) and the MOND-Shared-class-law companion role, correctly stated.
- **Registered per-galaxy per-definition sweep** (incl. TEST-09's slope) still **unrun** — this diagnostic does
  not discharge it.
- **Count/preprint questions** still gate on dp.
- **Bucket 0 unchanged (0); arc AT REST.**

## So what

Two turns ago I caught the ledger over-refuting the framework; one turn ago I caught myself over-defending its
no-go; today the critique caught me over-claiming a *cage* — and the mechanism is the one its author names: a
structural, self-deprecating reframe is exactly the profile a self-criticism-rewarding culture waves through. My
"the work was never in doubt" *felt* like rigorous humility and was in fact an over-claim — the work decided
which of two branches obtained, and one was a win. Under-claiming has the same failure mode as over-claiming: a
declared property ("cannot win," "f_DM,max is robust") that fails when you compute on it. The fix leaves the
verdict exactly where it was and makes the reasoning true: the framework's galaxy sector cannot *fit better* than
MOND, and on SPARC it *fails*, and the number that says so is 28 galaxies, not one. Bucket 0 stays 0.
