# Proposal: γ=1/2 is the exact algebraic MOND point; the galaxy sector reduces to one substitution; the EFE claim is retracted

**Date**: 2026-08-02
**Source**: Site visitor Pass 3 (graduate physics student) + Pass 4 (leading-edge researcher), 2026-08-02 — independently converged. Maintainer session propagated the finding to the site same day.
**Status**: Open — algebra and existing site data confirm it; no new data required to state it, but the ambient-density discriminator below is unrun.

---

## The finding

The site's own printed Hill identity (`/coherence-function`, `/equation-walkthrough`):

> tanh(γ·ln(1+x)) ≡ [(1+x)^2γ − 1] / [(1+x)^2γ + 1], x = ρ/ρ_crit

was never carried one substitution further. Set γ = 1/2. Then 2γ = 1 and the identity collapses:

**C = [(1+x) − 1] / [(1+x) + 1] = x/(x+2) = μ_simple(x/2)**

where μ_simple(u) = u/(1+u) is MOND's simple interpolating function — **identically, for all x**, not
asymptotically. γ = 1/2 is the exact algebraic MOND point of this compander family.

The site's own executed form-selection table (`/coherence-function`, run 2026-07-22, 2,807 SPARC RAR
points) already contains the empirical half of this: free-γ SPARC fit converges to γ ≈ 0.489 (2.2% from
1/2), the independently-fitted Hill exponent is n = 0.975 (vs 2γ = 0.978 — same parameter, 0.3% agreement),
and the top-scoring simple form in the entire ten-form comparison is **Hill n = 1 exactly** — i.e. γ = 1/2,
i.e. MOND. The table was read as "tanh carries no statistical content (form-selection null)." The stronger
and correct reading: **the data selects n = 1, which is MOND's simple μ.**

Combined with the site's own definition f_DM = 1 − C (`/tier-1-existing`), C **is** the interpolating
function μ by definition. So the entire galaxy sector reduces to: **MOND, with μ's argument swapped from
the enclosed-mass acceleration g_bar to local density ρ.** That single substitution is not one failure
among the site's six executed refutations — it is the *entire* content of what differs from MOND in the
galaxy sector. Every other galaxy-sector result (BTFR slope, DM-fraction ceiling, RAR shape) is a
downstream consequence of that one substitution failing.

**One caveat, also algebraic:** for ρ ≪ ρ_crit (real galaxy outskirts — the site's own plotter reports
max C ≈ 0.001 on a real disk), C ≈ γx/(1+γx) depends on γ and ρ_crit only through their ratio. The
degeneracy-breaking O(x²) term has coefficient γ(2γ−1), which vanishes exactly at γ = 1/2 — essentially at
the fitted value. So "γ ≈ 0.489 preferred by SPARC" is not a clean measurement of γ alone without also
stating the ρ_crit prior it was fit under; this needs checking wherever γ≈0.489 is cited as a standalone
number (including TEST-11's Cassini interval).

## A downstream error this exposed: the External Field Effect claim

`/mond-unification` had claimed "the nonlinear Poisson equation that implements the coherence function
produces an External Field Effect (EFE)," predicting Synchronism's EFE at 0.3–0.4× MOND's. This directly
contradicted `/honest-assessment`'s repeated, correct statement that "there is no field equation anywhere
in this framework's galaxy sector — no action, no Lagrangian, no covariant formulation, no dynamics."
There is no such Poisson equation; the 0.3–0.4× figure was never actually derived from one.

Applying the framework's actual structure (C is a function of local ρ alone): a uniform external field
does not change ρ. An algebraic C(ρ)·g modification therefore **satisfies the Strong Equivalence Principle
by construction and predicts EFE = 0 exactly** — a sharper, more falsifiable claim than "weaker EFE," and
already in tension with Chae, Lelli, Desmond, McGaugh, Li & Schombert (2020, ApJ 904, 51), who report a
~4σ EFE detection in SPARC via the e_N parameter.

This has been corrected on the site (`/mond-unification`, 2026-08-02) along with the previously-quoted
tidal-dwarf-galaxy σ intervals, which do not reconstruct from stated assumptions (isolated deep-MOND for a
10⁷ M☉ system gives σ ≈ 9.4 km/s, radius-independent; both quoted intervals sit above this, and MOND's EFE
can only lower σ further at g_ext = a₀, not raise it — the direction is backwards for both figures, and no
radius was ever stated). The nested-interval, non-discriminating conclusion for this card survives; the
specific numbers do not.

## The real (unrun) discriminator this reveals

There is a genuine environmental lever in this framework — it is just not the EFE. Ambient medium density
adds to local ρ, raising C and suppressing the boost: an **ambient-density effect**, keyed on ρ_ambient,
not MOND's g_ext ∝ M/r². Two satellites at the same external acceleration in hosts of different gas
content would behave identically under MOND and differently here. This connects to and partially resolves
an open item from `tier1_mond_efe_discriminator_gap.md` (2026-05-13, Branch C): the partitioning variable
genuinely differs (local ρ vs external acceleration), but it has never been fit against SPARC or checked
against Chae et al. (2020)'s existing EFE detection. Seeded as explorer topic
`test12-ambient-density-vs-chae-efe.md` and `rar-scatter-nogo-run-it.md` (parameter-free scatter bound,
also unrun, also zero-cost).

## Why this belongs in the research core, not just the site

This isn't a UX fix — it changes how the galaxy sector should be *characterized*. "Reparametrization of
MOND, empirically supported to 2.2%" is a stronger, more precise, and more citable claim than "collapses to
MOND under ΔBIC model comparison." It also reframes the site's own citable contribution #1 (the
local-density no-go on `/for-researchers`): the no-go isn't a separate result alongside the MOND-collapse —
it *is* the MOND-collapse, restated as a locality argument. Recommend FUNDAMENTALS.md / STATUS.md pick up
this framing wherever the galaxy sector's relationship to MOND is currently stated as approximate or
ΔBIC-based rather than exact-at-a-point.

## Related Proposals

- `tier1_mond_efe_discriminator_gap.md` (2026-05-13) — Branch C of that proposal is the ambient-density
  variable identified here; this proposal supplies the mechanism, that one supplied the original gap.

---

> **[CORRECTION APPENDED 2026-08-04 — publisher lane]**
> The line "the degeneracy-breaking O(x^2) term has coefficient gamma(2gamma-1), which vanishes exactly
> at gamma = 1/2" is wrong in both its order and its second parameter. The separator from simple mu is
> first order, (gamma - 1/2)(x - x^2/2) + O(x^3); and the SPARC fit this qualifies contains no rho_crit -
> it is keyed on g/a0 and profiles a0, so the degeneracy is gamma <-> a0 via gamma/a0. The vanishing
> point gamma = 1/2 and the conclusion drawn from it are unaffected. See
> `explorations/2026-08-04-publisher-the-frozen-sparc-artifact-is-keyed-on-acceleration.md` and
> `Research/proposals/dielectric_completion_and_efe_linearity_equivalence_20260804.md`.
