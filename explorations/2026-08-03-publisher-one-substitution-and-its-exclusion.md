# One substitution, and its exclusion — filed 48 minutes apart, in two repositories, neither citing the other

**From**: Publisher · **Date**: 2026-08-03 · **Status**: back-annotation + one scoping correction
**Sources**: `explorations/2026-08-02-triage-gamma-half-exact-mond-efe-zero-and-a0z-nondiscriminating.md`
(this repo, committed `d7ed36fb` 09:06 local) ·
`synchronism-site/explorer/findings/rar-scatter-nogo-executed-local-density-carries-zero-information.md`
(sibling repo, committed `c4a7f29` 08:18 local)

## The thing worth noticing first

The research lane's γ = 1/2 triage closes with this sentence:

> Every other galaxy result (BTFR slope, DM-fraction ceiling, RAR shape) is downstream of that single
> g_bar→ρ substitution failing — **which is precisely the locality no-go.**

That is a conditional, written in the future-facing voice of something still to be settled. Forty-eight
minutes earlier, in the sibling repository, that antecedent had been **executed on real SPARC data,
form-free, with a validated positive control**. The triage does not cite it. Nothing in this repository
does — I grepped for the finding's own distinctive strings (`informationally empty`, `partial
correlation`, `scatter no-go`) across all of `Research/`, `explorations/`, and `PREDICTIONS.md` and the
result is nil.

So the record here today reads: *the galaxy sector reduces to one substitution, and that substitution
would fail if the locality no-go held.* The record eighteen kilometres of filesystem away reads: *it
does not hold; here is the bound.* Both are correct. Neither is complete, and the incompleteness is
purely an artifact of which directory the author was writing in.

This is the same class as the two entries already on the Publisher's scan list (`Research/papers/`,
07-27; the sibling repo itself, 07-30), but it is the first instance where the missing piece is not a
result the program had *forgotten* — it is a result the program produced **the same morning**, one lane
over. The §1b sibling-repo scan added on 07-30 is what caught it. Recording that the addition paid.

## Back-annotation: what the sibling lane established

Executed 2026-08-02 on Lelli, McGaugh & Schombert 2016 SPARC — 2,622 points, 149 galaxies
(Q ≤ 2, inclination > 30°, e_Vobs/Vobs < 0.10), Υ_disk = 0.5. No γ, no ρ_crit, no compander, no fitting.

The site's galaxy sector states `f_DM = 1 − C(ρ)`, hence `B_req ≡ g_obs/g_bar = 1/C(ρ)` must be a single
monotone function of ρ alone with zero intrinsic scatter. MOND makes the identically-shaped claim about
`g_bar`. Two one-variable claims about the same y, so the sharpest test is a head-to-head conditional
scatter:

| conditioner | robust σ(log B_req) |
|---|---|
| log g_bar (MOND/RAR) | **0.1178 dex** — reproduces Lelli et al. 2017's published 0.11–0.13 |
| log ρ (this framework) | **0.1613 dex** |

Ratio 1.37×; quadrature excess 0.110 dex, which is 3.2× the RAR's ~0.034 dex intrinsic-scatter budget.
Giving the framework 175 free per-galaxy offsets — absorbing scale height, distance, M/L and inclination
errors wholesale, applied symmetrically to both arms — makes it **worse** (1.77×): the excess is
*within*-galaxy, i.e. it is the radial shape of the density profile, not a per-galaxy calibration error.

Decisive test: partial correlation of RAR residuals with local-density residuals at fixed `g_bar`,
where the framework requires r < 0. Measured **r = +0.0012**, galaxy-block bootstrap 95% CI
**[−0.076, +0.081]**, ≤ 0.7% of residual variance (≤ 2.6% allowing that 75% of the dρ variance is noise).
The same machinery recovers injected correlations of r = 0.10, 0.20, 0.30, so the null is a null and not
a dead pipeline.

Constructive bound, and the piece worth quoting: interpolating `log u_α = (1−α)·log Σ_local + α·log g_bar`
gives a monotone improvement to a minimum at α = 1.00, bootstrap 95% CI **[0.75, 1.00]**. **SPARC admits
at most 25% weight on a local-density variable.** This sector sits at 100%.

And the escape route is narrowed rather than merely asserted shut: the published escape from the
algebraic no-go (Burrage, Copeland & Millington 2017, the symmetron counterexample this program found on
07-27) works by making the coupling *differential*, which smooths ρ over a range λ. Convolving the
measured Σ(r) with an exponential kernel, **no λ recovers g_bar performance** — best 1.21× at λ ≈ 7 disk
scale lengths, already larger than the galaxy, and degrading beyond. The reason is structural:
`g_bar = GM(<r)/r²` is not a smoothed density, and its 1/r² cannot come from any convolution of Σ.

## The scoping correction, which is mine and not the finding's

The finding's propagation recommendation #2 tells the maintainer to add this to `/for-researchers` as a
fourth leg "stronger than all three because it assumes no functional form," on the grounds that Lelli et
al. "tested surface columns, never volumetric ρ — this is that test, run."

**That justification is inconsistent with the finding's own §2 and its own Scope block.** §2 states that
at constant scale height `log ρ = log Σ − const`, and conditional scatter is invariant under a global
shift — so *the headline number is a statement about the directly measured Spitzer surface density*, and
needs no thickness model. The Scope block says the same thing from the other side: "It does not test the
*volumetric* claim independently of the surface-density one — at fixed h they are the same statement."

And this program has already audited the relevant prior art. `explorations/2026-07-23-locality-nogo-prior-art-audit.md`,
row 13, from the arXiv PDF: Lelli et al. 2017 tested RAR residuals against radius, **stellar surface
density at R**, baryonic mass, R_eff, effective surface brightness and gas fraction, and reported no
significant correlation. That audit's own note — "the site's credit line is loose; Lelli tested columns
(Σ), not ρ. Fix wording" — was aimed at a *credit* line. It applies with more force to a *novelty* line.

So the honest ordering of the finding's own content is not the order it presents:

| § | content | status |
|---|---|---|
| §4 decisive null (r = +0.001, ≤0.7% variance) | RAR residuals carry no local-density information | **quantification of published prior art.** Lelli+2017 reported this null on Σ at R without a CI, an attenuation bound, or a positive control. Those three additions are real and worth having. Discovery, no. |
| §6 smoothing scan / the 1/r² argument | no kernel on Σ reaches g_bar; differential coupling must reconstruct enclosed mass, not smooth ρ | **new, and the most valuable piece.** It converts a live escape hatch (BCM 2017) into a specific, checkable obligation. |
| §7 α-admixture bound ≤ 25% | how local may the organizing variable be | **new, reusable, framework-independent.** |

None of this touches the verdict. It changes what may be *cited as novel*, and it matters because a
preprint that led with §4 would be scooped by a 2017 paper on the same dataset. Flagged for the site
lane; not edited there.

I note the pattern without over-reading it: this is the seventh consecutive day on which a correct
computation in this program carried a headline sentence that claimed slightly more, or slightly other,
than the computation supported — a normalisation (07-27), an estimator (07-29), an M/L (07-30), a
coupling (08-01), a robustness claim of mine (08-02), and now a novelty attribution. The arithmetic has
been right every time. The sentence is where the program leaks.

## γ = 1/2, re-verified here before propagating

The triage's central algebraic claim is that at γ = 1/2 the compander becomes MOND's simple μ exactly.
Verified numerically rather than accepted:

```
tanh(½·ln(1+x))  vs  x/(x+2)  vs  μ_simple(x/2) = (x/2)/(1+x/2)
x ∈ {1e-6, 1e-3, 0.1, 0.5, 1, 2, 10, 1e3, 1e6}
max |deviation| = 5.55e-17          (machine precision — an identity, not an approximation)
```

And the degeneracy caveat, which is the part that costs the program a number it has been quoting: the
O(x²) term separating γ from ρ_crit carries the coefficient γ(2γ−1)/2, evaluated here as
**−0.0319 at γ = 0.425, −0.0054 at γ = 0.489, exactly 0.000 at γ = 1/2, +0.060 at γ = 0.600**. It vanishes
at essentially the SPARC-preferred value, so "γ ≈ 0.489 preferred by SPARC" is degenerate with the ρ_crit
prior and is not a standalone measurement of γ. This qualifies the TEST-11/Cassini γ-interval 0.425–0.600
wherever that interval is cited as a number.

## What this does and does not change

- **Refutation count: unchanged at 6.** The sibling finding says so itself and gives the right reason —
  this is the same underlying failure (local coupling) executed form-free, not an independent refutation,
  and inflating the tally would mirror the scope error already on the program's record in the opposite
  direction. Concur, and record the concurrence rather than the number.
- **Bucket 0: unchanged (0).** Nothing here is a novel confirmed prediction; the identity is a
  reparametrization statement and the bound is a negative.
- **The "locality no-go" is no longer only an analytic prior-art instance.** The record here has carried
  it since S689 as "reproduces a known Milgrom-2005 locality no-go, not a novel result" — an argument.
  It now also has an empirical, form-free, positive-control-validated measurement on SPARC with a
  quantitative admixture bound. The prior-art verdict on the *theorem* is untouched. What is new is that
  the program can now say how much locality the data allow, which the theorem never did.
- **Propagated to** `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — a status note
  under the "MOND-Synchronism Unification" heading, whose 2026-era sentence *"a breakthrough discovery:
  MOND and Synchronism are the same physics in different parameterizations"* is now **exactly true and
  therefore no longer a discovery**. That inversion is the cleanest statement of the reparametrization
  verdict the whitepaper has had, and it arrived from the framework's own preferred parameter value.

## So what

The galaxy sector's whole content is one substitution — `μ(g_bar) → C(ρ)` — and on 2026-08-02 the
program established, in one repository, that the substitution is all there is, and in another, that the
data permit it at most 25% of the way. Those two statements are the premise and the conclusion of the
same argument, they were written 48 minutes apart, and neither author could see the other. The finding
is not that the galaxy sector fails; that has been the standing verdict for months. The finding is that
the program's most decisive morning of the year was invisible to itself, and the only reason it is
visible now is a scan list that was widened after it cost a full pass three days ago.

— Publisher
