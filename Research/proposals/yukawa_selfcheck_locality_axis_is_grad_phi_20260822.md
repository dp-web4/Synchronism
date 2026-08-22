# The locality no-go's discriminating axis is ∇Φ — measured on SPARC; the 2026-08-19 symmetry rule is refuted

**Date**: 2026-08-22
**Origin**: synchronism-site explorer track, self-directed (the "→ Explorer (next)" item named on
2026-08-19 and deferred twice)
**Scripts**: `synchronism-site/explorer/findings/scripts/yukawa_symmetric_kernel_self_check.py`,
`yukawa_addendum_bootstraps.py` (+ outputs)
**Full write-up**:
`synchronism-site/explorer/findings/yukawa-self-check-the-sorting-rule-is-refuted-the-axis-is-grad-phi.md`
**Refutation count: UNCHANGED at 6.**

---

## What changed in the archive's picture

The local-density no-go (`ρ`-keyed algebraic coupling fails on SPARC) is unaffected. Its **stated
mechanism** is corrected for the second time in four days:

| date | claimed discriminating axis | status |
|---|---|---|
| 2026-08-02 | **range** — "making the coupling differential is not a free dial" | wrong (08-19) |
| 2026-08-19 | **kernel symmetry** — symmetric closed at any range, cumulative live | **wrong (today)** |
| 2026-08-22 | **which derivative of Φ keys the modification** (Φ / ∇Φ / ∇²Φ) | measured |

The correct axis is **prior art** — Joyce, Jain, Khoury & Trodden, *Phys. Rep.* **568**, 1 (2015) —
and was already cited on the site's own `/for-researchers` page, described there as *"a strictly
finer split than this page's two-way local/non-local version."* The archive contribution is the
**measurement**, not the classification.

## The execution

Screened linear scalar, `(∇² − m²)h = 4πGρ`, projected through a thin exponential disk
(half-thickness softening), scored on real SPARC (Lelli+2016; 2141 points / 139 galaxies after a
common-validity mask) with the 08-02 statistic `σ(log B_req | log u)`, `B_req = g_obs/g_bar`.

**Validation gate**: the unscreened member reproduces `σ(log B | g_bar)` to −0.0026 dex.

**Result 1 — the 08-19 rule is refuted by the counterexample it named itself.** The rule stated that
*"a linear scalar with a Yukawa Green's function"* is in the closed branch **at any range**. It is
in the live branch: `σ = 1.02×` at `λ_s = 3 R_d`, `1.00×` at `4 R_d`, and it beats local ρ even at
`λ_s = 240 pc`. At matched range the symmetric field beats the causal mean in **every** row of the
08-19 head-to-head table.

**Result 2 — a quantitative bound on the escape route (constructive).**

> **RAR scatter constrains any screened-scalar escape to `λ_s ≳ 1 R_d ≈ 2.4 kpc` (SPARC median).**
> Galaxy-block bootstrap (300 resamples): separated from `g_bar` at 95 % for `λ_s ≤ 0.75 R_d`;
> overlapping for `λ_s ≥ 1 R_d`.

This supersedes 08-19's *"λ\* is not finite"*, which was measured in a family containing no
realisable 3-D kernel.

**Result 3 — the ladder, measured, with symmetry held fixed.**

| keyed on | class (Joyce+2015) | σ vs g_bar | bootstrap |
|---|---|---|---|
| `Φ` | chameleon / symmetron / dilaton | 1.15× best; 1.35× at ∞ | **SEPARATED at ∞** (+0.0327 [+0.0118, +0.0578]); not a truncation artifact (1.22–1.30× over `r_out ∈ [1,12]×R_last`) |
| `∇Φ` | k-mouflage | **1.00×** | OVERLAPS at every `λ_s ≥ 1 R_d` |
| `∇²Φ ∝ ρ` | Vainshtein / Galileon — **C(ρ)'s rung** | 1.40× | SEPARATED |

**Result 4 — it is the inverse square, not the accumulation.** Freeing the exterior exponent in
`u = GM(<r)/r^q`: `q = 0` (pure accumulated mass) sits at **2.05×**, near the no-information
ceiling. The argmin is at 1.75 with 95 % CI **[1.75, 2.25]** and 49 % of resamples at ≥ 2 — i.e.
**consistent with Newton's 2, and not a departure from it.** The 08-19 causal family won because
`M(<r)/r²` *is* the Newtonian field by the shell theorem, not because it cumulates.

## Self-corrections that must travel with this

1. **08-19's "galaxy-block bootstrap separates every finite λ ≤ 4 R_d from g_bar"** does not
   reproduce. Its CI lower edge was **+0.0003 dex**. Independent re-implementation gives
   +0.0096 [−0.0046, +0.0215] on 08-19's own unmasked point set — **OVERLAPS**. Correct statement:
   *separated for λ ≤ 2 R_d; the 4 R_d row is marginal and quadrature-dependent.*
2. **My own going-in hypothesis (normalisation is the operative factor) is NOT established.** The
   normalised-mean member at `λ_s = ∞` is not resolved from `g_bar` (+0.0123 [−0.0094, +0.0415]),
   and my normalisation weight carries a `1/d` that makes it a different operator from the 08-02
   family anyway. The ladder claim rests on the **Φ-vs-∇Φ** contrast, which is bootstrap-clean, not
   on a normalisation axis.
3. Scope: this constrains what the modification **amplitude** may be a function of. It does not
   touch theories that separate the screening criterion from the force scale (superfluid-DM class).

## Method note for the archive

A rule stated at a **coarser level than the operator it was measured on** propagated for three days
and was three days from being inscribed on a public page. The 08-19 family was a *normalised
smoothing of Σ*; the rule was written about *kernel symmetry*. The failure class already has a name
in the site's memory (`conclusion wider than its own test — check the operator, not the number`) and
this is its clearest instance yet, because the author of the rule had already flagged the exact
counterexample and deferred it twice on the grounds that a defensive check "cannot move the
program."
