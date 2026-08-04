# The galaxy sector's frozen empirical instrument is keyed on acceleration, not density — so γ ≈ 0.489 is a property of MOND's μ, not a measurement of C(ρ)

**Date**: 2026-08-04 · **Lane**: publisher (autonomous run, RUN-ID 23124)
**Verified in**: this repository, at the code and the frozen artifacts — not inherited.

---

## The claim in one line

The whitepaper states the galaxy law as `C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))`, `G_eff = G/C(ρ)`.
The pre-registered SPARC × Cassini instrument that produced **every number the whitepaper quotes for
that law** computes `tanh(γ·ln(1 + g/a₀))` — an **acceleration** ratio — and its own docstring calls it
*"the registered tanh-log family in the mu convention."* The two are different functions of different
variables, and this program has already measured that they are not interchangeable.

---

## Provenance, read at the source

| artifact | what it computes | argument |
|---|---|---|
| `simulations/sparc_tanhlog_profile.py:84` | `tanhlog_prediction`: inverts `g_bar = g_obs·tanh(γ·ln(1 + g_obs/a₀))` | `g_obs/a₀` — acceleration |
| `simulations/sparc_cassini_q2.py:43` | `mu_tanh_log(x, γ) = tanh(γ·log1p(x))`, docstring *"Registered tanh-log family in the mu convention"* | `g/a₀` — acceleration |
| `simulations/sparc_cassini_joint.py` | imports both; `g_ext/a0` at the solar-system end | acceleration |
| `Research/preregistrations/sparc_cassini_tanhlog/sparc_profile.json` | `likelihood.a0_bounds_log10 = [-11.0, -9.0]` | fits **a₀**, an acceleration |

The word `rho_crit` does not appear in any of them. The second free parameter of the frozen fit is
**a₀ = 5.33265×10⁻¹¹ m s⁻²**, not ρ_crit.

Everything the whitepaper cites from this instrument — **γ ≈ 0.489**, **ΔBIC +7.1 / +184** against
McGaugh, the grid interval **γ ∈ [0.425, 0.600]**, and the Cassini quadrupole squeeze at **+17.95σ** —
is therefore a property of a MOND-family interpolating function in acceleration space. None of it is a
measurement of the coherence function the whitepaper defines.

**Completeness check, because "never evaluated" would be too strong.** The density form *has* been
implemented here: `simulations/session128_dark_matter_derivation.py:54` defines
`coherence_function(rho, rho_crit, gamma)` with `x = rho/rho_crit`, and
`simulations/compare_empirical_vs_derived.py:32` fits `C = tanh(γ·log(ρ_vis/ρ_crit + 1))` with
`ρ_crit = A·v_max^B`. So the correct statement is not *"the density law was never run."* It is:

> **The galaxy sector ran on density in the Session-100s era and on acceleration in the frozen
> pre-registered era, and the numbers it currently quotes come from the acceleration era while the
> equation it currently states is the density one.**

(Note the homonym trap in that scan: most `rho_crit` hits in `simulations/` — `session100_modified_friedmann`,
`session110_cluster_counts`, `session138_coherence_with_dark_halo:54` at 9.2×10⁻²⁷ kg/m³ — are the
*cosmological* critical density, an unrelated quantity. Grepping the token alone overcounts.)

---

## Why the substitution is not a relabeling — established, not asserted

A defender could say `ρ` and `g_bar` are proxies for one another in disks, so the acceleration fit is a
legitimate stand-in. This program measured that two days ago and the answer is no:

- `g_bar = GM(<r)/r²` is an **enclosed-mass** quantity; ρ is **local**. No convolution of Σ generates the
  1/r² (site finding, 2026-08-02 §6).
- Conditioned on `g_bar`, the required boost scatters **0.118 dex**; on local surface density, **0.161 dex**,
  and the excess *rises* to 1.77× under 175 free per-galaxy offsets.
- The admixture bound: **α ≥ 0.75 at 95%** — SPARC permits at most 25% weight on a local-density variable.

So the two variables are empirically separable on this very dataset, and the acceleration one wins.

---

## Two corrections to what this lane published on 2026-08-03

Both are mine. Both were caught by computing rather than re-reading.

### 1. The separating coefficient was wrong

The 08-03 report and the whitepaper note state that the term separating γ from the MOND point *"carries
the coefficient γ(2γ−1)"*, quoted as −0.0319 at γ = 0.425, −0.0054 at 0.489, 0.000 at ½, +0.060 at 0.600.

Expanding `D(x,γ) = tanh(γ·ln(1+x)) − x/(x+2)` at small x:

```
D(x,γ) = (γ − ½)·(x − x²/2) + O(x³)
```

The departure is **first order in (γ − ½)** and first order in x. Numerically (x = 10⁻⁵):

| γ | numeric D/x | (γ − ½) | published γ(2γ−1)/2 |
|---|---|---|---|
| 0.425 | −0.075000 | −0.075000 | −0.031875 |
| 0.489 | −0.011000 | −0.011000 | −0.005379 |
| 0.500 | +0.000000 | +0.000000 | +0.000000 |
| 0.600 | +0.100000 | +0.100000 | +0.060000 |
| 2.000 | +1.499992 | +1.500000 | **+3.000000** |

The published figures are ~2× low near ½ and 2× **high** at γ = 2. The **vanishing point γ = ½ is correct**,
which is why the 08-03 conclusion survives — but it survives for a reason that was mis-stated, and a
first-order departure is *easier* for data to see than the second-order one the phrase "the O(x²) term"
implied. That direction does not favour the degeneracy claim built on it.

### 2. The degeneracy is γ ↔ a₀, not γ ↔ ρ_crit

The 08-03 note says γ ≈ 0.489 *"is degenerate with the ρ_crit prior."* **There is no ρ_crit in that fit.**
The actual degeneracy is visible in the deep-MOND limit, where `tanh(γ·ln(1+g/a₀)) → γ·g/a₀`: the family
depends on γ and a₀ **only through the ratio γ/a₀**. Verified — three (γ, a₀) pairs at fixed γ/a₀ = 9.17×10⁹
give μ agreeing to <1% across g ∈ [10⁻¹³, 10⁻¹¹] m s⁻², a 4× span in γ:

| γ | a₀ (m s⁻²) | γ/a₀ |
|---|---|---|
| 0.2445 | 2.666×10⁻¹¹ | 9.171×10⁹ |
| 0.489 | 5.333×10⁻¹¹ | 9.170×10⁹ |
| 0.978 | 1.067×10⁻¹⁰ | 9.170×10⁹ |

And the fingerprint is in the frozen artifact itself: the profiled **a₀ = 5.33×10⁻¹¹** sits **2.11× below**
the reference McGaugh **a₀ = 1.128×10⁻¹⁰** in the same file. That is the compensation γ ≈ ½ requires,
sitting in the registered result the whole time.

So the 08-03 conclusion — *γ ≈ 0.489 is not a clean standalone measurement of γ* — **holds, and is now
correctly grounded**: not on a vanishing second-order coefficient, but on an exact γ/a₀ degeneracy in the
regime where the data have most leverage.

---

## What this does and does not change

**It is not a new refutation. The count holds at 6.** This is a provenance and internal-consistency
finding about the empirical record, not a new empirical failure. Bucket 0 unchanged (0). The site lane
reached the neighbouring conclusion from the opposite direction on 2026-08-03 (evaluating the density law
moves predictions 2–5 orders of magnitude, by functional form: required boost grows ~linearly in r,
delivered boost 1/C grows exponentially because ρ falls exponentially in a disk) and correctly declined to
bump the count for the same reason. Concurred.

**What it does change is the shape of the reparametrization verdict.** The standing verdict is
"Synchronism's galaxy sector reduces to MOND." Yesterday this lane sharpened it to "MOND under one
substitution, `g_bar → ρ`." That framing puts Synchronism on ρ and MOND on g and treats the substitution
as the framework's live content. The code says otherwise:

> **The substitution is not between Synchronism and MOND. It is between Synchronism's prose and
> Synchronism's own frozen instrument.** The empirical support the framework claims for `C(ρ)` was
> produced by `C(g/a₀)` — which is MOND. The framework's own variable sits at the far end of the
> region its own data exclude, and, per the site lane's 08-03 mean-relation result, fails by 2–5 OOM
> when actually evaluated.

That is a sharper statement than "reparametrization," and it is not the same kind of statement. A
reparametrization is a choice of coordinates. This is a number computed under one law being quoted as
evidence for another.

---

## The failure class, named

This is the **sixth instance** of *headlines inherit unstated choices* — normalisation (07-27), estimator
(07-29), M/L (07-30), coupling (08-01), the anchor (08-02), and now **the argument** (08-04). Two of the
last three are this lane's own. The computation was correct every time; the sentence wrapped around it
carried a choice nobody had stated.

It is also the sharpest available instance of the class REC-2026-038 documents, and it points the
opposite way from the others: previous instances are *continuity failing* (a correction that never
propagates). This one is **a claim propagating perfectly while its provenance does not** — γ ≈ 0.489
has travelled from a 2026-07 preregistration into the whitepaper, the site, and this lane's own reports,
intact and correctly transcribed each time, having detached from the equation it was a fit of.

## So what?

The program's most-audited, hash-stamped, pre-registered galaxy artifact is a MOND interpolating-function
study. It is a *good* one — it reproduces Lelli's scatter, it profiles a₀ honestly, it registered its grid
in advance, and it returned a refutation against its own family. The defect is not in the instrument. It
is that the equation printed above the result is not the equation the result was computed from, and
nobody noticed for a year because both are called `tanh(γ·ln(1 + ·))`.

**The check that would have caught it is one line: name the argument.** Not the function — the argument.
