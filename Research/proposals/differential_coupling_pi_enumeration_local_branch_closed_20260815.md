# The local differential-coupling branch is closed; the 2026-07-27 scope demotion rests on a misclassification

**Date**: 2026-08-15
**Origin**: synchronism-site explorer track
**Site finding**: `explorer/findings/differential-coupling-pi-enumeration-local-branch-closed.md`
**Script**: `explorer/findings/scripts/differential_coupling_pi_enumeration_real_sparc.py` (+ output)
**Refutation count**: **UNCHANGED at 6** — constructive-lead closure, not a new refutation.

---

## Claim

Within the class of **local** couplings `F(ρ, ∇ρ, ∇²ρ, ∇³ρ; G, a₀)` entering an AQUAL-type field
equation `∇·[F ∇Φ] = 4πGρ`, **no member reproduces the radial acceleration relation.** The galaxy
sector's last un-eliminated constructive direction is closed.

Separately: the 2026-07-27 demotion of the locality no-go's scope was made on the authority of a
counterexample (BCM 2017) that does not belong to the class the no-go covers. The scope should be
**restored**, not merely re-asserted — it now has empirical backing it never had before.

---

## 1. The class is exactly two-dimensional (this is what makes it a closure)

Buckingham π on `F(ρ, s₁ = |∇ρ|, s₂ = ∇²ρ; G, a₀)`: five quantities, dimensional matrix of rank 3
over (M, L, T), so the null space is **exactly 2-dimensional**. A convenient basis:

- **scale group** `x_diff = G ρ² / (a₀ |∇ρ|)` — the differential analogue of `x = ρ/ρ_crit`
- **shape group** `q = ρ ∇²ρ / |∇ρ|²` — free of *both* G and a₀

Adding `∇³ρ` adds exactly one further group. This converts "which F?" from an open search into a
closed question, and lets the **entire class** be tested at once by 2-D conditioning rather than
form-by-form.

**Note in the framework's favour:** the differential branch needs **no new constant**. The algebraic
branch requires a free `ρ_crit` (degenerate with γ); `x_diff` uses only G and a₀ and `q` uses
neither. The branch was strictly less parameterised than the one it would have replaced, which is
why it warranted execution rather than dismissal.

## 2. The vacuum fork (analytic, on the framework's own baseline geometry)

For `ρ ∝ e^{−r/R_d}`: `x_diff = GρR_d/a₀ → 0` as `ρ → 0` (inherits the exact pathology that killed
the algebraic branch on 2026-08-03), while `q ≡ 1` at **every** radius of any exponential (constant
coupling = a renormalisation of G, not a modification of gravity).

Additionally, in the **midplane** of a thin disk the first derivative is purely radial but the second
derivative is dominated by the vertical term, giving `q ≈ f″(0)(R_d/h)² ≈ −2(R_d/h)²` — i.e. `q` is
set by the aspect ratio, the worst-constrained quantity in the problem. The empirical section below
uses radial-only groups, which is the steelman.

## 3. Execution on real SPARC (Lelli+2016; N = 2614, 145 galaxies; Q ≤ 2, i > 30°, e_V/V ≤ 0.10)

Symmetry integrates the operator once, so the **required** coupling is `F_req = g_bar/g_obs`. A
theory is viable iff `F_req` is a tight single-valued function of its π groups. No functional form,
no γ, no ρ_crit, no fitting.

| variable | σ(log F_req \| ·) dex | vs g_bar |
|---|---|---|
| no-information ceiling (unconditional) | 0.3090 | 2.63× |
| log g_bar (target) | **0.1174** | 1.00× |
| log ρ (algebraic) | 0.1607 | 1.36× |
| **log x_diff** | **0.1945** | **1.65×** |
| **q** | **0.3002** | **2.55×** |
| (x_diff, q) — **full class, joint** | **0.1798** | **1.53×** |
| (g_bar, q) — does q add to MOND? | 0.1171 | 0.99× |
| (g_bar, g_bar) — binning-cost control | 0.1210 | 1.03× |

**The decisive test.** Correlation of each group with the RAR residual `dB`, controlling for g_bar:

| group | r(dB, · \| g_bar) | variance explained |
|---|---|---|
| log x_diff | −0.0335 | 0.11 % |
| q | −0.0225 | 0.05 % |
| \|∇ln ρ\| | +0.0382 | 0.15 % |
| dln ρ/dln r | −0.0394 | 0.16 % |
| 3rd-order group c₃ | −0.0279 | 0.08 % |
| *(local ρ, 2026-08-02)* | — | *≤ 0.70 %* |

**Every group in the complete differential class carries less information about the missing physics
than the algebraic variable it was proposed to replace.**

## 4. The four objections, pre-empted

- **Permutation null** (200×, same estimator): `q` measured 0.3002 vs null 0.3067 ± 0.0024, **z = −2.7**
  — indistinguishable from its own destroyed version. "q is just noisy" is not available as a rescue,
  because the null was built from the same noisy variable.
- **Is x_diff even new?** `log x_diff = 0.912 log ρ + c`, r = 0.8534 (72.8 % of variance is ρ); the
  residual after projecting ρ out scores **0.2893**, at the 0.3090 ceiling. The genuinely differential
  part carries nothing.
- **Degeneracy** (the 07-28 kill criterion): honest result running *against* the tidy argument — in
  real galaxies `q` is 95.3 % within-galaxy and the density scale length L has 0.382 dex within-galaxy
  scatter, so the escape from the horn-2 degeneracy genuinely exists. **It is simply empty.** That is
  a stronger closure, because it does not lean on the exponential-disk idealisation.
- **ϒ systematic** — the dominant one, which dissolved the γ concordance on 2026-08-14: best
  differential group / g_bar is **flat at 1.34–1.40× across ϒ ∈ [0.30, 0.80]**. Also robust across
  16 combinations of derivative estimator × window × scale height × gas treatment.

This is a useful calibration datum for the marginalisation guardrail itself: the ϒ check removed a
false affirmation on 08-14 and leaves this closure standing. The guardrail is **discriminating, not
a universal solvent.**

## 5. The scope demotion should be reversed

BCM 2017 = Burrage, Copeland & Millington, *"Radial acceleration relation from symmetron fifth
forces"*, PRD 95, 064050 (erratum PRD 95, 129902). φ solves a **nonlinear screened field equation**
`∇²φ = V_eff′(φ, ρ)`; φ(r) depends on the whole surrounding matter configuration. BCM's force is the
gradient of a function of **φ**, not of the local **ρ** — and BCM's own closed form,

```
g_sym = g_bar / ( exp √(g_bar/g†) − 1 ),   g† ≈ 1.20 × 10⁻¹⁰ m s⁻²
```

is written in **g_bar**, an integral of ρ. A local function of ρ and its derivatives cannot be
written this way.

So the escape class the site files it under — *"differential **local**-density coupling
(symmetron-class)"* — is **empty** (§1–4) **and unpopulated by its own exemplar**. The 07-27
rescoping ("the real axis was never local-vs-non-local; it was algebraic-vs-differential") should be
withdrawn; the real axis *was* locality, which is what the site's own text says three sentences later
when it names Milgrom's non-locality theorem (astro-ph/0510117) as the root obstruction.

**Restoring the scope does not rescue C(ρ)** — algebraic stays dead, and the differential branch that
might have rescued it is now closed too. It strengthens the *negative-result* deliverable.

**And the surviving branch costs EFE = 0.** Screened scalars are environment-dependent by
construction and generically violate the SEP. The framework's EFE = 0 claim rests on the algebraic
C(ρ)·g law satisfying the SEP by construction; buying non-locality to survive the RAR buys EFE ≠ 0.

## 6. Scope and what would overturn this

**Not established:** anything about **non-local** couplings (screened scalars, integro-differential
couplings, auxiliary-field theories) — that branch is alive and is where BCM lives; anything about
couplings keyed on non-density variables (e.g. `∇Φ`, which is MOND and which the data selects). The
`F g_obs = g_bar` symmetry integration is the site's own working approximation for a disk, not exact.

**Would overturn:** a local `F(ρ, ∇ρ, …)` predicting `dB` at the ≥ 5 % level on this sample; or a
demonstration that the symmetry integration fails badly enough for a disk that `F_req` is the wrong
target.

**Open gate (dp).** Preprint statement available: *no local density-keyed coupling, algebraic or
differential at any derivative order, reproduces the RAR.* Gate item: read BCM's field equations
directly (abstract + closed form + screening literature suffice for the classification, but not for
publication). Credit line required — local-density instance of Milgrom's non-locality theorem, not a
new theorem.

## 7. Next question for the sector

The non-local branch is the only one left. The sharp question is no longer "can it fit the RAR"
(BCM shows a screened scalar can) but **"does anything in that branch keep EFE = 0?"** — because that
is the framework's only claimed discriminator, and §5 argues the branch forfeits it generically.
