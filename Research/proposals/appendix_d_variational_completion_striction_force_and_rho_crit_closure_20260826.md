# Appendix D §D.2 — the variational completion, the omitted striction force, and a new a-priori closure

**Back-annotated from synchronism-site explorer, 2026-08-26.**
Source finding: `synchronism-site/explorer/findings/l2-is-not-l3-for-a-disc-and-the-action-adds-a-force-the-tests-omit.md`
Scripts + full output: `synchronism-site/explorer/findings/scripts/l2_vs_l3_and_the_missing_striction_force.py`

Follows the 2026-08-25 back-annotation establishing that Appendix D's field equation **is** Refracted
Gravity (Matsakos & Diaferio 2016, arXiv:1603.04943), with `C_Ω ≡ ε` to 2.2e-16.

---

## 1. What was executed

M&D 2016 **§2.2.1** ("Assumptions for the permittivity"), **Eq. (2.6)**, wrote down the candidate
Lagrangian for this exact equation and deferred it:

> "The consequences of a variational approach applied to a possible RG Lagrangian of the form
> **L = ε/8πG (∇Φ)² + ρΦ** should also be investigated. However, in this paper we limit our study to
> assess whether this new idea appears promising, at least phenomenologically. If so, we will explore
> all these fundamental issues elsewhere."

> **[CITATION CORRECTED 2026-08-27 — Publisher, from the paper's own text.]** This proposal as filed
> gave the location as "§6". §6 of arXiv:1603.04943 is *Predictions* (§6.1 globular clusters … §6.5
> galaxy groups and clusters) and contains no Lagrangian. The quotation itself is verbatim and the
> substance is unaffected; only the pointer was wrong. Noted because it is the **third consecutive
> pass** to carry a citation-precision defect on this one paper family (08-25: the routed cite named
> a non-author, "Sanna, *Pipino*, Diaferio"; 08-26: `C_ρ`'s residual was quoted at an unstated γ),
> and because a wrong section pointer is the kind of error that survives review — a reader who
> follows it finds a plausible-looking section and no contradiction. Two sign conventions also
> differ and neither is an error: M&D write `L = +ε(∇Φ)²/8πG + ρΦ`, this proposal writes
> `L = −C|∇Φ|²/8πG − ρΦ`; the overall sign flip leaves the Euler–Lagrange equations identical. Use
> **M&D's** form when the result is offered upstream.

Ten years open. Executed here.

> **[PRIOR-ART GATE RUN 2026-08-27 — Publisher — result: CLEAN, and the reason is structural.]**
> Before this striction result travels outward it was screened against the one place it could
> already exist: the covariant completion, **Sanna, Matsakos & Diaferio 2023** (A&A 674 A209;
> arXiv:2109.11217), read at source. It is a scalar–tensor theory,
> `S = (1/16πG)∫d⁴x√g[φR + (W(φ)/φ)∇^α φ∇_α φ + 2V(φ)] + ∫d⁴x√g L_m` with `W(φ) = −1`,
> `V(φ) = −Ξφ`, and **matter minimally coupled to the metric alone**; the weak-field limit is
> `∇·(φ∇Φ) ≃ 8πGρ` with `φ = 2ε`. No electrostriction term appears — **and it cannot**, because in
> CRG `φ` is an independent dynamical field obeying its own equation of motion, *not* a constitutive
> function `ε(ρ)`. There is no `ε′(ρ)` to vary against, so there is no striction force. The paper
> also nowhere states that the 2016 non-covariant formulation had a third-law or momentum-conservation
> problem. **Two things follow.** (i) The striction result is **not** prior art — it is a genuine
> open contribution to an active external programme. (ii) It sharpens, rather than softens, this
> repository's standing observation that *the covariant escape drops locality* (4th instance): the
> covariant completion purchases momentum conservation precisely by giving up the density-keying
> that is the whole content of `C(ρ)`. Coverage stated: full text of both papers via arXiv HTML/PDF;
> the CRG action, `W`, `V`, `φ = 2ε` and the minimal matter coupling were read at source, the
> surrounding derivations were not re-derived here.

## 2. The result, analytically

With `L = −C(ρ)|∇Φ|²/(8πG) − ρΦ`:

- **Vary Φ** → `∇·[C(ρ)∇Φ] = 4πGρ`, i.e. Appendix D §D.2 exactly (site label **L2**).
- **Vary the matter** (`δρ = −∇·(ρξ)`; ρ enters through `C(ρ)` as well as `ρΦ`) →

```
g = −∇Φ − ∇Ψ ,      Ψ ≡ C′(ρ) |∇Φ|² / (8πG)                             (★)
```

`−∇Ψ` is the gravitational analogue of the **Korteweg–Helmholtz electrostriction force** in a
dielectric whose permittivity depends on mass density (Landau & Lifshitz, *ECM* §15). The physics is
textbook; RG is *defined* by the dielectric analogy; the term was simply never carried over.

**It is not optional.** For the naive law `g = −∇Φ`, splitting off a divergence gives, for an
isolated system, `∫ρ(−∇Φ)dV = −(1/8πG)∫|∇Φ|²∇C dV ≠ 0` — a net self-force. Integrating the
striction force by parts gives exactly `+(1/8πG)∫|∇Φ|²∇C dV`. Hence:

> **Either the theory carries `−∇Ψ`, or it violates Newton's third law.**

## 3. The result, numerically

Axisymmetric finite-volume solve. Validation: uniform sphere vs `−GM/r` to **7e-4**; exponential
disc vs **exact Hankel transform** to **0.13%**; Gauss check on the discrete solution 0.984–1.003;
`L2 ≡ L3` reproduced to **0.18%** *when the source is spherical*. All headline numbers converged
to <2% over resolution (200×220 → 300×340) and box (400 → 1500 kpc).

**(a) `g = g_bar/C` (L3) is not the solution of L2 for a disc.** Gauss's law sets the field by the
permittivity averaged over the *enclosing sphere*; a thin disc is measure-zero on its own sphere
(at r = 8 kpc, ⟨C⟩ = 0.148 vs C(midplane) = 0.923). Max **L2/L3 = 5.89**. Closed form:
`max(L2/L3) → 1/C_min = B_max`, measured **3.178** against `1/Ω_m = 3.175`.

**(b) The omitted force.** `0.72×` gravity radially; **up to 164× K_z vertically**, in a thin shell
at the height where ρ crosses ρ_c. Generic: 245× at q=1, 5.2× at q=2, 10× for the tanh-log form.

**(c) Third-law violation.** Two equal-mass unequal-density spheres, uniform grid, C=1 control =
0.001: naive law gives **−3.6365 ×** the pair's own mutual force; the identity reproduces it to
**0.4%**; the striction term cancels it to **0.0%**. Isolated pair self-accelerates at **0.047 a₀**.

## 4. The new closure — citable, needs no data

ρ_crit scan for the framework's own `C = C_min + (1−C_min)tanh(γ ln(1+ρ/ρ_crit))`:

| ρ_crit (M☉/pc³) | knee at R | max L2/L3 | vertical striction / K_z |
|---|---|---|---|
| 1e-3 | 21 kpc | 3.05 | 858 |
| 1e-2 | 14 kpc | 2.80 | 104 |
| 1e-1 | 8 kpc | 2.43 | 12.9 |
| 1 | 3 kpc | 1.87 | 0.51 |
| 652 (= 0.029·150²) | never | 1.004 | 0.0005 |

> **A density-keyed permittivity theory cannot place its knee inside a disc without introducing a
> vertical force that exceeds the vertical gravity.** The omitted term drops below K_z only for a
> knee inside R ≲ 3 kpc — inside the radius where no boost is needed.

This is a constraint on the class, derived, requiring no data, and it is a direct answer to M&D's
deferred question.

## 5. Consequences for PREDICTIONS.md and the archive

1. **Appendix D states four mutually inconsistent dynamical laws** (checked, not conjectured):

   | where | statement | law |
   |---|---|---|
   | §D.2 | `∇·[C∇Φ] = 4πGρ` | L2 |
   | §D.5 | `S = ∫[−m√(−g ẋẋ) − λU(x)]dτ`, `U ∝ −ln C(ρ)`, **λ "to be calibrated"** | L2 + a fifth force |
   | §D.6.1 | `(1/r²)d/dr[r²dΦ/dr] = 4πGρ/C` | **L1** |
   | §D.6.2 | effective metric sourced by `ρ_eff = ρ/C` | L1, covariantly |

   λ is never calibrated anywhere in the appendix, and no test in either repo uses §D.5's force.
   **And §D.5's force cannot be the missing conservation term**, by a short theorem: any force
   `−∇f(ρ)` gives **zero** net force on an isolated system, since `∫ρ∇f(ρ)dV = ∫∇F(ρ)dV = 0` with
   `F′ = ρf′`. So no purely density-dependent extra force can restore the third law — the restoring
   term must depend on the field, and (★) is uniquely it. Verified numerically: §D.5's net force is
   **1.9e−08** of the violation at 30 kpc separation, and the λ that would be needed disagrees by
   **12×** across three separations.
2. **All galaxy-sector results computed from `g = g_bar/C` are L3 results.** Their *verdicts* are
   unaffected at the archive's own `ρ_crit = 0.029·V²` (every term ≤ 0.4% there), but the
   *mechanism* recorded for the failure — "C falls outward, required boost rises" — is the L3
   mechanism. L2's boost saturates at the ceiling, giving a near-constant rescaling.
3. **At `ρ_crit = 0.029·V²` the galaxy sector is `G → G/Ω_m`** — a constant rescaling of the
   gravitational constant by 3.17, C pinned at the floor from 1 kpc outward. This is the reason the
   Lagrangian question was safely ignorable, and it should be recorded as such rather than as a
   result about Lagrangians.
4. **Third structural difference from MOND, and the first non-constant one:** AQUAL is variational
   by construction (Bekenstein & Milgrom 1984, for exactly this reason); `C(ρ)` gravity is not
   unless it carries `−∇Ψ`.
5. **q < 1/2 pathology (RG-specific, offered upstream):** `Ψ ∝ ρ^{2q−1}|∇Φ|²` diverges as ρ→0 for
   q < 1/2. Cesare et al. 2020's fitted **q = 0.47** sits 0.03 below that boundary, so RG's own
   best-fit striction potential does not decay into vacuum. Cheap to check; appears unchecked.

## 6. Open

- **The covariant branch.** Sanna, Matsakos & Diaferio 2023 (A&A 674 A209; arXiv:2109.11217) *[author list corrected 2026-08-29 — "Pipino" was carried over from the routed note; he is on no RG paper]* complete this to a
  scalar–tensor theory with `φ = 2ε` in the weak-field limit. There `φ` is an *independent dynamical
  field*, so **ε is not a local function of ρ** — momentum is conserved and (★) dissolves, at the
  price of the locality that is the framework's advertised distinguishing feature. Fourth time an
  escape from this obstruction has required dropping locality (BCM 2017; Verlinde `M_B(<r)`;
  MOND Σ; the 2026-08-19 ∇Φ axis). **Nobody in this programme has examined branch (C).** It is the
  highest-value open item in the galaxy sector.
- **Does a massless tracer feel `−∇Ψ`?** This is M&D's flagged open SEP question. If gas feels it
  and stars do not, that is itself a sharp, cheap prediction.
- **Multi-component discs.** A thicker gas layer broadens the transition shell. Predicted: peak
  ratio falls, §4's monotone trade-off survives (it is set by *where the knee is*). Untested.
