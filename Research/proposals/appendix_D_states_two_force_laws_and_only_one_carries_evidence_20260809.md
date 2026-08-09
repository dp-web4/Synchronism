# Proposal: Appendix D.2 states two force laws in one equation block; only one carries evidence, and the other is eliminable a priori

**Date**: 2026-08-09
**Lane**: Publisher (autonomous)
**Status**: FILED — §D.2/§D.6.1 equation change is a *physics* decision, gated on dp. Editorial corrections applied.
**Script**: `simulations/publisher_20260809_appendixD_coupling_fork.py` (6 checks, all re-derived in-repo)
**Upstream**: `synchronism-site/explorer/findings/the-archive-has-a-field-equation-and-it-is-not-the-one-the-site-uses.md` (2026-08-08), back-annotated at `Research/proposals/appendix_D_field_equation_is_not_the_site_force_law_20260808.md`

---

## 1. The defect

`manuscripts/Appendix_D_Synchronism_in_General_Relativistic_Form.md` §D.2 (committed **2025-12-01**, `4400d54f`) presents these as one statement, the second following "in this limit" from the first:

```
L1:  ∇²Φ = 4πG ρ/C(ρ)                    (§D.2 opening; restated verbatim in §D.6.1 and §D.7)
L3:  g_obs = g_bar/C(ρ),  V_obs = V_bar/√C   (§D.2 "in this limit")
```

**L3 does not follow from L1.** L3 is the spherical Gauss solution of a different field equation:

```
L2:  ∇·[C(ρ)∇Φ] = 4πGρ      ⇒  C·g·r² = G·M_bar(<r)  ⇒  g = g_bar/C
```

L1 instead makes `ρ/C` the *source*: `g(r) = (G/r²)∫₀^r 4πr'²(ρ/C)dr'`. The two coincide **iff C is spatially constant**, and C is the framework's entire content.

L2 was written into whitepaper §5.15 on **2026-08-04** as a "momentum-conserving completion," on the explicit premise that *"the galaxy sector has no field equation."* One had been on file in this repository for seven months. **The 08-04 amendment is therefore not a rival to Appendix D — it is Appendix D's repair**: it is the equation whose spherical solution is already §D.2's own second line.

## 2. Measurement

`simulations/publisher_20260809_appendixD_coupling_fork.py`. In-repo SPARC mass models (Lelli, McGaugh & Schombert 2016), 113 galaxies at Q ≤ 2 and i > 30°. Both laws are evaluated on the **same** density, built self-consistently as `ρ_sph(r) = (4πr²)⁻¹ d[V_bar²r/G]/dr`, which reproduces the observed baryonic curve exactly and introduces no scale-height or disk-geometry convention.

| check | result |
|---|---|
| (1) L2 ⇒ L3 in spherical symmetry | exact, 0.000e+00 relative |
| (2) L1 ≡ L3 iff C spatially constant | constant-C control: 2.6×10⁻¹⁵ dex; C(ρ) at γ=2: 0.92 dex |
| (3) fork amplitude, outermost point | **median +0.821 dex**, IQR [+0.515, +1.416], range [+0.023, +3.203]; **105/113 > 0.3 dex, 47/113 > 1 dex** |
| (4) γ-dependence | **max \|gap(γ=2) − gap(γ=0.489)\| = 0.0000 dex** |
| (4b) *why* γ-invariant | every point at ρ/ρ_crit ≤ 0.045 (median 6.1×10⁻⁶); C → γρ/ρ_crit to within 2.4%; 1/C ∝ 1/γ in **both** laws ⇒ γ cancels |
| (4c) which is closer to data | **L1 is closer in 113/113 galaxies** at both γ (median +4.64 vs +5.65 dex at γ=2; +5.25 vs +6.26 at γ=0.489) |
| (5) L1's vacuum limit | `ρ/C → ρ_crit/γ` exactly (<10⁻⁴ rel. for γ ∈ [0.25,5] to ρ/ρ_crit = 10⁻¹⁰) |
| (6) provenance | L1: **1** script (`session72_spherical_toy_model.py`, 2025-12-01), **0** whitepaper citations. L3: **26** scripts, every quoted number. L2: **0** scripts |

**Independent reproduction.** The sibling explorer lane reported median 0.81 dex on five galaxies the same morning. This is 113 galaxies with a different density construction and agrees to 0.011 dex on the median. **Confirmed, not accepted.**

**Correction to the upstream framing.** That lane reports γ-invariance as a bare fact. It is a **linear-regime artifact of ρ_crit = 0.029·V_flat²**: no SPARC galaxy reaches C's knee, which is the same fact as the already-banked 2–5 dex miss. Correct scope: *the fork survives every γ at this ρ_crit*. A ρ_crit recalibration reaching the knee would have to be re-run.

**Direction, stated so it cannot be misread.** The reading eliminated in §3 is the one that fits *better*: L1 is closer to the observed accelerations than L2 = L3 in **113/113** galaxies. Both are catastrophically wrong and neither is viable — the elimination is **structural, not goodness-of-fit**. The absolute offsets above are computed on `ρ_sph` and are **not** comparable to §5.15's "2–5 orders of magnitude" (which uses `ρ = Σ/2h`); the *fork* is a ratio of two laws on one ρ and is robust to that choice.

## 3. Two eliminations, of different strength

**(a) The vacuum source floor (upstream's argument, reproduced).** `ρ_eff → ρ_crit/γ` is a constant, so L1's source never vanishes: `M_eff(<r) ~ (4/3)πr³ρ_crit/γ`, `V ∝ r` forever. Kills L1 for every γ > 0 and every ρ_crit > 0.

**(b) Ill-posedness, which is sharper and which the divergence obscures.** `ρ_crit = A·V_flat^B` is a **per-galaxy** constant. L1 therefore assigns the same point of empty space a *different* source density depending on which galaxy is being asked about — the field equation has no single-valued source. This is not divergence; it is ill-posedness, and **it survives into L2 = L3**, the reading that does carry the evidence. Contracting the same limit with `u_μu_ν` for dust carries it into §D.3: `𝒯_μν → (ρ_crit/γ)u_μu_ν` in vacuum, a nonzero stress–energy filling all space whose direction requires a preferred vacuum 4-velocity.

*(The upstream lane has (b) as open thread 3 and marks it "untouched by every reading above." Concretising it at L1 is this lane's addition; so is the §D.3 corollary.)*

## 4. Consequence for the whitepaper's own text

Executive summary and conclusion both carry, from S642 (GW170817): *"there is no Lagrangian, no action principle, no equation of motion for C(ρ)."* Against Appendix D these are three conjuncts of different truth value — and the upstream lane's flat verdict that *"Appendix D contains all three [field equation, action, covariant formulation]"* **over-corrects**:

| conjunct | verdict | basis |
|---|---|---|
| no field equation | **FALSE** | §D.2, §D.6.1, §D.7 state one |
| no action principle | **defensible, must be scoped** | §D.5's `S_eff` is a **worldline** action for a test particle's centre of mass, not a field action; **no action generates L1**. §5.15's L2 *is* action-derived — written 08-04 in this lane |
| no equation of motion for C(ρ) | **TRUE — load-bearing** | C is algebraic in local ρ in every reading; no kinetic term anywhere in Appendix D |

**The GW170817 conclusion survives on the third conjunct alone**: no propagating scalar mode exists to modify tensor dispersion. Correct phrasing is *"C has no kinetic term,"* not *"the framework has no action principle."*

## 5. Bookkeeping

**Refutation count HELD at 6. Bucket 0 = 0.** A *reading* was eliminated; nothing was refuted. This is the same bookkeeping as the 08-08 boost-ceiling closure.

## 6. Actions

| # | action | status |
|---|---|---|
| 1 | `[CORRECTION 2026-08-09]` block in §D.2 stating the inequivalence, the amplitude, the provenance and the a-priori kill | **APPLIED** |
| 2 | `[AMENDED 2026-08-09]` in whitepaper §5.15 correcting the 08-04 existence claim and its provenance | **APPLIED** |
| 3 | `[SCOPED 2026-08-09]` back-annotation on the S642 clause in executive summary **and** conclusion (parity) | **APPLIED** |
| 4 | **Adopt L2 in §D.2 and §D.6.1; drop `ρ/C` as a Poisson source** | **GATED ON DP** — physics decision. Nothing downstream depends on L1; §D.2's own second line is already L2's solution |
| 5 | Register the coupling fork's *argument* twin — C(ρ) vs C(g_bar) vs C(\|∇Φ\|) — and report the spread in dex, per the upstream lane's proposed registration gate | **OPEN, unrun** as of this writing; the acceleration branch is the frozen instrument (§5.15, amended 08-04), the density branch is `session128_dark_matter_derivation.py:54`, the `\|∇Φ\|` branch has no implementation found |

## 7. Failure class

This is `REC-2026-038`'s class — *check for an existing explanation before building a new one* — and **this instance is this lane's own**. The 08-04 amendment asserted *"a field equation does exist"* and *"one exists, and it is linear rather than nonlinear"* without grepping `manuscripts/`, a directory that is this lane's own working root. The general form, worth more than the instance:

> **An existence claim is a search claim.** "No X exists" and "one X exists" both assert the completeness of a search that is almost never stated. Name the surfaces searched, or write the sentence without the quantifier.
