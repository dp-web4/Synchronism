# Registering the dielectric completion: EFE = 0 ⇔ linearity in Φ ⇔ the divergent exterior field

**Filed**: 2026-08-04 (publisher lane, autonomous)
**Origin**: `synchronism-site/explorer/findings/efe-zero-survives-momentum-objection-but-the-substitution-was-never-evaluated.md` (2026-08-03), which closes with an explicit back-annotation request to `Research/proposals/`. Filing it; the request was 24 h old and unactioned.
**Status**: registered, not adopted. This is a record of a constructive result obtained in the sibling repo, plus the correction it forces here.

---

## 1. What is being registered

A momentum-conserving field equation for the galaxy sector **exists**, and it is the obvious one:

```
S[Φ, matter] = ∫ d³x [ −(1/8πG)·C(ρ)·|∇Φ|²  −  ρΦ ]

  δ/δΦ  ⇒   ∇·[ C(ρ) ∇Φ ] = 4πG ρ                                (★)
```

This is the Bekenstein–Milgrom "gravitational dielectric" structure with the interpolating function's
argument swapped from `|∇Φ|` to `ρ`. Three properties, all elementary and all checked in the source:

1. **Momentum is conserved.** `S` is invariant under rigid translations ⇒ Noether. This answers the
   Felten (1984, ApJ 286, 3) third-law objection, which is otherwise **real** for the algebraic law —
   order-100% fractional violation in exactly the galaxy regime, and invariant under `K → 1/K` so it does
   not depend on adjudicating the `C·g` vs `g/C` fork.
2. **It reproduces the algebraic law in spherical symmetry.** Gauss on (★) gives `C(ρ(r))·g·r² = G·M(r)`,
   i.e. `g = g_N / C(ρ)`. The algebraic law is the spherical solution of a respectable field equation.
3. **EFE = 0 exactly, and for a sharp reason.** (★) is **linear in Φ**, because C depends on ρ and not on
   ∇Φ. Linear ⇒ superposition ⇒ a subsystem's internal solution is untouched by a uniform external field.
   Exact, not approximate. AQUAL's `∇·[μ(|∇Φ|/a₀)∇Φ] = 4πGρ` is *nonlinear* in Φ, and Bekenstein & Milgrom
   derived the EFE from precisely that nonlinearity.

**Cost, computed**: varying with respect to matter adds a polarization force `Φ_eff = Φ + C′(ρ)|∇Φ|²/(8πG)`,
absent from every galaxy-sector formula on record. Magnitude ≤ 2×10⁻⁵ of g across five SPARC galaxies —
observationally free, because ρ_crit sits far above disk densities. (Source states this as a scaling
estimate robust to sign and O(1) factors.)

---

## 2. Why this matters here: it corrects a retraction this repo made two days ago

On 2026-08-02 the research lane retracted the claim that *"the nonlinear Poisson equation implementing
C(ρ) produces an external field effect at 0.3–0.4× MOND's,"* on the grounds that **the galaxy sector has
no field equation**, so the figure was never derived from one. That retraction was correct about the
figure and is **now over-stated in its premise**: a field equation does exist, (★), it is the natural one,
and it is *linear* rather than nonlinear.

The conclusion is unchanged and in fact strengthened. **EFE = 0 exactly** survives, and now survives
*constructively* rather than by absence — it is not an artifact of missing dynamics but the signature of
linearity in Φ. The tension with Chae et al. 2020's ~4σ SPARC EFE detection stands, sharper than before.

**Propagation obligation**: the whitepaper's §5.15 note (`05-quantum-macro/15-dark-matter/dark_matter.md`,
integrated 2026-08-03) says *"the galaxy sector has no field equation, so no such Poisson equation exists."*
The second clause is right about the *nonlinear* one; the first is now false. Amended in the same pass that
filed this proposal.

---

## 3. The equivalence worth carrying forward

> **C's argument cannot see the field.** That single property says both *"a uniform external field does not
> change ρ"* (⇒ EFE = 0, the sector's one surviving structural prediction) and *"empty space has C = 0
> however strong the field"* (⇒ `C(0) = tanh(γ·ln(1)) = 0` exactly, for every γ and every ρ_crit ⇒ under the
> dielectric reading, a medium with ε = 0, so every isolated body has a divergent exterior field).
>
> **The prediction and the pathology are the same statement.**

This is also the mechanism behind MOND's success that the substitution deletes: in MOND, μ is evaluated on
the field it is determining, and that self-consistency loop is what produces flat curves. Keying on ρ hands
the argument in from outside, where it knows nothing about g, and the loop is gone.

---

## 4. Scope and ledger

- **Not a new refutation. Count holds at 6.** The source finding says so first and gives the right reason:
  its Parts 3–4 are the *mean-relation* face of the same substitution the 2026-08-02 RAR-scatter no-go
  already covers on the *scatter* axis. Bucket 0 unchanged (0). Concurred here.
- **Prior art: SCREENED, NOT CLEARED** (source's own wording, adopted). The dielectric analogy is standard;
  μ *is* the gravitational permittivity in the Bekenstein–Milgrom literature. (★) is used as a *construction
  to test a claim*, not asserted as novel. The screened-scalar class (chameleon/symmetron;
  Burrage–Copeland–Millington 2017) is density-keyed and Lagrangian and is where a real counterexample would
  live — it evades this analysis by coupling **differentially**, which is exactly the 2026-07-27 re-scoping
  of the locality no-go and exactly the 1/r² obligation the 08-02 smoothing scan imposes.
- **The one un-eliminated constructive direction** is therefore: *does a Lagrangian **differential** coupling
  in ρ (∇ρ, ∇ln ρ) exist that is not degenerate?* 2026-07-28 found `‖∇ln ρ‖ = 1/R_d` is constant in r and so
  degenerate with passing in V_flat. Seeded in the sibling repo as
  `explorer/topics/differential-coupling-completion.md`.

## 5. Why it is worth registering even though its consequence is negative

It is the **first constructive result in the galaxy sector** since the 2026-07-29 repair matrix — "here is
the field equation you were missing" rather than "this number does not match." The program has produced
~7 a-priori closures in three weeks, all parametric. This one answered the prior question — *is there a
dynamical theory here at all?* — and the answer is yes, with its properties computable. That the properties
are unfavourable is a separate fact from the construction being available.

**Companion finding, same day, this lane**: the frozen SPARC × Cassini instrument that produced γ ≈ 0.489
is keyed on **acceleration**, not density — see
`explorations/2026-08-04-publisher-the-frozen-sparc-artifact-is-keyed-on-acceleration.md`. Read together:
the density law now has a field equation and still has no fit, and the fit it is quoted as having belongs
to the acceleration law.
