# Proposal: The galaxy sector's internal locality ledger — it is non-local by construction, in the wrong variables

**Raised:** 2026-09-07, maintainer track, from the site's visitor researcher + graduate-physics personas
**Status:** proposal — data-free structural claim, no execution required
**Relates to:** PREDICTIONS.md Bucket 2 (local-density no-go), `explorations/` Milgrom-non-locality
instance, `Research/OPEN_QUESTION_*` MRH formalization, the Fisher unidentifiability artifact

---

## The claim

The framework's galaxy sector is standardly defended as **local** — C is a function of local volumetric
density ρ(r) — and the local-density no-go is standardly framed as the collision between that locality
and the RAR's non-local organizing variable g_bar. That framing is **half wrong in the framework's own
favour**, and correcting it produces a refutation cheaper than anything currently on the ledger.

**The sector is not local.** Two of the mechanism's inputs come from outside any local neighborhood:

1. **ρ_crit = A·V_flat².** V_flat is the asymptotic rotation speed — defined in the r → ∞ limit, and via
   the BTFR fixed by the galaxy's *total* baryonic mass. So the quantity actually evaluated at radius r is
   **C(ρ(r), M_total)**, not C(ρ(r)). A threshold on a local field has been keyed to a global label.

2. **B_max = 1/Ω_m = 3.17.** A per-galaxy boost ceiling set by a *cosmological* parameter. A galaxy's
   relevancy horizon does not contain Ω_m.

Both are violations of the MRH's own **Predictive Closure** criterion ("adding elements outside the MRH
does not materially improve prediction; if it does, the MRH was incorrectly specified"). Both are stated
in FUNDAMENTALS/MRH as *definitional*, not as approximations, which is what makes this a self-consistency
failure rather than a modelling choice.

**Why it is worth registering.** It is **data-free**. It needs no SPARC fit, no BTFR slope, no ΔBIC, no
Cassini bound. It would stand unchanged if every empirical test on the Tier-1 ledger had *passed*. Every
other galaxy-sector refutation in the ledger is contingent on borrowed data; this one is not, which makes
it the cheapest entry we have and the only one that could have been written before any test was run.

## The mechanical consequence — this is the unidentifiability result

In the small-x regime SPARC actually samples (median x = ρ/ρ_crit ≈ 7×10⁻⁵):

    C ≈ γ·x = γ·ρ(r) / (A·V_flat²)

γ and A enter **only as the ratio γ/A**. That is one free number per galaxy, not two. This is exactly what
the Fisher correlation ρ(ln γ, ln A) = +1.000000 measures. The density-keyed unidentifiability artifact and
this locality violation are therefore **the same defect observed from two directions**: the model is fit
per-galaxy against a global label, so its local parameters are never separately measurable. The artifact is
currently derived from small-x linearization and presented as a statistical property; it is better
understood as a *structural* consequence of the non-local keying, and it is visible in the exact closed form
with no limit taken.

This also explains why the Refracted Gravity contrast case (Matsakos & Diaferio 2016) escapes: RG uses the
same ρ-keyed algebraic form but with a **universal** knee density, not a per-object one. Universal knee ⇒
knee lands inside the sampled range ⇒ two measurable numbers ⇒ identifiable. The identifiability is not
about the functional form at all; it is about whether the normalizing scale is imported from outside the
system. That is a sharper and more general statement of artifact 5 than the one we publish.

## The shape of the whole problem

    Needs one non-local variable it does not have:   g_bar (enclosed-mass acceleration) — the no-go
    Uses two non-local variables it should not have: V_flat, Ω_m                        — this proposal

The sector is **non-local by construction in precisely the place it claims locality, and non-local in the
wrong variables.** Stating it this way unifies what the ledger currently carries as two unrelated rows
(the Milgrom-non-locality instance; the Fisher degeneracy) and adds a third that costs nothing to verify.

## What this is not

- It is **not** a claim that the empirical refutations are wrong or superseded. They stand.
- It is **not** a claim that ρ_crit ∝ V_flat² is *numerically* wrong — that is the separate, already-registered
  sign result (MOND-matching forces ρ_crit ∝ V^−2; the framework asserts V^+2; 240×–3×10⁵× magnitude error).
  **The sign error is downstream.** The category error is that a local field was keyed to a global label at
  all, and it would remain even if the exponent were right.
- It is **not** an argument that the escape is closed. Two repairs are outside this objection and, as far as
  the maintainer can determine, **untested rather than refuted**: (a) MRH-smoothed density — replace ρ(r)
  with ρ averaged over an intrinsically-defined horizon; (b) a ∇ρ-keyed or gradient/symmetron-class coupling.
  Both are explicitly outside the local-density no-go's stated scope. Neither appears anywhere in the
  archive as attempted. *Untested ≠ refuted*, and the ledger should say which.

## Requested actions

1. Register the internal-locality violation as its own Bucket-2 row with a named refutation criterion:
   *"exhibit a galaxy-sector realization in which every input to C is available within the MRH of the point
   at which C is evaluated."* Currently no realization on record meets it.
2. Rewrite artifact 5 (density-keyed unidentifiability) to derive the degeneracy from the non-local keying
   rather than from small-x linearization — the linearization is a symptom.
3. Resolve the *untested vs refuted* ambiguity on the two escape routes above. One line each: "not
   attempted" or "attempted, result X." This program polices that distinction well everywhere else and
   its absence here reads as a closed door that was never opened.
4. Note in MRH/FUNDAMENTALS that the framework's flagship application violates its own foundational
   horizon criterion. That is the kind of thing SPINE exists to keep from oscillating away.

## Site status

Published 2026-09-07 on `/mrh` ("Does the framework respect its own horizon?"), with cross-references from
`/galaxy-plotter` and `/for-researchers` artifact 1.
