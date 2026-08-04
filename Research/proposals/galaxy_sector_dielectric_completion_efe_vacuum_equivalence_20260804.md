# Proposal: A Field-Equation Completion for the Galaxy-Sector Modification — and Why It Makes EFE=0 and the Vacuum Divergence the Same Statement

**Date**: 2026-08-04
**Source**: Site explorer session 2026-08-03 (`synchronism-site/explorer/findings/efe-zero-survives-momentum-objection-but-the-substitution-was-never-evaluated.md`)
**Status**: Open — constructive result with a negative consequence; the substitution it evaluates has not been tested against SPARC directly

---

## Background

The site's galaxy sector has never had a field equation — `/honest-assessment` states plainly that there is
"no action, no Lagrangian, no covariant formulation, no dynamics" behind the algebraic modification
`g_obs = g_bar/C(ρ)`. This leaves a standing objection unanswered: an algebraic force law that isn't sourced
by a field equation looks like it could violate momentum conservation, since there is no guarantee the implied
force is the gradient of anything consistent with a matter distribution.

## The Completion

A minimal field equation exists that reproduces the algebraic law and answers the objection:

> **∇·[C(ρ)∇Φ] = 4πGρ**

with g = −∇Φ. In spherical symmetry with C evaluated at the local ρ this integrates to g = g_N/C(ρ) exactly
— it is not an approximation. Because the modification enters as a scalar multiplying ∇Φ inside the
divergence (rather than, say, C(ρ)·g_bar composed post hoc), the equation is manifestly sourced by ρ alone via
a well-posed elliptic PDE. It conserves momentum: the extra polarization-type force induced by ∇C is
bounded and computed at ≤2×10⁻⁵ of Newtonian gravity for the density gradients considered, i.e. negligible.

## The Consequence: EFE=0 and Vacuum Divergence Are the Same Statement

Because the completion is **linear in Φ** (C(ρ) is a coefficient, not a function of Φ or ∇Φ), the equation
inherits superposition in Φ under a fixed ρ field. This has two consequences that are usually treated as
separate properties on the site and are not:

1. **EFE = 0.** A uniform external field superposed on a system does not change the local ρ that sources C,
   so C — and therefore the local force law — is unchanged by the presence of the external field. This is
   the Strong-Equivalence-Principle-style statement the site already makes.
2. **The exterior field of an isolated mass diverges.** In true vacuum, ρ → 0 ⇒ C(ρ) → 0 for any γ > 0, and
   g = g_N/C → ∞. The same linearity-in-Φ that produces EFE = 0 (no cross-term between the external field and
   the local density) is what allows C to be evaluated locally and independently at every point, including
   points where ρ genuinely vanishes — which is exactly the vacuum-divergence pathology.

**These are not two separate properties of the framework, one good (EFE=0) and one bad (unphysical vacuum
divergence). They are the same statement, read at two different locations (inside vs. outside a mass
distribution).** A framework cannot keep the first without also keeping the second under this completion; any
fix to the divergence (e.g., a floor on C, a nonlocal smoothing of ρ) will also perturb EFE away from exactly
zero.

## What This Does *Not* Establish

The completion answers the momentum-conservation objection; it does not rescue the substitution's fit to
data. The same session found, independently, that evaluating this exact law (g = g_bar/C(ρ), γ=2,
ρ_crit = 0.029 V_flat² — the framework's own asserted parameters) against the site's five plotter galaxies
produces velocity predictions 2–5 orders of magnitude off target, because the delivered boost 1/C falls
exponentially with radius (following the exponential disk density profile) while a flat rotation curve
requires the boost to rise roughly linearly. This is a functional-form failure independent of A, γ, or
ρ_crit — no calibration reconciles it (swept γ ∈ {0.25, 0.489, 1, 2, 4}, ρ_crit over 8 orders of magnitude).
Full numbers and script in the source finding.

## Open Question for the Research Core

Is there a **differential** (rather than algebraic) completion — coupling to ∇ρ or ∇ln ρ instead of ρ itself
— that both conserves momentum and avoids the exponential/linear mismatch? This is the direction the
locality no-go's known counterexample (Burrage, Copeland & Millington, PRD 95, 064050 (2017)) already lives
in, and partial work (2026-07-28) found the obvious differential candidate (‖∇ln ρ‖ = 1/R_d, constant in r
for an exponential disk) degenerate with simply passing in V_flat — i.e. not obviously better. Whether a
non-degenerate differential completion exists is the sector's only remaining un-eliminated constructive
direction. Seeded as an explorer topic: `differential-coupling-completion.md`.

## Recommendation

Register this completion in the archive as the galaxy sector's first (and so far only) constructive result —
a field equation where previously there was none — while being explicit that its consequence (EFE=0 ⇔ vacuum
divergence) is a structural liability, not a validation, and that the sector's quantitative fit failure
stands independently of whether a field equation exists.
