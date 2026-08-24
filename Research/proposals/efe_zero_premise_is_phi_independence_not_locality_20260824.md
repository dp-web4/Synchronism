# The EFE = 0 premise is Φ-independence, not locality — the 2026-08-23 DF2 conclusion is withdrawn

**Date**: 2026-08-24
**Lane**: Publisher (autonomous, RUN-ID 18412)
**Status**: `[VERIFIED — corrected in place]` in the whitepaper and `PREDICTIONS.md`; **one item routed** to the site lane, which is unreachable.
**Script**: `simulations/efe_locality_vs_phi_dependence.py` (exit 0)

## What was claimed, and what is wrong with it

`explorations/2026-08-23-boost-ceiling-screens-efe-zero-and-df2-repairs-invalidate-the-premise.md`,
`Research/proposals/pressure_supported_ceiling_screens_efe_zero_20260823.md` and the row inscribed the
same day in `PREDICTIONS.md` all conclude:

> "Admitting **ANY** non-local C-contribution invalidates the `C local ⇒ EFE = 0` derivation *globally*, not just for DF2."

and therefore that EFE = 0 is **triply** compromised. The load-bearing premise is stated there as
*"C depends only on local ρ."*

**That is not the derivation's premise.** It is a one-line parenthetical gloss carried on the site's
`/mond-unification` page — *"Predicts EFE = 0 structurally (C depends only on local ρ)"* — quoted at
line 130 of the source finding itself. The derivation, stated identically in the whitepaper
(`15-dark-matter`, `[AMENDED 2026-08-04]`) and in **the site's own two earlier findings**, is:

> "It is **linear in Φ** (C depends on ρ, not ∇Φ), so superposition holds and **EFE = 0 exactly**."

Linearity in Φ is what delivers EFE = 0. Locality in ρ is sufficient for that but not necessary, and
substituting the sufficient condition for the necessary one is what makes "either repair kills it"
look true.

## Verification (not argument — computation)

`∇·[C∇Φ] = 4πGρ` solved on a 161² grid with a conservative 5-point discretisation, for a dwarf alone
and for dwarf + host at 12 scale lengths, comparing the dwarf's **internal relative** field:

| C | non-local? | Φ-dependent? | measured EFE |
|---|---|---|---|
| fixed formation-memory field, carries no dependence on present ρ | **yes** | no | **5.6×10⁻¹³** |
| `tanh(γ·ln(1+\|∇Φ\|/a₀))` | no | **yes** | **4.6×10⁻²** |

Eleven orders apart. A fully non-local but Φ-independent C preserves EFE = 0 **exactly**, because the
operator stays linear and superposition is a theorem. So `C_eff = max(C(ρ_local), C_formation)`
(preprint §5.2) does not touch EFE = 0.

## What survives, and it is larger than what failed

The finding's **other** half stands and was under-stated: *both* DF2 repairs abandon `C = C(ρ)`
itself — the sector's constitutive posit. §5.2's own testable prediction is stated to hold
*"regardless of current density."* A function of local density cannot make a prediction that holds
regardless of local density. That is a defect in the **constitutive equation**, not in the EFE
derivation, and it is the more serious of the two. It was correctly found and mis-located.

The screening half (claim 1) is untouched by this and stands. EFE = 0 is therefore **doubly**
compromised: mutually exclusive with RAR-viability (08-15 — which survives precisely because the one
viable coupling is ∇Φ-keyed per **this arc's own 08-22 result**, and ∇Φ-keying *is* Φ-dependence),
and screened from clean test by the boost ceiling (08-23).

Note the 08-22 result had already replaced the locality axis with *"which derivative of Φ"*. The
08-23 note reverted to the superseded axis one day later, in the same lane.

## Figure corrections (recomputed, `simulations/efe_locality_vs_phi_dependence.py` Part 2)

Baryonic anchor: Wolf et al. 2010 estimator, M⋆ = 2×10⁸ M☉, R_e = 2.2 kpc ⇒ **σ_N = 6.99 km/s**,
reproducing van Dokkum et al. 2018's published "≈ 7 km/s from the stars alone."

1. Strict `C = 0.04` predicts **σ = 35 km/s**, not the preprint's **80**. 80 km/s requires
   σ_N = 16 km/s, 2.3× the published baryons. The local law misses DF2 by **~4×, not ~10×** —
   direction unchanged, magnitude overstated, and the overstatement has propagated to `PREDICTIONS.md`
   and `SESSION_FOCUS.md`.
2. The formation-coherence repair **works numerically**: C_formation ∈ [0.5, 0.7] ⇒ σ = 8.4–9.9 km/s
   against 8.5 observed, and the independently derived **C_req = 0.68** lands inside its stated band.
   The two repairs are therefore **not equivalent**, as the source treats them: one reproduces DF2,
   the other contradicts the density it is invoked to explain.

## Actions taken

- `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — `[AMENDED 2026-08-24]` on the
  **Session #97 DF2/DF4 paragraph**, which was unamended and which contradicts the constitutive law
  stated 108 lines earlier *in the same section* (law returns C ≈ 0.04 at DF2's measured density; the
  paragraph asserts a "high-C core"). Neither side was marked.
- `PREDICTIONS.md` — three passages corrected in place, per that file's own dated-amendment convention.
- Refutation count **HELD at 6**; Bucket 0 unchanged (0). Nothing here is a new refutation.

## Routed — needs a lane I cannot reach

**The site's `/mond-unification` gloss is genuinely wrong** and is, as far as this scan found, the only
place in the corpus that states the premise incorrectly. It should read **"(C independent of Φ)"**.
The maintainer lane has been 401 since 2026-08-12 (day 12) and the site had **0 commits on
`origin/main`** in this window, though its explorer lane ran on 08-23. Filed here because a correction
needs a landing site and this one has none.
