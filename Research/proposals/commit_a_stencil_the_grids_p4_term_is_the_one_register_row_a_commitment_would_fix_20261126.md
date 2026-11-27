# Proposal: commit a stencil — the grid's p⁴ term is the one register row a commitment would turn into a bet

*Site maintainer, 2026-11-26. From a simulated-visitor researcher persona (site visitor log 2026-11-26, Pass 4),
a pre-registered maintainer run, and site explorer 2026-11-25. Count 6, range 4–8, Bucket 0 = 0; nothing here moves
a bucket.*

## What was found

1. **The grid's O(p⁴/M_Pl²) dispersion depends on direction at the archive's own tick.** This was pre-registered at site
   `823387d` before the script existed, and 4 of 5 predictions held. Script:
   `synchronism-site/maintainer/scripts/grid_dispersion_anisotropy_by_stencil.py` (+ `_PREREG.md`, `_output.txt`).
   Setup: a massless scalar, leapfrog with λ = c·tick/spacing, and any symmetric 27-point stencil.
   - The general form is η(n̂) = η₀ + η₄(Σn_i⁴ − 3/5), a cubic (ℓ = 4) pattern fixed in the substrate frame.
   - Face-neighbour stencil: η = (λ² − Σn_i⁴)/12. At its stability limit λ = 1/√3, that is −1/18 on an axis and 0 on
     a body diagonal.
   - The direction dependence can be removed by a stencil choice. Every such stencil then has η ≡ (λ² − 1)/12 < 0,
     and none is stable above λ = √3/2.
   - At FUNDAMENTALS' one cell per tick (λ = 1), all 651 stable stencils found depend on direction (|η₄| ≥ 1/24).
   - The miss (P4): no stable stencil in the class makes the p⁴ term vanish.
2. **The site's prior-art cell was wrong.** It cited causal sets for the dispersion, but causal sets are discreteness
   *without* Lorentz violation (Bombelli, Henson & Sorkin, gr-qc/0605006). Isotropic modified dispersion has LQG
   precedents (Gambini & Pullin, gr-qc/9809038, photons at linear order). The frame-locked ℓ = 4 pattern maps onto the
   anisotropic dimension-6 SME coefficients (Kostelecký & Mewes, arXiv:0905.0031). All three abstracts were checked
   on arXiv on 2026-11-26.
3. **The ontology differs from 't Hooft's CAI in two places, and both are costs** (explorer 2026-11-25). The stated
   rule is a contraction, and the phase-lock erases information. The particles are ontic lumps, not templates. See
   the PREDICTIONS 2026-11-26 block.

## Why this matters for the research direction

Across the register, every row that differs from MOND, ΛCDM or standard QM is excluded, belongs to someone else, or
is gated. This row is different in kind. Its magnitude (|η| ≈ 0.02–0.14) sits at the reach of ultra-high-energy
cosmic-ray threshold bounds. Its sign and its angular pattern are set by a choice the archive has never made. The
11-21 proposal said refutability needs a committed member; this is the cheapest place that applies. One sentence in
FUNDAMENTALS turns an uncommitted family into a number that can lose.

The frame question it raises: the anisotropy is forced by the archive's *own* commitment (one cell per tick), not by
an extra assumption. If the project lead keeps λ = 1, the framework predicts frame-locked cubic anisotropy at
p⁴/M_Pl². Whether existing bounds already exclude it is the open comparison. Either way, the verdict would belong to
the framework, not to someone else's theory. If λ < 1 with an isotropic stencil, the framework predicts isotropic subluminal η = (λ² − 1)/12 for that
species. Either way, the substrate stops being an interpretation in this one place.

## Recommendations (gate on the project lead)

1. **Commit an update rule:** λ and the stencil (or a stencil class) per field species. A FUNDAMENTALS §1 sentence is
   enough to make this row a Bucket 1 entry with a named refutation criterion.
2. Then compute direction-dependent reaction thresholds (photon decay, vacuum Čerenkov, GZK) for the committed
   stencils and compare them with the anisotropic d = 6 SME bounds and the Auger arrival-direction data. This is
   explorer work and is seeded as a topic.
3. **SPINE "The one move":** name CAI and the two differences. Offer the two readings: (a) the move is CAI plus a
   Klein–Gordon lattice rule, inheriting its open problems; (b) the move rests on ontic pattern-particles, whose
   single-particle interference is an unpaid burden.

## Caveats

- Scalar field and leapfrog only, nearest 27-point stencils only, O(k⁴) only. Fermions and gauge fields on a lattice
  have their own discretizations (doublers, Wilson terms), and their angular patterns differ.
- The λ = 1 statement comes from a numerical scan of the weight plane, not a proof.
- Nothing here is an external review.
