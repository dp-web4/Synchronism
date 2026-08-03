# Proposal: Register the chemistry-sector orientation flip against C(ρ)'s monotonicity axiom

**Origin**: `synchronism-site/explorer/findings/two-coherence-orientations-chemistry-uses-the-flipped-one.md`
(2026-07-29, EXECUTED, parameter-free). Back-annotated 2026-08-03 by the maintainer track after
propagating the site-side fixes (`/coherence-function`, `/sound-velocity`, `/electronegativity`,
`/dark-matter-failure`, `/galaxy-plotter`).

## The result

`C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))` is monotone **increasing** in ρ for every γ, ρ_crit > 0 — this is
stated as a design axiom ("no paradoxical inversions") on the site's `/coherence-function` page and
is true by construction of the equation.

The framework's single strongest empirical correlation — chemistry sound velocity vs. C, r = 0.982,
the largest "validated"-class cohort on the site (1,703 phenomena) — runs the **opposite** direction.
Sound speed v = √(K/ρ) has ρ in the denominator by definition; on 22 real elemental solids,
Spearman(C(ρ), sound velocity) = **−0.32 for every (γ, ρ_crit) tested** (parameter-free — the sign
of a monotone transform's rank correlation cannot depend on its parameters). Diamond
(ρ = 3.51 g/cm³) ranks *more* coherent than lead (ρ = 11.34 g/cm³) in the site's own worked example —
anti-monotone in density.

Independently, the galaxy sector's missing-gravity term is coded as `v² = v_b² + (V_flat·C)²` —
proportional to C — while dark matter (the phenomenon that term is supposed to explain) is
identified elsewhere on the site as **low-C**. Same direction of error, different sector, checkable
with no data (`/dark-matter-failure`, `/galaxy-plotter`).

A repair-matrix check (in the finding above) shows four of six documented sign problems across the
program (BEC/BCS ordering, radial dC/dr<0, CFD viscosity/Bullet Cluster, chemistry ordering) collapse
to **one** orientation error on the "level" axis (low-ρ ↔ high-C should be swapped), not four
independent bugs. Two remain genuinely separate (the γ=2/√N_corr sharpness axis, and the ρ_crit
V-exponent location axis, the latter itself found estimator-dependent on 2026-07-29 in a related
session).

Flipping the orientation (`C → C(ρ_crit/ρ)`) makes the sign ledger self-consistent but was tested
and found **not sufficient**: on the two best-measured SPARC spirals it still loses to zero-parameter
MOND while spending two free parameters per galaxy. So this is not a rescue — it is a bookkeeping
correction with no predictive payoff, on top of an already-dead galaxy sector.

## Why this belongs in the ledger

This is the first finding in the program where the framework's **best empirical result and its worst
sector agree with each other** on a repair — every prior sign-finding could be read as "the failed
sectors are telling us something." This one says the success sector votes the same way. That is a
stronger, more specific signal than another isolated failure, and it reframes "chemistry validates
the framework" as itself contingent on which orientation convention you read the same number under.

## Recommendation

1. Record in PREDICTIONS.md (or STATUS.md, whichever tracks internal-consistency findings) that the
   chemistry sector's coherence ordering is opposite the governing equation's stated monotonicity
   axiom, with the Spearman −0.32 figure and its parameter-independence.
2. Do **not** treat the orientation flip as a rescue of the galaxy sector — it fails the same
   zero-parameter-MOND comparison it was tested against, just with a corrected sign.
3. Leave open: whether the flip survives in the cluster sector (untested), and whether "C" carries a
   third, still-different orientation in the consciousness sector (D/S form, flagged separately).
4. This does not change Bucket 0 (still 0) or move any existing refuted/reparametrization
   classification — it is an internal-consistency finding about the framework's own equation vs. its
   own flagship result, not a new external-data test.

## Site-side propagation (already done, 2026-08-03)

- `/coherence-function`: monotonicity axiom now carries the contradiction inline.
- `/sound-velocity`, `/electronegativity`: demoted from deprecated `validated` badge to
  `reparametrization`, with both the circularity (2026-05-06 finding) and sign (this finding) caveats.
- `/dark-matter-failure`, `/galaxy-plotter`: cross-linked note on the C-proportional boost term vs.
  low-C dark matter identification.
