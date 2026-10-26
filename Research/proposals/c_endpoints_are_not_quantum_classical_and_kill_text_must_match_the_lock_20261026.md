# C's endpoints are not quantum/classical; the TEST-26 kill text must match the sign lock; the plotter's "observed" points were hand-drawn

*Site maintainer, 2026-10-26. Source: four simulated-visitor personas (enthusiast, technical writer, graduate physics,
researcher) and site explorer 2026-10-24. No bucket moves. Count 6 (range on file 3–8). Bucket 0 = 0.*

## 1. The quantum/classical labelling of C is backwards (correction, archive and site)

The site's ρ_crit page opened "C transitions from quantum (C → 0) to classical (C → 1)" as density rises. /gamma-regimes
labelled the γ = 2 end "strongly quantum (single electron)" and γ → 0 "classical (macroscopic)". The archive README:161
carries "γ ~ 1 = quantum-classical boundary".

A graduate-physics persona pointed out that this is backwards against statistical mechanics. Quantum degeneracy sets in
at **high** phase-space density (nλ_dB³ ≳ 1), and a dilute gas is the classical limit. The most quantum many-body
states, condensates and superconductors, have N_corr ≫ 1 and land at γ → 0, the end the labelling calls classical. This
is the audited sign inversion of γ = 2/√N_corr seen from the other side. It is not a new defect, but it was never
connected to the endpoint labels. The enthusiast persona's "why do superconductors score low?" is the lay form of it.

**Done (site):** C runs "uncorrelated (C → 0) → collective (C → 1)" and is stated not to be a degeneracy parameter. The
γ bands are now "few-body / boundary / many-body collective", with a note on why the old labels were withdrawn and the
Curie–Weiss √N argument imported. Lint rules added.
**Done (archive):** bracket note under README:161.
**Recommendation (gates on the project lead):** FUNDAMENTALS/STATUS should say what C *is* in one sentence that does
not use "quantum" or "classical". Suggested: "C is a proposed collective-response fraction. It is not quantum coherence
and not a degeneracy parameter."

## 2. The TEST-26 kill text had not caught up with the 2026-09-23 sign-lock restatement

The sign lock forbids **thawing and crossing**: w₀ > −1 with wₐ < 0, crossed or not. That is the 09-23 block, with the
kill made ΛCDM's per the 09-27 correction. Three site sentences and the drafted registration on /top-5-tests still read
"fires if DR3 robustly requires the crossing." A thawing-only DR3 result would then kill the family on the page's own
logic without firing the written kill. For a pre-registered criterion, that is the defect that matters most, and it has
to be fixed before DR3.

**Done (site):** all four now read "DR3 excludes Λ in the thawing or crossing direction (w₀ > −1 with wₐ < 0, whether
or not w = −1 is crossed)". A lint rule retires the old wording.
**For the project lead:** if TEST-26 is registered, use that wording, and name the horn (mean-density vs local fluid).
This repeats the 09-23 recommendation. What is new is that the drafted registration text itself had the narrow
wording.

## 3. Record: the Galaxy Curve Plotter's "observed SPARC" points were hand-drawn for ~4.5 months

A graduate persona read the plotter's JS and found UGC 128 entered as an "LSB dwarf", V_flat 55 km/s, R_d 4 kpc. SPARC
(Lelli+2016 Table 1 and rotmod) gives a giant LSB disc: V_flat 129.3, R_disk 5.95 kpc. This is NGC 2403's
low-surface-brightness twin from the classic de Blok & McGaugh test. Checked against the in-repo SPARC files, all five
plotter rows were approximations:
- NGC 2403 at 3 kpc: plotter 118, SPARC 100.
- NGC 7331 V_flat: 250 vs 239.
- NGC 2403 R_d: 2.7 vs 1.39.

**Not a ledger issue.** Every ledger test reads the SPARC files directly. The plotter is a display, and no Bucket-2 row
rests on it.

**Done:** the rows are now real SPARC points with error bars. Max C_ρ on the shipped discs moves from "≈ 0.001" to
0.0004–0.008. The floor still binds everywhere, so the display verdict is unchanged.

**Method note:** this is the fourth "compiled number stood for months against its primary" instance (07-10, 10-09,
10-10, now), and the first in a display tool. The tell was a galaxy type ("dwarf") that contradicted the galaxy's name.

## 4. ρ_crit inventory drained (explorer 2026-10-24)

- **Site:** the six-object table is on /critical-density#one-symbol-several-knees, now retitled "the Reference Density
  (Saturation Knee)" with nav and meta matching the glossary.
- **Archive:** back-annotated the inverted (V/c)² ratio in `rho_crit_derivation_calibration_vs_prediction.md` and added
  a cosmological corollary to `one_scale_insufficiency_theorem_cluster_gap.md`.
- **Open question, untested:** is there a rule that predicts each sector's knee without fitting? If one exists, the
  unifying object is a₀, not ρ_crit.

## 5. Frame (for the project lead)

All four personas rated trust high. The two expert personas say the self-audit is more thorough than they would have
done. The researcher persona's verdict: "reject the physics; publish the A2ACW calibration (Youden J = 0, PPV 0/9,
cross-vendor agreement 0.41–0.52) as a short methods paper."

Today's HIGH items were all **propagation failures**: a correction reached one sentence and not its sibling. A visitor
loop that mainly finds propagation failures has reached its steady state as a physics critic. Its remaining value is as
a lint-generator. That is useful, but it is not discovery.

The two live, unknown-outcome experiments on file are:
- the cross-vendor A2ACW control (CBP 10-24 block); **[CBP-Claude 2026-10-26: label fix. The 10-24 block is not an A2ACW control. It is a cross-vendor *mutual-discrimination* run: measure how often agent a accepts vendor-b's correct vs wrong answers, then one critique-revise round vs voting on the same answers. It tests agent ensembles, not the debate protocol's novelty detection, so its outcome says nothing about A2ACW specificity. Template: `explorations/2026-10-16-interacting-ensemble-positive-control.md`.]**
- the post-cutoff arm (10-12).

Both are outward-facing, and both gate on the project lead. Recommendation: authorise one of them, or rule on the
stopping question (10-04/10-24). Without either, the loop will keep moving labels.
