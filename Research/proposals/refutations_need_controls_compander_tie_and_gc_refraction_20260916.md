# Proposal: two refutations lacked the controls we demand of claims. The compander ties MOND once only the function is swapped, and the GC window omits refraction of the Galactic field.

**From:** site maintainer, 2026-09-16 (WAKE, written before the rest of the session's fixes were pushed)
**Trigger:** site visitor log 2026-09-16, graduate-student and researcher personas, independently
**Count:** stays 6. Bucket 0 stays 0. **No bucket moves.** One dp-gated recommendation is withdrawn.

## Why this is research, not site polish

The visitor researcher's verdict was: *"The live risk is no longer overclaiming. It's that the refutations are less
guarded than the claims were."* Today two of those refutations were checked, and both needed a qualifier. One of them
underpins a recommendation routed to dp (TEST-25 reclassification, proposal of 2026-09-10, item 4).

## 1. "The compander form fails on its own, 2.10×": withdrawn

**Claim (explorer 2026-09-09; site /parameter-derivations item 2; proposal 20260910 item 4).**
- Same argument (|∇Φ|), same floor, same 153 SPARC discs, same field-equation solver, "swapping only the function".
- MOND μ scores χ²/N 51.45; the tanh-log compander scores 108.10: 2.10×.

**What the source actually varied** (`synchronism-site/explorer/findings/scripts/ceiling_vs_likelihood_where_does_the_boost_deficit_hide.py`):
- the function;
- the knee: 0.32 a₀, where the exact γ = ½ identity tanh(½ ln(1+x)) = μ_simple(x/2) needs 0.5 a₀;
- the floor **form**: MOND μ used a clip, max(μ, 0.089); the compander used an affine floor, 0.089 + 0.911·tanh(…),
  which removes boost at every radius.

**Pre-registered controls** (site commit `04c4037`, before running;
`synchronism-site/maintainer/scripts/compander_form_isolation_controls.py` + `_output.txt`):

| run | function | γ | knee | floor | χ²/N | ×R0 |
|---|---|---|---|---|---|---|
| R0 | MOND μ | – | a₀ | clip | 51.45 | 1.000 |
| R1 | tanh-log | 0.5 | a₀/2 | clip | 51.45 | 1.000 (identity control; algebraic max diff 1.1×10⁻¹⁶) |
| **R2** | **tanh-log** | **0.489** | **a₀/2** | **clip** | **50.83** | **0.988** |
| R3 | tanh-log | 0.489 | 0.32 a₀ | clip | 77.98 | 1.516 |
| R4 | tanh-log | 0.5 | a₀/2 | affine | 67.69 | 1.316 |
| R5 | tanh-log | 0.489 | 0.32 a₀ | affine | 108.10 | 2.101 (source reproduced) |

- Prediction 1 (R1 = R0) held.
- Prediction 2 (R2/R0 < 1.15) held: the result is 0.988.
- Prediction 3 (the gap is carried by the knee and/or the floor form) held for both.

**Consequences.**
- The compander stays **Reparametrization**, now without the "and it loses" clause. Swapped honestly, it ties MOND's μ,
  which is what an exact identity predicts.
- **Proposal 20260910 item 4 ("reclassify TEST-25 as framework-specific; the 2.10× closes MOND's escape hatch") —
  recommendation withdrawn.**
  - What survives: the compander's only escape from Cassini is raising γ (return exponent q = 2γ), and TEST-25's own
    γ scan already shows the SPARC-retained interval fails.
  - The asymmetry premise is weaker than stated. Desmond, Hees & Famaey 2024 report the tension persists across all
    interpolating-function families they tested, so MOND is not obviously free to walk away either.
  - Recommendation now: **keep "inherited from MOND", question open.**
- The 51.45 baseline is MOND in that pipeline (algebraic MOND scores 52.21 there). The visitor's alternative reading —
  that the right reference is the 21.2 quoted elsewhere — is a cross-pipeline comparison and is **not** adopted.

## 2. "A density-keyed law has no external field to appeal to": true only algebraically

- The GC window (ρ_c ∈ 0.1–300 M☉/pc³ excluded at γ = 0.489) is computed as g = g_N/C(ρ). For an isolated spherical
  cluster that is exact under L2, ∇·[C∇Φ] = 4πGρ, by flux conservation.
- But L2 is linear in Φ, so the Galactic field **superposes**. Inside a body whose C varies, it is refracted and is
  **non-uniform**. The result is a differential acceleration linear in g_ext: not MOND's nonlinear EFE, and not zero.

**Order-of-magnitude check** (not a pre-registered test; expectation written first):
`synchronism-site/maintainer/scripts/gc_external_field_refraction_estimate.py` + `_output.txt`.
- Model: exact l = 1 mode, Plummer 2×10⁵ M☉, a = 3 pc, floor 0.315, g_ext = 1.57×10⁻¹⁰ m/s².
- Beyond the knee, the non-uniform residual is **0.29–0.37 g_ext** at every tested (γ, ρ_c).
- Its ratio to the cluster's own gravity is **0.25–0.67 at 20–30 pc** and **>1 by ~35–40 pc** (Newtonian Jacobi radius
  there ~80 pc).
- The term is dipolar, so it cancels at first order in angle-averaged dispersion profiles. The slope statistic may
  survive, but the untested signatures are anisotropy and truncation.

**Consequences.**
- "EFE = 0 exactly" (and PREDICTIONS' "EFE=0 is the observable signature of locality" thread) describes the
  **algebraic** reading. For L2 with ∇C ≠ 0 inside the system, the correct statement is **"no nonlinear EFE; linear
  refraction of g_ext present."**
- `simulations/efe_locality_vs_phi_dependence.py` (2026-08-24; EFE = 5.6×10⁻¹³ for a Φ-independent C) should be checked
  for *what it subtracted*. If the refracted external solution was removed, the zero restates superposition. It would
  not show that internal relative accelerations are absent.
- **A candidate framework-independent observable for the whole ε(ρ) class (Refracted Gravity included):** GC truncation
  radii scaling linearly with g_ext rather than as Jacobi's (M/M_gal)^{1/3} R_GC, with a dipolar outer-halo asymmetry
  aligned with the Galactic field. Seeded to the site explorer: `explorer/topics/refraction-efe-linear-permittivity-gravity-gc-truncation.md`.

## 3. Provenance (explorer 2026-09-15 P0, propagated)

(ε₀, Q, ρ_c) = (0.089, 0.47, 0.0083 M☉/pc³) is Refracted Gravity's **elliptical (E0)** mean, Cesare et al. 2022
(arXiv:2102.12499). It is not the Cesare+2020 DiskMass disc fit it was attributed to.
- Disc values: 0.56 / 0.92 / 7.4×10⁻⁴ (mean) and 0.661 / 1.79 / 4.3×10⁻³ (joint).
- Both disc knees also sit below the GC exclusion band, so no verdict changes.
- The archive proposals 20260909 and 20260915 use 0.0083 as "RG's knee"; read it as E0.

## Decisions for dp

1. Note the withdrawal of 20260910 item 4's recommendation (no action needed if it was not yet ruled).
2. Whether density-keyed per-object tests must declare, besides γ and D (09-15), **algebraic vs field-equation
   (L2/L3) treatment of external fields**.

## Method lesson

"Swap only X" needs a written list of everything else that differs between the two code paths, checked against the code
and not the prose. Here the floor had the same *value* in both runs but a different *form*. A value match hid a
functional-form mismatch, the same shape as the site's two C(ρ) floored/unfloored forms.
