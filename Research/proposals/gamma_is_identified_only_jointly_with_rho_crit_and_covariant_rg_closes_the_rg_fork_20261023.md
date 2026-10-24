# γ is identified only jointly with ρ_crit; covariant RG has no working dynamics (2026-10-23)

*Site maintainer, from simulated-visitor graduate and researcher personas (site visitor log 2026-10-23) and site explorer
2026-10-22. No bucket moves, count 6, Bucket 0 = 0.*

## 1. What γ does in C = tanh(γ ln(1+x)), x = ρ/ρ_crit

The archive and the site treat γ as a transition **sharpness** (Hill coefficient, "the γ = 2/√N_corr sharpness slot").
The visitor persona derived, and I re-checked numerically, that the function does not support that role:

- γ is the initial slope (C ≈ γx at small x) and half the saturation exponent (1 − C ≈ 2(1+x)^(−2γ)).
- On a log-density axis the steepest slope max dC/d ln x is bounded: 0.250, 0.322, 0.375, 0.408 at γ = ½, 1, 2, 4, and
  0.448 as γ → ∞ (max u sech²u at 2u tanh u = 1). The site printed 0.446 for the limit; that is the γ = 100 value.
- For γ ≳ 1 the change happens at small x, where ln(1+x) ≈ x, so C ≈ tanh(γx). This depends only on γρ/ρ_crit. Raising γ
  moves the knee to ≈ 0.77 ρ_crit/γ. It does not sharpen it.

**Identifiability, quantified** (`synchronism-site/maintainer/scripts/gamma_rho_crit_degeneracy.py`, estimate, not
pre-registered). Refitting ρ_crit for the alternative γ, the maximum |ΔC| over x ∈ [10⁻⁴, 10³] is:

| γ_ref | γ_alt | max \|ΔC\| |
|---|---|---|
| 1 | 2 | 0.036 |
| 2 | 4 | 0.020 |
| 4 | 16 | 0.015 |
| any | any, at x ≤ 10⁻³ | 2×10⁻⁷ (exact: only γ/ρ_crit enters) |

So a factor 16 in N_corr (γ = 4 → 1) is a 1.5% effect on C across seven decades of density once ρ_crit floats. And
ρ_crit does float: A in ρ_crit = A·V² is fitted and underived (S687, the 644× row).

## 2. Consequences

1. **γ = 2/√N_corr cannot be tested by a density-keyed fit unless ρ_crit is fixed independently.** This is a precondition
   on the Bucket 3 "N_corr-divergence test" note. N_corr can converge on one operational definition and still not be
   checkable against C(ρ), because what the data identify is a combination.
2. **This is the density-side twin of the a₀′/γ result** (2026-10-04 block): in the acceleration-keyed fits the identified
   scale is a₀′/γ; in the density-keyed form at SPARC densities (x ≈ 7×10⁻⁵) it is γ/ρ_crit. The free density fit's
   γ → 0.046 is therefore a statement about γ/ρ_crit, not about γ.
3. **The γ = 2/√N_corr sign-inversion argument survives in a narrower form.** It holds on a linear axis at fixed ρ_crit
   (mean-field finite-size scaling: rounding width ∝ N^(−1/2), so sharpness should grow as √N). It does not hold as "γ is a
   sharpness", because no γ makes a step on a log-density axis. That Bucket 2 row's executed content (r = +0.55 predicted,
   −0.55 measured) is unaffected; only the explanatory gloss changes.
4. **For the "step-function" reading in SPINE** (saturation makes emergence a step): the γ-family cannot produce a step on
   a log axis. Any step SPINE means must come from somewhere other than γ in this function. This is a reading note, not a
   refutation.

## 3. Covariant Refracted Gravity (back-annotating site explorer 2026-10-22)

The 2026-10-22 proposal (§2) and the 09-24 block of PREDICTIONS said the covariant RG completion's LLR status was
unchecked. The site explorer checked it (PREREG site `4e3898f` + `1244148`; not re-run by me). Sanna, Matsakos & Diaferio
2023 fix W = −1 and V = −Ξφ. The trace equation is exactly linear (□φ − 2Ξφ = 8πG T): ω = −1 Brans–Dicke with a 2.2 Gpc
range, so it cannot screen.
- Cassini: γ_PPN = 0, ~4×10⁴σ.
- LLR Nordtvedt: η = 1, ~2×10³σ.
- LLR Ġ/G: the paper's own FLRW system forces q − ζ = 2 and |Ġ/G| ≥ 0.99 H₀, ≥ 900σ at every admissible point.
- Weak field: φ varies ≤ 10⁻⁴ across a galaxy. CRG does not contain RG's ε(ρ).

So the 10-22 question "adopt RG's dynamics or keep the undeclared family" has no dynamics on the adopt side. RG is a
fitted galaxy-scale equation without a working covariant theory, and the 09-19 "a covariant completion exists" credit
reverses. **For the project lead:** the fork reduces to "build a covariant completion or keep the undeclared family";
and whether to tell the RG authors (outward-facing) is still open.

## 4. Also recorded (site, record-facing)

- The site called A2ACW's 47/3,308 (1.4%) survival rate "an upper bound on what in-distribution adversarial AI-AI
  collaboration can find". With calibration specificity 0/6 the filter removes real discoveries too, so the rate bounds
  nothing from either side until the post-cutoff arm runs. Corrected on the site, not in the archive (no archive
  instance found by grep).
- The site's refutation definition said "executed on external data against a registered criterion", which excluded
  the Bell/CHSH construction check and the γ = 2 pin (no TEST ID). It now matches the counting rule (stated failure
  condition, actually run, failed) and the footer's "5 on external data + 1 construction check" (cf. 10-07 (ii)).
- A visitor researcher persona asked why a 26-test registry with no win branch is still presented as a forward
  programme. This is the 10-04 stopping question from the reader's side. The enthusiast persona asked the same thing as
  "what would it take to declare a win?". The honest answer, "nothing currently registered", is not yet anyone's
  headline. **Gates on the project lead.**

## Recommendations (gate on the project lead)

1. Add the identifiability precondition to the Bucket 3 N_corr-divergence note.
2. Retire "sharpness" for γ in FUNDAMENTALS/STATUS wording, or state it as "sharpness on a linear axis at fixed ρ_crit".
3. Record RG's covariant result as part of the RG status line requested 10-22.
