# Proposal: the count measures registrations, not framework content; the density fit is weak, not switched off

*Site maintainer, 2026-11-04. Sources: the site visitor log for 2026-11-04 (simulated graduate and researcher personas)
and site explorer 2026-11-02 (`synchronism-site/explorer/findings/density-law-knee-sits-below-every-disc-and-the-formula-A-fits-as-newton.md`,
pre-registered at site `c62bc9d`; I re-read both `_output.txt` files and did not re-run them). No bucket moves. The count
stays at 6, the range at 3–8, and Bucket 0 at 0.*

## 1. The strongest data exclusions of the archive's own equations are outside the count

A researcher persona read the TEST-26 card on the site and saw two dark-energy predictions that don't fit together. One
is the substituted family, which sits at ΛCDM's corner. The other is the Brans–Dicke completion (B), which forces
wₐ = +0.23…+0.60. The persona asked for the DR2 likelihood to be run on completion B. That run already exists. On
DESI DR2 + CMB + SN, completion B is Δχ² ≥ +79 behind ΛCDM at every ω. Completion A, which is Appendix D §D.3 solved on
FRW exactly as written, gives pure Einstein–de Sitter and χ² ≈ 9,900 (PREDICTIONS 2026-08-12 block).

So the program has an executed, data-facing exclusion of an equation the archive itself wrote down, at χ² ≈ 9,900. It
does not appear among the 6. The reason is the printed rule ("a stated failure condition, actually executed, fired as
registered"). Nobody wrote a failure condition for §D.3 before solving it. The same rule also excludes:

- the uniform G/C growth reading (10-12);
- the dark-energy sector run without CDM (09-22);
- the density law at its own formula's A. At A = 4.6×10⁻⁵ the best fit is Newton's, ΔBIC +6593 against the acceleration law. With γ = 2 at that A it is +10 372 (explorer 11-02).

**What this shows about the count.** The count measures how disciplined the registrations were. It does not measure how
much of the framework's content has met data. Those are different quantities, and the headline uses one to stand for
the other. Two kinds of reader get two different wrong pictures from this:

- A reader who takes "6" as a framework-content score is over-counted on rows #4 (an amplitude the equation never
  produced) and #6 (a simulation).
- The same reader is under-counted on the covariant equation, which failed worst of all.

**Recommendation (gates on the project lead).** Keep the 6 as the registration count. Beside it, print a second,
descriptive line that does not change the count: "executed exclusions of archive-written equations with no prior
failure condition: §D.3 as written (EdS, χ² ≈ 9,900), completion B (Δχ² ≥ +79), the density law at its formula's A
(ΔBIC +6593)". That would answer the persona's "N framework-specific / M convention-dependent / K inherited / 1
theorem-in-simulation" request without the site choosing a counting rule. The site's TEST-26 card now says why
completions A and B are outside the 6. The descriptive line itself has not been added; adding it is the ruling asked for here.

## 2. Back-annotation: the γ = 0.046 density fit has a weak density response, not a switched-off one

The 08-24 block's phrase "γ→0.046 = density-dependence switched OFF" is wrong, and so is the site gloss "at fixed
A = 0.029; C ≈ γx, only γ/A identified". The +2843 fit has one global knee, fitted at **ρ_c = 4.4×10⁻⁶ M☉/pc³**:

- That knee sits below 97% of SPARC disc points by more than 100×, so the discs are in the log regime, C ≈ γ ln(ρ/ρ_c).
  There γ and ρ_c are separately identified (Hessian correlation +0.89).
- The boost scales as ρ^−0.10 (5–95%: 0.06–0.17).
- It beats a constant boost by ΔBIC +945 (≈ 6σ at galaxy level). The response is weak, but it is not off.

The "only γ/A identified" property belongs to the *fixed-A* family (γ/A = 1.26×10⁷ at its optimum), not to this fit.
**The verdict is unchanged and stronger:** the density law loses at every calibration tried:

| Calibration | ΔBIC |
|---|---|
| Free global knee | +2843 |
| Free A | +3309 |
| Formula A, γ free (the fit becomes Newton) | +6593 |
| Formula A, γ = 2 | +10 372 |
| Published A, γ = 2 | +16 486 |

Open, cheap, seeded to the site explorer: is ρ^−0.10 just regression dilution of g_bar? Fit log B against log g_bar and
log ρ jointly.

## 3. Record (site only)

- The illustrative surfaces had drifted from their sources. The glossary's "Kill Criterion Triggered" examples were
  TEST-03 R² = 0.14 (misattributed) and TEST-04a ≈ 0.55 (superseded). The Test Roadmap read "FAILED —" on TEST-04a.
  The Tier 1 badge-reading example lagged its card. Lint rules were added, and the positive control fires 9 times on HEAD.
- TEST-25's card state now says which kill it means. The published, marginalized 8.7σ holds. The in-house empty
  intersection holds at point level only.
- TEST-02: the wide/close ratio is null, but absolute dynamics are not. At the published knee (unfloored) the boost at
  the Sun is ~3.5×10⁴ in g, and the Oort limit excludes it. The site no longer lets "identically null" read as "untested".
