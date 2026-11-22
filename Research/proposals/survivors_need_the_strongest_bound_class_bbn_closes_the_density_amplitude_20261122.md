# Survivors need the strongest bound class: BBN closes the density-keyed amplitude, single-photon interference closes the slow scan

*Site maintainer, 2026-11-22. Source: the site's simulated-visitor log for 2026-11-22 (researcher and graduate-student
personas), site explorer finding 2026-11-21 (Bell + LIV pincer), and one maintainer computation. No bucket moves.
Count 6, range 4–8, Bucket 0 = 0.*

## What the feedback revealed (WAKE)

Three of the site's "survives" or "below reach" verdicts were each judged against the bound the archive happened to
name, not against the strongest bound class in the literature. Two of the three fall when the stronger class is applied.
The third has not been computed. This is a research-process gap, not a site gap: every verdict of the form "survives",
"untested, not refuted" or "below reach" should have a named strongest-bound check before it is printed.

## (i) BBN sees the density-keyed absolute boost (computed)

The site and the 11-20 back-annotation say no local measurement sees the density-keyed law's absolute boost, because
lab G is measured at the same smoothed density (G_lab = G_bare/C(ρ_s)). That is true locally. But in the first minutes
ρ ≈ 2.5×10¹⁸ M☉/pc³, far above every knee on file, so C = 1 and G_BBN = G_bare. The law therefore predicts

    G_BBN / G_lab = C(ρ_local),    ρ_local ≈ 0.09 M☉/pc³ (smoothed ≥ pc)

Bound, retrieved 2026-11-22: Alvey, Sabti, Escudero & Fairbairn, EPJC 80, 148 (2020), arXiv:1910.10730,
"G_BBN/G_0 = 0.99^{+0.06}_{-0.05} at 2σ". (Copi, Davis & Krauss, PRL 92, 171301 (2004): 1.01 +0.20/−0.16 at 68%,
weaker.) Script `synchronism-site/maintainer/scripts/bbn_g_ratio_density_keyed.py`, predictions committed before the run
(site `42e99f2`); output beside it.

| Knee (M☉/pc³) | γ | C form | predicted G_BBN/G_lab | z (Alvey) |
|---|---|---|---|---|
| 0.42 (archive 0.028·V^0.5, MW) | 2 | floored | 0.57 | 16.9 |
| 0.42 | 2 | unfloored | 0.37 | 24.8 |
| 0.42 | 0.489 | floored | 0.38 | 24.4 |
| 0.42 | 0.489 | unfloored | 0.09 | 35.8 |
| 0.161 (velocity-blind, measured) | 2 | floored | 0.80 | 7.5 |
| 0.161 | 2 | unfloored | 0.71 | 11.2 |
| 0.161 | 0.489 | floored / unfloored | 0.46 / 0.21 | 21.1 / 31.0 |
| 0.0083 (RG elliptical) | 0.489 | floored / unfloored | 0.89 / 0.84 | 4.1 / 6.1 |
| 0.0083 | 2 | either | 1.000 | −0.4 (passes) |

Scorecard: P1, P2, P4 held as written. P3 held numerically (4.1–6.1σ, inside the registered 4–8σ) but its text called
the site's +9.4% the floored form; it is the unfloored one. The 2σ lower half-width is treated as Gaussian.

Consequences:
- A **third root** for the density-keyed amplitude, independent of SPARC, and independent of LLR's assumption about how
  density varies along the Sun's path.
- It excludes the **floored archive-knee law that the Oort vertical-gradient test lets through** (11-20).
- At Refracted Gravity's knee, γ = 0.489 is excluded at 4–6σ; γ = 2 is Newtonian there and passes.
- **Condition:** applies only when density is smoothed over a scale that puts the laboratory and the solar neighbourhood at
  one density, which is the reading the galaxy fits and the site's amplitudes use. Pointwise, lab G sits at Earth's density
  (C = 1) and this says nothing; the ephemeris argument is the right one there.
- Not yet done: RG *as published* at its three calibrations (its own ε function and smoothing length D; pairs with 11-21
  (ii), where only DMS joint survives LLR at D ≳ 30 pc). Queued for the site explorer.

## (ii) The slow-scan QM survivor fails single-photon interference (argued, not computed)

The site's Two Reframes page kept "T_scan ≥ 1.6–9.4 µs" as what survives the g²(0) check (explorer 2026-10-09). A scan
slower than the photon puts each photon in one arm for its whole life: a which-path variable, so averaged over photons the
state is a mixture with no relative phase between arms, and no fringes. Single photons interfere. Grangier, Roger & Aspect
(Europhys. Lett. 1, 173, 1986) is titled "…a new light on single-photon interferences" and reports a heralded single-photon
Mach–Zehnder (visibility near 98% from recall; publisher page refused retrieval). The two halves of one experiment close the
window from both sides. The R0 instantaneous-sampling escape faces the same check: fringes need relative phases between
modes, and a list of modes with relative phases is the amplitude vector. So the literal sequential-scan reading is excluded at
every period; the reading that keeps the phases is QM renamed (Bucket 3 in spirit; the site keeps the claim Speculative as an
interpretation). Archive docs that describe superposition as "temporal scanning" (SPINE's CRT section, FUNDAMENTALS
§Witnessing, the Observer Synchronization Framework) should carry this objection where they state the analogy. **Gates on the
project lead** for SPINE wording.

## (iii) Dimension-≥5 LIV: "below reach" was judged on time-of-flight only (not computed)

The B7 / Phase-12 "structurally untestable" verdict for the even-k quadratic term used time-of-flight. Threshold and
reaction-kinematics bounds (UHECR near GZK, PeV photons) reach much closer to the Planck scale for n = 2, because at
10²⁰ eV the lattice term k⁴a² exceeds m_p² by ~10⁶. Whether they bite depends on the grid's sign per species (a
nearest-neighbour lattice is subluminal for every field, which partly cancels in reactions). Literature named by the persona
from recall: Maccione, Taylor, Mattingly & Liberati 2009; LHAASO PeV-photon decay. Queued for the site explorer; until then
the site prints "none by time-of-flight; threshold bounds unchecked".

## (iv) The joint requirement, sharpened by the site explorer (2026-11-21)

The 11-21 stake ("no nearest-neighbour rule meets it without a protecting symmetry") was half right. A one-loop lattice
Yukawa computation with tree-level-equal speeds gives δc² = 0 (to 10⁻¹³) when time and space share one Euclidean spacing
and a finite, non-decoupling δc² otherwise (−0.27 and +2.46 × y²/16π² at ξ = 2, ½); a real-time tick has no symmetric point;
Edwards, Joó & Lin (PRD 78, 054501, 2008) show the same effect in lattice QCD (bare anisotropies 4.3 vs 3.4 for one
renormalized 3.5). Bell's 1986 lattice model (Phys. Rep. 137, 49) hides its frame in a beable layer the matter dynamics do
not feel. And there are three realist Bell escapes, not one (Wharton & Argaman 2020, RMP 92, 021002): foliation,
retrocausal, superdeterminism; superdeterminism works with a nearest-neighbour rule but does nothing for the loops.

**Proposed replacement for the 11-21 SPINE sentence (gates on the project lead):** *"A substrate can carry a hidden
preferred frame only if matter is not made of the thing that ticks: the frame must sit in a layer the matter loops cannot see.
FUNDAMENTALS makes matter the recurring pattern of the ticking field, so meeting both Bell and Lorentz tests would mean giving
up that identity, or tuning every particle's speed of light to one part in 10¹⁸–10³⁰."* This names the specific commitment
at stake: entity = recurring pattern of the substrate.

## Recommendation

A standing step: before any verdict reads "survives", "untested, not refuted" or "below reach", name the strongest bound class
for that observable and say whether it was checked. Candidates to run through it next: B7 Umklapp (threshold kinematics, as
in (iii)); the a₀(z) row (whether any high-z dynamical probe beyond RC100/Ciocan constrains n); the EFE = 0 statement
(already screened, 08-23). The site explorer can do this as a sweep; it is cheaper than discovering each one from a persona.

## Gates on the project lead
1. SPINE: the sharpened joint-requirement sentence (iv), replacing the 11-21 draft.
2. SPINE / FUNDAMENTALS: carry the mixture objection where the CRT scanning analogy is stated (ii).
3. Whether the density-keyed exclusion (now three roots: SPARC ΔBIC, LLR, BBN) becomes a booked result (ledger governance,
   carried from 11-12/11-20).
4. Carried: a dated committed configuration before Gaia DR4 (2 December, 10 days), or a statement that the framework is a
   family with no committed member (11-21 (iii)). The site now says "the framework is a family" in its top paragraph.
