# Proposal: the "EFE-free corner" is not a defined model, and declining to write a light-deflection rule does not protect the boost cap

**Date**: 2026-11-11 · **From**: site maintainer (Claude), from the 2026-11-11 simulated-visitor researcher and graduate personas
**Status**: for the project lead. No bucket moves. Count 6, range 4–8, Bucket 0 = 0.
**Script**: `synchronism-site/maintainer/scripts/efe_free_corner_which_acceleration.py` (+ `_PREREG.md` at site commit `3d53a42`, `_output.txt`).
**Blind status**: not blind. I did the hand estimate before registering; the script checks it with four controls (all pass). 6/6 predictions held.

## Why this is a research item and not a site fix

Two of the framework's surviving positions exist only because a rule has not been written. In both cases the visitor asked
what the rule is, and in both cases working through the options shows that every way of writing it is either excluded or
turns the position into something already on the ledger. That is different from "untested". It means the positions are
kept alive by the absence of a definition.

## 1. Which acceleration is the argument of C?

The site has carried an "EFE-free corner" since 2026-10-16: the acceleration-keyed compander
g_N = g·tanh(γ ln(1 + g/a₀′)) at γ ≈ 0.87–0.95, "read algebraically with no external field effect". It was scored against
perihelia, globular clusters, dwarf spheroidals and (prospectively) Gaia DR4 wide binaries. An algebraic law has to say what
g is. There are two readings.

**(T) g is each body's total acceleration.** This is the only reading a local law can evaluate without being told which
bodies form "the system". Executed:

- The Sun's Galactic acceleration (1.6–2.2×10⁻¹⁰ m/s²) sits in the transition, so it is boosted by ν_e = 1.26–1.44.
  Saturn's total Newtonian field is 6.5×10⁻⁵ m/s², where the boost is 1 to parts in 10⁹. Saturn is therefore left behind:
  a uniform relative acceleration (ν_e − 1)·g_N,gal = **4.6–4.9×10⁻¹¹ m/s²**, pointing within 6° of Saturn's apsidal line.
- Integrated numerically (control: zero force gives 0.009 mas/cy; the secular formula is reproduced to 0.1%), Saturn's
  perihelion turns by **−81,000 to −88,000 mas per century**, against EPM2011's −0.32 ± 0.47. That is 1.7–1.9×10⁵σ at
  every grid point. One in 10⁵ force directions would bring the perihelion rate inside 3σ, and there the eccentricity drift
  is largest instead.
- The same reading **has an external field effect**. For a wide binary the orientation-averaged boost is ν_e(1 + L_e/3),
  the same expression as the quasi-linear QUMOND result (control: at γ = ½ it returns 1.375–1.525 against the 11-06 row's
  1.38–1.53). At γ = 0.87–0.95 it is **1.17–1.30**, flat in internal acceleration.
- Applied to atoms instead of bodies, total accelerations are thermal and enormous, so there is no boost anywhere and
  galaxies are Newtonian. (This is the classic objection to pointwise algebraic MOND; I am citing it from memory, with no
  locator. Bekenstein & Milgrom 1984 is the field-equation fix. Not re-read today.)

**(I) g is the internal relative acceleration of a chosen system.** This is what every site number used: the boost for a
wide binary rises from 1.48 at g_N = 1.0×10⁻¹⁰ to 2.69 at 0.2×10⁻¹⁰ and never levels off. To evaluate it, the law must be
told that the binary is "the system" and the Galaxy is not, and that Saturn-plus-Sun is "the system" and the Galaxy is not.
No document in the archive writes that rule.

**What follows.**
- The corner is not "alive, disfavoured, not novel". It is **not a defined model**. Its one evaluable reading is excluded by
  five orders of magnitude; its other reading needs a partition rule.
- A partition rule is what the MRH is supposed to be. This is the second place the framework would need the MRH to carry a
  number or a procedure. The first was the density-keyed law's smoothing length (2026-09-15: "a missing declared
  parameter, and the first place the MRH would have to carry a number"). The γ ladder ("γ = 2 for every resolved-member
  system", 2026-09-08) is a third request of the same kind. All three ask: what decides where a system ends?
- I do not think this is an accident of three unrelated tests. The framework's distinctive concept is a context boundary,
  and each time a calculation needs that boundary, the calculation is done by a person or an agent choosing it by hand.
  A written rule would be falsifiable immediately (it has to give the planets no anomalous push and still boost a wide
  binary), which is probably why it is more useful than another sweep.

**Recommendation (gates on the project lead):**
1. Either state the system rule for the algebraic law, or retire the corner as a model and keep only the sentence
   "isolated MOND without an external field effect is disfavoured at about 3σ by globular clusters and testable by DR4",
   which is a statement about MOND, not about this framework.
2. If a rule is stated, the first three checks are fixed: planets (no uniform push above ~10⁻¹⁵ m/s², loosely), Earth–Moon
   (same, against lunar laser ranging), and whether a binary inside the Galaxy counts as isolated.
3. Decide whether "what decides where a system ends" should be an explicit open question in FUNDAMENTALS/STATUS, with the
   three instances listed. I think it should; it is the most concrete thing the MRH has been asked to do.

## 2. Which potential does light follow?

The 11-10 block made the KiDS lensing exclusion of the boost cap conditional on a light-deflection rule the framework has
not written. A persona asked today whether gravitational waves and light follow the same potential, and pointed to prior
art I had not connected. Working through the readings on the explorer's transcription of the eight KiDS-1000 bins
(Brouwer et al. 2021; reliable range g_bar > 10⁻¹³ m/s², their Sec. 5.4, read verbatim today):

| reading of light | framework's own cap 1/Ω_m = 3.17 | most permissive cap Ω_m/Ω_b = 6.39 |
|---|---|---|
| (a) light follows the boosted potential | exceeded in 3 reliable bins at ×5 baryons (4.1–8.9σ), in 1 at ×10 (4.2σ) | 1 bin at ×5 (4.1σ, at the reliability edge); passes every reliable bin at ×10 (explorer 11-10) |
| (b) light follows the visible matter only | lensing would show no excess; the three reliable bins below 10⁻¹² show 24–46×, **7.2–10.6σ above even ×10** | same |
| (c) unseen matter does the lensing | cap not tested by lensing; this is the dark-matter fork (09-22), and then the boost is not what holds galaxies up either | same |
| partial slip, lensing boost = B(1+η)/2 | needs η ≥ 1.5 at ×10, η ≥ 3.9 at ×5 (η = 1 is GR) | not computed |

So **for the framework as written, no reading of light passes the reliable KiDS bins without dark matter doing the
lensing.** The 11-10 softening ("the leg is a statement about completions that have one") was right for the 6.39 variant
and too generous for the 3.17 cap. The 6.39 variant is the one whose cosmology fails (10-27 (ii)).

**Gravitational waves.** Boran, Desai, Kahya & Woodard 2018 (PRD 97, 041501; abstract read today): GW170817 "rules out a
class of modified gravity theories, termed dark matter emulators", in which ordinary matter and light feel the potential a
dark halo would give while gravitational waves feel only the visible matter; the Shapiro delay difference through the
Milky Way is "about 400 days" (the persona recalled ~10³ days). The framework has no radiative sector that survives
(Bucket 2, spin-0), so it has no rule here either. But the constraint on any future rule is now two-sided: light must
follow the boosted potential (or dark matter must exist), and gravitational waves must follow the same one as light.
TEST-15 on the site registered a null "at 10⁻¹⁶" that can never fire; the question behind it was answered in 2018.

**Recommendation (gates on the project lead):** record "light and gravitational waves follow the potential that moves
matter" as a *requirement on any completion*, not as a prediction, and retire TEST-15 as answered by others. Whether
reading (c) is the framework's dark-matter story is the 09-22 fork, still open.

## 3. Outward-facing items, carried a fourth day

A researcher persona's top item today: the site states, on single AI sessions, that named third parties' code or fits have
problems (Saad & Ting's separation-error setting; Chae et al.'s consolidated Γ; Refracted Gravity's χ²/N; covariant RG
against Cassini). The site now labels each "not confirmed by the authors, who have not been contacted" and says the RG
numbers are this site's solver on SPARC, not a reproduction of the published DiskMass fits. The label is a stopgap. The
actual fix is to tell the authors, which only the project lead can authorise. A draft, to make that a two-minute decision:

> **Draft issue for the Saad & Ting repository (not sent).**
> Title: `run_gamma_analysis.py` line 106: separation error is max(5% of s, |line-of-sight distance error|)?
>
> In the public code at commit `4333e9f`, line 106 of `run_gamma_analysis.py` sets the uncertainty on the projected
> separation to the larger of 5% of s and the absolute line-of-sight distance error. For the 36-binary sample that gives
> a median of roughly 190% of s. The paper's text describes the separation uncertainty as about 0.1%. Is the line as
> intended? In mock catalogues we built with a known boost, this setting returned about 0.65 of the injected value, while
> a physical separation error returned it unbiased. We may well be misreading the code, and would be glad to be corrected.
> (This reading was made by an AI research assistant working on a public notebook-style project and has not been
> independently reviewed.)

If the lead prefers not to contact anyone, the alternative is to remove the claims from the public site until DR4 analyses
appear. Leaving them published and unsent is the one option I would argue against.

## What I am not claiming

- Reading (T) was never asserted by the archive, so its exclusion is not a refutation of the framework and does not enter
  the count.
- The lensing table uses the explorer's bin transcription (2026-10-04), not re-read today, and the ×5 / ×10 allowances are
  flat factors; Brouwer et al.'s hot-gas test is an isothermal halo with M_gas = M⋆ inside 100 kpc (their Sec. 5.4), a
  profile, not a factor.
- The planetary margin (10⁵σ) is against a published perihelion residual, not an ephemeris refit. A uniform force cannot be
  absorbed into a mass or a single orbital element, but I have not shown what a full refit would leave.
