# Draft note to the Refracted Gravity authors: lunar laser ranging bounds the smoothing length (2026-11-27)

*Site maintainer, for the project lead. **Not sent.** Sending it is outward-facing and gates on dp. No bucket moves; count 6;
range 4–8; Bucket 0 = 0.*

## Why this exists

Every simulated-visitor researcher persona since 2026-11-17 has ended on the same point, and today's put it plainest: "the
binding constraint is no longer data. It is the twelve undecided forks, and the fact that no human outside the project has yet
been handed one of the citable results." Outreach has been routed to dp as a yes/no question on 11-08 (Saad & Ting), 11-17
(proposal `the_external_value_has_moved_to_other_peoples_theories_20261117.md`, decision 2), 11-21 (ii) and 11-24 (Chae et
al.). None has been acted on, and none came with a text to send. A decision is cheaper when the artifact already exists, so this
file is the artifact. dp can send it, edit it, or say no.

Why RG first, of the five "other people's theories" results:
- It is the one where the site's result is about **the authors' own theory as published**, at their own calibrations, with
  their own free parameter (the smoothing length D of Matsakos & Diaferio 2016 §2.2). The wide-binary and open-cluster results
  are about readings of MOND; the dSph result is about saturating boosts in general.
- It asks them a question only they can answer (what D is meant to be), so it invites a reply rather than a rebuttal.
- The researcher persona named it as "the cheapest route to the first external review".

Facts below are from `synchronism-site/explorer/findings/rg-as-published-has-a-smoothing-length-where-its-galaxy-fits-apply-llr-excludes-it-the-surviving-window-is-sub-parsec.md`
(PREREG site `4167dd5`, script `explorer/scripts/rg_smoothing_length_window.py`; 3 of 6 registered predictions held, 1 partial,
2 failed, all reported) and the 2026-09-23/10-27 LLR findings. Not re-run by me today. **Before sending, someone should re-run
the script** (it is minutes) and confirm the three readings of the permittivity below against the authors' text.

---

## The draft

> **Subject:** Lunar laser ranging and the smoothing length in Refracted Gravity: a question about D
>
> Dear Dr Matsakos, Prof. Diaferio and Dr Cesare,
>
> I am writing about a constraint on Refracted Gravity that we came across while auditing a different framework (Synchronism,
> https://synchronism-site.vercel.app/), whose density-keyed field equation turned out to be yours, ∇·[ε(ρ)∇Φ] = 4πGρ. We credit
> RG as the prior art. The analysis was done by AI research agents in our group and has not been reviewed by anyone outside it,
> which is the main reason for writing: we would like to know whether we have read your theory correctly.
>
> **The argument.** Matsakos & Diaferio (2016, §2.2) note that a density-dependent permittivity needs a smoothing length D,
> "tens of astronomical units or larger". If D ≳ 1 pc, every kernel centred in the Solar System has the same content to
> (r/D)² ≲ 10⁻⁷, so ε is uniform across the Solar System and G_eff = G/ε(ρ_s), with ρ_s the smoothed density around the Sun.
> The level of ε cancels against laboratory G. Its rate does not: the Sun moves through the disc and past nearby stars, so
> Ġ/G = −(d ln ε/d ln ρ_s)(d ln ρ_s/dt). Lunar laser ranging bounds |Ġ/G| ≤ 2.2×10⁻¹³ yr⁻¹ at 2σ (Hofmann & Müller 2018). We
> estimated d ln ρ_s/dt from Gaia DR3 neighbours for D ≥ 1 pc and from a hand list of nearby stars plus the Local Interstellar
> Cloud edge below that.
>
> **What we find,** for your three published calibrations (the DiskMass mean and joint fits of Cesare et al. 2020, and the
> elliptical E0 fit of Cesare et al. 2022), ratio to the LLR bound, ≥ 1 meaning excluded:
>
> | D (pc) | 0.1 | 0.3 | 0.5 | 1 | 3 | 10 | ≥ 30 |
> |---|---|---|---|---|---|---|---|
> | DMS mean | 5×10⁻⁶ | 0.07 | 130 | 7.8×10³ | 2.3×10³ | 62 | 1.2–9 |
> | DMS joint | 4×10⁻⁹ | 4×10⁻⁵ | 1.1 | 940 | 1.0×10³ | 30 | 0.25–6.2 |
> | E0 | 1.5 | 950 | 5×10⁵ | 7×10⁶ | 10⁶ | 2.5×10⁴ | ~10³ |
>
> (Permittivity written as in Cesare's 2024 review, Eq. 10; the ≥ 30 pc range spans smoothed solar densities 0.043–0.084
> M☉ pc⁻³.) So at the smoothing lengths where smooth-density galaxy fits make sense (a kernel holding ≳ 10 stars, D ≳ 2 pc near
> the Sun), lunar ranging excludes all three calibrations, with one conditional exception: the DMS joint fit at D ≳ 30 pc if the
> smoothed solar density is ≳ 0.06 M☉ pc⁻³. Below ~0.3 pc the local tests pass, but there ε follows the ISM phase structure
> (stars sit in ε ≈ 1 bubbles), and the smooth-density fits no longer compute the theory. The Solar System sets the other edge:
> D ≳ 10–70 AU depending on the kernel.
>
> **Our questions.**
> 1. Is D meant to be a fixed length, or is it tied to something (a kernel holding a given number of stars, the MRH of the
>    system being fitted)? The answer decides whether this is a constraint or a no-go.
> 2. We found three forms of ε(ρ) in circulation in our own code, which differ in whether the log is natural or base-10 and in
>    whether Q multiplies or divides. Which is yours? The exclusion holds for all three at 0.5–10 pc, but the windows differ.
> 3. Is there a reason the Solar System's ε should be evaluated differently (for example, a time-averaged ρ_s over a
>    relaxation time)? A relaxation time τ ≳ 80 Myr would evade the bound, but it is a different, orbit-averaged law.
>
> We would be glad of any correction. The scripts and pre-registrations are public in
> https://github.com/dp-web4/synchronism-site (`explorer/scripts/rg_smoothing_length_window.py`). The "cancels in the ratio,
> leaks into the rate" argument applies to any gravity law keyed on ambient density, so if it is right it may be of wider use.
>
> With thanks,
> [dp's name and affiliation]
> (Analysis: Claude, an AI research agent in the Synchronism project; one session, pre-registered, not externally reviewed.)

---

## Checks before sending (for whoever sends it)

1. Re-run `explorer/scripts/rg_smoothing_length_window.py` and confirm the table.
2. Confirm the three ε forms against Cesare 2024 Eq. 10 and Matsakos & Diaferio 2016 Eq. 4.1 (the site memory
   "grep every in-repo implementation of a rival equation" records that the forms diverged in this repo before).
3. Confirm the recipients' current addresses from their papers; do not guess.
4. Decide whether to attach the BBN result (G_BBN/G_lab, 2026-11-22). It is **not** included above: it was computed for the
   framework's tanh law at RG's knee value, not for RG's own ε, so it would overstate what we checked.

## Decision requested

Send (as is, or edited), or not. If not, say so, and the maintainer will stop routing outreach as an open item. Either answer
ends a loop that has run since 2026-11-08.
