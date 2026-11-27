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

---

## Pre-send check (explorer, 2026-11-27)

*PREREG site `378f2f2` (before the script); script `synchronism-site/explorer/scripts/rg_note_presend_check.py`, output saved
beside it. Finding: `synchronism-site/explorer/findings/rg-note-presend-check-numbers-hold-q2-is-answered-2021-llr-closes-the-exception-kernel-shape-is-the-real-escape.md`.*

**Verdict: send an edited version (v2 below), not the draft as written.** Every number in the table survived. Three things
in the text needed to change, and one sentence needed adding.

### What held
- **Re-run.** `rg_smoothing_length_window.py` reproduces its saved output byte-for-byte. (The saved file also carries a
  hand-appended pointwise block that the script does not print.) Section A of the new script is an identity control: it
  reprints the note's table exactly.
- **Hand check, independent of the script.** DMS mean at D = 1 pc gives 7.7×10³ and DMS joint gives 920. The table says
  7.8×10³ and 940. α Cen alone (1.34 pc, −22 km/s) supplies d ln ρ_s/dt ≈ 10⁻⁵ yr⁻¹ at D = 1 pc.
- **ε form, read from LaTeX source.** Cesare et al. 2020 (arXiv:2003.07377) and the 2024 review (2404.06538, Eq. 10) both
  write `tanh[ln(ρ/ρ_c)^Q]`, which is Q·ln. The review also says "the larger its value, the steeper the transition", which
  rules out /Q. M&D 2016 writes `\log`.
- **Calibrations, read verbatim from the review.**
  - DMS mean (0.56, 0.92, −25.30) and DMS joint (0.661, 1.79, −24.54): review l.402.
  - E0 (0.089, 0.47, −24.25): review l.588.
- **LLR bound.** Hofmann & Müller 2018, (7.1 ± 7.6)×10⁻¹⁴ yr⁻¹, is correct.
- **The topic's stake.** Under a Gaussian kernel, every D = 1–10 pc cell stays ≥ 10× the bound: the smallest is 30×.
  - The star-count escape makes things worse. With a kernel holding the N nearest stars, the kernel's edge moves with the
    N-th neighbour. That gives 4×10³ to 5×10⁷×, sharp edge included.
- **Time-averaging.** No RG paper mentions averaging ρ over time; the only averaging is M&D's spatial D (grep of all three
  sources). So the note's question 3 stays a question.

### What changes
1. **Drop question 2.** The authors' papers already answer it (Q·ln). The other two forms were our own code's slips, and
   asking would tell them we had not read their papers.
2. **Use the newer LLR bound, and keep the old one in view.** Biskupek, Müller & Torre 2021 (arXiv:2012.12032, abstract)
   give (−5.0 ± 9.6)×10⁻¹⁵ yr⁻¹, so |Ġ/G| ≤ 2.4×10⁻¹⁴ at 2σ. That is 9.2× tighter than the 2018 bound.
   - The note's one conditional exception closes: DMS joint at D ≥ 30 pc goes to 2.3–5.2× at ρ_s = 0.084.
   - For a smooth kernel, all three calibrations are then excluded at every D ≥ 0.5 pc.
   - The passing window barely moves: D ≲ 0.31 pc (DMS mean), ≲ 0.43 pc (DMS joint), ≲ 0.06 pc (E0).
3. **State the kernel, and give the sharp-edged case.** This escape is new, and it is real. With a top-hat kernel, a star
   inside the kernel contributes a constant, so between boundary crossings only gas changes ρ_s.
   - At D = 1 pc the Sun's own top-hat term dominates ρ_s, so DMS joint passes the 2018 bound (0.03–0.9×).
   - At D = 10 pc both DMS calibrations straddle the 2018 bound (0.04–2.7×), and the chance of a star crossing the kernel
     edge inside the LLR span is 8–23 %.
   - D = 2–5 pc stays excluded under either bound.
   - So the maintainer's "every D = 1–10 pc cell ≥ 10×" is true for smooth kernels only. The authors would find this
     within a day. The note should say it first.
4. **Add one sentence on covariant RG** (Sanna, Matsakos & Diaferio 2023). The review itself points to the covariant theory's
   weak-field limit as the possible fix for ε(ρ), so it is the first reply to expect. Our 10-22 result is that the covariant
   scalar tracks the potential, not the density, so this argument does not carry over (the Sun's orbit gives Ġ/G ~ 10⁻¹⁵).
   Our separate concern about the covariant theory stays out of this note: one question per note.

Also checked: no paper applying LLR or Ġ/G to RG turned up (one search; not exhaustive). The Torino group is still testing RG:
Pizzuti, …, Diaferio 2025 (A&A 698, A83) on cluster kinematics. Recipient addresses are still to be taken from the authors'
current papers.

### v2 of the note (send-ready apart from name, affiliation and addresses)

> **Subject:** Lunar laser ranging and the smoothing length in Refracted Gravity: a question about D
>
> Dear Dr Matsakos, Prof. Diaferio and Dr Cesare,
>
> I am writing about a constraint on Refracted Gravity that we came across while auditing a different framework (Synchronism,
> https://synchronism-site.vercel.app/). Its density-keyed field equation turned out to be yours, ∇·[ε(ρ)∇Φ] = 4πGρ, and we
> credit RG as the prior art. The analysis was done by AI research agents in our group and has not been reviewed by anyone
> outside it. That is the main reason for writing: we would like to know whether we have read your theory correctly.
>
> **The argument.** Matsakos & Diaferio (2016, §2.2) note that a density-dependent permittivity needs a smoothing length D.
> If D ≳ 1 pc, every kernel centred in the Solar System has the same content to (r/D)² ≲ 10⁻⁷. So ε is uniform across the
> Solar System, and G_eff = G/ε(ρ_s), where ρ_s is the smoothed density around the Sun.
> - The level of ε cancels against laboratory G. Its rate does not.
> - The Sun moves past nearby stars and gas, so Ġ/G = −(d ln ε/d ln ρ_s)(d ln ρ_s/dt).
> - Lunar laser ranging gives Ġ/G = (−5.0 ± 9.6)×10⁻¹⁵ yr⁻¹ (Biskupek, Müller & Torre 2021), or (7.1 ± 7.6)×10⁻¹⁴
>   (Hofmann & Müller 2018).
> - We estimated d ln ρ_s/dt from Gaia DR3 neighbours for D ≥ 1 pc, and below that from nearby stars plus the edge of the
>   Local Interstellar Cloud.
>
> **What we find.** We used your Eq. (ε = ε₀ + (1−ε₀)½{1 + tanh[ln(ρ/ρ_c)^Q]}) at your three published calibrations: the
> DiskMass mean and joint fits (Cesare et al. 2020) and the E0 elliptical fit (Cesare et al. 2022). The table gives the
> ratio to the 2021 2σ bound, with a Gaussian kernel of width D; ≥ 1 means excluded. Divide by 9.2 for the 2018 bound.
>
> | D (pc) | 0.1 | 0.3 | 0.5 | 1 | 3 | 10 | ≥ 30 |
> |---|---|---|---|---|---|---|---|
> | DMS mean | 5×10⁻⁵ | 0.6 | 1.2×10³ | 7×10⁴ | 2×10⁴ | 570 | 11–85 |
> | DMS joint | 4×10⁻⁸ | 4×10⁻⁴ | 10 | 9×10³ | 10⁴ | 280 | 2.3–57 |
> | E0 | 14 | 9×10³ | 5×10⁶ | 6×10⁷ | 9×10⁶ | 2×10⁵ | ~10⁴ |
>
> The ≥ 30 pc range spans smoothed solar densities of 0.043–0.084 M☉ pc⁻³. Galaxy fits on a smooth density presuppose
> roughly ten or more stars per kernel, which means D ≳ 2 pc near the Sun. At those smoothing lengths, lunar ranging
> excludes all three calibrations. The local tests pass only below ~0.3 pc. There ε follows the phase structure of the
> interstellar medium (stars sit in ε ≈ 1 bubbles), and the smooth-density fits no longer compute the theory. The Solar
> System sets the lower edge, D ≳ 10–70 AU depending on the kernel.
>
> **Two things weaken this, which you may want to weigh.**
> - **Kernel shape.** With a sharp-edged (top-hat) kernel, a star inside contributes a constant, so only gas moves ρ_s
>   between crossings.
>   - D = 2–5 pc stays excluded.
>   - Near D = 1 pc the Sun's own term dominates ρ_s. The DiskMass joint fit then passes the 2018 bound (0.03–0.9×) and
>     straddles the 2021 bound (0.3–8×).
>   - Near D = 10 pc both DiskMass fits straddle the 2018 bound (0.04–2.7×) and sit at 0.4–25× the 2021 bound.
>   - Under this kernel a star crossing the kernel edge is a step in G, which the ranging record would see. For a 10 pc
>     kernel the chance of one during the record is about 8–23 %.
> - **Covariant RG.** This argument concerns the non-relativistic ε(ρ) used in the galaxy fits. In the covariant formulation
>   (Sanna et al. 2023) the scalar tracks the potential rather than the density, and the same estimate gives only
>   Ġ/G ~ 10⁻¹⁵ yr⁻¹.
>
> **Our questions.**
> 1. Is D meant to be a fixed length, or is it tied to something, such as the kernel's shape or a given number of stars?
>    (Tying it to a fixed number of stars makes the bound stronger, because the kernel edge then moves with the
>    outermost star.) The answer decides whether this is a constraint or a no-go.
> 2. Is there a reason the Solar System's ε should be evaluated differently, for example from a ρ_s averaged over time? A
>    relaxation time τ ≳ 80 Myr would evade the bound, but that is a different, orbit-averaged law, and we found no such
>    averaging in your papers.
>
> We would be glad of any correction. The scripts and pre-registrations are public at https://github.com/dp-web4/synchronism-site
> (`explorer/scripts/rg_smoothing_length_window.py`, `rg_note_presend_check.py`). The argument that the level of ε cancels
> against laboratory G while its rate does not applies to any gravity law keyed on ambient density. If it is right, it may be
> of wider use.
>
> With thanks,
> [dp's name and affiliation]
> (Analysis: Claude, an AI research agent in the Synchronism project; pre-registered, re-checked once, not externally reviewed.)
