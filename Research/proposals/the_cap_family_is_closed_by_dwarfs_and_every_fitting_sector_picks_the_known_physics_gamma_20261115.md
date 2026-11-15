# The cap family is closed by dwarfs, and every fitting sector picks the known-physics γ

*Site maintainer, 2026-11-15. For dp. Status: proposal (two decisions, one candidate sentence).*

## What changed today
- **The 11-14 proposal said the bounded-boost family survives at its most permissive cap, Ω_m/Ω_b = 6.39, except
  conditionally on lensing. That was wrong on the archive's own data.** The explorer's 2026-08-23 pressure-supported table
  had Draco (B_req ≈ 101) and Crater II (≈ 60) on file. It was tested at 3.17 and 13.7 only. Joined at 6.39
  (`synchronism-site/maintainer/scripts/cap6p39_vs_dsph.py`, PREREG `c7d729c`, 5/5 held, identity control exact):
  - 4.9σ each at M/L = 4, surviving a six-system trials correction;
  - no gas census and no light rule involved;
  - escape only by non-equilibrium (σ tidally inflated ≥ 1.3×) or M/L ≥ 8–12.
  - Draco excludes every cap up to ≈ 14.7 after trials.
- A researcher persona raised it. The persona's other claim, that "a weaker light rule makes lensing worse", was already
  in `LENSING_LEG_ASSUMPTION`.

## Decision 1: retire the cap family as a class
- **The only cosmic ratio above 14.7 on file is 1/Ω_b = 20.3, the no-dark-matter cap.** Dwarfs exclude it at fiducial
  M/L (Draco 4.2σ) but not at M/L = 4 (2.8σ).
- Every Ω_m-built cap (3.17, 5.4, 6.39) is excluded by matter dynamics without a light rule.
- **Recommendation:** state on the ledger that "a boost ceiling tied to a cosmic ratio that includes dark matter" is
  excluded, and stop evaluating members one at a time. The site has spent four sessions on 6.39.
  - Escape to name: non-equilibrium dwarfs.
  - Hypothetical member to name: a cap ≥ 15 derived from something other than Ω.
- This does not change the count (TEST-10 is already in the floor at 3.17). It closes a question that keeps reopening.

## Decision 2 (carried, day 8): the third-party code claim
The /wide-binaries text now says only that our reproduction differs from the paper on one input setting. The specifics
stay in the linked finding until the authors are consulted, because of a professional-courtesy point from a researcher
persona. The finding's filename still contains "bug". Contacting Saad & Ting is yours; I have not.

## Candidate surviving-content sentence (from a graduate persona's question 4)
> Wherever the framework fits data, it fits as a function physics already had. In galaxies every γ is a MOND interpolating
> function. Under the algebraic wiring, the fitted γ (0.489 ± 0.11) lands on the value where it *is* MOND's simple μ
> (γ = ½). In cosmology the fitted γ (0.487, mean-density reading) lands on the value where the sector *is* Λ. No sector
> has measured γ away from ½ at significance except by pinning or refuting it.

Caveat, checked today: under the AQUAL action (L3) γ = ½ gives x(x+3)/(x+2)², not μ_simple, and the L3 refit sits at
γ ≈ 0.44. So "picks the known-physics γ" is an L2 statement, while "fits only as a function physics already had" holds
under both wirings. Adopt, edit, or reject.

## Also recorded
- **Which field equation each test used** (SPARC/plotter: algebraic μ-form; TEST-25: QUMOND from μ = C; TEST-09/10:
  floored g_bar-keyed division). The L3 action has never met Cassini, so it is **untested** there, not passed.
- **Scope correction to 10-04 (ii):** dwarfs close 20.3 at fiducial M/L only.
