# Proposal: the Cassini pass at γ ≳ 1.5–2 held a₀ fixed; TEST-25 leaves the count range; DR4 has no survival window

**From:** site maintainer, 2026-11-06. It started from a simulated-visitor researcher persona who asked for the Gaia DR4 γ window
to be pre-registered before December.
**Gates on the project lead:** (a) adopting the DR4 reading key; (b) the count-rule ruling, now with a floor of 4; (c) the stopping
question below.

## What was found
1. **One preparatory scan carried a survival window for seven weeks.** The scan is `simulations/sparc_cassini_q2.py`'s
   `target_scan`, which holds a₀ = 1.03×10⁻¹⁰ for every γ. Since 09-17 its output, "the compander passes Cassini at γ ≳ 1.5–2", has
   carried four things:
   - the 09-29 cascade ("TEST-25 needs point independence");
   - the 10-06 DR4 "one-γ squeeze with a narrow survival window";
   - the 10-24 range change (4–8 → 3–8);
   - six site pages.

   Re-profiling a₀ from the frozen SPARC surface at each γ gives Cassini failures of **14.5–19.3σ for every γ from 0.3 to 3**. The
   identity control at γ = 0.489 reproduces +17.96σ.
2. **Consequence for the count:** the range on file goes back to **4–8**. The floor is TEST-10, environment, TEST-25 and Bell/CHSH.
   Only TEST-09 and the γ = 2 pin move with the rule.
3. **Consequence for DR4:** no outcome lands on a surviving member. A reading key that maps boost to γ is committed before the data
   (`synchronism-site/maintainer/scripts/wb_dr4_reading_key_PREREG.md`). It uses the quasi-linear QUMOND point-mass boost
   ν_e(1 + L_e/3) in place of the 10-06 bracket, which was too wide to register.

## Why this matters beyond the number
This is the **under-refutation** direction: a survival window that was never there. It is the mirror of the over-refutations the
ledger has caught since July. Both came from the same failure: **a parameter held fixed in a scan that the downstream claim needed
refit**, the "fit value needs its partner parameter" pattern. The joint test's own design was right. It profiled a₀ per γ. The
window came from the side door, a scan labelled "not the registered test".

Recommendation: wherever a ledger sentence says "passes" or "fails" at some γ, record whether a₀ (or ρ_c) was profiled at that γ.
This is the a₀ twin of the 10-24 comparator recommendation (comparator, Δk, N, replication unit).

## The frame question (fourth time raised)
The persona's closing question was: "At what point does the program stop auditing the framework and declare the galaxy and cosmology
sectors closed?" Today's result removes the last named survival window in the acceleration-keyed galaxy sector. That sector is
excluded across γ 0.3–3 by SPARC × Cassini, and the density-keyed sector by SPARC (+2843) and LLR. Inherited or not, no sector has an
open window on a registered test.

The project lead should rule on the stopping question (10-04, 10-24) or name what would reopen it. Without a ruling, the daily loop
will keep auditing closed branches, and each audit finds bookkeeping, not physics.
