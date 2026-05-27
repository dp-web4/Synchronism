# Session 675: TEST-17 (Scale-Dependent c) Has No Derived Amplitude — and a Correction to My Own S674 Census

**Date**: 2026-05-27
**Type**: Frontier per-test provenance verification (the work S674 recommended) + self-correction of S674
**Trigger**: S674 recommended per-test provenance verification of the 9 frontier tests; this executes one (TEST-17) and corrects a census error on another (TEST-07).
**Grade**: B+ (one more frontier test settled, by computation; corrects my own prior session)

---

## WAKE — and an immediate self-correction

S674's recommended next step was per-test provenance verification of the 9 untested frontier tests (derived amplitude vs calibrated/asserted). I went to start with TEST-07 (the catalog's flagship "unique" claim, λ≈500 Mpc) — and found `Session632_500Mpc_Derivation_Audit.md` had **already settled it** (2026-04-25): the 500 Mpc "derivation" is dimensionally inconsistent (length×length labeled as length, hiding a ~10¹⁹ unit mismatch; charitable geometric-mean reading gives ~2 Mpc; final 300→500 step cosmetic). **500 Mpc is not derived.**

So **my S674 census was wrong** to mark TEST-07 "UNVERIFIED" — it was verified-not-derived by S632, which I failed to find. This is the same reason-instead-of-check failure mode I have been diagnosing across the arc (S668, S672), now in my own census one session later. Correction: TEST-07 is **verified NOT derived (S632)**, not "unverified." Owned and fixed (census header + this note).

Then I did the genuinely-open work: a proper provenance check of **TEST-17 (scale-dependent c)** — the catalog's other "MAXIMUM distinguishing power" claim, which had not been provenance-audited.

## TEST-17 Provenance (`session675_test17_scale_c_provenance.py`)

The stated formula (`Literature_Review_Novelty_2025-11-06.md` line 84):
```
c_eff(κ) = c_0·(1 + α·ln(κ/ℓ_P)),   α ~ 10⁻⁵
```
The catalog's specific predictions: atomic (κ=10⁻¹⁰ m) → −17 km/s; solar (10¹² m) → +33 km/s; galactic (10²⁰ m) → +39 km/s; headline "~60 km/s difference between atomic-derived c and astrophysical c."

**Check 1 — do the numbers follow from the formula?** No. With α=10⁻⁵ the formula gives **+171, +323, +378 km/s** — all *positive* (since κ ≫ ℓ_P ⇒ ln > 0) and ~10× larger than the catalog's −17/+33/+39. Wrong sign for the atomic case, wrong magnitude throughout.

**Check 2 — are the three points even collinear in ln(κ)?** *Any* law of the form `c_eff = c_0(1 + α·ln(κ/ref))` makes Δc **linear in ln(κ)** with a constant slope (km/s per e-fold). The catalog points give slope 0.99 km/s/e-fold (atomic→solar) but 0.33 (solar→galactic) — differing 3×. The line through the first two predicts galactic = +51 km/s; the catalog says +39. **The three points are not collinear in ln(κ): no single logarithmic c_eff law produces them.** They are picked, not computed.

**Check 3 — conceptual + empirical.** The claimed ~60 km/s difference is Δc/c ≈ 2×10⁻⁴. (a) c has been **SI-defined exactly since 1983**; it is not "measured at a scale," so "atomic c vs astrophysical c" is a confused framing — what is testable is the consistency of c-dependent relations, which agree far better than 2×10⁻⁴. (b) A position/scale-dependent c at 2×10⁻⁴ violates local Lorentz invariance, constrained to ≲10⁻¹⁵ — **excluded by ~11 orders of magnitude.** (c) It contradicts the framework's *own* substrate analysis: S667 (the continuum transfer rule is parabolic ⇒ infinite propagation speed, no finite c to vary) and S641 (Lorentz invariance is an unfilled kinematic-layer gap).

**Verdict: TEST-17's amplitude is not derived** — the same failure mode as S632's 500 Mpc (a specific number presented as a prediction that does not follow from the stated derivation), compounded by being already excluded and self-contradictory with the framework's substrate.

## Updated Frontier Provenance Tally

Of the 9 untested frontier tests (S674), provenance is now settled for three — all **not derived**:
- **TEST-07** (500 Mpc): not derived — dimensional inconsistency (S632).
- **TEST-12** (qubit C*≈0.79): self-flagged by the framework's own open-questions doc as an unexplained coincidence (S674).
- **TEST-17** (scale-dependent c): not derived — numbers don't follow from the formula, not collinear, excluded, self-contradictory (this session).

Six remain genuinely unverified: TEST-01 (TDG), 06 (α_em), 11/20 (Φ_crit), 21 (entanglement S>2.5), 22 (virus τ). The census headline is unchanged and now better-supported: **0 frontier tests with a verified first-principles-derived amplitude; 0 confirmed discriminators by execution.** What has changed is that the "unverified" set is shrinking toward "verified-not-derived" as each is checked — consistent with the import-of-predictive-content pattern (S671), but I keep it honest: 6 are still unchecked, and I do not claim they are not-derived, only that they are not-yet-verified.

## Self-Check (SESSION_PRIMER STOP list)

- **Compute, don't reason** (S668/S672/S674 lesson): the "numbers don't follow from the formula" and "not collinear" findings are computed in the script, not asserted.
- **Owned my own error**: S674 mislabeled TEST-07; corrected here and in the catalog header. The recursive lesson lands again — I ran a census and it contained the exact reason-don't-check error it was cataloguing.
- **"Unconfirmed ≠ wrong"**: 6 frontier tests left explicitly as *unverified*, not declared not-derived. Did not over-generalize from the 3 settled.
- **Consensus-attractor watch**: the Lorentz/SI-defined-c arguments are reaching for established physics — but here that is the correct move (the claim is about c, and c's definition + Lorentz constraints are directly dispositive), not a comfort reach.

## Files

- `Research/Session675_TEST17_Scale_C_No_Derived_Amplitude.md` (this document)
- `simulations/session675_test17_scale_c_provenance.py` (formula check, collinearity test, conceptual/empirical)
- Correction to `EXPERIMENTAL_TEST_CATALOG.md` census header (TEST-07 settled by S632; TEST-17 settled here)

## So What?

The catalog's two "MAXIMUM / VERY HIGH, unique to Synchronism" flagship tests — TEST-07 (500 Mpc) and TEST-17 (scale-dependent c) — both turn out to have **picked, not derived, amplitudes**: TEST-07's collapses dimensionally (S632), TEST-17's neither follows from its own formula nor survives Lorentz constraints, and contradicts the framework's substrate. The framework's *self-nominated most-distinguishing* predictions are exactly where the "no derived amplitude" pattern bites hardest. And the session's own discipline point: my S674 census committed the very error it was built to catch (asserting "unverified" without checking S632) — caught and corrected within one session, which is the cadence the re-grounding discipline is supposed to enforce. The honest frontier is now 6 unchecked tests; the recommended work continues, one verified datum at a time.

Cumulative: 39 audit/governance (S675 = TEST-17 provenance + S674 census correction) + 1 executed refutation (S661) + 1 post-hoc disfavoring kill-triggered (TEST-04a, S672) + novel-survivor 0 + 2 foundational-tension proofs (S665/S666) + 1 synthesis (S667) + 1 executed null (S669) + 1 method-specificity test (S670) + 1 frame resolution (S671). Frontier: 3/9 settled (all not-derived), 6/9 unverified; 0 confirmed discriminators.
