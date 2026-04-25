# Session 632: 500 Mpc Derivation Audit — Same Failure Mode as S631

**Date**: 2026-04-25
**Type**: Targeted archive audit responding to 2026-04-25 back-annotation
**Grade**: A- (clean disposition, dimensionally rigorous)

---

## Trigger

Site maintainer's back-annotation `Research/proposals/cosmic_interference_500mpc_derivation.md` (2026-04-25), responding to a four-persona visitor review:

> "If derivable from a coherence wavelength formula, this is the test that matters most for my field."

The site's TEST-07 claims oscillations at λ ~ 500 Mpc, links to a 404 page, and provides no derivation. The proposal asks: where does 500 Mpc come from?

## Archive Search Result

The derivation lives in `Research/Cosmic_Interference_Search_Protocol.md`, Session #4 Track C, 2025-11-08. The chain:

1. Cluster mass M ~ 10¹⁵ M☉
2. R_MRH ~ (GM/c²) · (c/H₀)
3. With numerical evaluation:
   ```
   R_MRH ~ (6.7×10⁻¹¹ · 2×10⁴⁵ / 9×10¹⁶) · (3×10⁸ / 2.3×10⁻¹⁸)
   R_MRH ~ 1.5×10¹⁸ · 1.3×10²⁶ m
   R_MRH ~ 2×10⁴⁴ m ~ 600 Mpc
   ```
4. Then λ_interference ~ R_MRH/2 ~ 300 Mpc
5. Then "in terms of observable structure: λ_obs ~ 500 Mpc"

## Dimensional Audit

The formula R_MRH = (GM/c²)·(c/H₀) is **length × length**. Its result has units of **m²** (area), not m (length).

Numerical check:
- GM/c² with M = 10¹⁵ M☉: half the Schwarzschild radius ≈ 1.5 × 10¹⁸ m ≈ 0.05 kpc
- c/H₀: Hubble length ≈ 1.3 × 10²⁶ m ≈ 4200 Mpc
- Product: 1.5×10¹⁸ × 1.3×10²⁶ = **2×10⁴⁴ m²**, **not** 2×10⁴⁴ m

The text labels this product "2×10⁴⁴ m ~ 600 Mpc". But 600 Mpc = 1.85×10²⁵ m. The "m" label hides a unit mismatch by ~10¹⁹.

If we charitably interpret the intended formula as a geometric mean (which would have units of length):
- √(r_s · R_H) for r_s = 3×10¹⁸ m, R_H = 1.3×10²⁶ m
- = √(3.9×10⁴⁴ m²) = 6.2×10²² m ≈ **2 Mpc**, not 500 Mpc

To make √(r_s · R_H) ~ 500 Mpc, we'd need r_s · R_H = (1.5×10²⁵ m)² = 2.4×10⁵⁰ m², requiring r_s = 1.85×10²⁴ m for the Hubble length we have. That r_s corresponds to M ≈ 6×10¹⁶ M☉ — about 60× the most massive known clusters.

No honest dimensional reading of the stated formula produces 500 Mpc for cluster masses.

## The R_MRH/2 → 500 Mpc Step

After the dimensionally inconsistent "R_MRH ~ 600 Mpc" claim:
- "λ_interference ~ R_MRH/2 ~ 300 Mpc"
- "in terms of observable structure: λ_obs ~ 500 Mpc"

The factor R_MRH/2 has no derivation. The leap from 300 Mpc to 500 Mpc has no explanation. Both factors are cosmetic adjustments.

## Verdict

**500 Mpc is not derived from the framework.** The Session #4 Track C "derivation" is:
- Dimensionally inconsistent (length × length labeled as length)
- Even with charitable correction (geometric mean), gives ~2 Mpc, not 500 Mpc
- Concludes via two unjustified factor adjustments (R_MRH → R_MRH/2 → 500/300 × 300 Mpc)

The number 500 Mpc was chosen, then a derivation chain was written that nominally connects to framework parameters but does not survive dimensional scrutiny.

## What This Adds

**Three site-derivation audits, all same failure mode:**

| Site claim | Archive source | Failure |
|------------|---------------|---------|
| BTFR n ≈ 2.2 (S631) | Session #48 Track B | Self-labeled "not rigorous"; observed n = 3.85 ± 0.09 |
| A = 4π/(α²GR₀²) (S631) | Session #66 | α = 1.0 (fiducial), not fine-structure constant |
| TEST-07 λ ~ 500 Mpc (S632) | Session #4 Track C | Dimensionally inconsistent; cosmetic factor adjustments |

S631 framed this as "public-claim failure" — distinct from S617–628's internal-physics failures. S632 strengthens that frame: a third site claim with archive content that doesn't support the public framing. The pattern is now empirically established, not anecdotal.

The visitor researcher persona has now triggered two productive audits in three days. The external-feedback channel works.

## What I Tried Before Concluding

- Searched archive for "500 Mpc", "coherence wavelength", "cosmic interference" — found Session #4 Track C and Literature_Review_Novelty_2025-11-06.md as primary sources.
- Attempted to verify the dimensional formula at face value: fails.
- Attempted geometric-mean reinterpretation (a charitable reading): gives ~2 Mpc, not 500.
- Tried the proposal's Approach 2 (ρ_crit = A·V² applied to clusters with σ ~ 1000 km/s): A·V² is itself dimensionless under the current parameter conventions (A in (km/s)⁻², V in km/s), so it produces no length scale. This route doesn't help either.

## Recommended Site-Side Actions (Not Taken — Site Is Reference-Only)

1. **TEST-07 → relabel "Speculative" or "Exploratory"** with explicit note that the 500 Mpc scale is not derived from framework parameters and the Session #4 Track C dimensional argument is inconsistent.
2. **Remove the 404 link** to /cosmic-interference, or replace it with an archive pointer to the (corrected) failure-mode documentation.
3. **Asymmetric falsifier label**: as the visitor noted, "any oscillation anywhere" is not a sharp prediction. Without a derived scale, this test does not currently distinguish Synchronism from a generic "something at sub-BAO scales" claim.

## Connection to Other Open Issues (from the proposal)

The proposal explicitly noted:
- BAO ~10⁻⁴ shift (TEST-04) has the same derivation gap
- Wide-binary density coefficient (TEST-02) lacks a predicted magnitude

If the operator wants to extend this audit, those are the next two targets. The methodology is reproducible: archive search → dimensional check → disposition (derived | semi-empirical | undeerived).

## So What?

S631 found two site claims unsupported by archive content. S632 finds a third. The external-feedback channel reliably surfaces these in days; internal review across 600+ sessions surfaced none of them. The governance gap is structural, not anecdotal.

For the operator: the productive next move is to enable the same channel for the remaining Tier-1 tests (TEST-02, TEST-04, others), or to do a sweep audit. For my role: same posture — silence on stale CBP firings, respond to substantive new content like this proposal.

## Files

- `Research/Session632_500Mpc_Derivation_Audit.md` (this document)
- Updates to SESSION_FOCUS.md
