# Proposal: A-from-Jeans Closure — The 5% Number Belongs to a Different Law

**Submitted:** 2026-06-07
**Source:** Synchronism site Explorer track (reconstructed from primary sources: Session 65,
Session 66 markdown + `session66_A_gap_investigation.py`, Session 631, Session 644)
**Priority:** HIGH — closes the framework's last surviving first-principles claim
**Supersedes the open framing in:** `a_from_jeans_r0_universality_flaw.md` (2026-06-07)

## Result

The "A = 4π/(β_J²·G·R₀²) ≈ 0.0294 vs empirical 0.028, 5% agreement" claim does not survive a
re-run of the original computation. Three findings, all numerically reproduced
(`synchronism-site/explorer/scripts/a_from_jeans_chain_of_custody.py`):

1. **The 5% computation derives the wrong scaling law.** The only calculation that yields 0.0294
   (`session66_A_gap_investigation.py`) uses α = 4.5 and R₀ = 0.07 kpc/(km/s)^0.75 — the
   *fitted slope* of the size-velocity relation R_half = R₀·V^0.75 — and produces
   **ρ_crit ∝ V^0.5** (Session 65 confirms B = 0.5 explicitly). The framework's operative law is
   **ρ_crit ∝ V²**. The derivation underpins a law the framework does not use.

2. **The stated formula does not reproduce its own number.** A = 4π/(β_J²·G·R₀²) with β_J = 1,
   R₀ = 8 kpc (the inputs in the S66 markdown, S644, and the public site) gives
   **A ≈ 4.6×10⁻⁵**, ~600× from 0.028. The S66 markdown's own arithmetic (4.57×10⁻⁵) reaches
   0.0294 only via an **unexplained 644× "unit conversion."** The 5% was imported from the
   V^0.5 script, which has different inputs and a different exponent.

3. **Three fitted/chosen inputs, not one.** α = 4.5 (fitted; varies 1.3–4.5 by galaxy type per
   S53), R₀ = 0.07 (fitted size-velocity slope), and 4π (grid-searched post-hoc from a list of
   constants near 12). "R₀ = 8 kpc, the Sun's galactocentric radius" is a *mislabel* introduced
   in the S66 markdown — it is not the quantity used to obtain 5%.

## Why this matters for the archive

This is the static analogue of the 2026-05-25 DESI epistemic-regression event. There a confident
correction overwrote a verified result; here a headline number (0.0294) detached from its
generating computation and propagated through ~600 sessions and onto the public site. S631 and
S644 both *re-read* the derivation and restated its inputs, but **neither re-ran the arithmetic**
— if they had, β_J=1, R₀=8 kpc would have returned 4.6×10⁻⁵, not 0.0294. The methodology lesson:
re-reading a derivation is not auditing it. Auditing means re-executing it.

## Disposition

- **A-from-Jeans → Reparametrization / Audited-Negative.** It is a calibration of A=0.028 (itself
  the V^0.5-law empirical coefficient) reused with mismatched units under a V² law it does not
  derive.
- **Framework status: zero first-principles predictions with an independent derivation.** a₀
  (Milgrom coincidence), Σ₀ (Freeman re-expression), R₀ (follows from a₀), Γ = γ²(1−c)
  (Palma-Suominen-Ekert 1996), RAR γ=2 (refuted ΔBIC=+184), A-from-Jeans (this closure). The
  physics audit reaches a clean terminal state.
- **Remaining live route (Session 644 Path C):** independently-measured σ → β_J → test whether
  ρ_crit = V²/(G β_J² R_half²) reproduces A across SPARC without circularity. Note it would
  predict the V^0.5 law, not V². Worth running on SPARC σ data ($0).

## Recommended site updates

See `synchronism-site/explorer/findings/a-from-jeans-chain-of-custody-failure.md` →
"Action: Maintainer". In brief: replace the /parameter-derivations open-question box with the
closure; downgrade the badge; update honest-assessment to "zero independently-derived
first-principles predictions"; reconcile the V^0.5 / V² exponent inconsistency.
