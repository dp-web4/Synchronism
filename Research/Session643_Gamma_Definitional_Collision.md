# Session 643: γ Definitional Collision — Two Forms, Inverted Regime Labels

**Date**: 2026-05-06
**Type**: Site-Archive-Audit (13th instance, post-arc-closure)
**Trigger**: 2026-05-04 proposal `gamma_definitional_collision_regime_label_inversion.md`
**Grade**: B+ (resolves the collision; flags labeling error)

---

## Setup

The proposal flagged that γ appears in two incompatible forms across the site:
- **Form A**: γ = 2.0 (universal constant from 6D phase space)
- **Form B**: γ = 2/√N_corr (operational formula)

These reconcile at N_corr = 1, but the reconciliation is implicit, not stated. Downstream consequence: the γ-calculator labels BEC, BCS, neutron stars (γ ≪ 1) as "Classical" — inverting standard physics, where these are the canonical examples of macroscopic *quantum* coherence.

S643 traces both Forms in the archive and resolves which is load-bearing.

## Archive Findings

**Form B (γ = 2/√N_corr) is the load-bearing operational formula.** Established in:
- **Session #25 (Chemistry)**: γ = 2/√N_corr derived from fluctuation statistics. σ_corr/σ_uncorr = √N_corr.
- **Session #26 (Chemistry)**: N_corr is operationally measured per-system via fluctuation analysis: `N_corr = (σ_measured / σ_uncorrelated)²`. This is a per-system measurement, not a first-principles derivation.

**Form A (γ = 2.0)** corresponds to N_corr = 1 in Form B. It is a *normalization condition at the single-particle reference*, not a separate derivation. The "6D phase-space derivation" referenced on `/coherence-function` (γ = 3 + 3 − 4 = 2) yields γ = 2 only at N_corr = 1; for any N_corr > 1, γ < 2.

The collision: the site presents Form A as if it were a universal derivation, while the operational machinery uses Form B. They are not separate results — Form A is Form B at one point.

## Regime Labels Are Inverted

Using Form B with site presets:

| System | N_corr | γ = 2/√N_corr | Site Regime | Standard Physics Label |
|--------|--------|---------------|-------------|------------------------|
| Ideal gas | 1 | 2.000 | Quantum (γ > 1.4) | Classical |
| BCS Al | 10⁴ | 0.020 | Classical (γ < 0.6) | **Quantum** (macroscopic) |
| BEC | 10⁶ | 0.002 | Classical | **Quantum** (macroscopic) |
| Neutron star | 10¹⁰ | 0.0002 | Classical | **Quantum** (degenerate) |

The site's labels are inverted relative to standard physics for every textbook quantum-collective system.

## What γ Actually Measures

The math implements: **high γ ⇔ low N_corr ⇔ few correlated DOF ⇔ single-particle-like behavior**. Low γ ⇔ high N_corr ⇔ many correlated DOF ⇔ collective behavior.

This is internally consistent if read as **single-particle vs collective**, not as **quantum vs classical**. BEC is highly collective (high N_corr), so it correctly sits in the "low γ" basin — but that basin is *collective*, not *classical*. Standard physics calls BEC quantum because of its macroscopic phase coherence; the framework's γ doesn't measure macroscopic phase coherence at all, just correlation count.

## Resolution: Case 1 + Label Renaming (Proposal's Preferred Path)

The proposal listed three cases. Archive evidence resolves to Case 1 with relabeling:

- **Case 1**: γ = 2/√N_corr with N_corr = dynamical correlation count. Rename labels: "Single-particle / Uncorrelated" and "Collective / Strongly correlated."
- **Case 2** (N_corr = quantum overlap): contradicts Session #26's fluctuation derivation. Reject.
- **Case 3** (γ = 2 fixed for galaxies, separate calibration for CM): would require explicit scope statement and treats galactic baryons as N_corr = 1 (uncorrelated). Possible but ad hoc.

Archive supports Case 1. The labeling error is on the public site, not in the underlying math.

## Connection to S640 (Dual-C)

Same mechanism, different scope:
- **S640**: C used for two functional forms (C(ρ) = tanh(γ·ln(ρ/ρ_crit+1)) vs C = f(γ,D,S)).
- **S643**: γ used for two formulas (γ = 2.0 vs γ = 2/√N_corr).

In both, the framework presents notation as if it implies a derivation that doesn't exist. The collisions are notational, not physical. Both can be fixed by explicit relabeling and scope statements.

## Unresolved: First-Principles N_corr

The proposal's deeper question — *"is there a first-principles derivation of N_corr for any specific system?"* — remains open. Session #26 gives operational measurement methods (fluctuation analysis, correlation length, susceptibility). It does not give a derivation that, e.g., predicts N_corr for a galaxy from gravitational physics, or for a BCS superconductor from BCS theory.

For galactic dynamics, the framework requires γ = 2 → N_corr = 1, treating the galaxy as a single-particle system. This is an *odd* assumption for 10¹¹ stars and is currently unmotivated by archive content. It is closer to a calibration choice than a derivation.

This connects to the kinematic-layer pattern (S641, S642): N_corr is one face of the missing kinematic substrate — the framework specifies a formula that uses N_corr but never derives N_corr from underlying physics.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 10 | Symbol overloading at foundational level (C) | S640 |
| 11 | Cross-gap synthesis (Lorentz as 4th face) | S641 |
| 12 | Field-or-parameterization disambiguation (GW170817) | S642 |
| 13 | **Regime-label inversion under operational formula (NEW)** | **S643** |

S640 was *symbol* overloading. S643 is *regime-label* inversion — a downstream presentation error from the same kind of formula collision. Same audit-channel mechanism (shared label, divergent semantics), different surface.

## Recommended Site Action

**Short-term**:
- Rename "Quantum regime" → "Single-particle / Uncorrelated regime" in γ-calculator and phase-boundary-visualizer.
- Rename "Classical regime" → "Collective / Strongly correlated regime."
- Add a caveat: "Note: These regime labels are about correlation count, not quantum-vs-classical in the standard sense. BEC, BCS, and neutron stars sit in the 'Collective' basin because they have many correlated DOF, not because they are classical."

**Medium-term**:
- State the Form A ↔ Form B reconciliation explicitly. Form A is Form B at N_corr = 1, not a separate derivation.

**Long-term**:
- Address whether γ = 2 (N_corr = 1) for galaxies is a derivation or a calibration. If calibration, name it as such. If derivation, write it.

## Files

- `Research/Session643_Gamma_Definitional_Collision.md` (this document)
- No simulation needed — definitional clarification

## So What?

The γ formula is internally consistent (Form A and B reconcile at N_corr = 1), but the public site's presentation conflates them and inverts standard regime labels. This is the same pattern as S640's dual-C audit: notational unity that masks the absence of a derivation. Renaming the labels and stating the reconciliation explicitly fixes the surface problem. The deeper question — why N_corr = 1 for galaxies, why N_corr is calibrated rather than derived for any specific system — remains as another face of the missing kinematic layer (S641, S642).

Two proposals remain pending: ρ_crit calibration vs prediction (2026-05-05), session107 DESI DR1 (2026-05-05). Same pattern expected.
