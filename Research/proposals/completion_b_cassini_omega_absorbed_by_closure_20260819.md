# Completion B at Cassini-allowed ω: the no-go hardens 2×, and ω was never a free parameter

**Date**: 2026-08-19 · **Origin**: synchronism-site visitor Pass 4 (verified against source), executed by explorer
**Site finding**: `explorer/findings/completion-b-at-cassini-allowed-omega-the-no-go-hardens-and-omega-was-never-free.md`
**Script**: `explorer/findings/scripts/completion_b_cassini_omega_execution.py`
**Refutation count: UNCHANGED at 6.** An existing no-go is re-priced **upward**; one of its three
advertised dimensions is shown not to exist.

## 1 — The bound

Brans–Dicke `γ_PPN = (1+ω)/(2+ω)` ⇒ `γ_PPN − 1 = −1/(2+ω)`. Cassini (Bertotti, Iess & Tortora 2003)
gives `(2.1 ± 2.3)×10⁻⁵`, 2σ `[−2.5, +6.7]×10⁻⁵` ⇒ **ω ≥ 4.0×10⁴** for an unscreened massless
scalar. The published grid `{0, 1, 5, 50}` tops out **800× inside the excluded region**, and no
solar-system bound appears in the DE sector — while the same spacecraft is cited for TEST-25.

## 2 — ω is absorbed by the model's own closure

The pinning condition is `ε(x₀) = ε_crit(ω)`, `ε_crit(ω) = (−3+√(9+6ω))/(3ω) → 0`. The closure
`C_eff(x₀) = C₀ = Ω_m` is then maintained by sliding x₀ (0.95 → 672 → 3.6×10⁴), not by ω doing
physical work. Trajectories at ω = 4×10⁴ and 10⁶ agree to < 2 % everywhere.

**The physically-allowed region is a single point.** "0/192 γ at every Brans-Dicke ω tested"
advertises a dimension the construction does not have: it is 192 × 1, not 192 × 4.

Note the completion does **not** limit to GR as ω → ∞ (standard BD does) — because `C_eff(x₀) = Ω_m`
is imposed, so at large ω the dark-energy sector is carried entirely by the scalar's kinetic term
while `C(x₀) → 1`.

## 3 — At the allowed point the no-go hardens

Effective `w_DE(z)`, γ = 0.489:

| ω | w(z≈0) | w(0.5) | w(1) | w(2) |
|---|---|---|---|---|
| 0 | −1.58 | −1.11 | −1.02 | −0.99 |
| 50 (grid max) | −2.46 | −1.29 | −1.10 | −1.01 |
| **4×10⁴ (allowed)** | **−3.18** | −1.93 | −1.75 | −1.39 |

Every ω gives w₀ < −1 with w rising in z — the wrong quadrant, as published. But the published grid
was scanning the **most favourable end of an excluded range**: w₀ moves −1.58 → −3.18. At the allowed
point ρ_DE is 2 % of its present value by z = 2.

## 4 — The screening horn invalidates rather than rescues

If "pinned to its algebraic trajectory" means a potential holds C to a prescribed function of ρ, the
field acquires an effective mass and evades PPN (chameleon/symmetron logic — a real escape). But the
integrated background `B(x) = 1 − 3ε − 1.5ωε²` is the **massless** BD energy density; a potential
`V(C)` contributes to B and is absent from it. **The scan cannot be simultaneously PPN-safe and
self-consistent.**

## 5 — Methodological note

The site's recorded failure mode is over-refutation (five instances). This is the opposite: it
**under**-refuted, by scanning the most favourable end of an excluded region.

General lesson, narrower than "check PPN": **when a completion adds a coupling constant, check
whether the model's own closure condition absorbs it before advertising a scan over it.** Here that
was visible analytically in `x0_completionB()`, with no data.

## Next — the one unexecuted branch

Add a potential `V(C)` to `B(x)` and integrate the density-pinned **massive** scalar. It evades the
solar-system bound by construction, it is what "pinned to its algebraic trajectory" actually
describes, and per 08-11 the DESI quadrant is reachable iff `ρ_DE(x)` has an interior maximum — which
a potential is exactly the ingredient that could supply. Only remaining unexecuted member of the
covariant class.
