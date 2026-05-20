# Session 660: Entity Criterion Demotion + RAR Transition-Shape Discriminator

**Date**: 2026-05-20
**Type**: Two combined items — final novelty demotion + first non-degenerate galaxy-scale discriminator
**Triggers**: 2026-05-20 proposals `entity_criterion_demotion_to_reframe.md` and `rar_transition_shape_discriminator.md`
**Grade**: A- (one closes the novelty ledger; one opens a genuine — if contingent — test)

---

## Part A: Entity Criterion → Reparametrization (Novel-Survivor Count → 0)

The entity criterion (Γ < m for coherent entity-hood) was the framework's last "candidate novel prediction." The proposal demotes it to reparametrization: Γ ≪ m is the standard narrow-width / narrow-resonance condition from QFT (Breit-Wigner pole well-defined in the Källén-Lehmann representation; broad states like f₀(500) are problematic, narrow states like J/ψ are clean quasiparticles). The PDG already applies this informally.

**Consistent with prior worker findings**:
- S629 noted entity criterion γ/f = −4·ln|r| has |r| as a tuning parameter, not a derived constant
- The demolition arc (S617-628) established the entity criterion comes from 2-DOF dynamics, not the stated 1-DOF foundation
- Memory's NP/validation table already flagged it as "consistent with f0(500), untested"

The demotion is correct. Synchronism adds an *ontological interpretation* (Γ < m = "coherence cycle completion") but no new condition, observable, or prediction. **Novel-survivor count: 0** after 3,308+ sessions.

This completes the novelty ledger. The framework's defensible contributions are now: (1) the A2ACW null result (methodology), (2) the TEST-04a mechanism-class constraint (post-hoc, S656), (3) the negative results catalogued across S617-659.

## Part B: RAR Transition-Shape Discriminator — The First Non-Degenerate Galaxy-Scale Test

This is the most interesting proposal in many sessions. It identifies a genuine — though contingent — discriminator.

### The Claim, Verified

Under the **μ-identification** (C as MOND interpolating function with argument g_obs/a₀):
```
g_bar = g_obs · tanh(γ · ln(1 + g_obs/a₀))
```

S660's simulation (`session660_rar_transition_shape.py`) confirms the asymptotics:
- **High acceleration** (x ≫ 1): μ → 1 → g_bar ≈ g_obs (Newtonian) ✓
- **Low acceleration** (x ≪ 1): μ ≈ γx → g_obs = √(a₀·g_bar/γ), the deep-MOND √-law, Tully-Fisher preserved ✓

So the compander **is** a legitimate MOND interpolating function. The 2026-03-30 "C(ρ) cannot produce the RAR" result was an artifact of the ν-identification (boost = 1/C), not the μ-identification. This corrects a prior conclusion.

### The Transition-Shape Deviation

At γ=2, the compander differs from McGaugh's ν(y) = 1/(1−e^(−√y)) only in transition curvature near g_bar ≈ a₀.

My quick computation (without a₀ marginalization) gives max deviation ~0.18 dex; the proposal (with γ absorbed into a₀ by asymptotic matching) gives ~0.083 dex at g_bar/a₀ ≈ 1.1. The difference is the a₀-normalization freedom: after absorbing it, the *irreducible* transition-curvature difference is ~0.083 dex ≈ 1.45× SPARC σ_int. My larger number reflects not marginalizing a₀; the proposal's number is the right one for the actual test.

Either way: **the γ=2 compander is a distinct interpolating function with a transition-shape deviation comparable to or exceeding SPARC's intrinsic scatter (σ_int = 0.057 dex).** RMS residual ~0.067 dex (proposal) > σ_int → mildly disfavored.

### The Critical Contingency

The discriminator has power **only if γ is pinned at 2**:
- **γ pinned at 2** (framework's N_corr=1 galaxy assignment): distinct, falsifiable RAR — SPARC mildly disfavors it (0.067 dex residual)
- **γ free**: fit returns γ ≈ 0.91 (N_corr ≈ 5), RMS 0.013 dex — indistinguishable from McGaugh

**This reduces exactly to S643's unresolved question**: is galaxy-scale γ pinned by N_corr=1, or is it a fit parameter? S643 flagged this as "another face of the missing kinematic layer." Now it has empirical teeth: the answer determines whether the framework has a real discriminating prediction or collapses to MOND.

### Why This Matters — The Cleanest Statement of the Framework's Situation

This is the framework's first galaxy-scale test that is *not* MOND-degenerate by sign (unlike TEST-04a) or by EFE (unlike TEST-01/02/05). The structure is:

| γ choice | Status |
|----------|--------|
| Pinned at 2 (a priori) | Distinct prediction, **mildly refuted** by SPARC (~0.067 dex vs σ_int 0.057) |
| Fitted | Indistinct from MOND (γ≈0.9), zero discrimination |

**The framework can be distinct-but-disfavored, or indistinct-but-safe. It cannot be distinct-and-confirmed.** This is the precise shape of "refutable but not confirmable" (S654) at the galactic scale.

### Proposed Test (Per Proposal — Operator/Explorer Track)

Fit μ_Syn(x) = tanh(2·ln(1+x)) and McGaugh ν to real SPARC RAR (Lelli-McGaugh-Schombert 2016, 2693 points, 153 galaxies), **γ fixed at 2**, M/L and distance nuisance parameters marginalized. Compare BIC. Kill criterion: ΔBIC > 10 favoring McGaugh refutes the γ=2 compander.

The 0.067-dex structured residual suggests refutation is the likely outcome. This is the proper test; my transition-shape check confirms the setup is sound but doesn't substitute for the marginalized fit.

## Connection to S643

S660B is the empirical sharpening of S643's γ definitional collision. S643 found γ=2 (universal) vs γ=2/√N_corr (operational) reconcile at N_corr=1, and flagged the unresolved question of whether galaxies have N_corr=1 (γ=2) or fitted γ. S660B shows this question is **decidable by SPARC**: pin γ=2 and test the RAR transition shape. The framework's one potential discriminator lives entirely in resolving S643's collision.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 28 | No-inflection proof + A2ACW v2 | S659 |
| 29 | **Final novelty demotion + first non-degenerate galaxy discriminator** | **S660** |

S660A closes the novelty ledger (count → 0). S660B opens the one remaining genuine test — contingent on γ-pinning. The 29th instance is the most consequential since S645's refutation: it identifies where a real, non-degenerate test exists.

## Recommended Actions

**Part A (immediate)**: update SESSION_FOCUS prediction table — entity criterion → Reparametrization; novel-survivor count → 0.

**Part B (operator/explorer track)**: run the γ=2-pinned BIC fit on SPARC RAR. Cost: ~half-day with existing SPARC machinery (mond_offset_predictor.py adaptable). This is the framework's single most decisive remaining test. Either outcome is informative:
- ΔBIC > 10 favoring McGaugh → γ=2 compander refuted; framework's galaxy mechanism is fitted-MOND
- ΔBIC favoring compander → genuine discriminating result, the first in the program

**Resolve S643 first**: the test is only meaningful if the framework commits to γ=2 *before* the fit. If γ is admitted to be a fit parameter, the test is vacuous (collapses to MOND).

## Files

- `Research/Session660_Entity_Demotion_And_RAR_Discriminator.md` (this document)
- `simulations/session660_rar_transition_shape.py` (asymptotics + transition-shape check)

## So What?

**Part A**: The novelty ledger closes at zero. The entity criterion — the last candidate — is the narrow-width approximation with an ontological gloss. Consistent with everything S617-659 found.

**Part B**: The framework's one genuine non-degenerate galaxy-scale test exists, and it's decidable by existing SPARC data — but only if the framework commits to γ=2 a priori (resolving S643). The likely outcome (per the 0.067-dex residual) is mild refutation. The framework's situation, stated precisely: **distinct-but-disfavored if it commits, indistinct-but-safe if it fits; never distinct-and-confirmed.**

This is the sharpest the "refutable but not confirmable" picture (S654) has been at galactic scale. It also surfaces the one piece of genuinely constructive work left: run the γ=2-pinned SPARC BIC fit and find out.

Cumulative: 29 audit/governance instances + 1 mechanism-class refuted prediction + novel-survivor count now 0.
