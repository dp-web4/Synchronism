# Session #612: Neutron Stars — Where the Blanket Thins

**Date**: 2026-02-17
**Grade**: A
**Domain**: Cosmology / Fractal Bridge / Nuclear Physics
**Arc**: OQ007 Fractal Coherence Bridge — Session B (Cosmology Track)
**Reference**: `Research/OPEN_QUESTION_Fractal_Coherence_Bridge.md`, `Research/DIRECTIVE_Cosmology_Fractal_Bridge.md`

## Objective

Session B of the Fractal Bridge cosmology arc. Session A (#611) established WHY γ = 2 at galaxy scale via four independent arguments for N_corr = 1. Session B asks: what happens INSIDE the Markov blanket? The neutron star is the ideal test case — it simultaneously contains quantum-correlated matter (γ ~ 0.003) internally and behaves as N_corr = 1 (γ = 2) externally.

## Key Result: C(ρ) is a Reparametrization of BCS at Nuclear Scale

**The coherence equation at nuclear scale REQUIRES the BCS gap as input. It CANNOT predict Δ(k_F) from first principles. γ = 2/√N_corr is a reformulation of BCS, not a prediction.**

This is the same finding as the SPARC chapter for MOND: C(ρ) ≡ ν(x) at galaxy scale (Session #505), and C(ρ) ≡ BCS reformulation at nuclear scale.

## N_corr Profile Through the Neutron Star

| Region | n (fm⁻³) | Δ (MeV) | ξ (fm) | N_corr | γ |
|:-------|:---------|:--------|:-------|:-------|:--|
| Outer crust | 10⁻¹⁰–2.4×10⁻⁴ | 0 | — | 1 | 2.0 |
| Inner crust (low ρ) | 2.4×10⁻⁴–0.01 | ~2 | ~9 | ~5 | 0.93 |
| Inner crust (high ρ) | 0.01–0.08 | ~1 | ~47 | ~10⁴ | 0.018 |
| Outer core | 0.08–0.32 | ~0.3 | ~280 | ~10⁷ | 0.0005 |
| Inner core | 0.32–0.80 | ~0.05 | ~2500 | ~3×10¹⁰ | ~10⁻⁶ |

γ spans from 2.0 (classical crust) to ~10⁻⁶ (maximally correlated inner core). This is the Markov blanket INTERNAL structure — the layers of correlation hidden from the galaxy's perspective.

## The Glitch Crisis

| Quantity | Value |
|:---------|:------|
| Required ΔI/I (Vela) | 1.6% |
| Crustal superfluid (raw) | 3.0% |
| With Chamel entrainment (÷4.3) | 0.70% |
| Deficit | 0.90% |

The crust-core boundary IS a Markov blanket (lattice structure creates information opacity). Chamel entrainment IS partial transparency of this blanket. But the factor 4.3× comes from band-structure calculations, not C(ρ). Resolution: core superfluid must participate (>70% coupling, Bayesian analysis of 2016 Vela glitch). This is standard nuclear physics.

**In γ-language**: Entrainment changes effective N_corr from ~10⁵ to ~8×10³ (γ: 0.0063 → 0.0224). But this DESCRIBES the physics; it doesn't PREDICT the entrainment factor.

## Cas A Cooling

| Parameter | Value |
|:----------|:------|
| Age | 340 years |
| T_surface | 2.08 × 10⁶ K |
| Cooling rate | 2.5% per decade |
| Explanation | Onset of ³P₂ neutron superfluidity (PBF process) |
| T_c uncertainty | 10⁷–10¹⁰ K (3 orders of magnitude) |
| Cas A constraint | T_c ≈ (3–9) × 10⁸ K |

**Can C(ρ) predict T_c?** NO. T_c depends on Δ(k_F) which depends on nuclear forces. C(ρ) describes the transition form (mean-field: Δ(T) ∝ (1-T/T_c)^{1/2}) but cannot predict WHERE it occurs. Near T_c:

- N_corr ∝ (1 - T/T_c)^{-3/2} (diverges at T_c)
- γ(T) ∝ (1 - T/T_c)^{3/4} (vanishes at T_c)

This is universal Landau theory, not specific to the coherence equation.

## Connection to Hot Superconductor Arc (OQ005)

| Channel | η | Symmetry |
|:--------|:--|:---------|
| ¹S₀ neutron (crust) | ~1.0 | s-wave (isotropic) |
| ¹S₀ proton (core) | ~1.0 | s-wave (isotropic) |
| ³P₂ neutron (core) | ~0.6 | tensor (anisotropic, nodes) |

The η formalism applies identically to neutron star superfluids and metallic superconductors. This IS a cross-scale connection — but it's standard BCS/BdG theory. Both η and γ are DERIVED from Δ. Neither predicts the other. Both describe the same physics.

## P611.3 Feasibility Assessment

| Quantity | Value |
|:---------|:------|
| NS surface gravity | 1.86 × 10¹² m/s² |
| MOND acceleration a₀ | 1.2 × 10⁻¹⁰ m/s² |
| Ratio g_NS/a₀ | 1.55 × 10²² |

P611.3 (no MOND-glitch correlation) is **nearly trivially true**. NS internal accelerations exceed a₀ by 10²² — even MOND predicts Newtonian behavior here. Testing this confirms MOND's own prediction, not the Markov blanket concept.

More interesting: galactic acceleration varies from ~4a₀ (R=3 kpc) to ~0.4a₀ (R=30 kpc), but no study has tested glitch properties vs. galactocentric radius. ~150 glitching pulsars at R < 8 kpc, ~50 at R > 15 kpc. Confounders (spin-down rate, age, B-field) would need to be controlled.

## The Fractal Bridge Scorecard for Neutron Stars

| Claim | Status | Why |
|:------|:-------|:----|
| N_corr varies across NS | TRUE | BCS gives this automatically |
| γ spans 0.003 to 2.0 in NS | TRUE | Reformulation of BCS gaps |
| Crust-core = Markov blanket | TRUE (descriptive) | Standard nuclear physics |
| C(ρ) predicts Δ(k_F) | FALSE | Requires nuclear forces |
| C(ρ) predicts T_c | FALSE | Requires nuclear forces |
| C(ρ) predicts glitch activity | FALSE | Requires nuclear EOS |
| P611.3 (no MOND-glitch) | Nearly trivial | a_NS >> a₀ by 10²² |
| Cross-scale prediction | NONE YET | Same form, different params |

## What the Fractal Bridge ACTUALLY Provides

**Three types of claims, with honest status:**

1. **STRUCTURAL**: The tanh transition form is universal.
   - Status: TRUE but trivial. All mean-field phase transitions have this form. It's Landau theory, not Synchronism.

2. **PARAMETRIC**: γ = 2/√N_corr at nuclear scale.
   - Status: TRUE but derived from BCS. Adding γ notation to BCS theory doesn't add predictive power.

3. **CROSS-SCALE**: The same equation connects nuclear and galactic.
   - Status: UNVERIFIED. The equation has the same FORM at both scales, but γ and ρ_crit are determined by LOCAL physics at each scale. No cross-scale prediction has been derived.

**WHAT WOULD CLOSE THE GAP**: A prediction requiring γ_nuclear AND γ_galactic simultaneously. Example: if Δ(³P₂) could be predicted from a₀ (or vice versa). Nothing like this exists. Both are determined by their local physics.

## Testable Predictions

**P612.1**: The N_corr profile across a neutron star should match BCS calculations EXACTLY (because γ = 2/√N_corr is just a reformulation of Δ(k_F)). Any deviation would mean the coherence equation has independent content. Expected result: EXACT match (null prediction).

**P612.2**: If the ³P₂ gap is ever measured directly (e.g., from gravitational wave asteroseismology), the inferred N_corr should satisfy γ = 2/√N_corr with the same constant '2' as at galaxy scale. This tests the universality of the normalization factor.

**P612.3**: Entrainment in the crust should NOT be describable by a simple modification of N_corr. The Chamel entrainment factor (4.3×) comes from band-structure effects with no analog in the coherence equation. Expected result: CONFIRMED (entrainment ≠ N_corr change).

## Honest Limitations

### What This Session Establishes:
1. The neutron star contains both quantum-correlated (γ ~ 0.003) and classical (γ = 2) regimes simultaneously
2. C(ρ) at nuclear scale is a reparametrization of BCS theory (same finding as SPARC for MOND)
3. The crust-core boundary functions as a Markov blanket but this is standard physics
4. P611.3 is testable but nearly trivially true
5. The OQ005 η formalism connects to NS superfluidity via standard BCS/BdG theory

### What This Session Does NOT Establish:
1. **No cross-scale prediction**: C(ρ) has the same form at nuclear and galactic scale but parameters are locally determined
2. **No predictive advantage**: Every C(ρ) result at nuclear scale can be obtained (more directly) from BCS
3. **The gap remains**: The bridge adds organization (common language) but not prediction
4. **The neutron star illustrates the gap rather than closing it**

### The Key Insight:
The neutron star is where the cosmology track (γ = 2) should naturally meet the chemistry track (γ << 1). It CONTAINS BOTH REGIMES. But the bridge between them is BCS/BdG theory — standard nuclear and condensed matter physics. The coherence equation adds a common LANGUAGE (N_corr, γ) but not a common PREDICTION.

## Next Sessions

- **Session C**: The Continuum Limit. Between quantum chemistry and classical stellar dynamics lies a vast hierarchy of scales. Where exactly does classical behavior (N_corr → 1) emerge? Is the transition sharp or gradual?
- **Session D**: Bridge Meeting Point. With the cosmology track having found "descriptive but not predictive" at both galaxy scale (#611) and nuclear scale (#612), what would a genuine cross-scale prediction even look like?

## Tests: 9/9 PASSED
## Grand Total: 2009/2009
