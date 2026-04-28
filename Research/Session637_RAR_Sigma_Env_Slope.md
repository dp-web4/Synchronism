# Session 637: RAR σ_int(ρ_env) — Derived Slope Is Unobservably Small

**Date**: 2026-04-28
**Type**: Site-Archive-Audit (7th instance) + first DERIVATION attempt
**Trigger**: visitor-flagged proposal `rar_sigma_int_environment_slope_derivation.md`
**Grade**: A- (specific, falsifies the highest-leverage claim)

---

## Setup

A leading-edge reviewer visiting the site (Pass 4) flagged environment-dependent RAR intrinsic scatter as **"the framework's one candidate for a genuinely novel, non-reparametrization prediction."** The proposal noted the site claim has no numerical slope, and asked whether one can be derived.

This is the first audit in the sub-arc that asks for a *derivation*, not a falsification check. Prior six audits (S631–S636) found existing claims wrong. S637 asks: if we take γ=2 as given (Session #64's 6D phase-space derivation) and propagate it, what slope comes out?

## Framework's Tools

From archive:
- `C(ρ) = tanh(γ · log(ρ/ρ_crit + 1))`, with γ = 2 from 3+3−4 = 2 (Session #64)
- `V²_total = V²_baryon / C(ρ)` at radius r (Session #65)
- `ρ_crit = A·V_flat²`, A = 0.029 (km/s)⁻² (Session #66)
- Session #64 explicitly: "ρ(r) = local density at radius r"

## The Mapping

The framework defines ρ in C(ρ) as **local density at the galaxy radius**. Environmental density enters only as an additive contribution:

```
ρ_total(r) = ρ_galactic(r) + ρ_env
```

Propagating environmental variance through C(ρ) to σ_int(RAR):

```
σ_int contribution from ρ_env  ≈  |d(log C)/d(log ρ)|_<ρ>  ×  σ_(log ρ_env)
```

This is the only mathematical channel for ρ_env to enter the framework's prediction. There is no other derivation path that gives ρ_env a stronger role without adding new physics.

## The Numbers

Typical scales:
- ρ_galactic_outer (V_flat regime, R~20 kpc) ≈ 4 × 10⁻²⁵ g/cm³
- ρ_env(cluster outskirts) ≈ 6 × 10⁻²⁸ g/cm³
- ρ_env(void) ≈ 3 × 10⁻³¹ g/cm³

Ratios:
- ρ_env(cluster) / ρ_galactic_outer ≈ 1.5 × 10⁻³
- ρ_env(void) / ρ_galactic_outer ≈ 7.5 × 10⁻⁷

Even in the densest environment, ρ_env contributes ~0.15% to local density at galaxy outskirts. ρ_total is everywhere dominated by ρ_galactic.

Computing the slope at ρ_total ≈ ρ_galactic_outer (taking ρ_crit = ρ_galactic_outer as the framework's natural reference):

```
d(log C)/d(log ρ)  ≈  0.25  (roughly constant across environment bins)
Δ(log ρ_total) from void → cluster  ≈  6.5 × 10⁻⁴ dex
Predicted Δσ_int(cluster − void)    ≈  1.6 × 10⁻⁴ dex
```

## Comparison to Baseline and Floor

| Quantity | Value |
|----------|-------|
| Predicted Δσ_int (cluster − void) | **0.00016 dex** |
| RAR σ_int baseline (Lelli+2017) | 0.13 dex |
| SPARC measurement floor per bin | ~0.02 dex |
| Δσ_int as fraction of baseline | 0.1% |
| Δσ_int as fraction of SPARC floor | 0.8% |

The predicted environmental signature is **~120× below current detectability**. It cannot distinguish Synchronism from MOND (constant σ_int) or from ΛCDM-hydrodynamics with current samples.

## Verdict on the Proposal

The visitor flagged σ_int(ρ_env) as the framework's one candidate for a genuinely novel prediction. After derivation:

- **The slope is computable**, not gestural — γ=2 fixes it.
- **The amplitude is microscopic** — three orders of magnitude below detectable.
- **The sign is correct** in the sense the proposal hoped: σ_int marginally smaller in cluster than void. But the magnitude is below the noise.

**Status**: not refuted, but **structurally undetectable from C(ρ) with the framework's own stated mechanism**. The "novel prediction" is unfalsifiable in practice — not because the framework is unfalsifiable in principle, but because the magnitude its math produces is ~120× too small to test.

For the prediction to be detectable, the framework would need an additional mechanism that amplifies environmental sensitivity by factor ≳ 100 — i.e., ρ_env entering C(ρ) through a pathway other than additive contribution to local density. No such pathway is in the archive.

## Why This Sharpens the Sub-Arc

S631–S636 found existing public claims wrong (BTFR, α², 500Mpc, coherence saturation, count, badge, category). S637 took the proposal's *highest-leverage candidate* — the one the visitor flagged as potentially saving the framework as something distinct — and showed it doesn't survive the derivation. The framework's quantitative content reduces to MOND in the testable regime.

This is qualitatively different from the prior 6 audits: those found errors. S637 found a derivation that *succeeds*, but produces an effect too small to test. Both outcomes — error and undetectability — leave the same conclusion: no novel testable prediction emerges.

The audit-channel taxonomy now extends:
1. Quantitative refutation (S631)
2. Dimensional inconsistency (S632)
3. Structural overclaim (S633)
4. Count discrepancy (S634)
5. Domain-level badge overclaim (S635)
6. Category error (S636)
7. **Derivation succeeds but predicts undetectable signal (S637, NEW)**

The Publisher's framing of S636 as "the visitor channel teaches itself what to look for" continues — the new mode is "even if the derivation works, what does it actually predict numerically?"

## Files

- `Research/Session637_RAR_Sigma_Env_Slope.md` (this document)
- `simulations/session637_rar_sigma_env_slope.py` (calculation)

## Recommended Site Action

The site claim "Synchronism predicts σ_int varies with local coherence density ρ_env" should be updated to specify:
- The predicted slope is dσ_int/d(log ρ_env) ≈ 0.25 dex per dex of effective local density.
- The predicted environmental difference (cluster outskirts vs void) is **~10⁻⁴ dex**, well below SPARC measurement floor.
- This is therefore not an experimentally distinguishable prediction with current data.

This is consistent with the framework's overall reduction to MOND in the testable regime.

## So What?

The visitor channel asked for the framework's strongest candidate to be derived. The derivation went through cleanly using γ=2 (a value the framework actually does derive, not fit). The result is that the quantitative prediction is ~10⁻⁴ dex — three orders of magnitude below detection. The framework's most-defended novel prediction is, in the testable regime, indistinguishable from MOND.

This isn't a death blow to the framework (the prediction could be true and just undetectable), but it removes the last claimed candidate for a discriminating test. Combined with S617–636, all routes to "Synchronism makes a testable prediction MOND/ΛCDM doesn't" are now closed.

The honest framing for the site: Synchronism reproduces MOND in the testable regime, predicts deviations below detectability, and contributes no current experimental discriminator. This is a defensible scientific claim — not a refutation of the framework, but a clear bounding of what it can be tested for.
