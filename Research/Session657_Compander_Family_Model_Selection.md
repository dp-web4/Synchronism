# Session 657: Compander-Family Model Selection (AIC/BIC) — Endorsement with Caveats

**Date**: 2026-05-17
**Type**: Methodology endorsement (not executed computation)
**Trigger**: 2026-05-17 proposal `compander_family_model_selection_aic_bic.md`
**Grade**: B (endorses; flags scope as operator-level)

---

## Setup

The proposal asks the framework to run an AIC/BIC comparison across the compander family — tanh, Hill/Naka-Rushton, logistic, erf, μ-law, Gompertz — on SPARC + chemistry + Tc datasets. This is the natural follow-up to the S653 compander commitment, and the site's own `/why-synchronism` page already concedes the comparison is needed.

S657 endorses the proposal, notes prior partial results, and flags the scope.

## The Proposal's Logic Is Sound

- The site admits ("any S-curve with the same saturation properties would fit the same data equally well") that tanh has no privileged status as a compander.
- C(ρ) is not a self-consistency equation (S636/S638/S640/S652); tanh's special status in mean-field Landau theory does not transfer.
- AIC/BIC across the compander family is the standard model-selection diagnostic.

The proposal correctly identifies this as the central methodological gap left by the post-Kimi "Findings vs Framings" discipline: the framework concedes adequacy but hasn't *demonstrated* it via comparison.

## Prior Finding (Per Proposal)

The 2026-03-27 explorer session ran Hill vs tanh on the coupling-coherence dataset. Initial result Hill won ΔAIC=4; after a baseline fix, tanh won ΔAIC=17.6. So at least *one* dataset has been compared, with tanh winning decisively.

This is one data point. SPARC + chemistry + Tc require separate fits. The proposal's request for the full panel is genuine new computation, not just retrieval.

## Caveats for the Comparison

If pursued, the comparison must respect:

1. **AIC/BIC limitations on non-nested models**: tanh and Hill aren't nested (Hill contains tanh as n→∞ on log axes, per proposal, but the parameter penalty matters). Standard AIC/BIC handles this but Bayesian model selection with explicit priors would be more rigorous.

2. **Different datasets may favor different forms**: galactic dynamics, chemistry correlations, and superconductor T_c live in very different regimes (different relevant ρ ranges, different noise characteristics, different baselines). One winner across all three would be surprising. The likely outcome is "different compander forms win in different regimes," which is itself informative.

3. **Baseline carries**: the 2026-03-27 result swung from Hill-wins-by-4 to tanh-wins-by-17.6 when the baseline was fixed. Any new comparison must explicitly state the baseline and verify it.

4. **The structural finding doesn't change**: whether tanh, Hill, or logistic wins, C(ρ) remains a compander with no field equation (S652), no governing dynamics (S638), no self-consistency loop (S636). Winning a compander family model-selection contest is local within the compander class; it does not promote C(ρ) to a derived object.

## What Each Outcome Means

- **Tanh wins by AIC/BIC on all datasets**: Confirms current choice; framework gains a defensible "we ran the model selection" statement; structural picture unchanged.
- **Hill wins on some datasets**: Important — means C(ρ) should switch to Hill for those domains. The current parametrization is wrong for those regimes. Reframes the "framework" as a *family* of compander parametrizations selected per-domain.
- **No significant difference**: Confirms the site's honest admission. Tanh is one adequate representative of the compander class, not the unique correct answer. The "ONE EQUATION" framing weakens further.

Each outcome is informative; none changes the underlying conclusion that the framework's predictive content is whatever the compander class predicts.

## Why I Haven't Run It

Full SPARC fitting requires the predictor machinery (mond_offset_predictor.py from S588) and would need adaptation to each compander variant. Chemistry data requires the 1,703-phenomenon corpus access (with S647's caveats about N_corr method specification). T_c data requires the superconductor dataset. Each fit is straightforward; setting them up consistently is ~half a day's work and is operator-track research, not a worker-session quick check.

The proposal estimates 2-4 hours. That's realistic for a focused executor session with all the machinery pre-existing. It's not realistic to fold into one back-annotation cycle alongside other audits.

**Recommendation**: this is an operator/explorer-track task. Worker sessions can confirm the methodology, flag prior results, and note caveats — which is what S657 does.

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 25 | Failure-as-contribution reframing | S656 |
| 26 | **Compander-family model-selection scope endorsement** | **S657** |

S657 is a scope/methodology endorsement, not an audit finding. The 26th "instance" continues the pattern of governance-adjacent work — confirming that operator-track proposals have sound methodology and audit findings support them.

## Recommended Operator Action

1. **If pursuing the comparison** (per proposal): assign to explorer track with explicit baseline specification, Bayesian model selection (or AIC/BIC with proper non-nested handling), and per-domain reporting.
2. **If not pursuing**: update `/why-synchronism` to reference the prior 2026-03-27 result (tanh wins coupling-coherence by ΔAIC=17.6) and acknowledge the broader comparison is open.
3. **Either way**: maintain S653's Frame B commitment. The compander framing stands regardless of which specific compander wins model selection.

## Files

- `Research/Session657_Compander_Family_Model_Selection.md` (this document)

## So What?

The proposal's methodology is sound and addresses a real gap. Prior work on coupling-coherence shows tanh wins by ΔAIC=17.6 on that dataset; the broader SPARC + chemistry + T_c comparison is open and operator-track.

The structural picture is stable across whatever outcome: tanh, Hill, or any other compander winning model selection doesn't promote C(ρ) from forward-map to derived field equation. Compander family model selection is local within the class S636/S638/S640/S652 already established as the framework's actual structural location.

Cumulative: 26 audit/governance instances + 1 mechanism-class refuted prediction. The audit channel continues to grow slowly through governance-adjacent recommendations rather than new structural findings.
