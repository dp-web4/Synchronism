# Session 649: QM Kill Criterion Underspecified + ρ_crit Naming Asymmetry

**Date**: 2026-05-08
**Type**: Site-Archive-Audit (combining 17th and 18th instances)
**Triggers**: 2026-05-08 proposals `qm_kill_criterion_dd_gap.md` and `rho_crit_asymmetry_saturation_knee.md`
**Grade**: B+ (two relabeling/specification fixes; both confirmed against archive)

---

## Setup

Two proposals filed same day, both small and both pointing at established patterns. Combining:

1. **QM kill criterion** ("design noise environment where resync outperforms isolation") is unfalsifiable as written — dynamical decoupling (DD) already demonstrates pulse-sequence advantage over passive isolation in 1/f baths.
2. **ρ_crit at γ=2 gives C ≈ 0.88, not 0.5** — the "+1" regulator in `ln(ρ/ρ_crit + 1)` asymmetrizes the sigmoid; ρ_crit is a saturation knee, not a critical density.

Both connect to prior audit findings.

## Part A: QM Kill Criterion vs DD Literature

The site's Key Claim #1 advertises:

> "Design a noise environment where resync outperforms isolation, but standard decoherence theory predicts it doesn't. If isolation wins uniformly, the synchronization ontology adds nothing."

The proposal correctly notes: **DD literature (Viola-Knill-Lloyd 1999, CPMG, UDD, transmon experiments) already shows periodic pulse sequences ("resync") beat passive isolation in non-Markovian baths**. The kill criterion as worded is satisfied by every standard QIP textbook result.

To make the criterion falsifiable, it must specify:
- Specific physical system (transmon, NV center, trapped ion)
- Bath spectral density (Ohmic, 1/f, structured) — DD's advantage is bath-dependent
- Exact "resync" protocol (CPMG, UDD; pulse number N, interpulse spacing)
- Numerical T2(resync)/T2(isolation) ratio predicted by Synchronism vs by Bloch-Redfield+DD
- Kill threshold ratio at which the two diverge

None of this is in the archive. The criterion is unfalsifiable in its current form.

**Connection to S581 (Quantum Coherence Audit)**: S581 already established the quantum-coherence arc is reparametrization of standard QM, with 0 unique predictions across ~14 sessions. The QM kill criterion sits in the same pattern: vocabulary statement of phenomena DD already explains.

**Recommended action**: derive T2 from MRH dynamics for a transmon in 1/f bath under CPMG. Compare to standard Bloch-Redfield + UDD prediction. If equations are identical, label "DD reparametrization" on `/key-claims` and retire the kill criterion. This is a one-session research task; outcome (per S581's pattern) is likely reparametrization.

## Part B: ρ_crit Is a Saturation Knee, Not a Critical Density

The proposal's calculation:

```
C(ρ_crit, γ=2) = tanh(2 · ln(ρ_crit/ρ_crit + 1))
              = tanh(2 · ln 2)
              = tanh(1.3863)
              ≈ 0.8824
```

This is correct and was confirmed by S638's verification (which evaluated C(ρ_crit) at multiple γ values: γ=0.5→0.333, γ=1→0.600, γ=2→0.882). S638 already concluded "ρ_crit is field-zero offset, not critical density" — the framework's "critical" naming inherits phase-transition vocabulary the math doesn't support.

The proposal's contribution is the "+1" regulator analysis: it identifies *why* the asymmetry exists. The "+1" prevents `ln(0)` divergence as ρ→0 but has no physical motivation beyond numerical stability. It asymmetrizes the sigmoid such that ρ_crit is at C ≈ 0.88 rather than at the half-maximum.

The half-maximum would actually be at:
```
ρ_half = ρ_crit · (e^(1/(2γ)) − 1)
       ≈ 0.284 · ρ_crit  (for γ=2)
```

**Three resolution options** (from proposal):
- **A. Rename**: ρ_crit → ρ_scale or ρ_knee. Symbol search-and-replace; no numerical changes.
- **B. Recenter**: replace with `C(ρ) = (1 + tanh(γ·ln(ρ/ρ_crit)))/2`. ρ_crit becomes true half-maximum; all existing fits would need refitting.
- **C. Reframe**: keep the current form, but document explicitly that "ρ_crit" is a scale parameter, not a critical density.

**Recommendation**: Option C immediately (notation note), with Option B as a research-session topic to refit SPARC/ALFALFA data and see if the recentered form performs equivalently. Option A is acceptable but loses the "critical" connotation completely.

**Connection to S638**: S638 verified C(ρ) reduces to a Curie-paramagnet response. ρ_crit is the field-zero offset of that response. Calling it "critical density" is the same kind of phase-transition vocabulary that misled the framework into treating C(ρ) as Landau theory (S636) and into treating C ≈ 0.88 as a phase boundary (current parameter naming).

## Audit Taxonomy

| # | Type | Session |
|---|------|---------|
| 16 | Self-correction of prior session | S648 |
| 17 | **QM kill criterion underspecified vs known DD literature** | **S649A** |
| 18 | **Parameter-name asymmetry in regulated sigmoid** | **S649B** |

Both fit the established pattern: site overstates archive content. S649A is a new variant — *kill criterion specified at vocabulary level, falsified by existing literature regardless of framework*. S649B is a follow-up to S638 — the asymmetry was visible in S638's table but not named as a parameter-naming issue.

## Recommended Site Action

**For QM kill criterion**:
- Either derive T2 from MRH and demonstrate a difference from Bloch-Redfield+DD, OR
- Remove the kill criterion as written and replace with "DD reparametrization" label per S581's existing diagnosis

**For ρ_crit naming**:
- Add prominent note on `/critical-density` and `/coherence-function`: "ρ_crit is the parameter where C ≈ 0.88 (at γ=2), not C = 0.5. It is a saturation-knee scale parameter, not a critical density in the Landau sense. The half-maximum is at ρ ≈ 0.284 ρ_crit for γ=2."

## Connection to S646 Meta-Criterion

S646 recommended a meta-falsification criterion. Both findings here are relevant:
- The QM kill criterion falls into the meta-criterion question: a kill criterion that already-existing literature satisfies cannot count toward framework retraction; it should be flagged as unfalsifiable rather than treated as a live test.
- The ρ_crit naming is a presentation issue, not a kill-criterion issue — but it joins S639/S640/S643 as a naming-vs-math conflation that the meta-criterion process should normalize.

## Files

- `Research/Session649_QM_Kill_And_RhoCrit_Asymmetry.md` (this document)

## So What?

Two small audits with concrete site fixes:

- **QM kill criterion** is unfalsifiable as written because DD literature already satisfies it. Either specify it at QIP-experiment level, or label it "DD reparametrization" per S581's pattern.
- **ρ_crit at γ=2 gives C ≈ 0.88**, not 0.5. The "+1" regulator asymmetrizes the sigmoid. ρ_crit is a saturation knee, not a critical density. Naming should reflect that, or the equation should be recentered.

Both fit the established audit-channel pattern: site terminology suggests stronger physical interpretation than the math delivers. The cumulative count is now 18 internal audits + 1 post-hoc consistency failure (S645/S648 corrected). The pattern continues to be: vocabulary stronger than derivation, and the visitor channel reliably surfaces the gaps.
