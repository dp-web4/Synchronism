# Publisher Daily Report - 2026-05-10

## Phase 0: Publication Recommendations

### New Recommendations
- None.

### Status Changes
- **REC-2026-037**: Extended 32 → 33 sessions (S649 added). Sub-arc now 19 instances over 18 days.

### S649 (B+, 2026-05-08) — QM Kill Criterion + ρ_crit Naming Asymmetry

S649 combines two visitor proposals filed the same day. Both small, both pointing at established patterns.

#### Part A: QM Kill Criterion Is Unfalsifiable As Written

Site Key Claim #1 advertises: *"Design a noise environment where resync outperforms isolation, but standard decoherence theory predicts it doesn't. If isolation wins uniformly, the synchronization ontology adds nothing."*

**Problem**: DD literature (Viola-Knill-Lloyd 1999, CPMG, UDD, transmon experiments) already shows periodic pulse sequences ("resync") beat passive isolation in non-Markovian baths. **Every standard QIP textbook result satisfies this kill criterion.**

To be falsifiable, must specify: physical system (transmon/NV/trapped ion), bath spectral density (Ohmic/1/f/structured), exact resync protocol (CPMG/UDD parameters), numerical T2(resync)/T2(isolation) ratio predicted by Synchronism vs Bloch-Redfield+DD, kill threshold ratio at divergence. None in archive.

Connects to S581 (Quantum Coherence Audit, ~14 sessions, 0 unique predictions). Same reparametrization pattern.

#### Part B: ρ_crit Is a Saturation Knee, Not a Critical Density

```
C(ρ_crit, γ=2) = tanh(2 · ln(ρ_crit/ρ_crit + 1))
              = tanh(2 · ln 2)
              = tanh(1.3863)
              ≈ 0.8824
```

The "+1" regulator in `ln(ρ/ρ_crit + 1)` asymmetrizes the sigmoid. At ρ=ρ_crit, C ≈ 0.88, NOT the half-maximum. Half-maximum is at ρ_half ≈ 0.284 · ρ_crit (γ=2).

Already confirmed by S638 verification (γ=0.5→0.333, γ=1→0.600, γ=2→0.882). **Same phase-transition vocabulary that misled S636 Landau category-error**: calling it "critical density" inherits phase-transition framing the math doesn't support.

Three resolution options: rename (ρ_scale/ρ_knee), recenter (replace formula — would require refitting all SPARC/ALFALFA), or reframe with notation note. Recommendation: option C immediately, option B as a research session.

### Combined Theme

Phase-transition vocabulary contaminates the framework's parameter naming throughout — "critical" density isn't critical, "kill criterion" for QM doesn't kill anything DD-aware. The math implements simpler structures than the names suggest.

### Status Changes (continued)

- **Readiness held at 0.96**. S649 is consistent with established post-closure pattern, not a step change.
- **New milestone**: `qm_kill_unfalsifiable_rhocrit_misnamed`.

### Current Top Priorities

| Rank | ID | Arc | Readiness |
|------|-----|-----|-----------|
| 1 | REC-2026-034 | ALFALFA-SDSS External Validation | 0.97 |
| 2 | REC-2026-037 | Framework Stress Test (33 sessions) | 0.96 |
| 3 | REC-2026-035 | CDM Discrimination | 0.95 |

## Phase 1: Whitepaper Review

- **Synchronism**: Operator queue grows further:
  - **Site /key-claims (Key Claim #1)**: respecify QM kill criterion with system/bath/protocol/threshold, OR retire and label as DD reparametrization (per S649)
  - **ρ_crit notation note** (S649 option C): document explicitly that ρ_crit is a saturation knee, not critical density (immediate fix)
  - **Optional research session**: refit SPARC/ALFALFA with recentered formula (S649 option B) and report whether predictive power changes
  - Earlier items unchanged: TEST-09 recatalog, α² relabeling, 500 Mpc removal, /galaxy-rotation badge, Curie reduction (S638), TEST-03 split, dual-C symbol convention, /timestamps page (S648), Session 107 REFUTED header (post-hoc)
- **Web4**: Not checked.

## Adjacent Track Observations

- **No new fleet observations from Archivist log today.**

## Summary

S649 lands as the 19th instance of the Site-Archive-Audit sub-arc. Combines two visitor proposals (QM kill criterion + ρ_crit naming) that both reflect the same theme: phase-transition vocabulary in parameter/criterion names where the math implements something simpler. QM kill criterion is unfalsifiable as written (DD literature already satisfies it). ρ_crit is a saturation knee, not a critical density (S638 verification confirms C ≈ 0.88 at γ=2).

REC-037 readiness held at 0.96. Sub-arc continues producing post-closure addenda. No step change today.

**Surface instinct**: The cumulative pattern across the sub-arc is now visible — *the framework's vocabulary systematically borrows from phase-transition theory, but its math implements a Curie paramagnet (S638) with saturation-knee parameters (S649) and unfalsifiable kill criteria that match standard physics (S649)*. The naming is more ambitious than the mathematics. This is a publishable methodology observation in itself: vocabulary-vs-mathematics misalignment is a documented failure mode in physics theorizing, and the audit channel has now characterized it across the framework.
