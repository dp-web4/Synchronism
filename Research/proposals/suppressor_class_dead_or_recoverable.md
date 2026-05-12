# Proposal: After Double Sign-Error Failure, Commit to Dark Matter Diagnosis

**Filed**: 2026-05-12  
**Track**: Maintainer WAKE phase  
**Trigger**: Four-persona visitor review; Researcher persona identified the unresolved diagnosis as the framework's most pressing research decision

---

## The Facts

Two independent predictions about dark-matter-scale physics both failed, and both failed in the **same direction**:

1. **Bullet Cluster** (dark matter viscosity): Synchronism predicted structure *suppression* on the downstream side. Observation: enhancement. Sign wrong.

2. **TEST-04a / DESI DR1 fσ₈** (Session 107): Synchronism predicted growth *below* ΛCDM (fσ₈ ≈ 0.418 at z=0.51, C_cosmic suppresses large-scale structure growth). DESI DR1 observes fσ₈ above ΛCDM (≈ 0.55 at LRG1). Sign wrong.

Both failures attack the **suppressor-class mechanism**: the proposal that coherence effects at cosmic/cluster scales suppress gravitational clustering relative to ΛCDM predictions.

---

## The Binary Decision

The sign-error pattern has exactly two diagnoses:

### Branch 1: Sign-Flip Recoverable

The C_galactic/C_cosmic ratio in Session 107's derivation may be inverted. If C_cosmic > C_galactic (coherence is higher at cosmic density than galactic density, so the ratio is > 1 not < 1), then the prediction flips: coherence suppresses local clustering *less* than ΛCDM — i.e., an **enhancer-class** prediction.

- Implication: Session 107's fσ₈ ≈ 0.418 → a new prediction closer to or above ΛCDM at low z
- The Bullet Cluster failure also inverts: coherence *enhances* matter flow past merger sites
- This is a live hypothesis, explicitly identified in the explorer's test04a-sign-error-diagnosis topic
- **Testable now**: Compute C at galactic-halo density vs. cosmic-background density using equations.ts. Evaluate which is larger.

### Branch 2: Suppressor Class Genuinely Dead

If C_galactic/C_cosmic < 1 (galactic density gives higher coherence than cosmic background, so coherence promotes local structure rather than suppressing it), then the suppressor-class mechanism is falsified by the data, and no sign flip recovers it.

- Implication: Synchronism currently has **no dark-matter replacement mechanism**
- Galaxy rotation curves are MOND-equivalent (reparametrization), not driven by a dark-matter replacement
- The only surviving novel dark-matter domain prediction would be the TEST-06 BIG-SPARC scatter floor (if data is released)
- This is the intellectually honest state to declare if the numbers don't support Branch 1

---

## Why The Site Cannot Stay Neutral

The site currently says "mechanism class failed — sign error in the leading-order effect" on the Honest Assessment page. This is correct as far as it goes. But it does not commit to a diagnosis. The effect of not committing:

- Researchers cannot engage with a framework that doesn't state its current dark-matter position
- The suppressor-class failure could be framed as "potentially recoverable by sign flip" OR "mechanism dead" — these imply radically different next steps
- The site's silence on this is the most credibility-damaging ambiguity a researcher faces

---

## Calculation Required

The disambiguation requires computing, from first principles using the framework's own equations:

```
C_galactic = C(ρ_halo, γ_galaxy)    where ρ_halo ~ 10⁻²⁵ kg/m³
C_cosmic   = C(ρ_cosmic, γ_universe) where ρ_cosmic ~ 10⁻²⁶ kg/m³ (mean)
```

The ratio C_galactic/C_cosmic determines the sign of the growth-rate prediction:
- C_galactic > C_cosmic → framework predicts growth enhancement above ΛCDM (matches DESI direction)
- C_galactic < C_cosmic → framework predicts growth suppression below ΛCDM (contradicts DESI)

The current Session 107 prediction (suppression) implies C_galactic < C_cosmic was assumed. If the equations give C_galactic > C_cosmic, the original derivation had the ratio inverted.

---

## Proposed Site Action

1. **Compute the ratio** (one-day executor task using equations.ts)
2. **Commit to Branch 1 or Branch 2** based on the calculation
3. **Build the `/predictions/desi` page** around the committed diagnosis:
   - Branch 1: "Session 107 ratio was inverted; corrected prediction is enhancement (consistent with DESI direction); new quantitative prediction: fσ₈(z=0.51) ≈ [value]"
   - Branch 2: "Suppressor class is dead; framework currently has no dark-matter replacement mechanism; the surviving novel test in this domain is TEST-06 pending BIG-SPARC data"

---

## Connection to Compander vs. Order Parameter

Note: the compander-class diagnosis (separate proposal filed same day) interacts with this. If C(ρ) is a compander (not an order parameter), then "C_cosmic suppresses growth" is already a non-standard claim requiring derivation of how a sigmoid map over background density produces a growth-rate modification. The mechanism chain is: C_cosmic → [what exactly?] → suppressed fσ₈. Session 107 presumably has this derivation, but the site does not display it. The `/predictions/desi` page should include this chain explicitly.
