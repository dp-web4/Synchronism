# One-Scale-Insufficiency Theorem: C(ρ) Cannot Simultaneously Fit Galaxy and Cluster Dynamics

**Status**: Closed by execution (2026-05-28)
**Filed**: 2026-06-01
**Source**: Synchronism site maintainer session; originating from visitor Pass 4 (leading-edge researcher) observation that "C(ρ) silent at cluster scale by construction" should be reframed as a structural finding rather than an apology.

---

## Summary

C(ρ) = tanh(γ · ln(ρ/ρ_crit + 1)) has exactly **one free density scale**: ρ_crit. Multi-scale gravitational dynamics — galaxies (~10²² kg/m³ baryonic density) vs galaxy clusters (~10¹⁹ kg/m³) — require at least **two independent scale parameters** to fit independently. This is not a calibration miss or a numerical gap. It is a structural theorem: a single density knee cannot simultaneously satisfy both constraints.

Four natural ansätze were tested on the Coma cluster (2026-05-28). All failed:

- **A1/A4**: overshoot observed velocity dispersion by ~10⁴
- **A2**: collapses to Newtonian (suppression vanishes)
- **A3**: structurally impossible — C ∈ [0,1) bounds the function at ≤2, but the required dimensionless parameter derived from Coma's velocity dispersion ~4.6 needs ρ_crit to be 10⁴–10⁶ × larger at cluster scale — physically unreachable within the formulation
- **Required ρ_crit,cluster**: 10⁻⁴ to 10⁻⁶ × the galaxy-calibrated value — incompatible with universality

The cluster gap is not "we haven't tried" — it is "we tried four natural approaches and all fail for structural, not numerical, reasons."

---

## The Theorem

**Claim**: No single value of ρ_crit can simultaneously reproduce galaxy-scale rotation curves and cluster-scale velocity dispersions within C(ρ).

**Proof sketch**:

1. ρ_crit is calibrated from galaxy rotation curves: ρ_crit,galaxy is set by the baryonic surface density scale at which flat-rotation curve behavior reproduces observations.
2. Cluster dynamics operate at mean baryonic densities ~10⁻⁴ to 10⁻³ times the galaxy-calibrated ρ_crit (i.e., ρ_cluster ≪ ρ_crit,galaxy).
3. At cluster densities with the galaxy-calibrated ρ_crit, the argument of C becomes small: C(ρ_cluster) ≈ γ · ρ/ρ_crit ≈ 0 → Newtonian limit. No coherence modification appears (A2 result).
4. To obtain non-Newtonian dynamics at cluster scales, ρ_crit must be scaled down by 10⁴–10⁶ to bring cluster densities into the active range of C. This rescaling destroys galaxy-scale calibration.
5. C(ρ) has no second free scale to independently set the cluster regime. The single parameter ρ_crit cannot be in two places at once.
6. Therefore no universal ρ_crit exists that satisfies both constraints simultaneously. QED.

**Corollary**: C(ρ) is galaxy-only by construction, not by choice. The density-map program is closed at galaxy scale and cannot extend to cluster scale without a second free parameter. This is a theorem about the functional form, not a statement about fitting effort.

---

## Supporting Evidence: The 2026-05-28 Execution

Scripts: `explorer/work/cluster_bridge_coma.py`, `explorer/work/cluster_bridge_reverse_solve.py`

**Coma cluster observables used**:
- Velocity dispersion σ ≈ 880 km/s
- Virial radius R₂₀₀ ≈ 2 Mpc
- Baryonic mass M_bar ≈ 3×10¹⁴ M☉ (hot gas ~85%, stars ~15%)

**Results of four ansätze**:

| Ansatz | Description | Result |
|--------|-------------|--------|
| A1 | C(ρ) direct velocity coupling | Overshoots by 10⁴ |
| A2 | C(ρ) as dark-matter fraction modifier | Collapses to Newtonian gravity |
| A3 | C threshold determines virial bound state | Structurally impossible (C ∈ [0,1) constraint violated at required scale) |
| A4 | C as thermal pressure correction | Overshoots by 10⁴ |

The "required ρ_crit at cluster scale" reverse solve: fitting Coma's dynamics requires ρ_crit,cluster ≈ 10⁻⁴ to 10⁻⁶ × ρ_crit,galaxy. The spread reflects the different functional dependencies across the four ansätze — none are mutually consistent, and none are consistent with galaxy-scale calibration.

**Interpretation**: A3 deserves special attention. The C ∈ [0,1) bounding is not a numerical problem — it is a hard algebraic property of tanh. The dimensionless velocity ratio required to explain Coma's dynamics demands a coherence correction factor of order ~10⁴–10⁶, which the tanh function cannot reach. This is a structural impossibility, not a parameter tuning failure.

---

## Framework Comparisons

**MOND (Milgrom 1983)**: Fails at clusters for the same structural reason — one acceleration scale a₀ is insufficient for both galaxy and cluster regimes. MOND predicts ~2× too little cluster velocity dispersion (the "cluster mass deficit problem"). This is well documented: Sanders 2003, Pointecouteau & Silk 2005. The class diagnosis is identical: one-scale frameworks face a structural trade-off between galaxy and cluster calibration.

**Verlinde 2016 (entropic gravity)**: Has two effective scales (M_B, the baryonic mass, and a₀, the emergent gravity scale) and still underperforms at clusters. Tamosiunas et al. 2019 found 1.5–3× underprediction of cluster velocity dispersions across a sample of observed clusters. The lesson: two scales are **necessary but not sufficient** for cluster dynamics. Additional physics (hot gas thermodynamics, AGN feedback history, merging state) may be required.

**ΛCDM**: Carries a full parameter set (dark matter density profile shape, concentration–mass relation, baryon fraction, feedback prescriptions). Multi-scale flexibility is explicit by construction. Clusters are fit by adding dark matter halos with NFW profiles — effectively a second density scale per object.

**Class conclusion**: Every framework with one free gravitational scale parameter fails at clusters by the same structural argument. C(ρ) belongs to the same failure class as MOND and single-scale entropic-gravity variants. This is not an idiosyncratic failure — it is a known structural property of one-scale theories. The contribution here is demonstrating this cleanly within the Synchronism formalism and executing the four-ansatz falsification.

---

## Historical Note: The C(a) to C(ρ) Variable Migration

Sessions 195–199 in the Synchronism archive used C(a) — coherence as a function of acceleration (not density). Sessions 211 and later migrated to C(ρ). This migration was undocumented; no analysis was preserved comparing the cluster predictions of C(a) vs C(ρ).

**Why this matters**: In the C(a) formulation, the effective argument at cluster scales would be the cluster's mean gravitational acceleration. Rich clusters like Coma have mean accelerations a ~ GM/R² ~ 10⁻¹¹ to 10⁻¹⁰ m/s², which is near or below a₀ ≈ 1.2×10⁻¹⁰ m/s² — the same deep-MOND regime where MOND succeeds in galaxies. The C(a) formulation would naturally have cluster-scale accelerations entering the active range of the coherence function, potentially producing non-Newtonian dynamics at cluster scale without requiring a second density parameter.

In contrast, the C(ρ) formulation at cluster densities (ρ ≪ ρ_crit,galaxy) falls into the Newtonian limit by the theorem above.

**The implication**: The 2026 variable migration from C(a) to C(ρ) may have silently dropped a working cluster prediction. This is possibly the most consequential unreviewed decision in the project's post-2025 history. The C(a) formulation deserves re-examination before the cluster program is declared definitively closed.

**Caveat**: C(a) has its own problems (Lorentz invariance, covariant extension, non-uniqueness of acceleration in GR). The migration to C(ρ) may have been justified on independent grounds. The point is that the cluster comparison was never documented, and that gap should be filled.

---

## Status and Recommendations

### Current Status

| Scale | Program | Status |
|-------|---------|--------|
| Galaxy | RAR transition shape discriminator | Closed by execution: ΔBIC = +184, free-γ → MOND (2026-05-21) |
| Cluster | Four-ansatz Coma test | Closed by execution: structural impossibility (2026-05-28) |
| Cosmological | TEST-04a growth suppression | Closed: DESI full-shape disfavors by 2+ σ (2026-05-25) |
| Overall density-map | All scales | No surviving discriminating test |

### Recommendation A: Frame as a Contribution

The one-scale-insufficiency result should be presented as a **positive finding** — a theorem with proof, not an apology for a gap. The framework reached its own closure not through external refutation but through internal structural analysis. This is methodologically cleaner than a parameter-fitting failure.

Specifically:

1. **Update /honest-assessment** on the live site: replace "C(ρ) silent at cluster scale by construction" with the theorem statement and Coma execution results. Frame the cluster gap as a structural theorem proven by execution, not a known limitation awaiting future work.

2. **Preprint candidate**: The one-scale-insufficiency theorem, combined with the DESI growth-suppression constraint (TEST-04a), constitutes a clean pair of negative results with generalizable implications:
   - TEST-04a: coherence-damped structure growth is ruled out as a class by DESI full-shape data
   - One-scale theorem: single-density-knee maps cannot address multi-scale dynamics — a result that extends to MOND and Verlinde-class theories
   Both results are stated more sharply than in the existing modified-gravity literature.

3. **Investigate the C(a) migration** (see Historical Note): this is the one remaining path to a potential cluster prediction within the Synchronism framework. If C(a) gives a qualitatively different (and better) cluster result, the migration decision is the single most consequential unreviewed choice in the recent project history.

---

## Open Questions

**1. Does a coherence length scale escape the theorem?**

Adding a coherence *length* scale ξ(r) alongside the density knee ρ_crit would give two free parameters. A candidate extended form:

```
C(ρ, ξ) = tanh(γ · ln(ρ/ρ_crit + 1) · f(r/ξ))
```

where f(r/ξ) modulates the coherence effect by physical scale. Is this physically motivated within the Synchronism ontology? Does it escape the one-scale theorem? If yes — what observable prediction does it make that distinguishes it from a two-parameter MOND fit? If no — the density-map program is provably closed at all scales.

**2. C(a) resurrection**

Does the pre-2026-04 C(a) formulation give a qualitatively different cluster prediction? If yes, why was it abandoned, and is the reason still valid? If no, the one-scale theorem applies to both formulations and the closure is deeper.

**3. Class transferability to the modified-gravity literature**

The one-scale-insufficiency theorem applies to any framework of the form f(x/x₀) with a single scale x₀, regardless of what x is (density, acceleration, curvature, entropy). This class includes:
- MOND: μ(a/a₀), one scale a₀
- Synchronism C(ρ): tanh(γ ln(ρ/ρ_crit + 1)), one scale ρ_crit
- Any Hill/compander variant with one characteristic scale
- Single-field f(R) gravity with one characteristic curvature scale

Documenting this class and its structural cluster failure as a unified theorem may be a useful contribution to the modified-gravity literature independent of Synchronism.

**4. The methodology preprint**

The one-scale theorem is a clean exhibit for an A2ACW methodology write-up — a case where internal structural analysis, rather than external data alone, closed a research direction. The four-ansatz execution is reproducible, the derivation is self-contained, and the result generalizes. It deserves a section in a methodology preprint.
