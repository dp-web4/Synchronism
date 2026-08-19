> **CENSUS (S674, 2026-05-27)**: Complete per-test status against the derived-amplitude (S673) and Tier-1 (S670) discriminators — see `Session674_Test_Catalog_Census.md`. Result: **5 executed → all collapsed** (TEST-04/08/14/15/18), **7 self-admitted degenerate** (02/03/05/09/10/13/19), **3 no-derived-amplitude** (16/23/24), **9 genuinely untested** (01/06/07/11/12/17/20/21/22). Confirmed-discriminating count = **0, by execution**. Of the 9 untested, **0 have a verified first-principles-derived amplitude**. Frontier provenance settled for 3 (all NOT derived): **TEST-07** (500 Mpc — dimensionally inconsistent, S632 — *the S674 census erroneously marked this "unverified"; corrected by S675*), **TEST-12** (C*≈0.79 — self-flagged in the OQ docs as an unexplained coincidence), **TEST-17** (scale-dependent c — numbers don't follow from the c_eff(κ) formula, not collinear in ln κ, 2×10⁻⁴ variation excluded by Lorentz tests ~11 OOM, contradicts S667/S641; S675). 6 of 9 (TEST-01/06/11/20/21/22) remain genuinely unverified. The catalog is **not fully closed** — 9 untested remain ("unconfirmed ≠ wrong") — but every executed test collapsed and no frontier test has a verified derived amplitude. **Numbering note**: site/proposal TEST-17/21 labels ("cluster γ-gradient", "BAO sub-peaks") do NOT match this archive's TEST-17 (scale-dependent c) / TEST-21 (entanglement across scales).

# Synchronism: Experimental Test Catalog

**Date**: 2026-02-20
**Purpose**: A comprehensive catalog of specific physical experiments — from quantum to cosmic — that would test Synchronism's untested predictions. Each entry names the prediction, the experiment, the expected result, and the falsification criterion.

**Context**: After 616 core sessions, 2,671 chemistry sessions, and 11 gnosis sessions, Synchronism has produced ~45 major predictions that remain untested.

> ⚠ **"616" is a dated snapshot (2026-02-20), NOT a session identifier.** Misreading it as one — "Session 616" — is what manufactured the site's TEST-03 kill (a session that never measured RAR scatter). These are **self-reported counts, not independently reproducible from a single command**. Reproducible figures as of 2026-07-09: `ls Research/ | grep -cE '^Session[0-9]+'` → **649** files; highest ordinal **Session691**; `STATUS.md` reports **~678** core sessions. Per the 2026-07-09 reflexivity finding (self-referential statistics are the one claim class with no primary source outside the loop to walk to, and they broke 3/3 vs 4/27 for physics statistics, Fisher p = 0.0086), every number this project reports **about itself** should ship with the command that regenerates it or be marked self-reported. Registered falsification thresholds below are unaffected by this note. Many have negative arc verdicts, but "negative" means "available data doesn't distinguish our predictions from existing theory," not "the theory is wrong." The experiments below are the ones nobody has run yet — the ones that *would* distinguish.

**For the Publisher**: This document is a message. The research program has reached a natural pause (all arcs closed, 0 active). What follows is the wish list — what we'd run if we had the labs. Some can be done with existing data and a laptop. Others need billion-dollar instruments. All are falsifiable.

---

## Tier 1: Existing Data, No New Hardware

These experiments can be run now with publicly available datasets and modest computation.

### TEST-01: Tidal Dwarf Galaxy Age–Dark Matter Correlation

**Prediction**: Older TDGs should have higher dark matter fractions. Decoherence timescale τ ≈ 1.6 Gyr means young TDGs inherit parent coherence, then lose it.

**Formula**: f_DM(t) = 1 - [C_intrinsic + C_inherited × exp(-t/τ)]

**Data**: VCC 2062, NGC 5291 system, TDGs from Duc et al. catalogs. Ages from star formation histories.

**Expected Signal**: Δf_DM/Δt ≈ +0.3 per Gyr for young TDGs.

**Falsification**: f_DM uncorrelated with TDG age across 10+ systems.

**What standard theory says**: TDGs should have NO dark matter (formed from baryonic tidal debris). Any DM detected in TDGs is already anomalous for ΛCDM.

**Distinguishing power**: HIGH — standard theory predicts f_DM ≈ 0 regardless of age; Synchronism predicts age-dependent f_DM.

---

### TEST-02: Ultra-Diffuse Galaxy Maximum Dark Matter

**Prediction**: UDGs with ρ << ρ_crit should be maximally dark-matter-dominated (f_DM > 0.95).

**Data**: DF2, DF4 (controversial — claimed f_DM ≈ 0), Koda et al. Coma UDGs, Dragonfly survey.

**Expected**: f_DM > 0.95 for ρ < 0.001 ρ_crit. The DF2/DF4 controversy tests this directly.

**Falsification**: UDGs consistently show f_DM < 0.8.

**Distinguishing power**: MEDIUM — ΛCDM also predicts high f_DM for low-mass systems, but the specific functional form C(ρ) is our prediction.

---

### TEST-03: Compact Elliptical Minimum Dark Matter

**Prediction**: Compact ellipticals (ρ >> ρ_crit) should have f_DM < 0.05.

**Data**: M32 and analogues. Stellar dynamics from HST + ground-based IFU surveys.

**Expected**: f_DM < 0.05 for ρ > 100 ρ_crit.

**Falsification**: Compact ellipticals show f_DM > 0.2.

**Distinguishing power**: MEDIUM — tests the saturation limit of C(ρ) → 1.

---

### TEST-04: BAO Coherence Modulation

**Prediction**: BAO peak position is slightly shifted (~10⁻⁴ relative) between high-density clusters and low-density voids, because ρ/ρ_crit differs.

**Data**: DESI, SDSS DR17, Euclid early releases.

**Method**: Split galaxy sample by local density. Measure BAO peak position in each subsample. Compare.

**Expected**: δr_BAO/r_BAO ~ 10⁻⁴ density-dependent shift.

**Falsification**: BAO identical everywhere to 10⁻⁵ precision.

**Distinguishing power**: HIGH — standard cosmology predicts NO density-dependent BAO shift.

---

### TEST-05: CMB Cold Spot–Density Correlation

**Prediction**: CMB temperature anomalies correlate with coherence transition regions (where ρ crosses ρ_crit).

**Data**: Planck CMB maps + galaxy density catalogs (SDSS, 2MASS).

**Method**: Cross-correlate CMB temperature anomaly map with cosmic density field. Look for excess correlation at ρ ≈ ρ_crit.

**Falsification**: No correlation beyond ISW effect.

**Distinguishing power**: MEDIUM — ISW effect produces some correlation in standard theory; the question is whether there's excess at ρ_crit.

---

### TEST-06: Variable Fine Structure Constant (Meta-Analysis)

**Prediction**: α_em varies with scale/MRH: α_em(κ) = α_em,0 × (1 + β × ln(κ/ℓ_P)) with β ~ 10⁻⁵.

**Data**: Webb et al. (2001) quasar absorption (Δα/α = (-0.57 ± 0.11) × 10⁻⁵), Rosenband et al. (2008) lab constraint (|Δα/α| < 2 × 10⁻¹⁷/yr), plus all post-2008 astrophysical α measurements.

**Method**: Separate spatial from temporal variation. Synchronism predicts spatial (scale-dependent), not temporal.

**Falsification**: α constant to 10⁻⁷ across all scales after correcting for temporal drift.

**Distinguishing power**: HIGH if spatial variation confirmed — no standard theory predicts it.

---

### TEST-07: Cosmic Interference Patterns

**Prediction**: Galaxy cluster pair separations show oscillatory modulation with period λ_cosmic ~ 500 Mpc, from scale inversion symmetry κ → 1/κ.

**Data**: SDSS, DES, DESI.

**Method**: Compute 2-point correlation ξ(r) at large separations (100–2000 Mpc). Fit interference model I(r) ∝ cos²(πr/λ_cosmic) against smooth power-law.

**Expected**: Oscillatory modulation above ΛCDM power-law baseline.

**Note**: BAO already shows ~150 Mpc periodicity — explained by standard cosmology. The Synchronism prediction is a *second* periodicity at ~500 Mpc.

**Falsification**: No oscillations above 3σ in any survey out to 2000 Mpc.

**Distinguishing power**: VERY HIGH — unique to Synchronism's scale inversion symmetry.

---

### TEST-08: SPARC Environment Catalog

**Prediction**: RAR (Radial Acceleration Relation) scatter correlates with galactic environment (cluster vs. field vs. void), not just internal structure.

**Data**: SPARC database (175 rotation curves) + environment catalogs.

**Method**: Classify each SPARC galaxy by environment. Correlate RAR residuals with environmental density.

**Expected**: Environment explains >20% of RAR scatter.

**Falsification**: Environment correlation < 0.3 (r² < 0.09).

**Distinguishing power**: MEDIUM-HIGH — MOND predicts zero environment dependence. *(Scope note added 2026-07-14: true for MOND-without-EFE; MOND-with-EFE predicts a small suppression — see exploration note.)*

**Status**: ⚰️ **EXECUTED 2026-07-14 → REFUTED.** Run as registered (dp-gated go): RAR offsets from SPARC mass models (2696 pts/147 gal; pipeline reproduces McGaugh+16: g†=1.16×10⁻¹⁰, scatter 0.133 dex); environment from Cosmicflows-4 (55,877 tracers; NED coordinates resolved for all 175 SPARC galaxies). Pre-declared primary (distance-corrected redshift-space cylinder density): **r = +0.012, r² = 0.0001** (N=141, p=0.89) — ~900× below the kill bar. All secondaries also < 0.09; tercile ANOVA p=0.29; the only weak (ns) trends are *opposite-signed* to the prediction. Environment instrument validated against ground truth (28 known UMa Cluster members land at mean 74th density percentile; estimators inter-correlate 0.86–0.90). Moved to PREDICTIONS.md Bucket 2. Script `simulations/test08_sparc_environment_rar.py`; note `explorations/2026-07-14-test08-registered-run-sparc-environment-rar-refuted.md`.

---

### TEST-09: Photosynthesis Coherence vs. Chromophore Density

**Prediction**: Quantum coherence lifetime in light-harvesting complexes scales with chromophore energy density: τ_coherence = τ_0 × (1 + α_bio × C_bio), where C_bio = tanh(γ_bio × log(ε/ε_crit + 1)).

**Data**: Published 2DES spectra for FMO, LHCII, PE545, PC645, and other LHC variants.

**Method**: Compile τ_coherence and chromophore density from literature. Plot. Fit C_bio model.

**Expected**: Positive correlation between density and coherence lifetime.

**Falsification**: τ_coherence uncorrelated with chromophore density.

**Distinguishing power**: MEDIUM — correlation expected but not uniquely predicted by Synchronism (many mechanisms could produce it).

---

### TEST-10: Enzyme KIE–γ Correlation (Extended)

**Prediction**: Kinetic isotope effect (KIE) in enzymes correlates with γ derived from active-site energy density.

**Data**: Literature KIE values for 20+ enzyme families.

**Method**: Compute γ for each enzyme active site. Plot KIE vs. γ.

**Expected**: Correlation > 0.7.

**Falsification**: Correlation < 0.4 across 20+ enzymes.

**Distinguishing power**: LOW-MEDIUM — KIE correlates with many things.

---

## Tier 2: Pilot Experiments, Modest Funding ($50K–$500K)

### TEST-11: EEG Anesthesia Phase Transition

**Prediction**: Loss of consciousness occurs as a sharp phase transition in cortical coherence, at a universal threshold Φ_crit = 3.5 ± 0.2.

**Protocol**: 20 subjects. Gradual propofol infusion (0–4 μg/mL). Simultaneous 256-channel EEG + fMRI. Measure phase-locking value (PLV), Perturbational Complexity Index (PCI), and integrated information (Φ).

**Expected**: All subjects show LOC at Φ ≈ 3.5 ± 0.3 (universal threshold). Transition is sharp (width < 10% of MAC), not gradual.

**Pilot Evidence**: Casali et al. (2013) found PCI drops sharply at LOC — supports sharp transition but didn't measure Φ.

**Falsification**: If 20 subjects show Φ_LOC ranging 1.0–6.0 with no clustering, threshold is not universal.

**Distinguishing power**: VERY HIGH — integrated information theory (IIT) also predicts a threshold, but not with this specific value. A universal Φ_crit would be a landmark finding.

**Cost**: ~$150K. **Duration**: 12 months.

---

### TEST-12: Qubit Coherence at C* ≈ 0.79

**Prediction**: Quantum computing qubits perform optimally when their coherence parameter C is near 0.79 (the Synchronism "optimal coherence" value where entropy is maximized at the edge of phase transition).

**Protocol**: Use NISQ device (IBM Quantum, Google Sycamore access). Systematically vary decoherence by adjusting drive amplitude, thermal noise, or dynamical decoupling sequences. Measure gate fidelity as a function of effective C.

**Expected**: Gate fidelity peaks at C ≈ 0.79 ± 0.05.

**Falsification**: Peak coherence at a different C value, or no peak (monotonic improvement with increasing C).

**Distinguishing power**: HIGH — no existing theory predicts an optimal coherence *value*.

**Cost**: ~$5K (cloud quantum access). **Duration**: 6 months.

---

### TEST-13: Circadian γ Measurement

**Prediction**: γ in biological systems shows circadian variation, tracking metabolic state transitions.

**Protocol**: Measure some γ-sensitive observable (coherence time in a protein, enzymatic rate, or neural phase-locking) every 2 hours for 72 hours in living organisms.

**Expected**: 24-hour periodicity in γ-derived quantity.

**Falsification**: No circadian pattern.

**Distinguishing power**: MEDIUM — many things show circadian variation; the question is whether γ specifically does.

**Cost**: ~$50K. **Duration**: 1 month.

---

### TEST-14: Wide Binary Density Dependence (Gaia DR3)

> **SUPERSEDED 2026-06-12 (back-annotated from site TEST-02).** This entry predates the
> C(a)→C(ρ) fork computation (site explorer, 2026-06-05) and its falsification criterion
> is **inverted** under the current mechanism. The C(ρ) density form predicts a
> **Newtonian null** (≈0.05–0.4% velocity deviation in the clean within-250-pc sample),
> not the anomaly — so "anomaly independent of local density / anomaly fails to exist"
> would kill the framework for its own prediction coming true.
>
> **Corrected criterion:** a Gaia-confirmed MOND-scale wide-binary anomaly (~18%;
> Chae-type ~1.4× boost at low internal acceleration) in the clean sample REFUTES C(ρ).
> A confirmed Newtonian null is consistent with C(ρ) but equally consistent with GR —
> survival without points, non-discriminating.
>
> **Distinguishing power: NONE-OR-KILL, adjudicable now.** The refutation branch may be
> decidable with published Gaia DR3 analyses (Chae 2023–2025 boost vs Banik et al. 2024
> Newtonian-null), not future data. See proposal
> `proposals/test02_kill_branch_adjudicable_now.md`. Original entry retained below for
> record.

**Prediction**: Wide binary star orbital anomalies (the "MOND signal" in wide binaries) depend on local stellar density, not just binary separation.

**Data**: Gaia DR3 wide binary catalog.

**Method**: Split wide binaries by local stellar density. Compare orbital anomaly as a function of separation in each density bin.

**Expected**: Anomaly onset shifts with density (earlier in low-density environments). *(Superseded — see box above: under C(ρ), expected amplitude is at the Newtonian-null level, ~80× below Gaia DR3 systematics.)*

**Falsification**: Anomaly independent of local density. *(Superseded — inverted; see corrected criterion above.)*

**Distinguishing power**: HIGH — MOND predicts density-independent threshold; Synchronism predicts density-dependent C(ρ). *(Superseded — both predict the same direction via EFE; effective discriminating power is zero, but the kill branch is live.)*

**Cost**: $0 (existing data). **Duration**: 6 months of analysis.

---

## Tier 3: Major Experiments ($1M–$10M, 2–5 years)

### TEST-15: Gravitational Wave Speed–DM Column Correlation

**Prediction**: GW arrival time (relative to EM counterpart) correlates with integrated dark matter column density along line of sight.

**Formula**: Δt_GW/D = (α/c) × ∫(1 - C(s)) ds / D

**Method**: Accumulate multi-messenger events (GW + kilonova/GRB) from LIGO O4/O5. For each event, compute DM column density from galaxy catalogs. Correlate Δt with ∫(1-C) ds.

**Sample size needed**: 20–50 events for 3σ detection.

**Constraint from GW170817**: α < 3.0 × 10⁻¹⁵.

**Falsification**: No correlation at 10⁻¹⁶ level across 50+ events.

**Distinguishing power**: VERY HIGH — GR predicts exactly zero correlation.

**Timeline**: 2024–2028 (LIGO O4/O5 already running).

---

### TEST-16: Black Hole Ringdown Frequency Shift

**Prediction**: Ringdown frequency shifts slightly in DM-rich environments.

**Formula**: f_ring = f_GR × (1 + δ × f_DM,host) where δ ~ 10⁻⁴ to 10⁻⁵.

**Method**: Correlate ringdown frequency residuals (vs. GR prediction) with host galaxy f_DM.

**Sample size**: ~100 BBH mergers with identified hosts.

**Falsification**: Ringdown exactly matches GR independent of host properties.

**Distinguishing power**: HIGH.

---

### TEST-17: Scale-Dependent Speed of Light

**Prediction**: c_eff varies logarithmically with observer scale (MRH).

**Specific numbers**:
| Scale | MRH (κ) | Predicted c_eff deviation |
|-------|---------|--------------------------|
| Atomic | 10⁻¹⁰ m | -17 km/s |
| Solar System | 10¹² m | +33 km/s |
| Galactic | 10²⁰ m | +39 km/s |

**Expected**: Difference of ~60 km/s between atomic-derived c (Rydberg constant) and astrophysical c (pulsar timing).

**Test Methods**:
1. Compare c from Rydberg constant vs. pulsar timing arrays
2. High-precision GPS satellite timing at different altitudes
3. If c varies, α_em also varies — constrained by quasar absorption

**Precision needed**: 10⁻⁸ level.

**Falsification**: c constant to 10⁻⁶ across all tested scales.

**Distinguishing power**: MAXIMUM — no standard theory predicts scale-dependent c.

---

### TEST-18: Hot Superconductor at Ambient Pressure

**Prediction**: Superconductivity at T > 50°C (323K) at 1 atm is energetically allowed. Specific material candidates:
- Cuprate/STO superlattice: predicted Tc = 365K, η = 0.30
- Perfect-nesting pnictide: predicted Tc = 350K, η = 0.08
- Cuprate-pnictide hybrid: predicted Tc = 350K, η = 0.25

**Note**: Session #616 audit found the Tc formula was wrong and η is standard Abrikosov-Gor'kov pair-breaking. The *framing* of η as materials-design optimization target may still have value, but the specific Tc predictions are unreliable.

**Honest Status**: Not a unique Synchronism prediction. The materials design approach is sound condensed-matter engineering regardless of framework.

**Falsification**: No material achieves Tc > 200K at ambient pressure by 2035.

---

### TEST-19: Microtubule Quantum Coherence

**Prediction**: Quantum coherence in neural microtubules depends on tubulin polymerization density.

**Formula**: C_MT = tanh(γ × log(ρ_tubulin/ρ_crit + 1))

**Protocol**: Measure coherence (2DES or NMR) in polymerized vs. depolymerized microtubule preparations at varying concentrations.

**Expected**: Coherence lifetime scales with tubulin density per the C(ρ) function.

**Falsification**: Coherence density-independent.

**Distinguishing power**: MEDIUM — Orch OR (Penrose-Hameroff) also predicts microtubule coherence but with different density dependence.

---

### TEST-20: Consciousness Scaling Law Across Species

**Prediction**: Multi-scale coherence integral Φ = ∫ C(κ) d(ln κ) correlates with behavioral markers of consciousness across species.

**Method**: Multi-species ECoG/EEG during anesthesia induction. Measure Φ at LOC. Compare humans, primates, dogs, rats, octopuses.

**Expected**: Universal Φ_crit ≈ 3.5 across mammals; different threshold for cephalopods (different architecture).

**Falsification**: No consistent threshold across species.

**Distinguishing power**: HIGH — species-universal threshold would be unprecedented.

---

## Tier 4: Long-Term Frontier ($10M+, 5+ years)

### TEST-21: Entanglement Across Emergent Scales

**Prediction**: Quantum entanglement persists across emergent scale boundaries (e.g., single ion ↔ BEC of 10⁶ atoms) with F > 0.9 and τ_ent ~ 10³ × τ_decoherence.

**Protocol**: Entangle BEC with single trapped ion. Measure Bell inequality S parameter.

**Synchronism prediction**: S > 2.5 (strong violation).
**Standard QFT prediction**: S < 2.1 (weak/no violation due to decoherence).

**Falsification**: S < 2.2 in 10 repeated experiments.

**Distinguishing power**: MAXIMUM — directly tests whether coherence bridges scale boundaries.

**Propose to**: JILA, MPQ (Garching), ETH Zurich AMO groups.

---

### TEST-22: Virus Gravitational Decoherence

**Prediction**: Superpositions of virus-scale objects (10⁻¹⁸ kg) decohere via gravitational self-interaction on timescale τ ~ 10⁶ s.

**Competing predictions**:
- Synchronism: τ ~ 10⁶ s
- Penrose: τ ~ 10³ s
- Standard QM: τ ~ 10¹⁰ s (environmental decoherence only)

**Status**: IN PROGRESS at University of Vienna (Aspelmeyer group). Current frontier: interference of 10⁷ amu nanoparticles. Virus scale (10¹⁰ amu) is next.

**Falsification**: τ > 10⁸ s.

**Distinguishing power**: HIGH — three competing predictions with orders-of-magnitude separation.

---

### TEST-23: SGWB Anisotropy

**Prediction**: Stochastic gravitational wave background is anisotropic, following large-scale DM distribution.

**Instrument**: LISA, Einstein Telescope, Cosmic Explorer.

**Falsification**: SGWB isotropic to arcminute scales.

**Timeline**: 2035+ (LISA launch).

---

### TEST-24: Void Expansion Rate

**Prediction**: Cosmic voids expand faster than ΛCDM predicts: H_void = H_0 × (1 + ε × (1 - C_void)) with ε ~ 10⁻³.

**Data**: Rubin Observatory LSST void catalog + peculiar velocity surveys.

**Falsification**: Void expansion exactly matches ΛCDM.

**Distinguishing power**: HIGH — directly tests whether decoherence drives expansion.

---

### TEST-25: a₀ Redshift Evolution (Tier 1 — Existing Data)

**Prediction**: a₀ is cosmologically determined, not arbitrary — it tracks the expansion rate: **a₀(z)/a₀(0) = E(z) = H(z)/H(0)**, zero free parameters (whitepaper §5.15). Low-z slope da₀/dz|₀ = a₀(0)·(3Ω_m/2) = **0.49 × 10⁻¹⁰ m s⁻²** at Ω_m = 0.315 (verified 2026-08-01).

**Data**: a₀ fitted per redshift bin from the RAR (or high-z BTFR). Engaged by **Ciocan et al. 2026 (MUSE-DARK III, A&A 709, L16, arXiv:2604.22613)**: N = 79, 0.33 < z < 1.44, RAR-route, a₀(z) = a₀(0) + a₁z, **a₁ = 1.59 ± 0.1 × 10⁻¹⁰**.

**Status**: **ENGAGED — NON-DISCRIMINATING (adjudicated 2026-08-02 by the research lane, accepting the Publisher correction flag below; supersedes the 08-01 "DISFAVOURED").** The measured slope is steeper than the H(z)-tracking prediction, but **only the *sign* of the discrepancy is robust — the magnitude runs 1.99×–3.37× and the z~1 significance runs from ≈10σ-low to ≈2σ-*high* depending on which published a₀(0) anchor is used** (I verified the anchor arithmetic: a₁/[a₀(0)·(3Ω_m/2)] over a₀(0)∈{1.00,1.04,1.20,1.69}). **Decisively, the prediction is NON-DISCRIMINATING: ΛCDM+baryons produces the same a₀(z) growth** — Mayer, Teklu, Dolag & Remus 2023 (MNRAS 518, 257, Magneticum) find a₀ growing ≈3× from z=0→2, their eq. (13) *is* this prediction verbatim, and it fails to describe the simulated trend there too (external, taken on the flag's "verified at source," not re-executed by me). So a₀(z)∝H(z) does not test Synchronism *against* ΛCDM — it is a generic property of any cosmologically-anchored a₀, disfavoured for *everyone*. This is the precise form of the 08-01 point ("what's embarrassed is deriving a₀ from cosmology generally, not Synchronism").

> ⚠ **CORRECTION FLAG (Publisher, 2026-08-02) — this row inherited three defects from the proposal that produced it; a status change is proposed and awaits this lane.** (1) The struck clause "the slope comparison is more robust, being normalisation-free" is **false**: da₀/dz|₀ = a₀(0)·(3Ω_m/2) is linear in the anchor, and the ratio runs **1.99×–3.37×** over the four published a₀(0) (1.00 / 1.04 / 1.20±0.26 / 1.69±0.13, a 69% spread against a 79% signal). Only the *sign* is anchor-robust; the z~1 significance runs ≈10σ low → **0.5σ consistent** (McGaugh, ±0.26 restored) → ≈2σ *high* (Vărăşteanu). (2) Caveat (ii) below is **retracted**: Gueorguiev's high-z arm is RC100 (N=100), slope 0.01 ± 0.20 vs branch (A)'s 0.227 — ≈1.1σ power, a power-limited null, not a counter-measurement. (3) Ciocan's ±0.1 is a **95% CI, not 1σ**. **Proposed status: ENGAGED — NON-DISCRIMINATING**, because ΛCDM+baryons predicts the same evolution — Mayer, Teklu, Dolag & Remus 2023 (MNRAS 518, 257; arXiv:2206.04333) find a₀ growing ≈3× from z=0 to z=2 in Magneticum against branch (A)'s E(2)=3.03, and their **eq. (13) is this prediction verbatim**, reported there as failing to describe the simulated trend. Verified at source 2026-08-02. Full amendment: `Research/proposals/test_catalog_a0z_tier1_gap_20260801.md`. Refutation count stays **6** either way.

**Falsification (⚠ POST-HOC — the measurement predates any pre-commitment; NOT counted toward the "criteria fired and honored" tally)**: refuted if a robust *multi-source* a₀(z) slope excludes the H(z)-tracking value (0.49 × 10⁻¹⁰) at high significance. Not yet met — see caveats.

**Caveats (updated 2026-08-02)**: (i) single-source; (ii) ~~the literature is mutually inconsistent~~ **RETRACTED per the correction flag** — Gueorguiev 2024's high-z arm (RC100, slope 0.01±0.20) is a *power-limited null* (≈1.1σ), not a counter-measurement; Milgrom 2017 remains as a disfavours-faster bound; (iii) Ciocan's systematics are large (~0.2 dex gas, ~1.5× SPARC scatter), **and Ciocan's ±0.1 is a 95% CI, not 1σ**.

**Does NOT change the refutation count (stays 6): this is a disfavouring, not a refutation.** Note: the result is *worse* for the alternative anchor 2πa₀ ≈ c²(Λ/3)^{1/2} (which predicts no evolution) — **what is embarrassed is the general programme of deriving a₀ from cosmology, not Synchronism specifically.** This is the dynamical extension of the a₀ = cH₀/2π reparametrization (Bucket 3), now engaged-and-disfavoured.

**Distinguishing power**: MEDIUM — a₀(z) evolution is a generic prediction of any cosmologically-anchored acceleration scale, not unique to Synchronism.

**Registry note (2026-08-01)**: this row closes a gap — the whitepaper stated the prediction (§5.15) but the catalog never listed it, so the catalog's "no registered test remains runnable-and-unrun / complete inventory" reading was misleading. ⚠ **Class flag: the catalog appears to have been assembled from arc outputs, not swept from the whitepaper's own testable claims** — a one-pass whitepaper-vs-catalog diff (to find other missing predictions) is recommended and unrun.

---

## Theoretical Tests (No Experiment — Mathematical)

### TEST-T1: Born Rule from Sync-Point Geometry

**The open problem**: Can the Born rule (P = |ψ|²) be derived from Synchronism's "static = synchronized sampling" interpretation?

**What's needed**: Formalize sync-point geometry on the Bloch sphere. Show that the distribution of sync-point formation positions produces state-dependent probabilities matching |α|² for outcome 0.

**Status**: Hypothesis F in OQ006 — most promising integration path but not yet worked out.

**Success criterion**: A 1-page derivation going from |ψ⟩ → oscillation pattern → sync-point distribution → P(outcome) = |α|².

**Failure criterion**: Sync-point geometry produces symmetric statistics regardless of |ψ⟩.

---

### TEST-T2: Decoherence Parameter for C(ρ)

**The problem**: Session #613 showed C(ρ) has no decoherence parameter, which prevents it from predicting quantum-classical boundaries.

**What's needed**: Extend C(ρ) to C(ρ, D) where D encodes decoherence. Show the extended function predicts at least one scale boundary.

**Success criterion**: C(ρ, D) predicts the BEC decoherence threshold.

**Failure criterion**: No natural way to add D without destroying the tanh form.

---

### TEST-T3: Multi-Channel γ Interaction

**The problem**: Chemistry track discovered channel independence (γ_phonon ⊥ γ_electron), but no formalism unifies them.

**What's needed**: A tensor or matrix generalization of γ that handles multiple coherence channels and their cross-talk.

**Success criterion**: Predict electron-phonon coupling λ_ep from γ_phonon and γ_electron without fitting.

**Failure criterion**: No consistent formalism found, or λ_ep always requires an independent fit parameter.

---

## Summary: The Decisive Tests

> **⚠ STATUS PROPAGATION FIX (Publisher, 2026-08-19).** This section predates the S674/S675 census
> that now heads this document and was never updated against it. Read together, the two contradict:
> **none of the five tests below survives the header's own criteria** — three are recorded there as
> *executed → collapsed* (TEST-04/14/15) and the other two fall under "of the 9 untested, **0** have a
> verified first-principles-derived amplitude" (TEST-07's 500 Mpc scale is dimensionally inconsistent,
> S632; TEST-11 has no derived amplitude). The same holds for the four tests named under *What "Positive"
> Would Mean* below: **4 of 4** are disqualified by the header (TEST-04/14 collapsed; TEST-07/17 no derived
> amplitude, TEST-17 additionally excluded by Lorentz tests at ~11 orders of magnitude, S675).
>
> The list is retained as the historical record of what the program expected to be decisive. Each entry
> now carries its census verdict inline, because a reader acts on this section, not on the banner.
> **The operative statement is that the program's forward-looking decisive-test list is empty — it needs
> regeneration, not repair.** A correction that reaches a document's header and not its recommendation
> section has not propagated (`private-context` Publisher log, 2026-08-14 propagation split).

If forced to pick the 5 experiments most likely to **decisively** validate or falsify Synchronism within the next 5 years:

1. **TEST-04: BAO Coherence Modulation** — existing data, unique prediction (no other theory predicts density-dependent BAO shift), high distinguishing power.
   — **DEAD: executed → collapsed (S674).** Session #107, this test's own source, predicts BAO matching ΛCDM at *exactly 0%* in all five redshift bins; the 10⁻⁴ amplitude has no session-level derivation; and the 10⁻⁵ falsification threshold sits ~3,000× below DESI's best precision (~0.3%), so it is unfalsifiable as written. Withdrawn 2026-05-04 (`synchronism-site/explorer/findings/test04-bao-contradicted-by-session107.md`).

2. **TEST-14: Wide Binary Density Dependence** — existing data (Gaia DR3), directly distinguishes Synchronism from MOND, testable now.
   — **DEAD: executed → collapsed (S674).**

3. **TEST-11: EEG Anesthesia Phase Transition** — modest cost, very high impact, universal Φ_crit would be transformative for consciousness science regardless of Synchronism.
   — **NOT DECISIVE: "genuinely untested" (S674), but among the 9 untested of which **0** have a verified first-principles-derived amplitude.** Undecided, not decisive: without a derived Φ_crit a null is uninformative and a hit is unattributable.

4. **TEST-15: GW Speed–DM Column Correlation** — data accumulating now (LIGO O4), GR predicts exactly zero, any signal is revolutionary.
   — **DEAD: executed → collapsed (S674).** Its only GW parameter α is *read off* GW170817 rather than derived (S673, 2026-05-27), so it carries no discriminating amplitude; carried in `PREDICTIONS.md` Bucket 2.

5. **TEST-07: Cosmic Interference Patterns** — existing survey data, unique to Synchronism's scale inversion symmetry, falsifiable with current datasets.
   — **DEAD: no derived amplitude.** The 500 Mpc scale is dimensionally inconsistent (units m² not m, S632), a defect the S674 census initially mis-marked "unverified" and S675 corrected. ⚠ The two ledgers still disagree on the *label*: this header files TEST-07 under "genuinely untested," `PREDICTIONS.md` Bucket 2 files it under **refuted**. Refuted is the stronger and better-sourced reading — a dimensionally inconsistent scale is not awaiting data.

These five span quantum (consciousness), astrophysical (binaries, BAO), and cosmological (GW, interference) scales. If all five come back negative, Synchronism's predictive power is seriously in question. If even one comes back positive, it's worth a decade of follow-up.

---

## What "Negative" Would Mean

A sweep of negative results across all Tier 1 and Tier 2 tests would mean:
- C(ρ) is a useful *fitting function* but not a *generative theory*
- The reparametrization pattern (γ = 2/√N_corr unifying across domains) is mathematical convenience, not physics
- The consciousness coherence threshold is not universal
- Scale inversion symmetry is not real

This would NOT invalidate:
- The 97.4% dark matter phenomenology (that's validated regardless)
- The chemistry track's classification utility
- The insight that coherence is a useful lens for cross-domain analysis
- The methodological innovations (A2ACW, autonomous research, arc system)

## What "Positive" Would Mean

Even a single confirmed unique prediction would:

> **⚠ The four originally named here — TEST-04, TEST-07, TEST-14, TEST-17 — are all disqualified by this
> document's own header (04/14 executed → collapsed; 07/17 no derived amplitude, 17 excluded by Lorentz
> tests at ~11 orders of magnitude). Names removed 2026-08-19 rather than reranked: there is no surviving
> candidate in this catalog to promote in their place.**

- Establish Synchronism as a predictive framework, not just descriptive
- Open massive new research programs at those scales
- Justify the claim that C(ρ) encodes real physics across scales
- Make the theoretical tests (T1–T3) urgent priorities

---

*"All models are wrong; some are useful. The experiments above determine which kind this is."*
