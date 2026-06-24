# Phase-18 (generative axis, QC) — the compatibility lens against REAL device data: scope, test, verdict (2026-06-24)

**Status:** `[ACTIVE-MRH]` — generative-axis follow-through on Phase-17. dp's brief: scope the publicly
available QC data that could discriminate the compatibility-floor design rule; if it exists, test
whether a compatibility metric ⟨C⟩ beats standard metrics (connectivity, mean error rate); if it
doesn't, specify the minimal experiment. Honesty + anti-overclaim paramount: a NULL or DATA-LIMITED
result is a perfectly good outcome. **Generative/applied; Bucket 0 = 0 (unchanged).**
**Sim:** [`simulations/phase18_qc_compatibility_vs_symmetry_alignment.py`](../simulations/phase18_qc_compatibility_vs_symmetry_alignment.py)
· result: `simulations/results/phase18_qc_compatibility_vs_symmetry_alignment_result.json`
**Author:** CBP-Claude (Opus 4.8).

---

## TL;DR verdict

**(ii) testable-but-data-limited for the strong claim, trending (iii) subsumed-by-standard-metrics for
the scalar form.** The QC field has *already, independently* established the load-bearing half of the
hypothesis — that correlation **structure**, not magnitude, gates logical performance at fixed
marginal error rate — but it expresses it with a **sharper, more specific** metric than a generic
scalar ⟨C⟩: **alignment of the correlated errors with the code's stabilizer symmetry**. The lens
**names the right design axis** (real synthesis/reframing value) but does **not** beat the field's own
structural metric; a scalar ⟨C⟩ is a coarsening of it. The one genuinely open piece — testing the
metric on **real-device correlation structure vs measured logical performance across configurations**
— is reachable with public data (Google Willow syndrome data is on Zenodo) but **nobody has run it
under a compatibility/⟨C⟩ framing**, and the more specific symmetry-alignment metric would likely win.

---

## 1. Data scope — what is publicly available?

The hypothesis needs three layers of data. Honest availability assessment:

| Layer | Availability | Detail |
|---|---|---|
| **(a) Device coupling/topology graphs** | **Available** | IBM (Qiskit `backend.coupling_map` / backend properties), Google (published device layouts), IonQ/Quantinuum (all-to-all by architecture). Trivial to obtain. |
| **(b) Correlated-error / crosstalk structure** | **Partially published, and *extractable*** | IBM calibration data does **not** ship crosstalk/correlated-error terms by default. BUT the standard **p_ij correlation matrix** (probability that detection events i,j fire together) is extractable directly from QEC syndrome data — Google's d=5 surface code yields a 600×600 p_ij matrix (arXiv:2207.06431). Simultaneous/correlated RB, ZZ-crosstalk maps, and `LayerFidelity` (Qiskit) characterize crosstalk but are reported per-paper, not as a uniform public corpus. |
| **(c) The DISCRIMINATING piece — logical/QEC performance correlated with noise-correlation STRUCTURE across configurations** | **Scarce (the real gap)** | Logical error rates are published per-experiment, but a **cross-configuration study that varies correlation structure at fixed marginal rate and measures the logical-performance response on real hardware** is essentially absent. The fixed-marginal results that exist are **simulation**, not device. |

**Key public asset:** Google released the **Willow surface-code syndrome datasets on Zenodo**
([10.5281/zenodo.13273331](https://doi.org/10.5281/zenodo.13273331), and dynamic-codes
[10.5281/zenodo.14238907](https://doi.org/10.5281/zenodo.14238907)). From the raw syndrome data one can
compute the p_ij correlation matrix and reconstruct a **detector error model (DEM)** including
correlated/hyperedge terms — an active 2025–26 line (arXiv:2504.14643, 2310.12448, 2606.16288;
analogous IBM `ibm_miami` runs). So the *correlation structure* of a real device is recoverable. What
the public data does **not** give you is the **controlled cross-configuration sweep** of structure-at-
fixed-marginal needed to cleanly attribute logical performance to ⟨C⟩ vs standard metrics.

## 2. The decisive literature finding (convergence — and a sharper metric)

The QC field has independently reached the core of the compatibility hypothesis, in simulation:

- **arXiv:2506.15490** ("Symmetry in Multi-Qubit Correlated Noise Errors Enhances Surface Code
  Thresholds"): at **fixed marginal error rate**, correlated errors whose symmetry **aligns with the
  code's stabilizer structure** *raise* the threshold (~15–25% in their sims); **misaligned**
  correlations *lower* it (~5–10%). A **sharp** transition, not gradual. Structure, not magnitude.
- **arXiv:2410.23779** ("Detrimental non-Markovian errors for surface code memory"): also fixes the
  marginal rate; finds some correlation classes (data-qubit, pairwise) benign and others ("streaky"
  multi-time on syndrome qubits / two-qubit gates) ruin the logical scaling. Here the degradation is
  **gradual** in a correlation-strength parameter (no hard floor).

**This is exactly the compatibility-lens prediction** — "compatible (collective/aligned) coupling
helps, incompatible (destructive/frustrating) coupling hurts, *independent of magnitude*." That the QC
field arrived here from its own angle is **corroboration that the frame captures something real**
(dp's standing correction: convergence ≠ redundancy). **But** the field's metric is **more specific**:
not a scalar "fraction of compatible couplings" but *which* correlations are compatible — those whose
support **matches a stabilizer**. The lens's ⟨C⟩ is a coarse scalar over that structure.

(Note on the floor: Phase-17's *hard* compatibility floor came from frustrated **Kuramoto** —
appropriate for collective phase-sync. In the **QEC** setting the relevant object is logical error vs
correlation structure; 2506.15490 shows a sharp threshold shift, 2410.23779 shows gradual degradation.
So "sharp floor" is regime-dependent in QC, **not** universal — an honest refinement of Phase-17's
"hard floor" language when transferred to QEC specifically.)

## 3. The test I could run on real-device-*shaped* data (Phase-18 sim)

I did **not** have a clean cross-configuration real-device dataset to fit (layer (c) is the gap), so I
built the **discriminating mechanism test** instead: a minimal Monte-Carlo toy (not a decoder, not
real data) of a logical-parity observable, **holding the marginal single-qubit error rate and the
correlated-pair fraction FIXED**, and sweeping only whether correlated pairs are **stabilizer-aligned**
(both flips inside one stabilizer support → cancel in the logical parity → benign) or **misaligned**
(straddle the logical cut → fatal weight-2 logical flip).

**Results** (`p_marginal=0.02`, `f_corr=0.4`, 400k shots/condition):
- **Structure-only logical-error gap at identical marginal rate:** misaligned `0.1325` → aligned
  `0.1264` (gap `0.0061`). Real, in the predicted direction; modest because the toy is deliberately
  weak (the literature's threshold effect is far larger near threshold). Reproduces the
  *fixed-marginal, structure-matters* finding qualitatively.
- **The discriminator:** `corr(⟨C⟩, LER) = −0.961` and `corr(aligned_fraction, LER) = −0.961` —
  **identical**, because the scalar ⟨C⟩ is an **affine relabel** of the alignment fraction. The scalar
  "compatibility" predicts logical error *only because* it is a coarsening of the structural metric the
  field already uses. It adds **no independent predictive power** over symmetry-alignment.
- **Control:** sweeping the marginal rate (the *standard* metric) at fixed structure moves LER
  `0.071→0.133→0.232` — confirming marginal rate is its own, **orthogonal**, axis (so ⟨C⟩ is *not*
  redundant with mean error rate; it is redundant with the *structural* metric).

**What this shows honestly:** ⟨C⟩ is **orthogonal to the standard scalar metrics** (mean error rate,
connectivity) — that part of the lens is real and useful — but it is **subsumed by the field's own
structural metric** (stabilizer-symmetry alignment), which says *which* correlations count as
compatible. The lens does not beat it; it is a less-specific synonym for it.

## 4. If the strong claim is to be earned: the minimal real-device experiment

Because layer (c) is the gap, here is the concrete, falsifiable experiment a QC group could run (and
that public Willow syndrome data partly enables for the *analysis* half):

**Design (compatibility-floor / structure-gating test):**
1. Fix **code distance** and **mean physical error rate** (the two standard metrics) across arms.
2. Prepare **two correlated-noise structures** with the **same marginal single-qubit error rate** and
   the **same total correlated fraction f**, differing **only** in ⟨C⟩ ≡ stabilizer-alignment of the
   correlations: Arm-HI (correlations engineered to lie within stabilizer supports, e.g. via tailored
   dynamical decoupling / pulse scheduling that makes residual ZZ collective) vs Arm-LO (correlations
   straddling logical cuts, e.g. crosstalk left structured across the boundary).
3. Measure **logical error rate** vs code distance in both arms; extract the **p_ij correlation
   matrix** / DEM in each to *verify* the achieved ⟨C⟩ (this analysis is doable today on Willow data).
4. **Prediction (compatibility lens):** a **sharp logical-performance gap** between arms **not**
   explained by code distance or mean physical error rate — and, if a floor exists in this regime, an
   Arm-LO ⟨C⟩ below which adding code distance / lowering gate error **does not** restore exponential
   suppression.
5. **Discriminating analysis:** regress LER on {distance, mean error rate, connectivity, ⟨C⟩} and on
   {…, full stabilizer-alignment structure}. **The lens earns the strong claim ONLY if** scalar ⟨C⟩
   adds predictive power *beyond mean-error-rate/connectivity* (expected: yes — orthogonal axis) **AND
   is not fully captured by the finer alignment metric** (expected from Phase-18: no — it is
   subsumed). If ⟨C⟩ tracks LER but the alignment metric tracks it strictly better, verdict is (iii)
   for the scalar; the lens's value is the **reframing**, not a new winning metric.

This is the door between "generative vocabulary" and "demonstrated generative novelty." It is
reachable; the analysis half (extract structure from real syndrome data, correlate to logical
outcomes) can start now with the public Willow dataset, the controlled-arms half needs hardware time.

## 5. Honest verdict (i / ii / iii)

- **NOT (i)** "demonstrably predictive on real data": no one (including me here) has shown a
  compatibility metric beating standard metrics on real-device cross-configuration logical data; the
  controlled data doesn't exist publicly and I did not manufacture it.
- **(ii) testable-but-data-limited** for the design rule itself: the discriminating layer-(c) data is
  scarce; the minimal experiment is specified above and is falsifiable.
- **Trending (iii) "reduces to / is subsumed by standard metrics"** for the **scalar** ⟨C⟩: the QC
  field already has the *structural* result (correlation structure gates logical performance at fixed
  marginal rate) with a **more specific** metric (stabilizer-symmetry alignment). Phase-18 shows scalar
  ⟨C⟩ is an affine coarsening of that metric — orthogonal to mean-error-rate (good) but not a new
  winner over the field's structural metric.

**The defensible claim:** the compatibility lens **reframes usefully** — it names correlation-
structure / compatible-vs-destructive coupling as the gating design axis, unifies DFS +
correlated-noise-aware QEC + symmetry-tailored codes under one frame, and the QC field's independent
convergence on this axis is real corroboration. **The indefensible claim** (which I am NOT making) is
that a Synchronism ⟨C⟩ metric **demonstrably beats** standard QC metrics on real data: that requires
the experiment in §4, and the more specific stabilizer-alignment metric would likely win even then.
Synthesis value: yes. Demonstrated metric novelty: not shown, and Phase-18 suggests subsumption.

**Bucket 0 = 0, unchanged.** This is the applied/generative axis; no physics bet was added.

---

## Sources

- arXiv:2408.13687 / Nature s41586-024-08449-y — "Quantum error correction below the surface code
  threshold" (Google Willow; Λ=2.14, d=7 @0.143%/cycle; rare correlated/cosmic-ray events).
- Zenodo [10.5281/zenodo.13273331](https://doi.org/10.5281/zenodo.13273331) (Willow syndrome data),
  [10.5281/zenodo.14238907](https://doi.org/10.5281/zenodo.14238907) (dynamic surface codes).
- arXiv:2207.06431 — surface-code scaling; p_ij detection-event correlation matrix (600×600, d=5).
- arXiv:2506.15490 — symmetry-aligned multi-qubit correlated noise **enhances** thresholds (fixed
  marginal rate; sharp transition; simulation).
- arXiv:2410.23779 — detrimental non-Markovian / "streaky" correlated errors (fixed marginal rate;
  gradual degradation; simulation).
- arXiv:2504.14643, 2310.12448, 2606.16288 — estimating detector error models / correlation structure
  from real syndrome data (Willow, IBM `ibm_miami`).
- arXiv:2003.02354 — correlated randomized benchmarking; Qiskit `LayerFidelity` (simultaneous RB) for
  crosstalk characterization; IBM calibration data omits crosstalk by default.
