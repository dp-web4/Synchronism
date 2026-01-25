# PRIMARY ROUND 3 — CONSENSUS OUTPUT

**Sources**: claude1-3.md, claude2-3.md, nova-3.md
**Arbiter**: claude2
**Scope**: Session #291 only

---

## 1. Why arcsine appears (derivation sketch)

Assume state is sinusoidal: s(t) = A·sin(ωt + φ)

**Velocity approach:**
- Velocity: ds/dt = Aω·cos(ωt)
- Time spent near position s ∝ 1/|ds/dt|

**Change of variables (equivalent):**
- Sample phase θ = ωt + φ uniformly on [0, 2π)
- s = A·sin(θ), so ds/dθ = A·cos(θ)
- |dθ/ds| = 1/(A|cos(θ)|) = 1/√(A² - s²)

**Result:**
Two solutions per s (except endpoints) gives factor of 2/(2π) = 1/π:

**P(s) = 1/(π√(A² - s²))  for |s| < A**

This is the arcsine distribution. Diverges at ±A (extremes), minimum at s = 0. Pure kinematics—no QM invoked.

---

## 2. Binary emergence

**Mechanism: Finite resolution binning (thresholding).**

The arcsine distribution is continuous. Binary outcomes require coarse-graining:
- Detector with resolution Δs bins outcomes into "near +A" vs "near −A"
- Because P(s) diverges at ±A, these bins capture dominant probability mass
- Without thresholding, weak measurement shows full continuous arcsine shape

**Yes, this requires binning/thresholding.** Binary emergence is a detector property, not intrinsic to the distribution.

Any additional "phase-locking selects extremes" story is interpretive—dynamics not provided. **Not specified in #291 — Bridge Missing.**

---

## 3. Three failure regimes ("breaks if...")

1. **Sampling is not uniform in phase/time**: If measurement times are phase-correlated or engineered (not uniformly distributed), the arcsine transform doesn't apply. Different sampling measure → different distribution.

2. **Measurement record does not monotonically map to s**: If the readout r(t) doesn't reflect s such that histogram over trials corresponds to distribution of s, the prediction fails. **Not specified in #291 — Bridge Missing.**

3. **Resolution too high or noise distorts**: If detector resolution is fine enough (or transfer function/noise distorts), "near ±A" bins don't dominate—outcomes remain broadly continuous or peak elsewhere.

*(Supplementary from claude2: Anharmonic oscillation or ξt < 4 also break the model, but these are input assumptions rather than experimental failure modes.)*

---

## 4. Does #291 reproduce QM measurement probabilities?

### **SUBSET ONLY.**

**#291 reproduces:**
- Existence of binary outcomes (two dominant regions)
- Why extremes are favored over midpoints
- Symmetric ~50/50 statistics when initial phase is unknown

**#291 does NOT reproduce:**
- Arbitrary Born rule probabilities |⟨ψ|basis⟩|²
- How superposition coefficients (α, β) map to outcome probabilities
- Basis-dependent statistics (X vs Z measurement)
- POVM structure / measurement operator formalism

The arcsine distribution is symmetric: P(near +A) ≈ P(near −A). There is no mechanism in #291 to generate asymmetric probabilities (e.g., P(|0⟩) = 0.7, P(|1⟩) = 0.3) from a prepared state.

The mapping from quantum state preparation → oscillation parameters (A, ω, φ) that would recover arbitrary probabilities is **Not specified in #291 — Bridge Missing.**

---

## Summary Table

| Question | Consensus Answer |
|----------|------------------|
| Arcsine origin | Uniform-phase sampling of sinusoid (standard math) |
| Binary emergence | Requires finite-resolution binning (detector property) |
| Failure regimes | Non-uniform sampling, bad readout mapping, fine resolution |
| QM reproduction | **SUBSET ONLY** — symmetric binary, not Born rule |

---

**STOP. Ready for CHALLENGER.**
