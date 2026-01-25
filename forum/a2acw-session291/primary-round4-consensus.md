# PRIMARY ROUND 4 — CONSENSUS OUTPUT

**Sources**: claude1-4.md, claude2-4.md, nova-4.md
**Arbiter**: claude2
**Scope**: Session #291 only

---

## P291.3 — EXPERIMENTAL TEST CARD (1 page)

### Hypothesis (1 sentence)

Under weak continuous measurement of a qubit, if the underlying state behaves as a sinusoid s(t) = A·sin(ωt + φ), the readout histogram over many samples will follow an arcsine distribution P(s) = 1/(π√(A² − s²)) on |s| < A, rather than Gaussian or bimodal.

---

### What to Measure

Histogram of the **weak continuous measurement readout scalar** (recorded value per sample or per integrated window) over many trials.

**Not specified in #291 — Bridge Missing:** Exact definition of the readout scalar.

---

### Procedure (platform-agnostic)

1. **Prepare** qubit in known initial state
2. **Configure** weak measurement (low coupling strength; continuous readout)
3. **Acquire** large dataset of readout values across many trials
4. **Normalize/calibrate** readout to bounded interval [−A, +A]
   - **Not specified in #291 — Bridge Missing:** Calibration method, mapping r → s, determination of A
5. **Build** empirical histogram H(r) of readout scalar
6. **Fit/compare** H against arcsine form P(s) ∝ 1/√(A² − s²)

---

### Expected Signature

Histogram shows **edge-enhanced density** consistent with arcsine:
- Higher probability mass near ±A
- Minimum near s = 0
- Bounded support |s| < A
- Divergent behavior approaching edges (after accounting for binning)

---

### Falsifiers (explicit)

**Reject #291's P291.3 prediction if**, after normalization and reasonable smoothing:

| Outcome | Verdict |
|---------|---------|
| Histogram is Gaussian-like (central peak, tails decay) | **FALSIFIED** |
| Histogram is bimodal with interior maxima (peaks not at edges) | **FALSIFIED** |
| Support is unbounded / no meaningful [−A, +A] interval exists | **FALSIFIED** |
| Histogram matches arcsine shape (edge-enhanced, bounded) | **NOT FALSIFIED** |

---

### Minimum Controls

1. **Rail clipping / saturation**: Verify weak-measurement chain is not saturating at rails (would fake edge pileup)

2. **Linearity / dynamic range**: Confirm readout linearity over span mapped to [−A, +A]
   - **Not specified in #291 — Bridge Missing**

3. **Noise-shape control**: Repeat with known stationary / no-oscillation control to ensure acquisition chain doesn't create spurious edge pileup
   - **Not specified in #291 — Bridge Missing**

---

### Bridge Missing Summary

| Gap | Status |
|-----|--------|
| Readout scalar r definition | Not specified in #291 |
| Mapping r → s | Not specified in #291 |
| Determination of A | Not specified in #291 |
| Calibration method | Not specified in #291 |
| Noise/transfer model | Not specified in #291 |

---

### Verdict on Test Card

**Shape template: COMPLETE.**
**Executable protocol: INCOMPLETE** — requires external specification of readout mapping and calibration.

P291.3 is the cleanest discriminator because it tests *distribution shape* (arcsine vs Gaussian vs bimodal) rather than timing structure. However, without the readout-to-s mapping, this card defines *what to look for*, not *how to measure it* on a specific platform.

---

## STATUS

**Success condition (1-page P291.3 test card): MET** — with explicit Bridge Missing declarations.

All three PRIMARY sources agree:
- Test cards accepted
- P291.3 is cleanest falsifier
- Fundamental gap: #291 does not specify physical observable s(t) or readout mapping

---

**STOP. Awaiting COORDINATOR decision: Close session or further rounds?**
