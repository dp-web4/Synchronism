[STANCE: Skeptical/Analytical] [ARCHITECTURE: compression/over-smoothing bias] [CONFIDENCE: High] [VERIFICATION: None] [CHALLENGE_COUNT: 0] A2ACW v0.1 round 4

## P291.3 — 1-page experimental test card (Weak measurement histogram)

**Hypothesis (1 sentence):**
 Under weak continuous measurement, the underlying qubit “state” treated as a sinusoid $s(t)=A\sin(\omega t+\phi)$ will yield an arcsine-distributed readout over many samples: $P(s)=\frac{1}{\pi\sqrt{A^2-s^2}}$ on $|s|, rather than Gaussian or bimodal. Session291

### What to measure

Histogram over many trials of the **weak continuous measurement readout scalar** (the recorded measurement value per sample or per integrated window). **Not specified in #291 — Bridge Missing.** round 4 Session291

### Procedure (platform-agnostic)

1. Prepare a qubit and configure **weak measurement** (low coupling strength; continuous readout). Session291
2. Acquire a large dataset of readout values (per sample or per fixed integration window) across many trials. Session291
3. Normalize/calibrate the readout to a bounded interval $[-A,+A]$. **Not specified in #291 — Bridge Missing.** round 4
4. Build the empirical histogram $H(r)$ of the readout scalar $r$. **Not specified in #291 — Bridge Missing.** round 4
5. Fit/compare $H$ against the arcsine form $P(s)\propto 1/\sqrt{A^2-s^2}$ over $|s|. Session291

### Expected signature

Histogram shows **edge-enhanced density** consistent with arcsine: higher mass near $\pm A$, minimum near 0, bounded support $|s|. Session291 Session291

### Falsifiers (explicit)

Fail #291’s P291.3 prediction if, after normalization and reasonable smoothing:

- Histogram is **Gaussian-like**, or
- Histogram is **bimodal with interior maxima** inconsistent with arcsine edge divergence, or
- Support is unbounded / no meaningful $|s| interval can be established. round 4 Session291

### Minimum controls / confounds

- **Rail clipping / saturation control:** verify the weak-measurement chain is not saturating/clipping at rails (would fake edge pileup). round 4
- **Linearity / dynamic range calibration:** confirm readout linearity over the span mapped to $[-A,+A]$. **Not specified in #291 — Bridge Missing.** round 4
- **Noise-shape control:** repeat with known stationary/no-oscillation control readout (or randomized trigger timing) to ensure the acquisition chain doesn’t itself create edge pileup. **Not specified in #291 — Bridge Missing.** round 4

### Bridge-missing declarations (must be filled externally; not in #291)

- Exact definition of the **readout scalar** $r$ to histogram. **Not specified in #291 — Bridge Missing.** round 4
- Mapping from $r$ to the model variable $s$ and determination of $A$. **Not specified in #291 — Bridge Missing.** round 4
- Noise/transfer model sufficient to justify arcsine detectability under weak measurement. **Not specified in #291 — Bridge Missing.** Session291