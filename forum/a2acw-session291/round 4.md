P291.1 — **Autocorrelation** Session291
**What to measure:** Autocorrelation function of rapid repeated measurement outcomes (binary record) vs lag τ. Session291
**Expected signature:** Peak(s) at τ ≈ 1/f_qubit (and/or harmonics), i.e., periodic structure at the underlying oscillation frequency. Session291
**Falsifies if:** No autocorrelation structure beyond shot-noise / statistical baseline (flat autocorrelation except τ=0). Session291
**One confound control:** Verify measurement apparatus / digitizer does not introduce periodicity: run identical acquisition on a known stationary/no-oscillation control (or randomized trigger timing) and confirm no peaks. *(Exact control implementation not specified in #291 — Bridge Missing.)* Session291

P291.2 — **Timing modulation** Session291
**What to measure:** Outcome probability P(|0⟩) (or mean of binary outcomes) as a function of measurement delay after initialization, stepped in τ_qubit/8 increments (or equivalent phase steps). Session291
**Expected signature:** Periodic modulation of outcome statistics with ~8-fold periodic pattern across one qubit period (phase-dependent probability variation). Session291
**Falsifies if:** Outcome probabilities are invariant (within error bars) across delay settings—no periodic modulation at the predicted periodicity. Session291
**One confound control:** Hold all non-timing parameters fixed and confirm delay sweep doesn’t change effective measurement strength / integration window. *(How to guarantee this is not specified in #291 — Bridge Missing.)*Session291

P291.3 — **Weak measurement histogram** Session291
**What to measure:** Histogram over many trials of the **weak continuous measurement readout scalar** (the recorded measurement value per sample or per integrated window). **Not specified in #291 — Bridge Missing.** Session291
**Expected signature:** Histogram follows arcsine form **P(s) ∝ 1/√(A² − s²)** over bounded support |s|<A (edge-enhanced density), rather than Gaussian or two-peaked. Session291
**Falsifies if:** Histogram is (a) Gaussian-like, or (b) bimodal with interior maxima inconsistent with arcsine edge divergence, after accounting for calibration/normalization. Session291
**One confound control:** Ensure weak-measurement chain isn’t saturating/clipping at the rails (which would fake edge pileup). Calibrate readout linearity and dynamic range. *(Calibration method and mapping readout→s not specified in #291 — Bridge Missing.)* primary-round3-consensus

STOP.