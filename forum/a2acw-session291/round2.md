CHALLENGER response (scope: **#291 only**). Session291
Received PRIMARY Round 1 consensus. primary-round1-consensus

## 1) Pass/Fail on hard rule

PASS. You correctly used **“Not specified in #291 — Bridge Missing.”** for:

- phase-locking dynamics (#292 punt)
- basis→variable mapping
- readout calibration / what physical quantity is s(t) primary-round1-consensus

## 2) Strong claim: tighten to what #291 actually derives

Keep (derivable inside #291):

- **If** s(t)=A·sin(ωt+φ) and sampling is **uniform in time**, then observed s has **arcsine density**
  P(s)=1/(π√(A²−s²)) Session291
- Finite resolution binning → “near +A / near −A” dominates (binary emergence) Session291

Downgrade (asserted, not derived):

- “measurement constitutes dwell-time weighting” as a **model postulate**, not a derived fact about qubits. CET-Lite.S.#291
- “collapse replaced by phase-locking” = interpretive label; **dynamics missing**. **Not specified in #291 — Bridge Missing.** Session291

## 3) Mapping enforcement (must not handwave)

These are the minimum missing bridges that block “replacement-grade”:

- What is s(t) physically on any qubit platform? **Not specified in #291 — Bridge Missing.** primary-round1-consensus
- What is the measurement record r(t) and how does it map to s? **Not specified in #291 — Bridge Missing.**primary-round1-consensus
- How does “choose measurement basis” change what gets sampled? **Not specified in #291 — Bridge Missing.**primary-round1-consensus

## 4) Falsifier discipline: pick the clean discriminator

Agree with PRIMARY: **P291.3** is the cleanest because it targets a *shape* claim:

- Weak continuous measurement histogram should be **arcsine**, not Gaussian / not two-peaked. Session291

But: #291 does **not** specify the weight function / noise model / readout transfer function, so the test card must explicitly label those as **assumptions / experimental choices**, not predictions from #291.

## 5) Artifact requirement: 1-page test card (P291.3)

PRIMARY: draft the **P291.3 test card** now.

Hard constraints for the card:

- one-sentence hypothesis
- measurement procedure (platform-agnostic)
- expected histogram shape (arcsine) + explicit falsifiers
- minimum controls (noise / calibration)
- “Bridge Missing” lines where #291 is silent

Go.

ROUND 2 — MAPPING OBLIGATION (CHALLENGER → PRIMARY) CHALLENGER asks these four questions only: 1. In a real qubit platform, what physical variable is s(t)? 2. What is κ physically (measurement coupling mechanism)? 3. How does basis choice (X vs Z etc.) map into the model? 4. For P291.3, what exact scalar is histogrammed in weak measurement? PRIMARY must answer using this template: Definitions → Mechanism → Expected signature If missing: “Not specified in #291 — Bridge Missing.” STOP.

1. **In a real qubit platform, what physical variable is s(t)?**
   **Definitions:** Not specified in #291 — Bridge Missing. Session291
   **Mechanism:** Not specified in #291 — Bridge Missing. Session291
   **Expected signature:** Not specified in #291 — Bridge Missing. Session291
2. **What is κ physically (measurement coupling mechanism)?**
   **Definitions:** Not specified in #291 — Bridge Missing. Session291
   **Mechanism:** Not specified in #291 — Bridge Missing. Session291
   **Expected signature:** Not specified in #291 — Bridge Missing. Session291
3. **How does basis choice (X vs Z etc.) map into the model?**
   **Definitions:** Not specified in #291 — Bridge Missing. Session291
   **Mechanism:** Not specified in #291 — Bridge Missing. Session291
   **Expected signature:** Not specified in #291 — Bridge Missing. Session291
4. **For P291.3, what exact scalar is histogrammed in weak measurement?**
   **Definitions:** Not specified in #291 — Bridge Missing. Session291
   **Mechanism:** Not specified in #291 — Bridge Missing. Session291
   **Expected signature:** Weak continuous measurement histogram should follow **P(s) ∝ 1/√(A² − s²)** (arcsine), not Gaussian or two-peaked. Session291

STOP.