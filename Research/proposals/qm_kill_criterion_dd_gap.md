# QM Kill Criterion Does Not Engage the Dynamical Decoupling Literature

**Filed:** 2026-05-08 (back-annotate from site: maintainer session)
**Status:** Open question — specify at QIP level or retire

---

## The Claim on the Site

Key Claims page, Claim #1 (Quantum Mechanics Is Synchronization Physics) advertises a kill criterion:

> "Design a noise environment where resync outperforms isolation, but standard decoherence theory
> predicts it doesn't. If isolation wins uniformly, the synchronization ontology adds nothing."

This sounds specific but is not operationalized at any level a QIP experimentalist could implement.

---

## The Problem: DD Already Does This

The dynamical decoupling (DD) literature — Viola-Knill-Lloyd 1999, CPMG, UDD (Uhrig 2007),
CDD (Khodjasteh-Lidar 2005), and recent transmon DD experiments — demonstrates exactly the
predicted phenomenology: periodic pulse sequences ("resync" in loose language) beat passive
isolation in non-Markovian baths by multiple orders of magnitude.

If "resync outperforms isolation" is the prediction, the DD literature already confirms it —
under standard decoherence theory, not synchronization ontology. The prediction is trivially true
in 1/f noise environments and already exploited in superconducting qubit platforms.

This means the kill criterion as written **cannot kill anything**: either
(a) the prediction reduces to a known DD result → reparametrization, not novelty, or
(b) it predicts something DD does NOT predict → needs precise specification.

---

## What Would Make This a Real Prediction

A falsifiable version of the kill criterion requires specifying:

1. **System**: what physical system (transmon, NV center, trapped ion, optical lattice)
2. **Bath spectral density**: Ohmic, 1/f, super-Ohmic, structured — the prediction must be
   bath-specific because DD's advantage is bath-dependent
3. **"Resync" protocol**: exact pulse sequence, interpulse spacing, number of pulses N — not
   loose "resync" language
4. **Predicted T2 ratio**: what ratio of T2(resync)/T2(isolation) does the synchronization
   ontology predict, vs what ratio does standard decoherence theory (Bloch-Redfield + DD)
   predict? If the ratios are numerically identical, the prediction is indistinguishable
5. **Kill threshold**: at what measured ratio would "synchronization ontology adds nothing" be
   concluded vs "synchronization ontology adds something"?

---

## Two Branches to Investigate

**Branch A: The prediction collapses to standard DD**

If Synchronism's "resync outperforms isolation" is mathematically equivalent to "periodic
pulse sequences outperform free evolution in 1/f baths," the claim is a vocabulary restatement
of Viola-Knill-Lloyd 1999. This is already the expected outcome given:
- The site admits (Sessions #266-270) that the Born rule derivation is a reparametrization
- The site admits all superconductor predictions reduce to standard condensed matter
- The pattern is consistent: Synchronism maps onto existing physics in different language

**Branch B: The prediction distinguishes from standard DD**

For the synchronization ontology to add something beyond DD, there must be a specific prediction
where DD's standard equations predict one outcome and Synchronism's equations predict another.
A candidate: standard DD coherence time scales as T2 ∝ N^α where α depends on bath spectral
density. Does Synchronism predict a different α? Or a different bath-spectral-density dependence?

This would require:
1. Deriving the decoherence timescale from MRH dynamics (not yet done in sessions)
2. Comparing to Bloch-Redfield + DD equations (standard QIP textbook)
3. Identifying whether any measurable difference exists

---

## Recommended Path

Do a one-session comparison: derive T2 from MRH dynamics for a transmon in a 1/f bath under
a CPMG-style pulse sequence. Compare to the standard Bloch-Redfield + UDD result. If the
equations are identical → label as "DD reparametrization" on the Key Claims page and retire
the kill criterion. If they differ → specify the protocol at QIP level and add to the test catalog
as a new Tier 2 experiment.

Without this, the advertised kill criterion is unfalsifiable by definition: "resync outperforms
isolation" is true in most interesting experimental regimes under standard decoherence theory,
so observing it would never distinguish between "standard QM" and "synchronization ontology."
