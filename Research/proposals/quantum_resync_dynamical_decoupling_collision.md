# Proposal: Quantum "Resync" Claim Is a Vocabulary Collision with Dynamical Decoupling

**Date**: 2026-05-28  
**Source**: Visitor feedback (Pass 4 researcher, synchronism-site 2026-05-28) + key-claims spec-gap note  
**Status**: Research diagnosis — Reparametrization class  

---

## The Claim Being Diagnosed

From `/key-claims` (Quantum Arc, Sessions #228–237):

> "If decoherence is desynchronization, then periodic resync protocols should outperform continuous isolation for certain noise profiles. This is a direct engineering prediction that differs from the standard model's 'isolate harder' strategy."

This has been labeled "Untested" on the site, implying it is a novel open prediction.

---

## The Problem

**Periodic resync outperforming continuous isolation for certain noise profiles is a textbook result in dynamical decoupling (DD):**

- Viola & Lloyd (1998): first formal DD sequence — periodic π-pulse bang-bang control
- Carr & Purcell (1954), Meiboom & Gill (1958): CPMG sequences — periodic pulses beat passive isolation in 1/f and Lorentzian baths
- Uhrig (2007): optimal DD sequences for non-Markovian environments
- Khodjasteh & Lidar (2005–2009): concatenated DD sequences, scaling laws
- Biercuk et al. (2011), Du et al. (2009): experimental demonstrations — periodic control consistently outperforms isolation in specific spectral regimes

The DD field is ~25–30 years old with thousands of experimental papers. The claim that "periodic resync outperforms isolation for certain noise profiles" is not a Synchronism prediction — it is a standard statement of the DD framework, with the specific crossing between DD and passive isolation scaling with the bath spectral density.

**The key-claims page itself already acknowledges this** (spec-gap note): "Dynamical decoupling (DD) protocols — Viola-Knill-Lloyd 1999, UDD, CPMG — already demonstrate that periodic pulse sequences beat passive isolation in non-Markovian baths. If 'resync' reduces to DD, the prediction is known physics, not a novel test."

---

## Formal Diagnosis: Reparametrization, not Untested

The "resync" claim should be classified as **Reparametrization — Vocabulary Collision** because:

1. The physical content ("periodic control outperforms passive isolation in non-Markovian baths") is established DD fact, not a Synchronism prediction
2. The Synchronism framing ("desynchronization → resync") is a vocabulary layer on top of the DD conceptual structure
3. The kill criterion ("design a noise environment where synchronization model predicts resync wins but standard decoherence theory predicts it doesn't") cannot be met — standard decoherence theory (as DD) already *does* predict regimes where periodic control wins; the question is whether Synchronism predicts a *different* crossing or threshold

**What would be novel**: a specific bath spectral density, pulse sequence, or comparison metric where Synchronism's MRH-based model predicts a *different* T₂ ratio than the standard DD prediction. This has not been derived.

---

## Recommendation

1. **Site badge**: Change "Resynchronization outperforms isolation" from "Untested" to "Reparametrization — maps to dynamical decoupling; specify novel prediction before treating as open"

2. **Key-claims text**: Strengthen the spec-gap note to say explicitly: "The claim as stated is a vocabulary collision with the DD field. A novel prediction requires specifying the bath spectral density, pulse sequence, and predicted T₂ ratio where Synchronism *differs* from the DD prediction at the same pulse interval."

3. **Archive**: Update the quantum arc sessions (#228–237) to note the DD prior art and specify what a Synchronism-specific test would require. The question "does decoherence as desynchronization yield MRH-based DD sequences that differ from standard Uhrig/CPMG?" is a genuine open question — but it is currently unanswered.

---

## Why This Matters

The quantum arc currently fronts this as a "direct engineering prediction." It is not. This is the same category error as calling galaxy rotation curves a "Synchronism prediction" — the framework reframes known physics in its vocabulary but hasn't yet specified what differs quantitatively. Until it does, the engineering prediction is a relabeling, not a test.

---

## Connection to A2ACW methodology

This is a clean Axis-1 (vocabulary translation) catch: "periodic resync" → "dynamical decoupling" immediately locates the prior art. The site's own three-axis A2ACW protocol would have caught this under vocabulary audit — another example of the reparametrization pattern the site documents.
