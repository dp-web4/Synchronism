# Back-annotate from site: A2ACW Specificity Null Baseline — The Reparametrization Detector Has Untested Specificity

**Filed:** 2026-05-22 (explorer track, synchronism-site)
**Source finding:** `synchronism-site/explorer/findings/a2acw-detector-false-positive-rate-null-baseline.md`
**Predecessor:** `a2acw_vocabulary_asymmetry_followup.md`, `a2acw_v2_three_axis_protocol.md` (2026-05-19)

## The claim being audited

The vocabulary-asymmetry experiment (2026-05-19) reported that modern-register translation catches **4/4** prior-art-class A2ACW demotions (vs. 0/6 for temporal asymmetry), and the three-axis protocol catches **6/6**. These numbers are cited as evidence that A2ACW is a *reparametrization detector* that generalizes beyond Synchronism — the project's strongest surviving contribution after the physics closed.

## The flaw

**All six test cases were already known to be reparametrizations** (post-hoc selected as the demotions). 4/4 and 6/6 are **sensitivity (TPR)** figures on a **positive-only** set. **Specificity (TNR / 1−FPR) was never measured.** A detector characterized only on positives is uncharacterized — the identical null-class flaw the project diagnosed for chemistry r=0.982 (high correlation forced by target monotonicity; control never run).

## The control, run

Held-out benchmark of 6 canonical genuine discoveries the detector must NOT flag (Dirac 1928, Bell 1964, BCS 1957, Higgs 1964, Hawking 1974, Noether 1918), each chosen because its modern-register translation *triggers* prior art for its ingredients. Plus 3 out-of-distribution true reparametrizations (Eddington α⁻¹, tired light, generic MOND-assuming TFR "derivation") to confirm sensitivity off the Synchronism sample.

Two decision rules tested:
- **R1 (literal protocol rule):** flag if a canonical prior-art reference is named within one round of the translated claim.
- **R2 (steelmanned):** flag if the claim *reduces to* prior art with nothing novel added.

### Results (best-case knowledgeable adversary, upper bound)

| | Reparametrization (n=3) | Genuine discovery (n=6) |
|---|:---:|:---:|
| **R1 flags** | 3 (TP) | **6 (FP)** |
| R1 silent | 0 | 0 |

R1 specificity = **0%**. Every genuine discovery names canonical antecedents under translation, because all non-trivial physics has antecedents. R1 detects *has-a-canonical-name*, not *reparametrization*.

| | Reparametrization (n=3) | Genuine discovery (n=6) |
|---|:---:|:---:|
| **R2 flags** | 3 (TP) | 0 (FP) |
| R2 silent | 0 | 6 (TN) |

R2 specificity = 100% — but **all** discrimination is supplied by a *novelty judgment* ("does the claim reduce to, or extend, the retrieved prior art?") that the protocol never operationalizes.

## Diagnosis: retrieval ≠ discrimination

Vocabulary translation contributes **retrieval** (puts the right prior art on the table — real, and the source of the 0/6→4/6 gain). It does **not** contribute **discrimination** (deciding the claim is *only* that prior art). Discrimination is a novelty judgment that is:
1. genuinely hard (BCS-vs-Cooper, Higgs-vs-Anderson divided human experts for decades);
2. exactly the step A2ACW failed — the temporal-asymmetry audit verified both agents *had* the prior art and didn't surface the reduction; the 6 demotions came from human audit, not the AI loop;
3. not automated by the protocol, so "A2ACW generalizes as a detector" is unsupported.

## Restated defensible claim (use this in any preprint)

> Forcing modern-register, framework-neutral restatement before adversarial review measurably improves **prior-art retrieval** (4/4 vs 0/6 on the demoted set). It is a **retrieval-augmentation / debiasing step for AI-assisted literature review**, not a reparametrization detector. The reparametrization-vs-discovery discrimination is an unautomated expert novelty judgment, and on a held-out benchmark that judgment carries 100% of the discriminating power; the protocol's literal rule has 0% specificity.

## Actions

1. **Gate the preprint.** The honest paper is *not* "we built a reparametrization detector." It is: "modern-register pre-translation improves prior-art retrieval in adversarial-AI review; the discovery/reparametrization discrimination remains an unautomated expert judgment that adversarial-AI loops fail." Null result about the detector + positive result about the retrieval aid, with the confusion-matrix figure above.
2. **Three-axis protocol caveat.** The 6/6 is sensitivity only; no axis has measured specificity; the prior-art axis has specificity 0% under its literal rule. Add to `a2acw_v2_three_axis_protocol.md`.
3. **Fresh-adversary test still needed** (`a2acw-v2-fresh-adversary-validation`): a weaker adversary degrades the R2 novelty judgment, strengthening the conclusion that the judgment (not retrieval) is load-bearing and fragile.

## Status

Productive failure. Eliminates "the project owns a validated transferable reparametrization detector"; replaces it with a narrower, honest, still-publishable retrieval-augmentation result. Consistent with the project's honest-assessment brand; overclaiming the detector would not be.
