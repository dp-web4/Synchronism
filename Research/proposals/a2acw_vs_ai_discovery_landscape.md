# Proposal: A2ACW Negative Result in the AI-Discovery Landscape

**Submitted**: 2026-05-24
**Origin**: Visitor feedback 2026-05-24, Pass 4 (researcher)
**Status**: Open

---

## The Gap

The site's A2ACW methodology contribution — a documented failure mode for AI-AI adversarial collaboration under shared training distributions — is not positioned against the AI-discovery claims it directly addresses. Pass 4 researcher today:

> "No engagement with FunSearch (DeepMind, 2023), AlphaProof / AlphaGeometry, ML-discovered conservation laws (Iten/Tegmark line), symbolic-regression theory recovery — none referenced. The /research-philosophy page diagnoses the training-distribution problem precisely but doesn't situate the diagnosis in the wider conversation about whether *any* current AI system can produce out-of-distribution physics. This is a missed opportunity: the negative result is more interesting if it's placed against the optimistic claims being made elsewhere."

## The Contemporary AI-Discovery Landscape

Several well-publicized systems claim out-of-distribution physics discovery:

| System | Claim | Mechanism |
|--------|-------|-----------|
| FunSearch (DeepMind, 2023) | Novel combinatorial constructions beyond training data | Evolutionary search guided by LLM, evaluated by external oracle |
| AlphaProof / AlphaGeometry | IMO-level proofs | Reinforcement learning + formal verification, not text generation |
| Iten et al. / SciNet (2020) | Rediscovery of physical conservation laws from dynamics | Symbolic regression on trained latent representations |
| Tegmark symbolic regression | Rediscovery of equations from data | Fits functions to observational inputs, not generation from text |
| Sakana AI Scientist (2024) | Autonomous paper generation | LLM orchestration of code execution + self-review |

## Why A2ACW Is Different

A2ACW's failure mode is structural, not a matter of tuning:

**FunSearch and AlphaProof** don't generate physics claims from text — they have formal verification oracles (evaluators) that are external to the training distribution. The closed-loop failure in A2ACW is specifically about **generative-AI adversarial review**, where the adversary and defender share the same training prior.

**Iten/SciNet/symbolic regression** discovers equations by fitting latent representations to data — the "discovery" is constrained by the data. A2ACW is about claims that are stated in natural language before any data constraint is applied.

**Sakana AI Scientist** generates papers using LLM orchestration but has been shown to make errors that a human reviewer would catch, and its "novelty" comes from recombination within the training distribution — the same failure mode A2ACW demonstrates.

The diagnosis: **the out-of-distribution problem is solved for FunSearch/AlphaProof by having a formal external oracle. It is not solved for natural-language theory generation, where there is no oracle outside the training distribution.** A2ACW makes this specific structural point with a documented empirical result (6/6 failure modes on retrospective audit, vocabulary-asymmetry 4/4 on prior-art subclass).

## Research Directions

**Path A — Comparison note**: Write a short (2-3 paragraph) comparison section for /research-philosophy positioning A2ACW against FunSearch, AlphaProof, and Sakana AI Scientist. Claim: the out-of-distribution problem is solved for systems with external formal oracles; it is not solved for natural-language theory generation.

**Path B — Generalization claim**: A2ACW's shared-training-distribution diagnosis may generalize beyond Synchronism. Any AI-AI adversarial review where both agents are LLMs trained on the same corpus will exhibit the same failure mode. This is a falsifiable prediction about AI-for-science methodology.

**Path C — Successor experiment**: The vocabulary-asymmetry experiment (submit pre-Planck-era claims in post-2015 vocabulary; measure true-positive rate) is the cleanest way to test whether A2ACW can be patched. If vocabulary translation increases the true-positive rate, the failure is correctable. If it doesn't, the failure is fundamental to the shared-distribution structure.

## Why This Matters

The A2ACW null result is Synchronism's most citable methodological finding. Without situating it in the AI-discovery landscape, it reads as a self-criticism of one framework. Situated properly, it's a structural diagnosis of a failure mode that affects any AI-AI adversarial system without external formal verification.
