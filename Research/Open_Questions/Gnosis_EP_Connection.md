# Gnosis: Self-Awareness Through Internal Coherence Detection

**Date**: January 11, 2026
**Status**: RESEARCH NOTE
**Origin**: Analysis of Gnosis repo (arXiv:2512.20578)

---

## What Gnosis Is

A lightweight (~5M param) "self-awareness head" that:
- Attaches to a frozen LLM backbone
- Reads hidden states + attention patterns + token confidence
- Predicts **correctness probability** before the answer is complete
- Achieves peak accuracy at ~40% generation completion

---

## Architecture (Three-Stream Introspection)

### 1. Hidden Circuit Encoder
```
last_hidden (B,S,D) → LayerNorm → Linear → d_tok
                    → Dilated Conv1D (d=1,2,4) with gated mixture
                    → SAB (self-attention) → PMA (pooling)
                    → D_HID features
```
**What it sees**: The trajectory of hidden states through generation - temporal patterns in the model's internal computation.

### 2. Attention Circuit Encoder
```
attention_maps (B,L,H,k,k) → Spectral stats (FFT, entropy, variance)
                          → CNN features on attention matrices
                          → Layer/Head embeddings + Set Transformer
                          → D_ATT features
```
**What it sees**: The structure of attention patterns - where the model is "looking" and how attention is distributed.

### 3. Confidence Circuit Encoder
```
token_probs (B,S) → [raw, derivative, logit] encoding
                  → Dilated Conv1D + SAB + PMA
                  → Rich stats (trend, volatility, drawdown)
                  → D_CONF features
```
**What it sees**: The evolution of model confidence over the generation - uncertainty patterns.

### Fusion: CorrectnessHeadLite
```
z_att, z_conf, z_hid → Learned gating (softmax weights)
                     → Concatenate weighted features
                     → MLP: D → 384 → 128 → 1
                     → Sigmoid → correctness probability
```

---

## Why This Matters for EP/SAGE

### This IS Coherence Backpropagation

From our Coherence_Backpropagation.md:
> "Memory is the functional deformation of the present caused by past instability."

Gnosis demonstrates exactly this:
- The hidden state trajectory contains signatures of impending failure
- These signatures are not explicit "I don't know" tokens
- They are **functional deformations** in the internal dynamics

### Error Signal Detection

| Gnosis Stream | Coherence Analog |
|---------------|------------------|
| Hidden trajectory | Phase coherence through time |
| Attention patterns | Spatial coherence in computation |
| Confidence evolution | Stability/instability signatures |

### Early Warning at 40%

The paper shows peak accuracy at ~40% completion. This means:
- Error signatures propagate through internal state before manifesting as wrong output
- The "decoherence" happens early and can be detected
- This is exactly "error signals bias future states toward stability" - but in reverse for detection

---

## Adaptation for SAGE

### What SAGE Could Learn From Gnosis

1. **Hidden State Monitoring**: Track hidden state trajectories during reasoning
2. **Attention Pattern Analysis**: Detect "scattered" vs "focused" attention
3. **Confidence Trajectory**: Watch for confidence drops during multi-step reasoning

### Key Differences

| Gnosis | SAGE EP |
|--------|---------|
| Correctness prediction | Epistemic state awareness |
| Binary (right/wrong) | Continuous calibration |
| Post-hoc analysis | Real-time monitoring |
| Frozen backbone | Integrated metacognition |

### Potential Integration

```
SAGE cogitation depth → Confidence trajectory encoder
SAGE attention over context → Attention circuit encoder
SAGE hidden state evolution → Hidden circuit encoder
                           ↓
                   Epistemic Calibration Score
```

---

## Key Findings from Code Analysis

### Feature Extraction
- **Spectral analysis**: FFT bands, spectral entropy on attention maps
- **Graph statistics**: Row/column variance, diagonal ratios, Laplacian trace
- **Temporal dynamics**: Total variation, slope, drawdown in confidence
- **Multi-scale processing**: Dilated convolutions at scales 1, 2, 4

### Training
- BCE loss on correctness labels
- Frozen backbone (only ~5M params trained)
- Transfer across model sizes (7B → 32B)

### Performance
- Outperforms 8B reward models for hallucination detection
- Generalizes across tasks (math, trivia, MMLU)
- Token-level probing validated internal signals exist

---

## Next Steps

1. **Experiment**: Try Gnosis on reasoning traces, not just answer correctness
2. **Adapt**: Build confidence trajectory encoder for SAGE cogitation
3. **Integrate**: Feed Gnosis-style features into SAGE's meta-learning
4. **Compare**: How does Gnosis score correlate with SAGE regret tracking?

---

## References

- [arXiv:2512.20578](https://arxiv.org/abs/2512.20578) - Original paper
- [dp-web4/Gnosis](https://github.com/dp-web4/Gnosis) - Fork of codebase
- [Coherence_Backpropagation.md](./Coherence_Backpropagation.md) - Related work on error-informed tick transitions
- SAGE/Web4 documentation - Target integration points

---

*"The error signal is already in the hidden state. Gnosis just learned to see it."*
