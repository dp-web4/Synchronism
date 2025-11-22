# 4. Complexity and Phase

## 4.1 Phase Modes

A useful class of solutions:

```math
\Psi_q(x, t, C) = A(x) e^{i(\omega t - kC)}
```

- A(x) — spatial pattern / envelope.
- ω — temporal frequency.
- k — "wave number" along the complexity dimension C.

## 4.2 Qubits as ΔC-Stable Modes

In this view:

- A **qubit** corresponds to a mode Ψ_q that:
  - remains stable within a band of C
  - can be phase-tuned (via ω, k) without breaking.

"Too weak" phase adjustments:

- k ≈ 0, negligible effect → no effective computation.

"Too strong" phase adjustments:

- |k| too large, pattern destroyed by diffusion in C (term with D_c).

## 4.3 Goldilocks Band

We can formalize:

```math
k_{\min} < |k| < k_{\max}
```

or, equivalently, a control parameter λ (phase tuning strength):

```math
\lambda_{\min} < \lambda < \lambda_{\max}
```

- Below λ_min → no meaningful gate.
- Above λ_max → decoherence.

Current quantum hardware can be interpreted as having:

- non-ideal material / environmental properties that make λ_{\max} quite small
- hence the perceived fragility of qubits.
