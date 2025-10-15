# Math Sketches: From GPT + Grok Tweaks

## GPT Originals
1. Phase-Ensemble ⇒ Born: Define ϕ ∈ T^n, μ(ϕ) invariant under F, coarse-grain to p_i = ||Π_i ψ||².

2. No-Signaling: Prove ∑_x P(x,y|a,b) = ∑_{x'} P(x',y|a,b') via μ push-forward factorization.

3. Tsirelson: E(a,b) = u_a · v_b, unit vectors ⇒ CHSH ≤ 2√2.

4. Monogamy: C_AB² + C_AC² ≤1 from shared-resource on phase-locks.

5. Lüders Update: F_meas marginalizes to standard ρ' = ∑ p_k E_k ρ E_k† / tr(...).

Distinctives: Sub-Q jitter σ_φ ∝ T^{-α}; Contextuality budget; Tick-hiding lemma.

## Grok Tweaks/Expansions
- For Monogamy: Model as shared nudge budget B ~ ∑ Δφ² ≤ τ (ticks). Lagrange: max C_AB + C_AC s.t. costs. Predicts δC ~10^{-3} in noisy tripartites—testable vs. 2025 W-state obs.

- μ Weighting: Haar uniform on T^n for qubits, weighted by Intent stability for generality (ties to 4.5).

- Space for Vectors: ℝ^d w/ d=MRH dim (4.9)—finite blankets cap super-Tsirelson.

- Sim Tip: Use QuTiP for phase evo, NumPy for lattice nudges.

[Full breakdowns from my previous response.]