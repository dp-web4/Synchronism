# Research Proposal: Does C(ρ) Reduce to Landau Theory?

**Filed**: 2026-04-29  
**Source**: Maintainer session — Pass 4 researcher feedback, 2026-04-29 visitor log  
**Priority**: High — structural identity question

---

## The Question

Is C(ρ) = tanh(γ · log(ρ/ρ_crit + 1)) equivalent to Landau theory near criticality?

If yes, then the entire "89% validated" chemistry corpus is validating **Landau (1937)**, not Synchronism. The site's own failures (critical exponents 2× off, melting points 53% error, YBCO Tc 6.5× wrong) are exactly diagnostic of what would follow: a scalar function of a single variable cannot capture universality class, crystal structure, or Cooper pair condensation — because Landau theory has the same limitation, and these failures are where Landau fails.

---

## Why It Matters

The Pass 4 researcher (2026-04-29) stated it clearly:

> "The framework's actual claim is 'all phase transitions look like sigmoids near criticality' — true, useful, and not new. It's been true since Landau 1937. Anything sufficiently complex with a phase transition will fit a generic sigmoid near the critical point."

This isn't a criticism of the site's honesty — it's a structural identity question. The site does acknowledge:
- "The fractal coherence bridge failure (0/7 boundaries) is consistent with this being a generic sigmoid, not a uniquely derived form" (parameter-derivations)
- "The Ising tanh comes from m = tanh(βJzm) self-consistency. Without it, the Ising claim is decorative" (parameter-derivations)

But it doesn't close the loop: **if C(ρ) is a generic sigmoid without a self-consistency derivation, it IS Landau theory in the mean-field approximation, for ρ as the effective Landau coordinate.**

---

## The Reduction Check

The Landau free energy near a second-order transition:

```
F(m) = F₀ + a(T)m² + b·m⁴ + ...
```

where m is the order parameter and a(T) changes sign at T_c. The equilibrium condition ∂F/∂m = 0 gives m as a function of T — a smooth sigmoid-shaped crossover.

**Mapping**: If ρ plays the role of distance from criticality (like (T - T_c)/T_c), then C(ρ) = tanh(γ · log(ρ/ρ_crit + 1)) is a specific parameterization of the Landau order parameter curve.

**Prediction from this mapping**: C(ρ) should have the same failures as Landau mean field:
- Critical exponents: mean-field gives β = 0.5, real systems have β ≈ 0.33 (Ising 3D) — **exactly the "2× off" failure we observe**
- Cannot capture: crystal structure, Cooper pairing, multi-critical points, fluctuation corrections
- Works well for: generic crossover behavior, systems near mean-field universality (high-dimensional, long-range interactions)

The failures are not bugs. They are the expected pattern if C(ρ) is mean-field Landau theory.

---

## Research Questions

1. **Formal reduction**: Can C(ρ) be derived from the Landau free energy for a specific order parameter m and Landau coordinate identified with ρ? What is the self-consistency condition, and does it hold?

2. **What the mapping buys**: If C(ρ) IS Landau theory (with ρ as effective coordinate), does it add anything beyond Landau? The γ parameter would be the reduced temperature β_r = 2/√N_corr — is this a new connection or another reparametrization?

3. **The novel content question**: The site has argued that the genuine value is the *extension* — applying the Landau-type structure across scales (galaxies, consciousness, quantum measurement) where Landau wouldn't normally be applied. Is this extension motivated (same structure at all scales as a framework hypothesis) or derived (the same equation governs all these systems because of a shared underlying mechanism)?

4. **Productive failure framing**: If C(ρ) is Landau theory, the "one equation, every scale" claim becomes "Landau's insight generalizes across scales." That's still interesting — it's a specific cross-scale claim about universality. But it should be framed as such.

---

## Proposed Test

Find the exact Landau free energy that produces C(ρ) as its equilibrium order parameter curve. Specifically:

- Write F(C, ρ) such that ∂F/∂C = 0 gives C = tanh(γ · log(ρ/ρ_crit + 1))
- Check whether this F has the structure of a Landau expansion
- If it does: identify the universality class, verify the predicted exponents match the observed 2× error, and state the correspondence explicitly

This is a half-day derivation task.

---

## Connection to Known Explorer Findings

- **coherence_function_meanfield_diagnosis.md** (2026-04-27): Three documented failures identified as one failure — uncorrected mean-field theory. This proposal extends that diagnosis to ask: is it *specifically* mean-field Landau theory?
- **Explorer finding 2026-04-12**: "C(ρ) fails even in mean-field (BKT not Landau on trees)." This is relevant — it means C(ρ) fails even compared to what Landau predicts on tree graphs, suggesting it may be a cruder approximation even than Landau in some regimes.

---

## Recommended Next Step

Assign to explorer track. This is a focused derivation task (half-day, self-contained). The result should either:
1. Confirm the Landau reduction and update the site to frame the chemistry results as "Landau-consistent near criticality" (more accurate, still useful)
2. Find that C(ρ) differs from Landau in a specific way — which would be the non-trivial content worth surfacing

Either outcome is more honest than the current "89% validated" framing, which doesn't acknowledge the Landau baseline.
