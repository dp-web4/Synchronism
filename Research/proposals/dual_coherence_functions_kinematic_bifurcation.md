# Proposal: Dual Coherence Functions — Kinematic Layer Bifurcation

**Filed**: 2026-05-07  
**Origin**: Site visitor feedback (4-persona pass), Pass 4 Researcher  
**Proposal type**: Research gap — structural

---

## The Problem

The Synchronism site currently uses two distinct coherence functions:

1. **C(ρ)** = tanh(γ · ln(ρ/ρ_crit + 1)) — density-based, used for galaxy dynamics, chemistry, phase transitions
2. **C = f(γ, D, S)** — D = decoherence rate, S = self-modeling capacity, used for consciousness and measurement

The landing page claims "one equation described reality from quantum to cosmic." By the site's own admission (on `/coherence-function`), this is false: there are at least two coherence functions, operating on different state variables, with no stated reduction from one to the other.

This is not a site-framing problem — it is a structural research gap.

---

## Three Cases

### Case A: C(γ, D, S) reduces to C(ρ) in some limit

If D and S are functions of ρ (e.g., decoherence rate ∝ environmental density, self-modeling capacity ∝ neural/computational density), then C(γ, D, S) → C(ρ) asymptotically. This would preserve the "one equation" claim and require showing:

- D(ρ) = explicit function
- S(ρ) = explicit function  
- C(γ, D(ρ), S(ρ)) ≡ C(ρ) in the relevant limits

**Status**: Not derived. No reduction written down in the archive.

### Case B: The framework has two distinct coherence functions with different domains

C(ρ) governs physical systems (below the level of self-modeling). C(γ, D, S) governs systems with cognitive architecture. They are genuinely different functions that happen to share the coherence variable. This would require:

- A regime map: when does each apply?
- A boundary condition: what triggers the switch from C(ρ) to C(γ, D, S)?
- A consistency check: in the "both-applicable" regime (e.g., simple neural systems), do they agree?

**Status**: No regime map exists. The site acknowledges the ambiguity without resolving it.

### Case C: C(γ, D, S) is the general form; C(ρ) is the specialization

If D = γ · ln(ρ/ρ_crit + 1) and S = 0 (for non-self-modeling systems), then C(ρ) could be a special case of the general form. This would require:

- Identifying what D represents in the physical (non-cognitive) domain
- Showing that S → 0 for inert matter in a physically meaningful sense  
- Writing the full C(γ, D, S) in explicit functional form (not just variables)

**Status**: Not derived. C(γ, D, S) is currently only stated as "a potential ambiguity problem requiring future research" (site's own language).

---

## Why This Matters

The "one equation" framing is the framework's central organizing claim. If there are two coherence functions:

1. The claim is false — and should be updated on the site
2. The framework's actual contribution is a *family* of coherence descriptions parameterized by system complexity, which is a different (and potentially more interesting) claim
3. The kinematic layer question (what does coherence operate *on*?) becomes acute: C(ρ) operates on density; C(γ, D, S) operates on what, exactly? A probability distribution over states? A configuration space?

---

## Proposed Research Task

For each of Cases A, B, C:
1. State whether the case is structurally possible given the current archive
2. If A or C: write down the explicit functional reduction
3. If B: write down the regime map and boundary conditions

If none of A/B/C can be written down coherently, the honest conclusion is that the dual-C problem is unresolved and the "one equation" framing should be withdrawn from the landing page.

---

## Connection to Kinematic Layer

The dual-C problem is one face of the kinematic-layer gap (MEMORY: "Born rule, dual-C, N_corr are three faces of one missing layer"). The framework has dynamics (how C evolves or applies) without a state space (what C operates on). Until the state space is written down, both C(ρ) and C(γ, D, S) are phenomenological fits, not equations of motion.

---

## Suggested Session Label

`Session_N: Dual Coherence Functions — Reduction or Bifurcation?`

Try: write down Case A explicitly. If it fails (S(ρ) cannot be derived from first principles), accept Case B and write the regime map. The result either closes the kinematic-layer gap or precisely locates it.
