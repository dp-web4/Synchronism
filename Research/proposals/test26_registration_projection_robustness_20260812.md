# Proposal: TEST-26's registration adjudicates an object that no longer exists — restate it class-level and projection-robust before adoption

**Date**: 2026-08-12 · **Origin**: site maintainer track, from visitor Pass 3 + Pass 4 (2026-08-12) crossed with the 2026-08-11 covariant-completion finding
**Gates on**: dp (adoption already does; this proposal changes *what* would be adopted)

## The problem, in two halves

**1. The draft criterion targets "the sign-locked locus" — which the archive's own covariant work killed the next day.**
The TEST-26 draft (2026-08-10) registers: *refuted if DESI DR3 excludes the sign-locked locus
(sign(w₀+1) = sign(wₐ)) at >3σ*. The 2026-08-11 completion scan
(`synchronism-site/explorer/findings/covariant-00-component-sign-lock-dies-desi-nogo-hardens.md`)
showed the literal sign lock is a property of the *substitution*, not the framework: completion A
(Appendix D as written) has no dark-energy sector at all, and completion B (Brans-Dicke-pinned)
populates mixed-sign (w₀, wₐ) pairs the lock forbade — while still never reaching the DESI quadrant
(0/192 γ at four ω). If DR3 is adjudicated against the *lock* and the lock is already dead by
internal derivation, the registration would be scoring a strawman — the same criterion-object drift
that produced the TEST-04a σ₈/fσ₈ swap, this time visible *before* adoption instead of after.

**The fix**: the registered object should be the class statement, which survived hardening:
*a slaved-dark-energy model (ρ_DE = ρ_m·F(x)) reaches the DESI quadrant iff ρ_DE(x) has an interior
maximum, and no completion of C = tanh(γ·ln(1+x)) — algebraic, substituted, or Brans-Dicke at any
scanned ω — produces one.* Pre-specify which completions are inside the registered class (all three),
and that the kill fires if DR3 robustly requires the crossing **all in-class members structurally
cannot produce**.

**2. The draft adjudicates in CPL (w₀, wₐ) space, and the projection-artifact literature — the
framework's best available defense — is uncited.**
Visitor Pass 4 (researcher persona, 2026-08-12): projecting a monotone non-CPL w(z) onto (w₀, wₐ)
is known to bias the recovered point, and there is an active literature arguing DESI's crossing
preference is partly a CPL-parameterization artifact — e.g. Shlivko & Steinhardt 2024
(arXiv:2405.03933), Cortês & Liddle 2024 (arXiv:2404.08056), Wolf, García-García & Ferreira
2024–25 (w₀wₐ reconstruction-bias papers). *(arXiv IDs from persona memory — verify against the
papers before external use; the standing from-memory flag applies.)* A site this compulsive about
steelmanning its executioners has so far missed the one live debate that cuts in its favor. Two
consequences for the registration:

- **The honest adjudication is a fit of the actual one-parameter w(z; γ) family** (substituted and
  completion-B forms; both are closed-form) **to the DR3 chains / compressed likelihoods**, not a
  quadrant read-off in a parameterization the model doesn't live in. If chain-level fitting is out
  of scope, the registration must at minimum pre-commit a projection-robustness check (does the
  crossing preference survive non-CPL parameterizations, per the Shlivko–Steinhardt method?) before
  the kill can fire.
- **The citation must go in both directions** — the same literature that could rescue the framework
  from a premature kill could also let a real crossing be dismissed as artifact. Pre-registering the
  robustness criterion *now*, before DR3, is what makes the eventual verdict clean either way.

## Also folded in (smaller)

- The draft's "the γ = 1/2 branch IS ΛCDM identically, so this test can only kill or tie" is now
  **substitution-conditional** (completion B has no ΛCDM member — the Möbius/Λ degeneracy is a
  property of the substitution). The kill-or-tie framing should carry that conditionality
  explicitly, as PREDICTIONS.md (2026-08-12 Publisher block) already does.
- Registration should state the freeze point: criterion frozen before DR3 data (with covariance
  source and SNe compilation named), per the draft's own discipline.

## Why this is a research-direction item, not site friction

The registration draft *is* a research artifact — it defines what the program's one live falsifiable
position means. Adopting it as drafted would re-inscribe an object (the literal lock) the archive
has already retired, and would leave the adjudication exposed to a parameterization bias the field
is actively litigating. Both fixes cost a paragraph now and an argument later.

## Suggested next actions

1. dp: decide whether TEST-26 adoption waits for the restated (class-level, projection-robust) form.
2. Explorer: fit the substituted and completion-B w(z; γ) families to the public DESI DR2 chains
   (or the compressed BAO+CMB+SN likelihoods) — quantifies today what the CPL quadrant only
   gestures at, and calibrates the projection bias for exactly this family. Topic being seeded
   site-side (`explorer/topics/fit-the-gamma-family-to-desi-chains.md`).
3. Verify the three projection-literature citations at the paper level before any appears on the
   site or in a registration.
