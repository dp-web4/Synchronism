# Proposal: Resolving the GW170817 / TEST-15 Question — Case 1 with No Derived Amplitude

**Filed:** 2026-05-26 (back-annotate from site: explorer)
**Resolves:** `gw170817_coherence_field_coupling_constraint.md` (2026-05-03, open three-case question)
**Source finding:** `synchronism-site/explorer/findings/test15-gw-dm-amplitude-closure.md`

---

## What the open question was

The 2026-05-03 proposal asked which of three cases describes the framework's GW
behavior: (1) derivative coupling → already falsified like TeVeS; (2) potential-only
coupling → survives but needs a dynamical mechanism; (3) parameterization, no action →
GW170817 doesn't apply. It left the classification open.

## Resolution: the archive already chose Case 1, and it closes TEST-15

Reading the primary sources (Session 59 §3.1 + Prediction 1; Session 642):

1. **Session 59 is a Case-1 propagation claim.** It writes c_g/c = 1 + α·(1−⟨C⟩_LOS) —
   the GW speed depends on the line-of-sight coherence field. This is exactly Case 1.

2. **The amplitude α is not derived.** Session 59 §3.1 sets "α ~ 10⁻¹⁵ (from GW170817
   constraint)." The only free parameter is read off the very experiment TEST-15 claims
   to perform. There is no independent prediction and no derived floor.

3. **GW170817 bounds α at near-maximal coupling, not vacuously.** Session 59's own path
   estimate is C_avg ≈ 0.1 → ⟨1−C⟩ ≈ 0.9. By C(ρ) = tanh(γ·ln(ρ/ρ_crit+1)), low-density
   extragalactic sightlines are low-coherence (C→0), so ⟨1−C⟩ ≈ O(1). Hence
   α·0.9 ≲ 10⁻¹⁵ → α ≲ 1.1×10⁻¹⁵, in the strong-coupling regime. The "GW170817 was a
   high-coherence path, low-coherence sightlines unconstrained" evasion conflates the
   strong-field *source* with the low-density *propagation path* and inverts C(ρ).

4. **The natural scale is dead.** The DM mechanism uses (1−C) ~ O(1) to produce O(1)
   mass discrepancies. If α tracked that coherence physics, α ~ O(1) — excluded by
   GW170817 at ~15 orders of magnitude. Survival requires the coherence field to couple
   to GW propagation 15 OOM more weakly than to galactic dynamics: an unstated,
   underived decoupling.

5. **Internal contradiction.** Session 59 (Case 1, propagation formula, α ≲ 10⁻¹⁵)
   and Session 642 (Case 3, "not a field theory, GW170817 doesn't apply") are mutually
   exclusive. You cannot keep TEST-15 *and* claim Case 3 — TEST-15 is the Case-1 claim.

## Verdict

TEST-15 is **not** a discriminating test. Either α ~ O(1) (refuted by GW170817) or α is
a free parameter bounded above only by the data it would test (no floor →
unfalsifiable as a positive prediction; GR-equivalent in the allowed range). It is the
GW-sector analogue of TEST-07: a domain where GR is distinct in principle, but with no
derived amplitude there is nothing to falsify. **Demote to exploratory (TEST-07
standard).** Novel discriminating-test count: unchanged at 0.

The honest Case-3 position ("we make no GW propagation claim; c_GW = c is assumed, not
derived, because we have no action") is *available* and defensible — but it requires
**removing TEST-15 from the test catalog**, because Case 3 is precisely the statement
that TEST-15's propagation prediction does not exist.

## Recommended archive actions

- Mark Session 59's TEST-15 / Prediction 1 as **superseded**: α is read off GW170817,
  not derived; the propagation claim is Case 1 and GR-equivalent in the surviving range.
- Note the Session 59 ↔ Session 642 contradiction explicitly in whichever doc serves as
  the GW canonical reference (Session 642 is the more honest; Session 59's "novel test"
  framing is the optimistic one that propagated to the site's Tier-3 page).
- Same audit recommended for TEST-17 (cluster γ-gradient) and TEST-21 (BAO sub-peaks):
  check for any derived amplitude; expected outcome is demotion to exploratory.
