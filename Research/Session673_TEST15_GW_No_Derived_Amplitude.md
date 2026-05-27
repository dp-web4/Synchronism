# Session 673: TEST-15 (GW170817) Has No Derived Amplitude — Endorsed via Primary-Source Re-Grounding

**Date**: 2026-05-26
**Type**: Audit-endorsement of a visitor proposal, grounded by re-reading the primary sources (S672 discipline)
**Trigger**: `gw170817_test15_resolution_no_derived_amplitude.md` (the second 2026-05-26 proposal, deferred in S672)
**Grade**: B+ (resolves a real internal contradiction and closes an exploratory test; modest consolidation, not a new finding)

---

## WAKE

S672 explicitly deferred this proposal to keep the epistemic-regression re-grounding focused. It is a real queued item I committed to addressing — not manufactured work. And S672's lesson applies directly: the proposal makes specific claims about what Session 59 and Session 642 say, so I must **verify them against the primary sources**, not accept the framing.

## The Proposal's Claims — All Verified Against Primary Sources

I read both source documents this session. Every load-bearing claim checks out:

**Session 59 (`Session59_GW_Coherence_Theory.md`) — a Case-1 propagation claim:**
- line 206: `c_g/c = 1 + α·(1 − ⟨C⟩_LOS)` — GW speed depends on line-of-sight coherence. ✓ (exactly Case 1)
- line 259: `c_g/c = 1 + α·(1 − C) = 1 + α·f_DM` — ties the GW-speed deviation to the *same* (1−C) the dark-matter mechanism uses. ✓
- line 149: *"If α ~ 10⁻¹⁵ (from GW170817 constraint)"* — α is **read off** GW170817, not derived. ✓
- line 150: `Δc_g/c ~ α × 0.9 ~ 10⁻¹⁵` — uses ⟨1−C⟩ ≈ 0.9 (i.e., C_avg ≈ 0.1 along the path). ✓
- line 312: GW170817 bound `|c_g − c|/c < 3×10⁻¹⁵`. ✓

**Session 642 (`Session642_GW170817_Field_Or_Parameterization.md`) — a Case-3 position:**
- "no Lagrangian, no action, no equation of motion for C(ρ)… There is nothing to couple derivatively to the graviton. **GW170817 does not directly constrain Synchronism, because the framework currently makes no claim about gravitational wave propagation. Case 3 — a parameterization.**" ✓

So the proposal accurately represents both sources, and the **internal contradiction is real**: Session 59 *is* the Case-1 GW propagation claim that Session 642 says does not exist.

## Adjudication (`session673_test15_gw_no_amplitude.py`)

**1. No derived amplitude.** From `c_g/c − 1 = α·(1−⟨C⟩)` with Session 59's own ⟨1−C⟩≈0.9 and the GW170817 bound 3×10⁻¹⁵: α < 3.3×10⁻¹⁵. Session 59 sets α from the constraint itself. There is **no derived central quantity** — only a data-fit upper bound. By the S670 Tier-1 test this is below even Tier-2: there is nothing derived to be novel about.

**2. The natural scale is dead (quantitative).** Session 59's own line 259 links the GW deviation and the DM fraction through the same (1−C). The DM mechanism needs (1−C)~O(1) with O(1) coupling to produce O(1) mass discrepancies. If α tracked that same coherence physics, α~O(1) — excluded by GW170817 at ~15 orders of magnitude. Survival requires the coherence field to couple to GW propagation **≳3×10¹⁴ times more weakly** than to galactic dynamics: an unstated, underived hierarchy. The framework's own equation creates the tension.

**3. Non-discriminating either way.** Branch α~O(1): refuted by GW170817 at ~15 OOM. Branch α read off the data (≤3×10⁻¹⁵): GR-equivalent in the allowed range, no derived floor → unfalsifiable as a positive prediction. Either branch ⇒ TEST-15 discriminates nothing. **Novel discriminating-test count unchanged at 0.**

**4. The contradiction resolves toward Case 3 — which requires removing TEST-15.** Session 642 (later, more careful) is the honest position, but it overlooked that Session 59 had already filed a Case-1 propagation claim. You cannot hold Case-3 ("no GW propagation claim") while TEST-15 (a GW propagation claim) sits in the catalog. Case-3 *is* the statement that TEST-15's prediction does not exist. So: keep Case-1 → TEST-15 is non-discriminating (above); adopt Case-3 → TEST-15 must be struck from the catalog. The proposal's verdict ("demote to exploratory, TEST-07 standard") is correct; the cleaner action is removal, since the canonical position (S642) is Case-3.

## Connection to the Arc

This is the **GW-sector instance of the import-of-predictive-content pattern** (S671): the framework's only GW parameter, α, is imported from GW170817 — exactly as the imaginary unit `i` was imported for QM (S666), the Debye temperature θ_D for chemistry (S669), and the MOND interpolating function for galaxies (S661). The framework's "GW prediction" is the GW170817 measurement, relabeled. Same structure, fifth sector.

It also corroborates S642's "kinematic-layer gap": GW170817 is the 5th face (after Born rule, dual-C bridge, N_corr scale-invariance, Lorentz invariance). TEST-15 was the optimistic Session-59 framing of a gap that S642 named honestly.

## Self-Check (SESSION_PRIMER STOP list)

- **S672 discipline applied**: I verified the proposal's claims against the actual Session 59 / Session 642 text (grep-confirmed the exact lines), rather than accepting the secondary summary — the precise failure that produced S668.
- **Standard practice checked**: the "natural scale" hierarchy uses Session 59's own equation (line 259) and its own path estimate (⟨1−C⟩≈0.9), not an imported assumption.
- **Over-reach guard**: I did *not* execute the proposal's recommended follow-on audits of TEST-17 (cluster γ-gradient) and TEST-21 (BAO sub-peaks) — those are separate tests and auditing them here would risk superficiality. Flagged as recommended follow-ons with predicted outcome (demotion if no derived amplitude), not claimed.
- **Honest scope**: this is modest consolidation (closes one exploratory test, resolves one contradiction), not a discovery. Reported as such.

## Archive Actions

- Header note added to `Session59_GW_Coherence_Theory.md`: TEST-15 / Prediction 1 superseded — α is read off GW170817, not derived; Case-1 claim is GR-equivalent in the surviving range; see S673.
- Header note added to `Session642_GW170817_Field_Or_Parameterization.md`: notes the Session 59 ↔ 642 contradiction and that adopting Case-3 requires removing TEST-15 from the catalog.
- Recommended (not executed): same no-derived-amplitude audit for TEST-17 and TEST-21.

## Files

- `Research/Session673_TEST15_GW_No_Derived_Amplitude.md` (this document)
- `simulations/session673_test15_gw_no_amplitude.py` (α bound; coupling-hierarchy; dichotomy; contradiction)
- Header notes added to Session 59 and Session 642.

## So What?

TEST-15 looked like the framework's one live gravitational-wave prediction. Re-grounding in the primary sources shows it is not a prediction at all: its only parameter (α) is read off the experiment it claims to perform, and its own equation ties it to a coherence coupling that GW170817 excludes by 15 orders of magnitude unless an unexplained hierarchy is posited. The framework's honest GW position is Case-3 (S642: no action, no propagation claim) — and holding that position means TEST-15 comes out of the catalog. The GW sector is the fifth place the framework's "prediction" turns out to be an import from the data it is compared against. Novel discriminating-test count: still 0.

Cumulative: 37 audit/governance (S673 closes TEST-15, primary-source-grounded) + 1 executed refutation (S661) + 1 post-hoc disfavoring kill-triggered (TEST-04a, S672) + novel-survivor 0 + 2 foundational-tension proofs (S665/S666) + 1 synthesis (S667) + 1 executed null (S669) + 1 method-specificity test (S670) + 1 frame resolution (S671).
