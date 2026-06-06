# Session 686: γ-Sign Flip Repairs S676 Ladder Under Reading A; Reading A vs B Is the Actual Fork `[ACTIVE-MRH]`

**Date**: 2026-06-06
**Type**: Verification of a site-back-annotated structural diagnosis (maintainer + explorer adjudication), connection to S676.
**Trigger**: `Research/proposals/gamma_ncorr_sign_inversion_sharpness.md` (site maintainer 2026-06-06 + explorer adjudication appended same day).
**Status**: `[ACTIVE-MRH]` — refines S676's "naming inversion" finding by surfacing the Reading A/B fork that S676 implicitly resolved into Reading A, substance preserved.

---

## §1 — What the proposal claims

The proposal identifies a sign inversion in `γ = 2/√N_corr` independent of the prefactor or scaling-law motivation:

- In statistical mechanics, `1/√N` is a *fluctuation width* — smaller width means a sharper transition.
- In `C(ρ) = tanh(γ·ln(ρ/ρ_crit+1))`, `γ` is in a *rate slot* — larger γ means a sharper transition.
- Placing `1/√N` (width-direction) into the γ slot (rate-direction) inverts the qualitative sign.

Result: the *most* correlated systems (BCS at `N_corr~10⁷`, BEC at `~10⁸`) get the *flattest* γ (`~6×10⁻⁴`) and therefore the most gradual C(ρ) — the opposite of physical reality, where these systems have sharp second-order transitions at a real `T_c`. The ideal gas (`N_corr=1`) gets the sharpest γ (`=2`) despite having no phase transition.

The explorer adjudication appended to the proposal makes a substantive addition:
1. Galaxy stars sit at `N_corr=1`, the swap-identity point of `2/√N ↔ 2√N`. Flipping the formula leaves every calibrated galaxy result unchanged.
2. The flip simultaneously repairs the sharpness inversion AND the coherence-magnitude inversion documented in S676 (`coherence_classicality_naming`). "One sign error, not two."
3. But neither sign has a fluctuation-width derivation. The deeper question is which *reading* of C the framework commits to:
   - **Reading A**: universal coherence scalar, presets comparable on one C axis. The visualizer and γ-calculator assume this.
   - **Reading B**: density-response function, γ as inverse effective T. The directional critique cannot even be posed.

## §2 — Verification (`session686_gamma_sign_flip_and_s676_ladder.py`)

**(1) Galaxy fixed point — trivial but worth stating.** At `N_corr=1`: `2/√1 = 2 = 2·√1`. The flip is the identity transformation at the galaxy regime. Every SPARC fit, `ρ_crit = A·V_flat²` calibration, and `γ=2` default carries through unchanged.

**(2) BCS shape under the flip.** At `N_BCS = 10⁷`:

| ρ/ρ_crit | C original (γ=6.3×10⁻⁴) | C flip (γ=6325) |
|---:|---:|---:|
| 10⁻⁶ | 0.000 | 0.006 |
| 10⁻⁵ | 0.000 | 0.063 |
| 10⁻⁴ | 0.000 | 0.560 |
| 10⁻³ | 0.000001 | 0.999994 |
| 1 | 0.0004 | 1.000 |
| 10 | 0.0015 | 1.000 |

Under the original formula BCS is pinned at C≈0 across the whole ρ range — the "BCS at C~0.0004" point S676 surfaced. Under the flip, C transitions sharply between `ρ/ρ_crit ~ 10⁻⁷` and `10⁻⁴`, saturating to 1 thereafter. Half-transition point: `ρ/ρ_crit ≈ 8.7×10⁻⁵`. The qualitative "sharp transition at a low ρ threshold then saturation" shape that BCS-with-a-real-`T_c` demands is recovered.

**(3) S676 ladder under the flip.** Evaluating C at `ρ=ρ_crit` for the S676 ladder:

| System | N_corr | C (original) | C (flip) |
|---|---:|---:|---:|
| Lone electron | 1 | 0.882 | 0.882 |
| Diatomic molecule | 2 | 0.753 | 0.961 |
| Small macromolecule | 100 | 0.138 | 1.000 |
| Mesoscale nanoparticle | 1000 | 0.044 | 1.000 |
| BCS superconductor | 10⁷ | 0.0004 | 1.000 |
| BEC | 10⁸ | 0.0001 | 1.000 |

Original column: monotone decreasing in `N_corr` (S676's "naming inversion"). Flip column: monotone increasing in `N_corr`. Under the flip, S676's "C anti-correlated with quantum coherence" becomes "C correlated with quantum coherence" — which is what fixing a sign error looks like — at zero cost to galaxy fits (the fixed point).

## §3 — The Reading A / Reading B fork

Both the proposal's sharpness claim and S676's magnitude claim require a single C axis across systems — Reading A. Under Reading B (C is a system-specific density-response, γ is inverse effective T), the ladder doesn't type-check at all, because comparing C(electron) to C(BCS) on the same scale presupposes the universal-scalar reading.

This is the actual research question, and the explorer's bottom line names it correctly. What the research-repo perspective adds: S676's "naming inversion" finding *implicitly* committed to Reading A. The verification I ran for S676 assumed I could ladder C across system identities. Under Reading B, that ladder isn't valid — and neither is this proposal's ladder, and neither is the visualizer that displays both. S676's finding therefore stands *conditional on Reading A*; it does not stand as a frame-independent fact.

So the foundational tension is genuinely structural, not parametric:

| Choice | Survives | Loses |
|---|---|---|
| Reading A + original γ=2/√N | Tooling (visualizer/calculator), `1/√N` motivation pointer | Two sign inversions (S676 magnitude, this proposal sharpness) on all preset systems |
| Reading A + flip γ=2√N | Tooling, both inversions repaired qualitatively, galaxy fixed point untouched | First-principles derivation (neither sign has one) |
| Reading B (either γ form) | Per-system formula directionality is defensible | Tooling (no cross-system C comparisons, presets invalid, visualizer collapses) |

There is **no choice that keeps both the framework's existing tooling and the original formula direction**. The framework is in a tools-vs-formula tension at the C *ontology* layer, not at the parameter layer.

The honest reading: the maintainer's site caveat (γ-calculator line 52, added 2026-06-06) flags the inversion as a problem, but the regime-description strings (lines 13/16) still assert the inverted C-magnitudes as fact. The tool's UI already contradicts its own caveat — which is the live-site manifestation of the deeper tools-vs-formula tension. Picking a reading lets the maintainer resolve the contradiction; not picking one keeps it visible.

## §4 — Implication for S676

S676 concluded: "C is anti-correlated with quantum coherence AND with synchronization." That conclusion stands *under Reading A*, with the qualifier that the original `γ=2/√N_corr` formula is in effect. Under the flip (still Reading A), the same C is *correlated* with quantum coherence — the inversion is repaired. Under Reading B, neither correlated nor anti-correlated is statable.

S676 therefore needs the same Reading-A qualifier the proposal's claim does. The substantive finding — *the framework's namesake variable behaves in a way that is regime-specific and ontology-dependent* — is preserved, but its sign is no longer load-bearing. Sign flips do not stabilize meaning. What does is committing to a reading. This refines S676 without overturning it: the inversion was real *under the reading my session assumed*, but my session did not surface the reading-choice.

## §5 — What this session does not output

- **No commitment to Reading A or B** — that is operator/coordinator work; my role is to make the fork explicit and verify the load-bearing computations.
- **No first-principles derivation of γ(N_corr)** — neither sign has one; the proposal's #1 conclusion ("the choice is purely about which qualitative inconsistency to carry") is not refuted by anything I did.
- **No site-text edits** — the proposal's "Proposed Site Action" and the explorer's "P1 live-site contradiction" surfacing are operator-channel work.
- **No retag of S676** — the substantive finding is preserved; this session adds the Reading-A qualifier. Tagging stays `[AUDITED-NEGATIVE]` on the old R(I) substrate; what shifts is the framing precision.
- **No cumulative tally.** Per the S679 discipline.

## §6 — Files

- `Research/Session686_Gamma_Sign_Flip_S676_Ladder_And_Reading_A_vs_B.md` (this document)
- `simulations/session686_gamma_sign_flip_and_s676_ladder.py` (galaxy fixed-point check; BCS shape under flip; S676 ladder under original vs flipped formula)

## §7 — So what (under the frame-doc disciplines)

A maintainer-filed structural diagnosis pairs with the explorer's deeper observation that the sign question is undecidable until the C ontology is committed. This session verifies the load-bearing numerical claims (galaxy fixed point; BCS shape recovery under flip; S676 ladder inversion under flip) and adds one research-repo refinement: S676's "naming inversion" finding implicitly committed to Reading A and stands *conditional* on that reading, with its sign reparable by a derivation-free formula flip that has the right qualitative effect on every preset and is free at the galaxy regime.

The actual research question — "commit C to Reading A or B" — is a frame question of the kind the autonomous prompt explicitly asks for, surfaced by the convergence of this proposal and S676. Both readings have a cost: Reading A keeps the tooling but requires a derivation-free sign flip (or carries the inversion); Reading B saves the formula's direction but invalidates every cross-system preset the framework's UI displays. The framework cannot honor both its tooling and its formula direction simultaneously.

The proposal's operator-level recommendations (site-text fix, derivation pursuit, tooling change) are not autonomous-session output.
