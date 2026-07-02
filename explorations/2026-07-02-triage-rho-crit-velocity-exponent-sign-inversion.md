# Triage — ρ_crit(V) is sign-inverted vs MOND-matching (V⁻² required, V⁺² asserted): sharpens the door-#1 no-go (2026-07-02)

**Status:** `[ACTIVE-MRH]` — triage of `Research/proposals/rho_crit_velocity_exponent_sign_inverted.md` (explorer, 2026-07-02), tripped the hold-checklist. In my Phase-11 / door-#1 (C(ρ)/MOND) territory. **Verdict: ACCEPT — the MOND-matching V⁻² requirement is derived correctly (I verified it independently, two ways), the magnitude checks out, and it is a genuine *sharpening* of the door-#1 locality no-go: a local-ρ(r)-keyed MOND mimic must have ρ_crit ∝ V⁻² (forced by BTFR + a₀, profile-independent), and C(ρ) asserts V⁺² — the wrong *sign*, which cannot be re-tuned away. Bucket 0 unchanged (0); this reinforces the already-refuted C(ρ)⇒MOND galaxy sector.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Independent verification

The load-bearing claim is the MOND-matching requirement ρ_crit ∝ V⁻². Checked both forms:
- BTFR `M = V⁴/(Ga₀)`; transition radius from `g_bar(r_t)=GM/r_t²=a₀` ⇒ `r_t = √(GM/a₀) = V²/a₀ ∝ V²`.
- mean density `ρ ~ M/r_t³ = (V⁴/Ga₀)/(V²/a₀)³ = a₀²/(GV²) ∝ V⁻²`; isothermal `ρ = V²/(4πGr_t²) = a₀²/(4πGV²) ∝ V⁻²`. Both **−2**. ✓
- Structural: `r_t ∝ V²` and `M ∝ V⁴` ⇒ `ρ ∝ V⁻²`, **independent of baryon profile**. ✓
- Magnitude: `a₀²/(GV²)` at V=150 km/s (a₀=1.2e-10) ≈ **0.14 M⊙/pc³** (isothermal ≈0.011) — in the proposal's 0.01–0.3 range for V=50–300; the framework's `0.029·V² ≈ 650 M⊙/pc³` is ~240–300,000× too high, gap growing with V. ✓ So with ρ_crit at inner-galaxy densities, the luminous disk sits below threshold and the C(ρ) knee is never crossed — the density-crossing is not what produces the rotation-curve transition.

The derivation and the magnitude are correct.

## Why this is stronger than the existing magnitude no-go

Prior door-#1 statements were **magnitude** ones (Phase-11: flat curve α=2 is outside the capacity family; S689: Milgrom-locality instance; ~1.7 dex cross-system ρ↔g_bar offset). This is a **sign** one, and sign errors are not fixable by re-tuning a coefficient:

> A modification keyed on local volumetric density must have a per-object threshold falling as **V⁻²**
> to track an a₀ acceleration threshold — forced by BTFR + the a₀ threshold, independent of the
> modification's functional form. C(ρ) asserts **V⁺²**: the wrong sign.

**Nuance (triage value) — the identification is forced, so the sign-inversion is not escapable.** One
could object "maybe ρ_crit tracks something other than the transition-radius baryonic density." But the
*mechanism's own claim* is that crossing ρ_crit **is** the rotation-curve (MOND) transition; that
identifies ρ_crit with the baryonic density where g_bar=a₀, which forces V⁻². To escape the V⁻²
requirement you must abandon "C(ρ) crossing produces the transition" — i.e. abandon the mechanism. So
the sign-inversion is structural, not a fixable normalization.

## Disposition

- **Sharpens, doesn't newly-refute:** the C(ρ)⇒MOND galaxy sector is already Bucket 2 (γ=2 refuted,
  ΔBIC=+184, S661). This adds an *independent, profile-free, sign-level* reason the C(ρ) galaxy
  mechanism fails, complementing Phase-11's family-exponent bound. Together: the capacity-family
  *profile* can't reach α=2 (Phase-11) **and** the ρ_crit *threshold V-scaling* has the wrong sign
  (this) — two different structural failures of the same door-#1 sector.
- **Transferable null:** it constrains *any* ρ(r)-keyed MOND mimic (per-object threshold must fall as
  V⁻²), so it belongs with the locality-triage null (preprint-strategy null 1) — a sharper, sign-level
  companion to the ~1.7 dex magnitude offset. Worth including if dp greenlights that preprint.
- **Ledger:** add the galaxy-sector self-audit line the proposal recommends (ρ_crit ∝ V² is not
  derived, and is anti-derived by MOND-matching which requires V⁻²). Done in PREDICTIONS.md.
- **Bucket 0 unchanged (0); arc AT REST.**

## So what

A real signal lifted the hold; the derivation is right and it tightens door #1 on a new axis. Door #1
was "MOND cage" by the profile-exponent bound (Phase-11); it is now *also* caged by the ρ_crit(V) sign
— the framework's local-density threshold scales the wrong way to track a₀, and no coefficient re-tune
fixes a sign. The genuinely-decisive door-#1 negative is now two-pronged (profile + threshold-sign),
both profile-independent. Recorded; arc returns to AT REST.
