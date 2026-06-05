# Session 685: TEST-04a S₈-Receding Strengthens Prior Disfavor + TEST-02 Wide-Binary Reduces to MOND+EFE Unconditionally `[ACTIVE-MRH]`

**Date**: 2026-06-05
**Type**: Verification of two site-back-annotated proposals' research questions, integration with prior archive work.
**Triggers**:
- `Research/proposals/test04a_s8_receding_baseline.md` (site maintainer, 2026-06-05)
- `Research/proposals/test02_triple_conditional_status.md` (site maintainer, 2026-06-05)
**Status**: `[ACTIVE-MRH]` — refines prior `[AUDITED-NEGATIVE]` findings (S672 for TEST-04a; S684 for the EFE/boost-ceiling sector now extended to wide-binary regime), substance preserved, scope-of-closure made more precise.

---

## §1 — What the proposals claim

Both proposals are site-side back-annotations — already applied to `synchronism-site/honest-assessment/` — and ask one research question each.

**Proposal 1 (TEST-04a S₈ receding)**: σ₈≈0.76 was calibrated to the S₈ lensing tension (KiDS/DES era, ~2-3σ from Planck). That tension itself is receding per DES Y3 6×2pt and KiDS-Legacy 2024-2025 reanalyses. The Synchronism σ₈≈0.76 is therefore "post-hoc against a moving baseline." Research Q: if S₈ → Planck (σ₈≈0.83), what is the status of the prediction?

**Proposal 2 (TEST-02 triple-conditional)**: TEST-02 (wide-binary density dependence) is triple-conditional, not a standing discriminator: (1) the anomaly is methodologically disputed; (2) MOND+EFE prediction not differentiated; (3) ~80× below Gaia DR3 reach. Research Q (focus on condition 2): does the Synchronism amplitude for σ_int(ρ_env) differ from MOND+EFE?

## §2 — Verification (`session685_sigma8_calibration_and_efe_degeneracy.py`)

### (1) Proposal 1: was σ₈≈0.76 derivable from C(ρ_cosmo) at galaxy-anchored ρ_crit?

Direct evaluation of the framework's coherence function `C(ρ) = tanh(γ·ln(ρ/ρ_crit+1))` at cosmological mean matter density:

- `ρ_crit (galaxy-anchored)` = 10⁻²³ kg/m³ (S678 lower edge)
- `ρ_cosmo (matter)` = Ω_m·ρ_crit_cosmo ≈ 3.0×10⁻²⁷ kg/m³
- ratio: `ρ_cosmo/ρ_crit ≈ 3×10⁻⁴`

At γ=2 (galaxy default): `C(ρ_cosmo) = tanh(2·ln(1 + 3×10⁻⁴)) ≈ 6×10⁻⁴`

At γ=2/√N_corr with cosmological N_corr~10⁶⁰: `C(ρ_cosmo) ≈ 6×10⁻³⁴`

Either way `C(ρ_cosmo) ≈ 0`. **There is no first-principles path from C(ρ) to a structure-growth modulation of order Δσ₈/σ₈ ≈ 0.07/0.83 ≈ 8% at the galaxy-anchored ρ_crit.** The σ₈≈0.76 "prediction" was a calibration to the S₈ tension, not a derivation — consistent with S668 (mechanism-class failure) and S672 (post-hoc disfavoring, kill triggered).

**Research question outcome (both options negative):**
- Option (a) — calibration anchor disappears: holds, because no first-principles σ₈ derivation exists to fall back on
- Option (b) — discrepancy with DESI grows: holds, because the calibrated σ₈≈0.76 vs DESI 0.841 was already 2.4σ disfavored, and S₈→Planck pushes the target up

Neither option offers a survival path. The S₈-receding observation strengthens S672's disfavor verdict by removing the calibration anchor itself.

### (2) Proposal 2: does the Synchronism wide-binary amplitude differ from MOND+EFE?

For 1 M_⊙ + 1 M_⊙ at typical wide-binary separations under McGaugh-Lelli-Schombert `ν_e(y) = 1/(1−exp(−√y))` and S684's boost-ceiling `B_max`:

| sep (AU) | y = g_int/a₀ | ν_MOND | σ/σ_MOND at B_max=3.17 | at B_max=5 | at B_max=∞ |
|---:|---:|---:|---:|---:|---:|
| 5000 | 1.977 | 1.32 | 1.000 | 1.000 | 1.000 |
| 10000 | 0.494 | 1.98 | 1.000 | 1.000 | 1.000 |
| 15000 | 0.220 | 2.67 | 1.000 | 1.000 | 1.000 |
| 20000 | 0.124 | 3.37 | 0.969 | 1.000 | 1.000 |
| 30000 | 0.055 | 4.79 | 0.814 | 1.000 | 1.000 |

The MOND boost at typical wide-binary separations (5000-15000 AU) is 1.3-2.7, **below the S684 fit-XOR-discriminate fork's `B_max=3.17` cap**. The boost-ceiling doesn't bite in the regime where the wide-binary anomaly was claimed.

**This is a STRONGER closure than I first expected.** S684's framing — "the fork's distinct branch (B_max~3.17) is RAR-refuted; the fork's fitting branch (B_max→∞) is MOND" — implicitly assumed the operating regime sits *above* the cap so both branches are reachable. For wide binaries, the operating regime sits *below* the cap on every branch, including the RAR-refuted distinct one. **The Synchronism wide-binary prediction reduces to MOND+EFE unconditionally across the S684 fork's model space**, not branch-dependent.

Condition (2) of TEST-02 fails not because of a parameter choice but because the wide-binary regime is *below where the framework's distinguishing knob even operates*. The triple-conditional collapses: condition (2) fails for the whole S684 fork model space, independent of conditions (1) (anomaly real?) and (3) (Gaia DR4 sensitivity) resolving.

## §3 — Integration with the existing pattern

The two proposals connect to prior `[AUDITED-NEGATIVE]` findings in compatible ways:

| Sector | Prior finding | This session adds |
|---|---|---|
| Cosmology (TEST-04 / TEST-04a) | S668 mechanism-class failure; S672 post-hoc disfavoring (2σ kill triggered) | The calibration anchor itself (S₈ tension) is receding, so the post-hoc indictment is "against a moving target" — strengthens the disfavor verdict by removing the anchor's stability |
| EFE/TDG/wide-binary (TEST-02) | S684 boost-ceiling fork: B_max=3.17 distinct but RAR-refuted; B_max→∞ ≡ MOND | The wide-binary regime sits BELOW B_max=3.17 on every branch, so the EFE-degeneracy is unconditional in that regime, not branch-dependent — generalizes S684 from "fork-determined" to "regime-determined" for the wide-binary case |

**A small refinement to S684's "structural fork across 3 sub-sectors" framing**: the cluster (S678/S683), galaxy RAR (S661), and EFE/TDG (S684) sectors all sit in the regime where the boost-ceiling cap is reachable, so the fork's two branches matter. The wide-binary sector sits BELOW the cap on every branch, so the fork doesn't divide into reachable choices — the result is the same on all branches. The fit-XOR-discriminate structure is *visible* in the 3 sub-sectors where the cap reaches the regime; in regimes below the cap, the framework's prediction *is* MOND directly, no choice involved.

The site maintainer's recommendation in Proposal 2 — retire "possibly TEST-02" as the standing novel prediction — is supported by this cross-sector argument independent of how conditions (1) and (3) resolve.

## §4 — What this session does not output

- **No verdict on broader sector closure**. The S668/S672 verdicts on TEST-04a stand as `[AUDITED-NEGATIVE]` on the old R(I) substrate; the active substrate reformulation (saturation reframe with independent **J**, Phase-1 sim) is upstream. The S684/this-session verdicts on TEST-02 stand on the same substrate, same caveat.
- **No first-principles σ₈ re-derivation**. Proposal 1 lists this as an open option (c); it would be operator/coordinator work and the verification above suggests it's not achievable with the current C(ρ) anchored at galaxy scale.
- **No MOND+EFE wide-binary amplitude computation in the same units as Synchronism**. Proposal 2 asks for this; the S684 cross-sector argument settles condition (2) without it because the wide-binary regime sits below where the fork's two branches diverge.
- **No engagement with conditions (1) (anomaly real?) or (3) (Gaia DR4 sensitivity)** of TEST-02. Those are observational and don't require framework-level adjudication.
- **No cumulative tally.** Per the S679 discipline.

## §5 — Files

- `Research/Session685_TEST04a_S8_Receding_and_TEST02_EFE_Degeneracy.md` (this document)
- `simulations/session685_sigma8_calibration_and_efe_degeneracy.py` (C(ρ_cosmo) evaluation at galaxy-anchored ρ_crit; wide-binary boost-ceiling scan across the anomaly's claimed separation range)

## §6 — So what (under the frame-doc disciplines)

Two site-side back-annotations resolved by cross-referencing prior `[AUDITED-NEGATIVE]` archive work plus one analytical check each. Both proposals' research questions reduce to negative conclusions already implicit in the archive, sharpened by this session:

1. **TEST-04a**: the S₈-receding observation makes "post-hoc against a moving target" precise; the σ₈≈0.76 prediction has no first-principles fallback at the galaxy-anchored ρ_crit (C(ρ_cosmo)~6×10⁻⁴), so both research-question options are negative.
2. **TEST-02 condition (2)**: the wide-binary regime sits below the S684 boost-ceiling cap on every branch of the fork; the framework's wide-binary prediction reduces to MOND+EFE unconditionally in that regime.

A small framing refinement to S684: the fit-XOR-discriminate fork's two branches are *reachable* only in regimes where the boost cap bites; in regimes below the cap (wide binaries) the prediction *is* MOND on every branch, no fork to choose between. The "structural fork across 3 sub-sectors" finding applies to regimes where the cap is reached; a 4th sub-sector (wide-binary) shows what happens below the cap — the framework predicts MOND directly.

The proposal-level operator/coordinator actions (re-derive σ₈ from first principles; compute MOND+EFE wide-binary amplitude in matching units; site-honest-assessment updates already applied) are not autonomous-session output.
