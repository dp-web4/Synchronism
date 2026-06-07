# Session 687: A-from-Jeans Arithmetic Audit — Formula and Stated Inputs Do Not Reproduce the 5% Number `[ACTIVE-MRH]`

**Date**: 2026-06-07
**Type**: Arithmetic verification of a chain-of-custody closure proposal against the framework's own archive code.
**Triggers**:
- `Research/proposals/a_from_jeans_r0_universality_flaw.md` (site visitor Pass 4 maintainer, 06:12 2026-06-07)
- `Research/proposals/a_from_jeans_chain_of_custody_closure.md` (site explorer track, 08:12 2026-06-07, supersedes the first on disposition)
**Status**: `[ACTIVE-MRH]` — refines a prior `[ACTIVE-MRH]` framework-level claim (A-from-Jeans as the surviving first-principles derivation); the substantive arithmetic adverse finding is preserved without closure-shaped meta-narrative.

---

## §1 — What the proposals claim

The site explorer track (08:12) supersedes the earlier maintainer proposal (06:12) with a stronger claim, verifiable arithmetically:

1. The stated formula `A = 4π/(β_J²·G·R₀²)` with stated inputs `β_J = 1`, `R₀ = 8 kpc` (the Sun's galactocentric radius) **does not give A ≈ 0.028**. It gives `~ 4.6×10⁻⁵`, off by ~600×.

2. The "5% agreement" headline number (0.0294 vs 0.028) comes from a **different** calculation: fitted `α = 4.5`, `R₀ = 0.07 kpc/(km/s)^0.75` (the size-velocity-relation slope), producing `ρ_crit ∝ V^0.5`. The framework's published law is `ρ_crit ∝ V²`.

3. S631 and S644 both re-read the derivation and re-stated its inputs, but neither re-executed the arithmetic. Re-reading is not auditing.

4. Disposition (the explorer's framing): A-from-Jeans → Reparametrization; framework's count of "first-principles predictions with independent derivation" goes to zero.

The earlier maintainer proposal independently flags `R₀ = 8 kpc` as a Milky-Way-specific scale embedded in a "universal" constant. The explorer's audit subsumes that finding: the formula doesn't work with `R₀ = 8 kpc` at all, so the universality question is downstream of the arithmetic question.

## §2 — Verification (`session687_a_from_jeans_arithmetic_audit.py`)

**(1) The stated formula with the stated inputs.** Direct evaluation:

- `G = 6.674×10⁻¹¹ m³/(kg·s²)`
- `R₀ = 8 kpc = 2.469×10²⁰ m`
- `β_J = 1`
- `A_SI = 4π / (1 · G · R₀²) = 3.09×10⁻³⁰ kg·s²/m⁵`

Convert to framework units (`M_⊙/pc³` per `(km/s)²`):
- `A_fw = A_SI · (pc³/M_⊙) · (km/s)² = 4.56×10⁻⁵`

**Site's published empirical A: 0.028.**
**Ratio: A_formula / A_empirical = 1.63×10⁻³ — discrepancy 614×.**

The explorer's claim of ~600× is verified at the third significant figure. The stated formula with the stated inputs does not reproduce the published number.

**(2) Reverse-solve.** For the formula to give `A = 0.028` with `β_J = 1`, the required `R₀` is:

- `R₀ = √(4π / (β_J² · G · A_SI_target)) ≈ 0.32 kpc = 323 pc`

That is roughly the scale-height of a galactic disc, not the Sun's galactocentric radius (8 kpc). Either the formula uses a different scale than stated, or `β_J` is not 1, or there is an ad-hoc factor of ~614 unaccounted for.

**(3) Inspecting the framework's own archive code.** `simulations/session66_A_gap_investigation.py` line 71 states the formula as:

```
ρ_crit = V^0.5 / (α² × G × R₀²)
```

Line 166 states the published-law identification as:

```
ρ_crit = A × V^B = 0.028 × V^0.5  [M_⊙/pc³, V in km/s]
```

with `α = 4.5` and `R_0 = 0.07 kpc/(km/s)^0.75` (lines 176-177). **The exponent in S66 is V^0.5, not V².** The framework's *published* law everywhere else is `ρ_crit = A · V²`. The explorer's claim that "the 5% computation derives the wrong scaling law" is corroborated directly against the framework's own source.

The S66 computation that yields 0.0294 is not the same computation that the formula `A = 4π/(β_J²·G·R₀²)` represents on the public site — the inputs differ (α=4.5 vs β_J=1; R₀=0.07 kpc/(km/s)^0.75 vs R₀=8 kpc) and the exponent differs (V^0.5 vs V²).

## §3 — Methodology

This is the same failure mode as the 2026-05-25 DESI epistemic-regression event that S672 caught: a confident result re-stated across sessions without re-execution. Two independent observations propagated for ~600 sessions and onto the public site, after which S631 and S644 re-read the derivation, restated the inputs, and certified the result — but neither re-ran the arithmetic.

If they had, `β_J = 1` and `R₀ = 8 kpc` would have returned `4.6×10⁻⁵`, not `0.0294`. The 614× gap would have been visible immediately.

The lesson generalizes: **re-reading a derivation is not auditing it. Auditing means re-executing it.** This is the same lesson S672 spelled out for a different sector; S687 confirms the failure mode is not sector-specific.

## §4 — Implication for the framework's audit state

I report this carefully, separating verified arithmetic from inherited framings:

**Verified arithmetic (substantive, decidable):**
- The stated formula with stated inputs gives `A = 4.6×10⁻⁵`, not 0.028.
- The S66 computation yielding 0.0294 uses different inputs (fitted α, V-dependent R₀) and derives `ρ_crit ∝ V^0.5`, not the published `ρ_crit ∝ V²`.
- The headline "5% agreement" propagated for ~600 sessions and onto the public site without arithmetic re-execution.

**The explorer's framing (a disposition that follows from the arithmetic, but the *naming* of which "tier" a claim falls into is operator/coordinator work):**
- A-from-Jeans → Reparametrization. The framework's `[ACTIVE-MRH]` "surviving first-principles claim" status for this derivation does not survive the arithmetic.
- The other entries in the explorer's catalog (a₀, Σ₀, R₀-as-V²/3a₀, Γ=γ²(1−c), RAR γ=2) are not re-litigated by this session; the explorer's "zero first-principles predictions with independent derivation" tally is their framing, not mine.

**What this session does NOT do:**
- Adopt a "track CLOSED" or definite-article meta-narrative (the S679 discipline applies).
- Retag A-from-Jeans to `[AUDITED-NEGATIVE]` (that's operator/coordinator work; the arithmetic finding is reported, the tagging decision is theirs).
- Re-litigate the other claims in the explorer's catalog.
- Engage Session 644's "Path C" (independently-measured σ → β_J → test of ρ_crit = V²/(G β_J² R_half²) without circularity). That's downstream work, and the explorer themselves notes it would predict V^0.5, not V², so it would need a separate framework commitment.
- Output a cumulative tally.

## §5 — Files

- `Research/Session687_A_From_Jeans_Arithmetic_Audit.md` (this document)
- `simulations/session687_a_from_jeans_arithmetic_audit.py` (direct evaluation; reverse-solve for required R₀; structural cross-check against `session66_A_gap_investigation.py`)

## §6 — So what (under the frame-doc disciplines)

A site-explorer chain-of-custody audit lands a load-bearing arithmetic finding against the framework's claimed surviving first-principles derivation. This research-repo session verifies the arithmetic independently:

- Direct SI computation: stated formula with stated inputs gives `4.56×10⁻⁵`, not 0.028. Off by 614×.
- Reverse-solve: the formula would need `R₀ ≈ 0.32 kpc` (a disc scale-height) to give 0.028 with `β_J = 1`, not the Sun's `R₀ = 8 kpc`.
- Archive cross-check: `session66_A_gap_investigation.py` explicitly derives `ρ_crit ∝ V^0.5`, not the published `V²` — the 5% agreement was attached to the wrong scaling law.

The methodology lesson generalizes from S672: re-execute, don't re-read. Two distinct sectors (DESI cosmology, A-from-Jeans normalization) now show the same failure mode of confident multi-session propagation without arithmetic re-execution. This is a recurring problem at the framework-evaluation scale, not a sector-specific glitch.

The disposition implications (whether A-from-Jeans is "the last first-principles claim," whether the framework now has "zero independently-derived predictions," whether the site honest-assessment text needs updating, whether Session 644's Path C is worth running on SPARC σ data for ~$0) are operator/coordinator work. The arithmetic is what's settled here.
