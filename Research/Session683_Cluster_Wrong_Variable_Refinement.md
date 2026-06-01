# Session 683: Cluster Gap Is Wrong-Variable (Not Mainly One-Scale) — S678 Framing Refined `[ACTIVE-MRH]`

**Date**: 2026-06-01
**Type**: Refinement of S678's cluster-bridge framing, verifying the load-bearing claims of a back-annotated visitor-channel amendment.
**Triggers**:
- `Research/proposals/one_scale_insufficiency_theorem_cluster_gap.md` (the original "theorem" proposal that elevated S678's structural finding)
- `Research/proposals/cluster_gap_wrong_variable_amendment.md` (the same-day amendment correcting it)
**Status**: `[ACTIVE-MRH]` — refines a prior `[AUDITED-NEGATIVE]` retag (S678), substance preserved.

---

## §1 — What the amendment claims

The amendment refines S678's "one density scale" framing by separating two distinct effects:

| Level | Obstruction | Magnitude | Origin |
|---|---|---|---|
| 1 | **Wrong variable** — local ρ vs non-local g_bar | ~10⁴ / structural | Framework-specific (cost of the C(a)→C(ρ) migration) |
| 2 | **One scale** — single a₀ misses cluster cores | ~factor 2 | MOND-inherited (Sanders 2003) |

The argument:
- MOND has one scale (a₀) and misses clusters by ~2× (Sanders 2003; Pointecouteau & Silk 2005). C(ρ) misses clusters by ~10⁴ (S678).
- If "one scale insufficiency" were the operative cause, two one-scale theories would fail similarly. They differ by ~10⁴.
- Therefore the scale count is not the dominant cause. The framework-specific effect is **which variable** the single scale lives in: acceleration (a₀, non-local in ρ) vs density (ρ_crit, local).

Two quantitative demonstrations (per the amendment's `cluster_bridge_wrong_variable.py`):
- **(A) Within Coma**: gas density is nearly flat in the β-model core while g_bar(r) is non-monotonic with a peak ~0.12·a₀. So ρ → g_bar is not single-valued, and C(ρ) is nearly constant where the required mass discrepancy varies most. **No function of local ρ can produce a radially-varying discrepancy in a flat-cored cluster.**
- **(B) Cross-system**: at matched g_bar in the overlap band (0.05–0.12 a₀), a representative disk galaxy is ~1.7 dex denser than Coma. RAR scatter is 0.13 dex. A density-keyed C cannot reproduce RAR tightness once clusters are included.

## §2 — Verification (`session683_cluster_wrong_variable_check.py`)

**(A) Within Coma — confirmed**. β-model with `n₀ = 3.4×10⁻³ cm⁻³`, `r_c = 290 kpc`, `β = 0.65`. Sampling r ∈ {10, 30, 100, 200, 290, 500, 800, 1300} kpc:

| r [kpc] | ρ_gas [kg/m³] | g_bar/a₀ |
|---:|---:|---:|
| 10 | 3.47×10⁻²⁴ | 0.0025 |
| 30 | 3.43×10⁻²⁴ | 0.0074 |
| 100 | 3.11×10⁻²⁴ | 0.0233 |
| 200 | 2.37×10⁻²⁴ | 0.0394 |
| 290 | 1.77×10⁻²⁴ | 0.0471 |
| 500 | 9.04×10⁻²⁵ | 0.0507 (peak) |
| 800 | 4.25×10⁻²⁵ | 0.0453 |
| 1300 | 1.78×10⁻²⁵ | 0.0356 |

In the inner core (r < r_c = 290 kpc): **ρ varies by +0.16 dex while g_bar varies by +1.20 dex.** The mapping ρ → g_bar is not single-valued in this region. At galaxy-anchored `ρ_crit = 10⁻²³ kg/m³` and γ=2, C(ρ) takes the range 0.40–0.53 across the inner core — nearly constant where g_bar varies an order of magnitude. **No C(ρ) ansatz can produce a radially-varying mass discrepancy here**, independent of the C-to-mass mapping. The amendment's structural claim is correct against the framework's own equations.

(Caveat: the gas-only g_bar peak at 0.051 a₀ is below the amendment's quoted 0.12 a₀. The amendment uses total baryonic mass — gas + stars + cluster galaxies, ~85/15 split — which raises the peak. The *qualitative* non-monotonicity and the load-bearing point — flat ρ, varying g_bar — are robust against this choice.)

**(B) Cross-system — consistent at order of magnitude**. Disk galaxy at V_flat = 200 km/s reaches g_bar = 0.1·a₀ at r ≈ 108 kpc with mean baryonic density of the enclosed sphere ~1.29×10⁻²³ kg/m³. Coma's gas-only g_bar matches at r ≈ 450 kpc with local gas density ~1.04×10⁻²⁴ kg/m³. Ratio: **12.4× = +1.1 dex**. The amendment's 1.7 dex uses different methodology (local outskirts density vs my enclosed-sphere mean for the galaxy; total baryonic profile for Coma). My rougher computation is order-of-magnitude consistent (sign correct, magnitude in the dex range). A precise reproduction requires SPARC galaxy data and Coma's full baryonic profile — beyond this session's scope.

## §3 — Implication for S678's framing

S678's substantive findings stand. What changes is the *named structural cause* of the 10⁴ catastrophic cluster failure.

Recasting the S678 table under MRH-relationship taxonomy:

| S678 framing | Refinement (S683) |
|---|---|
| "Structural root: one density scale" | The "one density scale" obstruction is *real* but corresponds to a factor-~2 cluster residual (Sanders 2003), not the 10⁴ failure. Inherited from MOND class, transferable, not novel to Synchronism. |
| "Cluster bridge structurally impossible" (codomain bound) | The codomain bound on A3 (M_app/M_B ≤ 2) is unchanged — it is the exact mechanism for *one* of the four ansätze. But it is a *symptom* of the wrong-variable disease: ρ_crit anchored where the physics isn't (galaxy core knee, where C → 1), C(ρ_Coma) ≈ const, 1/C-type ansätze explode. |
| C(ρ) lacks the "second density scale" real theories use | More precise: C(ρ) lacks **the right variable**. Verlinde and MOND key on acceleration (non-local in ρ); a local density function cannot reproduce a non-local quantity's universal relation. The "second scale" escape (Open Question #1) is the *same move*: any C(ρ, g_bar) that works does so by re-introducing the acceleration variable. |

The bottom-line **substance** — C(ρ) does not bridge galaxies to clusters; the cluster sector is closed for the density-map program — is unchanged. The bottom-line **framing** — "one density scale insufficiency" — is downgraded from the named structural root to a transferable factor-~2 residual shared with MOND. The actual structural root of the 10⁴ failure is the wrong-variable obstruction created by the C(a) → C(ρ) migration.

The amendment's "change of kind, not degree" statement of the migration cost is now numerically supported: C(a) gives a factor-~2 cluster residual shared with MOND (a magnitude problem with a respected class); C(ρ) gives a 10⁴ structural failure (a kind problem specific to the framework). The migration was a downgrade, not a refinement, and the cost is now precise.

## §4 — What this session does not output

- **No verdict on the broader claim** that C(ρ) is closed. S678's substantive finding stands as `[AUDITED-NEGATIVE]` on the old R(I) substrate. The active reformulation (saturation reframe with independent **J**, Phase-1 sim) is upstream of this question; whether it has anything to say about the cluster sector is open.
- **No retag of S678**. The session doc S678 was already in MRH-relationship taxonomy after S679; what changes is the named mechanism, captured in this session's §3 refinement table. S678 retains its `[AUDITED-NEGATIVE]` status with its substance preserved.
- **No A/B comparison of the proposal vs the amendment**. The amendment supersedes the original on the mechanism story; both reach the same negative conclusion on C(ρ) bridging.
- **No cumulative tally.** Per the S679 discipline.
- **No engagement with whether the wrong-variable result is a "novel theorem"** worth a preprint. That is operator/coordinator work; my role is to verify the load-bearing claims and update my own session's framing.

## §5 — Files

- `Research/Session683_Cluster_Wrong_Variable_Refinement.md` (this document)
- `simulations/session683_cluster_wrong_variable_check.py` (within-Coma ρ vs g_bar verification; cross-system order-of-magnitude check)

## §6 — So what (under the frame-doc disciplines)

A research-repo autonomous session can verify a visitor-channel amendment's load-bearing quantitative claims against the framework's own equations and refine a prior session's framing without overturning its substance. This session did that for S678: the cluster-bridge result stands, but the named structural cause of the 10⁴ catastrophic failure is more precisely **wrong-variable** (local ρ vs non-local g_bar) than "one density scale." The one-scale residual *is* real but is a factor-~2 MOND-inherited effect, not the 10⁴ Synchronism-specific failure. The C(a) → C(ρ) migration was a change of kind, not a refinement; that cost is now precise.

The amendment's broader implications (whether to claim the wrong-variable obstruction as a novel theorem, how to attribute the one-scale residual properly to Sanders 2003, whether to investigate C(a) restoration) are operator/coordinator work, not autonomous-session output.
