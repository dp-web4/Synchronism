# Session 635: Cosmology Domain Scorecard — Sixth Site-Archive Audit

**Date**: 2026-04-26
**Type**: Targeted archive sweep responding to 2026-04-26 back-annotation
**Scope**: Step 3 of `cosmology_claim_status_after_audits.md` only. Steps 1, 2, 4 explicitly out of scope (see end).
**Grade**: A- (clean classifications, narrow scope acknowledged)

---

## Trigger

`Research/proposals/cosmology_claim_status_after_audits.md` (filed 2026-04-26). Visitor (researcher persona) flagged a **propagation failure**: Bullet Cluster structural failure is documented at `/key-claims` and `/dark-matter-failure`, but `/galaxy-rotation` still carries a "Strongly Supported" tile. Both pages share the same coherence-viscosity foundation; the front-page badge contradicts what the documentation page acknowledges.

The proposal asks for four steps. This session does Step 3 (audit remaining cosmology predictions for novelty) using the existing archive. Steps 1, 2, and 4 require larger scopes — see closing section.

## The Scorecard

Classification scheme from the proposal:
- **Novel-unfalsified**: derivation present, prediction not shared with MOND/ΛCDM/other frameworks, not yet refuted
- **Reparametrization**: prediction follows from MOND/ΛCDM/dimensional analysis without unique Synchronism content
- **Refuted**: contradicted by data
- **Unanchored**: prediction stated without derivation

| Claim | Class | Archive provenance |
|-------|-------|-------------------|
| Dark matter = high-viscosity coherence (CFD viscosity mapping) | **Refuted** | Bullet Cluster sign error (Harvey 2015: σ/m < 0.47 cm²/g requires effectively inviscid; framework predicts MORE viscous). Documented in SESSION_FOCUS Validation State and `/dark-matter-failure`. |
| a₀ = cH₀/(2π) (NP1) | **Reparametrization (artifact)** | S574/S461. P(chance match) = 56%; α free gives α≈1, not 0.5. Shared with McCulloch 2007, Verlinde 2017, Smolin 2017 — common dimensional structure. |
| NP2 type-dependent scatter | **Reparametrization** | S574 Test 1: vanishes after 6-var M/L correction (variance ratio 0.823, p=0.437). Standard M/L physics. |
| NP3 a₀ redshift evolution | **Untestable with SPARC** | S574 documented; needs ν(z) data unavailable in archive. |
| NP4 V-shaped scatter at g† | **Reparametrization** | S574 Test 2: reproduced by noise-only Monte Carlo (r=0.939). Standard MOND consequence of nonlinear interpolation. |
| C(ρ) = tanh(γ·log(ρ/ρ_crit + 1)) | **Reparametrization + framing error** | S574 Test 4: equivalent to MOND ν(x). S633: saturates within ~1.6 decades for any γ; site's "80 orders of magnitude" framing is wrong. |
| 6-variable offset model | **Reparametrization** | S574 Test 6 + S526 + S528: All 6 coefficients MOND-derivable. V-L ratio = 4.03 = MOND's 4.0. No Synchronism variable adds to residuals. |
| TEST-04 BAO ~10⁻⁴ density-dependent shift | **Unanchored** | Catalog states magnitude (10⁻⁴) without derivation. Proposal explicitly flags this gap. |
| TEST-07 ~500 Mpc cosmic interference | **Unanchored** | S632: derivation in Session #4 Track C is dimensionally inconsistent (m² ≠ m); 500 Mpc not derived. |
| TEST-08 SPARC environment scatter | **Pending** | R²=0.14 (weak signal); ΔBIC vs MOND-only not yet computed. Cannot classify until Step 2 runs. |
| TEST-01 TDG age-DM | **Untested** | Catalog Tier 1, no SPARC equivalent in archive. |
| TEST-02 UDG max DM | **Untested** | Catalog Tier 1; DF2/DF4 controversy is the natural test. |
| TEST-03 Compact elliptical min DM | **Untested** | Catalog Tier 1. |
| TEST-05 CMB cold spot–density | **Untested** | Catalog Tier 1. |
| TEST-06 Variable α_em (spatial) | **Untested** | Catalog Tier 1. Webb+Rosenband data exist; not analyzed in archive. |

### Counts

- **Refuted**: 1 (CFD viscosity dark matter)
- **Reparametrization**: 5 (a₀, NP2, NP4, C(ρ), 6-var offset model)
- **Unanchored**: 2 (BAO shift magnitude, 500 Mpc scale)
- **Pending**: 1 (TEST-08 ΔBIC)
- **Untested**: 5 (TDG, UDG, compact ellipticals, CMB cold spot, variable α)
- **Untestable with current data**: 1 (NP3 a₀ redshift)

**Novel-unfalsified: 0.**

The cosmology domain currently has no claim that is both (a) derived from framework principles and (b) genuinely novel and (c) supported by data. The Untested-Tier-1 set could change this if the experiments are run; until then, the domain's net positive evidence over MOND/ΛCDM is zero.

## Direct Disposition: The /galaxy-rotation Page

The proposal's most concrete claim: `/galaxy-rotation` carries "Strongly Supported" while `/key-claims` documents the Bullet Cluster structural failure. Three observations from the scorecard:

1. The /galaxy-rotation page's claim depends on the McGaugh 2016 RAR interpolating function. That function is **MOND**. The badge on /galaxy-rotation should reflect that the underlying observational fit is MOND's, with Synchronism contributing an environmental scatter ansatz (R²=0.14, weak).
2. The dark matter mechanism (CFD viscosity) is refuted. /galaxy-rotation should not display badges that imply that mechanism works.
3. The strongest defensible badge for /galaxy-rotation is **"MOND Reparametrization"** or — if Step 2's ΔBIC test is favorable — **"MOND + Environmental Scatter (ΔBIC=X)."** Currently the badge overclaims by one or two tiers.

This corresponds to the proposal's "Step 4 minimum outcome" and is operator-side site action.

## Step 1 (Gradient Instead of Level) — Brief Note

The proposal asks: does coherence GRADIENT (rather than level) provide a dark-matter mechanism consistent with Bullet Cluster's σ/m < 0.47 cm²/g?

Quick analytical check: in the existing CFD mapping, viscosity μ_eff is a function of local density (or local C). Replacing this with a function of |∇C| has no obvious physical motivation in the existing framework — gradient terms appear in surface tension and pressure gradients, but the dark-matter "stickiness" claim depends on local viscosity, not surface effects.

Even granting the substitution: in observed dark-matter halos, density falls smoothly with radius, so |∇C| is smooth and tracks C-level closely. The sign of the inferred viscosity would not flip just by switching from level to gradient unless the gradient-viscosity relation is monotonically opposite to the level-viscosity relation. That requires a specific ad-hoc functional form.

**Tentative**: the gradient hypothesis doesn't immediately solve the sign error. A serious reformulation (Step 1) would need a derivation that justifies the gradient form from framework principles — not just "swap one functional dependence for another." Out of scope for this session.

## What This Adds

**Six site-archive audits, all same failure mode:**

| Claim | Source | Failure |
|-------|--------|---------|
| BTFR n≈2.2 (S631) | #48 | "not rigorous" / refuted |
| α² in A (S631) | #66 | α=1.0 fiducial |
| 500 Mpc (S632) | #4 | dimensionally inconsistent |
| 80 orders (S633) | site only | range vs smoothness |
| 47 contributions (S634) | #582 says 30 | 57% overcount |
| /galaxy-rotation "Strongly Supported" (S635) | scorecard | Depends on refuted DM mechanism + uncomputed ΔBIC |

S635 is slightly different from S631–634: rather than tracing one specific claim to its archive source, it consolidates the existing archive verdicts on the cosmology domain and shows that the public site's domain-level framing isn't supported by the archive's claim-level honesty.

The pattern of "archive is honest, site overclaims" now has six instances. The audit channel has been productive 6-for-6.

## Out of Scope (For Operator)

- **Step 1 (DM gradient mechanism)**: requires derivation work, not just archive search. Tentative archive read says gradient doesn't trivially save the sign. A focused theoretical session could probe this but isn't a quick audit.
- **Step 2 (ΔBIC for SPARC environment vs MOND-only)**: requires running fits on the 14,760-galaxy sample. Multi-hour computation. Best done as a dedicated session (S606 / S610 are nearby but not the same comparison).
- **Step 4 (site updates)**: operator-side; site is reference-only for workers. Five existing pending items now plus this scorecard's recommendations.

## Recommended Site-Side Actions (Operator)

1. `/galaxy-rotation` badge → "MOND Reparametrization" (or "MOND + Environmental Scatter (pending ΔBIC)")
2. `/galaxy-rotation` page → add prominent note linking to /dark-matter-failure
3. Cosmology domain → add scorecard summary: 0 novel-unfalsified, 1 refuted, 5 reparametrizations, 2 unanchored, 5 untested
4. Move TEST-04 (BAO shift), TEST-07 (500 Mpc) to "Speculative" or "Exploratory" until magnitudes derive
5. Run Step 2 (ΔBIC) when bandwidth allows; this is the single highest-leverage open computation

## So What?

The cosmology domain's public-facing posture is two tiers stronger than the archive supports. Six audits in six days have surfaced the same pattern: claim-by-claim, the archive is honest; in domain-level framing, the site overclaims. The fix is propagation, not new science.

For my role: protocol holding. Silence on stale CBP firings, response to substantive new content. The audit channel has settled into a productive rhythm; the bottleneck is operator action on the now-six pending items.

## Files

- `Research/Session635_Cosmology_Scorecard.md` (this document)
- Updates to SESSION_FOCUS.md
